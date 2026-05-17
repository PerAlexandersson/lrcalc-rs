//! Schur product and skew Schur expansion.
//!
//! Product expansion realizes a product as the skew Schur function of a
//! disconnected skew shape.  Skew expansion uses the variable-content beta
//! tableau enumerator, which accumulates all output contents in one search.

use crate::lrcoef::{beta_lr_content_expansion, LrCoefError};
use crate::partition::Partition;

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct SchurTerm {
    pub partition: Vec<i32>,
    pub coefficient: u128,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum SchurExpansionError {
    InvalidPartition,
    ArithmeticOverflow,
}

impl From<LrCoefError> for SchurExpansionError {
    fn from(value: LrCoefError) -> Self {
        match value {
            LrCoefError::InvalidPartition => Self::InvalidPartition,
            LrCoefError::ArithmeticOverflow => Self::ArithmeticOverflow,
        }
    }
}

pub fn schur_product_expansion(
    sh1: &[i32],
    sh2: &[i32],
    rows: i32,
    cols: i32,
) -> Result<Vec<SchurTerm>, SchurExpansionError> {
    validate_partition(sh1)?;
    validate_partition(sh2)?;
    let sh1 = trim_trailing_zeroes(sh1);
    let sh2 = trim_trailing_zeroes(sh2);
    let size = checked_partition_size(&sh1)?
        .checked_add(checked_partition_size(&sh2)?)
        .ok_or(SchurExpansionError::ArithmeticOverflow)?;
    let (outer, inner) = disconnected_product_skew_shape(&sh1, &sh2)?;
    let max_part = col_bound(cols, size)?;
    let mut terms = schur_skew_expansion(&outer, &inner, rows)?;
    terms.retain(|term| part_entry_i32(&term.partition, 0) <= max_part);
    Ok(terms)
}

pub fn schur_skew_expansion(
    outer: &[i32],
    inner: &[i32],
    rows: i32,
) -> Result<Vec<SchurTerm>, SchurExpansionError> {
    validate_partition(outer)?;
    validate_partition(inner)?;
    let outer = trim_trailing_zeroes(outer);
    let inner = trim_trailing_zeroes(inner);
    if !contains_partition(&outer, &inner) {
        return Ok(Vec::new());
    }
    let Some(optimized) = optimize_skew_shape(&outer, &inner, rows)? else {
        return Ok(Vec::new());
    };
    let max_labels = if optimized.requested_rows >= 0 {
        Some(
            usize::try_from(optimized.requested_rows)
                .map_err(|_| SchurExpansionError::ArithmeticOverflow)?,
        )
    } else {
        None
    };
    let mut terms = beta_lr_content_expansion(
        &optimized.outer,
        &optimized.inner,
        &optimized.fixed_content,
        max_labels,
    )?
    .into_iter()
    .map(|term| SchurTerm {
        partition: add_content_vectors(&optimized.fixed_content, &term.content),
        coefficient: term.coefficient,
    })
    .collect::<Vec<_>>();
    terms.sort_by(|left, right| right.partition.cmp(&left.partition));

    Ok(terms)
}

#[derive(Clone, Debug, Eq, PartialEq)]
struct OptimizedSkewShape {
    outer: Vec<i32>,
    inner: Vec<i32>,
    fixed_content: Vec<i32>,
    requested_rows: i32,
}

fn validate_partition(partition: &[i32]) -> Result<(), SchurExpansionError> {
    Partition::from_slice(partition)
        .map(|_| ())
        .ok_or(SchurExpansionError::InvalidPartition)
}

pub(crate) fn trim_trailing_zeroes(partition: &[i32]) -> Vec<i32> {
    let len = partition
        .iter()
        .rposition(|&part| part != 0)
        .map_or(0, |index| index + 1);
    partition[..len].to_vec()
}

fn checked_partition_size(partition: &[i32]) -> Result<i32, SchurExpansionError> {
    let mut sum = 0i32;
    for &part in partition {
        sum = sum
            .checked_add(part)
            .ok_or(SchurExpansionError::ArithmeticOverflow)?;
    }
    Ok(sum)
}

#[cfg(test)]
fn row_bound(rows: i32, size: i32) -> Result<usize, SchurExpansionError> {
    if rows >= 0 {
        usize::try_from(rows).map_err(|_| SchurExpansionError::ArithmeticOverflow)
    } else {
        usize::try_from(size).map_err(|_| SchurExpansionError::ArithmeticOverflow)
    }
}

fn col_bound(cols: i32, size: i32) -> Result<i32, SchurExpansionError> {
    if cols >= 0 {
        Ok(cols)
    } else {
        Ok(size)
    }
}

fn disconnected_product_skew_shape(
    top: &[i32],
    bottom: &[i32],
) -> Result<(Vec<i32>, Vec<i32>), SchurExpansionError> {
    let shift = part_entry_i32(bottom, 0);
    let mut outer = Vec::with_capacity(top.len() + bottom.len());
    let mut inner = Vec::with_capacity(top.len());

    for &part in top {
        outer.push(
            shift
                .checked_add(part)
                .ok_or(SchurExpansionError::ArithmeticOverflow)?,
        );
        inner.push(shift);
    }
    outer.extend_from_slice(bottom);

    Ok((trim_trailing_zeroes(&outer), trim_trailing_zeroes(&inner)))
}

fn contains_partition(outer: &[i32], inner: &[i32]) -> bool {
    inner
        .iter()
        .enumerate()
        .all(|(index, &part)| part <= outer.get(index).copied().unwrap_or(0))
}

fn optimize_skew_shape(
    outer: &[i32],
    inner: &[i32],
    requested_rows: i32,
) -> Result<Option<OptimizedSkewShape>, SchurExpansionError> {
    if !contains_partition(outer, inner) {
        return Ok(None);
    }

    let mut row_bound = outer.len();
    let mut inner_len = inner.len().min(row_bound);
    if row_bound <= inner.len() {
        while row_bound > 0 && part_entry_i32(inner, row_bound - 1) == outer[row_bound - 1] {
            row_bound -= 1;
        }
        inner_len = row_bound;
    }
    let mut row_first = 0usize;
    while row_first < inner_len && part_entry_i32(inner, row_first) == outer[row_first] {
        row_first += 1;
    }
    let row_span = row_bound.saturating_sub(row_first);
    if row_span == 0 {
        return Ok(Some(OptimizedSkewShape {
            outer: Vec::new(),
            inner: Vec::new(),
            fixed_content: Vec::new(),
            requested_rows,
        }));
    }

    let slen = row_span
        .checked_mul(2)
        .ok_or(SchurExpansionError::ArithmeticOverflow)?;
    let max_rows = if requested_rows >= 0 {
        usize::try_from(requested_rows).map_err(|_| SchurExpansionError::ArithmeticOverflow)?
    } else {
        slen.checked_add(1)
            .ok_or(SchurExpansionError::ArithmeticOverflow)?
    };

    let mut out = vec![0; slen];
    let mut inn = vec![0; slen];
    let mut content = vec![0; row_span.max(max_rows).max(1)];
    let mut content_len = 0usize;
    let mut content_size = 0i32;
    let mut full_cols = 0i32;

    let mut c2 = outer[row_first];
    let mut r2_top = row_first;
    let mut r2_bot = row_first;
    while r2_bot < row_bound && c2 <= outer[r2_bot] {
        r2_bot += 1;
    }

    let mut r0_top = r2_top;
    let mut r0_bot = r2_bot;
    let mut partial = PartialSkewShape {
        inn: &mut inn,
        out: &mut out,
        rows: max_rows,
        top: 0,
        bot: 0,
        col: 0,
    };

    for c1 in (0..c2).rev() {
        let r1_top = r0_top;
        let r1_bot = r0_bot;
        if c1 == 0 {
            r0_top = row_bound;
            r0_bot = row_bound;
        }
        while r0_bot < row_bound && c1 <= outer[r0_bot] {
            r0_bot += 1;
        }
        while r0_top < inner.len() && c1 <= inner[r0_top] {
            r0_top += 1;
        }

        let r0_to_r1_span = r0_bot.saturating_sub(r1_top);
        let r1_height = r1_bot.saturating_sub(r1_top);

        if r0_top < r1_bot && r0_to_r1_span < max_rows {
            continue;
        }

        if c1 == c2 - 1 && r1_height > max_rows {
            return Ok(None);
        }

        if c1 == c2 - 1 && r1_height == max_rows {
            full_cols += 1;
            c2 = c1;
            r2_top = r0_top;
            r2_bot = r0_bot;
            continue;
        }

        let mut component_size = 0i32;
        for row in r2_top..r1_bot {
            let mut left = part_entry_i32(inner, row);
            if left < c1 {
                left = c1;
            }
            let mut right = outer[row];
            if right > c2 {
                right = c2;
            }
            component_size += right - left;
        }

        if (r1_top == r2_top || r1_bot == r2_bot)
            && 0 < content_size
            && content_size < component_size
        {
            let mut equal_rows = 1usize;
            let width = content[0];
            while equal_rows < content_len && content[equal_rows] == width {
                equal_rows += 1;
            }
            add_component(
                &mut partial,
                &content[..content_len],
                None,
                ComponentWindow {
                    c0: 0,
                    r0_top: 0,
                    r0_bot: content_len,
                    c1: width,
                    r1_top: 0,
                    r1_bot: equal_rows,
                },
            );
        }

        if r1_top == r2_top && component_size > content_size {
            content_len = r1_bot - r1_top;
            ensure_len(&mut content, content_len);
            for row in r1_top..r2_bot {
                content[row - r1_top] = c2 - c1;
            }
            for row in r2_bot..r1_bot {
                content[row - r1_top] = outer[row] - c1;
            }
            content_size = component_size;
        } else if r1_bot == r2_bot && component_size > content_size {
            content_len = r2_bot - r2_top;
            ensure_len(&mut content, content_len);
            for row in (r1_top..r2_bot).rev() {
                content[r2_bot - 1 - row] = c2 - c1;
            }
            for row in (r2_top..r1_top).rev() {
                content[r2_bot - 1 - row] = c2 - part_entry_i32(inner, row);
            }
            content_size = component_size;
        } else if component_size > 0 {
            add_component(
                &mut partial,
                outer,
                Some(inner),
                ComponentWindow {
                    c0: c1,
                    r0_top: r1_top,
                    r0_bot: r1_bot,
                    c1: c2,
                    r1_top: r2_top,
                    r1_bot: r2_bot,
                },
            );
        }

        c2 = c1;
        r2_top = r0_top;
        r2_bot = r0_bot;
    }

    if full_cols > 0 {
        ensure_len(&mut content, max_rows);
        for value in content.iter_mut().take(content_len) {
            *value += full_cols;
        }
        for value in content.iter_mut().take(max_rows).skip(content_len) {
            *value = full_cols;
        }
        content_len = max_rows;
    }

    let bot = partial.bot;
    let col_shift = partial.col;
    drop(partial);

    out.truncate(bot);
    inn.truncate(bot);
    for value in &mut out {
        *value -= col_shift;
    }
    for value in &mut inn {
        *value -= col_shift;
    }
    out = trim_trailing_zeroes(&out);
    inn = trim_trailing_zeroes(&inn);
    content.truncate(content_len);
    content = trim_trailing_zeroes(&content);

    Ok(Some(OptimizedSkewShape {
        outer: out,
        inner: inn,
        fixed_content: content,
        requested_rows,
    }))
}

struct PartialSkewShape<'a> {
    inn: &'a mut [i32],
    out: &'a mut [i32],
    rows: usize,
    top: usize,
    bot: usize,
    col: i32,
}

#[derive(Clone, Copy)]
struct ComponentWindow {
    c0: i32,
    r0_top: usize,
    r0_bot: usize,
    c1: i32,
    r1_top: usize,
    r1_bot: usize,
}

fn add_component(
    partial: &mut PartialSkewShape<'_>,
    out0: &[i32],
    inn0: Option<&[i32]>,
    window: ComponentWindow,
) {
    let mut x = partial.top + partial.rows + window.r1_top - window.r1_bot;
    if x > partial.bot {
        x = partial.bot;
    }
    let y1 = x + window.r1_bot - window.r1_top;
    let z = y1 + window.r0_bot - window.r1_bot;

    for row in partial.bot..y1 {
        partial.out[row] = partial.col;
    }
    for row in y1..z {
        let source_row = row - x + window.r1_top;
        let c = part_entry_i32(out0, source_row);
        partial.out[row] = partial.col + c - window.c1;
    }

    let y0 = x + window.r0_top - window.r1_top;
    for row in x..y0 {
        let source_row = row - x + window.r1_top;
        let c = inn0.map_or(0, |inner| part_entry_i32(inner, source_row));
        partial.inn[row] = partial.col + c - window.c1;
    }
    for row in y0..z {
        partial.inn[row] = partial.col - window.c1 + window.c0;
    }

    partial.col -= window.c1 - window.c0;
    partial.top = y0;
    partial.bot = z;
}

fn ensure_len(values: &mut Vec<i32>, len: usize) {
    if values.len() < len {
        values.resize(len, 0);
    }
}

fn add_content_vectors(left: &[i32], right: &[i32]) -> Vec<i32> {
    let len = left.len().max(right.len());
    let mut out = Vec::with_capacity(len);
    for index in 0..len {
        out.push(part_entry_i32(left, index) + part_entry_i32(right, index));
    }
    trim_trailing_zeroes(&out)
}

fn part_entry_i32(partition: &[i32], index: usize) -> i32 {
    partition.get(index).copied().unwrap_or(0)
}

#[cfg(test)]
fn visit_partitions<F>(
    remaining: i32,
    rows_left: usize,
    max_part: i32,
    current: &mut Vec<i32>,
    visit: &mut F,
) -> Result<(), SchurExpansionError>
where
    F: FnMut(&[i32]) -> Result<(), SchurExpansionError>,
{
    if remaining == 0 {
        return visit(current);
    }
    if rows_left == 0 || max_part == 0 {
        return Ok(());
    }
    let rows_left_i32 =
        i32::try_from(rows_left).map_err(|_| SchurExpansionError::ArithmeticOverflow)?;
    if remaining > rows_left_i32.saturating_mul(max_part) {
        return Ok(());
    }

    let upper = remaining.min(max_part);
    for part in (1..=upper).rev() {
        current.push(part);
        visit_partitions(remaining - part, rows_left - 1, part, current, visit)?;
        current.pop();
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::lrcoef::lrcoef;

    fn term_map(terms: Vec<SchurTerm>) -> Vec<(Vec<i32>, u128)> {
        terms
            .into_iter()
            .map(|term| (term.partition, term.coefficient))
            .collect()
    }

    fn sorted_term_map(terms: Vec<SchurTerm>) -> Vec<(Vec<i32>, u128)> {
        let mut terms = term_map(terms);
        terms.sort();
        terms
    }

    fn scalar_skew_expansion(
        outer: &[i32],
        inner: &[i32],
        rows: i32,
    ) -> Result<Vec<(Vec<i32>, u128)>, SchurExpansionError> {
        if !contains_partition(outer, inner) {
            return Ok(Vec::new());
        }
        let size = checked_partition_size(outer)?
            .checked_sub(checked_partition_size(inner)?)
            .ok_or(SchurExpansionError::ArithmeticOverflow)?;
        let max_rows = row_bound(rows, size)?;
        let mut terms = Vec::new();
        visit_partitions(size, max_rows, size, &mut Vec::new(), &mut |content| {
            let coefficient = lrcoef(outer, inner, content)?;
            if coefficient != 0 {
                terms.push((content.to_vec(), coefficient));
            }
            Ok(())
        })?;
        terms.sort();
        Ok(terms)
    }

    fn scalar_product_expansion(
        sh1: &[i32],
        sh2: &[i32],
        rows: i32,
        cols: i32,
    ) -> Result<Vec<(Vec<i32>, u128)>, SchurExpansionError> {
        let size = checked_partition_size(sh1)?
            .checked_add(checked_partition_size(sh2)?)
            .ok_or(SchurExpansionError::ArithmeticOverflow)?;
        let max_rows = row_bound(rows, size)?;
        let max_part = col_bound(cols, size)?;
        let mut terms = Vec::new();

        visit_partitions(size, max_rows, max_part, &mut Vec::new(), &mut |lambda| {
            if !contains_partition(lambda, sh1) || !contains_partition(lambda, sh2) {
                return Ok(());
            }
            let coefficient = lrcoef(lambda, sh1, sh2)?;
            if coefficient != 0 {
                terms.push((lambda.to_vec(), coefficient));
            }
            Ok(())
        })?;
        terms.sort();
        Ok(terms)
    }

    #[test]
    fn product_of_single_boxes() {
        let terms = schur_product_expansion(&[1], &[1], -1, -1).unwrap();
        assert_eq!(term_map(terms), vec![(vec![2], 1), (vec![1, 1], 1)]);
    }

    #[test]
    fn product_respects_row_and_column_bounds() {
        let row_terms = schur_product_expansion(&[1], &[1], 1, -1).unwrap();
        assert_eq!(term_map(row_terms), vec![(vec![2], 1)]);

        let col_terms = schur_product_expansion(&[1], &[1], -1, 1).unwrap();
        assert_eq!(term_map(col_terms), vec![(vec![1, 1], 1)]);
    }

    #[test]
    fn product_respects_larger_row_and_column_bounds() {
        let row_terms = schur_product_expansion(&[2, 1], &[2, 1], 3, -1).unwrap();
        assert_eq!(
            sorted_term_map(row_terms),
            vec![
                (vec![2, 2, 2], 1),
                (vec![3, 2, 1], 2),
                (vec![3, 3], 1),
                (vec![4, 1, 1], 1),
                (vec![4, 2], 1),
            ]
        );

        let col_terms = schur_product_expansion(&[2, 1], &[2, 1], -1, 2).unwrap();
        assert_eq!(
            sorted_term_map(col_terms),
            vec![(vec![2, 2, 1, 1], 1), (vec![2, 2, 2], 1)]
        );
    }

    #[test]
    fn product_of_two_hooks_matches_upstream_example() {
        let terms = schur_product_expansion(&[2, 1], &[2, 1], -1, -1).unwrap();
        assert_eq!(
            sorted_term_map(terms),
            vec![
                (vec![2, 2, 1, 1], 1),
                (vec![2, 2, 2], 1),
                (vec![3, 1, 1, 1], 1),
                (vec![3, 2, 1], 2),
                (vec![3, 3], 1),
                (vec![4, 1, 1], 1),
                (vec![4, 2], 1),
            ]
        );
    }

    #[test]
    fn product_with_empty_partition_is_identity() {
        let terms = schur_product_expansion(&[3, 1], &[], -1, -1).unwrap();
        assert_eq!(term_map(terms), vec![(vec![3, 1], 1)]);
    }

    #[test]
    fn product_with_zero_row_bound_is_empty_unless_size_zero() {
        let nonempty = schur_product_expansion(&[1], &[1], 0, -1).unwrap();
        assert!(nonempty.is_empty());

        let empty = schur_product_expansion(&[], &[], 0, -1).unwrap();
        assert_eq!(term_map(empty), vec![(vec![], 1)]);
    }

    #[test]
    fn product_trims_trailing_zeroes() {
        let terms = schur_product_expansion(&[1, 0, 0], &[1, 0], -1, -1).unwrap();
        assert_eq!(term_map(terms), vec![(vec![2], 1), (vec![1, 1], 1)]);
    }

    #[test]
    fn product_rejects_invalid_partitions() {
        assert_eq!(
            schur_product_expansion(&[1, 2], &[1], -1, -1),
            Err(SchurExpansionError::InvalidPartition)
        );
        assert_eq!(
            schur_product_expansion(&[1], &[1, -1], -1, -1),
            Err(SchurExpansionError::InvalidPartition)
        );
    }

    #[test]
    fn disconnected_product_matches_scalar_for_small_shapes() {
        for size1 in 0..=5 {
            for sh1 in partitions_of(size1) {
                for size2 in 0..=5 {
                    for sh2 in partitions_of(size2) {
                        for rows in [-1, 0, 1, 2, 3, 4] {
                            for cols in [-1, 1, 2, 3, 4] {
                                assert_eq!(
                                    sorted_term_map(
                                        schur_product_expansion(&sh1, &sh2, rows, cols).unwrap()
                                    ),
                                    scalar_product_expansion(&sh1, &sh2, rows, cols).unwrap(),
                                    "sh1={sh1:?} sh2={sh2:?} rows={rows} cols={cols}"
                                );
                            }
                        }
                    }
                }
            }
        }
    }

    #[test]
    fn skew_two_cell_hook() {
        let terms = schur_skew_expansion(&[2, 1], &[1], -1).unwrap();
        assert_eq!(term_map(terms), vec![(vec![2], 1), (vec![1, 1], 1)]);
    }

    #[test]
    fn skew_respects_row_bound() {
        let terms = schur_skew_expansion(&[2, 1], &[1], 1).unwrap();
        assert_eq!(term_map(terms), vec![(vec![2], 1)]);
    }

    #[test]
    fn skew_larger_example_matches_upstream_example() {
        let terms = schur_skew_expansion(&[3, 2, 1], &[2, 1], -1).unwrap();
        assert_eq!(
            term_map(terms),
            vec![(vec![3], 1), (vec![2, 1], 2), (vec![1, 1, 1], 1)]
        );
    }

    #[test]
    fn skew_row_bound_filters_long_content_partitions() {
        let terms = schur_skew_expansion(&[3, 2, 1], &[2, 1], 2).unwrap();
        assert_eq!(term_map(terms), vec![(vec![3], 1), (vec![2, 1], 2)]);
    }

    #[test]
    fn skew_with_empty_inner_is_identity() {
        let terms = schur_skew_expansion(&[3, 1], &[], -1).unwrap();
        assert_eq!(term_map(terms), vec![(vec![3, 1], 1)]);
    }

    #[test]
    fn skew_with_equal_shapes_is_empty_partition_term() {
        let terms = schur_skew_expansion(&[3, 1], &[3, 1], -1).unwrap();
        assert_eq!(term_map(terms), vec![(vec![], 1)]);
    }

    #[test]
    fn skew_requires_inner_contained_in_outer() {
        let terms = schur_skew_expansion(&[2], &[2, 1], -1).unwrap();
        assert!(terms.is_empty());
    }

    #[test]
    fn skew_rejects_invalid_partitions() {
        assert_eq!(
            schur_skew_expansion(&[1, 2], &[1], -1),
            Err(SchurExpansionError::InvalidPartition)
        );
        assert_eq!(
            schur_skew_expansion(&[2], &[-1], -1),
            Err(SchurExpansionError::InvalidPartition)
        );
    }

    #[test]
    fn optimized_skew_matches_scalar_for_small_shapes() {
        for outer_size in 0..=7 {
            for outer in partitions_of(outer_size) {
                for inner_size in 0..=outer_size {
                    for inner in partitions_of(inner_size) {
                        if !contains_partition(&outer, &inner) {
                            continue;
                        }
                        for rows in [-1, 0, 1, 2, 3, 4] {
                            assert_eq!(
                                sorted_term_map(
                                    schur_skew_expansion(&outer, &inner, rows).unwrap()
                                ),
                                scalar_skew_expansion(&outer, &inner, rows).unwrap(),
                                "outer={outer:?} inner={inner:?} rows={rows}"
                            );
                        }
                    }
                }
            }
        }
    }

    fn partitions_of(n: i32) -> Vec<Vec<i32>> {
        fn go(remaining: i32, max_part: i32, current: &mut Vec<i32>, out: &mut Vec<Vec<i32>>) {
            if remaining == 0 {
                out.push(current.clone());
                return;
            }
            for part in (1..=remaining.min(max_part)).rev() {
                current.push(part);
                go(remaining - part, part, current, out);
                current.pop();
            }
        }

        let mut out = Vec::new();
        go(n, n, &mut Vec::new(), &mut out);
        out
    }
}
