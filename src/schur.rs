//! Schur product and skew Schur expansion.
//!
//! Product expansion realizes a product as the skew Schur function of a
//! disconnected skew shape.  Skew expansion uses the variable-content beta
//! tableau enumerator, which accumulates all output contents in one search.

use crate::lrcoef::{compact_skew_shape, visit_beta_lr_content_expansion_with_len, LrCoefError};
use crate::partition::Partition;
use std::collections::BTreeMap;

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct SchurTerm {
    pub partition: Vec<i32>,
    pub coefficient: u128,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct SignedSchurTerm {
    pub partition: Vec<i32>,
    pub coefficient: i128,
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
    let mut terms = Vec::new();
    visit_schur_skew_expansion(outer, inner, rows, |partition, coefficient| {
        terms.push(SchurTerm {
            partition,
            coefficient,
        });
    })?;
    terms.sort_by(|left, right| right.partition.cmp(&left.partition));

    Ok(terms)
}

pub(crate) fn visit_schur_skew_expansion<F>(
    outer: &[i32],
    inner: &[i32],
    rows: i32,
    visit: F,
) -> Result<(), SchurExpansionError>
where
    F: FnMut(Vec<i32>, u128),
{
    visit_schur_skew_expansion_with_len(outer, inner, rows, |_| {}, visit)
}

pub(crate) fn visit_schur_skew_expansion_with_len<B, F>(
    outer: &[i32],
    inner: &[i32],
    rows: i32,
    begin: B,
    mut visit: F,
) -> Result<(), SchurExpansionError>
where
    B: FnOnce(usize),
    F: FnMut(Vec<i32>, u128),
{
    validate_partition(outer)?;
    validate_partition(inner)?;
    let outer = trim_trailing_zeroes(outer);
    let inner = trim_trailing_zeroes(inner);
    if !contains_partition(&outer, &inner) {
        begin(0);
        return Ok(());
    }
    let (outer, inner) = compact_skew_shape(outer, inner);
    let Some(optimized) = optimize_skew_shape(&outer, &inner, rows)? else {
        begin(0);
        return Ok(());
    };
    let max_labels = if optimized.requested_rows >= 0 {
        Some(
            usize::try_from(optimized.requested_rows)
                .map_err(|_| SchurExpansionError::ArithmeticOverflow)?,
        )
    } else {
        None
    };
    let mut overflowed = false;
    visit_beta_lr_content_expansion_with_len(
        &optimized.outer,
        &optimized.inner,
        &optimized.fixed_content,
        max_labels,
        begin,
        |content, coefficient| {
            if overflowed {
                return;
            }
            let Some(partition) = add_content_vectors(&optimized.fixed_content, content) else {
                overflowed = true;
                return;
            };
            visit(partition, coefficient);
        },
    )?;
    if overflowed {
        return Err(SchurExpansionError::ArithmeticOverflow);
    }
    Ok(())
}

pub fn schur_coproduct_expansion(
    shape: &[i32],
    rows: i32,
    cols: i32,
    all: bool,
) -> Result<Vec<SchurTerm>, SchurExpansionError> {
    validate_partition(shape)?;
    if rows < 0 || cols < 0 {
        return Err(SchurExpansionError::InvalidPartition);
    }
    let rows = usize::try_from(rows).map_err(|_| SchurExpansionError::ArithmeticOverflow)?;
    let shape = trim_trailing_zeroes(shape);
    let rectangle = vec![cols; rows];
    let mut terms = schur_product_expansion(&shape, &rectangle, -1, -1)?;
    if !all {
        let mut filtered = Vec::with_capacity(terms.len());
        for term in terms {
            if !coproduct_is_redundant(&term.partition, rows, cols)? {
                filtered.push(term);
            }
        }
        terms = filtered;
    }
    Ok(terms)
}

pub fn fusion_reduce_partition(
    partition: &[i32],
    rows: usize,
    level: i32,
) -> Result<Option<(Vec<i32>, i32)>, SchurExpansionError> {
    validate_partition(partition)?;
    if partition
        .iter()
        .enumerate()
        .skip(rows)
        .any(|(_, &part)| part != 0)
    {
        return Ok(None);
    }

    let mut values = vec![0; rows];
    for (index, &part) in partition.iter().take(rows).enumerate() {
        values[index] = part;
    }
    let Some((reduced, sign)) = fusion_reduce_values(&values, level)? else {
        return Ok(None);
    };
    Ok(Some((trim_trailing_zeroes(&reduced), sign)))
}

pub fn schur_product_fusion_expansion(
    sh1: &[i32],
    sh2: &[i32],
    rows: i32,
    level: i32,
) -> Result<Vec<SignedSchurTerm>, SchurExpansionError> {
    validate_partition(sh1)?;
    validate_partition(sh2)?;
    if rows < 0 || level < 0 {
        return Err(SchurExpansionError::InvalidPartition);
    }
    if rows == 0 {
        if trim_trailing_zeroes(sh1).is_empty() && trim_trailing_zeroes(sh2).is_empty() {
            return Ok(vec![SignedSchurTerm {
                partition: Vec::new(),
                coefficient: 1,
            }]);
        }
        return Ok(Vec::new());
    }
    let rows = usize::try_from(rows).map_err(|_| SchurExpansionError::ArithmeticOverflow)?;
    if has_nonzero_entry_at_or_after(sh1, rows) || has_nonzero_entry_at_or_after(sh2, rows) {
        return Ok(Vec::new());
    }

    let ordinary_terms = schur_product_expansion(sh1, sh2, rows as i32, -1)?;
    let mut merged = BTreeMap::<Vec<i32>, i128>::new();
    for term in ordinary_terms {
        let Some((partition, sign)) = fusion_reduce_partition(&term.partition, rows, level)? else {
            continue;
        };
        let coefficient = i128::try_from(term.coefficient)
            .map_err(|_| SchurExpansionError::ArithmeticOverflow)?
            .checked_mul(i128::from(sign))
            .ok_or(SchurExpansionError::ArithmeticOverflow)?;
        let entry = merged.entry(partition.clone()).or_insert(0);
        *entry = entry
            .checked_add(coefficient)
            .ok_or(SchurExpansionError::ArithmeticOverflow)?;
        if *entry == 0 {
            merged.remove(&partition);
        }
    }

    let mut terms = merged
        .into_iter()
        .map(|(partition, coefficient)| SignedSchurTerm {
            partition,
            coefficient,
        })
        .collect::<Vec<_>>();
    terms.sort_by(|left, right| right.partition.cmp(&left.partition));
    Ok(terms)
}

pub(crate) fn fusion_reduce_values(
    values: &[i32],
    level: i32,
) -> Result<Option<(Vec<i32>, i32)>, SchurExpansionError> {
    let rows = values.len();
    if rows == 0 || level < 0 {
        return Err(SchurExpansionError::InvalidPartition);
    }
    let rows_i64 = i64::try_from(rows).map_err(|_| SchurExpansionError::ArithmeticOverflow)?;
    let level_i64 = i64::from(level);
    let n = rows_i64
        .checked_add(level_i64)
        .ok_or(SchurExpansionError::ArithmeticOverflow)?;
    if n <= 0 {
        return Err(SchurExpansionError::InvalidPartition);
    }

    let mut q = 0i64;
    let mut tmp = Vec::with_capacity(rows);
    for (index, &value) in values.iter().enumerate() {
        let index_i64 =
            i64::try_from(index).map_err(|_| SchurExpansionError::ArithmeticOverflow)?;
        let a = i64::from(value)
            .checked_add(rows_i64)
            .and_then(|a| a.checked_sub(index_i64))
            .and_then(|a| a.checked_sub(1))
            .ok_or(SchurExpansionError::ArithmeticOverflow)?;
        let b = floor_div_i64(a, n);
        q = q
            .checked_add(b)
            .ok_or(SchurExpansionError::ArithmeticOverflow)?;
        let shifted = a
            .checked_sub(
                b.checked_mul(n)
                    .ok_or(SchurExpansionError::ArithmeticOverflow)?,
            )
            .and_then(|a| a.checked_sub(rows_i64))
            .and_then(|a| a.checked_add(1))
            .ok_or(SchurExpansionError::ArithmeticOverflow)?;
        tmp.push(shifted);
    }

    let mut sign_parity = if rows % 2 == 1 { 0 } else { q.rem_euclid(2) };
    for index in 0..rows.saturating_sub(1) {
        let mut max_index = index;
        let mut max_value = tmp[max_index];
        for next in index + 1..rows {
            if max_value < tmp[next] {
                max_index = next;
                max_value = tmp[max_index];
            }
        }
        if max_index != index {
            tmp[max_index] = tmp[index];
            tmp[index] = max_value;
            sign_parity ^= 1;
        }
    }

    let mut reduced = vec![0; rows];
    for index in 0..rows {
        if index > 0 && tmp[index - 1] == tmp[index] {
            return Ok(None);
        }
        let index_i64 =
            i64::try_from(index).map_err(|_| SchurExpansionError::ArithmeticOverflow)?;
        let k = index_i64
            .checked_add(q)
            .ok_or(SchurExpansionError::ArithmeticOverflow)?;
        let quotient = k / rows_i64;
        let a = tmp[index]
            .checked_add(k)
            .and_then(|a| a.checked_add(quotient.checked_mul(level_i64)?))
            .ok_or(SchurExpansionError::ArithmeticOverflow)?;
        let target = usize::try_from(k.rem_euclid(rows_i64))
            .map_err(|_| SchurExpansionError::ArithmeticOverflow)?;
        reduced[target] = i32::try_from(a).map_err(|_| SchurExpansionError::ArithmeticOverflow)?;
    }

    let sign = if sign_parity == 0 { 1 } else { -1 };
    Ok(Some((reduced, sign)))
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub(crate) struct OptimizedSkewShape {
    pub(crate) outer: Vec<i32>,
    pub(crate) inner: Vec<i32>,
    pub(crate) fixed_content: Vec<i32>,
    pub(crate) requested_rows: i32,
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

fn has_nonzero_entry_at_or_after(partition: &[i32], index: usize) -> bool {
    partition.iter().skip(index).any(|&part| part != 0)
}

fn floor_div_i64(numerator: i64, denominator: i64) -> i64 {
    numerator.div_euclid(denominator)
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

fn coproduct_is_redundant(
    partition: &[i32],
    rows: usize,
    cols: i32,
) -> Result<bool, SchurExpansionError> {
    let rows_i64 = i64::try_from(rows).map_err(|_| SchurExpansionError::ArithmeticOverflow)?;
    let cols_i64 = i64::from(cols);
    let mut left_size = rows_i64
        .checked_mul(cols_i64)
        .and_then(|value| value.checked_neg())
        .ok_or(SchurExpansionError::ArithmeticOverflow)?;
    for row in 0..rows {
        left_size = left_size
            .checked_add(i64::from(part_entry_i32(partition, row)))
            .ok_or(SchurExpansionError::ArithmeticOverflow)?;
    }

    let mut right_size = 0i64;
    for &part in partition.iter().skip(rows) {
        right_size = right_size
            .checked_add(i64::from(part))
            .ok_or(SchurExpansionError::ArithmeticOverflow)?;
    }

    if left_size != right_size {
        return Ok(left_size < right_size);
    }

    for row in 0..rows {
        let diff = i64::from(part_entry_i32(partition, row))
            .checked_sub(cols_i64)
            .and_then(|value| value.checked_sub(i64::from(part_entry_i32(partition, rows + row))))
            .ok_or(SchurExpansionError::ArithmeticOverflow)?;
        if diff != 0 {
            return Ok(diff > 0);
        }
    }
    Ok(false)
}

fn contains_partition(outer: &[i32], inner: &[i32]) -> bool {
    inner
        .iter()
        .enumerate()
        .all(|(index, &part)| part <= outer.get(index).copied().unwrap_or(0))
}

pub(crate) fn optimize_skew_shape(
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
    let mut content_size = 0i64;
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
            full_cols = full_cols
                .checked_add(1)
                .ok_or(SchurExpansionError::ArithmeticOverflow)?;
            c2 = c1;
            r2_top = r0_top;
            r2_bot = r0_bot;
            continue;
        }

        let mut component_size = 0i64;
        for (row, &outer_row) in outer.iter().enumerate().take(r1_bot).skip(r2_top) {
            let mut left = part_entry_i32(inner, row);
            if left < c1 {
                left = c1;
            }
            let mut right = outer_row;
            if right > c2 {
                right = c2;
            }
            component_size = component_size
                .checked_add(i64::from(right) - i64::from(left))
                .ok_or(SchurExpansionError::ArithmeticOverflow)?;
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
            )?;
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
            )?;
        }

        c2 = c1;
        r2_top = r0_top;
        r2_bot = r0_bot;
    }

    if full_cols > 0 {
        ensure_len(&mut content, max_rows);
        for value in content.iter_mut().take(content_len) {
            *value = value
                .checked_add(full_cols)
                .ok_or(SchurExpansionError::ArithmeticOverflow)?;
        }
        for value in content.iter_mut().take(max_rows).skip(content_len) {
            *value = full_cols;
        }
        content_len = max_rows;
    }

    let bot = partial.bot;
    let col_shift = partial.col;

    out.truncate(bot);
    inn.truncate(bot);
    for value in &mut out {
        *value = value
            .checked_sub(col_shift)
            .ok_or(SchurExpansionError::ArithmeticOverflow)?;
    }
    for value in &mut inn {
        *value = value
            .checked_sub(col_shift)
            .ok_or(SchurExpansionError::ArithmeticOverflow)?;
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
) -> Result<(), SchurExpansionError> {
    let mut x = partial
        .top
        .checked_add(partial.rows)
        .and_then(|value| value.checked_add(window.r1_top))
        .and_then(|value| value.checked_sub(window.r1_bot))
        .ok_or(SchurExpansionError::ArithmeticOverflow)?;
    if x > partial.bot {
        x = partial.bot;
    }
    let y1 = x
        .checked_add(window.r1_bot)
        .and_then(|value| value.checked_sub(window.r1_top))
        .ok_or(SchurExpansionError::ArithmeticOverflow)?;
    let z = y1
        .checked_add(window.r0_bot)
        .and_then(|value| value.checked_sub(window.r1_bot))
        .ok_or(SchurExpansionError::ArithmeticOverflow)?;

    for row in partial.bot..y1 {
        partial.out[row] = partial.col;
    }
    for row in y1..z {
        let source_row = row - x + window.r1_top;
        let c = part_entry_i32(out0, source_row);
        partial.out[row] = partial
            .col
            .checked_add(c)
            .and_then(|value| value.checked_sub(window.c1))
            .ok_or(SchurExpansionError::ArithmeticOverflow)?;
    }

    let y0 = x + window.r0_top - window.r1_top;
    for row in x..y0 {
        let source_row = row - x + window.r1_top;
        let c = inn0.map_or(0, |inner| part_entry_i32(inner, source_row));
        partial.inn[row] = partial
            .col
            .checked_add(c)
            .and_then(|value| value.checked_sub(window.c1))
            .ok_or(SchurExpansionError::ArithmeticOverflow)?;
    }
    for row in y0..z {
        partial.inn[row] = partial
            .col
            .checked_sub(window.c1)
            .and_then(|value| value.checked_add(window.c0))
            .ok_or(SchurExpansionError::ArithmeticOverflow)?;
    }

    let delta = window
        .c1
        .checked_sub(window.c0)
        .ok_or(SchurExpansionError::ArithmeticOverflow)?;
    partial.col = partial
        .col
        .checked_sub(delta)
        .ok_or(SchurExpansionError::ArithmeticOverflow)?;
    partial.top = y0;
    partial.bot = z;
    Ok(())
}

fn ensure_len(values: &mut Vec<i32>, len: usize) {
    if values.len() < len {
        values.resize(len, 0);
    }
}

fn add_content_vectors(left: &[i32], right: &[i32]) -> Option<Vec<i32>> {
    let len = left.len().max(right.len());
    let mut out = Vec::with_capacity(len);
    for index in 0..len {
        out.push(part_entry_i32(left, index).checked_add(part_entry_i32(right, index))?);
    }
    Some(trim_trailing_zeroes(&out))
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
    use std::collections::BTreeMap;

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

    fn signed_term_map(terms: Vec<SignedSchurTerm>) -> Vec<(Vec<i32>, i128)> {
        terms
            .into_iter()
            .map(|term| (term.partition, term.coefficient))
            .collect()
    }

    fn sorted_signed_term_map(terms: Vec<SignedSchurTerm>) -> Vec<(Vec<i32>, i128)> {
        let mut terms = signed_term_map(terms);
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

    fn reference_coproduct_expansion(
        shape: &[i32],
        rows: i32,
        cols: i32,
        all: bool,
    ) -> Result<Vec<(Vec<i32>, u128)>, SchurExpansionError> {
        let rows_usize =
            usize::try_from(rows).map_err(|_| SchurExpansionError::ArithmeticOverflow)?;
        let rectangle = vec![cols; rows_usize];
        let mut terms = scalar_product_expansion(shape, &rectangle, -1, -1)?;
        if !all {
            let mut filtered = Vec::new();
            for (partition, coefficient) in terms {
                if !coproduct_is_redundant(&partition, rows_usize, cols)? {
                    filtered.push((partition, coefficient));
                }
            }
            terms = filtered;
        }
        terms.sort();
        Ok(terms)
    }

    fn scalar_fusion_expansion(
        sh1: &[i32],
        sh2: &[i32],
        rows: i32,
        level: i32,
    ) -> Result<Vec<(Vec<i32>, i128)>, SchurExpansionError> {
        let rows_usize =
            usize::try_from(rows).map_err(|_| SchurExpansionError::ArithmeticOverflow)?;
        let mut merged = BTreeMap::<Vec<i32>, i128>::new();
        for (partition, coefficient) in scalar_product_expansion(sh1, sh2, rows, -1)? {
            let Some((reduced, sign)) = fusion_reduce_partition(&partition, rows_usize, level)?
            else {
                continue;
            };
            let coefficient = i128::try_from(coefficient)
                .map_err(|_| SchurExpansionError::ArithmeticOverflow)?
                * i128::from(sign);
            let entry = merged.entry(reduced.clone()).or_insert(0);
            *entry += coefficient;
            if *entry == 0 {
                merged.remove(&reduced);
            }
        }
        Ok(merged.into_iter().collect())
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
    fn fusion_reduce_matches_upstream_spot_checks() {
        assert_eq!(fusion_reduce_partition(&[2], 2, 1).unwrap(), None);
        assert_eq!(
            fusion_reduce_partition(&[1, 1], 2, 1).unwrap(),
            Some((vec![1, 1], 1))
        );
        assert_eq!(
            fusion_reduce_partition(&[2], 2, 0).unwrap(),
            Some((vec![1, 1], -1))
        );
    }

    #[test]
    fn fusion_product_matches_upstream_examples() {
        let terms = schur_product_fusion_expansion(&[1], &[1], 2, 1).unwrap();
        assert_eq!(sorted_signed_term_map(terms), vec![(vec![1, 1], 1)]);

        let terms = schur_product_fusion_expansion(&[1], &[1], 2, 2).unwrap();
        assert_eq!(
            sorted_signed_term_map(terms),
            vec![(vec![1, 1], 1), (vec![2], 1)]
        );

        let terms = schur_product_fusion_expansion(&[2, 1], &[2, 1], 3, 2).unwrap();
        assert_eq!(
            sorted_signed_term_map(terms),
            vec![(vec![2, 2, 2], 1), (vec![3, 2, 1], 1)]
        );
    }

    #[test]
    fn fusion_product_matches_scalar_reference_for_small_shapes() {
        for size1 in 0..=4 {
            for sh1 in partitions_of(size1) {
                for size2 in 0..=4 {
                    for sh2 in partitions_of(size2) {
                        for rows in 1..=4 {
                            for level in 0..=4 {
                                assert_eq!(
                                    sorted_signed_term_map(
                                        schur_product_fusion_expansion(&sh1, &sh2, rows, level)
                                            .unwrap()
                                    ),
                                    scalar_fusion_expansion(&sh1, &sh2, rows, level).unwrap(),
                                    "sh1={sh1:?} sh2={sh2:?} rows={rows} level={level}"
                                );
                            }
                        }
                    }
                }
            }
        }
    }

    #[test]
    fn fusion_product_handles_invalid_or_too_long_inputs() {
        assert_eq!(
            schur_product_fusion_expansion(&[1], &[1], 0, 1),
            Ok(Vec::new())
        );
        assert_eq!(
            schur_product_fusion_expansion(&[0], &[0], 0, 1),
            Ok(vec![SignedSchurTerm {
                partition: Vec::new(),
                coefficient: 1,
            }])
        );
        assert_eq!(
            schur_product_fusion_expansion(&[1], &[1], 2, -1),
            Err(SchurExpansionError::InvalidPartition)
        );
        assert!(schur_product_fusion_expansion(&[1, 1], &[1], 1, 2)
            .unwrap()
            .is_empty());
    }

    #[test]
    fn coproduct_filters_redundant_terms() {
        let reduced = schur_coproduct_expansion(&[1], 1, 1, false).unwrap();
        assert_eq!(term_map(reduced), vec![(vec![2], 1)]);

        let all = schur_coproduct_expansion(&[1], 1, 1, true).unwrap();
        assert_eq!(term_map(all), vec![(vec![2], 1), (vec![1, 1], 1)]);
    }

    #[test]
    fn coproduct_matches_reference_for_small_shapes() {
        for size in 0..=6 {
            for shape in partitions_of(size) {
                for rows in 0..=4 {
                    for cols in 0..=4 {
                        for all in [false, true] {
                            assert_eq!(
                                sorted_term_map(
                                    schur_coproduct_expansion(&shape, rows, cols, all).unwrap()
                                ),
                                reference_coproduct_expansion(&shape, rows, cols, all).unwrap(),
                                "shape={shape:?} rows={rows} cols={cols} all={all}"
                            );
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
    fn skew_removes_large_empty_column_runs_before_expansion() {
        let terms = schur_skew_expansion(&[100_000], &[99_999], -1).unwrap();
        assert_eq!(term_map(terms), vec![(vec![1], 1)]);
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
