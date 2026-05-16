//! Schur product and skew Schur expansion.
//!
//! Product expansion is currently correctness-first: it enumerates candidate
//! output partitions and reuses the scalar Littlewood-Richardson engine.  Skew
//! expansion uses the variable-content beta tableau enumerator, which
//! accumulates all output contents in one search.

use crate::lrcoef::{beta_lr_content_expansion, lrcoef, LrCoefError};
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
    let max_rows = row_bound(rows, size)?;
    let max_part = col_bound(cols, size)?;
    let mut terms = Vec::new();

    visit_partitions(size, max_rows, max_part, &mut Vec::new(), &mut |lambda| {
        if !contains_partition(lambda, &sh1) || !contains_partition(lambda, &sh2) {
            return Ok(());
        }
        let coefficient = lrcoef(lambda, &sh1, &sh2)?;
        if coefficient != 0 {
            terms.push(SchurTerm {
                partition: lambda.to_vec(),
                coefficient,
            });
        }
        Ok(())
    })?;

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
    let size = checked_partition_size(&outer)?
        .checked_sub(checked_partition_size(&inner)?)
        .ok_or(SchurExpansionError::ArithmeticOverflow)?;
    let max_labels = if rows >= 0 {
        Some(row_bound(rows, size)?)
    } else {
        None
    };
    let mut terms = beta_lr_content_expansion(&outer, &inner, &[], max_labels)?
        .into_iter()
        .map(|term| SchurTerm {
            partition: term.content,
            coefficient: term.coefficient,
        })
        .collect::<Vec<_>>();
    terms.sort_by(|left, right| right.partition.cmp(&left.partition));

    Ok(terms)
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

fn contains_partition(outer: &[i32], inner: &[i32]) -> bool {
    inner
        .iter()
        .enumerate()
        .all(|(index, &part)| part <= outer.get(index).copied().unwrap_or(0))
}

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
}
