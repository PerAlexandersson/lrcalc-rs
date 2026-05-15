//! Kostka coefficients reduced to single Littlewood-Richardson coefficients.

use crate::lrcoef::{lrcoef, LrCoefError};

pub type LrTriple = (Vec<i32>, Vec<i32>, Vec<i32>);

/// Compute the ordinary Kostka coefficient `K_{shape, weight}` via LR.
pub fn kostka_via_lr(shape: &[i32], weight: &[i32]) -> Result<u128, LrCoefError> {
    let shape_sum = checked_sum_nonnegative(shape)?;
    let weight_sum = checked_sum_nonnegative(weight)?;
    if shape_sum != weight_sum {
        return Ok(0);
    }
    let (outer, inner, content) = kostka_lr_triple(shape, weight)?;
    lrcoef(&outer, &inner, &content)
}

/// Return `(tau, sigma, shape)` such that `K_{shape, weight} = c^tau_{sigma, shape}`.
pub fn kostka_lr_triple(shape: &[i32], weight: &[i32]) -> Result<LrTriple, LrCoefError> {
    if !valid_partition(shape) || !weight.iter().all(|&part| part >= 0) {
        return Err(LrCoefError::InvalidPartition);
    }
    let row_components: Vec<Vec<i32>> = weight
        .iter()
        .copied()
        .filter(|&part| part > 0)
        .map(|part| vec![part])
        .collect();
    let (outer, inner) = diagonal_concatenation(&row_components)?;
    Ok((outer, inner, trim_partition(shape)))
}

fn diagonal_concatenation(components: &[Vec<i32>]) -> Result<(Vec<i32>, Vec<i32>), LrCoefError> {
    let mut trimmed = Vec::new();
    for component in components {
        if !valid_partition(component) {
            return Err(LrCoefError::InvalidPartition);
        }
        let length = partition_length(component);
        if length != 0 {
            trimmed.push(&component[..length]);
        }
    }

    let mut suffix_widths = vec![0i32; trimmed.len() + 1];
    for i in (0..trimmed.len()).rev() {
        suffix_widths[i] = suffix_widths[i + 1]
            .checked_add(trimmed[i][0])
            .ok_or(LrCoefError::ArithmeticOverflow)?;
    }

    let row_count: usize = trimmed.iter().map(|component| component.len()).sum();
    let mut outer = Vec::with_capacity(row_count);
    let mut inner = Vec::with_capacity(row_count);
    for (i, component) in trimmed.iter().enumerate() {
        let offset = suffix_widths[i + 1];
        for &part in *component {
            outer.push(
                offset
                    .checked_add(part)
                    .ok_or(LrCoefError::ArithmeticOverflow)?,
            );
            inner.push(offset);
        }
    }

    trim_trailing_zeroes(&mut inner);
    Ok((outer, inner))
}

fn valid_partition(partition: &[i32]) -> bool {
    partition.iter().all(|&part| part >= 0)
        && partition.windows(2).all(|window| window[0] >= window[1])
}

fn partition_length(partition: &[i32]) -> usize {
    partition
        .iter()
        .rposition(|&part| part != 0)
        .map_or(0, |index| index + 1)
}

fn trim_partition(partition: &[i32]) -> Vec<i32> {
    partition[..partition_length(partition)].to_vec()
}

fn trim_trailing_zeroes(partition: &mut Vec<i32>) {
    while partition.last().is_some_and(|&part| part == 0) {
        partition.pop();
    }
}

fn checked_sum_nonnegative(values: &[i32]) -> Result<i32, LrCoefError> {
    let mut sum = 0i32;
    for &value in values {
        if value < 0 {
            return Err(LrCoefError::InvalidPartition);
        }
        sum = sum
            .checked_add(value)
            .ok_or(LrCoefError::ArithmeticOverflow)?;
    }
    Ok(sum)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn translates_weight_rows_to_skew_shape() {
        let (outer, inner, content) = kostka_lr_triple(&[2, 1], &[2, 1]).unwrap();
        assert_eq!(outer, vec![3, 1]);
        assert_eq!(inner, vec![1]);
        assert_eq!(content, vec![2, 1]);
    }

    #[test]
    fn computes_small_kostka_coefficients_via_lr() {
        assert_eq!(kostka_via_lr(&[], &[]), Ok(1));
        assert_eq!(kostka_via_lr(&[3], &[2, 1]), Ok(1));
        assert_eq!(kostka_via_lr(&[2, 1], &[2, 1]), Ok(1));
        assert_eq!(kostka_via_lr(&[1, 1, 1], &[2, 1]), Ok(0));
        assert_eq!(kostka_via_lr(&[2], &[1, 1]), Ok(1));
        assert_eq!(kostka_via_lr(&[1, 1], &[1, 1]), Ok(1));
        assert_eq!(kostka_via_lr(&[2, 1], &[1, 1, 1]), Ok(2));
    }

    #[test]
    fn weight_order_does_not_affect_result() {
        assert_eq!(kostka_via_lr(&[3, 1], &[2, 1, 1]), Ok(2));
        assert_eq!(kostka_via_lr(&[3, 1], &[1, 2, 1]), Ok(2));
    }

    #[test]
    fn size_mismatch_is_zero() {
        assert_eq!(kostka_via_lr(&[3], &[1, 1]), Ok(0));
    }
}
