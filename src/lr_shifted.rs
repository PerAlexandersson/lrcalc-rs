//! Shifted-Schur interval dynamic program for scalar LR coefficients.
//!
//! This is an experimental comparison engine.  It computes ordinary
//! Littlewood-Richardson coefficients through the recurrence for shifted
//! Schur structure constants, after the usual Buch scalar normalization.

use std::collections::HashMap;

use num_bigint::{BigInt, Sign};
use num_traits::{One, ToPrimitive, Zero};

use crate::lrcoef::{optim_coef, LrCoefError, OptimizedCoef};

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct ShiftedIntervalStats {
    pub value: u128,
    pub upper: Vec<i32>,
    pub lower: Vec<i32>,
    pub fixed: Vec<i32>,
    pub interval_partitions: usize,
    pub base_states: usize,
    pub total_pair_states: usize,
    pub peak_layer_states: usize,
    pub layer_states: Vec<usize>,
}

#[derive(Clone, Debug)]
struct IntervalData {
    parts: Vec<Vec<i32>>,
    by_rank: Vec<Vec<usize>>,
    add_neighbors: Vec<Vec<usize>>,
    remove_neighbors: Vec<Vec<usize>>,
    lower_index: usize,
    upper_index: usize,
}

pub fn lrcoef_shifted_interval_u128(
    outer: &[i32],
    inner1: &[i32],
    inner2: &[i32],
) -> Result<u128, LrCoefError> {
    Ok(lrcoef_shifted_interval_stats(outer, inner1, inner2)?.value)
}

pub fn lrcoef_shifted_interval_stats(
    outer: &[i32],
    inner1: &[i32],
    inner2: &[i32],
) -> Result<ShiftedIntervalStats, LrCoefError> {
    match optim_coef(outer, inner1, inner2)? {
        OptimizedCoef::Zero => Ok(ShiftedIntervalStats {
            value: 0,
            upper: Vec::new(),
            lower: Vec::new(),
            fixed: Vec::new(),
            interval_partitions: 0,
            base_states: 0,
            total_pair_states: 0,
            peak_layer_states: 0,
            layer_states: Vec::new(),
        }),
        OptimizedCoef::One => Ok(ShiftedIntervalStats {
            value: 1,
            upper: Vec::new(),
            lower: Vec::new(),
            fixed: Vec::new(),
            interval_partitions: 1,
            base_states: 1,
            total_pair_states: 1,
            peak_layer_states: 1,
            layer_states: vec![1],
        }),
        OptimizedCoef::Count(shape) => {
            let inner_size = partition_sum_usize(&shape.inner)?;
            let content_size = partition_sum_usize(&shape.content)?;
            let (fixed, lower) = if inner_size <= content_size {
                (shape.inner, shape.content)
            } else {
                (shape.content, shape.inner)
            };
            shifted_interval_stats_for_normalized(shape.outer, lower, fixed)
        }
    }
}

fn shifted_interval_stats_for_normalized(
    upper: Vec<i32>,
    lower: Vec<i32>,
    fixed: Vec<i32>,
) -> Result<ShiftedIntervalStats, LrCoefError> {
    if !partition_leq(&lower, &upper) {
        return Ok(ShiftedIntervalStats {
            value: 0,
            upper,
            lower,
            fixed,
            interval_partitions: 0,
            base_states: 0,
            total_pair_states: 0,
            peak_layer_states: 0,
            layer_states: Vec::new(),
        });
    }

    let fixed_size = partition_sum_usize(&fixed)?;
    let interval = build_interval(&lower, &upper)?;
    let base_states = interval.parts.len();
    let mut previous = HashMap::<u64, BigInt>::with_capacity(base_states);
    for (index, part) in interval.parts.iter().enumerate() {
        let value = shifted_schur_eval(&fixed, part)?;
        previous.insert(pair_key(index, index), value);
    }

    let mut layer_states = vec![base_states];
    let mut total_pair_states = base_states;
    let mut peak_layer_states = base_states;

    for gap in 1..=fixed_size {
        let mut current = HashMap::new();
        for upper_rank in gap..interval.by_rank.len() {
            let lower_rank = upper_rank - gap;
            for &upper_index in &interval.by_rank[upper_rank] {
                for &lower_index in &interval.by_rank[lower_rank] {
                    if !partition_leq(&interval.parts[lower_index], &interval.parts[upper_index]) {
                        continue;
                    }
                    let value = shifted_interval_next_value(
                        &interval,
                        &previous,
                        upper_index,
                        lower_index,
                        gap,
                    )?;
                    current.insert(pair_key(upper_index, lower_index), value);
                }
            }
        }
        let count = current.len();
        total_pair_states = total_pair_states
            .checked_add(count)
            .ok_or(LrCoefError::ArithmeticOverflow)?;
        peak_layer_states = peak_layer_states.max(count);
        layer_states.push(count);
        previous = current;
    }

    let value = previous
        .get(&pair_key(interval.upper_index, interval.lower_index))
        .cloned()
        .unwrap_or_else(BigInt::zero);
    let value = big_int_to_u128(value)?;

    Ok(ShiftedIntervalStats {
        value,
        upper,
        lower,
        fixed,
        interval_partitions: interval.parts.len(),
        base_states,
        total_pair_states,
        peak_layer_states,
        layer_states,
    })
}

fn shifted_interval_next_value(
    interval: &IntervalData,
    previous: &HashMap<u64, BigInt>,
    upper_index: usize,
    lower_index: usize,
    gap: usize,
) -> Result<BigInt, LrCoefError> {
    let mut numerator = BigInt::zero();
    for &lower_plus in &interval.add_neighbors[lower_index] {
        if partition_leq(&interval.parts[lower_plus], &interval.parts[upper_index]) {
            if let Some(value) = previous.get(&pair_key(upper_index, lower_plus)) {
                numerator += value;
            }
        }
    }
    for &upper_minus in &interval.remove_neighbors[upper_index] {
        if partition_leq(&interval.parts[lower_index], &interval.parts[upper_minus]) {
            if let Some(value) = previous.get(&pair_key(upper_minus, lower_index)) {
                numerator -= value;
            }
        }
    }

    let divisor = BigInt::from(gap);
    if &numerator % &divisor != BigInt::zero() {
        return Err(LrCoefError::ArithmeticOverflow);
    }
    Ok(numerator / divisor)
}

fn build_interval(lower: &[i32], upper: &[i32]) -> Result<IntervalData, LrCoefError> {
    let len = lower.len().max(upper.len());
    let lower = pad_partition(lower, len);
    let upper = pad_partition(upper, len);
    let lower_size = partition_sum_usize(&lower)?;
    let upper_size = partition_sum_usize(&upper)?;
    let max_rank = upper_size
        .checked_sub(lower_size)
        .ok_or(LrCoefError::ArithmeticOverflow)?;

    let mut parts = Vec::new();
    let mut current = vec![0; len];
    generate_interval_partitions(0, &lower, &upper, &mut current, &mut parts);
    let mut index_by_part = HashMap::with_capacity(parts.len());
    for (index, part) in parts.iter().enumerate() {
        index_by_part.insert(part.clone(), index);
    }

    let mut by_rank = vec![Vec::new(); max_rank + 1];
    for (index, part) in parts.iter().enumerate() {
        let rank = partition_sum_usize(part)?
            .checked_sub(lower_size)
            .ok_or(LrCoefError::ArithmeticOverflow)?;
        by_rank[rank].push(index);
    }

    let lower_index = *index_by_part
        .get(&lower)
        .ok_or(LrCoefError::ArithmeticOverflow)?;
    let upper_index = *index_by_part
        .get(&upper)
        .ok_or(LrCoefError::ArithmeticOverflow)?;

    let mut add_neighbors = vec![Vec::new(); parts.len()];
    let mut remove_neighbors = vec![Vec::new(); parts.len()];
    for (index, part) in parts.iter().enumerate() {
        for row in 0..len {
            let mut next = part.clone();
            next[row] = next[row]
                .checked_add(1)
                .ok_or(LrCoefError::ArithmeticOverflow)?;
            if next[row] <= upper[row] && is_partition(&next) {
                if let Some(&next_index) = index_by_part.get(&next) {
                    add_neighbors[index].push(next_index);
                }
            }

            if part[row] > lower[row] {
                let mut previous = part.clone();
                previous[row] -= 1;
                if is_partition(&previous) {
                    if let Some(&previous_index) = index_by_part.get(&previous) {
                        remove_neighbors[index].push(previous_index);
                    }
                }
            }
        }
    }

    Ok(IntervalData {
        parts,
        by_rank,
        add_neighbors,
        remove_neighbors,
        lower_index,
        upper_index,
    })
}

fn generate_interval_partitions(
    row: usize,
    lower: &[i32],
    upper: &[i32],
    current: &mut [i32],
    out: &mut Vec<Vec<i32>>,
) {
    if row == current.len() {
        out.push(current.to_vec());
        return;
    }

    let previous_bound = if row == 0 {
        upper[row]
    } else {
        current[row - 1].min(upper[row])
    };
    let next_lower = lower.get(row + 1).copied().unwrap_or(0);
    let minimum = lower[row].max(next_lower);
    for value in (minimum..=previous_bound).rev() {
        current[row] = value;
        generate_interval_partitions(row + 1, lower, upper, current, out);
    }
}

fn shifted_schur_eval(mu: &[i32], lambda: &[i32]) -> Result<BigInt, LrCoefError> {
    let n = mu.len().max(lambda.len());
    if n == 0 {
        return Ok(BigInt::one());
    }
    let mu = pad_partition(mu, n);
    let lambda = pad_partition(lambda, n);
    let mut numerator = vec![vec![BigInt::zero(); n]; n];
    let mut denominator = vec![vec![BigInt::zero(); n]; n];

    for i in 0..n {
        let x = i64::from(lambda[i])
            .checked_add(i64::try_from(n - i - 1).map_err(|_| LrCoefError::ArithmeticOverflow)?)
            .ok_or(LrCoefError::ArithmeticOverflow)?;
        for j in 0..n {
            let num_degree = usize::try_from(
                mu[j]
                    .checked_add(
                        i32::try_from(n - j - 1).map_err(|_| LrCoefError::ArithmeticOverflow)?,
                    )
                    .ok_or(LrCoefError::ArithmeticOverflow)?,
            )
            .map_err(|_| LrCoefError::ArithmeticOverflow)?;
            numerator[i][j] = falling_factorial(x, num_degree)?;
            denominator[i][j] = falling_factorial(x, n - j - 1)?;
        }
    }

    let denominator = determinant_bareiss(denominator)?;
    if denominator.is_zero() {
        return Err(LrCoefError::ArithmeticOverflow);
    }
    let numerator = determinant_bareiss(numerator)?;
    if &numerator % &denominator != BigInt::zero() {
        return Err(LrCoefError::ArithmeticOverflow);
    }
    Ok(numerator / denominator)
}

fn determinant_bareiss(mut matrix: Vec<Vec<BigInt>>) -> Result<BigInt, LrCoefError> {
    let n = matrix.len();
    if n == 0 {
        return Ok(BigInt::one());
    }
    if n == 1 {
        return Ok(matrix[0][0].clone());
    }

    let mut sign = BigInt::one();
    let mut previous_pivot = BigInt::one();
    for k in 0..n - 1 {
        let Some(pivot_row) = (k..n).find(|&row| !matrix[row][k].is_zero()) else {
            return Ok(BigInt::zero());
        };
        if pivot_row != k {
            matrix.swap(k, pivot_row);
            sign = -sign;
        }
        let pivot = matrix[k][k].clone();
        for i in k + 1..n {
            for j in k + 1..n {
                let value =
                    (&matrix[i][j] * &pivot - &matrix[i][k] * &matrix[k][j]) / &previous_pivot;
                matrix[i][j] = value;
            }
        }
        for row in matrix.iter_mut().take(n).skip(k + 1) {
            row[k] = BigInt::zero();
        }
        previous_pivot = pivot;
    }
    Ok(sign * matrix[n - 1][n - 1].clone())
}

fn falling_factorial(x: i64, degree: usize) -> Result<BigInt, LrCoefError> {
    let mut value = BigInt::one();
    for offset in 0..degree {
        let offset = i64::try_from(offset).map_err(|_| LrCoefError::ArithmeticOverflow)?;
        value *= x
            .checked_sub(offset)
            .ok_or(LrCoefError::ArithmeticOverflow)?;
    }
    Ok(value)
}

fn pair_key(upper_index: usize, lower_index: usize) -> u64 {
    ((upper_index as u64) << 32) | lower_index as u64
}

fn big_int_to_u128(value: BigInt) -> Result<u128, LrCoefError> {
    if value.sign() == Sign::Minus {
        return Err(LrCoefError::ArithmeticOverflow);
    }
    value.to_u128().ok_or(LrCoefError::ArithmeticOverflow)
}

fn pad_partition(partition: &[i32], len: usize) -> Vec<i32> {
    let mut result = vec![0; len];
    for (index, &part) in partition.iter().take(len).enumerate() {
        result[index] = part;
    }
    result
}

fn partition_sum_usize(partition: &[i32]) -> Result<usize, LrCoefError> {
    let mut sum = 0usize;
    for &part in partition {
        if part < 0 {
            return Err(LrCoefError::InvalidPartition);
        }
        sum = sum
            .checked_add(usize::try_from(part).map_err(|_| LrCoefError::InvalidPartition)?)
            .ok_or(LrCoefError::ArithmeticOverflow)?;
    }
    Ok(sum)
}

fn partition_leq(left: &[i32], right: &[i32]) -> bool {
    let len = left.len().max(right.len());
    (0..len).all(|index| part_entry(left, index) <= part_entry(right, index))
}

fn part_entry(partition: &[i32], index: usize) -> i32 {
    partition.get(index).copied().unwrap_or(0)
}

fn is_partition(partition: &[i32]) -> bool {
    partition.iter().all(|&part| part >= 0)
        && partition.windows(2).all(|window| window[0] >= window[1])
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::lrcoef::lrcoef;

    #[test]
    fn shifted_schur_base_values_for_one_box() {
        assert_eq!(shifted_schur_eval(&[1], &[1]).unwrap(), BigInt::from(1));
        assert_eq!(shifted_schur_eval(&[1], &[2]).unwrap(), BigInt::from(2));
        assert_eq!(shifted_schur_eval(&[1], &[1, 1]).unwrap(), BigInt::from(2));
    }

    #[test]
    fn computes_basic_lr_coefficients() {
        let cases = [
            (&[2][..], &[1][..], &[1][..]),
            (&[1, 1][..], &[1][..], &[1][..]),
            (&[3, 2, 1][..], &[2, 1][..], &[2, 1][..]),
            (&[4, 2][..], &[2, 1][..], &[2, 1][..]),
        ];
        for (outer, inner, content) in cases {
            assert_eq!(
                lrcoef_shifted_interval_u128(outer, inner, content).unwrap(),
                lrcoef(outer, inner, content).unwrap(),
                "outer={outer:?} inner={inner:?} content={content:?}"
            );
        }
    }

    #[test]
    fn agrees_with_buch_for_small_triples() {
        let parts = small_partitions(6, 4);
        for outer in &parts {
            let outer_size = partition_sum_usize(outer).unwrap();
            for inner in &parts {
                let inner_size = partition_sum_usize(inner).unwrap();
                if inner_size > outer_size || !partition_leq(inner, outer) {
                    continue;
                }
                for content in &parts {
                    if inner_size + partition_sum_usize(content).unwrap() != outer_size {
                        continue;
                    }
                    let shifted = lrcoef_shifted_interval_u128(outer, inner, content).unwrap();
                    let buch = lrcoef(outer, inner, content).unwrap();
                    assert_eq!(
                        shifted, buch,
                        "outer={outer:?} inner={inner:?} content={content:?}"
                    );
                }
            }
        }
    }

    fn small_partitions(max_size: i32, max_rows: usize) -> Vec<Vec<i32>> {
        let mut result = Vec::new();
        for size in 0..=max_size {
            let mut current = Vec::new();
            visit_partitions(size, max_rows, size, &mut current, &mut result);
        }
        result
    }

    fn visit_partitions(
        remaining: i32,
        rows_left: usize,
        max_part: i32,
        current: &mut Vec<i32>,
        out: &mut Vec<Vec<i32>>,
    ) {
        if remaining == 0 {
            out.push(current.clone());
            return;
        }
        if rows_left == 0 || max_part == 0 {
            return;
        }
        for part in (1..=remaining.min(max_part)).rev() {
            current.push(part);
            visit_partitions(remaining - part, rows_left - 1, part, current, out);
            current.pop();
        }
    }
}
