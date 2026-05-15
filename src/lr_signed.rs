//! Littlewood-Richardson coefficients by a signed Kostka expansion.

use crate::kostka_fast::{FastSkewKostkaEngine, KostkaFastError};
use std::collections::HashMap;

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum SignedLrError {
    InvalidInput,
    ArithmeticOverflow,
    StateTooWide,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum SignedLrMode {
    Complete,
    ElementaryConjugate,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct SignedLrStats {
    pub value: u128,
    pub mode: SignedLrMode,
    pub determinant_len: usize,
    pub valid_permutations: u64,
    pub aggregated_terms: usize,
    pub cached_transitions: usize,
}

/// Compute `c^outer_{inner, content}` using a signed sum of skew Kostka numbers.
pub fn lrcoef_signed_kostka(
    outer: &[i32],
    inner: &[i32],
    content: &[i32],
) -> Result<u128, SignedLrError> {
    Ok(lrcoef_signed_kostka_stats(outer, inner, content)?.value)
}

/// Compute `c^outer_{inner, content}` and return expansion statistics.
pub fn lrcoef_signed_kostka_stats(
    outer: &[i32],
    inner: &[i32],
    content: &[i32],
) -> Result<SignedLrStats, SignedLrError> {
    let outer = normalize_partition(outer)?;
    let inner = normalize_partition(inner)?;
    let content = normalize_partition(content)?;

    let outer_size = checked_sum(&outer)?;
    let inner_size = checked_sum(&inner)?;
    let content_size = checked_sum(&content)?;
    if inner_size > outer_size || !partition_less_equal(&inner, &outer) {
        return Ok(zero_stats(SignedLrMode::Complete, content.len()));
    }
    if outer_size - inner_size != content_size {
        return Ok(zero_stats(SignedLrMode::Complete, content.len()));
    }
    if content_size == 0 {
        return Ok(SignedLrStats {
            value: u128::from(trim_eq(&outer, &inner)),
            mode: SignedLrMode::Complete,
            determinant_len: 0,
            valid_permutations: 1,
            aggregated_terms: 1,
            cached_transitions: 0,
        });
    }

    let conjugate_content = conjugate_partition(&content);
    let (mode, det_partition, skew_outer, skew_inner) = if conjugate_content.len() < content.len() {
        (
            SignedLrMode::ElementaryConjugate,
            conjugate_content,
            conjugate_partition(&outer),
            conjugate_partition(&inner),
        )
    } else {
        (SignedLrMode::Complete, content, outer, inner)
    };

    let mut terms = HashMap::<Vec<i32>, i64>::new();
    let mut alpha = vec![0i32; det_partition.len()];
    let mut used = vec![false; det_partition.len()];
    let mut valid_permutations = 0u64;
    generate_signed_terms(
        &det_partition,
        0,
        false,
        &mut used,
        &mut alpha,
        &mut terms,
        &mut valid_permutations,
    )?;

    let mut total = 0i128;
    let skew_outer = to_i32_vec(&skew_outer)?;
    let skew_inner = to_i32_vec(&skew_inner)?;
    let mut engine = FastSkewKostkaEngine::new(&skew_outer, &skew_inner).map_err(map_fast_error)?;
    for (weight, signed_multiplicity) in &terms {
        if *signed_multiplicity == 0 {
            continue;
        }
        let k = engine.coefficient(weight).map_err(map_fast_error)?;
        let term = i128::try_from(k)
            .map_err(|_| SignedLrError::ArithmeticOverflow)?
            .checked_mul(*signed_multiplicity as i128)
            .ok_or(SignedLrError::ArithmeticOverflow)?;
        total = total
            .checked_add(term)
            .ok_or(SignedLrError::ArithmeticOverflow)?;
    }

    if total < 0 {
        return Err(SignedLrError::ArithmeticOverflow);
    }

    Ok(SignedLrStats {
        value: total as u128,
        mode,
        determinant_len: det_partition.len(),
        valid_permutations,
        aggregated_terms: terms.values().filter(|&&value| value != 0).count(),
        cached_transitions: engine.cached_transitions(),
    })
}

fn generate_signed_terms(
    det_partition: &[u32],
    row: usize,
    odd: bool,
    used: &mut [bool],
    alpha: &mut [i32],
    terms: &mut HashMap<Vec<i32>, i64>,
    valid_permutations: &mut u64,
) -> Result<(), SignedLrError> {
    if row == det_partition.len() {
        let mut weight = alpha.to_vec();
        weight.sort_unstable_by(|a, b| b.cmp(a));
        while weight.last().is_some_and(|&part| part == 0) {
            weight.pop();
        }
        let sign = if odd { -1 } else { 1 };
        let entry = terms.entry(weight).or_insert(0);
        *entry = entry
            .checked_add(sign)
            .ok_or(SignedLrError::ArithmeticOverflow)?;
        *valid_permutations = valid_permutations
            .checked_add(1)
            .ok_or(SignedLrError::ArithmeticOverflow)?;
        return Ok(());
    }

    for col in 0..det_partition.len() {
        if used[col] {
            continue;
        }
        let value = det_partition[row] as i32 + col as i32 - row as i32;
        if value < 0 {
            continue;
        }
        let inversions_added = used[col + 1..].iter().filter(|&&is_used| is_used).count();
        used[col] = true;
        alpha[row] = value;
        generate_signed_terms(
            det_partition,
            row + 1,
            odd ^ (inversions_added % 2 == 1),
            used,
            alpha,
            terms,
            valid_permutations,
        )?;
        used[col] = false;
    }
    Ok(())
}

fn normalize_partition(parts: &[i32]) -> Result<Vec<u32>, SignedLrError> {
    if !parts.iter().all(|&part| part >= 0) || !parts.windows(2).all(|w| w[0] >= w[1]) {
        return Err(SignedLrError::InvalidInput);
    }
    Ok(parts
        .iter()
        .copied()
        .take_while(|&part| part != 0)
        .map(|part| part as u32)
        .collect())
}

fn checked_sum(parts: &[u32]) -> Result<u32, SignedLrError> {
    let mut sum = 0u32;
    for &part in parts {
        sum = sum
            .checked_add(part)
            .ok_or(SignedLrError::ArithmeticOverflow)?;
    }
    Ok(sum)
}

fn partition_less_equal(inner: &[u32], outer: &[u32]) -> bool {
    let len = inner.len().max(outer.len());
    for i in 0..len {
        if *inner.get(i).unwrap_or(&0) > *outer.get(i).unwrap_or(&0) {
            return false;
        }
    }
    true
}

fn conjugate_partition(partition: &[u32]) -> Vec<u32> {
    let Some(&width) = partition.first() else {
        return Vec::new();
    };
    let mut conjugate = vec![0; width as usize];
    for &part in partition {
        for entry in conjugate.iter_mut().take(part as usize) {
            *entry += 1;
        }
    }
    conjugate
}

fn trim_eq(left: &[u32], right: &[u32]) -> bool {
    let left_len = left
        .iter()
        .rposition(|&part| part != 0)
        .map_or(0, |index| index + 1);
    let right_len = right
        .iter()
        .rposition(|&part| part != 0)
        .map_or(0, |index| index + 1);
    left[..left_len] == right[..right_len]
}

fn to_i32_vec(parts: &[u32]) -> Result<Vec<i32>, SignedLrError> {
    parts
        .iter()
        .copied()
        .map(|part| i32::try_from(part).map_err(|_| SignedLrError::ArithmeticOverflow))
        .collect()
}

fn map_fast_error(error: KostkaFastError) -> SignedLrError {
    match error {
        KostkaFastError::InvalidInput => SignedLrError::InvalidInput,
        KostkaFastError::ArithmeticOverflow => SignedLrError::ArithmeticOverflow,
        KostkaFastError::StateTooWide => SignedLrError::StateTooWide,
    }
}

fn zero_stats(mode: SignedLrMode, determinant_len: usize) -> SignedLrStats {
    SignedLrStats {
        value: 0,
        mode,
        determinant_len,
        valid_permutations: 0,
        aggregated_terms: 0,
        cached_transitions: 0,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::lrcoef::lrcoef;

    #[test]
    fn computes_basic_lr_coefficients() {
        assert_eq!(lrcoef_signed_kostka(&[3, 2, 1], &[2, 1], &[2, 1]), Ok(2));
        assert_eq!(lrcoef_signed_kostka(&[4, 2], &[2, 1], &[2, 1]), Ok(1));
        assert_eq!(lrcoef_signed_kostka(&[5, 1], &[2, 1], &[2, 1]), Ok(0));
    }

    #[test]
    fn uses_conjugate_for_tall_content() {
        let stats = lrcoef_signed_kostka_stats(
            &[7, 4, 3, 1],
            &[0],
            &[1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
        )
        .unwrap();
        assert_eq!(stats.mode, SignedLrMode::ElementaryConjugate);
        assert_eq!(stats.determinant_len, 1);
    }

    #[test]
    fn agrees_with_lr_port_for_representative_cases() {
        let cases = [
            (&[2, 1][..], &[2][..], &[1][..]),
            (&[3, 2, 1][..], &[2, 1][..], &[2, 1][..]),
            (&[5, 4, 3, 2, 1][..], &[3, 2, 1][..], &[4, 3, 1, 1][..]),
            (
                &[7, 6, 5, 4, 3, 2, 1][..],
                &[4, 4, 3, 2, 1][..],
                &[5, 4, 3, 2][..],
            ),
            (&[9, 7, 5, 3, 1][..], &[5, 3, 2, 1][..], &[6, 5, 3, 2][..]),
        ];
        for (outer, inner, content) in cases {
            assert_eq!(
                lrcoef_signed_kostka(outer, inner, content).unwrap(),
                lrcoef(outer, inner, content).unwrap()
            );
        }
    }
}
