//! Ehrhart interpolation for stretched Littlewood-Richardson coefficients.

use crate::lrcoef::{
    beta_lrcoef_buch_stretch_cache, beta_lrcoef_buch_stretched_counts_u128,
    lrcoef_buch_stretch_cache, lrcoef_buch_stretched_counts_u128, BetaLrBuchStretchCache,
    LrBuchCounts, LrBuchStretchCache, LrCoefError,
};
use num_bigint::{BigInt, ToBigInt};
use num_rational::BigRational;
use num_traits::{One, Zero};

#[derive(Clone, Debug, Eq, PartialEq)]
pub enum LrStretchError {
    InvalidInput,
    ArithmeticOverflow,
    StateTooWide,
    SingularInterpolation,
    NonIntegralValue,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct LrStretchPolynomial {
    pub dimension: usize,
    pub coefficients: Vec<BigRational>,
    pub h_vector: Vec<BigInt>,
    pub sample_points: Vec<(i64, BigInt)>,
}

/// Compute the Ehrhart h*-vector for
/// `t ↦ c^(t outer)_(t inner, t content)`.
pub fn lr_stretch_h_vector(
    outer: &[i32],
    inner: &[i32],
    content: &[i32],
) -> Result<LrStretchPolynomial, LrStretchError> {
    compute_lr_stretch_polynomial(outer, inner, content)
}

/// Evaluate the stretched coefficient at a large stretch factor using the
/// interpolated h*-vector.
pub fn lr_stretch_coefficient(
    outer: &[i32],
    inner: &[i32],
    content: &[i32],
    stretch: u64,
) -> Result<BigInt, LrStretchError> {
    let polynomial = compute_lr_stretch_polynomial(outer, inner, content)?;
    Ok(evaluate_h_vector(
        &polynomial.h_vector,
        polynomial.dimension,
        stretch,
    ))
}

/// Compute the Ehrhart h*-vector for the beta-shifted stretch
/// `t ↦ beta_lrcoef(t outer, t inner, t content, t beta)`.
pub fn beta_lr_stretch_h_vector(
    outer: &[i32],
    inner: &[i32],
    content: &[i32],
    beta: &[i32],
) -> Result<LrStretchPolynomial, LrStretchError> {
    compute_beta_lr_stretch_polynomial(outer, inner, content, beta)
}

/// Evaluate a beta-shifted stretched coefficient using interpolation.
pub fn beta_lr_stretch_coefficient(
    outer: &[i32],
    inner: &[i32],
    content: &[i32],
    beta: &[i32],
    stretch: u64,
) -> Result<BigInt, LrStretchError> {
    let polynomial = compute_beta_lr_stretch_polynomial(outer, inner, content, beta)?;
    Ok(evaluate_h_vector(
        &polynomial.h_vector,
        polynomial.dimension,
        stretch,
    ))
}

pub fn compute_lr_stretch_polynomial(
    outer: &[i32],
    inner: &[i32],
    content: &[i32],
) -> Result<LrStretchPolynomial, LrStretchError> {
    let Some(stretch_cache) =
        lrcoef_buch_stretch_cache(outer, inner, content).map_err(map_lrcoef_error)?
    else {
        return Ok(LrStretchPolynomial {
            dimension: 0,
            coefficients: vec![BigRational::zero()],
            h_vector: vec![BigInt::zero()],
            sample_points: Vec::new(),
        });
    };
    let dimension = stretch_cache.dimension();

    interpolate_stretch_polynomial(dimension, |stretch| {
        scaled_lr_counts(&stretch_cache, stretch)
    })
}

pub fn compute_beta_lr_stretch_polynomial(
    outer: &[i32],
    inner: &[i32],
    content: &[i32],
    beta: &[i32],
) -> Result<LrStretchPolynomial, LrStretchError> {
    let Some(stretch_cache) =
        beta_lrcoef_buch_stretch_cache(outer, inner, content, beta).map_err(map_lrcoef_error)?
    else {
        return Ok(LrStretchPolynomial {
            dimension: 0,
            coefficients: vec![BigRational::zero()],
            h_vector: vec![BigInt::zero()],
            sample_points: Vec::new(),
        });
    };
    let dimension = stretch_cache.dimension();

    interpolate_stretch_polynomial(dimension, |stretch| {
        scaled_beta_lr_counts(&stretch_cache, stretch)
    })
}

fn interpolate_stretch_polynomial<F>(
    dimension: usize,
    mut counts_at: F,
) -> Result<LrStretchPolynomial, LrStretchError>
where
    F: FnMut(u64) -> Result<LrBuchCounts, LrStretchError>,
{
    if dimension == 0 {
        return Ok(LrStretchPolynomial {
            dimension,
            coefficients: vec![BigRational::one()],
            h_vector: vec![BigInt::one()],
            sample_points: vec![(0, BigInt::one())],
        });
    }

    let mut points = Vec::<(i64, BigRational)>::with_capacity(dimension + 1);
    let mut sample_points = Vec::<(i64, BigInt)>::with_capacity(dimension + 1);
    points.push((0, BigRational::one()));
    sample_points.push((0, BigInt::one()));

    let sign_positive = dimension % 2 == 0;
    let mut stretch = 1u64;
    while points.len() <= dimension {
        let counts = counts_at(stretch)?;
        let full = counts
            .full
            .to_bigint()
            .ok_or(LrStretchError::ArithmeticOverflow)?;
        points.push((stretch as i64, BigRational::from_integer(full.clone())));
        sample_points.push((stretch as i64, full));

        if points.len() <= dimension {
            let count = counts
                .interior
                .to_bigint()
                .ok_or(LrStretchError::ArithmeticOverflow)?;
            let signed = if sign_positive {
                count.clone()
            } else {
                -count.clone()
            };
            points.push((-(stretch as i64), BigRational::from_integer(signed)));
            sample_points.push((-(stretch as i64), count));
        }
        stretch = stretch
            .checked_add(1)
            .ok_or(LrStretchError::ArithmeticOverflow)?;
    }

    let coefficients = interpolate(&points)?;
    let h_vector = compute_h_vector_from_coefficients(&coefficients, dimension)?;
    Ok(LrStretchPolynomial {
        dimension,
        coefficients,
        h_vector,
        sample_points,
    })
}

pub fn evaluate_h_vector(h_vector: &[BigInt], dimension: usize, stretch: u64) -> BigInt {
    let mut total = BigInt::zero();
    for (k, h_k) in h_vector.iter().enumerate().take(dimension + 1) {
        if h_k.is_zero() {
            continue;
        }
        let top = stretch
            .checked_add(dimension as u64)
            .and_then(|value| value.checked_sub(k as u64))
            .expect("dimension index should not exceed stretch + dimension");
        total += h_k * binomial_bigint(top, dimension as u64);
    }
    total
}

pub fn format_bigint_vector(vector: &[BigInt]) -> String {
    format!(
        "[{}]",
        vector
            .iter()
            .map(BigInt::to_string)
            .collect::<Vec<_>>()
            .join(", ")
    )
}

fn scaled_lr_counts(
    stretch_cache: &LrBuchStretchCache,
    stretch: u64,
) -> Result<LrBuchCounts, LrStretchError> {
    lrcoef_buch_stretched_counts_u128(stretch_cache, stretch).map_err(map_lrcoef_error)
}

fn scaled_beta_lr_counts(
    stretch_cache: &BetaLrBuchStretchCache,
    stretch: u64,
) -> Result<LrBuchCounts, LrStretchError> {
    beta_lrcoef_buch_stretched_counts_u128(stretch_cache, stretch).map_err(map_lrcoef_error)
}

fn compute_h_vector_from_coefficients(
    coefficients: &[BigRational],
    dimension: usize,
) -> Result<Vec<BigInt>, LrStretchError> {
    let mut h_vector = vec![BigInt::zero(); dimension + 1];
    for (k, h_entry) in h_vector.iter_mut().enumerate().take(dimension + 1) {
        let mut value = BigInt::zero();
        for j in 0..=k {
            let p_j = evaluate_coefficients(coefficients, j as i64);
            if !p_j.is_integer() {
                return Err(LrStretchError::NonIntegralValue);
            }
            let p_j = p_j.to_integer();
            let binom = binomial_bigint((dimension + 1) as u64, (k - j) as u64);
            if (k - j) % 2 == 0 {
                value += binom * p_j;
            } else {
                value -= binom * p_j;
            }
        }
        *h_entry = value;
    }
    Ok(h_vector)
}

fn evaluate_coefficients(coefficients: &[BigRational], x: i64) -> BigRational {
    let x = BigRational::from_integer(BigInt::from(x));
    let mut result = BigRational::zero();
    let mut power = BigRational::one();
    for coefficient in coefficients {
        result += coefficient * &power;
        power *= &x;
    }
    result
}

fn interpolate(points: &[(i64, BigRational)]) -> Result<Vec<BigRational>, LrStretchError> {
    let n = points.len();
    let mut matrix: Vec<Vec<BigRational>> = points
        .iter()
        .map(|&(x, ref y)| {
            let x = BigInt::from(x);
            let mut row = Vec::with_capacity(n + 1);
            let mut power = BigInt::one();
            for _ in 0..n {
                row.push(BigRational::from_integer(power.clone()));
                power *= &x;
            }
            row.push(y.clone());
            row
        })
        .collect();

    for col in 0..n {
        let Some(pivot_row) = (col..n).find(|&row| !matrix[row][col].is_zero()) else {
            return Err(LrStretchError::SingularInterpolation);
        };
        matrix.swap(col, pivot_row);
        let pivot = matrix[col][col].clone();
        for entry in matrix[col].iter_mut().skip(col) {
            *entry /= pivot.clone();
        }
        let pivot_tail = matrix[col][col..=n].to_vec();
        for (row, row_entries) in matrix.iter_mut().enumerate().take(n) {
            if row == col || row_entries[col].is_zero() {
                continue;
            }
            let factor = row_entries[col].clone();
            for (entry, pivot_entry) in row_entries[col..=n].iter_mut().zip(&pivot_tail) {
                let sub = factor.clone() * pivot_entry;
                *entry -= sub;
            }
        }
    }

    Ok(matrix.into_iter().map(|row| row[n].clone()).collect())
}

fn binomial_bigint(n: u64, k: u64) -> BigInt {
    if k > n {
        return BigInt::zero();
    }
    let k = k.min(n - k);
    let mut result = BigInt::one();
    for i in 1..=k {
        result = result * BigInt::from(n - k + i) / BigInt::from(i);
    }
    result
}

fn map_lrcoef_error(error: LrCoefError) -> LrStretchError {
    match error {
        LrCoefError::InvalidPartition => LrStretchError::InvalidInput,
        LrCoefError::ArithmeticOverflow => LrStretchError::ArithmeticOverflow,
    }
}

pub fn format_error(error: LrStretchError) -> String {
    match error {
        LrStretchError::InvalidInput => "invalid input".to_string(),
        LrStretchError::ArithmeticOverflow => "arithmetic overflow".to_string(),
        LrStretchError::StateTooWide => "state does not fit in u128".to_string(),
        LrStretchError::SingularInterpolation => "singular interpolation system".to_string(),
        LrStretchError::NonIntegralValue => "interpolation produced non-integral value".to_string(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::lrcoef::{beta_lrcoef_buch_stretch_cache, beta_lrcoef_buch_stretched_counts_u128};

    #[test]
    fn evaluates_h_vector_formula() {
        let h = vec![BigInt::one()];
        assert_eq!(evaluate_h_vector(&h, 0, 100), BigInt::one());
    }

    #[test]
    fn stretch_polynomial_matches_direct_counts() {
        let cases = [
            (&[3, 2, 1][..], &[2, 1][..], &[2, 1][..]),
            (&[4, 2][..], &[2, 1][..], &[2, 1][..]),
        ];
        for (outer, inner, content) in cases {
            let polynomial = compute_lr_stretch_polynomial(outer, inner, content).unwrap();
            let stretch_cache = lrcoef_buch_stretch_cache(outer, inner, content)
                .unwrap()
                .expect("test cases should have nonzero stretch families");
            assert_eq!(polynomial.h_vector.len(), polynomial.dimension + 1);
            for stretch in 1..=4 {
                let interpolated =
                    evaluate_h_vector(&polynomial.h_vector, polynomial.dimension, stretch);
                let direct = scaled_lr_counts(&stretch_cache, stretch)
                    .unwrap_or_else(|_| panic!("direct count failed at stretch {stretch}"));
                let direct = direct
                    .full
                    .to_bigint()
                    .expect("u128 should convert to BigInt");
                assert_eq!(
                    interpolated, direct,
                    "outer={outer:?} inner={inner:?} content={content:?} stretch={stretch}"
                );
            }
        }
    }

    #[test]
    fn beta_stretch_polynomial_matches_direct_counts() {
        let cases = [
            (&[3, 2, 1][..], &[2, 1][..], &[2, 1][..], &[][..]),
            (&[5, 3, 1][..], &[3, 2, 1][..], &[2, 1][..], &[2, 0][..]),
            (&[5, 3, 1][..], &[3, 2, 1][..], &[1, 2][..], &[3, 0][..]),
        ];
        for (outer, inner, content, beta) in cases {
            let polynomial =
                compute_beta_lr_stretch_polynomial(outer, inner, content, beta).unwrap();
            let stretch_cache = beta_lrcoef_buch_stretch_cache(outer, inner, content, beta)
                .unwrap()
                .expect("test cases should have nonzero stretch families");
            assert_eq!(polynomial.h_vector.len(), polynomial.dimension + 1);
            for stretch in 1..=4 {
                let interpolated =
                    evaluate_h_vector(&polynomial.h_vector, polynomial.dimension, stretch);
                let direct = beta_lrcoef_buch_stretched_counts_u128(&stretch_cache, stretch)
                    .unwrap_or_else(|_| panic!("direct count failed at stretch {stretch}"));
                let direct = direct
                    .full
                    .to_bigint()
                    .expect("u128 should convert to BigInt");
                assert_eq!(
                    interpolated, direct,
                    "outer={outer:?} inner={inner:?} content={content:?} beta={beta:?} stretch={stretch}"
                );
            }
        }
    }

    #[test]
    fn beta_stretch_with_empty_beta_specializes_to_lr_stretch() {
        let cases = [
            (&[3, 2, 1][..], &[2, 1][..], &[2, 1][..]),
            (&[4, 2][..], &[2, 1][..], &[2, 1][..]),
        ];
        for (outer, inner, content) in cases {
            assert_eq!(
                compute_beta_lr_stretch_polynomial(outer, inner, content, &[]).unwrap(),
                compute_lr_stretch_polynomial(outer, inner, content).unwrap(),
                "outer={outer:?} inner={inner:?} content={content:?}"
            );
        }
    }

    #[test]
    fn empty_positive_stretch_has_zero_h_vector() {
        let polynomial = compute_lr_stretch_polynomial(&[5, 1], &[2, 1], &[2, 1]).unwrap();
        assert_eq!(polynomial.h_vector, vec![BigInt::zero()]);
        assert_eq!(
            evaluate_h_vector(&polynomial.h_vector, polynomial.dimension, 100),
            BigInt::zero()
        );
    }
}
