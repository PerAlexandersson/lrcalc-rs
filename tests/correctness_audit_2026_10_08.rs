//! Regressions from the 2026-10-08 correctness audit of LR, beta-LR, and
//! packed Kostka dimensions, interiors, beta column strictness, and h-vector
//! evaluation.
//!
//! The expected values are mathematical answers.  Stretching polynomials are
//! checked against direct counts at fresh dilations that were not used for
//! interpolation, so agreement between two engines that share a dimension
//! routine is never the only evidence.

use lrcalc::{
    kostka::kostka_lr_triple,
    kostka_fast::{skew_kostka_fast_u128, skew_kostka_interior_u128},
    lr_ehrhart::{
        compute_beta_lr_stretch_polynomial, compute_lr_stretch_polynomial, evaluate_h_vector,
        lr_stretch_coefficient, LrStretchPolynomial,
    },
    lr_gt::{lrcoef_gt_dimension, lrcoef_gt_interior_u128, lrcoef_gt_u128},
    lrcoef::{
        beta_lr_content_expansion, beta_lrcoef, beta_lrcoef_buch_counts_u128,
        beta_lrcoef_buch_dimension, lrcoef, lrcoef_buch_dimension, lrcoef_buch_interior_u128,
    },
};
use num_bigint::BigInt;
use num_rational::BigRational;
use num_traits::{One, Zero};

fn scaled(parts: &[i32], dilation: i32) -> Vec<i32> {
    parts.iter().map(|part| part * dilation).collect()
}

fn evaluate(polynomial: &LrStretchPolynomial, dilation: i64) -> BigRational {
    let x = BigRational::from_integer(BigInt::from(dilation));
    polynomial
        .coefficients
        .iter()
        .rev()
        .fold(BigRational::zero(), |value, coefficient| {
            value * &x + coefficient
        })
}

fn rational(value: u128) -> BigRational {
    BigRational::from_integer(BigInt::from(value))
}

fn binomial(top: i64, k: i64) -> BigRational {
    (0..k).fold(BigRational::one(), |value, i| {
        value * BigRational::from_integer(BigInt::from(top - i))
            / BigRational::from_integer(BigInt::from(i + 1))
    })
}

fn partitions(total: i32, max_part: i32, max_parts: usize) -> Vec<Vec<i32>> {
    if total == 0 {
        return vec![Vec::new()];
    }
    if max_parts == 0 {
        return Vec::new();
    }
    let mut out = Vec::new();
    for first in (1..=total.min(max_part)).rev() {
        for mut rest in partitions(total - first, first, max_parts - 1) {
            rest.insert(0, first);
            out.push(rest);
        }
    }
    out
}

fn compositions(total: i32, max_parts: usize) -> Vec<Vec<i32>> {
    if total == 0 {
        return vec![Vec::new()];
    }
    if max_parts == 0 {
        return Vec::new();
    }
    let mut out = Vec::new();
    for first in 1..=total {
        for mut rest in compositions(total - first, max_parts - 1) {
            rest.insert(0, first);
            out.push(rest);
        }
    }
    out
}

/// `c^{(7,6,5,4,2)}_{(6,5,4,2),(3,2,2)}` is the translated Kostka number
/// `K_{(3,2,2),(1,1,1,2,2)}`.  Its five tableaux at dilation one omit label 3
/// from the third row, but dilation two has a tableau with one there, so the
/// affine dimension is four and the stretching polynomial is `binom(N+4, 4)`.
#[test]
fn ordinary_lr_stretch_predicts_fresh_counts() {
    let (outer, inner, content) = (&[7, 6, 5, 4, 2][..], &[6, 5, 4, 2][..], &[3, 2, 2][..]);
    assert_eq!(
        lrcoef_buch_dimension(outer, inner, content).unwrap(),
        Some(4)
    );
    assert_eq!(lrcoef_gt_dimension(outer, inner, content).unwrap(), Some(4));
    for dilation in 2..=3 {
        let (a, b, c) = (
            scaled(outer, dilation),
            scaled(inner, dilation),
            scaled(content, dilation),
        );
        assert_eq!(lrcoef_buch_dimension(&a, &b, &c).unwrap(), Some(4));
        assert_eq!(lrcoef_gt_dimension(&a, &b, &c).unwrap(), Some(4));
    }

    let polynomial = compute_lr_stretch_polynomial(outer, inner, content).unwrap();
    assert_eq!(polynomial.dimension, 4);
    assert_eq!(
        polynomial.h_vector,
        vec![BigInt::one(), 0.into(), 0.into(), 0.into(), 0.into()]
    );
    for dilation in 0..=6 {
        let expected = binomial(i64::from(dilation) + 4, 4);
        assert_eq!(evaluate(&polynomial, i64::from(dilation)), expected);
        if dilation > 0 {
            let (a, b, c) = (
                scaled(outer, dilation),
                scaled(inner, dilation),
                scaled(content, dilation),
            );
            assert_eq!(rational(lrcoef(&a, &b, &c).unwrap()), expected);
            assert_eq!(rational(lrcoef_gt_u128(&a, &b, &c).unwrap()), expected);
            // Reciprocity: the interior count is binom(N - 1, 4).
            let interior = binomial(i64::from(dilation) - 1, 4);
            assert_eq!(
                rational(lrcoef_buch_interior_u128(&a, &b, &c).unwrap()),
                interior
            );
            assert_eq!(
                rational(lrcoef_gt_interior_u128(&a, &b, &c).unwrap()),
                interior
            );
        }
    }
    assert_eq!(
        lr_stretch_coefficient(outer, inner, content, 10).unwrap(),
        BigInt::from(1001)
    );

    // The same family through empty-inner beta-LR with vacuous beta.
    let beta = [32, 24, 16, 8];
    let beta_polynomial =
        compute_beta_lr_stretch_polynomial(&[3, 2, 2], &[], &[1, 1, 1, 2, 2], &beta).unwrap();
    assert_eq!(beta_polynomial.dimension, 4);
    assert_eq!(evaluate(&beta_polynomial, 3), binomial(7, 4));
    assert_eq!(
        beta_lrcoef_buch_dimension(&[6, 4, 4], &[], &[2, 2, 2, 4, 4], &scaled(&beta, 2)).unwrap(),
        Some(4)
    );
}

/// With beta `(8, 4)` the Yamanouchi inequalities are vacuous, so the only
/// filling of a three-box column with content `(1,1,1)` is `1,2,3`.
#[test]
fn dominant_beta_counts_only_column_strict_tableaux() {
    assert_eq!(
        beta_lrcoef(&[1, 1, 1], &[], &[1, 1, 1], &[8, 4]).unwrap(),
        1
    );
    assert_eq!(
        beta_lrcoef_buch_counts_u128(&[1, 1, 1], &[], &[1, 1, 1], &[8, 4])
            .unwrap()
            .full,
        1
    );
    assert_eq!(
        skew_kostka_fast_u128(&[1, 1, 1], &[], &[1, 1, 1]).unwrap(),
        1
    );
    let expansion = beta_lr_content_expansion(&[1, 1, 1], &[], &[8, 4], Some(3)).unwrap();
    let column = expansion
        .iter()
        .find(|term| term.content == [1, 1, 1])
        .expect("content (1,1,1) occurs");
    assert_eq!(column.coefficient, 1);

    assert_eq!(
        skew_kostka_fast_u128(&[2, 2, 2], &[], &[1, 1, 1, 1, 2]).unwrap(),
        2
    );
    assert_eq!(
        beta_lrcoef(&[2, 2, 2], &[], &[1, 1, 1, 1, 2], &[28, 21, 14, 7]).unwrap(),
        2
    );
    assert_eq!(
        beta_lrcoef_buch_counts_u128(&[2, 2, 2], &[], &[1, 1, 1, 1, 2], &[28, 21, 14, 7])
            .unwrap()
            .full,
        2
    );
}

/// Straight shapes and all compositions, including every content permutation.
/// A dominant beta makes the beta-LR count a Kostka number, and the
/// diagonal-row translation makes it an ordinary LR coefficient.
#[test]
fn dominant_beta_and_translated_kostka_grid_matches_fresh_counts() {
    let mut checked = 0;
    for size in 1..=6 {
        for shape in partitions(size, size, 4) {
            for weight in compositions(size, 5) {
                let kostka = skew_kostka_fast_u128(&shape, &[], &weight).unwrap();
                let gap = size + 1;
                let beta = (0..weight.len())
                    .map(|index| (weight.len() - 1 - index) as i32 * gap)
                    .collect::<Vec<_>>();
                assert_eq!(beta_lrcoef(&shape, &[], &weight, &beta).unwrap(), kostka);
                assert_eq!(
                    beta_lrcoef_buch_counts_u128(&shape, &[], &weight, &beta)
                        .unwrap()
                        .full,
                    kostka
                );
                if kostka == 0 || kostka > 30 {
                    continue;
                }
                let beta_polynomial =
                    compute_beta_lr_stretch_polynomial(&shape, &[], &weight, &beta).unwrap();
                let dimension = beta_polynomial.dimension;
                if dimension > 4 {
                    continue;
                }
                let fresh = dimension as i32 + 2;
                let direct =
                    skew_kostka_fast_u128(&scaled(&shape, fresh), &[], &scaled(&weight, fresh))
                        .unwrap();
                assert_eq!(
                    evaluate(&beta_polynomial, i64::from(fresh)),
                    rational(direct),
                    "{shape:?} {weight:?}"
                );
                let sign = if dimension.is_multiple_of(2) { 1 } else { -1 };
                let interior =
                    skew_kostka_interior_u128(&scaled(&shape, 2), &[], &scaled(&weight, 2))
                        .unwrap();
                assert_eq!(
                    evaluate(&beta_polynomial, -2) * BigRational::from_integer(sign.into()),
                    rational(interior),
                    "{shape:?} {weight:?}"
                );
                assert_eq!(
                    beta_lrcoef_buch_dimension(
                        &scaled(&shape, 2),
                        &[],
                        &scaled(&weight, 2),
                        &scaled(&beta, 2)
                    )
                    .unwrap(),
                    Some(dimension)
                );

                let (outer, inner, content) = kostka_lr_triple(&shape, &weight).unwrap();
                let ordinary = compute_lr_stretch_polynomial(&outer, &inner, &content).unwrap();
                assert_eq!(ordinary.dimension, dimension, "{shape:?} {weight:?}");
                assert_eq!(
                    lrcoef_gt_dimension(&outer, &inner, &content).unwrap(),
                    Some(dimension)
                );
                assert_eq!(
                    evaluate(&ordinary, i64::from(fresh)),
                    rational(direct),
                    "{shape:?} {weight:?}"
                );
                checked += 1;
            }
        }
    }
    assert!(checked > 100);
}

/// Every LR triple with outer size at most seven and at most four rows.
#[test]
fn ordinary_lr_grid_matches_fresh_counts() {
    let mut checked = 0;
    for size in 1..=7 {
        for outer in partitions(size, size, 4) {
            for inner_size in 0..=size {
                for inner in partitions(inner_size, inner_size, 4) {
                    if inner.len() > outer.len()
                        || inner.iter().zip(&outer).any(|(small, large)| small > large)
                    {
                        continue;
                    }
                    for content in partitions(size - inner_size, size - inner_size, 4) {
                        if lrcoef(&outer, &inner, &content).unwrap() == 0 {
                            assert_eq!(
                                lrcoef_buch_dimension(&outer, &inner, &content).unwrap(),
                                None
                            );
                            continue;
                        }
                        let polynomial =
                            compute_lr_stretch_polynomial(&outer, &inner, &content).unwrap();
                        let dimension = polynomial.dimension;
                        assert_eq!(
                            lrcoef_gt_dimension(&outer, &inner, &content).unwrap(),
                            Some(dimension)
                        );
                        let (a, b, c) = (scaled(&outer, 2), scaled(&inner, 2), scaled(&content, 2));
                        assert_eq!(lrcoef_buch_dimension(&a, &b, &c).unwrap(), Some(dimension));
                        let fresh = dimension as i32 + 2;
                        let direct = lrcoef(
                            &scaled(&outer, fresh),
                            &scaled(&inner, fresh),
                            &scaled(&content, fresh),
                        )
                        .unwrap();
                        assert_eq!(
                            evaluate(&polynomial, i64::from(fresh)),
                            rational(direct),
                            "{outer:?} {inner:?} {content:?}"
                        );
                        let sign = if dimension.is_multiple_of(2) { 1 } else { -1 };
                        let expected =
                            evaluate(&polynomial, -2) * BigRational::from_integer(sign.into());
                        assert_eq!(
                            rational(lrcoef_buch_interior_u128(&a, &b, &c).unwrap()),
                            expected
                        );
                        assert_eq!(
                            rational(lrcoef_gt_interior_u128(&a, &b, &c).unwrap()),
                            expected
                        );
                        checked += 1;
                    }
                }
            }
        }
    }
    assert!(checked > 100);
}

/// Shape `(3,2,1)`, content `(1,1,3,1)` is a segment with `N + 1` points:
/// three upper bounds sum to a level total and are jointly forced.
#[test]
fn packed_kostka_interiors_use_exact_affine_hull() {
    for dilation in 1..=4 {
        let shape = scaled(&[3, 2, 1], dilation);
        let weight = scaled(&[1, 1, 3, 1], dilation);
        assert_eq!(
            skew_kostka_fast_u128(&shape, &[], &weight).unwrap(),
            (dilation + 1) as u128
        );
        assert_eq!(
            skew_kostka_interior_u128(&shape, &[], &weight).unwrap(),
            (dilation - 1) as u128
        );
    }
}

#[test]
fn h_vector_evaluation_supports_the_full_u64_domain() {
    let max = BigInt::from(u64::MAX);
    assert_eq!(evaluate_h_vector(&[BigInt::one()], 1, u64::MAX), &max + 1);
    // binom(N + 4, 4) at N = u64::MAX.
    let h = [BigInt::one(), 0.into(), 0.into(), 0.into(), 0.into()];
    let expected = (1..=4).fold(BigInt::one(), |value, i| value * (&max + i)) / 24;
    assert_eq!(evaluate_h_vector(&h, 4, u64::MAX), expected);
    assert_eq!(evaluate_h_vector(&h, 4, 0), BigInt::one());
}
