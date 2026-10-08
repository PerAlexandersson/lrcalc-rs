//! Broader deterministic checks, including non-dominant beta and zero entries.
use lrcalc::{
    kostka_fast::{skew_kostka_fast_u128, skew_kostka_interior_u128, KostkaFastError},
    lr_ehrhart::{compute_beta_lr_stretch_polynomial, evaluate_h_vector},
    lr_gt::{lrcoef_gt_dimension, lrcoef_gt_u128},
    lrcoef::{
        beta_lr_content_expansion, beta_lrcoef, beta_lrcoef_buch_dimension, lrcoef,
        lrcoef_buch_dimension,
    },
};
use num_bigint::BigInt;

fn partitions(total: i32, max: i32, rows: usize) -> Vec<Vec<i32>> {
    if total == 0 {
        return vec![vec![]];
    }
    if rows == 0 {
        return vec![];
    }
    let mut out = vec![];
    for first in 1..=max.min(total) {
        for mut tail in partitions(total - first, first, rows - 1) {
            tail.insert(0, first);
            out.push(tail);
        }
    }
    out
}

fn scaled(parts: &[i32], n: i32) -> Vec<i32> {
    parts.iter().map(|x| x * n).collect()
}

#[test]
fn lr_symmetries_and_zero_padding() {
    let mut checked = 0;
    for size in 0..=8 {
        for outer in partitions(size, size, 4) {
            for inner_size in 0..=size {
                for inner in partitions(inner_size, inner_size, 4) {
                    for content in partitions(size - inner_size, size, 4) {
                        let count = lrcoef(&outer, &inner, &content).unwrap();
                        assert_eq!(count, lrcoef(&outer, &content, &inner).unwrap());
                        // Tiny test oracle: no public conjugation API in this crate.
                        let conjugate = |v: &[i32]| {
                            (1..=v.first().copied().unwrap_or(0))
                                .map(|j| v.iter().filter(|&&x| x >= j).count() as i32)
                                .collect::<Vec<_>>()
                        };
                        assert_eq!(
                            count,
                            lrcoef(&conjugate(&outer), &conjugate(&inner), &conjugate(&content))
                                .unwrap()
                        );
                        let pad = |v: &[i32]| {
                            let mut p = v.to_vec();
                            p.extend([0, 0]);
                            p
                        };
                        assert_eq!(
                            count,
                            lrcoef(&pad(&outer), &pad(&inner), &pad(&content)).unwrap()
                        );
                        assert_eq!(count, lrcoef_gt_u128(&outer, &inner, &content).unwrap());
                        assert_eq!(
                            lrcoef_buch_dimension(&outer, &inner, &content).unwrap(),
                            lrcoef_gt_dimension(&outer, &inner, &content).unwrap()
                        );
                        checked += 1;
                    }
                }
            }
        }
    }
    eprintln!("LR symmetry/padding triples: {checked}");
}

#[test]
fn beta_with_zero_entries_matches_expansion_and_fresh_stretches() {
    let betas = [
        vec![],
        vec![1],
        vec![1, 1],
        vec![2, 1],
        vec![2, 2],
        vec![3, 1],
        vec![3, 2, 1],
    ];
    let mut checked = 0;
    for size in 0..=5 {
        for outer in partitions(size, size, 3) {
            for beta in &betas {
                let expansion = beta_lr_content_expansion(&outer, &[], beta, Some(3)).unwrap();
                for a in 0..=size {
                    for b in 0..=size - a {
                        let content = [a, b, size - a - b];
                        let scalar = beta_lrcoef(&outer, &[], &content, beta).unwrap();
                        let expected = expansion
                            .iter()
                            .find(|term| {
                                (0..3).all(|i| {
                                    term.content.get(i).copied().unwrap_or(0) == content[i]
                                })
                            })
                            .map_or(0, |term| term.coefficient);
                        assert_eq!(
                            scalar, expected,
                            "outer={outer:?} content={content:?} beta={beta:?}"
                        );
                        let polynomial =
                            compute_beta_lr_stretch_polynomial(&outer, &[], &content, beta)
                                .unwrap();
                        for n in [2, 3, 5] {
                            let direct = beta_lrcoef(
                                &scaled(&outer, n),
                                &[],
                                &scaled(&content, n),
                                &scaled(beta, n),
                            )
                            .unwrap();
                            assert_eq!(
                                evaluate_h_vector(
                                    &polynomial.h_vector,
                                    polynomial.dimension,
                                    n as u64
                                ),
                                BigInt::from(direct),
                                "outer={outer:?} content={content:?} beta={beta:?} n={n}"
                            );
                        }
                        assert_eq!(
                            beta_lrcoef_buch_dimension(&outer, &[], &content, beta).unwrap(),
                            beta_lrcoef_buch_dimension(
                                &scaled(&outer, 2),
                                &[],
                                &scaled(&content, 2),
                                &scaled(beta, 2)
                            )
                            .unwrap()
                        );
                        checked += 1;
                    }
                }
            }
        }
    }
    eprintln!("non-dominant beta/zero-content cases: {checked}");
}

#[test]
fn packed_bit_budget_is_checked_at_128_bits() {
    for (rows, width, fits) in [
        (8, 65535, true),
        (8, 65536, false),
        (14, 511, true),
        (15, 511, false),
    ] {
        let outer = vec![width; rows];
        let mut inner = outer.clone();
        inner[rows - 1] -= 1;
        let result = skew_kostka_fast_u128(&outer, &inner, &[1]);
        if fits {
            assert_eq!(result.unwrap(), 1);
        } else {
            assert_eq!(result, Err(KostkaFastError::StateTooWide));
        }
    }
}

#[test]
fn extreme_shape_totals_are_exact_or_checked_errors() {
    let m = i32::MAX;
    for rows in 1..=4 {
        let outer = vec![m; rows];
        let mut inner = outer.clone();
        inner[rows - 1] -= 1;
        if let Ok(count) = skew_kostka_fast_u128(&outer, &inner, &[1]) {
            assert_eq!(count, 1);
        }
        if let Ok(count) = lrcoef(&outer, &inner, &[1]) {
            assert_eq!(count, 1);
        }
        if let Ok(count) = lrcoef_gt_u128(&outer, &inner, &[1]) {
            assert_eq!(count, 1);
        }
    }
}

#[test]
fn kostka_zero_content_insertion_preserves_full_and_interior_counts() {
    for n in 1..=4 {
        let shape = scaled(&[3, 2, 1], n);
        for slot in 0..=4 {
            let mut content = scaled(&[1, 1, 3, 1], n);
            content.insert(slot, 0);
            assert_eq!(
                skew_kostka_fast_u128(&shape, &[], &content).unwrap(),
                (n + 1) as u128
            );
            assert_eq!(
                skew_kostka_interior_u128(&shape, &[], &content).unwrap(),
                (n - 1) as u128
            );
        }
    }
}
