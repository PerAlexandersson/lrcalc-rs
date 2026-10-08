//! Exact affine structure of LR, beta-LR, and skew Kostka tableau polytopes.
//!
//! The coordinates are the multiplicities `m[s][r]` of label `s + 1` in row
//! `r` of a skew shape `outer / inner`.  For Littlewood--Richardson tableaux
//! with a virtual Yamanouchi prefix `beta` (`beta = 0` gives ordinary LR
//! tableaux), the rational polytope is cut out by
//!
//! ```text
//! sum_s m[s][r] = outer_r - inner_r                      (row lengths)
//! sum_r m[s][r] = content_s                               (content)
//! m[s][r] >= 0                                            (lower)
//! inner_r + sum_{s' <= s} m[s'][r]
//!     <= inner_{r-1} + sum_{s' < s} m[s'][r-1],   r > 0   (diagonal)
//! beta_s + sum_{r' <= r} m[s][r']
//!     <= beta_{s-1} + sum_{r' < r} m[s-1][r'],    s > 0   (Yamanouchi)
//! ```
//!
//! The diagonal inequalities are column strictness, and the Yamanouchi
//! inequalities are the lattice-word conditions after reading row `r` from
//! right to left.  Skew Kostka polytopes omit the Yamanouchi family.
//!
//! These are exactly the inequality families whose tightness the Buch,
//! GT-chain, and packed Kostka counters record.  Earlier versions inferred
//! tightness from the lattice points at one dilation.  Those points need not
//! affinely span the rational polytope, so the inferred dimension could be
//! too small and the interior counts too large.  We instead compute the
//! implicit equalities exactly with [`crate::affine_hull`], whose result is
//! certified by an exact relative-interior point and dual certificate.  The
//! result is invariant under dilation.
//!
//! Dimension zero means a single rational point, not necessarily a lattice
//! point.  The stretching interpolation treats a zero-dimensional family as
//! the constant one.  If the polytope has a lattice point at dilation one,
//! that point is the whole polytope and the assumption holds.  This is so
//! for ordinary LR coefficients by saturation (Knutson--Tao) and for
//! straight-shape Kostka numbers, whose nonemptiness is the dominance
//! condition.  For other families (skew Kostka, general beta-LR) it is an
//! assumption of the interpolation, as is polynomiality.

use crate::affine_hull::{rational, AffineHull, RationalPolyhedron};

/// Implicit equalities and dimension of an LR or skew Kostka polytope.
///
/// `lower[s][r]`, `diagonal[s][r]` and `yamanouchi[s][r]` are true when the
/// corresponding inequality holds with equality on the whole polytope.
/// Entries without an inequality (`diagonal[s][0]`, `yamanouchi[0][r]`, and
/// all Yamanouchi entries of a Kostka polytope) are reported as tight, so a
/// counter never asks for strictness there.
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct ExactTightFlags {
    pub(crate) dimension: usize,
    pub(crate) lower: Vec<Vec<bool>>,
    pub(crate) diagonal: Vec<Vec<bool>>,
    pub(crate) yamanouchi: Vec<Vec<bool>>,
}

/// Exact implicit equalities of the tableau polytope, or `None` if it is
/// empty as a rational polytope.
///
/// `steps` is the number of labels.  Entries of `inner`, `content`, and `beta`
/// beyond their lengths are zero.  `beta = None` omits the Yamanouchi
/// inequalities (skew Kostka); `Some(&[])` gives ordinary LR tableaux.
pub(crate) fn exact_tight_flags(
    outer: &[i64],
    inner: &[i64],
    content: &[i64],
    steps: usize,
    beta: Option<&[i64]>,
) -> Option<ExactTightFlags> {
    let rows = outer.len();
    let entry = |values: &[i64], index: usize| values.get(index).copied().unwrap_or(0);
    let var = |step: usize, row: usize| step * rows + row;
    let mut system = RationalPolyhedron::new(rows * steps);

    for (row, &outer_part) in outer.iter().enumerate() {
        system.add_equality(
            (0..steps).map(|step| (var(step, row), rational(1))),
            rational(outer_part - entry(inner, row)),
        );
    }
    for step in 0..steps {
        system.add_equality(
            (0..rows).map(|row| (var(step, row), rational(1))),
            rational(entry(content, step)),
        );
    }
    if content.len() > steps && content[steps..].iter().any(|&part| part != 0) {
        return None;
    }

    let mut lower = vec![vec![None; rows]; steps];
    let mut diagonal = vec![vec![None; rows]; steps];
    let mut yamanouchi = vec![vec![None; rows]; steps];
    for step in 0..steps {
        for row in 0..rows {
            lower[step][row] =
                Some(system.add_inequality([(var(step, row), rational(-1))], rational(0)));
            if row > 0 {
                let terms = (0..=step)
                    .map(|earlier| (var(earlier, row), rational(1)))
                    .chain((0..step).map(|earlier| (var(earlier, row - 1), rational(-1))));
                diagonal[step][row] = Some(
                    system
                        .add_inequality(terms, rational(entry(inner, row - 1) - entry(inner, row))),
                );
            }
            if let Some(beta) = beta {
                if step > 0 {
                    let terms = (0..=row)
                        .map(|prefix_row| (var(step, prefix_row), rational(1)))
                        .chain(
                            (0..row).map(|prefix_row| (var(step - 1, prefix_row), rational(-1))),
                        );
                    yamanouchi[step][row] = Some(system.add_inequality(
                        terms,
                        rational(entry(beta, step - 1) - entry(beta, step)),
                    ));
                }
            }
        }
    }

    let AffineHull::Nonempty(hull) = system.affine_hull() else {
        return None;
    };
    let tight = |indices: Vec<Vec<Option<usize>>>| {
        indices
            .into_iter()
            .map(|row| {
                row.into_iter()
                    .map(|index| index.is_none_or(|index| hull.implicit_equalities[index]))
                    .collect()
            })
            .collect()
    };
    Some(ExactTightFlags {
        dimension: hull.dimension,
        lower: tight(lower),
        diagonal: tight(diagonal),
        yamanouchi: tight(yamanouchi),
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn translated_lr_witness_has_dimension_four() {
        // c^{(7,6,5,4,2)}_{(6,5,4,2),(3,2,2)}: every tableau at dilation one
        // omits 3 from the third row, but dilation two has one with a 3 there.
        let flags =
            exact_tight_flags(&[7, 6, 5, 4, 2], &[6, 5, 4, 2], &[3, 2, 2], 3, Some(&[])).unwrap();
        assert_eq!(flags.dimension, 4);
        assert!(!flags.lower[2][2]);
    }

    #[test]
    fn kostka_segment_has_jointly_forced_equalities() {
        // Shape (3,2,1), content (1,1,3,1): the polytope is a segment.
        let flags = exact_tight_flags(&[3, 2, 1], &[], &[1, 1, 3, 1], 4, None).unwrap();
        assert_eq!(flags.dimension, 1);
        let empty = exact_tight_flags(&[3, 3, 2, 1], &[], &[1, 1, 2, 4, 1], 5, None);
        assert_eq!(empty, None);
    }

    #[test]
    fn dimension_is_dilation_invariant() {
        // c^{(3,2,1)}_{(2,1),(2,1)} = 2.
        let base = exact_tight_flags(&[3, 2, 1], &[2, 1], &[2, 1], 2, Some(&[])).unwrap();
        let tripled = exact_tight_flags(&[9, 6, 3], &[6, 3], &[6, 3], 2, Some(&[])).unwrap();
        assert_eq!(base.dimension, 1);
        assert_eq!(base, tripled);
    }
}
