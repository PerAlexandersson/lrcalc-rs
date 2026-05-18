//! Schubert polynomial multiplication routines.
//!
//! These are direct Rust ports of the recursive routines in upstream
//! `schublib.c`, using owned vectors and maps internally instead of C
//! `ivlincomb` storage.

use std::collections::BTreeMap;

pub type LinearCombination = BTreeMap<Vec<i32>, i32>;
const MAX_STRING_CLASSES: usize = 1_000_000;

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum SchubertError {
    InvalidInput,
    ArithmeticOverflow,
}

pub fn trans_polynomial(w: &[i32], vars: i32) -> Result<LinearCombination, SchubertError> {
    if vars < 0 {
        return Err(SchubertError::InvalidInput);
    }
    let mut permutation = w.to_vec();
    trans_recursive(
        &mut permutation,
        usize::try_from(vars).map_err(|_| SchubertError::ArithmeticOverflow)?,
    )
}

pub fn monk_product(
    i: i32,
    slc: &LinearCombination,
    rank: i32,
) -> Result<LinearCombination, SchubertError> {
    if i <= 0 || rank < 0 {
        return Err(SchubertError::InvalidInput);
    }
    let mut out = LinearCombination::new();
    let rank = if rank == 0 { i32::MAX } else { rank };
    monk_add(
        usize::try_from(i).map_err(|_| SchubertError::ArithmeticOverflow)?,
        slc,
        rank,
        &mut out,
    )?;
    Ok(out)
}

pub fn multiply_poly_schubert(
    poly: &LinearCombination,
    perm: &[i32],
    rank: i32,
) -> Result<LinearCombination, SchubertError> {
    if rank < 0 {
        return Err(SchubertError::InvalidInput);
    }
    if poly.is_empty() {
        return Ok(LinearCombination::new());
    }
    let rank = if rank == 0 { i32::MAX } else { rank };
    let mut terms = Vec::with_capacity(poly.len());
    let mut maxvar = 0usize;
    for (monomial, &coefficient) in poly {
        let monomial = trim_trailing_zeroes(monomial);
        maxvar = maxvar.max(monomial.len());
        terms.push((monomial, coefficient));
    }
    let mut perm = perm.to_vec();
    perm.truncate(perm_group(&perm));
    let mut out = LinearCombination::new();
    multiply_poly_schubert_recursive(terms, maxvar, &perm, rank, &mut out)?;
    Ok(out)
}

pub fn multiply_schubert(
    w1: &[i32],
    w2: &[i32],
    rank: i32,
) -> Result<LinearCombination, SchubertError> {
    if rank < 0 || !valid_permutation(w1) || !valid_permutation(w2) {
        return Err(SchubertError::InvalidInput);
    }
    let mut left = w1.to_vec();
    let mut right = w2.to_vec();
    let mut left_len = perm_length(&left);
    let mut right_len = perm_length(&right);
    if left_len > right_len {
        std::mem::swap(&mut left, &mut right);
        std::mem::swap(&mut left_len, &mut right_len);
    }

    left.truncate(perm_group(&left));
    right.truncate(perm_group(&right));

    let rank = if rank == 0 { i32::MAX } else { rank };
    if rank != i32::MAX
        && (twice_sum_exceeds_rank(left_len, right_len, rank) || bruhat_zero(&left, &right, rank))
    {
        return Ok(LinearCombination::new());
    }

    let poly = trans_polynomial(&left, 0)?;
    multiply_poly_schubert(&poly, &right, rank)
}

pub fn multiply_schubert_strings(
    str1: &[i32],
    str2: &[i32],
) -> Result<LinearCombination, SchubertError> {
    if !strings_are_compatible(str1, str2) {
        return Err(SchubertError::InvalidInput);
    }
    let dimvec = string_dimension_vector(str1).ok_or(SchubertError::InvalidInput)?;
    let w1 = string_to_permutation(str1).ok_or(SchubertError::InvalidInput)?;
    let w2 = string_to_permutation(str2).ok_or(SchubertError::InvalidInput)?;
    let product = multiply_schubert(
        &w1,
        &w2,
        i32::try_from(w1.len()).map_err(|_| SchubertError::ArithmeticOverflow)?,
    )?;

    let mut out = LinearCombination::new();
    for (perm, coefficient) in product {
        let string = permutation_to_string(&perm, &dimvec).ok_or(SchubertError::InvalidInput)?;
        add_term(&mut out, string, coefficient)?;
    }
    Ok(out)
}

pub fn valid_permutation(w: &[i32]) -> bool {
    let n = w.len();
    let mut seen = vec![false; n];
    for &value in w {
        let Some(value) = value.checked_sub(1) else {
            return false;
        };
        let Ok(index) = usize::try_from(value) else {
            return false;
        };
        if index >= n || seen[index] {
            return false;
        }
        seen[index] = true;
    }
    true
}

pub fn strings_are_compatible(str1: &[i32], str2: &[i32]) -> bool {
    if str1.len() != str2.len() {
        return false;
    }
    string_dimension_vector(str1) == string_dimension_vector(str2)
}

fn trans_recursive(w: &mut Vec<i32>, mut vars: usize) -> Result<LinearCombination, SchubertError> {
    let n = perm_group(w);
    w.truncate(n);

    let mut r = n.saturating_sub(1);
    while r > 0 && w[r - 1] < w[r] {
        r -= 1;
    }
    if r == 0 {
        let key_len = vars.max(1);
        let mut out = LinearCombination::new();
        out.insert(vec![0; key_len], 1);
        return Ok(out);
    }
    vars = vars.max(r);

    let mut s = r + 1;
    while s < n && w[r - 1] > w[s] {
        s += 1;
    }

    let wr = w[r - 1];
    let ws = w[s - 1];
    let mut v = w.clone();
    v[s - 1] = wr;
    v[r - 1] = ws;

    let mut out = LinearCombination::new();
    for (mut monomial, coefficient) in trans_recursive(&mut v.clone(), vars)? {
        monomial[r - 1] = monomial[r - 1]
            .checked_add(1)
            .ok_or(SchubertError::ArithmeticOverflow)?;
        add_term(&mut out, monomial, coefficient)?;
    }

    let mut last = 0;
    let vr = v[r - 1];
    for i in (1..r).rev() {
        let vi = v[i - 1];
        if last < vi && vi < vr {
            last = vi;
            let mut next = v.clone();
            next[i - 1] = vr;
            next[r - 1] = vi;
            let tmp = trans_recursive(&mut next, vars)?;
            add_multiple(&mut out, &tmp, 1)?;
        }
    }

    Ok(out)
}

fn monk_add(
    i: usize,
    slc: &LinearCombination,
    rank: i32,
    out: &mut LinearCombination,
) -> Result<(), SchubertError> {
    for (w, &coefficient) in slc {
        let n = w.len();
        let wi = if i <= n {
            w[i - 1]
        } else {
            i32::try_from(i).map_err(|_| SchubertError::ArithmeticOverflow)?
        };

        if i <= n + 1 {
            let ulen = i.max(n);
            let mut last = 0;
            for j in (1..i).rev() {
                if last < w[j - 1] && w[j - 1] < wi {
                    last = w[j - 1];
                    let mut u = extend_permutation(w, ulen)?;
                    u[j - 1] = wi;
                    u[i - 1] = last;
                    add_term(
                        out,
                        u,
                        coefficient
                            .checked_neg()
                            .ok_or(SchubertError::ArithmeticOverflow)?,
                    )?;
                }
            }
        } else {
            let mut u = extend_permutation(w, i)?;
            u[i - 2] = i32::try_from(i).map_err(|_| SchubertError::ArithmeticOverflow)?;
            u[i - 1] = i32::try_from(i - 1).map_err(|_| SchubertError::ArithmeticOverflow)?;
            add_term(
                out,
                u,
                coefficient
                    .checked_neg()
                    .ok_or(SchubertError::ArithmeticOverflow)?,
            )?;
        }

        if i >= n + 1 {
            let mut u = extend_permutation(w, i + 1)?;
            u[i - 1] = i32::try_from(i + 1).map_err(|_| SchubertError::ArithmeticOverflow)?;
            u[i] = i32::try_from(i).map_err(|_| SchubertError::ArithmeticOverflow)?;
            add_term(out, u, coefficient)?;
        } else {
            let mut last = i32::MAX;
            for j in i + 1..=n {
                if wi < w[j - 1] && w[j - 1] < last {
                    last = w[j - 1];
                    let mut u = w.clone();
                    u[i - 1] = last;
                    u[j - 1] = wi;
                    add_term(out, u, coefficient)?;
                }
            }
            if last > i32::try_from(n).map_err(|_| SchubertError::ArithmeticOverflow)?
                && n < usize::try_from(rank).unwrap_or(usize::MAX)
            {
                let mut u = extend_permutation(w, n + 1)?;
                u[i - 1] = i32::try_from(n + 1).map_err(|_| SchubertError::ArithmeticOverflow)?;
                u[n] = wi;
                add_term(out, u, coefficient)?;
            }
        }
    }
    Ok(())
}

fn multiply_poly_schubert_recursive(
    terms: Vec<(Vec<i32>, i32)>,
    maxvar: usize,
    perm: &[i32],
    rank: i32,
    out: &mut LinearCombination,
) -> Result<(), SchubertError> {
    if terms.is_empty() {
        return Ok(());
    }
    if maxvar == 0 {
        return add_term(out, perm.to_vec(), terms[0].1);
    }

    let mut lower = Vec::new();
    let mut upper = Vec::new();
    let mut mv0 = 0usize;
    let mut mv1 = 0usize;
    for (mut monomial, coefficient) in terms {
        if monomial.len() < maxvar {
            mv0 = mv0.max(monomial.len());
            lower.push((monomial, coefficient));
        } else {
            monomial[maxvar - 1] -= 1;
            monomial = trim_trailing_zeroes(&monomial);
            mv1 = mv1.max(monomial.len());
            upper.push((monomial, coefficient));
        }
    }

    let mut res1 = LinearCombination::new();
    multiply_poly_schubert_recursive(upper, mv1, perm, rank, &mut res1)?;
    monk_add(maxvar, &res1, rank, out)?;
    if !lower.is_empty() {
        multiply_poly_schubert_recursive(lower, mv0, perm, rank, out)?;
    }
    Ok(())
}

fn extend_permutation(w: &[i32], len: usize) -> Result<Vec<i32>, SchubertError> {
    let mut out = Vec::with_capacity(len);
    out.extend_from_slice(w);
    for value in w.len() + 1..=len {
        out.push(i32::try_from(value).map_err(|_| SchubertError::ArithmeticOverflow)?);
    }
    Ok(out)
}

fn perm_length(w: &[i32]) -> i32 {
    let mut length = 0i32;
    for i in 0..w.len().saturating_sub(1) {
        for j in i + 1..w.len() {
            if w[i] > w[j] {
                length += 1;
            }
        }
    }
    length
}

fn perm_group(w: &[i32]) -> usize {
    let mut len = w.len();
    while len > 0 && w[len - 1] == len as i32 {
        len -= 1;
    }
    len
}

fn bruhat_zero(w1: &[i32], w2: &[i32], rank: i32) -> bool {
    let mut n1 = perm_group(w1);
    let n2 = perm_group(w2);
    if n1 > rank as usize || n2 > rank as usize {
        return true;
    }
    let (w1, w2) = if n1 > n2 {
        n1 = n2;
        (w2, w1)
    } else {
        (w1, w2)
    };
    for q in 1..n1 {
        let q2 = rank - q as i32;
        let mut r1 = 0;
        let mut r2 = 0;
        for p in 0..n1.saturating_sub(1) {
            if w1[p] <= q as i32 {
                r1 += 1;
            }
            if w2[p] > q2 {
                r2 += 1;
            }
            if r1 < r2 {
                return true;
            }
        }
    }
    false
}

fn twice_sum_exceeds_rank(left_len: i32, right_len: i32, rank: i32) -> bool {
    let lhs = (i64::from(left_len) + i64::from(right_len)) * 2;
    let rank = i64::from(rank);
    let rhs = rank * (rank - 1);
    lhs > rhs
}

fn string_dimension_vector(string: &[i32]) -> Option<Vec<i32>> {
    let mut classes = 0usize;
    for &value in string {
        let value = usize::try_from(value).ok()?;
        if value >= MAX_STRING_CLASSES {
            return None;
        }
        let next = value.checked_add(1)?;
        classes = classes.max(next);
    }
    let mut out = vec![0; classes];
    for &value in string {
        let index = usize::try_from(value).ok()?;
        out[index] += 1;
    }
    for index in 1..out.len() {
        out[index] += out[index - 1];
    }
    Some(out)
}

fn string_to_permutation(string: &[i32]) -> Option<Vec<i32>> {
    let mut dimvec = string_dimension_vector(string)?;
    let mut perm = vec![0; string.len()];
    for index in (0..string.len()).rev() {
        let class = usize::try_from(string[index]).ok()?;
        dimvec[class] -= 1;
        let target = usize::try_from(dimvec[class]).ok()?;
        perm[target] = i32::try_from(index + 1).ok()?;
    }
    Some(perm)
}

fn permutation_to_string(perm: &[i32], dimvec: &[i32]) -> Option<Vec<i32>> {
    let n = dimvec.last().copied().unwrap_or(0);
    let mut out = vec![0; usize::try_from(n).ok()?];
    let mut j = 0usize;
    for (class, &limit) in dimvec.iter().enumerate() {
        let limit = usize::try_from(limit).ok()?;
        while j < limit {
            let wj = if j < perm.len() {
                perm[j]
            } else {
                i32::try_from(j + 1).ok()?
            };
            let target = usize::try_from(wj.checked_sub(1)?).ok()?;
            if target >= out.len() {
                return None;
            }
            out[target] = i32::try_from(class).ok()?;
            j += 1;
        }
    }
    Some(out)
}

fn add_multiple(
    out: &mut LinearCombination,
    terms: &LinearCombination,
    scale: i32,
) -> Result<(), SchubertError> {
    for (key, &coefficient) in terms {
        let coefficient = coefficient
            .checked_mul(scale)
            .ok_or(SchubertError::ArithmeticOverflow)?;
        add_term(out, key.clone(), coefficient)?;
    }
    Ok(())
}

fn add_term(
    out: &mut LinearCombination,
    key: Vec<i32>,
    coefficient: i32,
) -> Result<(), SchubertError> {
    if coefficient == 0 {
        return Ok(());
    }
    let entry = out.entry(key.clone()).or_insert(0);
    *entry = entry
        .checked_add(coefficient)
        .ok_or(SchubertError::ArithmeticOverflow)?;
    if *entry == 0 {
        out.remove(&key);
    }
    Ok(())
}

fn trim_trailing_zeroes(values: &[i32]) -> Vec<i32> {
    let len = values
        .iter()
        .rposition(|&value| value != 0)
        .map_or(0, |index| index + 1);
    values[..len].to_vec()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn terms(out: LinearCombination) -> Vec<(Vec<i32>, i32)> {
        out.into_iter().collect()
    }

    fn abs_coeff_sum(out: &LinearCombination) -> i32 {
        out.values().map(|coefficient| coefficient.abs()).sum()
    }

    #[test]
    fn computes_basic_transitions() {
        assert_eq!(
            terms(trans_polynomial(&[1, 2], 0).unwrap()),
            vec![(vec![0], 1)]
        );
        assert_eq!(
            terms(trans_polynomial(&[2, 1], 0).unwrap()),
            vec![(vec![1], 1)]
        );
        assert_eq!(
            terms(trans_polynomial(&[2, 3, 1], 0).unwrap()),
            vec![(vec![1, 1], 1)]
        );
    }

    #[test]
    fn multiplies_small_schubert_polynomials() {
        assert_eq!(
            terms(multiply_schubert(&[1, 2], &[1, 2], 0).unwrap()),
            vec![(vec![], 1)]
        );
        assert_eq!(
            terms(multiply_schubert(&[2, 1], &[1, 2], 0).unwrap()),
            vec![(vec![2, 1], 1)]
        );
        assert_eq!(
            terms(multiply_schubert(&[2, 1], &[2, 1], 0).unwrap()),
            vec![(vec![3, 1, 2], 1)]
        );
    }

    #[test]
    fn multiplies_nontrivial_schubert_products() {
        assert_eq!(
            terms(multiply_schubert(&[2, 1, 3], &[1, 3, 2], 0).unwrap()),
            vec![(vec![2, 3, 1], 1), (vec![3, 1, 2], 1)]
        );
        assert_eq!(
            terms(multiply_schubert(&[2, 4, 1, 3], &[3, 1, 4, 2], 0).unwrap()),
            vec![
                (vec![4, 3, 2, 1], 1),
                (vec![4, 5, 1, 2, 3], 1),
                (vec![5, 2, 3, 1, 4], 1),
                (vec![5, 3, 1, 2, 4], 1),
            ]
        );
    }

    #[test]
    fn rank_bound_can_zero_schubert_product() {
        assert!(multiply_schubert(&[2, 1], &[2, 1], 2).unwrap().is_empty());
    }

    #[test]
    fn rejects_extreme_invalid_inputs_without_overflow() {
        assert!(!valid_permutation(&[i32::MIN]));
        assert!(!valid_permutation(&[i32::MAX]));
        assert!(string_dimension_vector(&[-1]).is_none());
        assert!(string_dimension_vector(&[i32::MAX]).is_none());
        assert_eq!(string_dimension_vector(&[0, 2]), Some(vec![1, 1, 2]));
        assert_eq!(string_dimension_vector(&[3]), Some(vec![0, 0, 0, 1]));
        assert!(permutation_to_string(&[i32::MIN], &[1]).is_none());
        assert!(!twice_sum_exceeds_rank(i32::MAX, i32::MAX, i32::MAX));
        let _ = twice_sum_exceeds_rank(i32::MAX, i32::MAX, i32::MIN);
    }

    #[test]
    fn multiplies_compatible_strings() {
        assert_eq!(
            terms(multiply_schubert_strings(&[0, 1], &[1, 0]).unwrap()),
            vec![(vec![1, 0], 1)]
        );
        assert_eq!(
            terms(multiply_schubert_strings(&[3], &[3]).unwrap()),
            vec![(vec![3], 1)]
        );
        let mixed = multiply_schubert_strings(&[0, 2, 1, 0, 2, 1], &[1, 0, 2, 1, 0, 2]).unwrap();
        assert_eq!(mixed.len(), 5);
        assert_eq!(abs_coeff_sum(&mixed), 6);
    }
}
