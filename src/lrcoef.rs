//! Single Littlewood-Richardson coefficient computation.
//!
//! This module ports the optimized single-coefficient path from upstream
//! `lrcalc`: first compact the triple with `optim_coef`, then count LR
//! tableaux with the branch-pruned search from `lrcoef_count`.

use std::collections::HashMap;

use num_rational::BigRational;
use num_traits::Zero;

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum LrCoefError {
    InvalidPartition,
    ArithmeticOverflow,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub(crate) struct SkewShape {
    pub(crate) outer: Vec<i32>,
    pub(crate) inner: Vec<i32>,
    pub(crate) content: Vec<i32>,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub(crate) enum OptimizedCoef {
    Zero,
    One,
    Count(SkewShape),
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct LrBuchInteriorStats {
    pub value: u128,
    pub weak_tableaux: u128,
    pub strict_tableaux: u128,
    pub strict_lower_constraints: usize,
    pub strict_diagonal_constraints: usize,
    pub strict_yamanouchi_constraints: usize,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct LrBuchInteriorMemoStats {
    pub value: u128,
    pub weak_tableaux: u128,
    pub memo_states: usize,
    pub cache_hits: usize,
    pub generated_transitions: usize,
    pub strict_lower_constraints: usize,
    pub strict_diagonal_constraints: usize,
    pub strict_yamanouchi_constraints: usize,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct LrBuchCounts {
    pub full: u128,
    pub interior: u128,
}

pub(crate) struct LrBuchStretchCache {
    outer: Vec<i32>,
    inner: Vec<i32>,
    content: Vec<i32>,
    tight: BuchTightFlags,
    dimension: usize,
}

#[derive(Clone, Copy, Debug, Default)]
struct LrCoefBox {
    value: i32,
    max: i32,
    row: usize,
    north: usize,
    east: usize,
    se_supply: i32,
    se_sz: i32,
    west_sz: i32,
}

#[derive(Clone, Copy, Debug, Default)]
struct LrCoefContent {
    cont: i32,
    supply: i32,
}

#[derive(Clone, Debug)]
struct BuchTightFlags {
    lower: Vec<Vec<bool>>,
    diagonal: Vec<Vec<bool>>,
    yamanouchi: Vec<Vec<bool>>,
}

#[derive(Clone, Debug)]
struct BuchTightData {
    tight: BuchTightFlags,
    weak_tableaux: u128,
}

#[derive(Clone, Debug, Eq, Hash, PartialEq)]
struct BuchMemoState {
    step: usize,
    shape: Vec<u32>,
    previous_prefix: Vec<u32>,
}

#[derive(Clone, Debug, Default, Eq, PartialEq)]
struct BuchMemoCounters {
    cache_hits: usize,
    generated_transitions: usize,
}

struct BuchMemoCounter<'a> {
    outer: Vec<u32>,
    content: Vec<u32>,
    tight: &'a BuchTightFlags,
    cache: HashMap<BuchMemoState, u128>,
    counters: BuchMemoCounters,
}

struct BuchInteriorPruner<'a> {
    tight: &'a BuchTightFlags,
    inner: Vec<u32>,
    increments: Vec<u32>,
    rows: usize,
    steps: usize,
}

/// Compute a single Littlewood-Richardson coefficient.
///
/// This returns `Err` only for invalid input partitions or integer overflow in
/// size arithmetic. A valid coefficient that vanishes is returned as `Ok(0)`.
pub fn lrcoef(outer: &[i32], inner1: &[i32], inner2: &[i32]) -> Result<u128, LrCoefError> {
    match optim_coef(outer, inner1, inner2)? {
        OptimizedCoef::Zero => Ok(0),
        OptimizedCoef::One => Ok(1),
        OptimizedCoef::Count(shape) => lrcoef_count(&shape.outer, &shape.inner, &shape.content),
    }
}

/// Compute a single coefficient and downcast to the upstream ABI return type.
pub fn lrcoef_i64(outer: &[i32], inner1: &[i32], inner2: &[i32]) -> Option<i64> {
    let coef = lrcoef(outer, inner1, inner2).ok()?;
    i64::try_from(coef).ok()
}

/// Count relative interior lattice points using Buch's tableau search order.
pub fn lrcoef_buch_interior_u128(
    outer: &[i32],
    inner: &[i32],
    content: &[i32],
) -> Result<u128, LrCoefError> {
    Ok(lrcoef_buch_interior_stats(outer, inner, content)?.value)
}

/// Count relative interior lattice points with Buch-search diagnostics.
pub fn lrcoef_buch_interior_stats(
    outer: &[i32],
    inner: &[i32],
    content: &[i32],
) -> Result<LrBuchInteriorStats, LrCoefError> {
    match optim_coef(outer, inner, content)? {
        OptimizedCoef::Zero => Ok(zero_buch_interior_stats()),
        OptimizedCoef::One => Ok(one_buch_interior_stats()),
        OptimizedCoef::Count(shape) => {
            lrcoef_buch_interior_stats_compacted(&shape.outer, &shape.inner, &shape.content)
        }
    }
}

/// Count relative interior lattice points using Buch tight facets and a
/// memoized GT-chain suffix counter.
pub fn lrcoef_buch_interior_memo_u128(
    outer: &[i32],
    inner: &[i32],
    content: &[i32],
) -> Result<u128, LrCoefError> {
    Ok(lrcoef_buch_interior_memo_stats(outer, inner, content)?.value)
}

/// Count relative interior lattice points with memoized suffix diagnostics.
pub fn lrcoef_buch_interior_memo_stats(
    outer: &[i32],
    inner: &[i32],
    content: &[i32],
) -> Result<LrBuchInteriorMemoStats, LrCoefError> {
    match optim_coef(outer, inner, content)? {
        OptimizedCoef::Zero => Ok(zero_buch_interior_memo_stats()),
        OptimizedCoef::One => Ok(one_buch_interior_memo_stats()),
        OptimizedCoef::Count(shape) => {
            lrcoef_buch_interior_memo_stats_compacted(&shape.outer, &shape.inner, &shape.content)
        }
    }
}

/// Count both the full LR polytope lattice points and its relative interior.
pub fn lrcoef_buch_counts_u128(
    outer: &[i32],
    inner: &[i32],
    content: &[i32],
) -> Result<LrBuchCounts, LrCoefError> {
    let stats = lrcoef_buch_interior_stats(outer, inner, content)?;
    Ok(LrBuchCounts {
        full: stats.weak_tableaux,
        interior: stats.strict_tableaux,
    })
}

/// Dimension of the LR GT/Yamanouchi polytope using Buch's tableau search.
pub fn lrcoef_buch_dimension(
    outer: &[i32],
    inner: &[i32],
    content: &[i32],
) -> Result<Option<usize>, LrCoefError> {
    match optim_coef(outer, inner, content)? {
        OptimizedCoef::Zero => Ok(None),
        OptimizedCoef::One => Ok(Some(0)),
        OptimizedCoef::Count(shape) => {
            lrcoef_buch_dimension_compacted(&shape.outer, &shape.inner, &shape.content)
        }
    }
}

pub(crate) fn lrcoef_buch_stretch_cache(
    outer: &[i32],
    inner: &[i32],
    content: &[i32],
) -> Result<Option<LrBuchStretchCache>, LrCoefError> {
    match optim_coef(outer, inner, content)? {
        OptimizedCoef::Zero => Ok(None),
        OptimizedCoef::One => Ok(Some(LrBuchStretchCache {
            outer: Vec::new(),
            inner: Vec::new(),
            content: Vec::new(),
            tight: BuchTightFlags::new(0, 0),
            dimension: 0,
        })),
        OptimizedCoef::Count(shape) => {
            let Some(tight_data) = buch_tight_data(&shape.outer, &shape.inner, &shape.content)?
            else {
                return Ok(None);
            };
            let dimension = buch_dimension_from_tight_flags(&tight_data.tight);
            Ok(Some(LrBuchStretchCache {
                outer: shape.outer,
                inner: shape.inner,
                content: shape.content,
                tight: tight_data.tight,
                dimension,
            }))
        }
    }
}

pub(crate) fn lrcoef_buch_stretched_counts_u128(
    cache: &LrBuchStretchCache,
    stretch: u64,
) -> Result<LrBuchCounts, LrCoefError> {
    if cache.dimension == 0 || stretch == 0 {
        return Ok(LrBuchCounts {
            full: 1,
            interior: 1,
        });
    }
    let outer = scale_partition_i32(&cache.outer, stretch)?;
    let inner = scale_partition_i32(&cache.inner, stretch)?;
    let content = scale_partition_i32(&cache.content, stretch)?;
    Ok(LrBuchCounts {
        full: lrcoef_count(&outer, &inner, &content)?,
        interior: lrcoef_count_strict_tableaux(&outer, &inner, &content, &cache.tight)?,
    })
}

impl LrBuchStretchCache {
    pub(crate) fn dimension(&self) -> usize {
        self.dimension
    }
}

fn valid_partition(part: &[i32]) -> bool {
    part.iter().all(|&x| x >= 0) && part.windows(2).all(|window| window[0] >= window[1])
}

fn part_length(part: &[i32]) -> usize {
    part.iter()
        .rposition(|&x| x != 0)
        .map_or(0, |index| index + 1)
}

fn part_entry(part: &[i32], index: usize) -> i32 {
    part.get(index).copied().unwrap_or(0)
}

fn part_sum(part: &[i32]) -> Result<i32, LrCoefError> {
    let mut sum = 0i32;
    for &part in part {
        sum = sum
            .checked_add(part)
            .ok_or(LrCoefError::ArithmeticOverflow)?;
    }
    Ok(sum)
}

fn scale_partition_i32(parts: &[i32], stretch: u64) -> Result<Vec<i32>, LrCoefError> {
    parts
        .iter()
        .copied()
        .map(|part| {
            if part < 0 {
                return Err(LrCoefError::InvalidPartition);
            }
            let value = (part as u64)
                .checked_mul(stretch)
                .ok_or(LrCoefError::ArithmeticOverflow)?;
            i32::try_from(value).map_err(|_| LrCoefError::ArithmeticOverflow)
        })
        .collect()
}

fn partition_to_u32_padded(parts: &[i32], len: usize) -> Result<Vec<u32>, LrCoefError> {
    if parts.len() > len && parts[len..].iter().any(|&part| part != 0) {
        return Err(LrCoefError::InvalidPartition);
    }
    let mut result = vec![0u32; len];
    for (index, &part) in parts.iter().take(len).enumerate() {
        if part < 0 {
            return Err(LrCoefError::InvalidPartition);
        }
        result[index] = part as u32;
    }
    Ok(result)
}

fn has_positive_part_beyond(part: &[i32], index: usize) -> bool {
    part.get(index).is_some_and(|&part| part > 0)
}

pub(crate) fn optim_coef(
    outer: &[i32],
    inner1: &[i32],
    inner2: &[i32],
) -> Result<OptimizedCoef, LrCoefError> {
    if !valid_partition(outer) || !valid_partition(inner1) || !valid_partition(inner2) {
        return Err(LrCoefError::InvalidPartition);
    }

    let mut n = part_length(outer);
    if has_positive_part_beyond(inner1, n) || has_positive_part_beyond(inner2, n) {
        return Ok(OptimizedCoef::Zero);
    }
    if n == 0 {
        return Ok(OptimizedCoef::One);
    }

    let mut nu = vec![0; n];
    let mut la = vec![0; n];
    let mut mu = vec![0; n];

    let mut sum = 0i32;
    for r in (0..n).rev() {
        nu[r] = outer[r];
        sum = sum
            .checked_add(nu[r])
            .ok_or(LrCoefError::ArithmeticOverflow)?;
    }

    let mut nla = n;
    while nla > inner1.len() {
        nla -= 1;
        la[nla] = 0;
    }
    while nla > 0 && inner1[nla - 1] == 0 {
        nla -= 1;
        la[nla] = 0;
    }
    for r in (0..nla).rev() {
        let x = inner1[r];
        la[r] = x;
        if nu[r] < x {
            return Ok(OptimizedCoef::Zero);
        }
        sum = sum
            .checked_sub(la[r])
            .ok_or(LrCoefError::ArithmeticOverflow)?;
    }

    let mut nmu = n;
    while nmu > inner2.len() {
        nmu -= 1;
        mu[nmu] = 0;
    }
    while nmu > 0 && inner2[nmu - 1] == 0 {
        nmu -= 1;
        mu[nmu] = 0;
    }
    for r in (0..nmu).rev() {
        let x = inner2[r];
        mu[r] = x;
        if nu[r] < x {
            return Ok(OptimizedCoef::Zero);
        }
        sum = sum
            .checked_sub(mu[r])
            .ok_or(LrCoefError::ArithmeticOverflow)?;
    }

    if sum != 0 {
        return Ok(OptimizedCoef::Zero);
    }

    let mut n0 = n + 1;
    let mut nu0 = 0;
    while n < n0 || nu[0] < nu0 {
        n0 = n;
        nu0 = nu[0];

        /* Horizontal compactification of nu/la. */
        let mu0 = mu[0];
        let mut lar1 = 0;
        let mut nur1 = 0;
        let mut r = n;
        while r > 0 {
            let i = r - 1;
            let lar = la[i];
            let nur = nu[i];
            if lar > nur1 || nur - lar1 > mu0 {
                break;
            }
            lar1 = lar;
            nur1 = nur;
            r -= 1;
        }
        let mut c = 0;
        while r > 0 {
            let i = r - 1;
            let lar = la[i];
            let nur = nu[i];
            if nur - lar > mu0 {
                return Ok(OptimizedCoef::Zero);
            }
            let mut ca = nur - lar1 - mu0;
            if ca < lar - nur1 {
                ca = lar - nur1;
            }
            if ca > 0 {
                c += ca;
            }
            if nur - c < mu[i] {
                return Ok(OptimizedCoef::Zero);
            }
            if nur == c {
                n = i;
                break;
            }
            la[i] = lar - c;
            nu[i] = nur - c;
            lar1 = lar;
            nur1 = nur;
            r -= 1;
        }

        /* Remove row of size mu[0] from nu/la. */
        let mu0 = mu[0];
        let mut r = 0;
        while r < n && nu[r] - la[r] < mu0 {
            r += 1;
        }
        if r < n {
            if nu[r] - la[r] > mu0 {
                return Ok(OptimizedCoef::Zero);
            }
            while r + 1 < n {
                la[r] = la[r + 1];
                nu[r] = nu[r + 1];
                r += 1;
            }
            for r in 0..n - 1 {
                mu[r] = mu[r + 1];
            }
            n -= 1;
        }

        /* Horizontal compactification of nu/mu. */
        let la0 = la[0];
        let mut mur1 = 0;
        let mut nur1 = 0;
        let mut r = n;
        while r > 0 {
            let i = r - 1;
            let mur = mu[i];
            let nur = nu[i];
            if mur > nur1 || nur - mur1 > la0 {
                break;
            }
            mur1 = mur;
            nur1 = nur;
            r -= 1;
        }
        let mut c = 0;
        while r > 0 {
            let i = r - 1;
            let mur = mu[i];
            let nur = nu[i];
            if nur - mur > la0 {
                return Ok(OptimizedCoef::Zero);
            }
            let mut ca = nur - mur1 - la0;
            if ca < mur - nur1 {
                ca = mur - nur1;
            }
            if ca > 0 {
                c += ca;
            }
            if nur - c < la[i] {
                return Ok(OptimizedCoef::Zero);
            }
            if nur == c {
                n = i;
                break;
            }
            mu[i] = mur - c;
            nu[i] = nur - c;
            mur1 = mur;
            nur1 = nur;
            r -= 1;
        }

        /* Remove row of size la[0] from nu/mu. */
        let la0 = la[0];
        let mut r = 0;
        while r < n && nu[r] - mu[r] < la0 {
            r += 1;
        }
        if r < n {
            if nu[r] - mu[r] > la0 {
                return Ok(OptimizedCoef::Zero);
            }
            while r + 1 < n {
                mu[r] = mu[r + 1];
                nu[r] = nu[r + 1];
                r += 1;
            }
            for r in 0..n - 1 {
                la[r] = la[r + 1];
            }
            n -= 1;
        }

        /* Vertical compactification of nu/la. */
        if n < nmu {
            nmu = n;
        }
        while nmu > 0 && mu[nmu - 1] == 0 {
            nmu -= 1;
        }
        if nmu == 0 {
            return Ok(OptimizedCoef::One);
        }
        let mut r = 0;
        while r < nmu && la[r] < nu[r] {
            r += 1;
        }
        while r < n && la[r] < nu[r] && nu[r] < la[r - nmu] {
            r += 1;
        }
        if r < n {
            let mut inu = r;
            let mut s = r.saturating_sub(nmu);
            let mut ila = s;
            while r < n && inu < nmu {
                if la[r] == nu[r] {
                    la[r] = -1;
                } else {
                    nu[inu] = nu[r];
                    if nu[inu] < mu[inu] {
                        return Ok(OptimizedCoef::Zero);
                    }
                    inu += 1;
                }
                r += 1;
            }
            while r < n {
                if la[r] == nu[r] {
                    la[r] = -1;
                    r += 1;
                    continue;
                }
                while la[s] == -1 {
                    s += 1;
                }
                if la[s] < nu[r] {
                    return Ok(OptimizedCoef::Zero);
                }
                if la[s] > nu[r] {
                    la[ila] = la[s];
                    ila += 1;
                    nu[inu] = nu[r];
                    if nu[inu] < mu[inu] {
                        return Ok(OptimizedCoef::Zero);
                    }
                    inu += 1;
                }
                r += 1;
                s += 1;
            }
            while s < n {
                if la[s] != -1 {
                    la[ila] = la[s];
                    ila += 1;
                }
                s += 1;
            }
            if inu < n && mu[inu] > 0 {
                return Ok(OptimizedCoef::Zero);
            }
            n = inu;
        }

        /* Remove column of size len(mu) from nu/la. */
        let mut r = nmu;
        while r <= n && nu[r - 1] <= la[r - nmu] {
            r += 1;
        }
        if r <= n {
            if r > nmu && nu[r - 1] > la[r - nmu - 1] {
                return Ok(OptimizedCoef::Zero);
            }
            if r < n && nu[r] > la[r - nmu] {
                return Ok(OptimizedCoef::Zero);
            }
            if r > nmu {
                let mut s = r - nmu;
                while s > 0 {
                    s -= 1;
                    la[s] -= 1;
                }
            }
            for s in (0..nmu).rev() {
                mu[s] -= 1;
            }
            for (s, entry) in nu.iter_mut().enumerate().take(r) {
                *entry -= 1;
                if *entry == 0 {
                    n = s;
                    break;
                }
            }
        }

        /* Vertical compactification of nu/mu. */
        if n < nla {
            nla = n;
        }
        while nla > 0 && la[nla - 1] == 0 {
            nla -= 1;
        }
        if nla == 0 {
            return Ok(OptimizedCoef::One);
        }
        let mut r = 0;
        while r < nla && mu[r] < nu[r] {
            r += 1;
        }
        while r < n && mu[r] < nu[r] && nu[r] < mu[r - nla] {
            r += 1;
        }
        if r < n {
            let mut inu = r;
            let mut s = r.saturating_sub(nla);
            let mut imu = s;
            while r < n && inu < nla {
                if mu[r] == nu[r] {
                    mu[r] = -1;
                } else {
                    nu[inu] = nu[r];
                    if nu[inu] < la[inu] {
                        return Ok(OptimizedCoef::Zero);
                    }
                    inu += 1;
                }
                r += 1;
            }
            while r < n {
                if mu[r] == nu[r] {
                    mu[r] = -1;
                    r += 1;
                    continue;
                }
                while mu[s] == -1 {
                    s += 1;
                }
                if mu[s] < nu[r] {
                    return Ok(OptimizedCoef::Zero);
                }
                if mu[s] > nu[r] {
                    mu[imu] = mu[s];
                    imu += 1;
                    nu[inu] = nu[r];
                    if nu[inu] < la[inu] {
                        return Ok(OptimizedCoef::Zero);
                    }
                    inu += 1;
                }
                r += 1;
                s += 1;
            }
            while s < n {
                if mu[s] != -1 {
                    mu[imu] = mu[s];
                    imu += 1;
                }
                s += 1;
            }
            if inu < n && la[inu] > 0 {
                return Ok(OptimizedCoef::Zero);
            }
            n = inu;
        }

        /* Remove column of size len(la) from nu/mu. */
        let mut r = nla;
        while r <= n && nu[r - 1] <= mu[r - nla] {
            r += 1;
        }
        if r <= n {
            if r > nla && nu[r - 1] > mu[r - nla - 1] {
                return Ok(OptimizedCoef::Zero);
            }
            if r < n && nu[r] > mu[r - nla] {
                return Ok(OptimizedCoef::Zero);
            }
            if r > nla {
                let mut s = r - nla;
                while s > 0 {
                    s -= 1;
                    mu[s] -= 1;
                }
            }
            for s in (0..nla).rev() {
                la[s] -= 1;
            }
            for (s, entry) in nu.iter_mut().enumerate().take(r) {
                *entry -= 1;
                if *entry == 0 {
                    n = s;
                    break;
                }
            }
        }
    }

    if n == 0 {
        return Ok(OptimizedCoef::One);
    }

    nu.truncate(n);
    if n < nla {
        nla = n;
    }
    while nla > 0 && la[nla - 1] == 0 {
        nla -= 1;
    }
    la.truncate(nla);
    if n < nmu {
        nmu = n;
    }
    while nmu > 0 && mu[nmu - 1] == 0 {
        nmu -= 1;
    }
    mu.truncate(nmu);

    Ok(OptimizedCoef::Count(SkewShape {
        outer: nu,
        inner: la,
        content: mu,
    }))
}

fn lrcoef_count(outer: &[i32], inner: &[i32], content: &[i32]) -> Result<u128, LrCoefError> {
    let outer_sum = part_sum(outer)?;
    let inner_sum = part_sum(inner)?;
    let content_sum = part_sum(content)?;
    let expected_outer_sum = inner_sum
        .checked_add(content_sum)
        .ok_or(LrCoefError::ArithmeticOverflow)?;
    debug_assert_eq!(outer_sum, expected_outer_sum);
    debug_assert!(content_sum > 1);

    let mut boxes = new_skewtab(outer, inner, part_length(content), content_sum)?;
    let mut counts = new_content(content);

    let n = content_sum;
    let mut pos = 0usize;
    let mut above = boxes[boxes[pos].north].value;
    let mut x = 1i32;
    let mut se_supply = n - counts[1].supply;
    let mut coef = 0u128;

    loop {
        while x > above
            && (counts[x as usize].cont == counts[x as usize].supply
                || counts[x as usize].cont == counts[(x - 1) as usize].cont)
        {
            se_supply += counts[x as usize].supply - counts[x as usize].cont;
            x -= 1;
        }

        if x == above || n - pos as i32 - se_supply <= boxes[pos].west_sz {
            if pos == 0 {
                break;
            }
            pos -= 1;
            se_supply = boxes[pos].se_supply;
            above = boxes[boxes[pos].north].value;
            x = boxes[pos].value;
            counts[x as usize].cont -= 1;
            se_supply += counts[x as usize].supply - counts[x as usize].cont;
            x -= 1;
        } else if pos + 1 < n as usize {
            boxes[pos].se_supply = se_supply;
            boxes[pos].value = x;
            counts[x as usize].cont += 1;
            pos += 1;
            se_supply = boxes[boxes[pos].east].se_supply;
            x = boxes[boxes[pos].east].value;
            above = boxes[boxes[pos].north].value;
            while x > boxes[pos].max {
                se_supply += counts[x as usize].supply - counts[x as usize].cont;
                x -= 1;
            }
            while x > above && se_supply < boxes[pos].se_sz {
                se_supply += counts[x as usize].supply - counts[x as usize].cont;
                x -= 1;
            }
        } else {
            coef = coef.checked_add(1).ok_or(LrCoefError::ArithmeticOverflow)?;
            pos -= 1;
            se_supply = boxes[pos].se_supply;
            above = boxes[boxes[pos].north].value;
            x = boxes[pos].value;
            counts[x as usize].cont -= 1;
            se_supply += counts[x as usize].supply - counts[x as usize].cont;
            x -= 1;
        }
    }

    Ok(coef)
}

fn lrcoef_buch_interior_stats_compacted(
    outer: &[i32],
    inner: &[i32],
    content: &[i32],
) -> Result<LrBuchInteriorStats, LrCoefError> {
    let rows = outer.len();
    if rows == 0 {
        return Ok(if part_sum(content)? == 0 {
            one_buch_interior_stats()
        } else {
            zero_buch_interior_stats()
        });
    }

    let Some(tight_data) = buch_tight_data(outer, inner, content)? else {
        return Ok(zero_buch_interior_stats());
    };
    let tight = tight_data.tight;
    let weak_tableaux = tight_data.weak_tableaux;

    let (strict_lower, strict_diagonal, strict_yamanouchi) = tight.strict_counts();
    let strict_tableaux = lrcoef_count_strict_tableaux(outer, inner, content, &tight)?;

    Ok(LrBuchInteriorStats {
        value: strict_tableaux,
        weak_tableaux,
        strict_tableaux,
        strict_lower_constraints: strict_lower,
        strict_diagonal_constraints: strict_diagonal,
        strict_yamanouchi_constraints: strict_yamanouchi,
    })
}

fn lrcoef_buch_interior_memo_stats_compacted(
    outer: &[i32],
    inner: &[i32],
    content: &[i32],
) -> Result<LrBuchInteriorMemoStats, LrCoefError> {
    let rows = outer.len();
    if rows == 0 {
        return Ok(if part_sum(content)? == 0 {
            one_buch_interior_memo_stats()
        } else {
            zero_buch_interior_memo_stats()
        });
    }

    let Some(tight_data) = buch_tight_data(outer, inner, content)? else {
        return Ok(zero_buch_interior_memo_stats());
    };
    let tight = tight_data.tight;
    let weak_tableaux = tight_data.weak_tableaux;

    let (strict_lower, strict_diagonal, strict_yamanouchi) = tight.strict_counts();
    let (value, memo_states, counters) =
        lrcoef_count_strict_tableaux_memo(outer, inner, content, &tight)?;

    Ok(LrBuchInteriorMemoStats {
        value,
        weak_tableaux,
        memo_states,
        cache_hits: counters.cache_hits,
        generated_transitions: counters.generated_transitions,
        strict_lower_constraints: strict_lower,
        strict_diagonal_constraints: strict_diagonal,
        strict_yamanouchi_constraints: strict_yamanouchi,
    })
}

fn lrcoef_buch_dimension_compacted(
    outer: &[i32],
    inner: &[i32],
    content: &[i32],
) -> Result<Option<usize>, LrCoefError> {
    if outer.is_empty() {
        return Ok((part_sum(content)? == 0).then_some(0));
    }
    let Some(tight_data) = buch_tight_data(outer, inner, content)? else {
        return Ok(None);
    };
    Ok(Some(buch_dimension_from_tight_flags(&tight_data.tight)))
}

fn buch_tight_data(
    outer: &[i32],
    inner: &[i32],
    content: &[i32],
) -> Result<Option<BuchTightData>, LrCoefError> {
    let rows = outer.len();
    let steps = part_length(content);
    let mut tight = BuchTightFlags::new(steps, rows);
    let mut weak_tableaux = 0u128;
    lrcoef_for_each_tableau(outer, inner, content, |boxes| {
        weak_tableaux = weak_tableaux
            .checked_add(1)
            .ok_or(LrCoefError::ArithmeticOverflow)?;
        let increments = tableau_increments(boxes, rows, steps)?;
        tight.observe_tableau(inner, &increments);
        Ok(())
    })?;

    if weak_tableaux == 0 {
        Ok(None)
    } else {
        Ok(Some(BuchTightData {
            tight,
            weak_tableaux,
        }))
    }
}

fn lrcoef_for_each_tableau<F>(
    outer: &[i32],
    inner: &[i32],
    content: &[i32],
    mut visit: F,
) -> Result<(), LrCoefError>
where
    F: FnMut(&[LrCoefBox]) -> Result<(), LrCoefError>,
{
    let inner_sum = part_sum(inner)?;
    let content_sum = part_sum(content)?;
    debug_assert_eq!(part_sum(outer)?, inner_sum + content_sum);

    if content_sum == 0 {
        return Ok(());
    }

    let mut boxes = new_skewtab(outer, inner, part_length(content), content_sum)?;
    let mut counts = new_content(content);

    let n = content_sum;
    let real_boxes = usize::try_from(n).map_err(|_| LrCoefError::ArithmeticOverflow)?;
    let mut pos = 0usize;
    let mut above = boxes[boxes[pos].north].value;
    let mut x = 1i32;
    let mut se_supply = n - counts[1].supply;

    loop {
        while x > above
            && (counts[x as usize].cont == counts[x as usize].supply
                || counts[x as usize].cont == counts[(x - 1) as usize].cont)
        {
            se_supply += counts[x as usize].supply - counts[x as usize].cont;
            x -= 1;
        }

        if x == above || n - pos as i32 - se_supply <= boxes[pos].west_sz {
            if pos == 0 {
                break;
            }
            pos -= 1;
            se_supply = boxes[pos].se_supply;
            above = boxes[boxes[pos].north].value;
            x = boxes[pos].value;
            counts[x as usize].cont -= 1;
            se_supply += counts[x as usize].supply - counts[x as usize].cont;
            x -= 1;
        } else if pos + 1 < real_boxes {
            boxes[pos].se_supply = se_supply;
            boxes[pos].value = x;
            counts[x as usize].cont += 1;
            pos += 1;
            se_supply = boxes[boxes[pos].east].se_supply;
            x = boxes[boxes[pos].east].value;
            above = boxes[boxes[pos].north].value;
            while x > boxes[pos].max {
                se_supply += counts[x as usize].supply - counts[x as usize].cont;
                x -= 1;
            }
            while x > above && se_supply < boxes[pos].se_sz {
                se_supply += counts[x as usize].supply - counts[x as usize].cont;
                x -= 1;
            }
        } else {
            boxes[pos].value = x;
            counts[x as usize].cont += 1;
            visit(&boxes[..real_boxes])?;
            counts[x as usize].cont -= 1;

            if pos == 0 {
                break;
            }
            pos -= 1;
            se_supply = boxes[pos].se_supply;
            above = boxes[boxes[pos].north].value;
            x = boxes[pos].value;
            counts[x as usize].cont -= 1;
            se_supply += counts[x as usize].supply - counts[x as usize].cont;
            x -= 1;
        }
    }

    Ok(())
}

fn lrcoef_count_strict_tableaux(
    outer: &[i32],
    inner: &[i32],
    content: &[i32],
    tight: &BuchTightFlags,
) -> Result<u128, LrCoefError> {
    let inner_sum = part_sum(inner)?;
    let content_sum = part_sum(content)?;
    debug_assert_eq!(part_sum(outer)?, inner_sum + content_sum);

    if content_sum == 0 {
        return Ok(1);
    }

    let mut boxes = new_skewtab(outer, inner, part_length(content), content_sum)?;
    let mut counts = new_content(content);
    let mut pruner = BuchInteriorPruner::new(tight, inner)?;

    let n = content_sum;
    let real_boxes = usize::try_from(n).map_err(|_| LrCoefError::ArithmeticOverflow)?;
    let mut pos = 0usize;
    let mut above = boxes[boxes[pos].north].value;
    let mut x = 1i32;
    let mut se_supply = n - counts[1].supply;
    let mut coef = 0u128;

    loop {
        while x > above
            && (counts[x as usize].cont == counts[x as usize].supply
                || counts[x as usize].cont == counts[(x - 1) as usize].cont)
        {
            se_supply += counts[x as usize].supply - counts[x as usize].cont;
            x -= 1;
        }

        if x == above || n - pos as i32 - se_supply <= boxes[pos].west_sz {
            if pos == 0 {
                break;
            }
            pos -= 1;
            pruner.unplace(&boxes[pos])?;
            se_supply = boxes[pos].se_supply;
            above = boxes[boxes[pos].north].value;
            x = boxes[pos].value;
            counts[x as usize].cont -= 1;
            se_supply += counts[x as usize].supply - counts[x as usize].cont;
            x -= 1;
        } else if pos + 1 < real_boxes {
            boxes[pos].se_supply = se_supply;
            boxes[pos].value = x;
            counts[x as usize].cont += 1;
            let row_complete = boxes[pos + 1].row != boxes[pos].row;
            if pruner.place(&boxes[pos], row_complete)? {
                pos += 1;
                se_supply = boxes[boxes[pos].east].se_supply;
                x = boxes[boxes[pos].east].value;
                above = boxes[boxes[pos].north].value;
                while x > boxes[pos].max {
                    se_supply += counts[x as usize].supply - counts[x as usize].cont;
                    x -= 1;
                }
                while x > above && se_supply < boxes[pos].se_sz {
                    se_supply += counts[x as usize].supply - counts[x as usize].cont;
                    x -= 1;
                }
            } else {
                pruner.unplace(&boxes[pos])?;
                counts[x as usize].cont -= 1;
                se_supply += counts[x as usize].supply - counts[x as usize].cont;
                x -= 1;
            }
        } else {
            boxes[pos].value = x;
            counts[x as usize].cont += 1;
            if pruner.place(&boxes[pos], true)? && pruner.is_complete_relative_interior() {
                coef = coef.checked_add(1).ok_or(LrCoefError::ArithmeticOverflow)?;
            }
            pruner.unplace(&boxes[pos])?;
            counts[x as usize].cont -= 1;

            if pos == 0 {
                break;
            }
            pos -= 1;
            pruner.unplace(&boxes[pos])?;
            se_supply = boxes[pos].se_supply;
            above = boxes[boxes[pos].north].value;
            x = boxes[pos].value;
            counts[x as usize].cont -= 1;
            se_supply += counts[x as usize].supply - counts[x as usize].cont;
            x -= 1;
        }
    }

    Ok(coef)
}

fn lrcoef_count_strict_tableaux_memo(
    outer: &[i32],
    inner: &[i32],
    content: &[i32],
    tight: &BuchTightFlags,
) -> Result<(u128, usize, BuchMemoCounters), LrCoefError> {
    let outer_sum = part_sum(outer)?;
    let inner_sum = part_sum(inner)?;
    let content_sum = part_sum(content)?;
    debug_assert_eq!(outer_sum, inner_sum + content_sum);

    if content_sum == 0 {
        return Ok((1, 0, BuchMemoCounters::default()));
    }

    let rows = outer.len();
    let steps = part_length(content);
    let outer = partition_to_u32_padded(outer, rows)?;
    let inner = partition_to_u32_padded(inner, rows)?;
    let content = partition_to_u32_padded(content, steps)?;
    let previous_prefix = vec![0u32; rows + 1];
    let mut counter = BuchMemoCounter::new(outer, content, tight);
    let value = counter.count(0, &inner, &previous_prefix)?;
    let memo_states = counter.cache.len();
    Ok((value, memo_states, counter.counters))
}

fn tableau_increments(
    boxes: &[LrCoefBox],
    rows: usize,
    steps: usize,
) -> Result<Vec<u32>, LrCoefError> {
    let mut increments = vec![0u32; rows.saturating_mul(steps)];
    for cell in boxes {
        let step = usize::try_from(cell.value - 1).map_err(|_| LrCoefError::ArithmeticOverflow)?;
        if step >= steps || cell.row >= rows {
            return Err(LrCoefError::ArithmeticOverflow);
        }
        let index = step
            .checked_mul(rows)
            .and_then(|value| value.checked_add(cell.row))
            .ok_or(LrCoefError::ArithmeticOverflow)?;
        increments[index] = increments[index]
            .checked_add(1)
            .ok_or(LrCoefError::ArithmeticOverflow)?;
    }
    Ok(increments)
}

fn new_content(content: &[i32]) -> Vec<LrCoefContent> {
    let n = part_length(content);
    debug_assert!(n > 0);
    let mut result = vec![LrCoefContent::default(); n + 1];
    result[0].cont = content[0];
    result[0].supply = content[0];
    for i in 0..n {
        result[i + 1].supply = content[i];
    }
    result
}

fn new_skewtab(
    outer: &[i32],
    inner: &[i32],
    max_value: usize,
    skew_size: i32,
) -> Result<Vec<LrCoefBox>, LrCoefError> {
    debug_assert!(valid_partition(outer));
    debug_assert!(valid_partition(inner));
    debug_assert!(outer[0] > 0);

    let n = usize::try_from(skew_size).map_err(|_| LrCoefError::ArithmeticOverflow)?;
    let max_value = i32::try_from(max_value).map_err(|_| LrCoefError::ArithmeticOverflow)?;
    let mut array = vec![LrCoefBox::default(); n + 2];
    let mut pos = n;

    for r in (0..outer.len()).rev() {
        let nu_0 = if r == 0 { outer[0] } else { outer[r - 1] };
        let la_0 = if r == 0 {
            outer[0]
        } else {
            part_entry(inner, r - 1)
        };
        let nu_r = outer[r];
        let la_r = part_entry(inner, r);
        let nu_1 = part_entry(outer, r + 1);
        for c in la_r..nu_r {
            pos -= 1;
            let north = if la_0 <= c && c < nu_0 {
                usize::try_from(pos as i32 - nu_r + la_0)
                    .map_err(|_| LrCoefError::ArithmeticOverflow)?
            } else {
                n
            };
            let east = if c + 1 < nu_r { pos - 1 } else { n + 1 };
            let west_sz = c - la_r;
            let (max, se_sz) = if c >= nu_1 {
                (max_value, 0)
            } else {
                let below = usize::try_from(pos as i32 + nu_1 - la_r)
                    .map_err(|_| LrCoefError::ArithmeticOverflow)?;
                (array[below].max - 1, array[below].se_sz + nu_1 - c)
            };
            array[pos] = LrCoefBox {
                value: 0,
                max,
                row: r,
                north,
                east,
                se_supply: 0,
                se_sz,
                west_sz,
            };
        }
    }

    debug_assert_eq!(pos, 0);
    array[n].value = 0;
    array[n + 1].value = max_value;
    array[n + 1].se_supply = 0;
    Ok(array)
}

impl BuchTightFlags {
    fn new(steps: usize, rows: usize) -> Self {
        Self {
            lower: vec![vec![true; rows]; steps],
            diagonal: vec![vec![true; rows]; steps],
            yamanouchi: vec![vec![true; rows]; steps],
        }
    }

    fn observe_tableau(&mut self, inner: &[i32], increments: &[u32]) {
        let rows = self.lower.first().map_or(0, Vec::len);
        let mut shape = (0..rows)
            .map(|row| part_entry(inner, row) as u32)
            .collect::<Vec<_>>();

        for step in 0..self.lower.len() {
            for row in 0..rows {
                let increment = increment_at(increments, rows, step, row);
                if increment > 0 {
                    self.lower[step][row] = false;
                }
                if row > 0 && shape[row] + increment < shape[row - 1] {
                    self.diagonal[step][row] = false;
                }
            }
            for (row, shape_row) in shape.iter_mut().enumerate().take(rows) {
                *shape_row += increment_at(increments, rows, step, row);
            }
            if step > 0 {
                for row in 0..rows {
                    let current_prefix = prefix_sum(increments, rows, step, row + 1);
                    let previous_prefix = prefix_sum(increments, rows, step - 1, row);
                    if current_prefix < previous_prefix {
                        self.yamanouchi[step][row] = false;
                    }
                }
            }
        }
    }

    fn strict_counts(&self) -> (usize, usize, usize) {
        let lower = self
            .lower
            .iter()
            .flat_map(|step| step.iter())
            .filter(|&&is_tight| !is_tight)
            .count();
        let diagonal = self
            .diagonal
            .iter()
            .flat_map(|step| {
                step.iter()
                    .enumerate()
                    .skip(1)
                    .map(|(_, is_tight)| is_tight)
            })
            .filter(|&&is_tight| !is_tight)
            .count();
        let yamanouchi = self
            .yamanouchi
            .iter()
            .enumerate()
            .skip(1)
            .flat_map(|(_, step)| step.iter())
            .filter(|&&is_tight| !is_tight)
            .count();
        (lower, diagonal, yamanouchi)
    }
}

impl<'a> BuchMemoCounter<'a> {
    fn new(outer: Vec<u32>, content: Vec<u32>, tight: &'a BuchTightFlags) -> Self {
        Self {
            outer,
            content,
            tight,
            cache: HashMap::new(),
            counters: BuchMemoCounters::default(),
        }
    }

    fn count(
        &mut self,
        step: usize,
        shape: &[u32],
        previous_prefix: &[u32],
    ) -> Result<u128, LrCoefError> {
        let state = BuchMemoState {
            step,
            shape: shape.to_vec(),
            previous_prefix: previous_prefix.to_vec(),
        };
        if let Some(&value) = self.cache.get(&state) {
            self.counters.cache_hits = self
                .counters
                .cache_hits
                .checked_add(1)
                .ok_or(LrCoefError::ArithmeticOverflow)?;
            return Ok(value);
        }

        let value = if step == self.content.len() {
            u128::from(shape == self.outer.as_slice())
        } else {
            let rows = self.outer.len();
            let mut next_shape = shape.to_vec();
            let mut next_prefix = vec![0u32; rows + 1];
            self.enumerate_step(
                step,
                0,
                self.content[step],
                shape,
                previous_prefix,
                &mut next_shape,
                &mut next_prefix,
            )?
        };
        self.cache.insert(state, value);
        Ok(value)
    }

    #[allow(clippy::too_many_arguments)]
    fn enumerate_step(
        &mut self,
        step: usize,
        row: usize,
        remaining: u32,
        shape: &[u32],
        previous_prefix: &[u32],
        next_shape: &mut [u32],
        next_prefix: &mut [u32],
    ) -> Result<u128, LrCoefError> {
        if row == self.outer.len() {
            if remaining == 0 {
                return self.count(step + 1, next_shape, next_prefix);
            }
            return Ok(0);
        }

        let lower_min = if self.tight.lower[step][row] { 0 } else { 1 };
        if remaining < lower_min {
            return Ok(0);
        }

        let base = shape[row];
        let max_from_outer = self.outer[row].saturating_sub(base);
        let max_from_horizontal_strip = if row == 0 {
            remaining
        } else {
            let weak_bound = shape[row - 1].saturating_sub(base);
            if self.tight.diagonal[step][row] {
                weak_bound
            } else {
                weak_bound.saturating_sub(1)
            }
        };
        let max_from_yamanouchi = if step == 0 {
            remaining
        } else {
            let prefix_sum = next_prefix[row];
            let weak_bound = previous_prefix[row].saturating_sub(prefix_sum);
            if self.tight.yamanouchi[step][row] {
                weak_bound
            } else {
                weak_bound.saturating_sub(1)
            }
        };
        let max_increment = remaining
            .min(max_from_outer)
            .min(max_from_horizontal_strip)
            .min(max_from_yamanouchi);
        if max_increment < lower_min {
            return Ok(0);
        }

        let mut total = 0u128;
        let prefix_sum = next_prefix[row];
        for increment in lower_min..=max_increment {
            self.counters.generated_transitions = self
                .counters
                .generated_transitions
                .checked_add(1)
                .ok_or(LrCoefError::ArithmeticOverflow)?;
            next_shape[row] = base
                .checked_add(increment)
                .ok_or(LrCoefError::ArithmeticOverflow)?;
            next_prefix[row + 1] = prefix_sum
                .checked_add(increment)
                .ok_or(LrCoefError::ArithmeticOverflow)?;
            total = total
                .checked_add(self.enumerate_step(
                    step,
                    row + 1,
                    remaining - increment,
                    shape,
                    previous_prefix,
                    next_shape,
                    next_prefix,
                )?)
                .ok_or(LrCoefError::ArithmeticOverflow)?;
        }
        next_shape[row] = base;
        next_prefix[row + 1] = 0;
        Ok(total)
    }
}

impl<'a> BuchInteriorPruner<'a> {
    fn new(tight: &'a BuchTightFlags, inner: &[i32]) -> Result<Self, LrCoefError> {
        let rows = tight.lower.first().map_or(0, Vec::len);
        let steps = tight.lower.len();
        let inner = (0..rows)
            .map(|row| {
                let part = part_entry(inner, row);
                if part < 0 {
                    return Err(LrCoefError::InvalidPartition);
                }
                Ok(part as u32)
            })
            .collect::<Result<Vec<_>, _>>()?;
        Ok(Self {
            tight,
            inner,
            increments: vec![0; rows.saturating_mul(steps)],
            rows,
            steps,
        })
    }

    fn place(&mut self, cell: &LrCoefBox, row_complete: bool) -> Result<bool, LrCoefError> {
        let step = self.step_for_cell(cell)?;
        let index = self.index(step, cell.row)?;
        self.increments[index] = self.increments[index]
            .checked_add(1)
            .ok_or(LrCoefError::ArithmeticOverflow)?;

        if !self.partial_row_is_possible(cell.row) {
            return Ok(false);
        }
        if row_complete && !self.completed_row_is_relative_interior(cell.row) {
            return Ok(false);
        }
        Ok(true)
    }

    fn unplace(&mut self, cell: &LrCoefBox) -> Result<(), LrCoefError> {
        let step = self.step_for_cell(cell)?;
        let index = self.index(step, cell.row)?;
        self.increments[index] = self.increments[index]
            .checked_sub(1)
            .ok_or(LrCoefError::ArithmeticOverflow)?;
        Ok(())
    }

    fn is_complete_relative_interior(&self) -> bool {
        (0..self.rows).all(|row| self.completed_row_is_relative_interior(row))
    }

    fn partial_row_is_possible(&self, row: usize) -> bool {
        if row > 0 {
            for step in 0..self.steps {
                if !self.tight.diagonal[step][row]
                    && self.shape_after_step(row, step) >= self.shape_before_step(row - 1, step)
                {
                    return false;
                }
            }
        }
        for step in 1..self.steps {
            if !self.tight.yamanouchi[step][row]
                && self.prefix_for_step(step, row + 1) >= self.prefix_for_step(step - 1, row)
            {
                return false;
            }
        }
        true
    }

    fn completed_row_is_relative_interior(&self, row: usize) -> bool {
        for step in 0..self.steps {
            if !self.tight.lower[step][row] && self.increment(step, row) == 0 {
                return false;
            }
            if row > 0
                && !self.tight.diagonal[step][row]
                && self.shape_after_step(row, step) >= self.shape_before_step(row - 1, step)
            {
                return false;
            }
            if step > 0
                && !self.tight.yamanouchi[step][row]
                && self.prefix_for_step(step, row + 1) >= self.prefix_for_step(step - 1, row)
            {
                return false;
            }
        }
        true
    }

    fn step_for_cell(&self, cell: &LrCoefBox) -> Result<usize, LrCoefError> {
        let step = usize::try_from(cell.value - 1).map_err(|_| LrCoefError::ArithmeticOverflow)?;
        if step >= self.steps || cell.row >= self.rows {
            return Err(LrCoefError::ArithmeticOverflow);
        }
        Ok(step)
    }

    fn index(&self, step: usize, row: usize) -> Result<usize, LrCoefError> {
        step.checked_mul(self.rows)
            .and_then(|value| value.checked_add(row))
            .ok_or(LrCoefError::ArithmeticOverflow)
    }

    fn increment(&self, step: usize, row: usize) -> u32 {
        self.increments[step * self.rows + row]
    }

    fn shape_after_step(&self, row: usize, step: usize) -> u32 {
        self.inner[row] + (0..=step).map(|s| self.increment(s, row)).sum::<u32>()
    }

    fn shape_before_step(&self, row: usize, step: usize) -> u32 {
        self.inner[row] + (0..step).map(|s| self.increment(s, row)).sum::<u32>()
    }

    fn prefix_for_step(&self, step: usize, row_count: usize) -> u32 {
        (0..row_count).map(|row| self.increment(step, row)).sum()
    }
}

fn buch_dimension_from_tight_flags(tight: &BuchTightFlags) -> usize {
    let rows = tight.lower.first().map_or(0, Vec::len);
    let steps = tight.lower.len();
    let variable_count = rows * steps;
    let mut equations = Vec::<Vec<i32>>::new();

    let var = |step: usize, row: usize| -> usize { step * rows + row };

    for step in 0..steps {
        let mut equation = vec![0; variable_count];
        for row in 0..rows {
            equation[var(step, row)] = 1;
        }
        equations.push(equation);
    }

    for row in 0..rows {
        let mut equation = vec![0; variable_count];
        for step in 0..steps {
            equation[var(step, row)] = 1;
        }
        equations.push(equation);
    }

    for step in 0..steps {
        for row in 0..rows {
            if tight.lower[step][row] {
                let mut equation = vec![0; variable_count];
                equation[var(step, row)] = 1;
                equations.push(equation);
            }

            if row > 0 && tight.diagonal[step][row] {
                let mut equation = vec![0; variable_count];
                for earlier in 0..=step {
                    equation[var(earlier, row)] += 1;
                }
                for earlier in 0..step {
                    equation[var(earlier, row - 1)] -= 1;
                }
                equations.push(equation);
            }

            if step > 0 && tight.yamanouchi[step][row] {
                let mut equation = vec![0; variable_count];
                for prefix_row in 0..=row {
                    equation[var(step, prefix_row)] += 1;
                }
                for prefix_row in 0..row {
                    equation[var(step - 1, prefix_row)] -= 1;
                }
                equations.push(equation);
            }
        }
    }

    variable_count.saturating_sub(rational_rank(&equations))
}

fn rational_rank(equations: &[Vec<i32>]) -> usize {
    let Some(first) = equations.first() else {
        return 0;
    };
    let column_count = first.len();
    let mut matrix: Vec<Vec<BigRational>> = equations
        .iter()
        .filter(|row| row.iter().any(|&entry| entry != 0))
        .map(|row| {
            row.iter()
                .map(|&entry| BigRational::from_integer(entry.into()))
                .collect()
        })
        .collect();

    let mut rank = 0usize;
    for col in 0..column_count {
        let Some(pivot_row) = (rank..matrix.len()).find(|&row| !matrix[row][col].is_zero()) else {
            continue;
        };
        matrix.swap(rank, pivot_row);
        let pivot = matrix[rank][col].clone();
        for entry in matrix[rank].iter_mut().skip(col) {
            *entry /= pivot.clone();
        }
        let pivot_tail = matrix[rank][col..column_count].to_vec();
        for (row, row_entries) in matrix.iter_mut().enumerate() {
            if row == rank || row_entries[col].is_zero() {
                continue;
            }
            let factor = row_entries[col].clone();
            for (entry, pivot_entry) in row_entries[col..column_count].iter_mut().zip(&pivot_tail) {
                let sub = factor.clone() * pivot_entry;
                *entry -= sub;
            }
        }
        rank += 1;
        if rank == matrix.len() {
            break;
        }
    }
    rank
}

fn increment_at(increments: &[u32], rows: usize, step: usize, row: usize) -> u32 {
    increments[step * rows + row]
}

fn prefix_sum(increments: &[u32], rows: usize, step: usize, row_count: usize) -> u32 {
    (0..row_count)
        .map(|row| increment_at(increments, rows, step, row))
        .sum()
}

fn zero_buch_interior_stats() -> LrBuchInteriorStats {
    LrBuchInteriorStats {
        value: 0,
        weak_tableaux: 0,
        strict_tableaux: 0,
        strict_lower_constraints: 0,
        strict_diagonal_constraints: 0,
        strict_yamanouchi_constraints: 0,
    }
}

fn one_buch_interior_stats() -> LrBuchInteriorStats {
    LrBuchInteriorStats {
        value: 1,
        weak_tableaux: 1,
        strict_tableaux: 1,
        strict_lower_constraints: 0,
        strict_diagonal_constraints: 0,
        strict_yamanouchi_constraints: 0,
    }
}

fn zero_buch_interior_memo_stats() -> LrBuchInteriorMemoStats {
    LrBuchInteriorMemoStats {
        value: 0,
        weak_tableaux: 0,
        memo_states: 0,
        cache_hits: 0,
        generated_transitions: 0,
        strict_lower_constraints: 0,
        strict_diagonal_constraints: 0,
        strict_yamanouchi_constraints: 0,
    }
}

fn one_buch_interior_memo_stats() -> LrBuchInteriorMemoStats {
    LrBuchInteriorMemoStats {
        value: 1,
        weak_tableaux: 1,
        memo_states: 1,
        cache_hits: 0,
        generated_transitions: 0,
        strict_lower_constraints: 0,
        strict_diagonal_constraints: 0,
        strict_yamanouchi_constraints: 0,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::lr_gt::{lrcoef_gt_dimension, lrcoef_gt_interior_dfs_u128};

    #[test]
    fn computes_basic_coefficients() {
        assert_eq!(lrcoef(&[], &[], &[]), Ok(1));
        assert_eq!(lrcoef(&[2, 1], &[2], &[1]), Ok(1));
        assert_eq!(lrcoef(&[3, 2, 1], &[2, 1], &[2, 1]), Ok(2));
        assert_eq!(lrcoef(&[4, 2], &[2, 1], &[2, 1]), Ok(1));
        assert_eq!(lrcoef(&[5, 1], &[2, 1], &[2, 1]), Ok(0));
    }

    #[test]
    fn detects_invalid_partitions() {
        assert_eq!(
            lrcoef(&[2, 3], &[1], &[1]),
            Err(LrCoefError::InvalidPartition)
        );
        assert_eq!(
            lrcoef(&[2, -1], &[1], &[1]),
            Err(LrCoefError::InvalidPartition)
        );
    }

    #[test]
    fn trims_trailing_zeroes_like_partitions() {
        assert_eq!(lrcoef(&[3, 2, 1, 0], &[2, 1, 0], &[2, 1, 0]), Ok(2));
    }

    #[test]
    fn buch_interior_matches_gt_dfs_for_representative_cases() {
        let cases = [
            (&[3, 2, 1][..], &[2, 1][..], &[2, 1][..]),
            (&[4, 2][..], &[2, 1][..], &[2, 1][..]),
            (
                &[20, 14, 9, 5, 2][..],
                &[14, 9, 5, 2][..],
                &[8, 6, 4, 2][..],
            ),
            (
                &[40, 28, 18, 10, 4][..],
                &[28, 18, 10, 4][..],
                &[16, 12, 8, 4][..],
            ),
            (
                &[60, 42, 27, 15, 6][..],
                &[42, 27, 15, 6][..],
                &[24, 18, 12, 6][..],
            ),
        ];
        for (outer, inner, content) in cases {
            assert_eq!(
                lrcoef_buch_interior_u128(outer, inner, content).unwrap(),
                lrcoef_gt_interior_dfs_u128(outer, inner, content).unwrap(),
                "outer={outer:?} inner={inner:?} content={content:?}"
            );
        }
    }

    #[test]
    fn buch_interior_matches_gt_dfs_for_all_small_triples() {
        for outer_size in 0..=7 {
            for outer in partitions_of(outer_size) {
                for inner_size in 0..=outer_size {
                    for inner in partitions_of(inner_size) {
                        if !partition_less_equal_i32(&inner, &outer) {
                            continue;
                        }
                        let content_size = outer_size - inner_size;
                        for content in partitions_of(content_size) {
                            assert_eq!(
                                lrcoef_buch_interior_u128(&outer, &inner, &content).unwrap(),
                                lrcoef_gt_interior_dfs_u128(&outer, &inner, &content).unwrap(),
                                "outer={outer:?} inner={inner:?} content={content:?}"
                            );
                        }
                    }
                }
            }
        }
    }

    #[test]
    fn buch_memo_interior_matches_buch_search_for_representative_cases() {
        let cases = [
            (&[3, 2, 1][..], &[2, 1][..], &[2, 1][..]),
            (&[4, 2][..], &[2, 1][..], &[2, 1][..]),
            (
                &[20, 14, 9, 5, 2][..],
                &[14, 9, 5, 2][..],
                &[8, 6, 4, 2][..],
            ),
            (
                &[40, 28, 18, 10, 4][..],
                &[28, 18, 10, 4][..],
                &[16, 12, 8, 4][..],
            ),
            (
                &[7, 6, 5, 4, 3, 2, 1][..],
                &[4, 4, 3, 2, 1][..],
                &[5, 4, 3, 2][..],
            ),
        ];
        for (outer, inner, content) in cases {
            assert_eq!(
                lrcoef_buch_interior_memo_u128(outer, inner, content).unwrap(),
                lrcoef_buch_interior_u128(outer, inner, content).unwrap(),
                "outer={outer:?} inner={inner:?} content={content:?}"
            );
        }
    }

    #[test]
    fn buch_memo_interior_matches_buch_search_for_all_small_triples() {
        for outer_size in 0..=7 {
            for outer in partitions_of(outer_size) {
                for inner_size in 0..=outer_size {
                    for inner in partitions_of(inner_size) {
                        if !partition_less_equal_i32(&inner, &outer) {
                            continue;
                        }
                        let content_size = outer_size - inner_size;
                        for content in partitions_of(content_size) {
                            assert_eq!(
                                lrcoef_buch_interior_memo_u128(&outer, &inner, &content).unwrap(),
                                lrcoef_buch_interior_u128(&outer, &inner, &content).unwrap(),
                                "outer={outer:?} inner={inner:?} content={content:?}"
                            );
                        }
                    }
                }
            }
        }
    }

    #[test]
    fn buch_dimension_matches_gt_dimension_for_representative_cases() {
        let cases = [
            (&[3, 2, 1][..], &[2, 1][..], &[2, 1][..]),
            (&[4, 2][..], &[2, 1][..], &[2, 1][..]),
            (
                &[7, 6, 5, 4, 3, 2, 1][..],
                &[4, 4, 3, 2, 1][..],
                &[5, 4, 3, 2][..],
            ),
            (
                &[60, 42, 27, 15, 6][..],
                &[42, 27, 15, 6][..],
                &[24, 18, 12, 6][..],
            ),
        ];
        for (outer, inner, content) in cases {
            assert_eq!(
                lrcoef_buch_dimension(outer, inner, content).unwrap(),
                lrcoef_gt_dimension(outer, inner, content).unwrap(),
                "outer={outer:?} inner={inner:?} content={content:?}"
            );
        }
    }

    #[test]
    fn buch_dimension_matches_gt_dimension_for_all_small_triples() {
        for outer_size in 0..=7 {
            for outer in partitions_of(outer_size) {
                for inner_size in 0..=outer_size {
                    for inner in partitions_of(inner_size) {
                        if !partition_less_equal_i32(&inner, &outer) {
                            continue;
                        }
                        let content_size = outer_size - inner_size;
                        for content in partitions_of(content_size) {
                            assert_eq!(
                                lrcoef_buch_dimension(&outer, &inner, &content).unwrap(),
                                lrcoef_gt_dimension(&outer, &inner, &content).unwrap(),
                                "outer={outer:?} inner={inner:?} content={content:?}"
                            );
                        }
                    }
                }
            }
        }
    }

    #[test]
    fn stretch_cache_strict_counts_match_uncached_counts() {
        let cases = [
            (&[6, 4, 2][..], &[4, 2][..], &[4, 2][..], 1..=5),
            (
                &[7, 6, 5, 4, 3, 2, 1][..],
                &[4, 4, 3, 2, 1][..],
                &[5, 4, 3, 2][..],
                1..=4,
            ),
        ];
        for (outer, inner, content, stretches) in cases {
            let cache = lrcoef_buch_stretch_cache(outer, inner, content)
                .unwrap()
                .expect("nonempty stretch cache");
            for stretch in stretches {
                let cached = lrcoef_buch_stretched_counts_u128(&cache, stretch)
                    .unwrap()
                    .interior;
                let scaled_outer = scale_partition_i32(outer, stretch).unwrap();
                let scaled_inner = scale_partition_i32(inner, stretch).unwrap();
                let scaled_content = scale_partition_i32(content, stretch).unwrap();
                let uncached =
                    lrcoef_buch_interior_u128(&scaled_outer, &scaled_inner, &scaled_content)
                        .unwrap();
                assert_eq!(
                    cached, uncached,
                    "outer={outer:?} inner={inner:?} content={content:?} stretch={stretch}"
                );
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

    fn partition_less_equal_i32(inner: &[i32], outer: &[i32]) -> bool {
        let len = inner.len().max(outer.len());
        for index in 0..len {
            if *inner.get(index).unwrap_or(&0) > *outer.get(index).unwrap_or(&0) {
                return false;
            }
        }
        true
    }
}
