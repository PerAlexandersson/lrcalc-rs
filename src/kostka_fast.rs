//! Fast ordinary Kostka coefficients.
//!
//! This is a specialized single-coefficient DP over the usual chain of
//! horizontal strips.  It uses `u128` counts, packs each intermediate partition
//! into one `u128`, streams horizontal-strip successors without allocating
//! `Partition` objects, and memoizes transitions `(state, strip_size)`.

use std::collections::hash_map::Entry;
use std::collections::{HashMap, HashSet};

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum KostkaFastError {
    InvalidInput,
    ArithmeticOverflow,
    StateTooWide,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct KostkaFastStats {
    pub value: u128,
    pub peak_states: usize,
    pub levels: Vec<usize>,
    pub cached_transitions: usize,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct KostkaInteriorStats {
    pub value: u128,
    pub peak_states: usize,
    pub levels: Vec<usize>,
    pub full_reachable_levels: Vec<usize>,
    pub strict_lower_constraints: usize,
    pub strict_diagonal_constraints: usize,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct KostkaCountsStats {
    pub full: u128,
    pub interior: u128,
    pub full_peak_states: usize,
    pub full_levels: Vec<usize>,
    pub interior_peak_states: usize,
    pub interior_levels: Vec<usize>,
    pub full_reachable_levels: Vec<usize>,
    pub strict_lower_constraints: usize,
    pub strict_diagonal_constraints: usize,
}

#[derive(Clone, Copy, Debug)]
struct PackedPartitions {
    bits: u32,
    mask: u128,
}

#[derive(Clone, Debug)]
struct StripTransition {
    target: u128,
    increments: Vec<u32>,
}

#[derive(Clone, Debug)]
struct KostkaTightFlags {
    lower: Vec<Vec<bool>>,
    diagonal: Vec<Vec<bool>>,
}

#[derive(Clone, Debug)]
pub struct FastSkewKostkaEngine {
    outer: Vec<u32>,
    inner: Vec<u32>,
    packer: PackedPartitions,
    lambda_key: u128,
    inner_key: u128,
    skew_size: u32,
    transition_cache: HashMap<(u128, u32), Vec<u128>>,
}

/// Compute the ordinary Kostka coefficient `K_{shape, weight}`.
pub fn kostka_fast_u128(shape: &[i32], weight: &[i32]) -> Result<u128, KostkaFastError> {
    Ok(kostka_fast_stats(shape, weight)?.value)
}

/// Compute the skew Kostka coefficient `K_{outer/inner, weight}`.
pub fn skew_kostka_fast_u128(
    outer: &[i32],
    inner: &[i32],
    weight: &[i32],
) -> Result<u128, KostkaFastError> {
    Ok(skew_kostka_fast_stats(outer, inner, weight)?.value)
}

/// Count relative interior lattice points of the ordinary GT/Kostka polytope.
pub fn kostka_interior_u128(shape: &[i32], weight: &[i32]) -> Result<u128, KostkaFastError> {
    Ok(skew_kostka_interior_stats(shape, &[], weight)?.value)
}

/// Count relative interior lattice points of the skew GT/Kostka polytope.
pub fn skew_kostka_interior_u128(
    outer: &[i32],
    inner: &[i32],
    weight: &[i32],
) -> Result<u128, KostkaFastError> {
    Ok(skew_kostka_interior_stats(outer, inner, weight)?.value)
}

/// Compute full and relative-interior counts for `GT(shape, weight)`.
pub fn kostka_counts_stats(
    shape: &[i32],
    weight: &[i32],
) -> Result<KostkaCountsStats, KostkaFastError> {
    skew_kostka_counts_stats(shape, &[], weight)
}

/// Compute full and relative-interior counts for `GT(outer/inner, weight)`.
pub fn skew_kostka_counts_stats(
    outer: &[i32],
    inner: &[i32],
    weight: &[i32],
) -> Result<KostkaCountsStats, KostkaFastError> {
    let outer = normalize_partition(outer)?;
    let inner = normalize_partition(inner)?;
    let mut weight = normalize_weight(weight)?;
    let outer_size = checked_sum(&outer)?;
    let inner_size = checked_sum(&inner)?;
    if inner_size > outer_size || !partition_less_equal(&inner, &outer) {
        return Ok(zero_counts_stats());
    }
    let skew_size = outer_size - inner_size;
    let weight_size = checked_sum(&weight)?;
    if skew_size != weight_size {
        return Ok(zero_counts_stats());
    }
    if skew_size == 0 {
        return Ok(if trim_eq(&outer, &inner) {
            one_counts_stats()
        } else {
            zero_counts_stats()
        });
    }
    if outer.is_empty() {
        return Ok(zero_counts_stats());
    }

    weight.sort_unstable_by(|a, b| b.cmp(a));
    if inner.is_empty() && !dominates(&outer, &weight) {
        return Ok(zero_counts_stats());
    }

    let packer = PackedPartitions::new(&outer)?;
    let inner_key = packer.pack_padded(&inner, outer.len())?;
    let outer_key = packer.pack(&outer)?;

    let mut dp: HashMap<u128, u128> = HashMap::new();
    dp.insert(inner_key, 1);
    let mut reachable = Vec::<HashSet<u128>>::with_capacity(weight.len() + 1);
    reachable.push(HashSet::from([inner_key]));
    let mut full_peak_states = dp.len();
    let mut full_levels = vec![dp.len()];

    for &strip_size in &weight {
        let mut next: HashMap<u128, u128> = HashMap::with_capacity(dp.len().saturating_mul(2));
        let mut next_reachable = HashSet::new();
        for (&state, &count) in &dp {
            enumerate_strip_transitions(packer, &outer, state, strip_size, |transition| {
                let entry = next.entry(transition.target).or_insert(0);
                *entry = entry
                    .checked_add(count)
                    .ok_or(KostkaFastError::ArithmeticOverflow)?;
                next_reachable.insert(transition.target);
                Ok(())
            })?;
        }
        full_peak_states = full_peak_states.max(next.len());
        full_levels.push(next.len());
        reachable.push(next_reachable);
        dp = next;
    }

    let full = dp.remove(&outer_key).unwrap_or(0);
    let full_reachable_levels = reachable.iter().map(HashSet::len).collect::<Vec<_>>();
    let coreachable = kostka_coreachable_levels(packer, &outer, &weight, outer_key, &reachable)?;
    if !coreachable[0].contains(&inner_key) {
        return Ok(KostkaCountsStats {
            full,
            full_peak_states,
            full_levels,
            full_reachable_levels,
            ..zero_counts_stats()
        });
    }

    let tight = kostka_tight_flags_on_complete_paths(packer, &outer, &weight, &coreachable)?;
    let (strict_lower, strict_diagonal) = tight.strict_counts();

    let mut strict_dp: HashMap<u128, u128> = HashMap::new();
    strict_dp.insert(inner_key, 1);
    let mut interior_peak_states = strict_dp.len();
    let mut interior_levels = vec![strict_dp.len()];

    for (step, &strip_size) in weight.iter().enumerate() {
        let mut next: HashMap<u128, u128> =
            HashMap::with_capacity(strict_dp.len().saturating_mul(2));
        for (&state, &count) in &strict_dp {
            enumerate_strip_transitions(packer, &outer, state, strip_size, |transition| {
                if !coreachable[step + 1].contains(&transition.target)
                    || !kostka_transition_is_relative_interior(
                        packer,
                        step,
                        state,
                        &transition,
                        &tight,
                    )
                {
                    return Ok(());
                }
                let entry = next.entry(transition.target).or_insert(0);
                *entry = entry
                    .checked_add(count)
                    .ok_or(KostkaFastError::ArithmeticOverflow)?;
                Ok(())
            })?;
        }
        interior_peak_states = interior_peak_states.max(next.len());
        interior_levels.push(next.len());
        strict_dp = next;
    }

    Ok(KostkaCountsStats {
        full,
        interior: strict_dp.remove(&outer_key).unwrap_or(0),
        full_peak_states,
        full_levels,
        interior_peak_states,
        interior_levels,
        full_reachable_levels,
        strict_lower_constraints: strict_lower,
        strict_diagonal_constraints: strict_diagonal,
    })
}

/// Compute `K_{shape, weight}` and return DP statistics useful for profiling.
pub fn kostka_fast_stats(
    shape: &[i32],
    weight: &[i32],
) -> Result<KostkaFastStats, KostkaFastError> {
    let shape = normalize_partition(shape)?;
    let mut weight = normalize_weight(weight)?;
    let shape_size = checked_sum(&shape)?;
    let weight_size = checked_sum(&weight)?;

    if shape_size != weight_size {
        return Ok(zero_stats());
    }
    if shape_size == 0 {
        return Ok(one_stats());
    }

    weight.sort_unstable_by(|a, b| b.cmp(a));
    if !dominates(&shape, &weight) {
        return Ok(zero_stats());
    }

    if shape.len() == 1 {
        return Ok(one_stats());
    }
    if shape.iter().all(|&part| part == 1) {
        return Ok(if weight.iter().all(|&part| part == 1) {
            one_stats()
        } else {
            zero_stats()
        });
    }
    if weight.iter().all(|&part| part == 1) {
        return Ok(KostkaFastStats {
            value: count_standard_tableaux(&shape)?,
            peak_states: 0,
            levels: Vec::new(),
            cached_transitions: 0,
        });
    }

    let mut engine = FastSkewKostkaEngine::from_normalized(shape, Vec::new())?;
    engine.stats_normalized(weight)
}

/// Compute `K_{outer/inner, weight}` and return DP statistics useful for profiling.
pub fn skew_kostka_fast_stats(
    outer: &[i32],
    inner: &[i32],
    weight: &[i32],
) -> Result<KostkaFastStats, KostkaFastError> {
    let outer = normalize_partition(outer)?;
    let inner = normalize_partition(inner)?;
    let mut weight = normalize_weight(weight)?;
    let outer_size = checked_sum(&outer)?;
    let inner_size = checked_sum(&inner)?;
    if inner_size > outer_size || !partition_less_equal(&inner, &outer) {
        return Ok(zero_stats());
    }
    let skew_size = outer_size - inner_size;
    let weight_size = checked_sum(&weight)?;

    if skew_size != weight_size {
        return Ok(zero_stats());
    }
    if skew_size == 0 {
        return Ok(if trim_eq(&outer, &inner) {
            one_stats()
        } else {
            zero_stats()
        });
    }

    weight.sort_unstable_by(|a, b| b.cmp(a));
    if inner.is_empty() {
        return kostka_fast_stats(
            &outer
                .iter()
                .copied()
                .map(|part| part as i32)
                .collect::<Vec<_>>(),
            &weight
                .iter()
                .copied()
                .map(|part| part as i32)
                .collect::<Vec<_>>(),
        );
    }

    let mut engine = FastSkewKostkaEngine::from_normalized(outer, inner)?;
    engine.stats_normalized(weight)
}

/// Count interior points of `GT(shape, weight)` with DP statistics.
pub fn kostka_interior_stats(
    shape: &[i32],
    weight: &[i32],
) -> Result<KostkaInteriorStats, KostkaFastError> {
    skew_kostka_interior_stats(shape, &[], weight)
}

/// Count interior points of `GT(outer/inner, weight)` with DP statistics.
pub fn skew_kostka_interior_stats(
    outer: &[i32],
    inner: &[i32],
    weight: &[i32],
) -> Result<KostkaInteriorStats, KostkaFastError> {
    let outer = normalize_partition(outer)?;
    let inner = normalize_partition(inner)?;
    let weight = normalize_weight(weight)?;
    let outer_size = checked_sum(&outer)?;
    let inner_size = checked_sum(&inner)?;
    if inner_size > outer_size || !partition_less_equal(&inner, &outer) {
        return Ok(zero_interior_stats());
    }
    let skew_size = outer_size - inner_size;
    let weight_size = checked_sum(&weight)?;
    if skew_size != weight_size {
        return Ok(zero_interior_stats());
    }
    if skew_size == 0 {
        return Ok(if trim_eq(&outer, &inner) {
            one_interior_stats()
        } else {
            zero_interior_stats()
        });
    }
    if outer.is_empty() {
        return Ok(zero_interior_stats());
    }

    let packer = PackedPartitions::new(&outer)?;
    let inner_key = packer.pack_padded(&inner, outer.len())?;
    let outer_key = packer.pack(&outer)?;

    let reachable = kostka_reachable_levels(packer, &outer, &weight, inner_key)?;
    let full_reachable_levels = reachable.iter().map(HashSet::len).collect::<Vec<_>>();
    let coreachable = kostka_coreachable_levels(packer, &outer, &weight, outer_key, &reachable)?;
    if !coreachable[0].contains(&inner_key) {
        return Ok(KostkaInteriorStats {
            full_reachable_levels,
            ..zero_interior_stats()
        });
    }

    let tight = kostka_tight_flags_on_complete_paths(packer, &outer, &weight, &coreachable)?;
    let (strict_lower, strict_diagonal) = tight.strict_counts();

    let mut dp: HashMap<u128, u128> = HashMap::new();
    dp.insert(inner_key, 1);
    let mut peak_states = dp.len();
    let mut levels = vec![dp.len()];

    for (step, &strip_size) in weight.iter().enumerate() {
        let mut next: HashMap<u128, u128> = HashMap::with_capacity(dp.len().saturating_mul(2));
        for (&state, &count) in &dp {
            enumerate_strip_transitions(packer, &outer, state, strip_size, |transition| {
                if !kostka_transition_is_relative_interior(packer, step, state, &transition, &tight)
                {
                    return Ok(());
                }
                let entry = next.entry(transition.target).or_insert(0);
                *entry = entry
                    .checked_add(count)
                    .ok_or(KostkaFastError::ArithmeticOverflow)?;
                Ok(())
            })?;
        }
        peak_states = peak_states.max(next.len());
        levels.push(next.len());
        dp = next;
    }

    Ok(KostkaInteriorStats {
        value: dp.remove(&outer_key).unwrap_or(0),
        peak_states,
        levels,
        full_reachable_levels,
        strict_lower_constraints: strict_lower,
        strict_diagonal_constraints: strict_diagonal,
    })
}

impl PackedPartitions {
    fn new(lambda: &[u32]) -> Result<Self, KostkaFastError> {
        let Some(&max_part) = lambda.first() else {
            return Ok(Self { bits: 1, mask: 1 });
        };
        let bits = bit_width(max_part);
        let width = (bits as usize)
            .checked_mul(lambda.len())
            .ok_or(KostkaFastError::ArithmeticOverflow)?;
        if width > 128 {
            return Err(KostkaFastError::StateTooWide);
        }
        let mask = (1u128 << bits) - 1;
        Ok(Self { bits, mask })
    }

    fn get(self, key: u128, row: usize) -> u32 {
        let shift = self.bits * row as u32;
        ((key >> shift) & self.mask) as u32
    }

    fn set(self, key: u128, row: usize, value: u32) -> u128 {
        let shift = self.bits * row as u32;
        let row_mask = self.mask << shift;
        (key & !row_mask) | ((value as u128) << shift)
    }

    fn pack(self, partition: &[u32]) -> Result<u128, KostkaFastError> {
        let mut key = 0u128;
        for (row, &part) in partition.iter().enumerate() {
            if part as u128 > self.mask {
                return Err(KostkaFastError::StateTooWide);
            }
            key = self.set(key, row, part);
        }
        Ok(key)
    }

    fn pack_padded(self, partition: &[u32], rows: usize) -> Result<u128, KostkaFastError> {
        if partition.len() > rows && partition[rows..].iter().any(|&part| part != 0) {
            return Err(KostkaFastError::InvalidInput);
        }
        self.pack(&partition[..partition.len().min(rows)])
    }
}

fn kostka_reachable_levels(
    packer: PackedPartitions,
    outer: &[u32],
    weight: &[u32],
    inner_key: u128,
) -> Result<Vec<HashSet<u128>>, KostkaFastError> {
    let mut levels = Vec::with_capacity(weight.len() + 1);
    let mut current = HashSet::new();
    current.insert(inner_key);
    levels.push(current);

    for &strip_size in weight {
        let mut next = HashSet::new();
        for &state in levels.last().expect("at least one level") {
            enumerate_strip_transitions(packer, outer, state, strip_size, |transition| {
                next.insert(transition.target);
                Ok(())
            })?;
        }
        levels.push(next);
    }

    Ok(levels)
}

fn kostka_coreachable_levels(
    packer: PackedPartitions,
    outer: &[u32],
    weight: &[u32],
    outer_key: u128,
    reachable: &[HashSet<u128>],
) -> Result<Vec<HashSet<u128>>, KostkaFastError> {
    let mut coreachable = vec![HashSet::new(); weight.len() + 1];
    if let Some(final_level) = reachable.last() {
        for &state in final_level {
            if state == outer_key {
                coreachable[weight.len()].insert(state);
            }
        }
    }

    for step in (0..weight.len()).rev() {
        for &state in &reachable[step] {
            let mut reaches_final = false;
            enumerate_strip_transitions(packer, outer, state, weight[step], |transition| {
                if coreachable[step + 1].contains(&transition.target) {
                    reaches_final = true;
                }
                Ok(())
            })?;
            if reaches_final {
                coreachable[step].insert(state);
            }
        }
    }

    Ok(coreachable)
}

fn kostka_tight_flags_on_complete_paths(
    packer: PackedPartitions,
    outer: &[u32],
    weight: &[u32],
    coreachable: &[HashSet<u128>],
) -> Result<KostkaTightFlags, KostkaFastError> {
    let mut tight = KostkaTightFlags::new(weight.len(), outer.len());
    for (step, &strip_size) in weight.iter().enumerate() {
        for &state in &coreachable[step] {
            enumerate_strip_transitions(packer, outer, state, strip_size, |transition| {
                if coreachable[step + 1].contains(&transition.target) {
                    tight.observe_transition(packer, step, state, &transition);
                }
                Ok(())
            })?;
        }
    }
    Ok(tight)
}

fn kostka_transition_is_relative_interior(
    packer: PackedPartitions,
    step: usize,
    source: u128,
    transition: &StripTransition,
    tight: &KostkaTightFlags,
) -> bool {
    for row in 0..transition.increments.len() {
        if !tight.lower[step][row] && transition.increments[row] == 0 {
            return false;
        }
        if row > 0 && !tight.diagonal[step][row] {
            let target_value = packer.get(transition.target, row);
            let diagonal_bound = packer.get(source, row - 1);
            if target_value == diagonal_bound {
                return false;
            }
        }
    }
    true
}

impl KostkaTightFlags {
    fn new(steps: usize, rows: usize) -> Self {
        Self {
            lower: vec![vec![true; rows]; steps],
            diagonal: vec![vec![true; rows]; steps],
        }
    }

    fn observe_transition(
        &mut self,
        packer: PackedPartitions,
        step: usize,
        source: u128,
        transition: &StripTransition,
    ) {
        for row in 0..transition.increments.len() {
            if transition.increments[row] > 0 {
                self.lower[step][row] = false;
            }
            if row > 0 {
                let target_value = packer.get(transition.target, row);
                let diagonal_bound = packer.get(source, row - 1);
                if target_value < diagonal_bound {
                    self.diagonal[step][row] = false;
                }
            }
        }
    }

    fn strict_counts(&self) -> (usize, usize) {
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
        (lower, diagonal)
    }
}

impl FastSkewKostkaEngine {
    pub fn new(outer: &[i32], inner: &[i32]) -> Result<Self, KostkaFastError> {
        let outer = normalize_partition(outer)?;
        let inner = normalize_partition(inner)?;
        let outer_size = checked_sum(&outer)?;
        let inner_size = checked_sum(&inner)?;
        if inner_size > outer_size || !partition_less_equal(&inner, &outer) {
            return Err(KostkaFastError::InvalidInput);
        }
        Self::from_normalized(outer, inner)
    }

    pub fn coefficient(&mut self, weight: &[i32]) -> Result<u128, KostkaFastError> {
        Ok(self.stats(weight)?.value)
    }

    pub fn stats(&mut self, weight: &[i32]) -> Result<KostkaFastStats, KostkaFastError> {
        let mut weight = normalize_weight(weight)?;
        weight.sort_unstable_by(|a, b| b.cmp(a));
        self.stats_normalized(weight)
    }

    pub fn cached_transitions(&self) -> usize {
        self.transition_cache.len()
    }

    fn from_normalized(outer: Vec<u32>, inner: Vec<u32>) -> Result<Self, KostkaFastError> {
        let packer = PackedPartitions::new(&outer)?;
        let lambda_key = packer.pack(&outer)?;
        let inner_key = packer.pack_padded(&inner, outer.len())?;
        let outer_size = checked_sum(&outer)?;
        let inner_size = checked_sum(&inner)?;
        Ok(Self {
            outer,
            inner,
            packer,
            lambda_key,
            inner_key,
            skew_size: outer_size - inner_size,
            transition_cache: HashMap::new(),
        })
    }

    fn stats_normalized(&mut self, weight: Vec<u32>) -> Result<KostkaFastStats, KostkaFastError> {
        let weight_size = checked_sum(&weight)?;
        if weight_size != self.skew_size {
            return Ok(zero_stats());
        }
        if self.skew_size == 0 {
            return Ok(if trim_eq(&self.outer, &self.inner) {
                one_stats()
            } else {
                zero_stats()
            });
        }
        if weight.len() == 1 {
            return Ok(KostkaFastStats {
                value: u128::from(self.is_horizontal_strip(weight[0])),
                peak_states: 0,
                levels: Vec::new(),
                cached_transitions: self.transition_cache.len(),
            });
        }
        if weight.len() == 2 {
            if let Some(value) = self.two_step_count(weight[0])? {
                return Ok(KostkaFastStats {
                    value,
                    peak_states: 0,
                    levels: Vec::new(),
                    cached_transitions: self.transition_cache.len(),
                });
            }
        }

        let mut dp: HashMap<u128, u128> = HashMap::new();
        dp.insert(self.inner_key, 1);

        let mut peak_states = dp.len();
        let mut levels = vec![dp.len()];

        for strip_size in weight {
            if strip_size == 0 {
                levels.push(dp.len());
                continue;
            }

            let mut next: HashMap<u128, u128> = HashMap::with_capacity(dp.len().saturating_mul(2));
            for (&state, &count) in &dp {
                let successors = match self.transition_cache.entry((state, strip_size)) {
                    Entry::Occupied(entry) => entry.into_mut(),
                    Entry::Vacant(entry) => {
                        let successors = successors(&self.packer, &self.outer, state, strip_size)?;
                        entry.insert(successors)
                    }
                };
                for &target in successors.iter() {
                    let entry = next.entry(target).or_insert(0u128);
                    *entry = (*entry)
                        .checked_add(count)
                        .ok_or(KostkaFastError::ArithmeticOverflow)?;
                }
            }

            peak_states = peak_states.max(next.len());
            levels.push(next.len());
            dp = next;
        }

        Ok(KostkaFastStats {
            value: dp.remove(&self.lambda_key).unwrap_or(0),
            peak_states,
            levels,
            cached_transitions: self.transition_cache.len(),
        })
    }

    fn is_horizontal_strip(&self, strip_size: u32) -> bool {
        if strip_size != self.skew_size {
            return false;
        }
        for row in 0..self.outer.len() {
            let outer = self.outer[row];
            let inner = self.part_inner(row);
            if outer < inner {
                return false;
            }
            let increment = outer - inner;
            if row > 0 {
                let max_increment = self.part_inner(row - 1).saturating_sub(inner);
                if increment > max_increment {
                    return false;
                }
            }
        }
        true
    }

    fn two_step_count(&self, first_strip_size: u32) -> Result<Option<u128>, KostkaFastError> {
        if self.outer.len() > 20 {
            return Ok(None);
        }

        let inner_size = checked_sum(&self.inner)?;
        let target_sum = inner_size
            .checked_add(first_strip_size)
            .ok_or(KostkaFastError::ArithmeticOverflow)?;
        let mut lower_sum = 0u32;
        let mut caps = Vec::with_capacity(self.outer.len());

        for row in 0..self.outer.len() {
            let lower = self.part_inner(row).max(self.part_outer(row + 1));
            let upper = self.outer[row].min(if row == 0 {
                self.outer[row]
            } else {
                self.part_inner(row - 1)
            });
            if lower > upper {
                return Ok(Some(0));
            }
            lower_sum = lower_sum
                .checked_add(lower)
                .ok_or(KostkaFastError::ArithmeticOverflow)?;
            caps.push(upper - lower);
        }

        if target_sum < lower_sum {
            return Ok(Some(0));
        }
        let target = target_sum - lower_sum;
        Ok(Some(count_bounded_sum(&caps, target)?))
    }

    fn part_outer(&self, row: usize) -> u32 {
        *self.outer.get(row).unwrap_or(&0)
    }

    fn part_inner(&self, row: usize) -> u32 {
        *self.inner.get(row).unwrap_or(&0)
    }
}

fn successors(
    packer: &PackedPartitions,
    lambda: &[u32],
    state: u128,
    strip_size: u32,
) -> Result<Vec<u128>, KostkaFastError> {
    let mut out = Vec::new();
    enumerate_successors(*packer, lambda, state, state, strip_size, 0, &mut out)?;
    Ok(out)
}

fn enumerate_successors(
    packer: PackedPartitions,
    lambda: &[u32],
    alpha_key: u128,
    partial_key: u128,
    remaining: u32,
    row: usize,
    out: &mut Vec<u128>,
) -> Result<(), KostkaFastError> {
    if row == lambda.len() {
        if remaining == 0 {
            out.push(partial_key);
        }
        return Ok(());
    }

    let base = packer.get(alpha_key, row);
    let max_from_lambda = lambda[row].saturating_sub(base);
    let max_from_strip = if row == 0 {
        remaining
    } else {
        packer.get(alpha_key, row - 1).saturating_sub(base)
    };
    let max_c = remaining.min(max_from_lambda).min(max_from_strip);

    for c in 0..=max_c {
        let value = base
            .checked_add(c)
            .ok_or(KostkaFastError::ArithmeticOverflow)?;
        let next_key = packer.set(partial_key, row, value);
        enumerate_successors(
            packer,
            lambda,
            alpha_key,
            next_key,
            remaining - c,
            row + 1,
            out,
        )?;
    }
    Ok(())
}

fn enumerate_strip_transitions<F>(
    packer: PackedPartitions,
    lambda: &[u32],
    alpha_key: u128,
    strip_size: u32,
    mut visit: F,
) -> Result<(), KostkaFastError>
where
    F: FnMut(StripTransition) -> Result<(), KostkaFastError>,
{
    let mut increments = vec![0; lambda.len()];
    enumerate_strip_transitions_recursive(
        packer,
        lambda,
        alpha_key,
        alpha_key,
        strip_size,
        0,
        &mut increments,
        &mut visit,
    )
}

#[allow(clippy::too_many_arguments)]
fn enumerate_strip_transitions_recursive<F>(
    packer: PackedPartitions,
    lambda: &[u32],
    alpha_key: u128,
    partial_key: u128,
    remaining: u32,
    row: usize,
    increments: &mut [u32],
    visit: &mut F,
) -> Result<(), KostkaFastError>
where
    F: FnMut(StripTransition) -> Result<(), KostkaFastError>,
{
    if row == lambda.len() {
        if remaining == 0 {
            visit(StripTransition {
                target: partial_key,
                increments: increments.to_vec(),
            })?;
        }
        return Ok(());
    }

    let base = packer.get(alpha_key, row);
    let max_from_lambda = lambda[row].saturating_sub(base);
    let max_from_strip = if row == 0 {
        remaining
    } else {
        packer.get(alpha_key, row - 1).saturating_sub(base)
    };
    let max_c = remaining.min(max_from_lambda).min(max_from_strip);

    for c in 0..=max_c {
        increments[row] = c;
        let value = base
            .checked_add(c)
            .ok_or(KostkaFastError::ArithmeticOverflow)?;
        let next_key = packer.set(partial_key, row, value);
        enumerate_strip_transitions_recursive(
            packer,
            lambda,
            alpha_key,
            next_key,
            remaining - c,
            row + 1,
            increments,
            visit,
        )?;
    }
    increments[row] = 0;
    Ok(())
}

fn normalize_partition(parts: &[i32]) -> Result<Vec<u32>, KostkaFastError> {
    if !parts.iter().all(|&part| part >= 0) || !parts.windows(2).all(|w| w[0] >= w[1]) {
        return Err(KostkaFastError::InvalidInput);
    }
    Ok(parts
        .iter()
        .copied()
        .take_while(|&part| part != 0)
        .map(|part| part as u32)
        .collect())
}

fn normalize_weight(parts: &[i32]) -> Result<Vec<u32>, KostkaFastError> {
    if !parts.iter().all(|&part| part >= 0) {
        return Err(KostkaFastError::InvalidInput);
    }
    Ok(parts
        .iter()
        .copied()
        .filter(|&part| part > 0)
        .map(|part| part as u32)
        .collect())
}

fn checked_sum(parts: &[u32]) -> Result<u32, KostkaFastError> {
    let mut sum = 0u32;
    for &part in parts {
        sum = sum
            .checked_add(part)
            .ok_or(KostkaFastError::ArithmeticOverflow)?;
    }
    Ok(sum)
}

fn dominates(shape: &[u32], weight: &[u32]) -> bool {
    let len = shape.len().max(weight.len());
    let mut shape_sum = 0u32;
    let mut weight_sum = 0u32;
    for i in 0..len {
        shape_sum = shape_sum.saturating_add(*shape.get(i).unwrap_or(&0));
        weight_sum = weight_sum.saturating_add(*weight.get(i).unwrap_or(&0));
        if shape_sum < weight_sum {
            return false;
        }
    }
    true
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

fn count_standard_tableaux(shape: &[u32]) -> Result<u128, KostkaFastError> {
    let n = checked_sum(shape)?;
    let mut factors: Vec<u128> = (2..=n).map(u128::from).collect();

    for hook in hook_lengths(shape) {
        let mut divisor = hook as u128;
        for factor in &mut factors {
            let g = gcd(*factor, divisor);
            if g > 1 {
                *factor /= g;
                divisor /= g;
                if divisor == 1 {
                    break;
                }
            }
        }
        if divisor != 1 {
            return Err(KostkaFastError::ArithmeticOverflow);
        }
    }

    let mut result = 1u128;
    for factor in factors {
        result = result
            .checked_mul(factor)
            .ok_or(KostkaFastError::ArithmeticOverflow)?;
    }
    Ok(result)
}

fn hook_lengths(shape: &[u32]) -> Vec<u32> {
    let mut hooks = Vec::new();
    for (row, &row_len) in shape.iter().enumerate() {
        for col in 0..row_len {
            let arm = row_len - col - 1;
            let leg = shape[row + 1..]
                .iter()
                .filter(|&&below_len| below_len > col)
                .count() as u32;
            hooks.push(arm + leg + 1);
        }
    }
    hooks
}

fn gcd(mut a: u128, mut b: u128) -> u128 {
    while b != 0 {
        let r = a % b;
        a = b;
        b = r;
    }
    a
}

fn bit_width(value: u32) -> u32 {
    (u32::BITS - value.leading_zeros()).max(1)
}

fn count_bounded_sum(caps: &[u32], target: u32) -> Result<u128, KostkaFastError> {
    if caps.is_empty() {
        return Ok(u128::from(target == 0));
    }

    let subsets = 1usize
        .checked_shl(caps.len() as u32)
        .ok_or(KostkaFastError::ArithmeticOverflow)?;
    let mut positive = 0u128;
    let mut negative = 0u128;
    for mask in 0..subsets {
        let mut shift = 0u32;
        let mut odd = false;
        for (index, &cap) in caps.iter().enumerate() {
            if (mask >> index) & 1 == 1 {
                shift = shift
                    .checked_add(cap)
                    .and_then(|value| value.checked_add(1))
                    .ok_or(KostkaFastError::ArithmeticOverflow)?;
                odd = !odd;
            }
        }
        if shift > target {
            continue;
        }
        let remaining = target - shift;
        let ways = binomial_u128(
            remaining as u64 + caps.len() as u64 - 1,
            caps.len() as u64 - 1,
        )?;
        if odd {
            negative = negative
                .checked_add(ways)
                .ok_or(KostkaFastError::ArithmeticOverflow)?;
        } else {
            positive = positive
                .checked_add(ways)
                .ok_or(KostkaFastError::ArithmeticOverflow)?;
        }
    }
    positive
        .checked_sub(negative)
        .ok_or(KostkaFastError::ArithmeticOverflow)
}

fn binomial_u128(n: u64, k: u64) -> Result<u128, KostkaFastError> {
    let k = k.min(n - k);
    let mut result = 1u128;
    for i in 1..=k {
        result = result
            .checked_mul((n - k + i) as u128)
            .ok_or(KostkaFastError::ArithmeticOverflow)?
            / i as u128;
    }
    Ok(result)
}

fn zero_stats() -> KostkaFastStats {
    KostkaFastStats {
        value: 0,
        peak_states: 0,
        levels: Vec::new(),
        cached_transitions: 0,
    }
}

fn one_stats() -> KostkaFastStats {
    KostkaFastStats {
        value: 1,
        peak_states: 1,
        levels: Vec::new(),
        cached_transitions: 0,
    }
}

fn zero_interior_stats() -> KostkaInteriorStats {
    KostkaInteriorStats {
        value: 0,
        peak_states: 0,
        levels: Vec::new(),
        full_reachable_levels: Vec::new(),
        strict_lower_constraints: 0,
        strict_diagonal_constraints: 0,
    }
}

fn one_interior_stats() -> KostkaInteriorStats {
    KostkaInteriorStats {
        value: 1,
        peak_states: 1,
        levels: Vec::new(),
        full_reachable_levels: Vec::new(),
        strict_lower_constraints: 0,
        strict_diagonal_constraints: 0,
    }
}

fn zero_counts_stats() -> KostkaCountsStats {
    KostkaCountsStats {
        full: 0,
        interior: 0,
        full_peak_states: 0,
        full_levels: Vec::new(),
        interior_peak_states: 0,
        interior_levels: Vec::new(),
        full_reachable_levels: Vec::new(),
        strict_lower_constraints: 0,
        strict_diagonal_constraints: 0,
    }
}

fn one_counts_stats() -> KostkaCountsStats {
    KostkaCountsStats {
        full: 1,
        interior: 1,
        full_peak_states: 1,
        full_levels: vec![1],
        interior_peak_states: 1,
        interior_levels: vec![1],
        full_reachable_levels: vec![1],
        strict_lower_constraints: 0,
        strict_diagonal_constraints: 0,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::kostka::{kostka_lr_triple, kostka_via_lr};
    use crate::lr_gt::lrcoef_gt_interior_u128;

    #[test]
    fn computes_small_coefficients() {
        assert_eq!(kostka_fast_u128(&[], &[]), Ok(1));
        assert_eq!(kostka_fast_u128(&[3], &[2, 1]), Ok(1));
        assert_eq!(kostka_fast_u128(&[2, 1], &[2, 1]), Ok(1));
        assert_eq!(kostka_fast_u128(&[2, 1], &[1, 1, 1]), Ok(2));
        assert_eq!(kostka_fast_u128(&[1, 1, 1], &[2, 1]), Ok(0));
    }

    #[test]
    fn computes_skew_coefficients() {
        assert_eq!(skew_kostka_fast_u128(&[3], &[1], &[1, 1]), Ok(1));
        assert_eq!(skew_kostka_fast_u128(&[1, 1], &[], &[1, 1]), Ok(1));
        assert_eq!(skew_kostka_fast_u128(&[2, 1], &[1], &[1, 1]), Ok(2));
    }

    #[test]
    fn weight_order_does_not_matter() {
        assert_eq!(kostka_fast_u128(&[3, 1], &[2, 1, 1]), Ok(2));
        assert_eq!(kostka_fast_u128(&[3, 1], &[1, 2, 1]), Ok(2));
    }

    #[test]
    fn handles_standard_tableaux_special_case() {
        assert_eq!(kostka_fast_u128(&[3, 2, 1], &[1, 1, 1, 1, 1, 1]), Ok(16));
        assert_eq!(
            kostka_fast_u128(&[4, 3, 2, 1], &[1, 1, 1, 1, 1, 1, 1, 1, 1, 1]),
            Ok(768)
        );
    }

    #[test]
    fn agrees_with_lr_translation_for_representative_cases() {
        let cases = [
            (&[3, 2, 1][..], &[2, 2, 2][..]),
            (&[4, 2, 1][..], &[3, 2, 1, 1][..]),
            (&[6, 4, 2][..], &[4, 4, 4][..]),
            (&[5, 4, 2, 1][..], &[4, 3, 2, 2, 1][..]),
            (&[8, 6, 4, 2][..], &[6, 5, 4, 3, 2][..]),
        ];
        for (shape, weight) in cases {
            assert_eq!(
                kostka_fast_u128(shape, weight).unwrap(),
                kostka_via_lr(shape, weight).unwrap()
            );
        }
    }

    #[test]
    fn computes_basic_interior_coefficients() {
        assert_eq!(kostka_interior_u128(&[], &[]), Ok(1));
        assert_eq!(kostka_interior_u128(&[2, 1], &[1, 1, 1]), Ok(0));
        assert_eq!(kostka_interior_u128(&[4, 2], &[2, 2, 1, 1]), Ok(0));
        assert_eq!(kostka_interior_u128(&[6, 4, 2], &[4, 4, 4]), Ok(1));
    }

    #[test]
    fn paired_counts_match_separate_kostka_paths() {
        let cases = [
            (&[3, 2, 1][..], &[2, 2, 2][..]),
            (&[6, 4, 2][..], &[4, 4, 4][..]),
            (&[8, 6, 4, 2][..], &[6, 5, 4, 3, 2][..]),
        ];
        for (shape, weight) in cases {
            let counts = kostka_counts_stats(shape, weight).unwrap();
            assert_eq!(counts.full, kostka_fast_u128(shape, weight).unwrap());
            assert_eq!(
                counts.interior,
                kostka_interior_u128(shape, weight).unwrap()
            );
        }
    }

    #[test]
    fn interior_matches_translated_lr_in_representative_cases() {
        let cases = [
            (&[2, 1][..], &[1, 1, 1][..]),
            (&[3, 1][..], &[2, 1, 1][..]),
            (&[4, 2][..], &[2, 2, 1, 1][..]),
            (&[3, 2, 1][..], &[2, 2, 2][..]),
            (&[5, 3, 1][..], &[3, 2, 2, 1, 1][..]),
            (&[6, 4, 2][..], &[4, 4, 4][..]),
        ];
        for (shape, weight) in cases {
            let (outer, inner, content) = kostka_lr_triple(shape, weight).unwrap();
            assert_eq!(
                kostka_interior_u128(shape, weight).unwrap(),
                lrcoef_gt_interior_u128(&outer, &inner, &content).unwrap(),
                "shape={shape:?} weight={weight:?} triple=({outer:?}, {inner:?}, {content:?})"
            );
        }
    }

    #[test]
    fn interior_matches_translated_lr_for_all_small_partition_weights() {
        for size in 0..=7 {
            for shape in partitions_of(size) {
                for weight in partitions_of(size) {
                    let (outer, inner, content) = kostka_lr_triple(&shape, &weight).unwrap();
                    assert_eq!(
                        kostka_interior_u128(&shape, &weight).unwrap(),
                        lrcoef_gt_interior_u128(&outer, &inner, &content).unwrap(),
                        "shape={shape:?} weight={weight:?} triple=({outer:?}, {inner:?}, {content:?})"
                    );
                }
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
}
