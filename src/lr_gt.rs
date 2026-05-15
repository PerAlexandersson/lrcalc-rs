//! Littlewood-Richardson coefficients by a `u128` GT-chain DP.
//!
//! This mirrors the GT/Kostka dynamic program: a tableau of shape
//! `outer / inner` and content `content` is a chain of partitions obtained by
//! adding horizontal strips.  The LR condition is enforced by carrying the
//! prefix sums of the previous strip and bounding the next strip by those
//! prefixes.

use std::collections::{HashMap, HashSet};

use crate::lrcoef::{optim_coef, LrCoefError, OptimizedCoef};
use num_rational::BigRational;
use num_traits::Zero;

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum LrGtError {
    InvalidInput,
    ArithmeticOverflow,
    StateTooWide,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct LrGtStats {
    pub value: u128,
    pub peak_states: usize,
    pub levels: Vec<usize>,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct LrGtInteriorStats {
    pub value: u128,
    pub peak_states: usize,
    pub levels: Vec<usize>,
    pub full_reachable_levels: Vec<usize>,
    pub strict_lower_constraints: usize,
    pub strict_diagonal_constraints: usize,
    pub strict_yamanouchi_constraints: usize,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct LrGtInteriorDfsStats {
    pub value: u128,
    pub weak_nodes: usize,
    pub strict_nodes: usize,
    pub strict_lower_constraints: usize,
    pub strict_diagonal_constraints: usize,
    pub strict_yamanouchi_constraints: usize,
}

#[derive(Clone, Copy, Debug)]
struct Packer {
    bits: u32,
    mask: u128,
    len: usize,
}

type State = (u128, u128);

#[derive(Clone, Debug)]
struct Transition {
    target: State,
    increments: Vec<u32>,
}

#[derive(Clone, Debug)]
struct TightFlags {
    lower: Vec<Vec<bool>>,
    diagonal: Vec<Vec<bool>>,
    yamanouchi: Vec<Vec<bool>>,
}

/// Compute the Littlewood-Richardson coefficient `c^outer_{inner, content}`.
pub fn lrcoef_gt_u128(outer: &[i32], inner: &[i32], content: &[i32]) -> Result<u128, LrGtError> {
    Ok(lrcoef_gt_stats(outer, inner, content)?.value)
}

/// Count relative interior lattice points of the LR GT/Yamanouchi polytope.
pub fn lrcoef_gt_interior_u128(
    outer: &[i32],
    inner: &[i32],
    content: &[i32],
) -> Result<u128, LrGtError> {
    lrcoef_gt_interior_dfs_u128(outer, inner, content)
}

/// Count relative interior points by a strict DFS rather than level DP.
pub fn lrcoef_gt_interior_dfs_u128(
    outer: &[i32],
    inner: &[i32],
    content: &[i32],
) -> Result<u128, LrGtError> {
    Ok(lrcoef_gt_interior_dfs_stats(outer, inner, content)?.value)
}

/// Dimension of the LR GT/Yamanouchi polytope, or `None` if it is empty.
pub fn lrcoef_gt_dimension(
    outer: &[i32],
    inner: &[i32],
    content: &[i32],
) -> Result<Option<usize>, LrGtError> {
    let outer = normalize_partition(outer)?;
    let inner = normalize_partition(inner)?;
    let content = normalize_partition(content)?;

    let outer_size = checked_sum(&outer)?;
    let inner_size = checked_sum(&inner)?;
    let content_size = checked_sum(&content)?;
    if inner_size > outer_size || !partition_less_equal(&inner, &outer) {
        return Ok(None);
    }
    if outer_size - inner_size != content_size {
        return Ok(None);
    }
    if content_size == 0 {
        return Ok(trim_eq(&outer, &inner).then_some(0));
    }
    if outer.is_empty() {
        return Ok(None);
    }

    let rows = outer.len();
    let partition_packer = Packer::new(bit_width(*outer.first().unwrap_or(&0)), rows)?;
    let prefix_packer = Packer::new(bit_width(content_size), rows + 1)?;
    let inner_key = partition_packer.pack_padded(&inner)?;
    let outer_key = partition_packer.pack(&outer)?;
    let start = (inner_key, 0);

    let reachable = reachable_levels(partition_packer, prefix_packer, &outer, &content, start)?;
    let coreachable = coreachable_levels(
        partition_packer,
        prefix_packer,
        &outer,
        &content,
        outer_key,
        &reachable,
    )?;
    if !coreachable[0].contains(&start) {
        return Ok(None);
    }

    let tight = tight_flags_on_complete_paths(
        partition_packer,
        prefix_packer,
        &outer,
        &content,
        &coreachable,
    )?;
    Ok(Some(lr_dimension_from_tight_flags(
        &outer, &content, &tight,
    )))
}

/// Compute `c^outer_{inner, content}` with DP state-count statistics.
pub fn lrcoef_gt_stats(
    outer: &[i32],
    inner: &[i32],
    content: &[i32],
) -> Result<LrGtStats, LrGtError> {
    match optim_coef(outer, inner, content).map_err(map_lrcoef_error)? {
        OptimizedCoef::Zero => Ok(zero_stats()),
        OptimizedCoef::One => Ok(one_stats()),
        OptimizedCoef::Count(shape) => {
            lrcoef_gt_stats_compacted(&shape.outer, &shape.inner, &shape.content)
        }
    }
}

/// Count relative interior lattice points with DP state-count statistics.
pub fn lrcoef_gt_interior_stats(
    outer: &[i32],
    inner: &[i32],
    content: &[i32],
) -> Result<LrGtInteriorStats, LrGtError> {
    match optim_coef(outer, inner, content).map_err(map_lrcoef_error)? {
        OptimizedCoef::Zero => Ok(zero_interior_stats()),
        OptimizedCoef::One => Ok(one_interior_stats()),
        OptimizedCoef::Count(shape) => {
            lrcoef_gt_interior_stats_compacted(&shape.outer, &shape.inner, &shape.content)
        }
    }
}

/// Count relative interior lattice points by first discovering the active
/// weak LR facets on complete paths, then traversing only strict transitions.
pub fn lrcoef_gt_interior_dfs_stats(
    outer: &[i32],
    inner: &[i32],
    content: &[i32],
) -> Result<LrGtInteriorDfsStats, LrGtError> {
    match optim_coef(outer, inner, content).map_err(map_lrcoef_error)? {
        OptimizedCoef::Zero => Ok(zero_interior_dfs_stats()),
        OptimizedCoef::One => Ok(one_interior_dfs_stats()),
        OptimizedCoef::Count(shape) => {
            lrcoef_gt_interior_dfs_stats_compacted(&shape.outer, &shape.inner, &shape.content)
        }
    }
}

fn lrcoef_gt_interior_dfs_stats_compacted(
    outer: &[i32],
    inner: &[i32],
    content: &[i32],
) -> Result<LrGtInteriorDfsStats, LrGtError> {
    let outer = normalize_partition(outer)?;
    let inner = normalize_partition(inner)?;
    let content = normalize_partition(content)?;

    let outer_size = checked_sum(&outer)?;
    let inner_size = checked_sum(&inner)?;
    let content_size = checked_sum(&content)?;
    if inner_size > outer_size || !partition_less_equal(&inner, &outer) {
        return Ok(zero_interior_dfs_stats());
    }
    if outer_size - inner_size != content_size {
        return Ok(zero_interior_dfs_stats());
    }
    if content_size == 0 {
        return Ok(if trim_eq(&outer, &inner) {
            one_interior_dfs_stats()
        } else {
            zero_interior_dfs_stats()
        });
    }
    if outer.is_empty() {
        return Ok(zero_interior_dfs_stats());
    }

    let rows = outer.len();
    let partition_packer = Packer::new(bit_width(*outer.first().unwrap_or(&0)), rows)?;
    let prefix_packer = Packer::new(bit_width(content_size), rows + 1)?;
    let inner_key = partition_packer.pack_padded(&inner)?;
    let outer_key = partition_packer.pack(&outer)?;
    let start = (inner_key, 0);

    let mut tight = TightFlags::new(content.len(), rows);
    let mut weak_nodes = 0usize;
    let has_complete_path = mark_tight_flags_by_weak_dfs(
        partition_packer,
        prefix_packer,
        &outer,
        &content,
        outer_key,
        0,
        start,
        &mut tight,
        &mut weak_nodes,
    )?;
    if !has_complete_path {
        return Ok(LrGtInteriorDfsStats {
            weak_nodes,
            ..zero_interior_dfs_stats()
        });
    }

    let (strict_lower, strict_diagonal, strict_yamanouchi) = tight.strict_counts();
    let mut strict_nodes = 0usize;
    let value = count_relative_interior_by_dfs(
        partition_packer,
        prefix_packer,
        &outer,
        &content,
        outer_key,
        0,
        start,
        &tight,
        &mut strict_nodes,
    )?;

    Ok(LrGtInteriorDfsStats {
        value,
        weak_nodes,
        strict_nodes,
        strict_lower_constraints: strict_lower,
        strict_diagonal_constraints: strict_diagonal,
        strict_yamanouchi_constraints: strict_yamanouchi,
    })
}

fn lrcoef_gt_interior_stats_compacted(
    outer: &[i32],
    inner: &[i32],
    content: &[i32],
) -> Result<LrGtInteriorStats, LrGtError> {
    let outer = normalize_partition(outer)?;
    let inner = normalize_partition(inner)?;
    let content = normalize_partition(content)?;

    let outer_size = checked_sum(&outer)?;
    let inner_size = checked_sum(&inner)?;
    let content_size = checked_sum(&content)?;
    if inner_size > outer_size || !partition_less_equal(&inner, &outer) {
        return Ok(zero_interior_stats());
    }
    if outer_size - inner_size != content_size {
        return Ok(zero_interior_stats());
    }
    if content_size == 0 {
        return Ok(if trim_eq(&outer, &inner) {
            one_interior_stats()
        } else {
            zero_interior_stats()
        });
    }

    let rows = outer.len();
    if rows == 0 {
        return Ok(zero_interior_stats());
    }

    let partition_packer = Packer::new(bit_width(*outer.first().unwrap_or(&0)), rows)?;
    let prefix_packer = Packer::new(bit_width(content_size), rows + 1)?;
    let inner_key = partition_packer.pack_padded(&inner)?;
    let outer_key = partition_packer.pack(&outer)?;
    let start = (inner_key, 0);

    let reachable = reachable_levels(partition_packer, prefix_packer, &outer, &content, start)?;
    let full_reachable_levels = reachable.iter().map(HashSet::len).collect::<Vec<_>>();
    let coreachable = coreachable_levels(
        partition_packer,
        prefix_packer,
        &outer,
        &content,
        outer_key,
        &reachable,
    )?;
    if !coreachable[0].contains(&start) {
        return Ok(LrGtInteriorStats {
            full_reachable_levels,
            ..zero_interior_stats()
        });
    }

    let tight = tight_flags_on_complete_paths(
        partition_packer,
        prefix_packer,
        &outer,
        &content,
        &coreachable,
    )?;
    let (strict_lower, strict_diagonal, strict_yamanouchi) = tight.strict_counts();

    let mut dp: HashMap<State, u128> = HashMap::new();
    dp.insert(start, 1);
    let mut peak_states = dp.len();
    let mut levels = vec![dp.len()];

    for (step, &strip_size) in content.iter().enumerate() {
        let mut next: HashMap<State, u128> = HashMap::with_capacity(dp.len().saturating_mul(2));
        for (&state, &count) in &dp {
            enumerate_relative_interior_transitions(
                partition_packer,
                prefix_packer,
                &outer,
                state.0,
                strip_size,
                previous_prefix_key(step, state),
                step,
                &tight,
                |target| {
                    if !coreachable[step + 1].contains(&target) {
                        return Ok(());
                    }
                    let entry = next.entry(target).or_insert(0);
                    *entry = entry
                        .checked_add(count)
                        .ok_or(LrGtError::ArithmeticOverflow)?;
                    Ok(())
                },
            )?;
        }
        peak_states = peak_states.max(next.len());
        levels.push(next.len());
        dp = next;
    }

    let mut value = 0u128;
    for ((partition_key, _), count) in dp {
        if partition_key == outer_key {
            value = value
                .checked_add(count)
                .ok_or(LrGtError::ArithmeticOverflow)?;
        }
    }

    Ok(LrGtInteriorStats {
        value,
        peak_states,
        levels,
        full_reachable_levels,
        strict_lower_constraints: strict_lower,
        strict_diagonal_constraints: strict_diagonal,
        strict_yamanouchi_constraints: strict_yamanouchi,
    })
}

fn lrcoef_gt_stats_compacted(
    outer: &[i32],
    inner: &[i32],
    content: &[i32],
) -> Result<LrGtStats, LrGtError> {
    let outer = normalize_partition(outer)?;
    let inner = normalize_partition(inner)?;
    let content = normalize_partition(content)?;

    let outer_size = checked_sum(&outer)?;
    let inner_size = checked_sum(&inner)?;
    let content_size = checked_sum(&content)?;
    if inner_size > outer_size || !partition_less_equal(&inner, &outer) {
        return Ok(zero_stats());
    }
    if outer_size - inner_size != content_size {
        return Ok(zero_stats());
    }
    if content_size == 0 {
        return Ok(if trim_eq(&outer, &inner) {
            one_stats()
        } else {
            zero_stats()
        });
    }

    let rows = outer.len();
    if rows == 0 {
        return Ok(zero_stats());
    }

    let partition_packer = Packer::new(bit_width(*outer.first().unwrap_or(&0)), rows)?;
    let prefix_packer = Packer::new(bit_width(content_size), rows + 1)?;
    let inner_key = partition_packer.pack_padded(&inner)?;
    let outer_key = partition_packer.pack(&outer)?;

    let mut dp: HashMap<(u128, u128), u128> = HashMap::new();
    enumerate_initial_extensions(
        partition_packer,
        prefix_packer,
        &outer,
        inner_key,
        content[0],
        |partition_key, prefix_key| {
            let entry = dp.entry((partition_key, prefix_key)).or_insert(0);
            *entry = entry.checked_add(1).ok_or(LrGtError::ArithmeticOverflow)?;
            Ok(())
        },
    )?;

    let mut peak_states = dp.len();
    let mut levels = vec![1, dp.len()];

    for &strip_size in &content[1..] {
        let mut next: HashMap<(u128, u128), u128> =
            HashMap::with_capacity(dp.len().saturating_mul(2));

        for (&(alpha_key, prefix_key), &count) in &dp {
            enumerate_yamanouchi_extensions(
                partition_packer,
                prefix_packer,
                &outer,
                alpha_key,
                strip_size,
                prefix_key,
                |partition_key, new_prefix_key| {
                    let entry = next.entry((partition_key, new_prefix_key)).or_insert(0);
                    *entry = entry
                        .checked_add(count)
                        .ok_or(LrGtError::ArithmeticOverflow)?;
                    Ok(())
                },
            )?;
        }

        peak_states = peak_states.max(next.len());
        levels.push(next.len());
        dp = next;
    }

    let mut value = 0u128;
    for ((partition_key, _), count) in dp {
        if partition_key == outer_key {
            value = value
                .checked_add(count)
                .ok_or(LrGtError::ArithmeticOverflow)?;
        }
    }

    Ok(LrGtStats {
        value,
        peak_states,
        levels,
    })
}

fn map_lrcoef_error(error: LrCoefError) -> LrGtError {
    match error {
        LrCoefError::InvalidPartition => LrGtError::InvalidInput,
        LrCoefError::ArithmeticOverflow => LrGtError::ArithmeticOverflow,
    }
}

fn reachable_levels(
    partition_packer: Packer,
    prefix_packer: Packer,
    outer: &[u32],
    content: &[u32],
    start: State,
) -> Result<Vec<HashSet<State>>, LrGtError> {
    let mut levels = Vec::with_capacity(content.len() + 1);
    let mut current = HashSet::new();
    current.insert(start);
    levels.push(current);

    for (step, &strip_size) in content.iter().enumerate() {
        let mut next = HashSet::new();
        for &state in levels.last().expect("at least one level") {
            enumerate_transition_details(
                partition_packer,
                prefix_packer,
                outer,
                state.0,
                strip_size,
                previous_prefix_key(step, state),
                |transition| {
                    next.insert(transition.target);
                    Ok(())
                },
            )?;
        }
        levels.push(next);
    }

    Ok(levels)
}

fn coreachable_levels(
    partition_packer: Packer,
    prefix_packer: Packer,
    outer: &[u32],
    content: &[u32],
    outer_key: u128,
    reachable: &[HashSet<State>],
) -> Result<Vec<HashSet<State>>, LrGtError> {
    let mut coreachable = vec![HashSet::new(); content.len() + 1];
    for &state in &reachable[content.len()] {
        if state.0 == outer_key {
            coreachable[content.len()].insert(state);
        }
    }

    for step in (0..content.len()).rev() {
        for &state in &reachable[step] {
            let mut reaches_final = false;
            enumerate_transition_details(
                partition_packer,
                prefix_packer,
                outer,
                state.0,
                content[step],
                previous_prefix_key(step, state),
                |transition| {
                    if coreachable[step + 1].contains(&transition.target) {
                        reaches_final = true;
                    }
                    Ok(())
                },
            )?;
            if reaches_final {
                coreachable[step].insert(state);
            }
        }
    }

    Ok(coreachable)
}

fn tight_flags_on_complete_paths(
    partition_packer: Packer,
    prefix_packer: Packer,
    outer: &[u32],
    content: &[u32],
    coreachable: &[HashSet<State>],
) -> Result<TightFlags, LrGtError> {
    let rows = outer.len();
    let mut tight = TightFlags::new(content.len(), rows);

    for (step, &strip_size) in content.iter().enumerate() {
        for &state in &coreachable[step] {
            enumerate_transition_details(
                partition_packer,
                prefix_packer,
                outer,
                state.0,
                strip_size,
                previous_prefix_key(step, state),
                |transition| {
                    if coreachable[step + 1].contains(&transition.target) {
                        tight.observe_transition(
                            partition_packer,
                            prefix_packer,
                            step,
                            state,
                            &transition,
                        );
                    }
                    Ok(())
                },
            )?;
        }
    }

    Ok(tight)
}

fn previous_prefix_key(step: usize, state: State) -> Option<u128> {
    (step > 0).then_some(state.1)
}

#[allow(clippy::too_many_arguments)]
fn enumerate_relative_interior_transitions<F>(
    partition_packer: Packer,
    prefix_packer: Packer,
    outer: &[u32],
    alpha_key: u128,
    strip_size: u32,
    previous_prefix_key: Option<u128>,
    step: usize,
    tight: &TightFlags,
    mut visit: F,
) -> Result<(), LrGtError>
where
    F: FnMut(State) -> Result<(), LrGtError>,
{
    enumerate_relative_interior_transitions_recursive(
        partition_packer,
        prefix_packer,
        outer,
        alpha_key,
        strip_size,
        previous_prefix_key,
        step,
        tight,
        0,
        0,
        alpha_key,
        0,
        &mut visit,
    )
}

#[allow(clippy::too_many_arguments)]
fn enumerate_relative_interior_transitions_recursive<F>(
    partition_packer: Packer,
    prefix_packer: Packer,
    outer: &[u32],
    alpha_key: u128,
    remaining: u32,
    previous_prefix_key: Option<u128>,
    step: usize,
    tight: &TightFlags,
    row: usize,
    prefix_sum: u32,
    partial_partition_key: u128,
    partial_prefix_key: u128,
    visit: &mut F,
) -> Result<(), LrGtError>
where
    F: FnMut(State) -> Result<(), LrGtError>,
{
    if row == outer.len() {
        if remaining == 0 {
            let prefix_key = prefix_packer.set(partial_prefix_key, row, prefix_sum);
            visit((partial_partition_key, prefix_key))?;
        }
        return Ok(());
    }

    let lower_min = if tight.lower[step][row] { 0 } else { 1 };
    if remaining < lower_min {
        return Ok(());
    }

    let base = partition_packer.get(alpha_key, row);
    let max_from_outer = outer[row].saturating_sub(base);
    let max_from_horizontal_strip = if row == 0 {
        remaining
    } else {
        let weak_bound = partition_packer
            .get(alpha_key, row - 1)
            .saturating_sub(base);
        if tight.diagonal[step][row] {
            weak_bound
        } else {
            weak_bound.saturating_sub(1)
        }
    };
    let max_from_yamanouchi = previous_prefix_key.map_or(remaining, |prefix_key| {
        let weak_bound = prefix_packer
            .get(prefix_key, row)
            .saturating_sub(prefix_sum);
        if tight.yamanouchi[step][row] {
            weak_bound
        } else {
            weak_bound.saturating_sub(1)
        }
    });
    let max_increment = remaining
        .min(max_from_outer)
        .min(max_from_horizontal_strip)
        .min(max_from_yamanouchi);
    if max_increment < lower_min {
        return Ok(());
    }

    let next_prefix_key = prefix_packer.set(partial_prefix_key, row, prefix_sum);
    for increment in lower_min..=max_increment {
        let value = base
            .checked_add(increment)
            .ok_or(LrGtError::ArithmeticOverflow)?;
        let next_partition_key = partition_packer.set(partial_partition_key, row, value);
        let next_prefix_sum = prefix_sum
            .checked_add(increment)
            .ok_or(LrGtError::ArithmeticOverflow)?;
        enumerate_relative_interior_transitions_recursive(
            partition_packer,
            prefix_packer,
            outer,
            alpha_key,
            remaining - increment,
            previous_prefix_key,
            step,
            tight,
            row + 1,
            next_prefix_sum,
            next_partition_key,
            next_prefix_key,
            visit,
        )?;
    }

    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn mark_tight_flags_by_weak_dfs(
    partition_packer: Packer,
    prefix_packer: Packer,
    outer: &[u32],
    content: &[u32],
    outer_key: u128,
    step: usize,
    state: State,
    tight: &mut TightFlags,
    weak_nodes: &mut usize,
) -> Result<bool, LrGtError> {
    *weak_nodes = weak_nodes
        .checked_add(1)
        .ok_or(LrGtError::ArithmeticOverflow)?;
    if step == content.len() {
        return Ok(state.0 == outer_key);
    }

    let mut has_completion = false;
    enumerate_transition_details(
        partition_packer,
        prefix_packer,
        outer,
        state.0,
        content[step],
        previous_prefix_key(step, state),
        |transition| {
            if mark_tight_flags_by_weak_dfs(
                partition_packer,
                prefix_packer,
                outer,
                content,
                outer_key,
                step + 1,
                transition.target,
                tight,
                weak_nodes,
            )? {
                tight.observe_transition(partition_packer, prefix_packer, step, state, &transition);
                has_completion = true;
            }
            Ok(())
        },
    )?;
    Ok(has_completion)
}

#[allow(clippy::too_many_arguments)]
fn count_relative_interior_by_dfs(
    partition_packer: Packer,
    prefix_packer: Packer,
    outer: &[u32],
    content: &[u32],
    outer_key: u128,
    step: usize,
    state: State,
    tight: &TightFlags,
    strict_nodes: &mut usize,
) -> Result<u128, LrGtError> {
    *strict_nodes = strict_nodes
        .checked_add(1)
        .ok_or(LrGtError::ArithmeticOverflow)?;
    if step == content.len() {
        return Ok(u128::from(state.0 == outer_key));
    }

    let mut value = 0u128;
    enumerate_relative_interior_transitions(
        partition_packer,
        prefix_packer,
        outer,
        state.0,
        content[step],
        previous_prefix_key(step, state),
        step,
        tight,
        |target| {
            value = value
                .checked_add(count_relative_interior_by_dfs(
                    partition_packer,
                    prefix_packer,
                    outer,
                    content,
                    outer_key,
                    step + 1,
                    target,
                    tight,
                    strict_nodes,
                )?)
                .ok_or(LrGtError::ArithmeticOverflow)?;
            Ok(())
        },
    )?;
    Ok(value)
}

impl TightFlags {
    fn new(steps: usize, rows: usize) -> Self {
        Self {
            lower: vec![vec![true; rows]; steps],
            diagonal: vec![vec![true; rows]; steps],
            yamanouchi: vec![vec![true; rows]; steps],
        }
    }

    fn observe_transition(
        &mut self,
        partition_packer: Packer,
        prefix_packer: Packer,
        step: usize,
        source: State,
        transition: &Transition,
    ) {
        for row in 0..transition.increments.len() {
            if transition.increments[row] > 0 {
                self.lower[step][row] = false;
            }

            if row > 0 {
                let target_value = partition_packer.get(transition.target.0, row);
                let diagonal_bound = partition_packer.get(source.0, row - 1);
                if target_value < diagonal_bound {
                    self.diagonal[step][row] = false;
                }
            }

            if step > 0 {
                let current_prefix = prefix_packer.get(transition.target.1, row + 1);
                let previous_prefix = prefix_packer.get(source.1, row);
                if current_prefix < previous_prefix {
                    self.yamanouchi[step][row] = false;
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

fn lr_dimension_from_tight_flags(outer: &[u32], content: &[u32], tight: &TightFlags) -> usize {
    let rows = outer.len();
    let steps = content.len();
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

    let rank = rational_rank(&equations);
    variable_count.saturating_sub(rank)
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

impl Packer {
    fn new(bits: u32, len: usize) -> Result<Self, LrGtError> {
        let bits = bits.max(1);
        let width = (bits as usize)
            .checked_mul(len)
            .ok_or(LrGtError::ArithmeticOverflow)?;
        if width > 128 {
            return Err(LrGtError::StateTooWide);
        }
        Ok(Self {
            bits,
            mask: (1u128 << bits) - 1,
            len,
        })
    }

    fn get(self, key: u128, index: usize) -> u32 {
        debug_assert!(index < self.len);
        let shift = self.bits * index as u32;
        ((key >> shift) & self.mask) as u32
    }

    fn set(self, key: u128, index: usize, value: u32) -> u128 {
        debug_assert!(index < self.len);
        let shift = self.bits * index as u32;
        let index_mask = self.mask << shift;
        (key & !index_mask) | ((value as u128) << shift)
    }

    fn pack(self, values: &[u32]) -> Result<u128, LrGtError> {
        if values.len() > self.len {
            return Err(LrGtError::InvalidInput);
        }
        let mut key = 0u128;
        for (index, &value) in values.iter().enumerate() {
            if value as u128 > self.mask {
                return Err(LrGtError::StateTooWide);
            }
            key = self.set(key, index, value);
        }
        Ok(key)
    }

    fn pack_padded(self, values: &[u32]) -> Result<u128, LrGtError> {
        if values.len() > self.len && values[self.len..].iter().any(|&part| part != 0) {
            return Err(LrGtError::InvalidInput);
        }
        self.pack(&values[..values.len().min(self.len)])
    }
}

fn enumerate_initial_extensions<F>(
    partition_packer: Packer,
    prefix_packer: Packer,
    outer: &[u32],
    alpha_key: u128,
    strip_size: u32,
    mut visit: F,
) -> Result<(), LrGtError>
where
    F: FnMut(u128, u128) -> Result<(), LrGtError>,
{
    enumerate_extensions(
        partition_packer,
        prefix_packer,
        outer,
        alpha_key,
        strip_size,
        None,
        0,
        0,
        alpha_key,
        0,
        &mut visit,
    )
}

fn enumerate_yamanouchi_extensions<F>(
    partition_packer: Packer,
    prefix_packer: Packer,
    outer: &[u32],
    alpha_key: u128,
    strip_size: u32,
    previous_prefix_key: u128,
    mut visit: F,
) -> Result<(), LrGtError>
where
    F: FnMut(u128, u128) -> Result<(), LrGtError>,
{
    enumerate_extensions(
        partition_packer,
        prefix_packer,
        outer,
        alpha_key,
        strip_size,
        Some(previous_prefix_key),
        0,
        0,
        alpha_key,
        0,
        &mut visit,
    )
}

fn enumerate_transition_details<F>(
    partition_packer: Packer,
    prefix_packer: Packer,
    outer: &[u32],
    alpha_key: u128,
    strip_size: u32,
    previous_prefix_key: Option<u128>,
    mut visit: F,
) -> Result<(), LrGtError>
where
    F: FnMut(Transition) -> Result<(), LrGtError>,
{
    let mut increments = vec![0; outer.len()];
    enumerate_transition_details_recursive(
        partition_packer,
        prefix_packer,
        outer,
        alpha_key,
        strip_size,
        previous_prefix_key,
        0,
        0,
        alpha_key,
        0,
        &mut increments,
        &mut visit,
    )
}

#[allow(clippy::too_many_arguments)]
fn enumerate_transition_details_recursive<F>(
    partition_packer: Packer,
    prefix_packer: Packer,
    outer: &[u32],
    alpha_key: u128,
    remaining: u32,
    previous_prefix_key: Option<u128>,
    row: usize,
    prefix_sum: u32,
    partial_partition_key: u128,
    partial_prefix_key: u128,
    increments: &mut [u32],
    visit: &mut F,
) -> Result<(), LrGtError>
where
    F: FnMut(Transition) -> Result<(), LrGtError>,
{
    if row == outer.len() {
        if remaining == 0 {
            let prefix_key = prefix_packer.set(partial_prefix_key, row, prefix_sum);
            visit(Transition {
                target: (partial_partition_key, prefix_key),
                increments: increments.to_vec(),
            })?;
        }
        return Ok(());
    }

    let base = partition_packer.get(alpha_key, row);
    let max_from_outer = outer[row].saturating_sub(base);
    let max_from_horizontal_strip = if row == 0 {
        remaining
    } else {
        partition_packer
            .get(alpha_key, row - 1)
            .saturating_sub(base)
    };
    let max_from_yamanouchi = previous_prefix_key.map_or(remaining, |prefix_key| {
        prefix_packer
            .get(prefix_key, row)
            .saturating_sub(prefix_sum)
    });
    let max_increment = remaining
        .min(max_from_outer)
        .min(max_from_horizontal_strip)
        .min(max_from_yamanouchi);
    let next_prefix_key = prefix_packer.set(partial_prefix_key, row, prefix_sum);

    for increment in 0..=max_increment {
        increments[row] = increment;
        let value = base
            .checked_add(increment)
            .ok_or(LrGtError::ArithmeticOverflow)?;
        let next_partition_key = partition_packer.set(partial_partition_key, row, value);
        let next_prefix_sum = prefix_sum
            .checked_add(increment)
            .ok_or(LrGtError::ArithmeticOverflow)?;
        enumerate_transition_details_recursive(
            partition_packer,
            prefix_packer,
            outer,
            alpha_key,
            remaining - increment,
            previous_prefix_key,
            row + 1,
            next_prefix_sum,
            next_partition_key,
            next_prefix_key,
            increments,
            visit,
        )?;
    }
    increments[row] = 0;

    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn enumerate_extensions<F>(
    partition_packer: Packer,
    prefix_packer: Packer,
    outer: &[u32],
    alpha_key: u128,
    remaining: u32,
    previous_prefix_key: Option<u128>,
    row: usize,
    prefix_sum: u32,
    partial_partition_key: u128,
    partial_prefix_key: u128,
    visit: &mut F,
) -> Result<(), LrGtError>
where
    F: FnMut(u128, u128) -> Result<(), LrGtError>,
{
    if row == outer.len() {
        if remaining == 0 {
            let prefix_key = prefix_packer.set(partial_prefix_key, row, prefix_sum);
            visit(partial_partition_key, prefix_key)?;
        }
        return Ok(());
    }

    let base = partition_packer.get(alpha_key, row);
    let max_from_outer = outer[row].saturating_sub(base);
    let max_from_horizontal_strip = if row == 0 {
        remaining
    } else {
        partition_packer
            .get(alpha_key, row - 1)
            .saturating_sub(base)
    };
    let max_from_yamanouchi = previous_prefix_key.map_or(remaining, |prefix_key| {
        prefix_packer
            .get(prefix_key, row)
            .saturating_sub(prefix_sum)
    });
    let max_increment = remaining
        .min(max_from_outer)
        .min(max_from_horizontal_strip)
        .min(max_from_yamanouchi);
    let next_prefix_key = prefix_packer.set(partial_prefix_key, row, prefix_sum);

    for increment in 0..=max_increment {
        let value = base
            .checked_add(increment)
            .ok_or(LrGtError::ArithmeticOverflow)?;
        let next_partition_key = partition_packer.set(partial_partition_key, row, value);
        let next_prefix_sum = prefix_sum
            .checked_add(increment)
            .ok_or(LrGtError::ArithmeticOverflow)?;
        enumerate_extensions(
            partition_packer,
            prefix_packer,
            outer,
            alpha_key,
            remaining - increment,
            previous_prefix_key,
            row + 1,
            next_prefix_sum,
            next_partition_key,
            next_prefix_key,
            visit,
        )?;
    }

    Ok(())
}

fn normalize_partition(parts: &[i32]) -> Result<Vec<u32>, LrGtError> {
    if !parts.iter().all(|&part| part >= 0) || !parts.windows(2).all(|w| w[0] >= w[1]) {
        return Err(LrGtError::InvalidInput);
    }
    Ok(parts
        .iter()
        .copied()
        .take_while(|&part| part != 0)
        .map(|part| part as u32)
        .collect())
}

fn checked_sum(parts: &[u32]) -> Result<u32, LrGtError> {
    let mut sum = 0u32;
    for &part in parts {
        sum = sum.checked_add(part).ok_or(LrGtError::ArithmeticOverflow)?;
    }
    Ok(sum)
}

fn partition_less_equal(inner: &[u32], outer: &[u32]) -> bool {
    let len = inner.len().max(outer.len());
    for index in 0..len {
        if *inner.get(index).unwrap_or(&0) > *outer.get(index).unwrap_or(&0) {
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

fn bit_width(value: u32) -> u32 {
    (u32::BITS - value.leading_zeros()).max(1)
}

fn zero_stats() -> LrGtStats {
    LrGtStats {
        value: 0,
        peak_states: 0,
        levels: Vec::new(),
    }
}

fn one_stats() -> LrGtStats {
    LrGtStats {
        value: 1,
        peak_states: 1,
        levels: Vec::new(),
    }
}

fn zero_interior_stats() -> LrGtInteriorStats {
    LrGtInteriorStats {
        value: 0,
        peak_states: 0,
        levels: Vec::new(),
        full_reachable_levels: Vec::new(),
        strict_lower_constraints: 0,
        strict_diagonal_constraints: 0,
        strict_yamanouchi_constraints: 0,
    }
}

fn zero_interior_dfs_stats() -> LrGtInteriorDfsStats {
    LrGtInteriorDfsStats {
        value: 0,
        weak_nodes: 0,
        strict_nodes: 0,
        strict_lower_constraints: 0,
        strict_diagonal_constraints: 0,
        strict_yamanouchi_constraints: 0,
    }
}

fn one_interior_stats() -> LrGtInteriorStats {
    LrGtInteriorStats {
        value: 1,
        peak_states: 1,
        levels: Vec::new(),
        full_reachable_levels: Vec::new(),
        strict_lower_constraints: 0,
        strict_diagonal_constraints: 0,
        strict_yamanouchi_constraints: 0,
    }
}

fn one_interior_dfs_stats() -> LrGtInteriorDfsStats {
    LrGtInteriorDfsStats {
        value: 1,
        weak_nodes: 1,
        strict_nodes: 1,
        strict_lower_constraints: 0,
        strict_diagonal_constraints: 0,
        strict_yamanouchi_constraints: 0,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::lrcoef::lrcoef;

    #[test]
    fn computes_basic_coefficients() {
        assert_eq!(lrcoef_gt_u128(&[], &[], &[]), Ok(1));
        assert_eq!(lrcoef_gt_u128(&[1], &[], &[1]), Ok(1));
        assert_eq!(lrcoef_gt_u128(&[3, 2, 1], &[2, 1], &[2, 1]), Ok(2));
        assert_eq!(lrcoef_gt_u128(&[4, 2], &[2, 1], &[2, 1]), Ok(1));
        assert_eq!(lrcoef_gt_u128(&[5, 1], &[2, 1], &[2, 1]), Ok(0));
    }

    #[test]
    fn agrees_with_lrcoef_port_for_representative_cases() {
        let cases = [
            (&[3, 2, 1][..], &[2, 1][..], &[2, 1][..]),
            (&[4, 3, 2, 1][..], &[3, 2, 1][..], &[1, 1, 1][..]),
            (
                &[7, 6, 5, 4, 3, 2, 1][..],
                &[4, 4, 3, 2, 1][..],
                &[5, 4, 3, 2][..],
            ),
            (&[30, 20, 10][..], &[20, 10][..], &[20, 10][..]),
            (
                &[12, 10, 8, 5, 3, 1][..],
                &[7, 6, 4, 2][..],
                &[9, 8, 5, 3, 1][..],
            ),
        ];
        for (outer, inner, content) in cases {
            assert_eq!(
                lrcoef_gt_u128(outer, inner, content).unwrap(),
                lrcoef(outer, inner, content).unwrap()
            );
        }
    }

    #[test]
    fn agrees_with_lrcoef_port_for_all_small_triples() {
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
                                lrcoef_gt_u128(&outer, &inner, &content).unwrap(),
                                lrcoef(&outer, &inner, &content).unwrap(),
                                "outer={outer:?} inner={inner:?} content={content:?}"
                            );
                        }
                    }
                }
            }
        }
    }

    #[test]
    fn computes_basic_interior_counts() {
        assert_eq!(lrcoef_gt_interior_u128(&[], &[], &[]), Ok(1));
        assert_eq!(lrcoef_gt_interior_u128(&[1], &[], &[1]), Ok(1));
        assert_eq!(lrcoef_gt_interior_u128(&[3, 2, 1], &[2, 1], &[2, 1]), Ok(0));
        assert_eq!(lrcoef_gt_interior_u128(&[4, 2], &[2, 1], &[2, 1]), Ok(1));
        assert_eq!(
            lrcoef_gt_interior_u128(&[5, 4, 3, 2, 1], &[3, 2, 1], &[4, 3, 1, 1]),
            Ok(0)
        );
    }

    #[test]
    fn optimized_interior_agrees_with_raw_path_for_small_triples() {
        for outer_size in 0..=6 {
            for outer in partitions_of(outer_size) {
                for inner_size in 0..=outer_size {
                    for inner in partitions_of(inner_size) {
                        if !partition_less_equal_i32(&inner, &outer) {
                            continue;
                        }
                        let content_size = outer_size - inner_size;
                        for content in partitions_of(content_size) {
                            let optimized =
                                lrcoef_gt_interior_u128(&outer, &inner, &content).unwrap();
                            let raw = lrcoef_gt_interior_stats_compacted(&outer, &inner, &content)
                                .unwrap()
                                .value;
                            assert_eq!(
                                optimized, raw,
                                "outer={outer:?} inner={inner:?} content={content:?}"
                            );
                        }
                    }
                }
            }
        }
    }

    #[test]
    fn strict_dfs_interior_agrees_with_dp_for_small_triples() {
        for outer_size in 0..=6 {
            for outer in partitions_of(outer_size) {
                for inner_size in 0..=outer_size {
                    for inner in partitions_of(inner_size) {
                        if !partition_less_equal_i32(&inner, &outer) {
                            continue;
                        }
                        let content_size = outer_size - inner_size;
                        for content in partitions_of(content_size) {
                            let dp = lrcoef_gt_interior_u128(&outer, &inner, &content).unwrap();
                            let dfs =
                                lrcoef_gt_interior_dfs_u128(&outer, &inner, &content).unwrap();
                            assert_eq!(
                                dfs, dp,
                                "outer={outer:?} inner={inner:?} content={content:?}"
                            );
                        }
                    }
                }
            }
        }
    }

    #[test]
    fn interior_count_is_at_most_full_count_for_small_triples() {
        for outer_size in 0..=6 {
            for outer in partitions_of(outer_size) {
                for inner_size in 0..=outer_size {
                    for inner in partitions_of(inner_size) {
                        if !partition_less_equal_i32(&inner, &outer) {
                            continue;
                        }
                        let content_size = outer_size - inner_size;
                        for content in partitions_of(content_size) {
                            let interior =
                                lrcoef_gt_interior_u128(&outer, &inner, &content).unwrap();
                            let full = lrcoef_gt_u128(&outer, &inner, &content).unwrap();
                            assert!(
                                interior <= full,
                                "outer={outer:?} inner={inner:?} content={content:?} interior={interior} full={full}"
                            );
                        }
                    }
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
