//! Littlewood-Richardson coefficients by a `u128` GT-chain DP.
//!
//! This mirrors the GT/Kostka dynamic program: a tableau of shape
//! `outer / inner` and content `content` is a chain of partitions obtained by
//! adding horizontal strips.  The LR condition is enforced by carrying the
//! prefix sums of the previous strip and bounding the next strip by those
//! prefixes.

use std::collections::{HashMap, HashSet};

use crate::kostka_fast::{
    kostka_counts_stats, kostka_fast_stats, KostkaCountsStats, KostkaFastError, KostkaFastStats,
};
use crate::lr_polytope::exact_tight_flags;
use crate::lrcoef::{lrcoef, lrcoef_buch_counts_u128, optim_coef, LrCoefError, OptimizedCoef};

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
pub struct LrGtYamanouchiMaskStats {
    pub value: u128,
    pub peak_states: usize,
    pub levels: Vec<usize>,
    pub enforced_rows: Vec<usize>,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct LrGtYamanouchiMaskCertifiedStats {
    pub stats: LrGtYamanouchiMaskStats,
    pub omitted_constraints_forced: bool,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum LrGtPartialCollapseMode {
    CertifiedMask,
    GtChainFallback,
}

impl LrGtPartialCollapseMode {
    pub const fn label(self) -> &'static str {
        match self {
            Self::CertifiedMask => "certified-mask",
            Self::GtChainFallback => "gt-fallback",
        }
    }
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct LrGtPartialCollapseStats {
    pub mode: LrGtPartialCollapseMode,
    pub enforced_rows: Vec<usize>,
    pub stats: LrGtStats,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum LrGtHybridMode {
    KostkaTranslation,
    CertifiedPartialCollapse,
    GtChainFallback,
}

impl LrGtHybridMode {
    pub const fn label(self) -> &'static str {
        match self {
            Self::KostkaTranslation => "kostka",
            Self::CertifiedPartialCollapse => "certified-mask",
            Self::GtChainFallback => "gt-fallback",
        }
    }
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct LrGtHybridStats {
    pub mode: LrGtHybridMode,
    pub enforced_rows: Vec<usize>,
    pub stats: LrGtStats,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum LrTableauHybridMode {
    KostkaTranslation,
    CertifiedPartialCollapse,
    BuchFallback,
}

impl LrTableauHybridMode {
    pub const fn label(self) -> &'static str {
        match self {
            Self::KostkaTranslation => "kostka",
            Self::CertifiedPartialCollapse => "certified-mask",
            Self::BuchFallback => "buch-fallback",
        }
    }
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct LrTableauHybridStats {
    pub mode: LrTableauHybridMode,
    pub enforced_rows: Vec<usize>,
    pub value: u128,
    pub peak_states: Option<usize>,
    pub levels: Option<Vec<usize>>,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum LrTableauHybridCountsMode {
    KostkaTranslation,
    BuchFallback,
}

impl LrTableauHybridCountsMode {
    pub const fn label(self) -> &'static str {
        match self {
            Self::KostkaTranslation => "kostka",
            Self::BuchFallback => "buch-fallback",
        }
    }
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct LrTableauHybridCountsStats {
    pub mode: LrTableauHybridCountsMode,
    pub full: u128,
    pub interior: u128,
    pub full_peak_states: Option<usize>,
    pub interior_peak_states: Option<usize>,
    pub full_levels: Option<Vec<usize>>,
    pub interior_levels: Option<Vec<usize>>,
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

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct LrGtCountsStats {
    pub full: u128,
    pub interior: u128,
    pub full_peak_states: usize,
    pub full_levels: Vec<usize>,
    pub interior_peak_states: usize,
    pub interior_levels: Vec<usize>,
    pub full_reachable_levels: Vec<usize>,
    pub strict_lower_constraints: usize,
    pub strict_diagonal_constraints: usize,
    pub strict_yamanouchi_constraints: usize,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum LrHybridCountsMode {
    KostkaTranslation,
    GtChain,
}

impl LrHybridCountsMode {
    pub const fn label(self) -> &'static str {
        match self {
            Self::KostkaTranslation => "kostka",
            Self::GtChain => "gt",
        }
    }
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct LrHybridCountsStats {
    pub mode: LrHybridCountsMode,
    pub counts: LrGtCountsStats,
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

#[derive(Clone, Debug)]
struct MaskDpValue {
    count: u128,
    min_prefixes: Vec<u32>,
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

    Ok(exact_gt_tight_flags(&outer, &inner, &content).map(|(dimension, _)| dimension))
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

/// Count with Yamanouchi inequalities enforced only at selected rows.
///
/// This is an exploration hook for partial Kostka collapse.  Passing every row
/// gives the usual LR GT count; passing fewer rows gives a controlled
/// relaxation with a smaller prefix state.
pub fn lrcoef_gt_yamanouchi_mask_stats(
    outer: &[i32],
    inner: &[i32],
    content: &[i32],
    enforced_rows: &[usize],
) -> Result<LrGtYamanouchiMaskStats, LrGtError> {
    lrcoef_gt_yamanouchi_mask_stats_raw(outer, inner, content, enforced_rows)
}

/// Count with a row-masked Yamanouchi DP and report a sufficient certificate.
///
/// If `omitted_constraints_forced` is true, then every omitted Yamanouchi
/// inequality held automatically on all relaxed reachable transitions, so the
/// masked count is the ordinary LR count.
pub fn lrcoef_gt_yamanouchi_mask_certified_stats(
    outer: &[i32],
    inner: &[i32],
    content: &[i32],
    enforced_rows: &[usize],
) -> Result<LrGtYamanouchiMaskCertifiedStats, LrGtError> {
    lrcoef_gt_yamanouchi_mask_certified_stats_raw(outer, inner, content, enforced_rows)
}

/// Candidate rows for partial Kostka collapse.
///
/// Empty output means either an exact Kostka translation, where the direct
/// Kostka fast path is better, or a case where this simple defect heuristic has
/// no useful row to enforce.
pub fn lrcoef_gt_partial_collapse_rows(
    outer: &[i32],
    inner: &[i32],
) -> Result<Vec<usize>, LrGtError> {
    let outer = normalize_partition(outer)?;
    let inner = normalize_partition(inner)?;
    Ok(partial_collapse_rows_normalized(&outer, &inner))
}

/// Safe full-count partial-collapse path.
///
/// The candidate masked DP is used only when the omitted Yamanouchi
/// inequalities are certified as forced.  Otherwise this falls back to the
/// ordinary GT-chain LR count.
pub fn lrcoef_gt_partial_collapse_stats(
    outer: &[i32],
    inner: &[i32],
    content: &[i32],
) -> Result<LrGtPartialCollapseStats, LrGtError> {
    let candidate = lrcoef_gt_partial_collapse_rows(outer, inner)?;
    let certified = lrcoef_gt_yamanouchi_mask_certified_stats(outer, inner, content, &candidate)?;
    if certified.omitted_constraints_forced {
        return Ok(LrGtPartialCollapseStats {
            mode: LrGtPartialCollapseMode::CertifiedMask,
            enforced_rows: candidate,
            stats: LrGtStats {
                value: certified.stats.value,
                peak_states: certified.stats.peak_states,
                levels: certified.stats.levels,
            },
        });
    }

    Ok(LrGtPartialCollapseStats {
        mode: LrGtPartialCollapseMode::GtChainFallback,
        enforced_rows: candidate,
        stats: lrcoef_gt_stats(outer, inner, content)?,
    })
}

/// Full-count hybrid selector for LR coefficients.
///
/// This is still a tableau/GT-chain strategy:
/// exact row-diagonal Kostka translations use the packed Kostka tableau DP,
/// near-Kostka cases may use a certified row-masked LR/GT DP, and all other
/// cases fall back to the ordinary GT-chain LR DP.
pub fn lrcoef_gt_hybrid_stats(
    outer: &[i32],
    inner: &[i32],
    content: &[i32],
) -> Result<LrGtHybridStats, LrGtError> {
    if let Some(weight) = kostka_translation_weight(outer, inner) {
        let stats =
            kostka_fast_to_lr_stats(kostka_fast_stats(content, &weight).map_err(map_kostka_error)?);
        return Ok(LrGtHybridStats {
            mode: LrGtHybridMode::KostkaTranslation,
            enforced_rows: Vec::new(),
            stats,
        });
    }

    let partial = lrcoef_gt_partial_collapse_stats(outer, inner, content)?;
    let mode = match partial.mode {
        LrGtPartialCollapseMode::CertifiedMask => LrGtHybridMode::CertifiedPartialCollapse,
        LrGtPartialCollapseMode::GtChainFallback => LrGtHybridMode::GtChainFallback,
    };
    Ok(LrGtHybridStats {
        mode,
        enforced_rows: partial.enforced_rows,
        stats: partial.stats,
    })
}

/// Production-oriented full-count hybrid selector.
///
/// This keeps the fast tableau DP paths for exact Kostka translations and
/// certified partial collapses, but uses Buch's branch-pruned tableau engine as
/// the general fallback.
pub fn lrcoef_tableau_hybrid_stats(
    outer: &[i32],
    inner: &[i32],
    content: &[i32],
) -> Result<LrTableauHybridStats, LrGtError> {
    if let Some(weight) = kostka_translation_weight(outer, inner) {
        let stats = kostka_fast_stats(content, &weight).map_err(map_kostka_error)?;
        return Ok(LrTableauHybridStats {
            mode: LrTableauHybridMode::KostkaTranslation,
            enforced_rows: Vec::new(),
            value: stats.value,
            peak_states: Some(stats.peak_states),
            levels: Some(stats.levels),
        });
    }

    let candidate = lrcoef_gt_partial_collapse_rows(outer, inner)?;
    if !candidate.is_empty() {
        let certified =
            lrcoef_gt_yamanouchi_mask_certified_stats(outer, inner, content, &candidate)?;
        if certified.omitted_constraints_forced {
            return Ok(LrTableauHybridStats {
                mode: LrTableauHybridMode::CertifiedPartialCollapse,
                enforced_rows: candidate,
                value: certified.stats.value,
                peak_states: Some(certified.stats.peak_states),
                levels: Some(certified.stats.levels),
            });
        }
    }

    Ok(LrTableauHybridStats {
        mode: LrTableauHybridMode::BuchFallback,
        enforced_rows: candidate,
        value: lrcoef(outer, inner, content).map_err(map_lrcoef_error)?,
        peak_states: None,
        levels: None,
    })
}

/// Production-oriented paired full/interior selector.
///
/// Exact row-diagonal Kostka translations use the packed Kostka DP.  General
/// LR shapes use Buch's tableau search, which computes the full count while
/// discovering the tight facets needed for the relative-interior count.
pub fn lrcoef_tableau_hybrid_counts_stats(
    outer: &[i32],
    inner: &[i32],
    content: &[i32],
) -> Result<LrTableauHybridCountsStats, LrGtError> {
    if let Some(weight) = kostka_translation_weight(outer, inner) {
        let stats = kostka_counts_stats(content, &weight).map_err(map_kostka_error)?;
        return Ok(LrTableauHybridCountsStats {
            mode: LrTableauHybridCountsMode::KostkaTranslation,
            full: stats.full,
            interior: stats.interior,
            full_peak_states: Some(stats.full_peak_states),
            interior_peak_states: Some(stats.interior_peak_states),
            full_levels: Some(stats.full_levels),
            interior_levels: Some(stats.interior_levels),
        });
    }

    let counts = lrcoef_buch_counts_u128(outer, inner, content).map_err(map_lrcoef_error)?;
    Ok(LrTableauHybridCountsStats {
        mode: LrTableauHybridCountsMode::BuchFallback,
        full: counts.full,
        interior: counts.interior,
        full_peak_states: None,
        interior_peak_states: None,
        full_levels: None,
        interior_levels: None,
    })
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

/// Compute full and relative-interior LR counts in one GT-chain DP setup.
pub fn lrcoef_gt_counts_stats(
    outer: &[i32],
    inner: &[i32],
    content: &[i32],
) -> Result<LrGtCountsStats, LrGtError> {
    match optim_coef(outer, inner, content).map_err(map_lrcoef_error)? {
        OptimizedCoef::Zero => Ok(zero_counts_stats()),
        OptimizedCoef::One => Ok(one_counts_stats()),
        OptimizedCoef::Count(shape) => {
            lrcoef_gt_counts_stats_compacted(&shape.outer, &shape.inner, &shape.content)
        }
    }
}

/// Compute full and relative-interior LR counts, using the packed Kostka DP
/// when the skew LR shape is the row-diagonal translation of a Kostka problem.
pub fn lrcoef_hybrid_counts_stats(
    outer: &[i32],
    inner: &[i32],
    content: &[i32],
) -> Result<LrHybridCountsStats, LrGtError> {
    if let Some(weight) = kostka_translation_weight(outer, inner) {
        let counts =
            kostka_counts_to_lr(kostka_counts_stats(content, &weight).map_err(map_kostka_error)?);
        return Ok(LrHybridCountsStats {
            mode: LrHybridCountsMode::KostkaTranslation,
            counts,
        });
    }

    Ok(LrHybridCountsStats {
        mode: LrHybridCountsMode::GtChain,
        counts: lrcoef_gt_counts_stats(outer, inner, content)?,
    })
}

fn lrcoef_gt_counts_stats_compacted(
    outer: &[i32],
    inner: &[i32],
    content: &[i32],
) -> Result<LrGtCountsStats, LrGtError> {
    let outer = normalize_partition(outer)?;
    let inner = normalize_partition(inner)?;
    let content = normalize_partition(content)?;

    let outer_size = checked_sum(&outer)?;
    let inner_size = checked_sum(&inner)?;
    let content_size = checked_sum(&content)?;
    if inner_size > outer_size || !partition_less_equal(&inner, &outer) {
        return Ok(zero_counts_stats());
    }
    if outer_size - inner_size != content_size {
        return Ok(zero_counts_stats());
    }
    if content_size == 0 {
        return Ok(if trim_eq(&outer, &inner) {
            one_counts_stats()
        } else {
            zero_counts_stats()
        });
    }

    let rows = outer.len();
    if rows == 0 {
        return Ok(zero_counts_stats());
    }

    let partition_packer = Packer::new(bit_width(*outer.first().unwrap_or(&0)), rows)?;
    let prefix_packer = Packer::new(bit_width(content_size), rows + 1)?;
    let inner_key = partition_packer.pack_padded(&inner)?;
    let outer_key = partition_packer.pack(&outer)?;
    let start = (inner_key, 0);

    let mut dp: HashMap<State, u128> = HashMap::new();
    dp.insert(start, 1);
    let mut reachable = Vec::<HashSet<State>>::with_capacity(content.len() + 1);
    reachable.push(HashSet::from([start]));
    let mut full_peak_states = dp.len();
    let mut full_levels = vec![dp.len()];

    for (step, &strip_size) in content.iter().enumerate() {
        let mut next: HashMap<State, u128> = HashMap::with_capacity(dp.len().saturating_mul(2));
        let mut next_reachable = HashSet::new();
        for (&state, &count) in &dp {
            enumerate_transition_details(
                partition_packer,
                prefix_packer,
                &outer,
                state.0,
                strip_size,
                previous_prefix_key(step, state),
                |transition| {
                    let entry = next.entry(transition.target).or_insert(0);
                    *entry = entry
                        .checked_add(count)
                        .ok_or(LrGtError::ArithmeticOverflow)?;
                    next_reachable.insert(transition.target);
                    Ok(())
                },
            )?;
        }
        full_peak_states = full_peak_states.max(next.len());
        full_levels.push(next.len());
        reachable.push(next_reachable);
        dp = next;
    }

    let mut full = 0u128;
    for ((partition_key, _), count) in &dp {
        if *partition_key == outer_key {
            full = full
                .checked_add(*count)
                .ok_or(LrGtError::ArithmeticOverflow)?;
        }
    }

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
        return Ok(LrGtCountsStats {
            full,
            full_peak_states,
            full_levels,
            full_reachable_levels,
            ..zero_counts_stats()
        });
    }

    let (_, tight) = exact_gt_tight_flags(&outer, &inner, &content)
        .expect("a polytope with a lattice point is nonempty");
    let (strict_lower, strict_diagonal, strict_yamanouchi) = tight.strict_counts();

    let mut strict_dp: HashMap<State, u128> = HashMap::new();
    strict_dp.insert(start, 1);
    let mut interior_peak_states = strict_dp.len();
    let mut interior_levels = vec![strict_dp.len()];

    for (step, &strip_size) in content.iter().enumerate() {
        let mut next: HashMap<State, u128> =
            HashMap::with_capacity(strict_dp.len().saturating_mul(2));
        for (&state, &count) in &strict_dp {
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
        interior_peak_states = interior_peak_states.max(next.len());
        interior_levels.push(next.len());
        strict_dp = next;
    }

    let mut interior = 0u128;
    for ((partition_key, _), count) in strict_dp {
        if partition_key == outer_key {
            interior = interior
                .checked_add(count)
                .ok_or(LrGtError::ArithmeticOverflow)?;
        }
    }

    Ok(LrGtCountsStats {
        full,
        interior,
        full_peak_states,
        full_levels,
        interior_peak_states,
        interior_levels,
        full_reachable_levels,
        strict_lower_constraints: strict_lower,
        strict_diagonal_constraints: strict_diagonal,
        strict_yamanouchi_constraints: strict_yamanouchi,
    })
}

fn lrcoef_gt_yamanouchi_mask_stats_raw(
    outer: &[i32],
    inner: &[i32],
    content: &[i32],
    enforced_rows: &[usize],
) -> Result<LrGtYamanouchiMaskStats, LrGtError> {
    let outer = normalize_partition(outer)?;
    let inner = normalize_partition(inner)?;
    let content = normalize_partition(content)?;

    let outer_size = checked_sum(&outer)?;
    let inner_size = checked_sum(&inner)?;
    let content_size = checked_sum(&content)?;
    let rows = outer.len();
    let enforced_rows = normalize_enforced_rows(enforced_rows, rows)?;

    if inner_size > outer_size || !partition_less_equal(&inner, &outer) {
        return Ok(zero_yamanouchi_mask_stats(enforced_rows));
    }
    if outer_size - inner_size != content_size {
        return Ok(zero_yamanouchi_mask_stats(enforced_rows));
    }
    if content_size == 0 {
        return Ok(if trim_eq(&outer, &inner) {
            one_yamanouchi_mask_stats(enforced_rows)
        } else {
            zero_yamanouchi_mask_stats(enforced_rows)
        });
    }
    if rows == 0 {
        return Ok(zero_yamanouchi_mask_stats(enforced_rows));
    }

    let row_to_slot = row_to_prefix_slot(rows, &enforced_rows);
    let partition_packer = Packer::new(bit_width(*outer.first().unwrap_or(&0)), rows)?;
    let prefix_packer = Packer::new(bit_width(content_size), enforced_rows.len())?;
    let inner_key = partition_packer.pack_padded(&inner)?;
    let outer_key = partition_packer.pack(&outer)?;

    let mut dp: HashMap<State, u128> = HashMap::new();
    dp.insert((inner_key, 0), 1);
    let mut peak_states = dp.len();
    let mut levels = vec![dp.len()];

    for (step, &strip_size) in content.iter().enumerate() {
        let mut next: HashMap<State, u128> = HashMap::with_capacity(dp.len().saturating_mul(2));
        for (&state, &count) in &dp {
            enumerate_masked_yamanouchi_extensions(
                partition_packer,
                prefix_packer,
                &outer,
                &row_to_slot,
                state.0,
                strip_size,
                previous_prefix_key(step, state),
                |target| {
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

    Ok(LrGtYamanouchiMaskStats {
        value,
        peak_states,
        levels,
        enforced_rows,
    })
}

fn lrcoef_gt_yamanouchi_mask_certified_stats_raw(
    outer: &[i32],
    inner: &[i32],
    content: &[i32],
    enforced_rows: &[usize],
) -> Result<LrGtYamanouchiMaskCertifiedStats, LrGtError> {
    let outer = normalize_partition(outer)?;
    let inner = normalize_partition(inner)?;
    let content = normalize_partition(content)?;

    let outer_size = checked_sum(&outer)?;
    let inner_size = checked_sum(&inner)?;
    let content_size = checked_sum(&content)?;
    let rows = outer.len();
    let enforced_rows = normalize_enforced_rows(enforced_rows, rows)?;

    if inner_size > outer_size || !partition_less_equal(&inner, &outer) {
        return Ok(certified_stats(
            zero_yamanouchi_mask_stats(enforced_rows),
            true,
        ));
    }
    if outer_size - inner_size != content_size {
        return Ok(certified_stats(
            zero_yamanouchi_mask_stats(enforced_rows),
            true,
        ));
    }
    if content_size == 0 {
        return Ok(if trim_eq(&outer, &inner) {
            certified_stats(one_yamanouchi_mask_stats(enforced_rows), true)
        } else {
            certified_stats(zero_yamanouchi_mask_stats(enforced_rows), true)
        });
    }
    if rows == 0 {
        return Ok(certified_stats(
            zero_yamanouchi_mask_stats(enforced_rows),
            true,
        ));
    }

    let row_to_slot = row_to_prefix_slot(rows, &enforced_rows);
    let enforced = enforced_row_flags(rows, &enforced_rows);
    let partition_packer = Packer::new(bit_width(*outer.first().unwrap_or(&0)), rows)?;
    let prefix_packer = Packer::new(bit_width(content_size), enforced_rows.len())?;
    let inner_key = partition_packer.pack_padded(&inner)?;
    let outer_key = partition_packer.pack(&outer)?;

    let mut dp: HashMap<State, MaskDpValue> = HashMap::new();
    dp.insert(
        (inner_key, 0),
        MaskDpValue {
            count: 1,
            min_prefixes: vec![0; rows],
        },
    );
    let mut peak_states = dp.len();
    let mut levels = vec![dp.len()];
    let mut omitted_constraints_forced = true;

    for (step, &strip_size) in content.iter().enumerate() {
        let mut next: HashMap<State, MaskDpValue> =
            HashMap::with_capacity(dp.len().saturating_mul(2));
        for (&state, value) in &dp {
            enumerate_masked_yamanouchi_extensions(
                partition_packer,
                prefix_packer,
                &outer,
                &row_to_slot,
                state.0,
                strip_size,
                previous_prefix_key(step, state),
                |target| {
                    let (prefixes_before, prefixes_after) =
                        transition_prefix_bounds(partition_packer, state.0, target.0, rows)?;
                    if step > 0
                        && omitted_constraints_forced
                        && omitted_yamanouchi_can_fail(
                            &enforced,
                            &prefixes_after,
                            &value.min_prefixes,
                        )
                    {
                        omitted_constraints_forced = false;
                    }

                    match next.entry(target) {
                        std::collections::hash_map::Entry::Occupied(mut entry) => {
                            let entry_value = entry.get_mut();
                            entry_value.count = entry_value
                                .count
                                .checked_add(value.count)
                                .ok_or(LrGtError::ArithmeticOverflow)?;
                            minimize_prefixes(&mut entry_value.min_prefixes, &prefixes_before);
                        }
                        std::collections::hash_map::Entry::Vacant(entry) => {
                            entry.insert(MaskDpValue {
                                count: value.count,
                                min_prefixes: prefixes_before,
                            });
                        }
                    }
                    Ok(())
                },
            )?;
        }
        peak_states = peak_states.max(next.len());
        levels.push(next.len());
        dp = next;
    }

    let mut result = 0u128;
    for ((partition_key, _), value) in dp {
        if partition_key == outer_key {
            result = result
                .checked_add(value.count)
                .ok_or(LrGtError::ArithmeticOverflow)?;
        }
    }

    Ok(certified_stats(
        LrGtYamanouchiMaskStats {
            value: result,
            peak_states,
            levels,
            enforced_rows,
        },
        omitted_constraints_forced,
    ))
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

    // The weak search above only certifies a lattice point.  Tightness on
    // the lattice points at one dilation does not determine the affine hull.
    let (_, tight) = exact_gt_tight_flags(&outer, &inner, &content)
        .expect("a polytope with a lattice point is nonempty");
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

    let (_, tight) = exact_gt_tight_flags(&outer, &inner, &content)
        .expect("a polytope with a lattice point is nonempty");
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

fn map_kostka_error(error: KostkaFastError) -> LrGtError {
    match error {
        KostkaFastError::InvalidInput => LrGtError::InvalidInput,
        KostkaFastError::ArithmeticOverflow => LrGtError::ArithmeticOverflow,
        KostkaFastError::StateTooWide => LrGtError::StateTooWide,
    }
}

fn kostka_counts_to_lr(stats: KostkaCountsStats) -> LrGtCountsStats {
    LrGtCountsStats {
        full: stats.full,
        interior: stats.interior,
        full_peak_states: stats.full_peak_states,
        full_levels: stats.full_levels,
        interior_peak_states: stats.interior_peak_states,
        interior_levels: stats.interior_levels,
        full_reachable_levels: stats.full_reachable_levels,
        strict_lower_constraints: stats.strict_lower_constraints,
        strict_diagonal_constraints: stats.strict_diagonal_constraints,
        strict_yamanouchi_constraints: 0,
    }
}

fn kostka_fast_to_lr_stats(stats: KostkaFastStats) -> LrGtStats {
    LrGtStats {
        value: stats.value,
        peak_states: stats.peak_states,
        levels: stats.levels,
    }
}

fn kostka_translation_weight(outer: &[i32], inner: &[i32]) -> Option<Vec<i32>> {
    if !valid_raw_partition(outer) || !valid_raw_partition(inner) {
        return None;
    }

    let outer_len = raw_partition_length(outer);
    if inner.iter().skip(outer_len).any(|&part| part != 0) {
        return None;
    }

    let mut weight = Vec::new();
    for row in 0..outer_len {
        let inner_row = raw_part(outer, row + 1);
        if raw_part(inner, row) != inner_row {
            return None;
        }

        let width = raw_part(outer, row).checked_sub(inner_row)?;
        if width > 0 {
            weight.push(width);
        }
    }

    Some(weight)
}

fn valid_raw_partition(parts: &[i32]) -> bool {
    parts.iter().all(|&part| part >= 0) && parts.windows(2).all(|w| w[0] >= w[1])
}

fn raw_partition_length(parts: &[i32]) -> usize {
    parts
        .iter()
        .rposition(|&part| part != 0)
        .map_or(0, |index| index + 1)
}

fn raw_part(parts: &[i32], index: usize) -> i32 {
    parts.get(index).copied().unwrap_or(0)
}

fn normalize_enforced_rows(enforced_rows: &[usize], rows: usize) -> Result<Vec<usize>, LrGtError> {
    if enforced_rows.iter().any(|&row| row >= rows) {
        return Err(LrGtError::InvalidInput);
    }
    let mut rows = enforced_rows.to_vec();
    rows.sort_unstable();
    rows.dedup();
    Ok(rows)
}

fn row_to_prefix_slot(rows: usize, enforced_rows: &[usize]) -> Vec<Option<usize>> {
    let mut slots = vec![None; rows];
    for (slot, &row) in enforced_rows.iter().enumerate() {
        slots[row] = Some(slot);
    }
    slots
}

fn enforced_row_flags(rows: usize, enforced_rows: &[usize]) -> Vec<bool> {
    let mut flags = vec![false; rows];
    for &row in enforced_rows {
        flags[row] = true;
    }
    flags
}

fn partial_collapse_rows_normalized(outer: &[u32], inner: &[u32]) -> Vec<usize> {
    let rows = outer.len();
    let Some(last_defect) =
        (0..rows).rfind(|&row| part_u32(inner, row) != part_u32(outer, row + 1))
    else {
        return Vec::new();
    };

    if part_u32(inner, last_defect) >= part_u32(outer, last_defect + 1) {
        return Vec::new();
    }

    if last_defect.saturating_mul(2) > rows {
        return Vec::new();
    }

    (0..last_defect).collect()
}

fn part_u32(parts: &[u32], index: usize) -> u32 {
    parts.get(index).copied().unwrap_or(0)
}

fn certified_stats(
    stats: LrGtYamanouchiMaskStats,
    omitted_constraints_forced: bool,
) -> LrGtYamanouchiMaskCertifiedStats {
    LrGtYamanouchiMaskCertifiedStats {
        stats,
        omitted_constraints_forced,
    }
}

fn transition_prefix_bounds(
    partition_packer: Packer,
    source_key: u128,
    target_key: u128,
    rows: usize,
) -> Result<(Vec<u32>, Vec<u32>), LrGtError> {
    let mut before = Vec::with_capacity(rows);
    let mut after = Vec::with_capacity(rows);
    let mut prefix = 0u32;
    for row in 0..rows {
        before.push(prefix);
        let source = partition_packer.get(source_key, row);
        let target = partition_packer.get(target_key, row);
        let increment = target.checked_sub(source).ok_or(LrGtError::InvalidInput)?;
        prefix = prefix
            .checked_add(increment)
            .ok_or(LrGtError::ArithmeticOverflow)?;
        after.push(prefix);
    }
    Ok((before, after))
}

fn omitted_yamanouchi_can_fail(
    enforced: &[bool],
    current_prefixes: &[u32],
    previous_min_prefixes: &[u32],
) -> bool {
    enforced
        .iter()
        .zip(current_prefixes)
        .zip(previous_min_prefixes)
        .any(|((&is_enforced, &current), &previous_min)| !is_enforced && current > previous_min)
}

fn minimize_prefixes(target: &mut [u32], source: &[u32]) {
    for (target, &source) in target.iter_mut().zip(source) {
        *target = (*target).min(source);
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

/// Exact dimension and implicit equalities of the LR GT/Yamanouchi polytope.
///
/// `None` means that the rational polytope is empty.  See
/// [`crate::lr_polytope`] for the inequality system and the certificate.
fn exact_gt_tight_flags(
    outer: &[u32],
    inner: &[u32],
    content: &[u32],
) -> Option<(usize, TightFlags)> {
    let widen = |values: &[u32]| {
        values
            .iter()
            .map(|&value| i64::from(value))
            .collect::<Vec<_>>()
    };
    let exact = exact_tight_flags(
        &widen(outer),
        &widen(inner),
        &widen(content),
        content.len(),
        Some(&[]),
    )?;
    Some((
        exact.dimension,
        TightFlags {
            lower: exact.lower,
            diagonal: exact.diagonal,
            yamanouchi: exact.yamanouchi,
        },
    ))
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

#[allow(clippy::too_many_arguments)]
fn enumerate_masked_yamanouchi_extensions<F>(
    partition_packer: Packer,
    prefix_packer: Packer,
    outer: &[u32],
    row_to_slot: &[Option<usize>],
    alpha_key: u128,
    strip_size: u32,
    previous_prefix_key: Option<u128>,
    mut visit: F,
) -> Result<(), LrGtError>
where
    F: FnMut(State) -> Result<(), LrGtError>,
{
    enumerate_masked_yamanouchi_extensions_recursive(
        partition_packer,
        prefix_packer,
        outer,
        row_to_slot,
        alpha_key,
        strip_size,
        previous_prefix_key,
        0,
        0,
        alpha_key,
        0,
        &mut visit,
    )
}

#[allow(clippy::too_many_arguments)]
fn enumerate_masked_yamanouchi_extensions_recursive<F>(
    partition_packer: Packer,
    prefix_packer: Packer,
    outer: &[u32],
    row_to_slot: &[Option<usize>],
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
    F: FnMut(State) -> Result<(), LrGtError>,
{
    if row == outer.len() {
        if remaining == 0 {
            visit((partial_partition_key, partial_prefix_key))?;
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
    let max_from_yamanouchi = match (previous_prefix_key, row_to_slot[row]) {
        (Some(prefix_key), Some(slot)) => prefix_packer
            .get(prefix_key, slot)
            .saturating_sub(prefix_sum),
        _ => remaining,
    };
    let max_increment = remaining
        .min(max_from_outer)
        .min(max_from_horizontal_strip)
        .min(max_from_yamanouchi);
    let next_prefix_key = if let Some(slot) = row_to_slot[row] {
        prefix_packer.set(partial_prefix_key, slot, prefix_sum)
    } else {
        partial_prefix_key
    };

    for increment in 0..=max_increment {
        let value = base
            .checked_add(increment)
            .ok_or(LrGtError::ArithmeticOverflow)?;
        let next_partition_key = partition_packer.set(partial_partition_key, row, value);
        let next_prefix_sum = prefix_sum
            .checked_add(increment)
            .ok_or(LrGtError::ArithmeticOverflow)?;
        enumerate_masked_yamanouchi_extensions_recursive(
            partition_packer,
            prefix_packer,
            outer,
            row_to_slot,
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

fn zero_yamanouchi_mask_stats(enforced_rows: Vec<usize>) -> LrGtYamanouchiMaskStats {
    LrGtYamanouchiMaskStats {
        value: 0,
        peak_states: 0,
        levels: Vec::new(),
        enforced_rows,
    }
}

fn one_yamanouchi_mask_stats(enforced_rows: Vec<usize>) -> LrGtYamanouchiMaskStats {
    LrGtYamanouchiMaskStats {
        value: 1,
        peak_states: 1,
        levels: vec![1],
        enforced_rows,
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

fn zero_counts_stats() -> LrGtCountsStats {
    LrGtCountsStats {
        full: 0,
        interior: 0,
        full_peak_states: 0,
        full_levels: Vec::new(),
        interior_peak_states: 0,
        interior_levels: Vec::new(),
        full_reachable_levels: Vec::new(),
        strict_lower_constraints: 0,
        strict_diagonal_constraints: 0,
        strict_yamanouchi_constraints: 0,
    }
}

fn one_counts_stats() -> LrGtCountsStats {
    LrGtCountsStats {
        full: 1,
        interior: 1,
        full_peak_states: 1,
        full_levels: vec![1],
        interior_peak_states: 1,
        interior_levels: vec![1],
        full_reachable_levels: vec![1],
        strict_lower_constraints: 0,
        strict_diagonal_constraints: 0,
        strict_yamanouchi_constraints: 0,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::kostka::kostka_lr_triple;
    use crate::kostka_fast::kostka_counts_stats;
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
    fn paired_counts_match_separate_gt_paths() {
        let cases = [
            (&[3, 2, 1][..], &[2, 1][..], &[2, 1][..]),
            (&[4, 2][..], &[2, 1][..], &[2, 1][..]),
            (
                &[7, 6, 5, 4, 3, 2, 1][..],
                &[4, 4, 3, 2, 1][..],
                &[5, 4, 3, 2][..],
            ),
        ];
        for (outer, inner, content) in cases {
            let counts = lrcoef_gt_counts_stats(outer, inner, content).unwrap();
            assert_eq!(counts.full, lrcoef_gt_u128(outer, inner, content).unwrap());
            assert_eq!(
                counts.interior,
                lrcoef_gt_interior_u128(outer, inner, content).unwrap()
            );
        }
    }

    #[test]
    fn hybrid_uses_kostka_translation_when_shape_matches() {
        let shape = [5, 4, 2, 1];
        let weight = [4, 3, 2, 2, 1];
        let (outer, inner, content) = kostka_lr_triple(&shape, &weight).unwrap();

        let kostka = kostka_counts_stats(&shape, &weight).unwrap();
        let hybrid = lrcoef_hybrid_counts_stats(&outer, &inner, &content).unwrap();

        assert_eq!(hybrid.mode, LrHybridCountsMode::KostkaTranslation);
        assert_eq!(hybrid.counts.full, kostka.full);
        assert_eq!(hybrid.counts.interior, kostka.interior);
        assert_eq!(hybrid.counts.full_peak_states, kostka.full_peak_states);
        assert_eq!(
            hybrid.counts.interior_peak_states,
            kostka.interior_peak_states
        );
    }

    #[test]
    fn hybrid_falls_back_for_general_lr_shape() {
        let outer = [4, 2];
        let inner = [2, 1];
        let content = [2, 1];

        let hybrid = lrcoef_hybrid_counts_stats(&outer, &inner, &content).unwrap();
        let gt = lrcoef_gt_counts_stats(&outer, &inner, &content).unwrap();

        assert_eq!(hybrid.mode, LrHybridCountsMode::GtChain);
        assert_eq!(hybrid.counts, gt);
    }

    #[test]
    fn hybrid_matches_gt_counts_for_all_small_triples() {
        for outer_size in 0..=6 {
            for outer in partitions_of(outer_size) {
                for inner_size in 0..=outer_size {
                    for inner in partitions_of(inner_size) {
                        if !partition_less_equal_i32(&inner, &outer) {
                            continue;
                        }
                        let content_size = outer_size - inner_size;
                        for content in partitions_of(content_size) {
                            let hybrid = lrcoef_hybrid_counts_stats(&outer, &inner, &content)
                                .unwrap()
                                .counts;
                            let gt = lrcoef_gt_counts_stats(&outer, &inner, &content).unwrap();
                            assert_eq!(
                                hybrid.full, gt.full,
                                "full mismatch: outer={outer:?} inner={inner:?} content={content:?}"
                            );
                            assert_eq!(
                                hybrid.interior, gt.interior,
                                "interior mismatch: outer={outer:?} inner={inner:?} content={content:?}"
                            );
                        }
                    }
                }
            }
        }
    }

    #[test]
    fn all_row_yamanouchi_mask_matches_gt_count_for_small_triples() {
        for outer_size in 0..=6 {
            for outer in partitions_of(outer_size) {
                let all_rows = (0..outer.len()).collect::<Vec<_>>();
                for inner_size in 0..=outer_size {
                    for inner in partitions_of(inner_size) {
                        if !partition_less_equal_i32(&inner, &outer) {
                            continue;
                        }
                        let content_size = outer_size - inner_size;
                        for content in partitions_of(content_size) {
                            let masked = lrcoef_gt_yamanouchi_mask_stats(
                                &outer, &inner, &content, &all_rows,
                            )
                            .unwrap()
                            .value;
                            let gt = lrcoef_gt_u128(&outer, &inner, &content).unwrap();
                            assert_eq!(
                                masked, gt,
                                "outer={outer:?} inner={inner:?} content={content:?}"
                            );
                        }
                    }
                }
            }
        }
    }

    #[test]
    fn partial_collapse_candidate_rows_match_probe_examples() {
        assert_eq!(
            lrcoef_gt_partial_collapse_rows(&[3, 2, 1], &[2, 1]).unwrap(),
            Vec::<usize>::new()
        );
        assert_eq!(
            lrcoef_gt_partial_collapse_rows(&[4, 2], &[2, 1]).unwrap(),
            Vec::<usize>::new()
        );
        assert_eq!(
            lrcoef_gt_partial_collapse_rows(&[7, 4, 2, 1], &[4, 2]).unwrap(),
            vec![0, 1]
        );
        assert_eq!(
            lrcoef_gt_partial_collapse_rows(&[7, 4, 2, 1], &[4, 3, 1]).unwrap(),
            Vec::<usize>::new()
        );
        assert_eq!(
            lrcoef_gt_partial_collapse_rows(&[7, 6, 5, 4, 3, 2, 1], &[4, 4, 3, 2, 1]).unwrap(),
            Vec::<usize>::new()
        );
    }

    #[test]
    fn certified_partial_collapse_candidate_is_exact_for_small_triples() {
        for outer_size in 0..=6 {
            for outer in partitions_of(outer_size) {
                for inner_size in 0..=outer_size {
                    for inner in partitions_of(inner_size) {
                        if !partition_less_equal_i32(&inner, &outer) {
                            continue;
                        }
                        let candidate = lrcoef_gt_partial_collapse_rows(&outer, &inner).unwrap();
                        let content_size = outer_size - inner_size;
                        for content in partitions_of(content_size) {
                            let certified = lrcoef_gt_yamanouchi_mask_certified_stats(
                                &outer, &inner, &content, &candidate,
                            )
                            .unwrap();
                            if certified.omitted_constraints_forced {
                                assert_eq!(
                                    certified.stats.value,
                                    lrcoef_gt_u128(&outer, &inner, &content).unwrap(),
                                    "outer={outer:?} inner={inner:?} content={content:?} candidate={candidate:?}"
                                );
                            }
                        }
                    }
                }
            }
        }
    }

    #[test]
    fn safe_partial_collapse_path_matches_gt_for_small_triples() {
        for outer_size in 0..=6 {
            for outer in partitions_of(outer_size) {
                for inner_size in 0..=outer_size {
                    for inner in partitions_of(inner_size) {
                        if !partition_less_equal_i32(&inner, &outer) {
                            continue;
                        }
                        let content_size = outer_size - inner_size;
                        for content in partitions_of(content_size) {
                            let partial =
                                lrcoef_gt_partial_collapse_stats(&outer, &inner, &content).unwrap();
                            let gt = lrcoef_gt_u128(&outer, &inner, &content).unwrap();
                            assert_eq!(
                                partial.stats.value, gt,
                                "outer={outer:?} inner={inner:?} content={content:?}"
                            );
                        }
                    }
                }
            }
        }
    }

    #[test]
    fn full_hybrid_uses_kostka_translation_first() {
        let shape = [5, 4, 2, 1];
        let weight = [4, 3, 2, 2, 1];
        let (outer, inner, content) = kostka_lr_triple(&shape, &weight).unwrap();
        let kostka = kostka_fast_stats(&shape, &weight).unwrap();
        let hybrid = lrcoef_gt_hybrid_stats(&outer, &inner, &content).unwrap();

        assert_eq!(hybrid.mode, LrGtHybridMode::KostkaTranslation);
        assert_eq!(hybrid.stats.value, kostka.value);
        assert_eq!(hybrid.stats.peak_states, kostka.peak_states);
        assert!(hybrid.enforced_rows.is_empty());
    }

    #[test]
    fn full_hybrid_uses_certified_partial_collapse() {
        let outer = [7, 4, 2, 1];
        let inner = [4, 2];
        let content = [5, 2, 1];
        let hybrid = lrcoef_gt_hybrid_stats(&outer, &inner, &content).unwrap();
        let gt = lrcoef_gt_u128(&outer, &inner, &content).unwrap();

        assert_eq!(hybrid.mode, LrGtHybridMode::CertifiedPartialCollapse);
        assert_eq!(hybrid.enforced_rows, vec![0, 1]);
        assert_eq!(hybrid.stats.value, gt);
    }

    #[test]
    fn full_hybrid_matches_gt_for_small_triples() {
        for outer_size in 0..=6 {
            for outer in partitions_of(outer_size) {
                for inner_size in 0..=outer_size {
                    for inner in partitions_of(inner_size) {
                        if !partition_less_equal_i32(&inner, &outer) {
                            continue;
                        }
                        let content_size = outer_size - inner_size;
                        for content in partitions_of(content_size) {
                            let hybrid = lrcoef_gt_hybrid_stats(&outer, &inner, &content)
                                .unwrap()
                                .stats
                                .value;
                            let gt = lrcoef_gt_u128(&outer, &inner, &content).unwrap();
                            assert_eq!(
                                hybrid, gt,
                                "outer={outer:?} inner={inner:?} content={content:?}"
                            );
                        }
                    }
                }
            }
        }
    }

    #[test]
    fn tableau_hybrid_uses_buch_fallback_for_noncertified_case() {
        let outer = [7, 4, 2, 1];
        let inner = [4, 3, 1];
        let content = [3, 2, 1];
        let hybrid = lrcoef_tableau_hybrid_stats(&outer, &inner, &content).unwrap();

        assert_eq!(hybrid.mode, LrTableauHybridMode::BuchFallback);
        assert_eq!(hybrid.value, lrcoef(&outer, &inner, &content).unwrap());
        assert_eq!(hybrid.peak_states, None);
    }

    #[test]
    fn tableau_hybrid_matches_buch_for_small_triples() {
        for outer_size in 0..=6 {
            for outer in partitions_of(outer_size) {
                for inner_size in 0..=outer_size {
                    for inner in partitions_of(inner_size) {
                        if !partition_less_equal_i32(&inner, &outer) {
                            continue;
                        }
                        let content_size = outer_size - inner_size;
                        for content in partitions_of(content_size) {
                            let hybrid = lrcoef_tableau_hybrid_stats(&outer, &inner, &content)
                                .unwrap()
                                .value;
                            let buch = lrcoef(&outer, &inner, &content).unwrap();
                            assert_eq!(
                                hybrid, buch,
                                "outer={outer:?} inner={inner:?} content={content:?}"
                            );
                        }
                    }
                }
            }
        }
    }

    #[test]
    fn tableau_counts_hybrid_uses_kostka_translation() {
        let shape = [5, 4, 2, 1];
        let weight = [4, 3, 2, 2, 1];
        let (outer, inner, content) = kostka_lr_triple(&shape, &weight).unwrap();
        let kostka = kostka_counts_stats(&shape, &weight).unwrap();
        let hybrid = lrcoef_tableau_hybrid_counts_stats(&outer, &inner, &content).unwrap();

        assert_eq!(hybrid.mode, LrTableauHybridCountsMode::KostkaTranslation);
        assert_eq!(hybrid.full, kostka.full);
        assert_eq!(hybrid.interior, kostka.interior);
        assert_eq!(hybrid.full_peak_states, Some(kostka.full_peak_states));
        assert_eq!(
            hybrid.interior_peak_states,
            Some(kostka.interior_peak_states)
        );
    }

    #[test]
    fn tableau_counts_hybrid_uses_buch_fallback_for_general_case() {
        let outer = [7, 4, 2, 1];
        let inner = [4, 3, 1];
        let content = [3, 2, 1];
        let hybrid = lrcoef_tableau_hybrid_counts_stats(&outer, &inner, &content).unwrap();
        let buch = lrcoef_buch_counts_u128(&outer, &inner, &content).unwrap();

        assert_eq!(hybrid.mode, LrTableauHybridCountsMode::BuchFallback);
        assert_eq!(hybrid.full, buch.full);
        assert_eq!(hybrid.interior, buch.interior);
        assert_eq!(hybrid.full_peak_states, None);
        assert_eq!(hybrid.interior_peak_states, None);
    }

    #[test]
    fn tableau_counts_hybrid_matches_gt_for_small_triples() {
        for outer_size in 0..=5 {
            for outer in partitions_of(outer_size) {
                for inner_size in 0..=outer_size {
                    for inner in partitions_of(inner_size) {
                        if !partition_less_equal_i32(&inner, &outer) {
                            continue;
                        }
                        let content_size = outer_size - inner_size;
                        for content in partitions_of(content_size) {
                            let hybrid =
                                lrcoef_tableau_hybrid_counts_stats(&outer, &inner, &content)
                                    .unwrap();
                            let gt = lrcoef_gt_counts_stats(&outer, &inner, &content).unwrap();
                            assert_eq!(
                                hybrid.full, gt.full,
                                "full mismatch: outer={outer:?} inner={inner:?} content={content:?}"
                            );
                            assert_eq!(
                                hybrid.interior, gt.interior,
                                "interior mismatch: outer={outer:?} inner={inner:?} content={content:?}"
                            );
                        }
                    }
                }
            }
        }
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
