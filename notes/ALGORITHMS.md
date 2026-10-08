# Algorithm Notes

This note summarizes the current implementation shape.  It replaces older
exploratory benchmark notes.

## Upstream Model

Classic `lrcalc` separates two tasks.

- Single LR coefficient: compact the triple, then count one coefficient with a
  Buch/Chaput-style tableau backtracker.
- Full Schur expansions: optimize the skew shape, iterate LR tableaux, and
  accumulate an `ivlincomb`.

`lrcalc-rs` follows that split.  The compatibility ABI stays C-shaped, while
the counting engines use native Rust or small C kernels internally.

## Main Engines

- `src/native/lrcoef_count.c`: native Buch-style single LR coefficient kernel
  with upstream-style compactification and pruning.
- `src/lrcoef.rs`: Rust wrapper and beta-prefix LR utilities, including
  full/interior counts and stretch helpers.
- `src/schur.rs`: Schur product, skew Schur expansion, coproduct, fusion, and
  direct ABI accumulation paths.
- `src/kostka_fast.rs`: packed ordinary and skew Kostka dynamic programs,
  including full/interior counts.
- `src/lr_gt.rs`: independent GT/Yamanouchi DP, mainly useful as an oracle,
  for dimensions, and for specialized stretched/count experiments.
- `src/lr_ehrhart.rs`: h-vector interpolation for stretched LR and beta-prefix
  LR families.
- `src/schubert.rs`: Schubert transitions, Monk multiplication, and Schubert
  product expansion.

## Selection Rules

The production path currently uses these principles.

- Use the Buch-style kernel for one ordinary LR coefficient.
- Use direct LR-tableau expansion for skew Schur, Schur product, coproduct, and
  fusion outputs.
- Detect exact Kostka translations in hybrid full/interior count routines and
  dispatch to the packed Kostka DP.
- Use explicit h-vector interpolation for stretched LR computations when the
  caller asks for the stretched/Ehrhart path.
- Keep GT-chain counting as an independent check and as a source of
  interior-count data, not as the default general LR coefficient engine.
- Compute LR, beta-LR, and skew Kostka dimensions and relative-interior masks
  from the exact affine hull of the rational polytope (`src/lr_polytope.rs`,
  `src/affine_hull.rs`).  One homogenized linear program (Freund, Roundy and
  Todd, 1985) identifies every implicit equality, and its answer is checked
  by an exact relative-interior point and an exact dual certificate.  Do not
  infer tightness from the tableaux at a fixed dilation: they need not span
  the polytope.  For example, the five tableaux of
  `c^{(7,6,5,4,2)}_{(6,5,4,2),(3,2,2)}` have no 3 in the third row, but the
  polytope has dimension four.  The Buch, GT-chain, and packed Kostka engines
  share this exact model, so their agreement is not independent evidence;
  regressions compare stretching polynomials with direct counts at fresh
  dilations instead.  Each affine hull costs one exact simplex solve, which
  is cheap compared with tableau enumeration on the tested shapes.

## Current Performance Picture

The current public benchmark table is
[LRCOEF_BENCHMARK_REPORT.md](LRCOEF_BENCHMARK_REPORT.md).  The important
high-level conclusions are:

- ordinary single LR coefficients are near upstream C `lrcalc` parity, often
  faster on setup-heavy zeros and compact shapes;
- large stretched LR examples can be much faster through explicit h-vector
  interpolation than through direct upstream counting;
- exact Kostka translations are best handled by the packed Kostka DP;
- full/interior stretched Kostka-as-LR counts use the hybrid Kostka dispatch;
- Sage wrapper overhead compresses differences, but the Rust library remains
  close to the packaged Sage `liblrcalc` on the measured suite.

## Release Hardening

Remaining compatibility work should stay focused.

- Broaden C ABI smoke tests for iterator, printing, ownership-transfer, and
  unusual empty-input cases.
- Add a full Sage rebuild/install recipe in addition to the current preload
  smoke test.
- Decide whether command-line term order must exactly match upstream hash-table
  iteration, or whether sorted-value equality is the supported guarantee.
- Keep benchmark reports generated from scripts rather than hand-edited timing
  notes.
