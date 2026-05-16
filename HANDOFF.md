# lrcalc-new Handoff

## Current Goal

Build a standalone Rust replacement for Anders Buch's `lrcalc`, with compatible
C ABI and CLI behavior where possible, while improving native coefficient
engines.

## Current State

- Rust crate builds as `liblrcalc` (`cdylib`, `staticlib`, `rlib`) plus a
  `lrcalc` binary.
- `src/abi.rs` exports the `ivector` allocation/copy/hash/sum subset and
  the core `ivlincomb` allocation/insertion/lookup/iteration/free surface.
  It also exports `schur_lrcoef`, `schur_mult`, and `schur_skew`.
- `src/schur.rs` contains Schur product and skew Schur expansion.  Product
  still enumerates bounded output partitions and reuses scalar `lrcoef`.
  Skew expansion uses the variable-content beta tableau enumerator with
  `beta=[]`, so all output contents are accumulated in one search.
- `src/lrcoef.rs` contains the primary Buch-style LR coefficient engine:
  compactification, pruned tableau search, interior counts, dimension, and
  stretch-cache helpers.  It also exposes beta-prefix LR counts for
  `outer/inner`, `content`, and virtual Yamanouchi prefix `beta`; the beta path
  supports paired full/interior counts and dimension.
- `src/lr_gt.rs` contains an independent GT-chain LR DP with stats, dimension,
  and relative-interior variants.  It also has a hybrid full-count selector
  that tries exact Kostka translation, certified partial Kostka collapse, then
  Buch fallback.  Its production paired full/interior selector tries exact
  Kostka translation, then Buch full/interior counts.
- `src/kostka_fast.rs` contains packed `u128` ordinary/skew Kostka DP and
  interior counts.
- `src/lr_signed.rs` contains the signed Kostka expansion for LR coefficients.
- `src/lr_ehrhart.rs` interpolates h-vectors for pure stretched LR families
  and beta-prefix stretched families.
- CLI commands expose the coefficient engines, beta-prefix LR counts,
  diagnostic/stat modes, and ordinary upstream-style `mult` and `skew`.
  `coprod`, fusion/quantum, `tab`, and `schubmult` are still missing.
- Benchmark scripts compare selected Rust paths against upstream C when an
  upstream binary is available.
- `src/bin/stretched_dp_bench.rs` compares paired full/interior counts for
  stretched Kostka DP and the equivalent stretched LR GT-chain DP.
- `src/bin/lr_hybrid_bench.rs` compares raw GT full and paired counts with the
  combined tableau hybrid selectors on exact and near-Kostka families.
- `src/bin/skew_kostka_ffi_bench.rs` compares direct Rust skew Kostka counts
  with upstream C repeated Schur multiplication.
- `src/bin/partial_collapse_probe.rs` explores row-masked Yamanouchi DPs for
  partial Kostka collapse candidates.

## Verified

`timeout 60s nice -n 10 cargo test` passed on 2026-05-16:
72 library tests, all benchmark-bin test targets, and doc-tests.

Product/skew Schur expansion sanity checks passed on 2026-05-16.  The Rust
CLI agrees with upstream C after sorting output lines for:
`mult 2 1 - 2 1` and `skew 3 2 1 / 2 1`.  CLI line order is not yet treated
as a compatibility guarantee.

Variable-content beta expansion passed on 2026-05-16.  With `beta=[]`, it
matches the Schur expansion coefficients for a representative skew shape.  With
a uniformly dominant finite `beta`, it matches skew Kostka counts over all
weights of a fixed label bound.

Beta-prefix sanity checks passed on 2026-05-16.  With `beta=[]`, beta counts
match ordinary LR full/interior counts on small triples.  With a strictly
dominating partition beta, beta counts match skew Kostka full/interior counts
on representative cases.

Beta-prefix stretch interpolation passed on 2026-05-16.  The beta stretch cache
compacts the base skew diagram, stores tight facets and dimension, then samples
`t*outer/t*inner`, `t*content`, and `t*beta` together.  `beta=[]` specializes to
the ordinary LR stretch polynomial in tests.  Dimension-only stretch commands
are available for LR and beta-prefix LR, avoiding h*-sampling when only the
degree is needed.

`timeout 60s nice -n 10 cargo run --release --bin stretched_dp_bench -- 5 3`
passed on 2026-05-16.  The run showed the equivalent LR GT-chain DP is about
`911x` slower than packed Kostka DP on the scaled Kostka-translation suite.
The hybrid LR path recognized every case as `kostka` mode and ran at about
`1.01x` the packed Kostka time, with matching full and interior counts.

`timeout 60s nice -n 10 cargo run --release --bin lr_hybrid_bench -- 5 3`
passed on 2026-05-16.  The production tableau full-count hybrid was about
`401x` faster than raw GT-chain full counts and about `1.10x` faster than Buch
full counts on the mixed exact/near-Kostka suite.  The paired full/interior
tableau selector was about `152x` faster than raw GT-chain paired counts.

`timeout 60s nice -n 10 scripts/skew_kostka_ffi_bench.sh 1000` passed on
2026-05-16 against `/tmp/lrcalc-upstream/src/.libs/liblrcalc.so`.  Direct Rust
skew Kostka DP took `0.024s`; upstream C repeated Schur multiplication took
`4.332s`, about `180x` slower.

## Main Gaps

- Replace the scalar-loop Schur product implementation with a shared expansion
  iterator/DP when performance matters.
- Implement C ABI functions for Schur coproduct, fusion/quantum, LR tableau
  iteration, and Schubert products.
- Add compatibility headers and C smoke tests for struct layout and exported
  symbols.
- Add Python/Sage rebuild tests against the Rust `liblrcalc`.
- Decide which native LR engine should serve each workload class after broader
  benchmarks.

## Next Useful Work

- Add C/Python smoke tests for `ivlincomb`, `schur_mult`, and `schur_skew`.
- Add an output-order decision for the CLI: either document unordered output or
  mimic upstream hash iteration more closely.
- Add an `nm`-based exported-symbol check for P0/P1 ABI coverage.
- Run `scripts/lrcoef_timed_suite.sh` and the Kostka/LR benchmark scripts
  against a freshly built upstream C binary.
- Add a broader skew Kostka corpus, especially larger sparse skew shapes and
  repeated fixed-shape weights that can reuse `FastSkewKostkaEngine`.
- Profile the Buch-port, GT-chain, and signed-Kostka paths on the same corpus
  before adding new optimizations.
- Benchmark the production tableau paired selector on a broader non-Kostka
  corpus, since it now avoids the GT-chain fallback.
- Use `notes/PARTIAL_KOSTKA_COLLAPSE.md` to guide larger certified partial-mask
  experiments for full counts.
- Extend certified partial collapse to paired full/interior counts only after a
  safe interior certificate is available; for now paired counts use Buch
  fallback outside exact Kostka translations.
- Benchmark beta-prefix stretch interpolation on larger skew Kostka-like
  families, especially cases where empty-row/empty-column compactification is
  substantial.

## Notes

- ABI inventory: `notes/ABI_SYMBOL_INVENTORY_WORKER.md`.
- Upstream ABI summary: `notes/UPSTREAM_ABI.md`.
- Algorithm survey and Rust roadmap: `notes/UPSTREAM_SOURCE_ALGORITHM_SURVEY.md`.
- Partial Kostka-collapse probe: `notes/PARTIAL_KOSTKA_COLLAPSE.md`.
