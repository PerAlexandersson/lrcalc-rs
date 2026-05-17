# lrcalc-new Handoff

## Current Goal

Build a standalone Rust replacement for Anders Buch's `lrcalc`, with compatible
C ABI and CLI behavior where possible, while improving native coefficient
engines.

## Current State

- Rust crate builds as `liblrcalc` (`cdylib`, `staticlib`, `rlib`) plus a
  `lrcalc` binary and a `schubmult` binary.
- `src/abi.rs` exports the upstream C ABI symbol surface: `ivector`,
  `ivlincomb`, `ilist`, `ivlist`, partition iterators, partition helpers,
  permutation/string helpers, LR-tableau iterators, Schur/fusion functions,
  and Schubert functions.  A release `nm` diff against upstream currently has
  no missing symbols; the only extra export is `lrcalc_new_abi_version`.
- `src/schur.rs` contains Schur product and skew Schur expansion.  Product
  expansion now realizes `s_mu s_nu` as the skew Schur function of a
  disconnected skew shape, so it shares the skew expansion backend instead of
  looping over scalar LR coefficients.  Coproduct expands against the
  containing rectangle and applies Buch's redundancy filter.  Skew expansion
  first applies an upstream-style `optim_skew` shape reduction, folds forced
  components into a beta prefix, and then uses the variable-content beta
  tableau enumerator so all residual contents are accumulated in one search.
  Fusion products currently use the upstream test identity: compute the
  row-bounded ordinary product, affine-reduce every term, and merge signed
  collisions.
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
- `src/schubert.rs` ports upstream `schublib.c` recursion for Schubert
  transition polynomials, Monk multiplication, permutation products, and string
  products.
- CLI commands expose the coefficient engines, beta-prefix LR counts,
  diagnostic/stat modes, and ordinary upstream-style `mult`, `skew`, and
  `coprod`.  `mult -f rows,level` and `mult -q rows,level` are implemented.
  CLI `tab` is implemented.  Standalone `schubmult` is implemented.
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

`timeout 60s nice -n 10 cargo test` passed on 2026-05-17:
124 library tests, all benchmark-bin test targets, and doc-tests.

Product/skew Schur expansion sanity checks passed on 2026-05-16.  The Rust
CLI agrees with upstream C after sorting output lines for:
`mult 2 1 - 2 1` and `skew 3 2 1 / 2 1`.  CLI line order is not yet treated
as a compatibility guarantee.

Variable-content beta expansion passed on 2026-05-16.  With `beta=[]`, it
matches the Schur expansion coefficients for a representative skew shape.  With
a uniformly dominant finite `beta`, it matches skew Kostka counts over all
weights of a fixed label bound.

Optimized skew expansion passed on 2026-05-16.  The regression test exhausts
all contained skew shapes of outer size at most `7` and row bounds
`[-1, 0, 1, 2, 3, 4]`, comparing the optimized path with scalar LR expansion.
On `skew 20 18 16 14 12 / 10 8 6 4 2`, Rust and upstream C both took about
`0.003s`.  On `skew 30 27 24 21 18 15 / 15 12 9 6 3`, Rust improved from
about `7.2s` to `4.1s` after `optim_skew`, then to `3.18s` after packed
content accumulation, then to `1.85s` after maintaining packed keys
incrementally during tableau search, then to `1.69s` after replacing the packed
content `HashMap<u128, u128>` with a specialized open-addressing table;
upstream C took about `1.37s`.

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

Schur product now uses disconnected skew Schur expansion.  On 2026-05-17,
`mult 8 6 4 2 - 8 6 4 2` ran about `1.75x` faster than the previous scalar
product loop, and `mult 12 9 6 3 - 12 9 6 3` ran about `2.60x` faster.

Fusion product support was added on 2026-05-17.  Focused Rust and ABI tests pass
for `fusion_reduce`, `fusion_reduce_lc`, and `schur_mult_fusion`.  CLI checks
against upstream passed for `mult -f 2,1 1 - 1`,
`mult -f 3,2 2 1 - 2 1`, and Maple quantum output for
`mult -m -q 3,2 2 1 - 2 1`, modulo output order.  A small sorted-output
comparison also matched upstream on eight fusion and six quantum examples.
`timeout 60s nice -n 10 cargo build --release` passed, and `nm -D` shows the
three fusion symbols plus `part_qdegree` and `part_qentry` exported from
`target/release/liblrcalc.so`.

The LR tableau iterator ABI (`lrit_new`, `lrit_good`, `lrit_next`,
`lrit_free`) was added on 2026-05-17.  Unit tests cover the small skew shape
`(2,1)/(1)`, row-bounded iteration, and an empty non-contained skew shape.
`nm -D` shows all four `lrit_*` symbols exported from the release shared
library.

CLI `tab` was added on 2026-05-17 on top of the same `lrit_*` iterator.
Output matched upstream exactly on five small cases, including row-bound and
weight-filtered examples.

Schubert support was added on 2026-05-17.  Focused Rust and ABI tests pass for
`trans`, `monk`, `mult_poly_schubert`, `mult_schubert`, and
`mult_schubert_str`.  The standalone `schubmult` binary matched upstream after
sorting output lines on all 576 `S_4` permutation products and 20 compatible
binary-string products of length 3.  `timeout 60s nice -n 10 cargo build
--release` passed, and `nm -D` shows the five Schubert symbols exported from
`target/release/liblrcalc.so`.

The remaining upstream-exported ABI helpers were added on 2026-05-17:
`iv_*`, `il_*`, `ivl_*`, `ivlc_equals/print/print_stat`, `part_*`,
`perm_*`, `pitr_*`, `lrit_count/expand/print/dump`, `lrcoef_count`,
`maple_*`, `optim_*`, `sksh_*`, `all_strings`, `all_perms`, and
`get_vect_arg`.  `timeout 60s nice -n 10 cargo test` and `timeout 60s nice -n
10 cargo build --release` pass.  A release `nm -D` comparison against upstream
shows no missing exported symbols and only the intentional extra
`lrcalc_new_abi_version`.
Upstream's `iv_new_init`, `il_new_init`, and `ivl_new_init` are C-variadic;
stable Rust cannot define true variadic exports, so the current symbols are
fixed-argument compatibility shims covering the first eight initializer values.

`scripts/schur_schubert_ffi_bench.sh 500` was added and passed on 2026-05-17
against `/tmp/lrcalc-upstream/src/.libs/liblrcalc.so`.  It verifies and times
Schur product, skew Schur, coproduct, fusion product, Schubert permutation
products, and Schubert string products in-process through upstream C FFI.
Current ratios, reported as Rust/upstream C, were: product `1.07x`, skew
`5.34x`, coproduct `0.80x`, fusion `1.09x`, Schubert permutations `0.29x`,
and Schubert strings `0.25x`.  The Schubert permutation diagnostic now includes
18 cases up to a selected S7 product; Rust ranges from about `0.04x`--`0.10x`
on tiny fixed-overhead products to `0.32x` on the S7 medium case.  This makes
skew Schur the main remaining performance gap in the Schur/Schubert surface.

## Main Gaps

- Continue low-level skew expansion tuning.  The high-level algorithm now
  matches upstream more closely, but dense cases still trail upstream C because
  the Rust path lacks the exact tight `lrit_next`-style iterator and packed
  output accumulator.
- Optimize fusion products by using the newly exposed `optim_fusion` path in
  `schur_mult_fusion` rather than reducing a full row-bounded product.
- Add C smoke tests for the ABI surface, including Schubert and `schubmult`.
- Add compatibility headers and C smoke tests for struct layout.
- Add Python/Sage rebuild tests against the Rust `liblrcalc`.
- Decide which native LR engine should serve each workload class after broader
  benchmarks.

## Next Useful Work

- Add C/Python smoke tests for `ivlincomb`, `schur_mult`, `schur_skew`,
  `schur_coprod`, Schubert products, and the low-level helper APIs.
- Add an output-order decision for the CLI: either document unordered output or
  mimic upstream hash iteration more closely.
- Add an `nm`-based exported-symbol check to CI so full ABI coverage stays
  visible.
- Run `scripts/lrcoef_timed_suite.sh` and the Kostka/LR benchmark scripts
  against a freshly built upstream C binary.
- Add a broader skew expansion corpus, especially shapes where `optim_skew`
  removes large forced components and dense shapes where accumulator overhead
  dominates.
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
