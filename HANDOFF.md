# lrcalc-rs Handoff

## Current Goal

Build a standalone Rust replacement for Anders Buch's `lrcalc`, with compatible
C ABI and CLI behavior where possible, while improving native coefficient
engines.

## Current State

- Rust crate builds as `liblrcalc` (`cdylib`, `staticlib`, `rlib`) plus a
  `lrcalc` binary and a `schubmult` binary.  Linux release builds now carry
  upstream-compatible SONAME `liblrcalc.so.2`.
- Compatibility headers for the upstream Python/Sage Cython surface and common
  public C headers live under `include/lrcalc/`.
  `scripts/stage_liblrcalc_prefix.sh` stages those headers plus
  `liblrcalc.so`, `liblrcalc.so.2`, `liblrcalc.so.2.0.0`, `liblrcalc.a`,
  `bin/lrcalc`, and `bin/schubmult` under `target/lrcalc-rs-prefix`.
- `src/abi.rs` exports the upstream C ABI symbol surface: `ivector`,
  `ivlincomb`, `ilist`, `ivlist`, partition iterators, partition helpers,
  permutation/string helpers, LR-tableau iterators, Schur/fusion functions,
  and Schubert functions.  A release `nm` diff against upstream currently has
  no missing symbols; the only extra export is `lrcalc_new_abi_version`.
- `src/schur.rs` contains native Schur product, skew Schur, and coproduct
  expansion.  Product expansion realizes `s_mu s_nu` as the skew Schur function
  of a disconnected skew shape, and coproduct expands against the containing
  rectangle with Buch's redundancy filter.  The native Rust skew path still uses
  the variable-content beta tableau enumerator after `optim_skew`.  The C ABI
  `schur_skew` path now special-cases beta-empty skew Schur expansion through
  upstream-style `optim_skew`, direct `lrit_expand`, and a specialized
  `IvLinComb` insertion path.  Fusion products currently use the upstream test
  identity: compute the row-bounded ordinary product, affine-reduce every term,
  and merge signed collisions.
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

`timeout 60s nice -n 10 cargo fmt --check` passed on 2026-05-18.

`timeout 60s nice -n 10 cargo clippy --all-targets -- -D warnings` passed on
2026-05-18.

`timeout 60s nice -n 10 cargo test -q` passed on 2026-05-18:
134 library tests, all benchmark-bin test targets, and doc-tests.

`timeout 60s nice -n 10 scripts/c_abi_smoke.sh` passed on 2026-05-18.  It
stages the Rust install prefix, compiles `tests/c_abi_smoke.c` against the
installed headers with `-Werror`, checks `ivlc_iter` layout against `size_t`
fields, links both shared and static `liblrcalc`, and exercises arbitrary-arity
source-level and exported-symbol initializer calls, low-level containers, Schur
helpers, LR-tableau iteration, Schubert multiplication, `optim_skew`,
`lrcoef_count`, and the staged `lrcalc`/`schubmult` binaries.

`timeout 120s nice -n 10 scripts/python_bindings_smoke.sh` passed on
2026-05-18.  It builds the upstream `python/lrcalc.pyx` Cython module against
the staged Rust `liblrcalc` prefix and checks `lrcoef`, `mult`, `skew`,
`coprod`, `mult_fusion`, Schubert multiplication, string Schubert
multiplication, and LR-tableau iteration from Python.

Sage 10.8 was installed from conda-forge on 2026-05-18:

- Miniforge: `/workspace/.tools/miniforge3`
- Sage environment: `/workspace/.conda-envs/sage`

`timeout 120s nice -n 10 scripts/sage_bindings_smoke.sh` passed on
2026-05-18.  It runs Sage's `sage.libs.lrcalc.lrcalc` wrapper with
`LD_PRELOAD=target/lrcalc-rs-sage-prefix/lib/liblrcalc.so.2`, verifies through
`/proc/self/maps` that the Rust `liblrcalc.so.2.0.0` is loaded, and checks
`lrcoef`, `mult`, `skew`, `coprod`, Schubert multiplication, and LR-tableau
iteration.

`timeout 120s nice -n 10 scripts/sage_lrcalc_bench.sh` passed on
2026-05-18 and wrote the tracked `notes/SAGE_LRCALC_BENCHMARK.md` report.
The script now writes to `target/SAGE_LRCALC_BENCHMARK.md` by default; set
`OUT=notes/SAGE_LRCALC_BENCHMARK.md` when intentionally refreshing the tracked
note.  It compares Sage's
wrapper using conda-forge C `liblrcalc` against the same wrapper using the Rust
`liblrcalc` via `LD_PRELOAD`.  Current geometric mean is `0.965x` Rust/Sage-C,
median is `1.034x`, and correctness signatures match on all 11 cases.

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
The ABI `schur_skew` path now uses optimized shape reduction plus direct
`lrit_expand` content accumulation.  In the current Sage wrapper benchmark,
`skew small` is `1.006x` Rust/Sage-C and `skew optimized` is `1.067x`.  A
larger local 23729-term probe for `30 27 24 21 18 15 / 15 12 9 6 3` had a best
Rust run around `1.40s` versus upstream C around `1.43s`, with run-to-run
noise; treat the dense skew ABI path as near parity, not a solved exact tie.

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
Upstream's `iv_new_init`, `il_new_init`, and `ivl_new_init` are C-variadic.
The exported symbols are now implemented by `src/native/abi_variadic.c`; the
staged C headers also provide source-level static inline variadic constructors,
and the C smoke tests cover both paths with more than eight initializer values.

`ivlc_add_multiple(..., LC_FREE_KEY)` deliberately clears the source table when
moving owned keys into a distinct destination.  This avoids the upstream
dangling-key/double-free hazard, but differs if downstream code observes
`ivlc_card(src)` after a transfer.

`scripts/schur_schubert_ffi_bench.sh 500` was added and passed on 2026-05-17
against `/tmp/lrcalc-upstream/src/.libs/liblrcalc.so`.  It verifies and times
Schur product, skew Schur, coproduct, fusion product, Schubert permutation
products, and Schubert string products in-process through upstream C FFI.  After
the direct `lrit_expand` ABI path, recent raw skew suite runs are much closer
than the old `5.34x` baseline but still noisy, roughly `1.3x`--`1.5x`
Rust/upstream C on selected optimized cases.  The Sage benchmark rows are
closer (`1.006x` and `1.067x`).  Schubert remains clearly faster in Rust in
both raw diagnostics and Sage wrapper timing.

Basic LR coefficient counting now uses a compact count-only tableau box with
32-bit indices, matching upstream's 32-byte `lrcoef_box` shape more closely.
The row-aware box remains for tableau/interior paths.  On
`lr_gt_vs_buch_bench -- 50000 mixed`, this changed Buch full-count time from
`0.128s` at commit `3107ad1` to `0.114s` after commit `af87dd4`.

`scripts/lrcoef_ffi_bench.sh 5000 mixed` compares the Rust `lrcoef_i64` ABI path
directly against upstream C `schur_lrcoef` through `dlopen`.  After adding the
native C compactification/count path, it passed on 2026-05-17 against
`/tmp/lrcalc-upstream/src/.libs/liblrcalc.so`; current overall Rust/upstream C
ratio is `0.95x` on 24 mixed cases.  Most setup-heavy zero/early-exit cases are
now faster than upstream, while count-heavy cases are near parity.  On
`scripts/lrcoef_ffi_bench.sh 10000 large-few-parts`, Rust is faster overall at
`0.96x` upstream C.

The native compactification in `src/native/lrcoef_count.c` intentionally mirrors
Buch's GPL-3.0-or-later `optshape.c` coefficient path for ABI performance, while
the public Rust path keeps a Rust fallback for errors and overflow.

## Main Gaps

- Broaden the skew benchmark corpus.  The ABI path now uses the same direct
  `lrit_expand` shape as upstream; remaining skew work is mostly small/medium
  constant factors and Rust ABI table overhead.
- Optimize fusion products by using the newly exposed `optim_fusion` path in
  `schur_mult_fusion` rather than reducing a full row-bounded product.
- Broaden C smoke tests for the ABI surface, including more iterator and
  printing cases.
- Audit staged headers against the complete upstream installed-header surface,
  especially allocator/template-header expectations.
- Expand Sage tests beyond the current `LD_PRELOAD` smoke test, or rebuild
  Sage's lrcalc package against the Rust install prefix.
- Decide whether native Rust skew expansion should keep the beta enumerator as
  its default or route beta-empty cases through the ABI-style lrit accumulator.

## Next Useful Work

- Broaden C/Python smoke tests for `ivlincomb`, Schur products, Schubert
  products, printing helpers, and low-level ownership-transfer APIs.
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
