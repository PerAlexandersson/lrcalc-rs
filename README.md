# lrcalc-rs

`lrcalc-rs` is a standalone Rust implementation targeting compatibility with
Anders Buch's `lrcalc` library and command-line tool.

The project/repository name is `lrcalc-rs`; the installed compatibility surface
keeps the classic names: `liblrcalc.so`, `liblrcalc.a`, `lrcalc`, `schubmult`,
and eventually `include/lrcalc/...`.

The current crate has a working single Littlewood-Richardson coefficient path,
several independent LR/Kostka counting engines, Ehrhart interpolation
experiments for stretched LR coefficients, and benchmark harnesses against
upstream C.

## Implemented

- `liblrcalc` library target configured as `cdylib`, `staticlib`, and `rlib`.
- C-facing `ivector`, `ivlincomb`, `ilist`, `ivlist`, partition iterator,
  partition, permutation/string, LR-tableau iterator, Schur, fusion, and
  Schubert ABI surfaces.  A release `nm` check currently has no missing
  upstream exported symbols; the only extra export is `lrcalc_new_abi_version`.
- `lrcalc coef` / `lrcalc lrcoef`, `lrcalc mult`, `lrcalc skew`,
  `lrcalc coprod`, and `lrcalc tab` CLI commands.  `mult` supports ordinary,
  fusion `-f`, and quantum-printing `-q` products.
- `schubmult` CLI and Schubert ABI exports `trans`, `monk`,
  `mult_poly_schubert`, `mult_schubert`, and `mult_schubert_str`.
- Native Buch-style single LR coefficient counter with upstream-style
  compactification and branch pruning.
- Beta-prefix LR counts for skew shape `outer/inner`, content `content`, and a
  virtual Yamanouchi prefix `beta`, including paired full/interior counts.
- Variable-content beta-prefix expansion, used by skew Schur expansion with
  `beta=[]` and by skew Kostka weight expansion with dominant finite `beta`.
- Upstream-style skew-shape optimization for skew Schur expansion, folding
  forced components into the beta prefix before the shared content expansion.
- Schur product and coproduct expansion via the same shared skew Schur backend,
  using disconnected skew shapes and upstream-style coproduct filtering.
- Fusion product expansion by row-bounded Schur multiplication followed by
  upstream-style affine fusion reduction.
- GT-chain LR counter, including relative-interior and dimension variants.
- Hybrid LR full-count stats that dispatch exact Kostka translations to the
  packed Kostka DP, use certified partial-collapse masks near Kostka shapes,
  and otherwise fall back to Buch's tableau search.
- Hybrid LR full/interior count stats that use packed Kostka on exact
  translations and Buch full/interior counts otherwise.
- Signed Kostka expansion for LR coefficients.
- Fast ordinary and skew Kostka dynamic programs.
- Ehrhart h-vector interpolation for pure stretched LR coefficients and the
  beta-prefix generalization.
- Benchmark binaries and shell scripts comparing Rust paths with upstream C.

## Not Yet Implemented

- Installed compatibility headers.
- Full C smoke tests and Python/Sage rebuild tests against the Rust library.

## Build And Test

```bash
timeout 60s nice -n 10 cargo test
timeout 60s nice -n 10 cargo build --release
```

The Rust library target is named `lrcalc`, so release builds should produce
`liblrcalc.so` and `liblrcalc.a` on Linux.

## Useful Commands

```bash
timeout 60s nice -n 10 cargo run --bin lrcalc -- coef 3 2 1 - 2 1 - 2 1
timeout 60s nice -n 10 cargo run --bin lrcalc -- mult 2 1 - 2 1
timeout 60s nice -n 10 cargo run --bin lrcalc -- mult -f 3,2 2 1 - 2 1
timeout 60s nice -n 10 cargo run --bin lrcalc -- mult -q 3,2 2 1 - 2 1
timeout 60s nice -n 10 cargo run --bin lrcalc -- skew 3 2 1 / 2 1
timeout 60s nice -n 10 cargo run --bin lrcalc -- tab 2 1 / 1
timeout 60s nice -n 10 cargo run --bin schubmult -- 2 1 - 2 1
timeout 60s nice -n 10 cargo run --bin lrcalc -- lr-gt 3 2 1 - 2 1 - 2 1
timeout 60s nice -n 10 cargo run --bin lrcalc -- lr-gt-hybrid-stats 7 4 2 1 - 4 2 - 5 2 1
timeout 60s nice -n 10 cargo run --bin lrcalc -- lr-tableau-hybrid-stats 7 4 2 1 - 4 2 - 5 2 1
timeout 60s nice -n 10 cargo run --bin lrcalc -- lr-tableau-hybrid-counts-stats 7 4 2 1 - 4 2 - 5 2 1
timeout 60s nice -n 10 cargo run --bin lrcalc -- lr-buch-counts 3 2 1 - 2 1 - 2 1
timeout 60s nice -n 10 cargo run --bin lrcalc -- beta-lr-buch-counts 5 3 1 - 3 2 1 - 2 1 - 2 0
timeout 60s nice -n 10 cargo run --bin lrcalc -- kostka-fast 3 2 1 - 2 2 2
timeout 60s nice -n 10 cargo run --bin lrcalc -- skew-kostka-fast 5 3 1 - 3 2 1 - 2 1
timeout 60s nice -n 10 cargo run --bin lrcalc -- lr-stretch-hvector 3 2 1 - 2 1 - 2 1
timeout 60s nice -n 10 cargo run --bin lrcalc -- beta-lr-stretch-dimension 3 3 2 1 1 - 1 1 - 2 2 1 1 1 1 - 4 3 2 1
timeout 60s nice -n 10 cargo run --bin lrcalc -- beta-lr-stretch-hvector 5 3 1 - 3 2 1 - 2 1 - 2 0
timeout 60s nice -n 10 cargo run --release --bin lr_hybrid_bench -- 5 3
timeout 60s nice -n 10 cargo run --release --bin stretched_dp_bench -- 5 3
timeout 60s nice -n 10 cargo run --release --bin partial_collapse_probe
timeout 60s nice -n 10 scripts/lrcoef_ffi_bench.sh 5000 mixed
timeout 60s nice -n 10 scripts/lrcoef_benchmark_report.sh
timeout 60s nice -n 10 scripts/skew_kostka_ffi_bench.sh 1000
timeout 60s nice -n 10 scripts/schur_schubert_ffi_bench.sh 500
```

Current `lr_hybrid_bench -- 5 3` result: the production tableau full-count
hybrid was about `401x` faster than raw GT-chain full counts and about `1.10x`
faster than Buch full counts.  The paired full/interior tableau selector was
about `152x` faster than raw GT-chain paired counts.

Current `stretched_dp_bench -- 5 3` result: pure LR GT-chain counts are about
`911x` slower than packed Kostka counts on the translated suite, while the
hybrid LR path is about `1.01x` the Kostka time and keeps matching full and
interior counts.

Current `skew_kostka_ffi_bench.sh 1000` result: direct Rust skew Kostka DP was
about `180x` faster than an upstream-C repeated Schur multiplication baseline.

Current `lrcoef_ffi_bench.sh 5000 mixed` result: the Rust `lrcoef_i64` ABI path
with native C compactification and counting matches upstream C on 24 mixed LR
cases and takes `0.95x` upstream C time overall.  Most setup-heavy zero and
early-exit cases are now faster than upstream, while count-heavy cases are near
parity.  On `lrcoef_ffi_bench.sh 10000 large-few-parts`, Rust is faster overall
at `0.96x` upstream C.

Current `schur_schubert_ffi_bench.sh 500` result: Schur product, coproduct, and
fusion are near upstream C parity (`1.07x`, `0.80x`, and `1.09x`
Rust/upstream respectively), skew Schur is slower on the included optimized
skew cases (`5.34x`).  On the broader Schubert diagnostic suite, Rust is faster
overall (`0.29x` for permutations and `0.25x` for strings), with the advantage
narrowing from tiny fixed-overhead cases to the larger S6/S7 examples.

Current optimized skew expansion snapshot: on
`skew 20 18 16 14 12 / 10 8 6 4 2`, Rust and upstream C both take about
`0.003s`.  On `skew 30 27 24 21 18 15 / 15 12 9 6 3`, Rust takes about
`1.69s` versus upstream C at about `1.37s`; before the optimizer, packed
accumulator, incremental packed keys, and custom packed-content table the Rust
path took about `7.2s`.

Benchmark scripts live under `scripts/`.  Most expect an upstream `lrcalc`
binary via `UPSTREAM_BIN`; see each script's help text before running.

## Reference

Upstream source:
<https://bitbucket.org/asbuch/lrcalc/src/master/>

Current notes use upstream commit
`8705a16e1575351684ed692f88552f97da3724f4`
from `2025-04-16`.
