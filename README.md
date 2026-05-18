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

## Upstream, Credits, and Ecosystem

`lrcalc-rs` is a compatibility-oriented Rust implementation of the classic
[Littlewood-Richardson Calculator](https://sites.math.rutgers.edu/~asbuch/lrcalc/)
by Anders S. Buch.  The upstream source repository is
[asbuch/lrcalc on Bitbucket](https://bitbucket.org/asbuch/lrcalc/src/master/).

The upstream project and ecosystem deserve explicit credit.  The original
`lrcalc` package is by Anders S. Buch; its upstream page also credits
Nicolas M. Thiery and Jean-Pierre Flori for the GNU automake system, and
Pierre-Emmanuel Chaput for suggestions that led to a large speedup in single
LR coefficient computation.  Sage's `lrcalc` interface credits Mike Hansen for
the core interface, and Anne Schilling, Nicolas M. Thiery, and Anders Buch for
fusion products, LR-tableau iteration, finalization, and documentation.

Classic `lrcalc` is used by several downstream systems and package ecosystems:

- [SageMath package `lrcalc`](https://doc.sagemath.org/html/en/reference/spkg/lrcalc.html)
- [SageMath `sage.libs.lrcalc` interface](https://doc.sagemath.org/html/en/reference/libs/sage/libs/lrcalc/lrcalc.html)
- [Python bindings on PyPI](https://pypi.org/project/lrcalc/)
- [Repology package overview](https://repology.org/project/lrcalc/versions)

## Mathematical Background

For background on the main objects computed by this project, see SymCat:

- [Schur polynomials](https://www.symmetricfunctions.com/schur.htm)
- [Littlewood--Richardson coefficients](https://www.symmetricfunctions.com/littlewoodRichardson.htm)
- [Kostka coefficients and Kostka--Foulkes polynomials](https://www.symmetricfunctions.com/kostkaFoulkes.htm)
- [Schubert polynomials](https://www.symmetricfunctions.com/schubert.htm)

## Support

To support SymCat and related symmetric-functions resources, see
[Ko-fi: symmetricfunctions](https://ko-fi.com/symmetricfunctions).

## Implemented

- `liblrcalc` library target configured as `cdylib`, `staticlib`, and `rlib`.
- Linux release builds use upstream-compatible SONAME `liblrcalc.so.2`.
- Compatibility headers for the upstream Python/Sage Cython surface under
  `include/lrcalc/`.
- Initial C ABI smoke coverage for staged headers, shared and static library
  links, `ivlc_iter` layout, low-level containers, Schur/LR-tableau helpers,
  Schubert helpers, staged CLI binaries, and the native LR coefficient kernel.
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
- Variable-content beta-prefix expansion, used by native Rust skew/content
  experiments and by skew Kostka weight expansion with dominant finite `beta`.
- Upstream-style skew-shape optimization for skew Schur expansion.  The C ABI
  `schur_skew` path uses optimized shape reduction plus direct `lrit_expand`
  accumulation for the beta-empty case.
- Schur product and coproduct expansion via skew Schur backends, using
  disconnected skew shapes and upstream-style coproduct filtering.
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

- Broader C layout/link smoke tests against the full upstream header inventory.
- Sage rebuild tests against the Rust library; the current Sage check is an
  `LD_PRELOAD` wrapper smoke test.
- Complete installed-header audit against the full upstream header inventory.

## License

`lrcalc-rs` is distributed under the GNU General Public License, version 3 or
any later version (`GPL-3.0-or-later`).  This matches the upstream
Littlewood-Richardson Calculator license and is compatible with SageMath's GPL
distribution model.

See [LICENSE](LICENSE) for the project notice and attribution notes, and
[COPYING](COPYING) for the full GPLv3 text.

## Repository Status

For contribution workflow notes, see [CONTRIBUTING.md](CONTRIBUTING.md).
For the public-repository and drop-in release checklist, see
[PUBLISHING.md](PUBLISHING.md).

## Build And Test

```bash
cargo test
cargo clippy --all-targets -- -D warnings
cargo build --release
scripts/c_abi_smoke.sh
scripts/python_bindings_smoke.sh
scripts/sage_bindings_smoke.sh
scripts/sage_lrcalc_bench.sh
```

The Rust library target is named `lrcalc`, so release builds should produce
`liblrcalc.so` and `liblrcalc.a` on Linux.  The staging helper
`scripts/stage_liblrcalc_prefix.sh` creates a local prefix under
`target/lrcalc-rs-prefix` with `include/lrcalc/`, `liblrcalc.a`,
`liblrcalc.so`, `liblrcalc.so.2`, `liblrcalc.so.2.0.0`, `bin/lrcalc`, and
`bin/schubmult`.

The Sage scripts use `sage -python` when `sage` is on `PATH`.  Alternatively,
set `SAGE_PYTHON=/path/to/sage/python`.

## Useful Commands

```bash
cargo run --bin lrcalc -- coef 3 2 1 - 2 1 - 2 1
cargo run --bin lrcalc -- mult 2 1 - 2 1
cargo run --bin lrcalc -- mult -f 3,2 2 1 - 2 1
cargo run --bin lrcalc -- mult -q 3,2 2 1 - 2 1
cargo run --bin lrcalc -- skew 3 2 1 / 2 1
cargo run --bin lrcalc -- tab 2 1 / 1
cargo run --bin schubmult -- 2 1 - 2 1
cargo run --bin lrcalc -- lr-gt 3 2 1 - 2 1 - 2 1
cargo run --bin lrcalc -- lr-gt-hybrid-stats 7 4 2 1 - 4 2 - 5 2 1
cargo run --bin lrcalc -- lr-tableau-hybrid-stats 7 4 2 1 - 4 2 - 5 2 1
cargo run --bin lrcalc -- lr-tableau-hybrid-counts-stats 7 4 2 1 - 4 2 - 5 2 1
cargo run --bin lrcalc -- lr-buch-counts 3 2 1 - 2 1 - 2 1
cargo run --bin lrcalc -- beta-lr-buch-counts 5 3 1 - 3 2 1 - 2 1 - 2 0
cargo run --bin lrcalc -- kostka-fast 3 2 1 - 2 2 2
cargo run --bin lrcalc -- skew-kostka-fast 5 3 1 - 3 2 1 - 2 1
cargo run --bin lrcalc -- lr-stretch-hvector 3 2 1 - 2 1 - 2 1
cargo run --bin lrcalc -- beta-lr-stretch-dimension 3 3 2 1 1 - 1 1 - 2 2 1 1 1 1 - 4 3 2 1
cargo run --bin lrcalc -- beta-lr-stretch-hvector 5 3 1 - 3 2 1 - 2 1 - 2 0
cargo run --release --bin lr_hybrid_bench -- 5 3
cargo run --release --bin stretched_dp_bench -- 5 3
cargo run --release --bin partial_collapse_probe
scripts/lrcoef_ffi_bench.sh 5000 mixed
scripts/lrcoef_benchmark_report.sh
scripts/skew_kostka_ffi_bench.sh 1000
scripts/schur_schubert_ffi_bench.sh 500
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
fusion are near upstream C parity.  After the direct `lrit_expand` ABI path,
raw in-process skew Schur runs are much closer than the old `5.34x` baseline,
with recent selected runs around `1.3x`--`1.5x` Rust/upstream C.  On the broader
Schubert diagnostic suite, Rust is faster overall, with the advantage narrowing
from tiny fixed-overhead cases to the larger S6/S7 examples.

Current Sage wrapper benchmark result: geometric mean `0.965x` Rust/Sage-C
across 11 cases, with matching correctness signatures.  The skew rows are near
parity (`1.006x` and `1.067x`), and Schubert multiplication is faster in Rust
through Sage (`0.445x`).

Benchmark scripts live under `scripts/`.  CLI-oracle scripts usually take
`UPSTREAM_BIN`; FFI benchmark scripts usually take `UPSTREAM_LIB` or
`UPSTREAM_LIB_DIR`.  Set `OUT=notes/SAGE_LRCALC_BENCHMARK.md` when intentionally
refreshing the tracked Sage benchmark note.

## Reference

Upstream source:
<https://bitbucket.org/asbuch/lrcalc/src/master/>

Current notes use upstream commit
`8705a16e1575351684ed692f88552f97da3724f4`
from `2025-04-16`.
