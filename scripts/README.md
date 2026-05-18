# Scripts

The scripts are grouped by purpose.  Most write build products, staged
prefixes, virtual environments, or benchmark reports under `target/`.

## Release And Smoke Checks

| Script | Purpose |
|---|---|
| `stage_liblrcalc_prefix.sh [PREFIX]` | Build release artifacts and stage headers, libraries, and binaries into a local prefix. |
| `c_abi_smoke.sh` | Compile and run C smoke tests against the staged prefix, including shared and static library checks. |
| `python_bindings_smoke.sh` | Rebuild upstream Python bindings against the staged Rust prefix and run API checks. |
| `sage_bindings_smoke.sh` | Run Sage's `sage.libs.lrcalc` wrapper against the staged Rust library via `LD_PRELOAD`. |

## Benchmark Reports

These scripts compare against an upstream `lrcalc` build.  Override
`UPSTREAM_BIN`, `UPSTREAM_LIB`, or `UPSTREAM_LIB_DIR` when the upstream build is
not in the default location.

| Script | Purpose |
|---|---|
| `lrcoef_benchmark_report.sh [OUT]` | Generate the main markdown LR/Kostka benchmark report. |
| `lrcoef_ffi_bench.sh [ARGS...]` | Run the in-process single-LR-coefficient ABI benchmark. |
| `schur_schubert_ffi_bench.sh [REPEAT]` | Benchmark Schur, fusion, coproduct, skew, and Schubert ABI paths. |
| `skew_kostka_ffi_bench.sh [REPEAT]` | Benchmark direct skew Kostka DP against upstream Schur-expansion paths. |
| `kostka_ffi_bench.sh [REPEAT]` | Benchmark direct Kostka DP against upstream C through a static FFI harness. |
| `lr_signed_ffi_bench.sh [REPEAT] [SUITE]` | Benchmark signed-Kostka LR expansion against upstream C through a static FFI harness. |
| `sage_lrcalc_bench.sh` | Compare Sage's packaged `liblrcalc` with the staged Rust library and write a markdown report. |
| `sage_lrcalc_bench.py` | Python helper used by `sage_lrcalc_bench.sh`; run inside Sage Python. |

## One-Off Diagnostics

| Script | Purpose |
|---|---|
| `kostka_compare_lrcalc.sh SHAPE - WEIGHT` | Translate one Kostka coefficient to LR form and compare Rust with upstream `lrcalc`. |
| `kostka_fast_vs_lrcalc_suite.sh [REPEAT]` | CLI-level benchmark for the fast Kostka path against upstream `lrcalc`. |
| `lrcoef_timed_suite.sh [REPEAT]` | CLI-level benchmark for single LR coefficients against upstream `lrcalc`. |
| `profile_skew_host.sh [options] [-- lrcalc-args...]` | Host-side Samply profiler helper for skew Schur expansion. |

## Upstream Build

Several scripts expect an upstream checkout with either a shared or static
library build:

```bash
git clone https://bitbucket.org/asbuch/lrcalc.git /tmp/lrcalc-upstream
cd /tmp/lrcalc-upstream
autoreconf -i
./configure --enable-shared --disable-static
make -j2
```

For static FFI benchmark scripts, configure upstream with
`--disable-shared --enable-static` instead.
