# lrcalc-new

`lrcalc-new` is a standalone Rust implementation targeting compatibility with
Anders Buch's `lrcalc` library and command-line tool.

The current crate has a working single Littlewood-Richardson coefficient path,
several independent LR/Kostka counting engines, Ehrhart interpolation
experiments for stretched LR coefficients, and benchmark harnesses against
upstream C.

## Implemented

- `liblrcalc` library target configured as `cdylib`, `staticlib`, and `rlib`.
- C-facing `ivector` allocation/copy/hash/sum functions.
- C ABI export for `schur_lrcoef`.
- `lrcalc coef` / `lrcalc lrcoef` CLI command.
- Native Buch-style single LR coefficient counter with upstream-style
  compactification and branch pruning.
- GT-chain LR counter, including relative-interior and dimension variants.
- Signed Kostka expansion for LR coefficients.
- Fast ordinary and skew Kostka dynamic programs.
- Ehrhart h-vector interpolation for pure stretched LR coefficients.
- Benchmark binaries and shell scripts comparing Rust paths with upstream C.

## Not Yet Implemented

- `ivlincomb` storage, iteration, and ownership-compatible C ABI.
- Schur product, skew Schur expansion, coproduct, fusion, and quantum product
  ABI functions.
- LR tableau iterator ABI.
- Schubert polynomial ABI and `schubmult` compatibility.
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
timeout 60s nice -n 10 cargo run -- coef 3 2 1 - 2 1 - 2 1
timeout 60s nice -n 10 cargo run -- lr-gt 3 2 1 - 2 1 - 2 1
timeout 60s nice -n 10 cargo run -- lr-buch-counts 3 2 1 - 2 1 - 2 1
timeout 60s nice -n 10 cargo run -- kostka-fast 3 2 1 - 2 2 2
timeout 60s nice -n 10 cargo run -- lr-stretch-hvector 3 2 1 - 2 1 - 2 1
timeout 60s nice -n 10 cargo run --release --bin stretched_dp_bench -- 5 3
```

Benchmark scripts live under `scripts/`.  Most expect an upstream `lrcalc`
binary via `UPSTREAM_BIN`; see each script's help text before running.

## Reference

Upstream source:
<https://bitbucket.org/asbuch/lrcalc/src/master/>

Current notes use upstream commit
`8705a16e1575351684ed692f88552f97da3724f4`
from `2025-04-16`.
