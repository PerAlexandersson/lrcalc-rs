# lrcalc-rs 0.1.0-rc.1

First release candidate for the Rust implementation of lrcalc, with the
`lrcalc` and `schubmult` commands and shared/static C-compatible libraries.

## Included

- LR, Kostka, Schur and Schubert computations with the existing C ABI.
- Exact affine-hull dimensions and relative-interior constraints for
  ordinary LR, beta-LR and skew Kostka computations.
- Corrected beta-LR column strictness and extreme h-vector evaluation.
- Boundary coverage for packing limits, translated tiny shapes, unused
  content coordinates, LR symmetries and stretched-count identities.
- Packed-state hashing and a dominant-beta Kostka shortcut, with fallback
  to the original LR counter when the packed engine cannot represent a case.
- A consistent unwind strategy for release libraries and Rust tests,
  avoiding abort/unwind output collisions on repeated builds. Panics may
  not unwind across the non-unwinding C ABI boundary.
- Debug/release CI and repeated release-library/test checks.

## Build and check

```sh
cargo build --release --locked
cargo test --locked
cargo test --release --locked
```

This GitHub prerelease supplies source archives; no prebuilt binaries or
crates.io packages are published with it. Build with a current stable Rust
toolchain and C compiler. The Linux shared-library SONAME remains
`liblrcalc.so.2`. See README.md for installation and ABI details, and
notes/BROAD_CORNER_AUDIT.md for the boundary review and repair follow-up.
