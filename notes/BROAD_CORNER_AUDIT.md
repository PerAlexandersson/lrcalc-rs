# Broader corner-case review, 2026-10-08

Repair follow-up: release builds now use `panic = "unwind"`, matching Rust
integration tests. Rebuilding a release library and rerunning release tests
in the previously affected cache no longer emits the artifact-collision
warnings or the BigInt identity error. CI now exercises both release tests
and the repeated library-build/test sequence. The historical observations
below describe the pre-repair state; no counting algorithm changed here.
Full release tests and strict clippy pass. Both C smoke programs also pass
against shared and static release libraries (four independently linked
binaries); no ABI boundary was changed.

Baseline: 4022bd3, isolated from the canonical checkout's uncommitted
optimization edits. This is a test-only checkpoint. No new lrcalc-rs failure
was found in this bounded review; it is not an exhaustive certification.

`tests/broad_corner_cases.rs` adds five deterministic tests:

- 4,147 LR triples, outer size at most eight and at most four rows:
  interchange the two lower partitions, simultaneously conjugate all three,
  append zero parts, compare full Buch/GT counts and dimensions.
- 1,477 beta-LR cases, outer size at most five, three content coordinates
  including zeros, and seven beta choices including non-dominant beta:
  compare scalar counts with content expansions, test the stretch polynomial
  at dilations 2, 3, and 5 against fresh direct counts, and check dimension
  invariance under dilation by two.
- Packed state sizes just below, at, and above the 128-bit limit, using
  large near-rectangular shapes whose skew difference is just one box.
- Partition totals beyond `i32::MAX`, again with a one-box skew difference:
  the public APIs must return one or a checked error, not panic/wrong count.
- Insert an unused content coordinate at every position in a one-dimensional
  Kostka fiber. At dilation n=1..4, full and interior counts remain n+1 and
  n-1, respectively.

All five tests pass in debug and release. The conjugation oracle is a tiny
test-only column count; no public conjugation API is available here. The
other tests use the existing library. No new production computation or
polynomial generator was introduced. Cross-backend agreement alone does
not exclude a shared mathematical error; the symmetry and closed-form
checks provide additional checks that do not depend on dimension code.

```sh
# Set CARGO_TARGET_DIR to an external build cache.
cargo test --locked --test broad_corner_cases -- --nocapture
cargo test --release --locked --test broad_corner_cases -- --nocapture
```

The full debug suite and strict clippy also pass. Release Cargo emits a
library-output filename-collision warning for the rlib/staticlib/cdylib
targets. The first release run passed; a later rebuild in the reused cache
failed with incompatible `BigInt` crate identities despite the same version.
A fresh external target directory (`/cargo-target/lrcalc-broad-audit-20261008`)
then built successfully and all five tests passed again. The build-artifact
collision remains a separate reproducibility concern; this review does not
fix it or claim its cause is fully diagnosed. The release profile uses
`panic = "abort"`, while test builds also require unwind artifacts.

This review does not exhaustively audit the C ABI, allocation failure, very
large state spaces, or all optional features.
