# lrcalc-new Handoff

## Current Goal

Build a standalone Rust replacement for Anders Buch's `lrcalc`, with compatible
C ABI and CLI behavior where possible, while improving native coefficient
engines.

## Current State

- Rust crate builds as `liblrcalc` (`cdylib`, `staticlib`, `rlib`) plus a
  `lrcalc` binary.
- `src/abi.rs` exports the `ivector` allocation/copy/hash/sum subset and
  `schur_lrcoef`.
- `src/lrcoef.rs` contains the primary Buch-style LR coefficient engine:
  compactification, pruned tableau search, interior counts, dimension, and
  stretch-cache helpers.
- `src/lr_gt.rs` contains an independent GT-chain LR DP with stats, dimension,
  and relative-interior variants.  It also has a hybrid full-count selector
  that tries exact Kostka translation, certified partial Kostka collapse, then
  Buch fallback.  Its production paired full/interior selector tries exact
  Kostka translation, then Buch full/interior counts.
- `src/kostka_fast.rs` contains packed `u128` ordinary/skew Kostka DP and
  interior counts.
- `src/lr_signed.rs` contains the signed Kostka expansion for LR coefficients.
- `src/lr_ehrhart.rs` interpolates h-vectors for pure stretched LR families.
- CLI commands expose the coefficient engines and diagnostic/stat modes, but
  not upstream `mult`, `skew`, `coprod`, `tab`, or `schubmult` behavior.
- Benchmark scripts compare selected Rust paths against upstream C when an
  upstream binary is available.
- `src/bin/stretched_dp_bench.rs` compares paired full/interior counts for
  stretched Kostka DP and the equivalent stretched LR GT-chain DP.
- `src/bin/lr_hybrid_bench.rs` compares raw GT full and paired counts with the
  combined tableau hybrid selectors on exact and near-Kostka families.
- `src/bin/partial_collapse_probe.rs` explores row-masked Yamanouchi DPs for
  partial Kostka collapse candidates.

## Verified

`timeout 60s nice -n 10 cargo test` passed on 2026-05-16:
55 library tests, all benchmark-bin test targets, and doc-tests.

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

## Main Gaps

- Implement `ivlincomb` ABI storage and iterators.
- Implement C ABI functions for Schur product, skew, coproduct, fusion/quantum,
  LR tableau iteration, and Schubert products.
- Add compatibility headers and C smoke tests for struct layout and exported
  symbols.
- Add Python/Sage rebuild tests against the Rust `liblrcalc`.
- Decide which native LR engine should serve each workload class after broader
  benchmarks.

## Next Useful Work

- Build the minimal `ivlincomb` ABI needed by Python: allocation, insertion,
  iteration, and `ivlc_free_all`.
- Implement `schur_skew` or `schur_mult` next, using the existing LR/Kostka
  engines and returning `ivlincomb`.
- Add an `nm`-based exported-symbol check for P0/P1 ABI coverage.
- Run `scripts/lrcoef_timed_suite.sh` and the Kostka/LR benchmark scripts
  against a freshly built upstream C binary.
- Profile the Buch-port, GT-chain, and signed-Kostka paths on the same corpus
  before adding new optimizations.
- Benchmark the production tableau paired selector on a broader non-Kostka
  corpus, since it now avoids the GT-chain fallback.
- Use `notes/PARTIAL_KOSTKA_COLLAPSE.md` to guide larger certified partial-mask
  experiments for full counts.
- Extend certified partial collapse to paired full/interior counts only after a
  safe interior certificate is available; for now paired counts use Buch
  fallback outside exact Kostka translations.

## Notes

- ABI inventory: `notes/ABI_SYMBOL_INVENTORY_WORKER.md`.
- Upstream ABI summary: `notes/UPSTREAM_ABI.md`.
- Algorithm survey and Rust roadmap: `notes/UPSTREAM_SOURCE_ALGORITHM_SURVEY.md`.
- Partial Kostka-collapse probe: `notes/PARTIAL_KOSTKA_COLLAPSE.md`.
