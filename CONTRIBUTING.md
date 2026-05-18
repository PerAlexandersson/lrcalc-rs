# Contributing

`lrcalc-rs` aims to be a drop-in Rust implementation of Anders Buch's
`lrcalc` library and command-line tools.  Changes should preserve the C ABI and
the CLI behavior unless the difference is documented as intentional.

## Development Setup

Required tools:

- Rust stable with Cargo
- a C compiler available as `cc`
- `ar`

Useful optional tools:

- an upstream `lrcalc` checkout for ABI and timing comparisons
- `nm` for exported-symbol checks
- `perf` or Samply for profiling

## Checks

Run these before opening a substantial pull request:

```bash
timeout 60s nice -n 10 cargo fmt --check
timeout 60s nice -n 10 cargo clippy --all-targets -- -D warnings
timeout 60s nice -n 10 cargo test
timeout 60s nice -n 10 cargo build --release
timeout 60s nice -n 10 scripts/c_abi_smoke.sh
```

If an upstream checkout is available at `/tmp/lrcalc-upstream`, also run:

```bash
timeout 120s nice -n 10 scripts/python_bindings_smoke.sh
timeout 120s nice -n 10 scripts/sage_bindings_smoke.sh
timeout 120s nice -n 10 scripts/sage_lrcalc_bench.sh
```

When comparing against upstream C, point benchmark scripts at an upstream build
with `UPSTREAM_LIB`, `UPSTREAM_LIB_DIR`, or `UPSTREAM_BIN` as appropriate.
The scripts under `scripts/` print the variables they expect; see
[scripts/README.md](scripts/README.md) for the script index.

## Compatibility Rules

- Keep ABI-facing allocation and iterator types isolated in `src/abi.rs`.
- Keep native Rust data structures separate from C-facing ownership rules.
- Do not remove an exported symbol without first documenting the compatibility
  impact.
- Keep output-order differences explicit; value equality is not always enough
  for drop-in command-line compatibility.
- Prefer small, focused benchmark cases when changing a hot counting loop.

## Notes

Current public algorithm, ABI, and benchmark notes live under `notes/`.
One-off profiling logs and local benchmark output should stay outside the
tracked tree, usually under `target/`.
