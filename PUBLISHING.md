# Public Repository Checklist

This checklist separates public Git readiness from a drop-in replacement
release.  The current repository is useful to publish as `lrcalc-rs`; the
compatibility release still needs a few packaging checks.

## Ready Before First Public Push

- Choose the public remote URL, then add it as `origin`.
- Add the same URL to `Cargo.toml` as `repository = "..."`
  once it is known.
- Confirm that `README.md`, `LICENSE`, `COPYING`, and `CONTRIBUTING.md` are
  acceptable as the first public-facing documentation.
- Decide whether to keep `notes/` in the public repo.  The notes are useful for
  development history, but they contain local benchmark paths and should not be
  treated as user documentation.
- Run the local pre-push checks:

```bash
timeout 60s nice -n 10 cargo fmt --check
timeout 60s nice -n 10 cargo test
timeout 60s nice -n 10 cargo build --release
timeout 120s nice -n 10 scripts/python_bindings_smoke.sh
```

## Drop-in Compatibility Release Blockers

- Add C smoke tests that link against the built `liblrcalc`.
- Audit and complete installed headers beyond the Python/Sage Cython surface,
  including the current non-variadic initializer shim caveat.
- Add Sage rebuild tests or a documented manual recipe.
- Re-run the upstream exported-symbol inventory against a release build.
- Re-run the benchmark report against a fresh upstream C checkout.
- Decide whether output order must exactly match upstream for each CLI command.

## Suggested First GitHub Settings

- Enable issues.
- Protect the default branch after the initial CI run passes.
- Add topics such as `lrcalc`, `littlewood-richardson`, `schur-functions`,
  `kostka`, `schubert-polynomials`, and `rust`.
- Add an initial release only after the C headers and smoke tests exist.
