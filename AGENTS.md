# lrcalc-rs Agent Guide

This is the `lrcalc-rs` standalone Rust project, currently checked out under
`/workspace/documents/lrcalc-new`.

## Scope

- Build a Rust implementation of Anders Buch's `lrcalc`.
- Preserve the public C ABI closely enough to act as a drop-in `liblrcalc`
  replacement.
- Preserve the command-line surface of the `lrcalc` executable.
- Prefer performance-oriented exact combinatorial algorithms.

## References

- Upstream source: <https://bitbucket.org/asbuch/lrcalc/src/master/>
- Local throwaway clone may be made under `/tmp/lrcalc-upstream`.
- Local reusable Rust code lives under `/workspace/rust`; check
  `/workspace/rust/AGENTS.md` before copying or depending on it.

## Working Rules

- Use Rust and Markdown as the main work surfaces.
- Do not vendor large chunks of upstream C code unless the licensing and purpose
  are explicit in notes.
- Keep ABI-facing types and allocation functions isolated from native Rust data
  structures.
- Use `timeout 60s nice -n 10 cargo ...` for Rust commands.
- Coordinate through `HANDOFF.md` and short notes under `notes/`.
- Avoid editing the same file as another worker.
