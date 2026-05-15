#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
RUST_BIN="${RUST_BIN:-$ROOT_DIR/target/release/lrcalc}"
UPSTREAM_BIN="${UPSTREAM_BIN:-/workspace/references/lrcalc-upstream/src/lrcalc}"

if [[ $# -eq 0 || "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
  cat <<'EOF'
Usage: scripts/kostka_compare_lrcalc.sh SHAPE - WEIGHT

Example:
  scripts/kostka_compare_lrcalc.sh 2 1 - 1 1 1

The script translates K_{SHAPE,WEIGHT} to a single LR coefficient and compares
the Rust result with upstream C lrcalc.
EOF
  exit 0
fi

if [[ ! -x "$UPSTREAM_BIN" ]]; then
  echo "upstream binary not found or not executable: $UPSTREAM_BIN" >&2
  exit 2
fi

(cd "$ROOT_DIR" && timeout 60s nice -n 10 cargo build --release >/dev/null)

if [[ ! -x "$RUST_BIN" ]]; then
  echo "Rust binary not found or not executable: $RUST_BIN" >&2
  exit 2
fi

triple="$("$RUST_BIN" kostka-lr-triple "$@")"
rust_out="$("$RUST_BIN" kostka-lr "$@")"
read -r -a triple_args <<< "$triple"
upstream_out="$("$UPSTREAM_BIN" coef "${triple_args[@]}")"

echo "Kostka input: $*"
echo "LR triple:    $triple"
echo "Rust:         $rust_out"
echo "upstream:     $upstream_out"

if [[ "$rust_out" != "$upstream_out" ]]; then
  echo "mismatch" >&2
  exit 1
fi
echo "comparison: ok"
