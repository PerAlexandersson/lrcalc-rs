#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
REPEAT="${1:-20000}"

cd "$ROOT_DIR"
timeout 60s nice -n 10 cargo run --release --bin kostka_bench -- "$REPEAT"
