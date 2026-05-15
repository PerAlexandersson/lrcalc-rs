#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
REPEAT="${1:-1}"

"$ROOT_DIR/scripts/lr_signed_ffi_bench.sh" "$REPEAT" large-three-part
