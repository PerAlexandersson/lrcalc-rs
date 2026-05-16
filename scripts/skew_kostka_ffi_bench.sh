#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
UPSTREAM_LIB="${UPSTREAM_LIB:-/tmp/lrcalc-upstream/src/.libs/liblrcalc.so}"
REPEAT="${1:-2000}"

if [[ ! -f "$UPSTREAM_LIB" ]]; then
  cat >&2 <<EOF
upstream static library not found:
  $UPSTREAM_LIB

Build it first:
  cd /tmp/lrcalc-upstream
  autoreconf -i
  ./configure --enable-shared --disable-static
  make -j2
EOF
  exit 2
fi

cd "$ROOT_DIR"
UPSTREAM_LIB="$UPSTREAM_LIB" \
  timeout 60s nice -n 10 cargo run --release --bin skew_kostka_ffi_bench -- "$REPEAT"
