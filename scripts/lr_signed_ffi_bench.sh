#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
UPSTREAM_LIB_DIR="${UPSTREAM_LIB_DIR:-/tmp/lrcalc-upstream/src/.libs}"
UPSTREAM_LIB="$UPSTREAM_LIB_DIR/liblrcalc.a"
if [[ $# -eq 0 ]]; then
  ARGS=(5000 mixed)
else
  ARGS=("$@")
fi

if [[ ! -f "$UPSTREAM_LIB" ]]; then
  cat >&2 <<EOF
upstream static library not found:
  $UPSTREAM_LIB

Build it first:
  cd /tmp/lrcalc-upstream
  autoreconf -i
  ./configure --disable-shared --enable-static
  make -j2
EOF
  exit 2
fi

cd "$ROOT_DIR"
RUSTFLAGS="--cfg upstream_lrcalc_ffi -L native=$UPSTREAM_LIB_DIR" \
  timeout 60s nice -n 10 cargo run --release --bin lr_signed_ffi_bench -- "${ARGS[@]}"
