#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PREFIX="${1:-"$ROOT/target/lrcalc-rs-prefix"}"

case "$PREFIX" in
  ""|"/"|"/usr"|"/usr/"|"/usr/local"|"/usr/local/")
    if [[ "${LRCALC_RS_STAGE_FORCE:-}" != "1" ]]; then
      echo "refusing to stage into unsafe prefix: $PREFIX" >&2
      echo "set LRCALC_RS_STAGE_FORCE=1 to override" >&2
      exit 2
    fi
    ;;
esac

cd "$ROOT"
timeout 60s nice -n 10 cargo build --release

rm -rf "$PREFIX/include/lrcalc"
mkdir -p "$PREFIX/include" "$PREFIX/lib" "$PREFIX/bin"
rm -f \
  "$PREFIX/lib/liblrcalc.a" \
  "$PREFIX/lib/liblrcalc.so" \
  "$PREFIX/lib/liblrcalc.so.2" \
  "$PREFIX/lib/liblrcalc.so.2.0.0" \
  "$PREFIX/bin/lrcalc" \
  "$PREFIX/bin/schubmult"

cp -R "$ROOT/include/lrcalc" "$PREFIX/include/"
cp "$ROOT/target/release/liblrcalc.so" "$PREFIX/lib/liblrcalc.so.2.0.0"
cp "$ROOT/target/release/liblrcalc.a" "$PREFIX/lib/liblrcalc.a"
cp "$ROOT/target/release/lrcalc" "$PREFIX/bin/lrcalc"
cp "$ROOT/target/release/schubmult" "$PREFIX/bin/schubmult"

ln -sfn liblrcalc.so.2.0.0 "$PREFIX/lib/liblrcalc.so.2"
ln -sfn liblrcalc.so.2 "$PREFIX/lib/liblrcalc.so"

printf '%s\n' "$PREFIX"
