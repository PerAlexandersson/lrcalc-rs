#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PREFIX="${1:-"$ROOT/target/lrcalc-rs-prefix"}"

cd "$ROOT"
timeout 60s nice -n 10 cargo build --release

rm -rf "$PREFIX/include/lrcalc" "$PREFIX/lib"
mkdir -p "$PREFIX/include" "$PREFIX/lib"

cp -R "$ROOT/include/lrcalc" "$PREFIX/include/"
cp "$ROOT/target/release/liblrcalc.so" "$PREFIX/lib/liblrcalc.so.2.0.0"
cp "$ROOT/target/release/liblrcalc.a" "$PREFIX/lib/liblrcalc.a"

ln -sfn liblrcalc.so.2.0.0 "$PREFIX/lib/liblrcalc.so.2"
ln -sfn liblrcalc.so.2 "$PREFIX/lib/liblrcalc.so"

printf '%s\n' "$PREFIX"
