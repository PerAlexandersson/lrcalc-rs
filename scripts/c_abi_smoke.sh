#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PREFIX="${PREFIX:-"$ROOT/target/lrcalc-rs-c-abi-prefix"}"
CC_BIN="${CC:-cc}"

"$ROOT/scripts/stage_liblrcalc_prefix.sh" "$PREFIX" >/dev/null

mkdir -p "$ROOT/target/c-abi-smoke"
"$CC_BIN" \
  -std=c11 \
  -Wall \
  -Wextra \
  -Werror \
  -I"$PREFIX/include" \
  "$ROOT/tests/c_abi_smoke.c" \
  -L"$PREFIX/lib" \
  -Wl,-rpath,"$PREFIX/lib" \
  -llrcalc \
  -o "$ROOT/target/c-abi-smoke/c_abi_smoke"

LD_LIBRARY_PATH="$PREFIX/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}" \
  "$ROOT/target/c-abi-smoke/c_abi_smoke"

"$CC_BIN" \
  -std=c11 \
  -Wall \
  -Wextra \
  -Werror \
  "$ROOT/tests/c_abi_symbol_smoke.c" \
  -L"$PREFIX/lib" \
  -Wl,-rpath,"$PREFIX/lib" \
  -llrcalc \
  -o "$ROOT/target/c-abi-smoke/c_abi_symbol_smoke"

LD_LIBRARY_PATH="$PREFIX/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}" \
  "$ROOT/target/c-abi-smoke/c_abi_symbol_smoke"

"$CC_BIN" \
  -std=c11 \
  -Wall \
  -Wextra \
  -Werror \
  -I"$PREFIX/include" \
  "$ROOT/tests/c_abi_smoke.c" \
  "$PREFIX/lib/liblrcalc.a" \
  -ldl \
  -lpthread \
  -lm \
  -o "$ROOT/target/c-abi-smoke/c_abi_smoke_static"

"$ROOT/target/c-abi-smoke/c_abi_smoke_static"

"$CC_BIN" \
  -std=c11 \
  -Wall \
  -Wextra \
  -Werror \
  "$ROOT/tests/c_abi_symbol_smoke.c" \
  "$PREFIX/lib/liblrcalc.a" \
  -ldl \
  -lpthread \
  -lm \
  -o "$ROOT/target/c-abi-smoke/c_abi_symbol_smoke_static"

"$ROOT/target/c-abi-smoke/c_abi_symbol_smoke_static"

"$PREFIX/bin/lrcalc" mult 1 - 1 >/dev/null
"$PREFIX/bin/schubmult" 1 2 - 1 2 >/dev/null
