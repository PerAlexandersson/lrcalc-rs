#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
RUST_BIN="${RUST_BIN:-$ROOT_DIR/target/release/lrcalc}"
UPSTREAM_BIN="${UPSTREAM_BIN:-/tmp/lrcalc-upstream/src/lrcalc}"
REPEAT="${REPEAT:-120}"

case "${1:-}" in
  -h|--help)
    cat <<'EOF'
Usage: scripts/kostka_fast_vs_lrcalc_suite.sh [repeat]

Compare lrcalc-rs's fast u128 Kostka DP against upstream C lrcalc,
using precomputed Kostka-to-LR triples for the lrcalc side.
EOF
    exit 0
    ;;
  "")
    ;;
  *)
    REPEAT="$1"
    ;;
esac

if ! [[ "$REPEAT" =~ ^[0-9]+$ ]] || (( REPEAT == 0 )); then
  echo "repeat must be a positive integer" >&2
  exit 2
fi
if [[ ! -x "$UPSTREAM_BIN" ]]; then
  echo "upstream binary not found or not executable: $UPSTREAM_BIN" >&2
  exit 2
fi

echo "building lrcalc-rs release binary..."
(cd "$ROOT_DIR" && timeout 60s nice -n 10 cargo build --release >/dev/null)

CASE_LABELS=()
CASE_SHAPES=()
CASE_WEIGHTS=()
CASE_TRIPLES=()

add_case() {
  CASE_LABELS+=("$1")
  CASE_SHAPES+=("$2")
  CASE_WEIGHTS+=("$3")
}

add_case "row shape small" "3" "2,1"
add_case "standard 21" "2,1" "1,1,1"
add_case "dominance zero" "1,1,1" "2,1"
add_case "single row all ones" "6" "1,1,1,1,1,1"
add_case "single column heavy zero" "1,1,1,1" "3,1"
add_case "balanced 321" "3,2,1" "2,2,2"
add_case "example 421" "4,2,1" "3,2,1,1"
add_case "standard 321" "3,2,1" "1,1,1,1,1,1"
add_case "standard 4321" "4,3,2,1" "1,1,1,1,1,1,1,1,1,1"
add_case "stretched 321 x2" "6,4,2" "4,4,4"
add_case "stretched 321 x3" "9,6,3" "6,6,6"
add_case "stretched 421 x2" "8,4,2" "6,4,2,2"
add_case "stretched 421 x3" "12,6,3" "9,6,3"
add_case "irregular A" "5,4,2,1" "4,3,2,2,1"
add_case "irregular B" "6,4,3,1" "5,3,3,2,1"
add_case "irregular C" "7,5,3,2" "6,4,3,2,2"
add_case "irregular D" "8,6,4,2" "6,5,4,3,2"
add_case "irregular E" "9,7,4,2,1" "7,5,4,3,2,2"
add_case "many labels 42" "4,2" "1,1,1,1,1,1"
add_case "many labels 5321" "5,3,2,1" "1,1,1,1,1,1,1,1,1,1,1"
add_case "many labels 7431" "7,4,3,1" "1,1,1,1,1,1,1,1,1,1,1,1,1,1,1"

csv_to_cli_parts() {
  printf '%s' "${1//,/ }"
}

fast_value() {
  local shape_parts weight_parts
  shape_parts="$(csv_to_cli_parts "$1")"
  weight_parts="$(csv_to_cli_parts "$2")"
  "$RUST_BIN" kostka-fast $shape_parts - $weight_parts
}

lr_triple() {
  local shape_parts weight_parts
  shape_parts="$(csv_to_cli_parts "$1")"
  weight_parts="$(csv_to_cli_parts "$2")"
  "$RUST_BIN" kostka-lr-triple $shape_parts - $weight_parts
}

lrcalc_value_from_triple() {
  local triple="$1"
  read -r -a triple_args <<< "$triple"
  "$UPSTREAM_BIN" coef "${triple_args[@]}"
}

now_ns() {
  date +%s%N
}

format_ns() {
  local ns="$1"
  local ms=$(( ns / 1000000 ))
  printf '%d.%03ds' "$(( ms / 1000 ))" "$(( ms % 1000 ))"
}

ratio() {
  awk -v num="$1" -v den="$2" 'BEGIN { printf "%.3f", num / den }'
}

echo "correctness pass: ${#CASE_LABELS[@]} cases"
for i in "${!CASE_LABELS[@]}"; do
  triple="$(lr_triple "${CASE_SHAPES[$i]}" "${CASE_WEIGHTS[$i]}")"
  CASE_TRIPLES+=("$triple")
  fast_out="$(fast_value "${CASE_SHAPES[$i]}" "${CASE_WEIGHTS[$i]}")"
  upstream_out="$(lrcalc_value_from_triple "$triple")"
  if [[ "$fast_out" != "$upstream_out" ]]; then
    echo "mismatch: ${CASE_LABELS[$i]}" >&2
    echo "shape: ${CASE_SHAPES[$i]}" >&2
    echo "weight: ${CASE_WEIGHTS[$i]}" >&2
    echo "LR triple: $triple" >&2
    echo "fast: $fast_out" >&2
    echo "lrcalc: $upstream_out" >&2
    exit 1
  fi
done
echo "correctness: ok"

time_fast() {
  local start end elapsed
  start="$(now_ns)"
  for (( pass = 0; pass < REPEAT; pass++ )); do
    for i in "${!CASE_LABELS[@]}"; do
      fast_value "${CASE_SHAPES[$i]}" "${CASE_WEIGHTS[$i]}" >/dev/null
    done
  done
  end="$(now_ns)"
  elapsed=$(( end - start ))
  LAST_ELAPSED_NS="$elapsed"
  printf '%-10s %s  (%d evals)\n' "fast:" "$(format_ns "$elapsed")" \
    "$(( REPEAT * ${#CASE_LABELS[@]} ))"
}

time_lrcalc() {
  local start end elapsed
  start="$(now_ns)"
  for (( pass = 0; pass < REPEAT; pass++ )); do
    for i in "${!CASE_LABELS[@]}"; do
      lrcalc_value_from_triple "${CASE_TRIPLES[$i]}" >/dev/null
    done
  done
  end="$(now_ns)"
  elapsed=$(( end - start ))
  LAST_ELAPSED_NS="$elapsed"
  printf '%-10s %s  (%d evals)\n' "lrcalc:" "$(format_ns "$elapsed")" \
    "$(( REPEAT * ${#CASE_LABELS[@]} ))"
}

echo "timed pass: repeat=$REPEAT, cases=${#CASE_LABELS[@]}"
suite_start="$(now_ns)"
time_fast
fast_ns="$LAST_ELAPSED_NS"
time_lrcalc
lrcalc_ns="$LAST_ELAPSED_NS"
suite_end="$(now_ns)"
echo "comparison: fast/lrcalc = $(ratio "$fast_ns" "$lrcalc_ns")x"
echo "total timed: $(format_ns "$(( suite_end - suite_start ))")"
