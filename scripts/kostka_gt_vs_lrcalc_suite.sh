#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
RUST_BIN="${RUST_BIN:-$ROOT_DIR/target/release/lrcalc}"
GT_BIN="${GT_BIN:-/workspace/rust/target/release/kostka}"
UPSTREAM_BIN="${UPSTREAM_BIN:-/workspace/references/lrcalc-upstream/src/lrcalc}"
REPEAT="${REPEAT:-120}"

case "${1:-}" in
  -h|--help)
    cat <<'EOF'
Usage: scripts/kostka_gt_vs_lrcalc_suite.sh [repeat]

Compare the existing GT-DP Kostka implementation in /workspace/rust/kostka
against upstream C lrcalc, using the Kostka-to-LR translation as the oracle
bridge. The timed section repeats a mixed fixed suite.

Environment:
  RUST_BIN       lrcalc-rs binary used only for Kostka-to-LR triples
  GT_BIN         /workspace/rust/kostka binary to test
  UPSTREAM_BIN   upstream C lrcalc binary used as oracle
  REPEAT         timed repetitions of the case list, default 120
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

echo "building GT-DP kostka release binary..."
(cd /workspace/rust && timeout 60s nice -n 10 cargo build --release -p kostka >/dev/null)

if [[ ! -x "$RUST_BIN" ]]; then
  echo "lrcalc-rs binary not found or not executable: $RUST_BIN" >&2
  exit 2
fi
if [[ ! -x "$GT_BIN" ]]; then
  echo "GT-DP kostka binary not found or not executable: $GT_BIN" >&2
  exit 2
fi

CASE_LABELS=()
CASE_SHAPES=()
CASE_WEIGHTS=()
CASE_TRIPLES=()

add_case() {
  CASE_LABELS+=("$1")
  CASE_SHAPES+=("$2")
  CASE_WEIGHTS+=("$3")
}

# Small and degenerate cases.
add_case "row shape small" "3" "2,1"
add_case "standard 21" "2,1" "1,1,1"
add_case "dominance zero" "1,1,1" "2,1"
add_case "single row all ones" "6" "1,1,1,1,1,1"
add_case "single column heavy zero" "1,1,1,1" "3,1"

# Classical matrix-size examples.
add_case "balanced 321" "3,2,1" "2,2,2"
add_case "example 421" "4,2,1" "3,2,1,1"
add_case "standard 321" "3,2,1" "1,1,1,1,1,1"
add_case "standard 4321" "4,3,2,1" "1,1,1,1,1,1,1,1,1,1"

# Stretched-looking balanced weights.
add_case "stretched 321 x2" "6,4,2" "4,4,4"
add_case "stretched 321 x3" "9,6,3" "6,6,6"
add_case "stretched 421 x2" "8,4,2" "6,4,2,2"
add_case "stretched 421 x3" "12,6,3" "9,6,3"

# Irregular shapes and weights.
add_case "irregular A" "5,4,2,1" "4,3,2,2,1"
add_case "irregular B" "6,4,3,1" "5,3,3,2,1"
add_case "irregular C" "7,5,3,2" "6,4,3,2,2"
add_case "irregular D" "8,6,4,2" "6,5,4,3,2"
add_case "irregular E" "9,7,4,2,1" "7,5,4,3,2,2"

# Many-label cases stress the GT chain depth and the LR row concatenation.
add_case "many labels 42" "4,2" "1,1,1,1,1,1"
add_case "many labels 5321" "5,3,2,1" "1,1,1,1,1,1,1,1,1,1,1"
add_case "many labels 7431" "7,4,3,1" "1,1,1,1,1,1,1,1,1,1,1,1,1,1,1"

csv_to_cli_parts() {
  local csv="$1"
  local spaced="${csv//,/ }"
  if [[ -z "$spaced" ]]; then
    printf '0'
  else
    printf '%s' "$spaced"
  fi
}

gt_value() {
  local shape="$1"
  local weight="$2"
  local out
  out="$("$GT_BIN" kostka --lambda "$shape" --weight "$weight" --format json)"
  sed -n 's/.*"value":"\([^"]*\)".*/\1/p' <<< "$out"
}

lr_triple() {
  local shape="$1"
  local weight="$2"
  local shape_parts weight_parts triple
  shape_parts="$(csv_to_cli_parts "$shape")"
  weight_parts="$(csv_to_cli_parts "$weight")"
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
  gt_out="$(gt_value "${CASE_SHAPES[$i]}" "${CASE_WEIGHTS[$i]}")"
  upstream_out="$(lrcalc_value_from_triple "$triple")"
  if [[ -z "$gt_out" || "$gt_out" != "$upstream_out" ]]; then
    echo "mismatch: ${CASE_LABELS[$i]}" >&2
    echo "shape:  ${CASE_SHAPES[$i]}" >&2
    echo "weight: ${CASE_WEIGHTS[$i]}" >&2
    echo "LR triple: $triple" >&2
    echo "GT-DP:    $gt_out" >&2
    echo "lrcalc:   $upstream_out" >&2
    exit 1
  fi
done
echo "correctness: ok"

time_gt() {
  local start end elapsed
  start="$(now_ns)"
  for (( pass = 0; pass < REPEAT; pass++ )); do
    for i in "${!CASE_LABELS[@]}"; do
      gt_value "${CASE_SHAPES[$i]}" "${CASE_WEIGHTS[$i]}" >/dev/null
    done
  done
  end="$(now_ns)"
  elapsed=$(( end - start ))
  LAST_ELAPSED_NS="$elapsed"
  printf '%-10s %s  (%d evals)\n' "GT-DP:" "$(format_ns "$elapsed")" \
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
echo "GT-DP binary:  $GT_BIN"
echo "lrcalc binary: $UPSTREAM_BIN"
suite_start="$(now_ns)"
time_gt
gt_ns="$LAST_ELAPSED_NS"
time_lrcalc
lrcalc_ns="$LAST_ELAPSED_NS"
suite_end="$(now_ns)"
echo "comparison: GT-DP/lrcalc = $(ratio "$gt_ns" "$lrcalc_ns")x"
echo "total timed: $(format_ns "$(( suite_end - suite_start ))")"
