#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
RUST_BIN="${RUST_BIN:-$ROOT_DIR/target/release/lrcalc}"
UPSTREAM_BIN="${UPSTREAM_BIN:-/workspace/references/lrcalc-upstream/src/lrcalc}"
REPEAT="${REPEAT:-250}"

case "${1:-}" in
  -h|--help)
    cat <<'EOF'
Usage: scripts/lrcoef_timed_suite.sh [repeat]

Compare the Rust single LR coefficient CLI against upstream C lrcalc, then
time repeated runs of both binaries on a fixed mixed suite.

Environment:
  RUST_BIN       Rust lrcalc binary to test
  UPSTREAM_BIN   upstream C lrcalc binary used as oracle
  REPEAT         timed repetitions of the case list, default 250
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
  cat >&2 <<EOF
upstream binary not found or not executable:
  $UPSTREAM_BIN

Build it first, for example:
  cd /workspace/references/lrcalc-upstream
  autoreconf -i
  ./configure --disable-shared --enable-static
  make -j2
EOF
  exit 2
fi

echo "building Rust release binary..."
(cd "$ROOT_DIR" && timeout 60s nice -n 10 cargo build --release >/dev/null)

if [[ ! -x "$RUST_BIN" ]]; then
  echo "Rust binary not found or not executable: $RUST_BIN" >&2
  exit 2
fi

CASE_LABELS=()
CASE_ARGS=()

add_case() {
  CASE_LABELS+=("$1")
  shift
  CASE_ARGS+=("$*")
}

# Empty, zero, containment, and size-mismatch cases.
add_case "empty coefficient" 0 - 0 - 0
add_case "empty outer nonempty content" 0 - 0 - 1
add_case "size mismatch zero" 1 - 0 - 0
add_case "inner not contained zero" 2 2 1 - 3 - 2

# Small Schur products, including the standard s_21 * s_21 expansion.
add_case "single box Pieri" 2 1 - 2 - 1
add_case "s21 s21 coefficient 2" 3 2 1 - 2 1 - 2 1
add_case "s21 s21 coefficient 1 A" 4 2 - 2 1 - 2 1
add_case "s21 s21 coefficient 1 B" 3 3 - 2 1 - 2 1
add_case "s21 s21 absent shape" 5 1 - 2 1 - 2 1

# Pieri-like rows and columns.
add_case "horizontal strip short" 5 3 1 - 3 2 1 - 2 1
add_case "vertical strip short" 4 3 2 1 - 3 2 1 - 1 1 1
add_case "row times row" 9 - 5 - 4
add_case "column times column" 2 2 2 1 1 - 1 1 1 1 - 1 1 1

# Upstream OOM-suite coefficient examples.
add_case "upstream oom medium" 5 4 3 2 1 - 3 2 1 - 4 3 1 1
add_case "upstream oom larger" 7 6 5 4 3 2 1 - 4 4 3 2 1 - 5 4 3 2

# Compactification-heavy stretched-looking cases.
add_case "stretched 321 scale 2" 6 4 2 - 4 2 - 4 2
add_case "stretched 321 scale 3" 9 6 3 - 6 3 - 6 3
add_case "wide compact small" 20 14 8 2 - 12 8 4 - 14 6 2
add_case "wide compact medium" 30 20 10 - 20 10 - 20 10
add_case "wide compact skewed" 28 21 15 6 - 18 12 6 - 20 14

# More irregular partitions for branch-pruned tableau counting.
add_case "irregular A" 9 7 5 3 1 - 5 3 2 1 - 6 5 3 2
add_case "irregular B" 10 8 6 4 2 - 6 4 3 1 - 7 6 4 3
add_case "irregular C" 11 9 7 4 2 1 - 7 5 3 2 - 8 7 4 2 1
add_case "irregular D" 12 10 8 5 3 1 - 7 6 4 2 - 9 8 5 3 1

run_coef() {
  local bin="$1"
  shift
  "$bin" coef "$@"
}

split_args() {
  local -n out_ref="$1"
  local arg_string="$2"
  read -r -a out_ref <<< "$arg_string"
}

now_ns() {
  date +%s%N
}

format_ns() {
  local ns="$1"
  local ms=$(( ns / 1000000 ))
  printf '%d.%03ds' "$(( ms / 1000 ))" "$(( ms % 1000 ))"
}

echo "correctness pass: ${#CASE_ARGS[@]} cases"
for i in "${!CASE_ARGS[@]}"; do
  args=()
  split_args args "${CASE_ARGS[$i]}"
  rust_out="$(run_coef "$RUST_BIN" "${args[@]}")"
  upstream_out="$(run_coef "$UPSTREAM_BIN" "${args[@]}")"
  if [[ "$rust_out" != "$upstream_out" ]]; then
    echo "mismatch: ${CASE_LABELS[$i]}" >&2
    echo "args: ${CASE_ARGS[$i]}" >&2
    echo "rust:     $rust_out" >&2
    echo "upstream: $upstream_out" >&2
    exit 1
  fi
done
echo "correctness: ok"

time_binary() {
  local name="$1"
  local bin="$2"
  local start end elapsed
  start="$(now_ns)"
  for (( pass = 0; pass < REPEAT; pass++ )); do
    for arg_string in "${CASE_ARGS[@]}"; do
      args=()
      split_args args "$arg_string"
      run_coef "$bin" "${args[@]}" >/dev/null
    done
  done
  end="$(now_ns)"
  elapsed=$(( end - start ))
  LAST_ELAPSED_NS="$elapsed"
  printf '%-10s %s  (%d evals)\n' "$name:" "$(format_ns "$elapsed")" \
    "$(( REPEAT * ${#CASE_ARGS[@]} ))"
}

ratio() {
  awk -v num="$1" -v den="$2" 'BEGIN { printf "%.3f", num / den }'
}

echo "timed pass: repeat=$REPEAT, cases=${#CASE_ARGS[@]}"
echo "upstream binary: $UPSTREAM_BIN"
echo "rust binary:     $RUST_BIN"
suite_start="$(now_ns)"
time_binary "upstream" "$UPSTREAM_BIN"
upstream_ns="$LAST_ELAPSED_NS"
time_binary "rust" "$RUST_BIN"
rust_ns="$LAST_ELAPSED_NS"
suite_end="$(now_ns)"
echo "comparison: rust/upstream = $(ratio "$rust_ns" "$upstream_ns")x"
echo "total timed: $(format_ns "$(( suite_end - suite_start ))")"
