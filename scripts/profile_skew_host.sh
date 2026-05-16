#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DEFAULT_OUTPUT="/tmp/lrcalc-skew-large-$(date +%Y%m%d-%H%M%S).json.gz"
OUTPUT="$DEFAULT_OUTPUT"
RATE="${SAMPLY_RATE:-1000}"
OPEN=0
BUILD=1
TARGET_ARGS=()

usage() {
  cat <<'EOF'
Usage: scripts/profile_skew_host.sh [options] [-- lrcalc-args...]

Build lrcalc-new with debug info and frame pointers, then profile it with
samply.  This script is intended for a non-Docker host terminal unless the
Docker container has perf capabilities.

Default profiled command:
  target/release/lrcalc skew 30 27 24 21 18 15 / 15 12 9 6 3

Options:
  -o, --output PATH   Write samply profile to PATH
  -r, --rate HZ       Sampling rate, default 1000
      --open          Open the samply UI instead of save-only mode
      --no-build      Reuse the existing release binary
  -h, --help          Show this help

Environment:
  RUST_BIN            Binary to profile, default target/release/lrcalc
  RUSTFLAGS           Defaults to "-C debuginfo=1 -C force-frame-pointers=yes"
  SAMPLY_RATE         Default sampling rate when --rate is not passed
  ALLOW_DOCKER=1      Permit running inside Docker after capabilities are fixed
EOF
}

while (($#)); do
  case "$1" in
    -h|--help)
      usage
      exit 0
      ;;
    -o|--output)
      if [[ $# -lt 2 ]]; then
        echo "--output requires a path" >&2
        exit 2
      fi
      OUTPUT="$2"
      shift 2
      ;;
    -r|--rate)
      if [[ $# -lt 2 ]]; then
        echo "--rate requires a positive integer" >&2
        exit 2
      fi
      RATE="$2"
      shift 2
      ;;
    --open)
      OPEN=1
      shift
      ;;
    --no-build)
      BUILD=0
      shift
      ;;
    --)
      shift
      TARGET_ARGS=("$@")
      break
      ;;
    *)
      echo "unknown option: $1" >&2
      usage >&2
      exit 2
      ;;
  esac
done

if ! [[ "$RATE" =~ ^[0-9]+$ ]] || ((RATE == 0)); then
  echo "sampling rate must be a positive integer" >&2
  exit 2
fi

if ((${#TARGET_ARGS[@]} == 0)); then
  TARGET_ARGS=(skew 30 27 24 21 18 15 / 15 12 9 6 3)
fi

if [[ -f /.dockerenv && "${ALLOW_DOCKER:-0}" != 1 ]]; then
  cat >&2 <<'EOF'
This looks like Docker.  The current Docker profile blocks perf sampling unless
the container is recreated with the profiling capabilities recorded in:
  /workspace/DOCKER_PROFILING_NEXT_BAKE.md

Run this from a non-Docker host terminal, or set ALLOW_DOCKER=1 after fixing
the container capabilities.
EOF
  exit 2
fi

if ! command -v samply >/dev/null 2>&1; then
  cat >&2 <<'EOF'
samply is not on PATH.

Install it on the host with:
  cargo install samply --locked
EOF
  exit 2
fi

if ! command -v cargo >/dev/null 2>&1; then
  echo "cargo is not on PATH" >&2
  exit 2
fi

if [[ -r /proc/sys/kernel/perf_event_paranoid ]]; then
  paranoid="$(cat /proc/sys/kernel/perf_event_paranoid)"
  if [[ "$paranoid" =~ ^[0-9]+$ ]] && ((paranoid > 1)); then
    cat >&2 <<EOF
/proc/sys/kernel/perf_event_paranoid is $paranoid.

For non-root sampling, lower it on the host:
  echo 1 | sudo tee /proc/sys/kernel/perf_event_paranoid
EOF
    exit 2
  fi
fi

cd "$ROOT_DIR"

if ((BUILD)); then
  export RUSTFLAGS="${RUSTFLAGS:--C debuginfo=1 -C force-frame-pointers=yes}"
  echo "building release binary with RUSTFLAGS=$RUSTFLAGS"
  if command -v timeout >/dev/null 2>&1 && command -v nice >/dev/null 2>&1; then
    timeout "${CARGO_BUILD_TIMEOUT:-120s}" nice -n "${NICE_LEVEL:-10}" cargo build --release
  else
    cargo build --release
  fi
fi

RUST_BIN="${RUST_BIN:-$ROOT_DIR/target/release/lrcalc}"
if [[ ! -x "$RUST_BIN" ]]; then
  echo "Rust binary not found or not executable: $RUST_BIN" >&2
  exit 2
fi

samply_args=(
  record
  --rate "$RATE"
  --output "$OUTPUT"
  --profile-name "lrcalc-new ${TARGET_ARGS[*]}"
)
if ((OPEN == 0)); then
  samply_args+=(--save-only)
fi

echo "profiling: $RUST_BIN ${TARGET_ARGS[*]}"
echo "output:    $OUTPUT"
exec samply "${samply_args[@]}" "$RUST_BIN" "${TARGET_ARGS[@]}"
