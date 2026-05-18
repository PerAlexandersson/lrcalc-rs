#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
UPSTREAM_PY="${UPSTREAM_PY:-/tmp/lrcalc-upstream/python}"
PREFIX="${PREFIX:-"$ROOT/target/lrcalc-rs-python-prefix"}"
VENV="${VENV:-/tmp/lrcalc-rs-python-smoke}"

if [[ ! -f "$UPSTREAM_PY/setup.py" || ! -f "$UPSTREAM_PY/lrcalc.pyx" ]]; then
  cat >&2 <<EOF
Upstream Python bindings not found at:
  $UPSTREAM_PY

Set UPSTREAM_PY to the upstream lrcalc/python directory.
EOF
  exit 2
fi

"$ROOT/scripts/stage_liblrcalc_prefix.sh" "$PREFIX" >/dev/null

python3 -m venv "$VENV"
# shellcheck disable=SC1091
source "$VENV/bin/activate"

python -m pip install --upgrade pip setuptools wheel cython >/dev/null

export CFLAGS="-I$PREFIX/include ${CFLAGS:-}"
export LDFLAGS="-L$PREFIX/lib -Wl,-rpath,$PREFIX/lib ${LDFLAGS:-}"

python -m pip install \
  --force-reinstall \
  --no-build-isolation \
  --no-cache-dir \
  "$UPSTREAM_PY" >/dev/null

LD_LIBRARY_PATH="$PREFIX/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}" python - <<'PY'
import lrcalc

def require(actual, expected, label):
    if actual != expected:
        raise AssertionError(f"{label}: expected {expected!r}, got {actual!r}")

require(lrcalc.lrcoef((3, 2, 1), (2, 1), (2, 1)), 2, "lrcoef")

require(
    lrcalc.mult((1,), (1,)),
    {(2,): 1, (1, 1): 1},
    "mult",
)

require(
    lrcalc.skew((2, 1), (1,)),
    {(2,): 1, (1, 1): 1},
    "skew",
)

require(
    lrcalc.coprod((2, 1)),
    {((2, 1), ()): 1, ((2,), (1,)): 1, ((1, 1), (1,)): 1},
    "coprod",
)

require(
    lrcalc.mult_fusion((1,), (1,), 3, 2),
    {(2,): 1, (1, 1): 1},
    "mult_fusion",
)

require(
    lrcalc.schubmult((2, 1), (2, 1)),
    {(3, 1, 2): 1},
    "schubmult",
)

require(
    lrcalc.schubmult_str((0, 1), (1, 0)),
    {(1, 0): 1},
    "schubmult_str",
)

words = list(lrcalc.lr_iterator((2, 1), (1,)))
if len(words) != 2 or any(len(word) != 2 for word in words):
    raise AssertionError(f"lr_iterator returned unexpected words: {words!r}")

print("python binding smoke passed")
print(f"module: {lrcalc.__file__}")
PY
