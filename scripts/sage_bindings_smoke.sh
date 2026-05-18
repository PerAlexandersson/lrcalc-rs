#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PREFIX="${PREFIX:-"$ROOT/target/lrcalc-rs-prefix"}"
SAGE_PYTHON="${SAGE_PYTHON:-/workspace/.conda-envs/sage/bin/python}"
PRELOAD_LIB="$PREFIX/lib/liblrcalc.so.2"

if [[ ! -x "$SAGE_PYTHON" ]]; then
  cat >&2 <<EOF
Sage Python not found:
  $SAGE_PYTHON

Set SAGE_PYTHON to a Sage-enabled Python executable.
EOF
  exit 2
fi

"$ROOT/scripts/stage_liblrcalc_prefix.sh" "$PREFIX" >/dev/null

LD_PRELOAD="$PRELOAD_LIB${LD_PRELOAD:+:$LD_PRELOAD}" "$SAGE_PYTHON" - <<PY
from sage.libs.lrcalc import lrcalc

expected_lib = "$PREFIX/lib/liblrcalc.so.2.0.0"

def part_dict(d):
    return {tuple(k): int(v) for k, v in d.items()}

def pair_dict(d):
    return {(tuple(a), tuple(b)): int(v) for (a, b), v in d.items()}

with open("/proc/self/maps", "r", encoding="utf-8") as maps:
    loaded = sorted({line.rsplit(None, 1)[-1] for line in maps if "liblrcalc" in line})

if loaded != [expected_lib]:
    raise AssertionError(f"expected {expected_lib!r} to be loaded, got {loaded!r}")

assert lrcalc.lrcoef([3, 2, 1], [2, 1], [2, 1]) == 2
assert part_dict(lrcalc.mult([1], [1])) == {(2,): 1, (1, 1): 1}
assert part_dict(lrcalc.skew([2, 1], [1])) == {(2,): 1, (1, 1): 1}
assert pair_dict(lrcalc.coprod([2, 1])) == {
    ((2, 1), ()): 1,
    ((2,), (1,)): 1,
    ((1, 1), (1,)): 1,
}
assert part_dict(lrcalc.mult_schubert([2, 1], [2, 1])) == {(3, 1, 2): 1}
words = list(lrcalc.lrskew([2, 1], [1]))
assert len(words) == 2, words

print("sage lrcalc smoke passed")
print(f"loaded liblrcalc: {loaded[0]}")
PY
