# Upstream ABI Notes

Reference source:
<https://bitbucket.org/asbuch/lrcalc/src/master/>

Snapshot inspected locally:

- commit: `8705a16e1575351684ed692f88552f97da3724f4`
- date: `2025-04-16 17:28:16 -0400`
- subject: `lrcalc tab: added weight argument.`

## Public Header Surface

The C headers in `src/` expose these main groups.

- `ivector.h`: flexible-array integer vector with `uint32_t length` and
  `int32_t array[1]`.
- `ivlincomb.h`: hash-table style integer linear combination keyed by
  `ivector *`, with iterator support.
- `schur.h`: Schur multiplication, fusion multiplication, skew Schur expansion,
  coproduct, and single LR coefficient.
- `lrcoef.h`: low-level LR coefficient count.
- `lriter.h`: LR tableau iterator.
- `part.h`: partition validation, conjugation, printing, and quantum/fusion
  partition helpers.
- `schublib.h`: Schubert polynomial transitions and multiplication.

The Python binding `python/liblrcalc.pxd` currently declares the ABI most
important for Python drop-in compatibility:

- `iv_new`, `iv_free`
- `part_qdegree`, `part_qentry`
- `ivlc_free_all`, `ivlc_first`, `ivlc_good`, `ivlc_next`,
  `ivlc_key`, `ivlc_value`
- `schur_mult`, `schur_mult_fusion`, `schur_skew`, `schur_coprod`,
  `schur_lrcoef`
- `trans`, `monk`, `mult_schubert`, `mult_schubert_str`
- `lrit_new`, `lrit_good`, `lrit_next`, `lrit_free`

## ABI Constraints

- `ivector` must remain `#[repr(C)]` compatible with:
  `uint32_t length; int32_t array[1];`.
- `iv_free` cannot rely on the current `length` field for deallocation because
  upstream code mutates `length` in routines such as partition chopping.
- ABI allocation should therefore use C `malloc`/`calloc`/`free` or another
  allocation scheme that does not need the original layout at free time.
- Native Rust data structures should not leak into the ABI boundary.

## Current ABI State

Implemented in Rust:

1. `ivector`, `ilist`, `ivlist`, and `ivlincomb` layouts and exported helper
   functions, with C allocation/free semantics at the ABI boundary.
2. `schur_lrcoef`, `lrcoef_count`, Schur product/fusion/skew/coproduct
   functions returning `ivlincomb`, and LR tableau iterator functions.
3. Partition helpers, quantum printing helpers, permutation/string helpers,
   Schubert functions, Maple printing helpers, and `optim_*`/`sksh_*` helpers.
4. Compatibility headers under `include/lrcalc/`, including the public
   allocator and template headers.  The `*_new_init` constructors are true
   exported C-variadic symbols, and the staged concrete headers also expose
   source-level variadic wrappers for rebuilt C callers.
5. `lrcalc` and `schubmult` CLI commands.  The staged prefix installs both
   binaries next to the library.

Verified smoke coverage:

1. `scripts/c_abi_smoke.sh` stages headers, shared/static libraries, and CLI
   binaries, then compiles and runs C programs against the staged prefix,
   covering both header-inline and exported-symbol variadic constructors.
2. `scripts/python_bindings_smoke.sh` rebuilds the upstream Python Cython
   module against the staged Rust prefix.
3. `scripts/sage_bindings_smoke.sh` validates Sage's wrapper through
   `LD_PRELOAD`.

Remaining drop-in compatibility work:

1. Full Sage rebuild/install testing rather than `LD_PRELOAD` smoke testing.
2. A fresh exported-symbol inventory check against a release build.
3. A final decision on exact CLI output order; current correctness checks
   compare sorted terms where upstream hash iteration order differs.
4. Wider C smoke coverage for printing helpers, iterator edge cases, and
   unusual ownership-transfer calls.
