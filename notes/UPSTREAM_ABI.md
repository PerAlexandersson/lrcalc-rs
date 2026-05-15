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

1. `ivector` layout plus `iv_new`, `iv_new_zero`, `iv_new_copy`, `iv_free`,
   `iv_set_zero`, `iv_cmp`, `iv_hash`, and `iv_sum`.
2. `schur_lrcoef` backed by the native Buch-style coefficient engine.
3. `lrcalc coef` / `lrcalc lrcoef` CLI output for single coefficients.

Still missing for Python/Sage drop-in compatibility:

1. `ivlincomb` allocation, insertion, iteration, and `ivlc_free_all`.
2. Schur expansion functions returning `ivlincomb`: `schur_mult`,
   `schur_mult_fusion`, `schur_skew`, and `schur_coprod`.
3. `part_qdegree` and `part_qentry` ABI exports for quantum output.
4. LR tableau iterator ABI: `lrit_new`, `lrit_good`, `lrit_next`,
   and `lrit_free`.
5. Schubert ABI: `trans`, `monk`, `mult_schubert`,
   and `mult_schubert_str`.
6. Compatibility headers and C/Python smoke tests.
