# ABI Symbol Inventory Worker Note

Role: `abi-worker`.

Scope checked:

- upstream clone: `/tmp/lrcalc-upstream`
- upstream commit: `8705a16e1575351684ed692f88552f97da3724f4`
- local Rust workspace guide: `/workspace/rust/AGENTS.md`

## Upstream Header Inventory

`src/Makefile.am` installs all `src/*.h` under `include/lrcalc`.
Template headers are therefore public too, but their concrete ABI names come
from `ilist.h`, `ivector.h`, `ivlist.h`, and `ivlincomb.h`.

### `alloc.h`

Release builds map these to libc/no-op macros:

- `ml_malloc`, `ml_calloc`, `ml_realloc`, `ml_free`
- `alloc_reset`, `alloc_getenv`, `alloc_report`
- `alloc_set_print`, `alloc_set_trap`, `alloc_set_fail`
- `alloc_set_trap_number`, `alloc_test_oom`

With `DEBUG_MEMORY`, the same names are real exported functions.

### `ivector.h`

Types/macros:

- `ivector { uint32_t length; int32_t array[1]; }`
- `iv_length(v)`
- `iv_elem(v,i)`

Concrete functions from `vector.tpl.h`:

- `iv_new`
- `iv_new_zero`
- `iv_new_init`
- `iv_free`
- `iv_new_copy`
- `iv_set_zero`
- `iv_copy`
- `iv_cmp`
- `iv_hash`
- `iv_sum`
- `iv_lesseq`
- `iv_mult`
- `iv_div`
- `iv_max`
- `iv_min`
- `iv_reverse`
- `iv_gcd`
- `iv_print`
- `iv_printnl`

Current Rust export check after the ABI hardening pass: release `nm -D` has no
missing upstream exported symbols; the only intentional extra symbol is
`lrcalc_new_abi_version`.  The `*_new_init` symbols are implemented as true
C-variadic exports through `src/native/abi_variadic.c`, and the concrete staged
headers also provide source-level variadic wrappers for rebuilt C callers.

### `ilist.h`

Types/macros:

- `ilist { int *array; size_t allocated; size_t length; }`
- `il_length(lst)`
- `il_elem(lst,i)`

Concrete functions from `list.tpl.h`:

- `il_init`
- `il_new`
- `il_new_init`
- `il_dealloc`
- `il_free`
- `il_reset`
- `il__realloc_array`
- `il_makeroom`
- `il_append`
- `il_poplast`
- `il_insert`
- `il_delete`
- `il_fastdelete`
- `il_extend`
- `il_copy`
- `il_new_copy`
- `il_reverse`

### `ivlist.h`

Types/macros:

- `ivlist { ivector **array; size_t allocated; size_t length; }`
- `ivl_length(lst)`
- `ivl_elem(lst,i)`

Concrete functions:

- `ivl_init`
- `ivl_new`
- `ivl_new_init`
- `ivl_dealloc`
- `ivl_free`
- `ivl_reset`
- `ivl__realloc_array`
- `ivl_makeroom`
- `ivl_append`
- `ivl_poplast`
- `ivl_insert`
- `ivl_delete`
- `ivl_fastdelete`
- `ivl_extend`
- `ivl_copy`
- `ivl_new_copy`
- `ivl_reverse`
- `ivl_free_all`

### `ivlincomb.h`

Types/constants:

- `ivlc_keyval_t { ivector *key; int32_t value; uint32_t hash; uint32_t next; }`
- `ivlincomb { uint32_t *table; ivlc_keyval_t *elts; uint32_t card;`
  `uint32_t free_elts; uint32_t elts_len; uint32_t elts_sz;`
  `uint32_t table_sz; }`
- `ivlc_iter { ivlincomb *ht; size_t index; size_t i; }`
- `IVLC_HASHTABLE_SZ`
- `IVLC_ARRAY_SZ`
- `LC_COPY_KEY`, `LC_FREE_KEY`, `LC_FREE_ZERO`, `LC_KEEP_ZERO`

Concrete functions from `hashtab.tpl.h`:

- `ivlc_card`
- `ivlc_init`
- `ivlc_new`
- `ivlc_dealloc`
- `ivlc_free`
- `ivlc_reset`
- `ivlc__grow_table`
- `ivlc__grow_elts`
- `ivlc_makeroom`
- `ivlc_lookup`
- `ivlc_insert`
- `ivlc_remove`
- `ivlc_equals`
- `ivlc_print_stat`
- `ivlc_good`
- `ivlc_first`
- `ivlc_next`
- `ivlc_key`
- `ivlc_value`
- `ivlc_keyval`
- `ivlc_dealloc_refs`
- `ivlc_dealloc_all`
- `ivlc_free_all`
- `ivlc_add_element`
- `ivlc_add_multiple`
- `ivlc_print`

ABI note: `python/liblrcalc.pxd` writes `ivlc_iter.index` and `ivlc_iter.i`
as `uint32_t`, but the installed C header uses `size_t`. Match the C header
layout for Rust ABI and generated compatibility headers.

### `part.h`

Types/constants:

- `part_iter { ivector *part; ivector *outer; ivector *inner; int length;`
  `int rows; int opt; }`
- `PITR_USE_OUTER`, `PITR_USE_INNER`, `PITR_USE_SIZE`

Functions:

- `part_valid`
- `part_decr`
- `part_length`
- `part_entry`
- `part_chop`
- `part_unchop`
- `part_leq`
- `part_conj`
- `part_print`
- `part_printnl`
- `part_print_lincomb`
- `part_qdegree`
- `part_qentry`
- `part_qprint`
- `part_qprintnl`
- `part_qprint_lincomb`
- `pitr_good`
- `pitr_first`
- `pitr_box_first`
- `pitr_box_sz_first`
- `pitr_sub_first`
- `pitr_sub_sz_first`
- `pitr_between_first`
- `pitr_between_sz_first`
- `pitr_next`

### `lrcoef.h`

- `lrcoef_count`

### `lriter.h`

Types:

- `lrit_box { int value; int max; int above; int right; }`
- `lrtab_iter { ivector *cont; int size; int array_len; lrit_box array[1]; }`

Functions:

- `lrit_new`
- `lrit_free`
- `lrit_print_skewtab`
- `lrit_dump`
- `lrit_dump_skew`
- `lrit_good`
- `lrit_next`
- `lrit_count`
- `lrit_expand`

### `schur.h`

- `schur_mult`
- `fusion_reduce`
- `fusion_reduce_lc`
- `schur_mult_fusion`
- `schur_skew`
- `schur_coprod`
- `schur_lrcoef`

### `schublib.h`

- `trans`
- `monk`
- `mult_poly_schubert`
- `mult_schubert`
- `mult_schubert_str`

### `perm.h`

- `perm_valid`
- `perm_length`
- `perm_group`
- `dimvec_valid`
- `bruhat_leq`
- `bruhat_zero`
- `str_iscompat`
- `all_strings`
- `all_perms`
- `string2perm`
- `str2dimvec`
- `perm2string`

### Other Headers

- `maple.h`: `maple_print_lincomb`, `maple_qprint_lincomb`
- `optshape.h`: `skew_shape`, `sksh_print`, `optim_mult`, `optim_fusion`,
  `optim_skew`, `optim_coef`, `sksh_dealloc`
- `vectarg.h`: `get_vect_arg`
- `list.tpl.h`, `vector.tpl.h`, `hashtab.tpl.h`: public templates, not
  standalone concrete ABI names without the instantiating macros.

## Python `.pxd` Inventory

`python/liblrcalc.pxd` exposes this compatibility-critical subset:

### `lrcalc/ivector.h`

- `ivector { uint32_t length; int32_t array[1]; }`
- `iv_new(uint32_t length) -> ivector *`
- `iv_free(ivector *v) -> void`

### `lrcalc/part.h`

- `part_qdegree(ivector *p, int level) -> int`
- `part_qentry(ivector *p, int i, int d, int level) -> int`

### `lrcalc/ivlincomb.h`

- opaque `ivlincomb`
- `ivlc_iter { ivlincomb *ht; uint32_t index; uint32_t i; }`
- `ivlc_free_all(ivlincomb *lc) -> void`
- `ivlc_first(ivlincomb *lc, ivlc_iter *itr) -> void`
- `ivlc_good(ivlc_iter *itr) -> bint`
- `ivlc_next(ivlc_iter *itr) -> void`
- `ivlc_key(ivlc_iter *itr) -> ivector *`
- `ivlc_value(ivlc_iter *itr) -> int32_t`

### `lrcalc/schur.h`

- `schur_mult(ivector *sh1, ivector *sh2, int rows, int cols, int partsz)`
  `-> ivlincomb *`
- `schur_mult_fusion(ivector *sh1, ivector *sh2, int rows, int level)`
  `-> ivlincomb *`
- `schur_skew(ivector *outer, ivector *inner, int rows, int partsz)`
  `-> ivlincomb *`
- `schur_coprod(ivector *sh, int rows, int cols, int partsz, int all)`
  `-> ivlincomb *`
- `schur_lrcoef(ivector *outer, ivector *inner1, ivector *inner2)`
  `-> long long`

### `lrcalc/schublib.h`

- `trans(ivector *w, int vars) -> ivlincomb *`
- `monk(int i, ivlincomb *slc, int rank) -> ivlincomb *`
- `mult_schubert(ivector *w1, ivector *w2, int rank) -> ivlincomb *`
- `mult_schubert_str(ivector *str1, ivector *str2) -> ivlincomb *`

### `lrcalc/lriter.h`

- `lrit_box { int value; int max; int above; int right; }`
- `lrtab_iter { ivector *cont; int size; int array_len; lrit_box array[1]; }`
- `lrit_new(ivector *outer, ivector *inner, ivector *content,`
  `int maxrows, int maxcols, int partsz) -> lrtab_iter *`
- `lrit_good(lrtab_iter *lrit) -> int`
- `lrit_next(lrtab_iter *lrit) -> void`
- `lrit_free(lrtab_iter *lrit) -> void`

## ABI Priority Tiers

P0, implemented coefficient target:

- exact `ivector` layout and allocation/free behavior
- `iv_new`, `iv_new_zero`, `iv_new_copy`, `iv_free`
- `iv_set_zero`, `iv_cmp`, `iv_hash`, `iv_sum`
- `schur_lrcoef`
- `lrcalc coef` output compatibility

P1, Python/Sage Schur compatibility:

- all `.pxd` Schur functions: `schur_mult`, `schur_mult_fusion`,
  `schur_skew`, `schur_coprod`, `schur_lrcoef`
- `ivlincomb` layout enough for iteration and ownership
- `ivlc_free_all`, `ivlc_first`, `ivlc_good`, `ivlc_next`,
  `ivlc_key`, `ivlc_value`
- `part_qdegree`, `part_qentry`
- `lrit_new`, `lrit_good`, `lrit_next`, `lrit_free`

P2, Python/Sage Schubert compatibility:

- `trans`
- `monk`
- `mult_schubert`
- `mult_schubert_str`
- likely also `mult_poly_schubert` for C-header users

P3, full installed-header C compatibility:

- remaining `iv_*`, `il_*`, `ivl_*`, `ivlc_*`
- `part_*`, `pitr_*`, `perm_*`
- `lrcoef_count`, `lrit_count`, `lrit_expand`
- `maple_*`, `optim_*`, `get_vect_arg`

P4, diagnostics and debug-only compatibility:

- debug-memory `alloc_*` and `ml_*` functions
- print/dump helpers not needed by Python import or normal Sage use

## Current LR Coefficient Strategy

The first coefficient target is implemented in `src/lrcoef.rs`.

- `schur_lrcoef` converts ABI `ivector` inputs to slices and calls
  `lrcoef_i64`.
- `lrcoef` uses a Rust version of upstream-style `optim_coef`
  compactification and then a Buch/Chaput-style branch-pruned tableau search.
- Counts use checked `u128` internally and downcast to the upstream
  `long long` ABI return type at `schur_lrcoef`.
- Additional native paths exist for comparison and research:
  `src/lr_gt.rs`, `src/lr_signed.rs`, `src/kostka_fast.rs`,
  and `src/lr_ehrhart.rs`.

Current ABI follow-up: keep the exported-symbol inventory in CI, broaden the C
smoke corpus for iterator and printing edge cases, and decide whether exact CLI
output order must follow upstream hash iteration.

## Benchmark And Test Plan

Build references:

- upstream C: build `/tmp/lrcalc-upstream` in `/tmp`, keep its `lrcalc` binary
  and `liblrcalc.so`
- Rust: use `timeout 60s nice -n 10 cargo test` and
  `timeout 60s nice -n 10 cargo build --release`

ABI checks:

- `nm -D --defined-only` on both shared libraries; compare against P0/P1/P2
  required names first, then full header inventory.
- C layout probe compiled against compatibility headers:
  `sizeof(ivector)`, `offsetof(ivector,array)`, `sizeof(ivlc_iter)`,
  `offsetof(ivlc_iter,index)`, `sizeof(lrit_box)`, `sizeof(lrtab_iter)`.
- C smoke program for `iv_new`, element writes, `iv_hash`, `iv_sum`,
  `iv_new_copy`, `iv_cmp`, `iv_free`, and `schur_lrcoef`.

Functional comparisons:

- Exhaust all partitions up to total size 8 first:
  compare `lrcalc coef outer - inner1 - inner2` stdout and exit status.
- Add random valid triples up to size 20 and selected large examples from the
  upstream README.
- When P1 is implemented, compare `mult`, `skew`, `coprod`, `mult_fusion`,
  `mult_quantum`, and `lr_iterator` through the upstream Python `lrcalc.pyx`
  rebuilt against the Rust library.
- Keep output-order checks separate from dictionary-value checks; upstream
  hashtable iteration order is not a mathematical contract.

Benchmarks:

- Record CSV rows: operation, outer, inner1, inner2, coefficient, C time,
  Rust time, peak states if available.
- Start with `coef` triples by size and by coefficient magnitude.
- Add CLI-level timing and direct C-ABI timing; CLI startup can hide small
  coefficient costs.
- Use the same generated input corpus for correctness and timing.

## Local Rust Reuse Candidates

Checked `/workspace/rust` under the local guide.

Good candidates:

- `/workspace/rust/kostka/src/lr.rs`
  - Contains `lr_dp(lambda, mu, nu, max_states) -> BigUint`.
  - Uses a horizontal-strip DP with Yamanouchi constraints.
  - Good as a reference algorithm and test oracle.
- `/workspace/rust/kostka/src/kostka_dp.rs`
  - Contains horizontal-strip enumeration and skew Kostka DP.
  - Useful for later `schur_skew` and multiplication.
- `/workspace/rust/kostka/src/partition.rs`
  - Small MIT partition type with containment, conjugation, hook lengths, and
    parsing.
- `/workspace/rust/combinatoric-core/src/partition.rs`
  - Richer MIT partition type with skew boxes, strip removal, dominance, and
    `kostka_to_lr`.
- `/workspace/rust/sym-poly/core/src/tableau.rs`
  - Has `Tableau`, `SkewTableau`, standard skew tableau iterators, and word
    utilities.
- `/workspace/rust/sym-poly/sym/src/symmetric_function.rs`
  - Has skew Schur via Jacobi-Trudi and basis conversion machinery.

Reasons not to depend on them immediately:

- `lrcalc-rs` is intended as standalone ABI crate; path dependencies on
  `/workspace/rust` would make packaging and replacement builds fragile.
- `kostka` pulls broad dependencies and CLI/database-adjacent features not
  needed at the ABI boundary.
- Local types use `u32`/`BigUint` and owned Rust structures; ABI code needs
  narrow `int32_t`/`long long` behavior and C allocation boundaries.
- `sym-poly-sym` is useful for algebraic verification, but Jacobi-Trudi plus
  basis conversion is not the simplest first implementation of one LR
  coefficient.

Recommendation: copy or reimplement the small partition and LR-DP ideas inside
`lrcalc-rs` with attribution if needed, then consider a feature-gated local
dependency only after ABI compatibility is stable.
