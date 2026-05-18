#ifndef LRCALC_IVLINCOMB_H
#define LRCALC_IVLINCOMB_H

#include <stddef.h>
#include <stdint.h>

#include "ivector.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct ivlc_keyval_t {
  ivector *key;
  int32_t value;
  uint32_t hash;
  uint32_t next;
} ivlc_keyval_t;

typedef struct ivlincomb {
  uint32_t *table;
  ivlc_keyval_t *elts;
  uint32_t card;
  uint32_t free_elts;
  uint32_t elts_len;
  uint32_t elts_sz;
  uint32_t table_sz;
} ivlincomb;

typedef struct ivlc_iter {
  ivlincomb *ht;
  size_t index;
  size_t i;
} ivlc_iter;

#define IVLC_HASHTABLE_SZ 2003
#define IVLC_ARRAY_SZ 100

#define LC_COPY_KEY 1
#define LC_FREE_KEY 0
#define LC_FREE_ZERO 2
#define LC_KEEP_ZERO 0

int ivlc_init(ivlincomb *ht, uint32_t tabsz, uint32_t eltsz);
ivlincomb *ivlc_new(uint32_t tabsz, uint32_t eltsz);
uint32_t ivlc_card(ivlincomb *ht);
void ivlc_dealloc(ivlincomb *ht);
void ivlc_free(ivlincomb *ht);
void ivlc_reset(ivlincomb *ht);
int ivlc__grow_table(ivlincomb *ht, uint32_t sz);
int ivlc__grow_elts(ivlincomb *ht, uint32_t sz);
int ivlc_makeroom(ivlincomb *ht, uint32_t sz);
ivlc_keyval_t *ivlc_lookup(ivlincomb *ht, ivector *key, uint32_t hash);
ivlc_keyval_t *ivlc_insert(
    ivlincomb *ht, ivector *key, uint32_t hash, int32_t value);
ivlc_keyval_t *ivlc_remove(ivlincomb *ht, ivector *key, uint32_t hash);
int ivlc_good(ivlc_iter *itr);
void ivlc_first(ivlincomb *ht, ivlc_iter *itr);
void ivlc_next(ivlc_iter *itr);
ivector *ivlc_key(ivlc_iter *itr);
int32_t ivlc_value(ivlc_iter *itr);
ivlc_keyval_t *ivlc_keyval(ivlc_iter *itr);
void ivlc_dealloc_refs(ivlincomb *ht);
void ivlc_dealloc_all(ivlincomb *ht);
void ivlc_free_all(ivlincomb *ht);
int ivlc_add_element(
    ivlincomb *ht, int32_t c, ivector *key, uint32_t hash, int opt);
int ivlc_add_multiple(ivlincomb *dst, int32_t c, ivlincomb *src, int opt);
int ivlc_equals(ivlincomb *ht1, ivlincomb *ht2, int opt_zero);
void ivlc_print(ivlincomb *ht, int opt_zero);
void ivlc_print_stat(ivlincomb *ht);

#ifdef __cplusplus
}
#endif

#endif
