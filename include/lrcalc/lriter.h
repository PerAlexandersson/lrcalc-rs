#ifndef LRCALC_LRITER_H
#define LRCALC_LRITER_H

#include "ivector.h"
#include "ivlincomb.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct lrit_box {
  int value;
  int max;
  int above;
  int right;
} lrit_box;

typedef struct lrtab_iter {
  ivector *cont;
  int size;
  int array_len;
  lrit_box array[1];
} lrtab_iter;

lrtab_iter *lrit_new(
    ivector *outer, ivector *inner, ivector *content,
    int maxrows, int maxcols, int partsz);
void lrit_free(lrtab_iter *lrit);
int lrit_good(lrtab_iter *lrit);
void lrit_next(lrtab_iter *lrit);
ivlincomb *lrit_count(lrtab_iter *lrit);
ivlincomb *lrit_expand(
    ivector *outer, ivector *inner, ivector *content,
    int maxrows, int maxcols, int partsz);
void lrit_print_skewtab(lrtab_iter *lrit, ivector *outer, ivector *inner);
void lrit_dump(lrtab_iter *lrit);
void lrit_dump_skew(lrtab_iter *lrit, ivector *outer, ivector *inner);

#ifdef __cplusplus
}
#endif

#endif
