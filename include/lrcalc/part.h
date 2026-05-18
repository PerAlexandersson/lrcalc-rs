#ifndef LRCALC_PART_H
#define LRCALC_PART_H

#include "ivector.h"
#include "ivlincomb.h"

#ifdef __cplusplus
extern "C" {
#endif

int part_valid(ivector *p);
int part_decr(ivector *p);
int part_length(ivector *p);
int part_entry(ivector *p, int i);
void part_chop(ivector *p);
void part_unchop(ivector *p, int len);
int part_leq(ivector *p1, ivector *p2);
ivector *part_conj(ivector *p);
void part_print(ivector *p);
void part_printnl(ivector *p);
void part_print_lincomb(ivlincomb *lc);
int part_qdegree(ivector *p, int level);
int part_qentry(ivector *p, int i, int d, int level);
void part_qprint(ivector *p, int level);
void part_qprintnl(ivector *p, int level);
void part_qprint_lincomb(ivlincomb *lc, int level);

typedef struct part_iter {
  ivector *part;
  ivector *outer;
  ivector *inner;
  int length;
  int rows;
  int opt;
} part_iter;

#define PITR_USE_OUTER 1
#define PITR_USE_INNER 2
#define PITR_USE_SIZE 4

int pitr_good(part_iter *itr);
int pitr_first(
    part_iter *itr, ivector *p, int rows, int cols,
    ivector *outer, ivector *inner, int size, int opt);
void pitr_box_first(part_iter *itr, ivector *p, int rows, int cols);
void pitr_box_sz_first(part_iter *itr, ivector *p, int rows, int cols, int size);
void pitr_sub_first(part_iter *itr, ivector *p, ivector *outer);
void pitr_sub_sz_first(part_iter *itr, ivector *p, ivector *outer, int size);
void pitr_between_first(part_iter *itr, ivector *p, ivector *outer, ivector *inner);
void pitr_between_sz_first(
    part_iter *itr, ivector *p, ivector *outer, ivector *inner, int size);
void pitr_next(part_iter *itr);

#ifdef __cplusplus
}
#endif

#endif
