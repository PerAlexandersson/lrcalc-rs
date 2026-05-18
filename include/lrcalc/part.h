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
int part_qdegree(ivector *p, int level);
int part_qentry(ivector *p, int i, int d, int level);
void part_qprint(ivector *p, int level);
void part_qprintnl(ivector *p, int level);

#ifdef __cplusplus
}
#endif

#endif
