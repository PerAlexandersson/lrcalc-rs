#ifndef LRCALC_SCHUR_H
#define LRCALC_SCHUR_H

#include "ivector.h"
#include "ivlincomb.h"

#ifdef __cplusplus
extern "C" {
#endif

ivlincomb *schur_mult(
    ivector *sh1, ivector *sh2, int rows, int cols, int partsz);
int fusion_reduce(ivector *la, int level, ivector *tmp);
int fusion_reduce_lc(ivlincomb *lc, int level);
ivlincomb *schur_mult_fusion(ivector *sh1, ivector *sh2, int rows, int level);
ivlincomb *schur_skew(ivector *outer, ivector *inner, int rows, int partsz);
ivlincomb *schur_coprod(
    ivector *sh, int rows, int cols, int partsz, int all);
long long schur_lrcoef(ivector *outer, ivector *inner1, ivector *inner2);

#ifdef __cplusplus
}
#endif

#endif
