#ifndef LRCALC_SCHUBLIB_H
#define LRCALC_SCHUBLIB_H

#include "ivector.h"
#include "ivlincomb.h"

#ifdef __cplusplus
extern "C" {
#endif

ivlincomb *trans(ivector *w, int vars);
ivlincomb *monk(int i, ivlincomb *slc, int rank);
ivlincomb *mult_poly_schubert(ivlincomb *poly, ivector *perm, int rank);
ivlincomb *mult_schubert(ivector *w1, ivector *w2, int rank);
ivlincomb *mult_schubert_str(ivector *str1, ivector *str2);

#ifdef __cplusplus
}
#endif

#endif
