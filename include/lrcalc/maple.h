#ifndef LRCALC_MAPLE_H
#define LRCALC_MAPLE_H

#include "ivlincomb.h"

#ifdef __cplusplus
extern "C" {
#endif

void maple_print_lincomb(ivlincomb *ht, char *letter, int nz);
void maple_qprint_lincomb(ivlincomb *lc, int level, char *letter);

#ifdef __cplusplus
}
#endif

#endif
