#ifndef LRCALC_PERM_H
#define LRCALC_PERM_H

#include "ivector.h"
#include "ivlist.h"

#ifdef __cplusplus
extern "C" {
#endif

int perm_valid(ivector *w);
int perm_length(ivector *w);
int perm_group(ivector *w);
int dimvec_valid(ivector *dv);
int bruhat_leq(ivector *w1, ivector *w2);
int bruhat_zero(ivector *w1, ivector *w2, int rank);

int str_iscompat(ivector *str1, ivector *str2);
ivlist *all_strings(ivector *dimvec);
ivlist *all_perms(int n);

ivector *string2perm(ivector *str);
ivector *str2dimvec(ivector *str);
ivector *perm2string(ivector *perm, ivector *dimvec);

#ifdef __cplusplus
}
#endif

#endif
