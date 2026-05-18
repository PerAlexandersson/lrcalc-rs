#ifndef LRCALC_OPTSHAPE_H
#define LRCALC_OPTSHAPE_H

#include "ivector.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct skew_shape {
  ivector *outer;
  ivector *inner;
  ivector *cont;
  int sign;
} skew_shape;

void sksh_dealloc(skew_shape *ss);
void sksh_print(ivector *outer, ivector *inner, ivector *cont);

int optim_mult(skew_shape *ss, ivector *sh1, ivector *sh2, int maxrows, int maxcols);
int optim_fusion(skew_shape *ss, ivector *sh1, ivector *sh2, int maxrows, int level);
int optim_skew(skew_shape *ss, ivector *outer, ivector *inner, ivector *content, int maxrows);
int optim_coef(skew_shape *ss, ivector *out, ivector *sh1, ivector *sh2);

#ifdef __cplusplus
}
#endif

#endif
