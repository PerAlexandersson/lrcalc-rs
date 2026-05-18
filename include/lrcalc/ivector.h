#ifndef LRCALC_IVECTOR_H
#define LRCALC_IVECTOR_H

#include <stdarg.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct ivector {
  uint32_t length;
  int32_t array[1];
} ivector;

#define iv_length(v) ((v)->length)
#define iv_elem(v, i) ((v)->array[(i)])

ivector *iv_new(uint32_t length);
ivector *iv_new_zero(uint32_t length);
ivector *iv_new_copy(ivector *v);
void iv_free(ivector *v);

static inline ivector *iv_new_init(uint32_t length, ...)
{
  ivector *v = iv_new(length);
  if (v == 0)
    return 0;

  va_list ap;
  va_start(ap, length);
  for (uint32_t i = 0; i < length; i++)
    v->array[i] = va_arg(ap, int);
  va_end(ap);

  return v;
}

void iv_set_zero(ivector *v);
int iv_cmp(ivector *v1, ivector *v2);
int32_t iv_hash(ivector *v);
int32_t iv_sum(ivector *v);
void iv_copy(ivector *dst, ivector *src);
int iv_lesseq(ivector *v1, ivector *v2);
void iv_mult(ivector *dst, int32_t c, ivector *src);
void iv_div(ivector *dst, ivector *src, int32_t c);
int32_t iv_max(ivector *v);
int32_t iv_min(ivector *v);
void iv_reverse(ivector *dst, ivector *src);
int32_t iv_gcd(ivector *v);
void iv_print(ivector *v);
void iv_printnl(ivector *v);

#ifdef __cplusplus
}
#endif

#endif
