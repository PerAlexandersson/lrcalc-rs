#include <stdlib.h>
#include <string.h>

#include "alloc.h"


typedef struct {
  SIZE_T length;
  VALUE_T array[1];
} VECTOR;

#ifdef DEBUG
INLINE VALUE_T * PREFIX(pelem)(VECTOR *v, SIZE_T i)
{
  claim(i >= 0);
  claim(i < v->length);
  return v->array + i;
}
#endif

INLINE VECTOR * PREFIX(new) (SIZE_T length)
{
  VECTOR *v = (VECTOR *)
    ml_malloc(sizeof(VECTOR) + length * sizeof(VALUE_T) - sizeof(VALUE_T));
  if (v == NULL)
    return NULL;
  v->length = length;
  return v;
}

INLINE VECTOR * PREFIX(new_zero) (SIZE_T length)
{
  VECTOR *v = (VECTOR *)
    ml_calloc(1, sizeof(VECTOR) + length * sizeof(VALUE_T) - sizeof(VALUE_T));
  if (v == NULL)
    return NULL;
  v->length = length;
  return v;
}

VECTOR * PREFIX(new_init) (SIZE_T length, ...);

INLINE void PREFIX(free) (VECTOR *v)
{
  ml_free(v);
}

INLINE VECTOR * PREFIX(new_copy) (VECTOR *v)
{
  VECTOR *vc = PREFIX(new)(v->length);
  if (vc == NULL)
    return NULL;
  memcpy(vc->array, v->array, v->length * sizeof(VALUE_T));
  return vc;
}

INLINE void PREFIX(set_zero) (VECTOR *v)
{
  memset(v->array, 0, v->length * sizeof(VALUE_T));
}

INLINE void PREFIX(copy) (VECTOR *d, VECTOR *s)
{
  claim(d->length == s->length);
  memcpy(d->array, s->array, d->length * sizeof(VALUE_T));
}

INLINE int PREFIX(cmp) (VECTOR *v1, VECTOR *v2)
{
  SIZE_T i;
  if (v1->length != v2->length)
    return v1->length - v2->length;
  for (i = 0; i < v1->length; i++)
    if (v1->array[i] != v2->array[i])
      return v1->array[i] - v2->array[i];
  return 0;
}

#ifdef INTEGER_VALUE
INLINE int32_t PREFIX(hash) (VECTOR *v)
{
  SIZE_T i;
  uint32_t h;
  h = v->length;
  for (i = 0; i < v->length; i++)
    h = ((h << 5) ^ (h >> 27)) + v->array[i];
  return h;
}
#endif

VALUE_T PREFIX(sum) (VECTOR *v);
int PREFIX(lesseq) (VECTOR *v1, VECTOR *v2);
void PREFIX(mult) (VECTOR *dst, VALUE_T c, VECTOR *src);
void PREFIX(div) (VECTOR *dst, VECTOR *src, VALUE_T c);
VALUE_T PREFIX(max) (VECTOR *v);
VALUE_T PREFIX(min) (VECTOR *v);
void PREFIX(reverse) (VECTOR *dst, VECTOR *src);

#ifdef INTEGER_VALUE
VALUE_T PREFIX(gcd) (VECTOR *v);
#endif

void PREFIX(print) (VECTOR *v);
void PREFIX(printnl) (VECTOR *v);
