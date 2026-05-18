#include <stdarg.h>
#include <stddef.h>
#include <stdint.h>

typedef struct ivector {
  uint32_t length;
  int32_t array[1];
} ivector;

typedef struct ilist {
  int32_t *array;
  size_t allocated;
  size_t length;
} ilist;

typedef struct ivlist {
  ivector **array;
  size_t allocated;
  size_t length;
} ivlist;

extern ivector *iv_new(uint32_t length);
extern ilist *il_new(size_t sz);
extern void il_free(ilist *lst);
extern int il_append(ilist *lst, int32_t x);
extern ivlist *ivl_new(size_t sz);
extern void ivl_free(ivlist *lst);
extern int ivl_append(ivlist *lst, ivector *x);

ivector *lrcalc_c_iv_new_init(uint32_t length, ...)
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

ilist *lrcalc_c_il_new_init(size_t sz, size_t count, ...)
{
  ilist *lst = il_new(sz);
  if (lst == 0)
    return 0;

  va_list ap;
  va_start(ap, count);
  for (size_t i = 0; i < count; i++) {
    int value = va_arg(ap, int);
    if (il_append(lst, value) != 0) {
      va_end(ap);
      il_free(lst);
      return 0;
    }
  }
  va_end(ap);

  return lst;
}

ivlist *lrcalc_c_ivl_new_init(size_t sz, size_t count, ...)
{
  ivlist *lst = ivl_new(sz);
  if (lst == 0)
    return 0;

  va_list ap;
  va_start(ap, count);
  for (size_t i = 0; i < count; i++) {
    ivector *value = va_arg(ap, ivector *);
    if (ivl_append(lst, value) != 0) {
      va_end(ap);
      ivl_free(lst);
      return 0;
    }
  }
  va_end(ap);

  return lst;
}
