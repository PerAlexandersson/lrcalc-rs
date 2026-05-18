#ifndef LRCALC_ILIST_H
#define LRCALC_ILIST_H

#include <stdarg.h>
#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct ilist {
  int32_t *array;
  size_t allocated;
  size_t length;
} ilist;

#define il_length(lst) ((lst)->length)
#define il_elem(lst, i) ((lst)->array[(i)])

int il_init(ilist *lst, size_t sz);
ilist *il_new(size_t sz);
void il_dealloc(ilist *lst);
void il_free(ilist *lst);
void il_reset(ilist *lst);
int il__realloc_array(ilist *lst, size_t sz);
int il_makeroom(ilist *lst, size_t sz);
int il_append(ilist *lst, int32_t x);
static inline ilist *il_new_init(size_t sz, size_t count, ...)
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
int32_t il_poplast(ilist *lst);
int il_insert(ilist *lst, size_t i, int32_t x);
int32_t il_delete(ilist *lst, size_t i);
int32_t il_fastdelete(ilist *lst, size_t i);
int il_extend(ilist *dst, ilist *src);
int il_copy(ilist *dst, ilist *src);
ilist *il_new_copy(ilist *lst);
int il_reverse(ilist *dst, ilist *src);

#ifdef __cplusplus
}
#endif

#endif
