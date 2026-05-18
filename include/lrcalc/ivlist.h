#ifndef LRCALC_IVLIST_H
#define LRCALC_IVLIST_H

#include <stdarg.h>
#include <stddef.h>

#include "ivector.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct ivlist {
  ivector **array;
  size_t allocated;
  size_t length;
} ivlist;

#define ivl_length(lst) ((lst)->length)
#define ivl_elem(lst, i) ((lst)->array[(i)])

int ivl_init(ivlist *lst, size_t sz);
ivlist *ivl_new(size_t sz);
void ivl_dealloc(ivlist *lst);
void ivl_free(ivlist *lst);
void ivl_free_all(ivlist *lst);
void ivl_reset(ivlist *lst);
int ivl__realloc_array(ivlist *lst, size_t sz);
int ivl_makeroom(ivlist *lst, size_t sz);
int ivl_append(ivlist *lst, ivector *x);
static inline ivlist *ivl_new_init(size_t sz, size_t count, ...)
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
ivector *ivl_poplast(ivlist *lst);
int ivl_insert(ivlist *lst, size_t i, ivector *x);
ivector *ivl_delete(ivlist *lst, size_t i);
ivector *ivl_fastdelete(ivlist *lst, size_t i);
int ivl_extend(ivlist *dst, ivlist *src);
int ivl_copy(ivlist *dst, ivlist *src);
ivlist *ivl_new_copy(ivlist *lst);
int ivl_reverse(ivlist *dst, ivlist *src);

#ifdef __cplusplus
}
#endif

#endif
