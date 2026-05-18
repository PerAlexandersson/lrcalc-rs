/*
 * List type for elements of user defined types.
 *
 * Macros used:
 *
 * LIST           :  Name of list type
 * PREFIX(name)   :  Returns e.g. l_name
 * VALUE_T        :  Type of list elements, e.g. void *
 * SIZE_T         :  Type for lengths, e.g. size_t
 */

#include <stdlib.h>
#include <string.h>

#include "alloc.h"


typedef struct {
  VALUE_T *array;
  SIZE_T allocated;
  SIZE_T length;
} LIST;

#ifdef DEBUG
INLINE VALUE_T * PREFIX(pelem)(LIST *lst, SIZE_T i)
{
  claim(i >= 0);
  claim(i < lst->length);
  return lst->array + i;
}
#endif

/* Initialize list structure. */
INLINE int PREFIX(init) (LIST *lst, SIZE_T sz)
{
  lst->array = (VALUE_T *) ml_malloc(sz * sizeof(VALUE_T));
  if (lst->array == NULL)
    return -1;
  lst->allocated = sz;
  lst->length = 0;
  return 0;
}

INLINE LIST * PREFIX(new) (SIZE_T sz)
{
  LIST *lst = (LIST *) ml_malloc(sizeof(LIST));
  if (lst == NULL)
    return NULL;
  if (PREFIX(init) (lst, sz) != 0)
    {
      ml_free(lst);
      return NULL;
    }
  return lst;
}

LIST * PREFIX(new_init)(SIZE_T sz, SIZE_T count, ...);

INLINE void PREFIX(dealloc) (LIST *v)
{
  ml_free(v->array);
}

INLINE void PREFIX(free) (LIST *v)
{
  ml_free(v->array);
  ml_free(v);
}


INLINE void PREFIX(reset) (LIST *lst)
{
  lst->length = 0;
}

int PREFIX(_realloc_array) (LIST *lst, SIZE_T sz);

INLINE int PREFIX(makeroom) (LIST *lst, SIZE_T sz)
{
  if (sz <= lst->allocated)
    return 0;
  else
    return PREFIX(_realloc_array) (lst, sz);
}

INLINE int PREFIX(append) (LIST *lst, VALUE_T x)
{
  if (PREFIX(makeroom) (lst, lst->length + 1) != 0)
    return -1;
  lst->array[(lst->length)++] = x;
  return 0;
}

INLINE VALUE_T PREFIX(poplast) (LIST *lst)
{
  claim(lst->length > 0);
  return lst->array[--(lst->length)];
}

INLINE int PREFIX(insert) (LIST *lst, SIZE_T i, VALUE_T x)
{
  SIZE_T n;
  claim(0 <= i && i <= lst->length);
  if (PREFIX(makeroom)(lst, lst->length + 1) != 0)
    return -1;
  n = lst->length - i;
  lst->length++;
  memmove(lst->array + i+1, lst->array + i, n * sizeof(VALUE_T));
  lst->array[i] = x;
  return 0;
}

INLINE VALUE_T PREFIX(delete) (LIST *lst, SIZE_T i)
{
  VALUE_T x;
  SIZE_T n;
  claim(0 <= i && i < lst->length);
  x = lst->array[i];
  lst->length--;
  n = lst->length - i;
  memmove(lst->array + i, lst->array + i+1, n * sizeof(VALUE_T));
  return x;
}

INLINE VALUE_T PREFIX(fastdelete) (LIST *lst, SIZE_T i)
{
  VALUE_T x;
  claim(0 <= i && i < lst->length);
  x = lst->array[i];
  lst->array[i] = lst->array[lst->length - 1];
  (lst->length)--;
  return x;
}

INLINE int PREFIX(extend) (LIST * dst, LIST * src)
{
  SIZE_T dlen = dst->length;
  SIZE_T slen = src->length;
  if (PREFIX(makeroom)(dst, dst->length + src->length) != 0)
    return -1;
  memmove(dst->array + dlen, src->array, slen * sizeof(VALUE_T));
  return 0;
}

INLINE int PREFIX(copy) (LIST *dst, LIST *src)
{
  if (PREFIX(makeroom)(dst, src->length) != 0)
    return -1;
  dst->length = src->length;
  memcpy(dst->array, src->array, dst->length * sizeof(VALUE_T));
  return 0;
}

INLINE LIST * PREFIX(new_copy) (LIST *lst)
{
  LIST *res = PREFIX(new) (lst->length);
  if (res == NULL)
    return NULL;
  res->length = lst->length;
  memcpy(res->array, lst->array, res->length * sizeof(VALUE_T));
  return res;
}

int PREFIX(reverse) (LIST *dst, LIST *src);

/* eof */
