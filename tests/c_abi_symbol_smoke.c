#include <assert.h>
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

ivector *iv_new_init(uint32_t length, ...);
void iv_free(ivector *v);
ilist *il_new_init(size_t sz, size_t count, ...);
void il_free(ilist *lst);
ivlist *ivl_new_init(size_t sz, size_t count, ...);
void ivl_free_all(ivlist *lst);

int main(void)
{
  ivector *v = iv_new_init(10, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10);
  assert(v != 0);
  assert(v->length == 10);
  assert(v->array[0] == 1);
  assert(v->array[8] == 9);
  assert(v->array[9] == 10);

  ilist *lst = il_new_init(1, 9, 1, 2, 3, 4, 5, 6, 7, 8, 9);
  assert(lst != 0);
  assert(lst->length == 9);
  assert(lst->array[8] == 9);

  ivlist *vl = ivl_new_init(1, 9,
      iv_new_init(1, 1), iv_new_init(1, 2), iv_new_init(1, 3),
      iv_new_init(1, 4), iv_new_init(1, 5), iv_new_init(1, 6),
      iv_new_init(1, 7), iv_new_init(1, 8), iv_new_init(1, 9));
  assert(vl != 0);
  assert(vl->length == 9);
  assert(vl->array[8]->array[0] == 9);

  ivl_free_all(vl);
  il_free(lst);
  iv_free(v);
  return 0;
}
