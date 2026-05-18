#include <assert.h>
#include <stddef.h>
#include <stdint.h>

#include "lrcalc/ilist.h"
#include "lrcalc/ivector.h"
#include "lrcalc/ivlincomb.h"
#include "lrcalc/ivlist.h"
#include "lrcalc/lrcoef.h"
#include "lrcalc/lriter.h"
#include "lrcalc/maple.h"
#include "lrcalc/optshape.h"
#include "lrcalc/part.h"
#include "lrcalc/perm.h"
#include "lrcalc/schublib.h"
#include "lrcalc/schur.h"
#include "lrcalc/vectarg.h"

typedef struct {
  ivlincomb *ht;
  size_t index;
  size_t i;
} expected_ivlc_iter;

int main(void)
{
  _Static_assert(sizeof(ivlc_iter) == sizeof(expected_ivlc_iter), "ivlc_iter size");
  _Static_assert(offsetof(ivlc_iter, i) == offsetof(expected_ivlc_iter, i), "ivlc_iter i");

  ivector *one = iv_new_init(1, 1);
  ivector *two_one = iv_new_init(2, 2, 1);
  ivector *three_two_one = iv_new_init(3, 3, 2, 1);
  ivector *inner = iv_new_init(1, 1);
  ivector *perm = iv_new_init(3, 2, 1, 3);
  assert(one != 0 && two_one != 0 && three_two_one != 0 && inner != 0);
  assert(part_valid(two_one) == 1);
  assert(perm != 0 && perm_valid(perm) == 1);

  ilist *left_list = il_new_init(1, 1, 7);
  ilist *right_list = il_new_init(1, 2, 11, 13);
  assert(left_list != 0 && right_list != 0);
  assert(il_extend(left_list, right_list) == 0);
  assert(left_list->length == 3);
  assert(left_list->array[0] == 7 && left_list->array[1] == 11 && left_list->array[2] == 13);
  il_free(right_list);
  il_free(left_list);

  ivlist *vl1 = ivl_new_init(1, 1, iv_new_init(1, 3));
  ivlist *vl2 = ivl_new_init(1, 1, iv_new_init(1, 4));
  assert(vl1 != 0 && vl2 != 0);
  assert(ivl_extend(vl1, vl2) == 0);
  assert(vl1->length == 2);
  ivl_free(vl2);
  ivl_free_all(vl1);

  ivlincomb *mult = schur_mult(one, one, -1, -1, -1);
  assert(mult != 0);
  assert(ivlc_card(mult) > 0);
  ivlc_iter itr;
  ivlc_first(mult, &itr);
  assert(ivlc_good(&itr) == 1);
  ivlc_free_all(mult);

  ivlincomb *skew = schur_skew(three_two_one, inner, -1, -1);
  assert(skew != 0);
  assert(ivlc_card(skew) > 0);
  ivlc_free_all(skew);

  ivlincomb *coprod = schur_coprod(two_one, 1, 1, -1, 0);
  assert(coprod != 0);
  ivlc_free_all(coprod);

  lrtab_iter *lrit = lrit_new(two_one, inner, 0, -1, -1, -1);
  assert(lrit != 0);
  assert(lrit_good(lrit) == 1);
  lrit_free(lrit);

  ivector *schub_left = iv_new_init(2, 2, 1);
  ivector *schub_right = iv_new_init(2, 2, 1);
  assert(schub_left != 0 && schub_right != 0);
  ivlincomb *schub = mult_schubert(schub_left, schub_right, 0);
  assert(schub != 0);
  ivlc_free_all(schub);
  iv_free(schub_right);
  iv_free(schub_left);

  skew_shape ss = {0, 0, 0, 0};
  assert(optim_skew(&ss, three_two_one, inner, 0, -1) == 0);
  sksh_dealloc(&ss);

  assert(lrcoef_count(three_two_one, inner, one) >= 0);

  iv_free(inner);
  iv_free(perm);
  iv_free(three_two_one);
  iv_free(two_one);
  iv_free(one);
  return 0;
}
