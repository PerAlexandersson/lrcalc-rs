#include <stdint.h>
#include <stdlib.h>

/* Count-only LR coefficient kernel, kept C-shaped to match Buch's hot loop. */

typedef struct {
  int32_t value;
  int32_t max;
  int32_t north;
  int32_t east;
  int32_t se_supply;
  int32_t se_sz;
  int32_t west_sz;
  int32_t padding;
} lrcoef_box;

typedef struct {
  int32_t cont;
  int32_t supply;
} lrcoef_content;

static inline int32_t part_entry(const int32_t *part, uintptr_t len, uintptr_t index)
{
  return (index < len) ? part[index] : 0;
}

static lrcoef_content *new_content(const int32_t *content, uintptr_t content_len)
{
  lrcoef_content *counts;
  uintptr_t i;

  counts = (lrcoef_content *)malloc((content_len + 1) * sizeof(lrcoef_content));
  if (counts == NULL)
    return NULL;

  counts[0].cont = content[0];
  counts[0].supply = content[0];
  for (i = 0; i < content_len; i++) {
    counts[i + 1].cont = 0;
    counts[i + 1].supply = content[i];
  }

  return counts;
}

static lrcoef_box *new_skewtab(
    const int32_t *outer,
    uintptr_t outer_len,
    const int32_t *inner,
    uintptr_t inner_len,
    int32_t max_value,
    int32_t skew_size)
{
  lrcoef_box *array;
  uintptr_t n;
  int32_t pos, rr;

  if (skew_size < 0 || outer_len > (uintptr_t)INT32_MAX)
    return NULL;
  n = (uintptr_t)skew_size;
  if (n > (UINTPTR_MAX / sizeof(lrcoef_box)) - 2)
    return NULL;

  array = (lrcoef_box *)malloc((n + 2) * sizeof(lrcoef_box));
  if (array == NULL)
    return NULL;

  pos = skew_size;
  for (rr = (int32_t)outer_len; rr-- > 0;) {
    int32_t nu_0 = (rr == 0) ? outer[0] : outer[rr - 1];
    int32_t la_0 = (rr == 0) ? outer[0] : part_entry(inner, inner_len, rr - 1);
    int32_t nu_r = outer[rr];
    int32_t la_r = part_entry(inner, inner_len, rr);
    int32_t nu_1 = part_entry(outer, outer_len, rr + 1);
    int32_t c;

    for (c = la_r; c < nu_r; c++) {
      lrcoef_box *box;
      int32_t north;
      int32_t east;
      int32_t west_sz;

      if (pos == 0) {
        free(array);
        return NULL;
      }
      pos--;
      box = array + pos;

      if (la_0 <= c && c < nu_0) {
        north = pos - nu_r + la_0;
      } else {
        north = skew_size;
      }

      east = (c + 1 < nu_r) ? pos - 1 : skew_size + 1;

      west_sz = c - la_r;
      box->north = north;
      box->east = east;
      box->west_sz = west_sz;

      if (c >= nu_1) {
        box->max = max_value;
        box->se_sz = 0;
      } else {
        int32_t below = pos + nu_1 - la_r;
        box->max = array[below].max - 1;
        box->se_sz = array[below].se_sz + nu_1 - c;
      }
    }
  }

  if (pos != 0) {
    free(array);
    return NULL;
  }
  array[n].value = 0;
  array[n + 1].value = max_value;
  array[n + 1].se_supply = 0;

  return array;
}

int64_t lrcalc_native_lrcoef_count_i64(
    const int32_t *outer,
    uintptr_t outer_len,
    const int32_t *inner,
    uintptr_t inner_len,
    const int32_t *content,
    uintptr_t content_len,
    int32_t content_sum)
{
  lrcoef_box *T;
  lrcoef_box *box;
  lrcoef_content *C;
  int32_t N;
  int32_t pos;
  int32_t x;
  int32_t above;
  int32_t se_supply;
  int64_t coef;

  if (outer == NULL || content == NULL || outer_len == 0 || content_len == 0)
    return -1;
  if (content_len > (uintptr_t)INT32_MAX)
    return -1;
  if (content_sum < 0)
    return -1;

  T = new_skewtab(outer, outer_len, inner, inner_len, (int32_t)content_len, content_sum);
  if (T == NULL)
    return -1;
  C = new_content(content, content_len);
  if (C == NULL) {
    free(T);
    return -1;
  }

  N = content_sum;
  pos = 0;
  box = T;
  above = T[box->north].value;
  x = 1;
  se_supply = N - C[1].supply;
  coef = 0;

  while (1) {
    while (x > above &&
           (C[x].cont == C[x].supply || C[x].cont == C[x - 1].cont)) {
      se_supply += C[x].supply - C[x].cont;
      x--;
    }

    if (x == above || N - pos - se_supply <= box->west_sz) {
      pos--;
      if (pos < 0)
        break;
      box--;
      se_supply = box->se_supply;
      above = T[box->north].value;
      x = box->value;
      C[x].cont--;
      se_supply += C[x].supply - C[x].cont;
      x--;
    } else if (pos + 1 < N) {
      box->se_supply = se_supply;
      box->value = x;
      C[x].cont++;
      pos++;
      box++;
      se_supply = T[box->east].se_supply;
      x = T[box->east].value;
      above = T[box->north].value;
      while (x > box->max) {
        se_supply += C[x].supply - C[x].cont;
        x--;
      }
      while (x > above && se_supply < box->se_sz) {
        se_supply += C[x].supply - C[x].cont;
        x--;
      }
    } else {
      if (coef == INT64_MAX) {
        coef = -1;
        break;
      }
      coef++;
      pos--;
      if (pos < 0)
        break;
      box--;
      se_supply = box->se_supply;
      above = T[box->north].value;
      x = box->value;
      C[x].cont--;
      se_supply += C[x].supply - C[x].cont;
      x--;
    }
  }

  free(T);
  free(C);
  return coef;
}
