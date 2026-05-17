#include <stdint.h>
#include <stdlib.h>

/* Native LR coefficient kernels, kept C-shaped to match Buch's hot paths. */

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

typedef struct {
  int32_t *storage;
  int32_t *outer;
  int32_t *inner;
  int32_t *content;
  uintptr_t outer_len;
  uintptr_t inner_len;
  uintptr_t content_len;
} lrcoef_optimized_shape;

static int valid_partition(const int32_t *part, uintptr_t len)
{
  uintptr_t i;
  int32_t previous;

  previous = 0;
  for (i = len; i-- > 0;) {
    int32_t value = part[i];
    if (value < previous || value < 0)
      return 0;
    previous = value;
  }
  return 1;
}

static uintptr_t part_length(const int32_t *part, uintptr_t len)
{
  while (len > 0 && part[len - 1] == 0)
    len--;
  return len;
}

static int part_length_i32(const int32_t *part, uintptr_t len, int32_t *out)
{
  uintptr_t trimmed = part_length(part, len);
  if (trimmed > (uintptr_t)INT32_MAX)
    return -1;
  *out = (int32_t)trimmed;
  return 0;
}

static int part_sum_i32(const int32_t *part, uintptr_t len, int32_t *out)
{
  uintptr_t i;
  int64_t sum;

  sum = 0;
  for (i = 0; i < len; i++) {
    if (part[i] < 0)
      return -1;
    sum += part[i];
    if (sum > INT32_MAX)
      return -1;
  }
  *out = (int32_t)sum;
  return 0;
}

static void free_optimized_shape(lrcoef_optimized_shape *shape)
{
  free(shape->storage);
  shape->storage = NULL;
  shape->outer = NULL;
  shape->inner = NULL;
  shape->content = NULL;
}

/*
 * Native coefficient compactification follows Buch's GPL lrcalc optshape.c,
 * but uses raw slices so lrcoef_i64 can avoid Rust-side Vec setup.
 */
static int optim_coef_raw(
    lrcoef_optimized_shape *shape,
    const int32_t *out,
    uintptr_t out_len,
    const int32_t *sh1,
    uintptr_t sh1_len,
    const int32_t *sh2,
    uintptr_t sh2_len)
{
  int32_t *storage, *la, *mu, *nu;
  int32_t N, Nla, Nmu, r, s, N0, nu0, la0, mu0, nur, lar, mur;
  int32_t lar1, mur1, nur1, c, ca, Inu, Ila, Imu;
  int64_t sum;

  shape->storage = NULL;
  shape->outer = NULL;
  shape->inner = NULL;
  shape->content = NULL;
  shape->outer_len = 0;
  shape->inner_len = 0;
  shape->content_len = 0;

  if (!valid_partition(out, out_len) ||
      !valid_partition(sh1, sh1_len) ||
      !valid_partition(sh2, sh2_len))
    return -1;
  if (out_len > (uintptr_t)INT32_MAX ||
      sh1_len > (uintptr_t)INT32_MAX ||
      sh2_len > (uintptr_t)INT32_MAX)
    return -1;
  if (part_length_i32(out, out_len, &N) != 0)
    return -1;
  if ((uintptr_t)N < sh1_len && sh1[N] > 0)
    return 0;
  if ((uintptr_t)N < sh2_len && sh2[N] > 0)
    return 0;
  if (N == 0)
    return 1;

  if ((uintptr_t)N > UINTPTR_MAX / (3 * sizeof(int32_t)))
    return -1;
  storage = (int32_t *)malloc((uintptr_t)N * 3 * sizeof(int32_t));
  if (storage == NULL)
    return -1;
  nu = storage;
  la = storage + N;
  mu = la + N;

  sum = 0;
  for (r = N - 1; r >= 0; r--) {
    nu[r] = out[r];
    sum += nu[r];
    if (sum > INT32_MAX)
      goto coef_error;
  }

  for (Nla = N; Nla > (int32_t)sh1_len; Nla--)
    la[Nla - 1] = 0;
  for (; Nla > 0 && sh1[Nla - 1] == 0; Nla--)
    la[Nla - 1] = 0;
  for (r = Nla - 1; r >= 0; r--) {
    int32_t x = sh1[r];
    la[r] = x;
    if (nu[r] < x)
      goto coef_zero;
    sum -= la[r];
    if (sum < -((int64_t)INT32_MAX) - 1)
      goto coef_error;
  }

  for (Nmu = N; Nmu > (int32_t)sh2_len; Nmu--)
    mu[Nmu - 1] = 0;
  for (; Nmu > 0 && sh2[Nmu - 1] == 0; Nmu--)
    mu[Nmu - 1] = 0;
  for (r = Nmu - 1; r >= 0; r--) {
    int32_t x = sh2[r];
    mu[r] = x;
    if (nu[r] < x)
      goto coef_zero;
    sum -= mu[r];
    if (sum < -((int64_t)INT32_MAX) - 1)
      goto coef_error;
  }

  if (sum != 0)
    goto coef_zero;

  N0 = N + 1;
  nu0 = 0;
  while (N < N0 || nu[0] < nu0) {
    N0 = N;
    nu0 = nu[0];

    /* Horizontal compactification of nu/la. */
    mu0 = mu[0];
    lar1 = 0;
    nur1 = 0;
    for (r = N - 1; r >= 0; r--) {
      lar = la[r];
      nur = nu[r];
      if (lar > nur1 || nur - lar1 > mu0)
        break;
      lar1 = lar;
      nur1 = nur;
    }
    c = 0;
    for (; r >= 0; r--) {
      lar = la[r];
      nur = nu[r];
      if (nur - lar > mu0)
        goto coef_zero;
      ca = nur - lar1 - mu0;
      if (ca < lar - nur1)
        ca = lar - nur1;
      if (ca > 0)
        c += ca;
      if (nur - c < mu[r])
        goto coef_zero;
      if (nur == c) {
        N = r;
        break;
      }
      la[r] = lar - c;
      nu[r] = nur - c;
      lar1 = lar;
      nur1 = nur;
    }

    /* Remove row of size mu[0] from nu/la. */
    r = 0;
    while (r < N && nu[r] - la[r] < mu0)
      r++;
    if (r < N) {
      if (nu[r] - la[r] > mu0)
        goto coef_zero;
      for (; r < N - 1; r++) {
        la[r] = la[r + 1];
        nu[r] = nu[r + 1];
      }
      for (r = 0; r < N - 1; r++)
        mu[r] = mu[r + 1];
      N -= 1;
    }

    /* Horizontal compactification of nu/mu. */
    la0 = la[0];
    mur1 = 0;
    nur1 = 0;
    for (r = N - 1; r >= 0; r--) {
      mur = mu[r];
      nur = nu[r];
      if (mur > nur1 || nur - mur1 > la0)
        break;
      mur1 = mur;
      nur1 = nur;
    }
    c = 0;
    for (; r >= 0; r--) {
      mur = mu[r];
      nur = nu[r];
      if (nur - mur > la0)
        goto coef_zero;
      ca = nur - mur1 - la0;
      if (ca < mur - nur1)
        ca = mur - nur1;
      if (ca > 0)
        c += ca;
      if (nur - c < la[r])
        goto coef_zero;
      if (nur == c) {
        N = r;
        break;
      }
      mu[r] = mur - c;
      nu[r] = nur - c;
      mur1 = mur;
      nur1 = nur;
    }

    /* Remove row of size la[0] from nu/mu. */
    r = 0;
    while (r < N && nu[r] - mu[r] < la0)
      r++;
    if (r < N) {
      if (nu[r] - mu[r] > la0)
        goto coef_zero;
      for (; r < N - 1; r++) {
        mu[r] = mu[r + 1];
        nu[r] = nu[r + 1];
      }
      for (r = 0; r < N - 1; r++)
        la[r] = la[r + 1];
      N -= 1;
    }

    /* Vertical compactification of nu/la. */
    if (N < Nmu)
      Nmu = N;
    while (Nmu > 0 && mu[Nmu - 1] == 0)
      Nmu--;
    if (Nmu == 0)
      goto coef_one;
    r = 0;
    while (r < Nmu && la[r] < nu[r])
      r++;
    while (r < N && la[r] < nu[r] && nu[r] < la[r - Nmu])
      r++;
    if (r < N) {
      Inu = r;
      s = (r > Nmu) ? (r - Nmu) : 0;
      Ila = s;
      for (; r < N && Inu < Nmu; r++) {
        if (la[r] == nu[r]) {
          la[r] = -1;
        } else {
          nu[Inu] = nu[r];
          if (nu[Inu] < mu[Inu])
            goto coef_zero;
          Inu++;
        }
      }
      while (r < N) {
        if (la[r] == nu[r]) {
          la[r] = -1;
          r++;
          continue;
        }
        while (la[s] == -1)
          s++;
        if (la[s] < nu[r])
          goto coef_zero;
        if (la[s] > nu[r]) {
          la[Ila] = la[s];
          Ila++;
          nu[Inu] = nu[r];
          if (nu[Inu] < mu[Inu])
            goto coef_zero;
          Inu++;
        }
        r++;
        s++;
      }
      while (s < N) {
        if (la[s] != -1) {
          la[Ila] = la[s];
          Ila++;
        }
        s++;
      }
      if (Inu < N && mu[Inu] > 0)
        goto coef_zero;
      N = Inu;
    }

    /* Remove column of size len(mu) from nu/la. */
    r = Nmu;
    while (r <= N && nu[r - 1] <= la[r - Nmu])
      r++;
    if (r <= N) {
      if (r > Nmu && nu[r - 1] > la[r - Nmu - 1])
        goto coef_zero;
      if (r < N && nu[r] > la[r - Nmu])
        goto coef_zero;
      for (s = r - Nmu - 1; s >= 0; s--)
        la[s]--;
      for (s = Nmu - 1; s >= 0; s--)
        mu[s]--;
      for (s = 0; s < r; s++) {
        nu[s]--;
        if (nu[s] == 0) {
          N = s;
          break;
        }
      }
    }

    /* Vertical compactification of nu/mu. */
    if (N < Nla)
      Nla = N;
    while (Nla > 0 && la[Nla - 1] == 0)
      Nla--;
    if (Nla == 0)
      goto coef_one;
    r = 0;
    while (r < Nla && mu[r] < nu[r])
      r++;
    while (r < N && mu[r] < nu[r] && nu[r] < mu[r - Nla])
      r++;
    if (r < N) {
      Inu = r;
      s = (r > Nla) ? (r - Nla) : 0;
      Imu = s;
      for (; r < N && Inu < Nla; r++) {
        if (mu[r] == nu[r]) {
          mu[r] = -1;
        } else {
          nu[Inu] = nu[r];
          if (nu[Inu] < la[Inu])
            goto coef_zero;
          Inu++;
        }
      }
      while (r < N) {
        if (mu[r] == nu[r]) {
          mu[r] = -1;
          r++;
          continue;
        }
        while (mu[s] == -1)
          s++;
        if (mu[s] < nu[r])
          goto coef_zero;
        if (mu[s] > nu[r]) {
          mu[Imu] = mu[s];
          Imu++;
          nu[Inu] = nu[r];
          if (nu[Inu] < la[Inu])
            goto coef_zero;
          Inu++;
        }
        r++;
        s++;
      }
      while (s < N) {
        if (mu[s] != -1) {
          mu[Imu] = mu[s];
          Imu++;
        }
        s++;
      }
      if (Inu < N && la[Inu] > 0)
        goto coef_zero;
      N = Inu;
    }

    /* Remove column of size len(la) from nu/mu. */
    r = Nla;
    while (r <= N && nu[r - 1] <= mu[r - Nla])
      r++;
    if (r <= N) {
      if (r > Nla && nu[r - 1] > mu[r - Nla - 1])
        goto coef_zero;
      if (r < N && nu[r] > mu[r - Nla])
        goto coef_zero;
      for (s = r - Nla - 1; s >= 0; s--)
        mu[s]--;
      for (s = Nla - 1; s >= 0; s--)
        la[s]--;
      for (s = 0; s < r; s++) {
        nu[s]--;
        if (nu[s] == 0) {
          N = s;
          break;
        }
      }
    }
  }

  if (N == 0)
    goto coef_one;

  if (N < Nla)
    Nla = N;
  while (Nla > 0 && la[Nla - 1] == 0)
    Nla--;
  if (N < Nmu)
    Nmu = N;
  while (Nmu > 0 && mu[Nmu - 1] == 0)
    Nmu--;

  shape->outer = nu;
  shape->inner = la;
  shape->content = mu;
  shape->storage = storage;
  shape->outer_len = (uintptr_t)N;
  shape->inner_len = (uintptr_t)Nla;
  shape->content_len = (uintptr_t)Nmu;
  return 2;

coef_one:
  free(storage);
  return 1;

coef_zero:
  free(storage);
  return 0;

coef_error:
  free(storage);
  return -1;
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

int64_t lrcalc_native_lrcoef_i64(
    const int32_t *outer,
    uintptr_t outer_len,
    const int32_t *inner1,
    uintptr_t inner1_len,
    const int32_t *inner2,
    uintptr_t inner2_len)
{
  lrcoef_optimized_shape shape;
  int32_t content_sum;
  int64_t coef;
  int sign;

  if ((outer == NULL && outer_len != 0) ||
      (inner1 == NULL && inner1_len != 0) ||
      (inner2 == NULL && inner2_len != 0))
    return -1;

  sign = optim_coef_raw(&shape, outer, outer_len, inner1, inner1_len, inner2, inner2_len);
  if (sign <= 1)
    return sign;

  if (part_sum_i32(shape.content, shape.content_len, &content_sum) != 0) {
    free_optimized_shape(&shape);
    return -1;
  }
  if (content_sum <= 1) {
    free_optimized_shape(&shape);
    return 1;
  }

  coef = lrcalc_native_lrcoef_count_i64(
      shape.outer,
      shape.outer_len,
      shape.inner,
      shape.inner_len,
      shape.content,
      shape.content_len,
      content_sum);
  free_optimized_shape(&shape);
  return coef;
}
