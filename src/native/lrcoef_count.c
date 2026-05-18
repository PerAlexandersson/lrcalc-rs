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
  int32_t value;
  int32_t max;
  int32_t above;
  int32_t right;
} lrit_content_box;

typedef struct {
  int32_t cont;
  int32_t supply;
} lrcoef_content;

typedef int32_t (*lrcalc_content_begin_fn)(void *ctx, uintptr_t term_count);
typedef int32_t (*lrcalc_content_emit_fn)(
    void *ctx,
    const int32_t *content,
    uintptr_t content_len,
    int64_t coefficient);

static inline int32_t part_entry(const int32_t *part, uintptr_t len, uintptr_t index)
{
  return (index < len) ? part[index] : 0;
}

static uintptr_t part_length(const int32_t *part, uintptr_t len);

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

typedef __uint128_t packed_key;

typedef struct {
  packed_key *keys;
  uint64_t *values;
  uintptr_t capacity;
  uintptr_t len;
  uintptr_t resize_at;
} packed_content_table;

typedef struct {
  uint32_t bits;
  uint32_t len_bits;
  packed_key len_mask;
  packed_key value_mask;
  uintptr_t max_len;
  uintptr_t len;
  uintptr_t overflow_labels;
  packed_key key;
} packed_content_state;

static inline uint64_t mix_u64(uint64_t value)
{
  value ^= value >> 33;
  value *= UINT64_C(0xff51afd7ed558ccd);
  value ^= value >> 33;
  value *= UINT64_C(0xc4ceb9fe1a85ec53);
  return value ^ (value >> 33);
}

static inline uint64_t mix_key(packed_key value)
{
  uint64_t low = (uint64_t)value;
  uint64_t high = (uint64_t)(value >> 64);
  return mix_u64(low ^ ((high << 32) | (high >> 32)));
}

static uint32_t bits_needed_u64(uint64_t value)
{
  uint32_t bits = 0;

  do {
    bits++;
    value >>= 1;
  } while (value != 0);
  return bits;
}

static packed_key mask_bits_u128(uint32_t bits)
{
  if (bits >= 128)
    return ~(packed_key)0;
  return (((packed_key)1) << bits) - 1;
}

static uintptr_t resize_threshold(uintptr_t capacity)
{
  uintptr_t threshold = capacity / 2;
  return threshold == 0 ? 1 : threshold;
}

static uintptr_t next_power_of_two(uintptr_t value)
{
  uintptr_t power = 2;

  while (power < value) {
    if (power > UINTPTR_MAX / 2)
      return 0;
    power *= 2;
  }
  return power;
}

static uintptr_t packed_initial_capacity(uintptr_t skew_size, uintptr_t label_count)
{
  uintptr_t cap;

  if (skew_size != 0 && label_count > UINTPTR_MAX / skew_size)
    cap = 65536;
  else {
    cap = skew_size * label_count;
    if (cap > UINTPTR_MAX / 8)
      cap = 65536;
    else
      cap *= 8;
  }

  if (cap < 64)
    cap = 64;
  if (cap > 65536)
    cap = 65536;
  return next_power_of_two(cap);
}

static int packed_table_init(
    packed_content_table *table,
    uintptr_t skew_size,
    uintptr_t label_count)
{
  uintptr_t capacity = packed_initial_capacity(skew_size, label_count);

  if (capacity == 0)
    return -1;
  table->keys = (packed_key *)calloc(capacity, sizeof(packed_key));
  table->values = (uint64_t *)calloc(capacity, sizeof(uint64_t));
  if (table->keys == NULL || table->values == NULL) {
    free(table->keys);
    free(table->values);
    table->keys = NULL;
    table->values = NULL;
    return -1;
  }
  table->capacity = capacity;
  table->len = 0;
  table->resize_at = resize_threshold(capacity);
  return 0;
}

static void packed_table_dealloc(packed_content_table *table)
{
  free(table->keys);
  free(table->values);
  table->keys = NULL;
  table->values = NULL;
  table->capacity = 0;
  table->len = 0;
  table->resize_at = 0;
}

static void packed_table_insert_existing(
    packed_content_table *table,
    packed_key key,
    uint64_t value)
{
  uintptr_t mask = table->capacity - 1;
  uintptr_t index = (uintptr_t)mix_key(key) & mask;

  while (table->values[index] != 0)
    index = (index + 1) & mask;

  table->keys[index] = key;
  table->values[index] = value;
  table->len++;
}

static int packed_table_grow(packed_content_table *table)
{
  packed_key *old_keys = table->keys;
  uint64_t *old_values = table->values;
  uintptr_t old_capacity = table->capacity;
  uintptr_t new_capacity;
  uintptr_t i;

  if (old_capacity > UINTPTR_MAX / 2)
    return -1;
  new_capacity = old_capacity * 2;
  table->keys = (packed_key *)calloc(new_capacity, sizeof(packed_key));
  table->values = (uint64_t *)calloc(new_capacity, sizeof(uint64_t));
  if (table->keys == NULL || table->values == NULL) {
    free(table->keys);
    free(table->values);
    table->keys = old_keys;
    table->values = old_values;
    return -1;
  }

  table->capacity = new_capacity;
  table->len = 0;
  table->resize_at = resize_threshold(new_capacity);
  for (i = 0; i < old_capacity; i++) {
    if (old_values[i] != 0)
      packed_table_insert_existing(table, old_keys[i], old_values[i]);
  }

  free(old_keys);
  free(old_values);
  return 0;
}

static int packed_table_add(packed_content_table *table, packed_key key)
{
  while (1) {
    uintptr_t mask = table->capacity - 1;
    uintptr_t index = (uintptr_t)mix_key(key) & mask;

    while (1) {
      uint64_t value = table->values[index];
      if (value == 0) {
        if (table->len >= table->resize_at) {
          if (packed_table_grow(table) != 0)
            return -1;
          break;
        }
        table->keys[index] = key;
        table->values[index] = 1;
        table->len++;
        return 0;
      }
      if (table->keys[index] == key) {
        if (value == (uint64_t)INT64_MAX)
          return -2;
        table->values[index] = value + 1;
        return 0;
      }
      index = (index + 1) & mask;
    }
  }
}

static int packed_state_init(
    packed_content_state *state,
    int32_t skew_size,
    uintptr_t label_count)
{
  uint32_t bits;
  uint32_t len_bits;
  uintptr_t max_len;

  if (skew_size < 0)
    return -1;
  bits = bits_needed_u64((uint64_t)skew_size);
  len_bits = bits_needed_u64((uint64_t)label_count);
  if (len_bits >= 128)
    return -2;
  max_len = (uintptr_t)((128 - len_bits) / bits);
  if (max_len == 0)
    return -2;

  state->bits = bits;
  state->len_bits = len_bits;
  state->len_mask = mask_bits_u128(len_bits);
  state->value_mask = mask_bits_u128(bits);
  state->max_len = max_len;
  state->len = 0;
  state->overflow_labels = 0;
  state->key = 0;
  return 0;
}

static inline packed_key packed_state_get(const packed_content_state *state, uintptr_t label)
{
  uint32_t shift = state->len_bits + state->bits * (uint32_t)(label - 1);
  return (state->key >> shift) & state->value_mask;
}

static inline void packed_state_set(
    packed_content_state *state,
    uintptr_t label,
    packed_key value)
{
  uint32_t shift = state->len_bits + state->bits * (uint32_t)(label - 1);
  packed_key mask = state->value_mask << shift;
  state->key = (state->key & ~mask) | (value << shift);
}

static inline void packed_state_write_len(packed_content_state *state)
{
  state->key = (state->key & ~state->len_mask) | (packed_key)state->len;
}

static inline void packed_state_place(packed_content_state *state, uintptr_t label)
{
  packed_key old;

  if (label > state->max_len) {
    state->overflow_labels++;
    return;
  }
  old = packed_state_get(state, label);
  packed_state_set(state, label, old + 1);
  if (label > state->len) {
    state->len = label;
    packed_state_write_len(state);
  }
}

static inline void packed_state_unplace(packed_content_state *state, uintptr_t label)
{
  packed_key old;

  if (label > state->max_len) {
    state->overflow_labels--;
    return;
  }
  old = packed_state_get(state, label);
  packed_state_set(state, label, old - 1);
  if (label == state->len && old == 1) {
    while (state->len > 0 && packed_state_get(state, state->len) == 0)
      state->len--;
    packed_state_write_len(state);
  }
}

static int unpack_packed_key(
    packed_key key,
    const packed_content_state *state,
    int32_t *content,
    uintptr_t *content_len);

static lrit_content_box *new_lrit_content_skewtab(
    const int32_t *outer,
    uintptr_t outer_len,
    const int32_t *inner,
    uintptr_t inner_len,
    uintptr_t beta_len,
    uintptr_t label_count,
    int32_t skew_size)
{
  lrit_content_box *array;
  uintptr_t trimmed_outer;
  uintptr_t trimmed_inner;
  uintptr_t array_len;
  int32_t out0;
  int32_t inn0;
  int32_t out1;
  int32_t inn1;
  int32_t out2;
  int32_t maxdepth;
  int32_t s;
  int32_t r;

  trimmed_outer = part_length(outer, outer_len);
  trimmed_inner = part_length(inner, inner_len);
  if (trimmed_inner > trimmed_outer)
    trimmed_inner = trimmed_outer;
  if (trimmed_outer > (uintptr_t)INT32_MAX ||
      trimmed_inner > (uintptr_t)INT32_MAX ||
      beta_len > (uintptr_t)INT32_MAX ||
      label_count > (uintptr_t)INT32_MAX)
    return NULL;
  if (skew_size < 0 || (uintptr_t)skew_size > (UINTPTR_MAX / sizeof(lrit_content_box)) - 2)
    return NULL;

  array_len = (uintptr_t)skew_size + 2;
  array = (lrit_content_box *)malloc(array_len * sizeof(lrit_content_box));
  if (array == NULL)
    return NULL;

  maxdepth = (int32_t)beta_len;
  for (r = 0; r < (int32_t)trimmed_outer; r++) {
    int32_t rowsz = outer[r] - part_entry(inner, trimmed_inner, (uintptr_t)r);
    if (rowsz > 0)
      maxdepth++;
  }
  if ((int32_t)label_count > maxdepth)
    label_count = (uintptr_t)maxdepth;

  s = 0;
  out1 = 0;
  out0 = (trimmed_outer == 0) ? 0 : outer[0];
  inn0 = (trimmed_outer == 0) ? out0 :
      (trimmed_outer <= trimmed_inner ? inner[trimmed_outer - 1] : 0);
  for (r = (int32_t)trimmed_outer; r-- > 0;) {
    int32_t c;

    out2 = out1;
    inn1 = inn0;
    out1 = outer[r];
    inn0 = (r == 0) ? out0 :
        ((uintptr_t)r <= trimmed_inner ? inner[(uintptr_t)r - 1] : 0);
    if (inn1 < out1)
      maxdepth--;
    for (c = inn1; c < out1; c++) {
      lrit_content_box *box = array + s;
      int32_t max_value;

      box->right = (c + 1 < out1) ? s + 1 : skew_size + 1;
      box->above = (c >= inn0) ? s + out1 - inn0 : skew_size;
      max_value = (c < out2) ? array[s - out2 + inn1].max - 1 : (int32_t)label_count - 1;
      box->max = (max_value < maxdepth) ? max_value : maxdepth;
      box->value = 0;
      s++;
    }
  }
  if (s != skew_size) {
    free(array);
    return NULL;
  }
  array[skew_size].value = -1;
  array[skew_size + 1].value = (int32_t)label_count - 1;
  return array;
}

static int lrit_content_minimal_fill(
    lrit_content_box *array,
    int32_t *content,
    packed_content_state *state,
    int32_t skew_size)
{
  int32_t s;

  for (s = skew_size; s-- > 0;) {
    lrit_content_box *box = array + s;
    int32_t x = array[box->above].value + 1;
    if (x > box->max)
      return 0;
    box->value = x;
    content[x]++;
    packed_state_place(state, (uintptr_t)x + 1);
  }
  return 1;
}

static int lrit_content_next(
    lrit_content_box *array,
    int32_t *content,
    packed_content_state *state,
    int32_t skew_size)
{
  lrit_content_box *box;
  lrit_content_box *box_bound = array + skew_size;

  for (box = array; box != box_bound; box++) {
    int32_t max_value = array[box->right].value;
    int32_t x;

    if (max_value > box->max)
      max_value = box->max;
    x = box->value;
    content[x]--;
    packed_state_unplace(state, (uintptr_t)x + 1);
    x++;
    while (x <= max_value && content[x] == content[x - 1])
      x++;
    if (x > max_value)
      continue;

    box->value = x;
    content[x]++;
    packed_state_place(state, (uintptr_t)x + 1);
    while (box != array) {
      box--;
      x = array[box->above].value + 1;
      box->value = x;
      content[x]++;
      packed_state_place(state, (uintptr_t)x + 1);
    }
    return 1;
  }

  return 0;
}

static int32_t lrcalc_native_beta_content_expand_lrit_i64(
    const int32_t *outer,
    uintptr_t outer_len,
    const int32_t *inner,
    uintptr_t inner_len,
    const int32_t *beta,
    uintptr_t beta_len,
    uintptr_t label_count,
    int32_t skew_size,
    lrcalc_content_begin_fn begin,
    lrcalc_content_emit_fn emit,
    void *ctx)
{
  uintptr_t trimmed_beta_len;
  uintptr_t effective_labels;
  uintptr_t i;
  lrit_content_box *array = NULL;
  int32_t *prefix_content = NULL;
  int32_t *content = NULL;
  packed_content_table terms;
  packed_content_state state;
  int32_t status;

  terms.keys = NULL;
  terms.values = NULL;
  terms.capacity = 0;
  terms.len = 0;
  terms.resize_at = 0;

  if (skew_size > 64)
    return -3;

  trimmed_beta_len = part_length(beta, beta_len);
  if (trimmed_beta_len > label_count)
    return -3;

  effective_labels = label_count;
  if (effective_labels < trimmed_beta_len)
    return -3;
  if (effective_labels == 0)
    return -3;

  if (packed_state_init(&state, skew_size, effective_labels) != 0)
    return -2;
  if (packed_table_init(&terms, (uintptr_t)skew_size, effective_labels) != 0)
    return -1;

  prefix_content = (int32_t *)calloc(effective_labels, sizeof(int32_t));
  if (prefix_content == NULL) {
    packed_table_dealloc(&terms);
    return -1;
  }
  for (i = 0; i < effective_labels; i++)
    prefix_content[i] = part_entry(beta, beta_len, i);

  array = new_lrit_content_skewtab(
      outer,
      outer_len,
      inner,
      inner_len,
      trimmed_beta_len,
      effective_labels,
      skew_size);
  if (array == NULL) {
    status = -1;
    goto cleanup;
  }

  if (lrit_content_minimal_fill(array, prefix_content, &state, skew_size)) {
    do {
      if (state.overflow_labels != 0) {
        status = -2;
        goto cleanup;
      }
      status = packed_table_add(&terms, state.key);
      if (status != 0)
        goto cleanup;
    } while (lrit_content_next(array, prefix_content, &state, skew_size));
  }

  content = (int32_t *)malloc(state.max_len * sizeof(int32_t));
  if (content == NULL && state.max_len != 0) {
    status = -1;
    goto cleanup;
  }
  if (begin(ctx, terms.len) != 0) {
    status = -1;
    goto cleanup;
  }
  for (i = 0; i < terms.capacity; i++) {
    uintptr_t content_len;

    if (terms.values[i] == 0)
      continue;
    if (unpack_packed_key(terms.keys[i], &state, content, &content_len) != 0) {
      status = -1;
      goto cleanup;
    }
    if (emit(ctx, content, content_len, (int64_t)terms.values[i]) != 0) {
      status = -1;
      goto cleanup;
    }
  }
  status = 0;

cleanup:
  free(content);
  free(array);
  free(prefix_content);
  packed_table_dealloc(&terms);
  return status;
}

static int unpack_packed_key(
    packed_key key,
    const packed_content_state *state,
    int32_t *content,
    uintptr_t *content_len)
{
  uintptr_t len = (uintptr_t)(key & state->len_mask);
  uintptr_t i;

  if (len > state->max_len)
    return -1;
  for (i = 0; i < len; i++) {
    uint32_t shift = state->len_bits + state->bits * (uint32_t)i;
    content[i] = (int32_t)((key >> shift) & state->value_mask);
  }
  *content_len = len;
  return 0;
}

static inline int beta_label_allowed(uintptr_t label, const int32_t *slack)
{
  return label == 1 || slack[label] > 0;
}

static inline void beta_place_label(uintptr_t label, uintptr_t label_count, int32_t *slack)
{
  if (label > 1)
    slack[label]--;
  if (label < label_count)
    slack[label + 1]++;
}

static inline void beta_unplace_label(uintptr_t label, uintptr_t label_count, int32_t *slack)
{
  if (label > 1)
    slack[label]++;
  if (label < label_count)
    slack[label + 1]--;
}

int32_t lrcalc_native_beta_content_expand_i64(
    const int32_t *outer,
    uintptr_t outer_len,
    const int32_t *inner,
    uintptr_t inner_len,
    const int32_t *beta,
    uintptr_t beta_len,
    uintptr_t label_count,
    int32_t skew_size,
    lrcalc_content_begin_fn begin,
    lrcalc_content_emit_fn emit,
    void *ctx)
{
  lrcoef_box *T = NULL;
  int32_t *slack = NULL;
  int32_t *content = NULL;
  packed_content_table terms;
  packed_content_state state;
  uintptr_t i;
  uintptr_t n;
  uintptr_t pos;
  int32_t x;
  int32_t above;
  int32_t status;

  terms.keys = NULL;
  terms.values = NULL;
  terms.capacity = 0;
  terms.len = 0;
  terms.resize_at = 0;

  if (outer == NULL || begin == NULL || emit == NULL)
    return -1;
  if (skew_size < 0 || label_count > (uintptr_t)INT32_MAX)
    return -2;
  if (skew_size == 0) {
    if (begin(ctx, 1) != 0)
      return -1;
    return emit(ctx, NULL, 0, 1);
  }
  if (outer_len == 0 || label_count == 0)
    return -1;

  status = lrcalc_native_beta_content_expand_lrit_i64(
      outer,
      outer_len,
      inner,
      inner_len,
      beta,
      beta_len,
      label_count,
      skew_size,
      begin,
      emit,
      ctx);
  if (status != -3)
    return status;

  status = packed_state_init(&state, skew_size, label_count);
  if (status != 0)
    return status;
  if (packed_table_init(&terms, (uintptr_t)skew_size, label_count) != 0)
    return -1;

  slack = (int32_t *)calloc(label_count + 2, sizeof(int32_t));
  if (slack == NULL) {
    packed_table_dealloc(&terms);
    return -1;
  }
  for (i = 2; i <= label_count; i++)
    slack[i] = part_entry(beta, beta_len, i - 2) - part_entry(beta, beta_len, i - 1);

  T = new_skewtab(outer, outer_len, inner, inner_len, (int32_t)label_count, skew_size);
  if (T == NULL) {
    free(slack);
    packed_table_dealloc(&terms);
    return -1;
  }

  n = (uintptr_t)skew_size;
  pos = 0;
  above = T[T[pos].north].value;
  x = (int32_t)label_count;

  while (1) {
    while (x > T[pos].max)
      x--;
    while (x > 0 && x > above && !beta_label_allowed((uintptr_t)x, slack))
      x--;

    if (x <= above) {
      uintptr_t label;

      if (pos == 0)
        break;
      pos--;
      above = T[T[pos].north].value;
      x = T[pos].value;
      label = (uintptr_t)x;
      beta_unplace_label(label, label_count, slack);
      packed_state_unplace(&state, label);
      x--;
    } else if (pos + 1 < n) {
      uintptr_t label = (uintptr_t)x;

      T[pos].value = x;
      beta_place_label(label, label_count, slack);
      packed_state_place(&state, label);
      pos++;
      x = T[T[pos].east].value;
      above = T[T[pos].north].value;
    } else {
      uintptr_t label = (uintptr_t)x;

      T[pos].value = x;
      beta_place_label(label, label_count, slack);
      packed_state_place(&state, label);
      if (state.overflow_labels != 0) {
        status = -2;
        goto cleanup;
      }
      status = packed_table_add(&terms, state.key);
      if (status != 0)
        goto cleanup;
      beta_unplace_label(label, label_count, slack);
      packed_state_unplace(&state, label);
      x--;
    }
  }

  content = (int32_t *)malloc(state.max_len * sizeof(int32_t));
  if (content == NULL && state.max_len != 0) {
    status = -1;
    goto cleanup;
  }
  if (begin(ctx, terms.len) != 0) {
    status = -1;
    goto cleanup;
  }
  for (i = 0; i < terms.capacity; i++) {
    uintptr_t content_len;

    if (terms.values[i] == 0)
      continue;
    if (unpack_packed_key(terms.keys[i], &state, content, &content_len) != 0) {
      status = -1;
      goto cleanup;
    }
    if (emit(ctx, content, content_len, (int64_t)terms.values[i]) != 0) {
      status = -1;
      goto cleanup;
    }
  }
  status = 0;

cleanup:
  free(content);
  free(T);
  free(slack);
  packed_table_dealloc(&terms);
  return status;
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
