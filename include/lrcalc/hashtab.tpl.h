#include <stdio.h>

#include "alloc.h"


/* This structure allows fast lookup, insertion, and reset.
   Remove is less optimized.  Iterating through the table will be
   slow if lots of elements are first inserted and then removed.
*/

typedef struct {
  KEY_T key;
  VALUE_T value;
  HASH_T hash;
  SIZE_T next;
} PREFIX(keyval_t);

typedef struct {
  SIZE_T *table;    /* Hash table. */
  PREFIX(keyval_t) *elts; /* List of elements indexed by table. */
  SIZE_T card;      /* Number of (key,value) pairs. */
  SIZE_T free_elts; /* First free keyval item, 0 if none. */
  SIZE_T elts_len;  /* Number of itmes used or on free_elts list. */
  SIZE_T elts_sz;   /* Allocated items on elts. */
  SIZE_T table_sz;  /* Allocated hash table size. */
} HASHTAB;


/* Minimal number of table entries for each element. */
#ifndef USE_FACTOR
#define USE_FACTOR 2
#endif

#ifndef INIT_TABLE_SIZE
#define INIT_TABLE_SIZE 2003
#endif

#ifndef INIT_ELT_SIZE
#define INIT_ELT_SIZE 100
#endif


INLINE SIZE_T PREFIX(card) (HASHTAB *ht)
{
  return ht->card;
}


/* Initialize hash table structure. */
INLINE int PREFIX(init) (HASHTAB *ht, SIZE_T tabsz, SIZE_T eltsz)
{
  ht->card = 0;
  ht->free_elts = 0;
  ht->elts_len = 1;
  ht->table_sz = tabsz;
  ht->table = ml_calloc(tabsz, sizeof(SIZE_T));
  if (ht->table == NULL)
    return -1;
  ht->elts_sz = eltsz;
  ht->elts = ml_malloc(eltsz * sizeof(PREFIX(keyval_t)));
  if (ht->elts == NULL)
    {
      ml_free(ht->table);
      return -1;
    }
  return 0;
}

INLINE HASHTAB * PREFIX(new) (SIZE_T tabsz, SIZE_T eltsz)
{
  HASHTAB *ht = (HASHTAB *) ml_malloc(sizeof(HASHTAB));
  if (ht == NULL)
    return NULL;
  if (PREFIX(init) (ht, tabsz, eltsz) != 0)
    {
      ml_free(ht);
      return NULL;
    }
  return ht;
}

INLINE void PREFIX(dealloc) (HASHTAB *ht)
{
  ml_free(ht->table);
  ml_free(ht->elts);
}

INLINE void PREFIX(free) (HASHTAB *ht)
{
  ml_free(ht->table);
  ml_free(ht->elts);
  ml_free(ht);
}


INLINE void PREFIX(reset) (HASHTAB *ht)
{
  memset(ht->table, 0, ht->table_sz * sizeof(SIZE_T));
  ht->card = 0;
  ht->free_elts = 0;
  ht->elts_len = 1;
}


int PREFIX(_grow_table) (HASHTAB *ht, SIZE_T sz);
int PREFIX(_grow_elts) (HASHTAB *ht, SIZE_T sz);

INLINE int PREFIX(makeroom) (HASHTAB *ht, SIZE_T sz)
{
  if (USE_FACTOR * sz > ht->table_sz)
    {
      if (PREFIX(_grow_table) (ht, sz) != 0)
        return -1;
    }
  /* First entry of ht->elts not used. */
  if (sz + 1 > ht->elts_sz)
    {
      if (PREFIX(_grow_elts) (ht, sz + 1) != 0)
        return -1;
    }
  return 0;
}


/* Return pointer to keyval_t, NULL if key not in table. */
INLINE PREFIX(keyval_t) * PREFIX(lookup) (HASHTAB *ht, KEY_T key, HASH_T hash)
{
  PREFIX(keyval_t) *elts = ht->elts;
  SIZE_T index, i;
  index = hash % ht->table_sz;
  i = ht->table[index];
  while (i != 0 && KEY_CMP(key, elts[i].key) != 0)
    i = ht->elts[i].next;
  return (i == 0) ? NULL : elts + i;
}


/* Call only if key is not in table.  Insert key into table and return
   a pointer to new value variable, NULL if memory allocation
   error. */
INLINE PREFIX(keyval_t) * PREFIX(insert)
     (HASHTAB *ht, KEY_T key, HASH_T hash, VALUE_T value)
{
  PREFIX(keyval_t) *elts, *kvs;
  SIZE_T index, i;
  if (PREFIX(makeroom) (ht, ht->card + 1) != 0)
    return NULL;
  ht->card++;
  elts = ht->elts;
  if (ht->free_elts != 0)
    {
      i = ht->free_elts;
      ht->free_elts = elts[i].next;
    }
  else
    {
      i = ht->elts_len++;
    }
  kvs = elts + i;
  kvs->key = key;
  kvs->hash = hash;
  kvs->value = value;
  index = hash % ht->table_sz;
  kvs->next = ht->table[index];
  ht->table[index] = i;
  return kvs;
}


/* Remove key from hashtable; return pointer to removed keyval_t, or NULL. */
INLINE PREFIX(keyval_t) * PREFIX(remove) (HASHTAB *ht, KEY_T key, HASH_T hash)
{
  PREFIX(keyval_t) * elts = ht->elts;
  SIZE_T i, *pi;
  pi = ht->table + (hash % ht->table_sz);
  i = *pi;
  while (i != 0 && KEY_CMP(key, elts[i].key) != 0)
    {
      pi = & elts[i].next;
      i = *pi;
    }
  if (i == 0)
    return NULL;
  ht->card--;
  *pi = elts[i].next;
  elts[i].next = ht->free_elts;
  ht->free_elts = i;
  return elts + i;
}


/* Return 1 if equal; ignore zero values unless opt_zero != 0. */
int PREFIX(equals)(HASHTAB *ht1, HASHTAB *ht2, int opt_zero);


void PREFIX(print_stat)(HASHTAB *ht);


typedef struct {
  HASHTAB *ht;
  size_t index;
  size_t i;
} PREFIX(iter);

INLINE int PREFIX(good) (PREFIX(iter) *itr)
{
  return (itr->i != 0);
}

INLINE void PREFIX(first) (HASHTAB *ht, PREFIX(iter) *itr)
{
  SIZE_T index;
  itr->ht = ht;
  index = 0;
  while (index < ht->table_sz && ht->table[index] == 0)
    index++;
  if (index == ht->table_sz)
    {
      itr->i = 0;
      return;
    }
  itr->index = index;
  itr->i = ht->table[index];
}

INLINE void PREFIX(next) (PREFIX(iter) *itr)
{
  HASHTAB *ht = itr->ht;
  PREFIX(keyval_t) *elts = ht->elts;
  SIZE_T index;
  if (elts[itr->i].next != 0)
    {
      itr->i = elts[itr->i].next;
      return;
    }
  index = itr->index + 1;
  while (index < ht->table_sz && ht->table[index] == 0)
    index ++;
  if (index == ht->table_sz)
    {
      itr->i = 0;
      return;
    }
  itr->index = index;
  itr->i = ht->table[index];
}

INLINE KEY_T PREFIX(key) (PREFIX(iter) *itr)
{
  return itr->ht->elts[itr->i].key;
}

INLINE VALUE_T PREFIX(value) (PREFIX(iter) *itr)
{
  return itr->ht->elts[itr->i].value;
}

INLINE PREFIX(keyval_t) * PREFIX(keyval) (PREFIX(iter) *itr)
{
  return itr->ht->elts + itr->i;
}

INLINE void PREFIX(dealloc_refs) (HASHTAB *ht)
{
  PREFIX(iter) itr;
  for (PREFIX(first)(ht, &itr); PREFIX(good)(&itr); PREFIX(next)(&itr))
    {
      PREFIX(keyval_t) *kv = PREFIX(keyval)(&itr);
      KEY_DEALLOC(kv->key);
      VALUE_DEALLOC(kv->value);
    }
}

INLINE void PREFIX(dealloc_all) (HASHTAB *ht)
{
  PREFIX(dealloc_refs)(ht);
  PREFIX(dealloc)(ht);
}

INLINE void PREFIX(free_all) (HASHTAB *ht)
{
  PREFIX(dealloc_refs)(ht);
  PREFIX(free)(ht);
}


#ifdef HASHTAB_LINCOMB

#define LC_COPY_KEY 1
#define LC_FREE_KEY 0

#define LC_FREE_ZERO 2
#define LC_KEEP_ZERO 0

INLINE int PREFIX(add_element)
     (HASHTAB *ht, VALUE_T c, KEY_T key, HASH_T hash, int opt)
{
  PREFIX(keyval_t) *kv;
  if (c == 0)
    {
      if (!(opt & LC_COPY_KEY))
        KEY_DEALLOC(key);
      return 0;
    }
  kv = PREFIX(lookup)(ht, key, hash);
  if (kv != NULL)
    {
      if (!(opt & LC_COPY_KEY))
        KEY_DEALLOC(key);
      kv->value += c;
      if (kv->value == 0 && (opt & LC_FREE_ZERO))
        {
          PREFIX(remove)(ht, kv->key, hash);
          KEY_DEALLOC(kv->key);
        }
      return 0;
    }
  if (PREFIX(makeroom)(ht, ht->card + 1) != 0)
    {
      if (!(opt & LC_COPY_KEY))
        KEY_DEALLOC(key);
      return -1;
    }
  if (opt & LC_COPY_KEY)
    {
      key = KEY_COPY(key);
      if (key == NULL)
        return -1;
    }
  kv = PREFIX(insert)(ht, key, hash, c);
  return 0;
}

INLINE int PREFIX(add_multiple)(HASHTAB *dst, VALUE_T c, HASHTAB *src, int opt)
{
  PREFIX(iter) itr;
  for (PREFIX(first)(src,&itr); PREFIX(good)(&itr); PREFIX(next)(&itr))
    {
      PREFIX(keyval_t) *kv = PREFIX(keyval)(&itr);
      if (PREFIX(add_element)(dst, c * kv->value, kv->key, kv->hash, opt) != 0)
        return -1;
    }
  return 0;
}

void PREFIX(print)(HASHTAB *ht, int opt_zero);

#endif
