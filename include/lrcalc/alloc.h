#ifndef _ALLOC_H
#define _ALLOC_H

#include <stdlib.h>
#include <stdio.h>


/* PINLINE: No compiled copy.  Can be used in many compilation units, but
            a compiled must be created with CINLINE in one unit.
   CINLINE: Include compiled copy.  Should be used in exactly one unit.
*/
#if __GNUC__ && !__GNUC_STDC_INLINE__
#define PINLINE extern inline
#define CINLINE inline
#else
#define PINLINE inline
#define CINLINE extern inline
#endif
#define INLINE PINLINE


#ifdef DEBUG
#define CLAIM
#endif

/* The macro claim works like the standard C macro assert except that when
   it fails it provokes a segmentation fault, which makes it possible to
   use a gdb.
*/
#ifdef CLAIM
#define claim(ex) ((void) ((ex) || (fprintf(stderr, "%s:%d: Claim `%s' failed\n", __FILE__, __LINE__, #ex), (*(char *)-1)=7)))
#else
#define claim(ex) ((void) 0)
#endif


#ifdef DEBUG
#define DEBUG_MEMORY
#endif

#ifndef DEBUG_MEMORY

#define ml_malloc malloc
#define ml_calloc calloc
#define ml_realloc realloc
#define ml_free free

#define alloc_reset()
#define alloc_getenv()
#define alloc_report()
#define alloc_set_print(prt)
#define alloc_set_trap(index)
#define alloc_set_fail(N, sf)
#define alloc_set_trap_number(N)
#define alloc_test_oom(f,arg)

#else /* DEBUG_MEMORY */

void *ml_malloc(size_t size);
void *ml_calloc(size_t num, size_t size);
void *ml_realloc(void *p, size_t size);
void ml_free(void *);

void alloc_reset();
void alloc_getenv();
void alloc_report();
void alloc_set_print(int prt);
void alloc_set_trap(size_t index);
void alloc_set_fail(size_t N, int segfault);
void alloc_set_trap_number(size_t N);
void alloc_test_oom(void (f)(void *), void *arg);

#endif /* DEBUG_MEMORY */

#endif /* _ALLOC_H */
