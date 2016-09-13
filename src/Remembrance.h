
#ifndef _H_REMEMBRANCE
#define _H_REMEMBRANCE

#include "stdlib.h"
#include "stdio.h"

#include "Options.h"

/* Packing these enables sorting on the remembrance structures */
struct remembrance_item {
    unsigned long index;
    unsigned long id;
    unsigned long parent;
};

/* A self-growing vector struct used for
   remembering where the empty separator
   blocks relate to. */
struct empty_cuts {
    unsigned long size;
    unsigned long alloc;
    struct remembrance_item * vector;
};

/* A self-growing vector struct used for
   remembering where the separator blocks
   end up. */
struct remembrance {
    unsigned long free_id;
    unsigned long size;
    unsigned long alloc;
    struct remembrance_item * vector;
    struct empty_cuts empty;
};

/** Maintenance functions and return functions. Note that
    the real work is in the remembrance functions defined
    and called in Permute.c/.h */
struct remembrance* remembrance_init(void);
struct remembrance_item * remembrance_get(struct remembrance * const m);
int remembrance_destroy(struct remembrance * m);
unsigned long remembrance_getFreeID(struct remembrance * const m);
int remembrance_relateHierarchy(struct remembrance * m);
int remembrance_write(struct remembrance * m, FILE * file);
int remembrance_item_compare( const void *one, const void *two );

#endif

