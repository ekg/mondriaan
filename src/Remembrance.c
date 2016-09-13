
#include "Remembrance.h"

/** Allocates and returns a remembrance struct initialised to hold
    two elements in the vector, the first of which is set to zero.
    The vector is pre-allocated to be able to hold four elements.
    Caller should set the second element to the appropiate dimension. */
struct remembrance* remembrance_init() {
    struct remembrance *ret = malloc( sizeof( struct remembrance ) );
    if( ret == NULL )
         fprintf(stderr, "remembrance_init(): Error during initialisation!\n");
    else {
        ret->size  = 0;
        ret->alloc = 2;
        ret->vector = malloc( ret->alloc*sizeof( struct remembrance_item ) );
        ret->empty.size = 0;
        ret->empty.alloc = 1;
        ret->empty.vector = malloc( ret->empty.alloc*sizeof( struct remembrance_item ) );
        ret->free_id = 1;
    }
    return ret;
}

/** Used for sorting remembrance vectors on index */
int remembrance_item_compare( const void *one, const void *two ) {
    const unsigned long a = ((struct remembrance_item*)one)->index;
    const unsigned long b = ((struct remembrance_item*)two)->index;
    if( a>b ) return  1;
    if( a<b ) return -1;
    return 0;
}


/** Will free a remembrance struct. */
int remembrance_destroy( struct remembrance * m ) {
    if( m == NULL ) {
        fprintf( stderr, "remembrance_destroy(): Null pointer supplied!\n" );
        return FALSE;
    }
    if( m->vector == NULL || m->empty.vector == NULL ) {
        fprintf( stderr, "remembrance_destroy(): Invalid remembrance structure!\n" );
        return FALSE;
    }

    free( m->vector );
    free( m->empty.vector );
    free( m );
    return TRUE;
}

struct remembrance_item * remembrance_get( struct remembrance * const m ) {
    if( m == NULL ) {
        fprintf( stderr, "remembrance_get(): Null pointer supplied!\n" );
        return NULL;
    }
    if( m->vector == NULL ) {
        fprintf( stderr, "remembrance_get(): Invalid remembrance struct!\n" );
        return NULL;
    }
    if( !remembrance_relateHierarchy( m ) ) {
        fprintf( stderr, "remembrance_get(): Unable to get relative hierarchy!\n" );
        return NULL;
    }
    return m->vector;
}

unsigned long remembrance_getFreeID( struct remembrance * const m ) {
    return m->free_id++;
}

/* This function invalidates the data in the empty cuts struct! (m->empty) */
int remembrance_relateHierarchy(struct remembrance * m) {
    long i;
    const long int size = m->size;
    struct remembrance_item * const vector = m->vector;

    /* extend and abuse empty cuts vector for index translation */
    /* worst case memory usage: O(p), where p is the number of parts requested;
       this is due to the number of IDs given out by remembrance is of O(p) */
    const long int tr_size = remembrance_getFreeID( m );
    unsigned long *translate;
    m->empty.vector = realloc( m->empty.vector, tr_size * sizeof( unsigned long ) );
    translate = ((unsigned long*)(m->empty.vector));
    if( translate == NULL ) {
        fprintf(stderr, "remembrance_relateHierarchy(): Realloc failed!\n");
        return FALSE;
    }
    for (i = 0; i < tr_size; i++ ) translate[ i ] = ULONG_MAX;
    for (i = 0; i < size; i++) {
        if( vector[i].id == ULONG_MAX ) continue;
        translate[ vector[ i ].id ] = i;
    }
    for (i = 0; i < tr_size; i++ )
        if( translate[ i ] == ULONG_MAX ) {
            fprintf( stderr, "remembrance_relateHierarchy(): ID %ld unset!\n", i );
            return FALSE;
        }
    for (i = 0; i < size; i++) {
        vector[ i ].parent = ( vector[i].parent == ULONG_MAX ? ULONG_MAX : translate[ vector[ i ].parent ] );
    }
    return TRUE;
}

int remembrance_write(struct remembrance * m, FILE *File) {
    long i;

    /*
    if (File == NULL || m == NULL) {
        fprintf(stderr, "remembrance_write(): Null arguments!\n");
        return FALSE;
    }
    */

    const struct remembrance_item * const vector = m->vector;
    const long int size = m->size;

    for (i = 0; i < size-1; i++) {
        /* Write only indices and relative hierarchy; not the internal IDs */
        fprintf(File, "%lu\t%lu\n", vector[i].index + 1,
                      (vector[i].parent == ULONG_MAX ? 0 : vector[i].parent + 1)
               );
    }
    fprintf(File, "%lu\n", vector[size-1].index+1);

    return TRUE;
}

