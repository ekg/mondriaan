#ifndef __Permute_h__
#define __Permute_h__

#include "Options.h"
#include "SparseMatrix.h"
#include "Remembrance.h"

struct ordertreerange {
    long PiLo, PiHi;

    long LeftCount, CutCount, RightCount;
    long LeftStart, CutStart, RightStart;
    long LeftChild, CutChild, RightChild;
    int WasTagged;

    /* For use with remembrance,
       to get partial tree hierarchy
       and separator indices. */
    int Meta; /* 0 block, 1 cut, 2 multiple cuts */
    unsigned long id, parent;    
};

struct ordertree {
    int *NetAssign;
    long *NetToRange;
    long *Pi;
    
    struct ordertreerange *Ranges;
    long NrRanges, NrNets;
    
    struct ordertreerange **TagList;
    long NrTags;
};

int CreateOrderTree(struct ordertree *Tree, const long NrNets, const long Lo, const long Hi, long *Net);
int GeneralCreateOrderTree(struct ordertree *Tree, const long NrNets, const long Lo, const long Hi, long * const * const Nets, const long NrNetArrays);
int DestroyOrderTree(struct ordertree *Tree);
int OrderTreeSplit(struct ordertree *Tree, const long Lo, const long Mid, const long Hi, long *Net, const struct opts *pOptions, struct remembrance * const m );
int GeneralOrderTreeSplit(struct ordertree *Tree, const long Lo, const long Mid, const long Hi, long * const * const Nets, const long NrNets, const struct opts *pOptions, struct remembrance * const m );

/* Remembrance data modifier functions */
int remembrance_add( struct remembrance * const m, const unsigned long x, const unsigned long id, const unsigned long parent );
int remembrance_empty( struct remembrance * const m, const unsigned long id, const unsigned long index, const unsigned long parent );
int remembrance_combine(struct remembrance * const m, struct ordertreerange * const r, const long r_size, const char empty_first);
/* Comparator functions */
int remembrance_ordertreerange_compare( const void *one, const void *two );
#endif

