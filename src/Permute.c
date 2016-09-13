#include <stdlib.h>
#include <stdio.h>

#include "Permute.h"

char mondriaan_warn_empty_blocks = 0;

int CreateOrderTree(struct ordertree *Tree, const long NrNets, const long Lo, const long Hi, long *Net) {
    long **Nets = &Net;
    return GeneralCreateOrderTree(Tree, NrNets, Lo, Hi, Nets, 1l);
}

int GeneralCreateOrderTree(struct ordertree *Tree, const long NrNets, const long Lo, const long Hi, long * const * const Nets, const long NrNetArrays) {
    long t, n, NrEmpty;
    
    if (Tree == NULL) {
        fprintf(stderr, "CreateOrderTree(): Null argument!\n");
        return FALSE;
    }
    
    Tree->NetAssign = (int *)malloc(NrNets*sizeof(int));
    Tree->NetToRange = (long *)malloc(NrNets*sizeof(long));
    Tree->Pi = (long *)malloc(NrNets*sizeof(long));
    Tree->Ranges = (struct ordertreerange *)malloc(NrNets*sizeof(struct ordertreerange));
    Tree->TagList = (struct ordertreerange **)malloc(NrNets*sizeof(struct ordertreerange *));
    
    if (Tree->NetAssign == NULL || Tree->NetToRange == NULL || Tree->Pi == NULL || Tree->Ranges == NULL || Tree->TagList == NULL) {
        fprintf(stderr, "CreateOrderTree(): Not enough memory!\n");
        return FALSE;
    }
    
    Tree->NrRanges = 1;
    
    Tree->Ranges[0].PiLo = 0;
    Tree->Ranges[0].PiHi = NrNets;

    /* Remembrance specific. */
    Tree->Ranges[0].Meta   = 0;
    Tree->Ranges[0].parent = ULONG_MAX;
    Tree->Ranges[0].id     = 0;
    
    Tree->NrNets = NrNets;
    Tree->NrTags = 0;
    
    for (t = 0; t < NrNets; ++t) {
        Tree->NetToRange[t] = 0;
        Tree->Pi[t] = t;
    }
    
    /* Find empty nets. */
    for (t = 0; t < NrNets; ++t) Tree->NetAssign[t] = 0;
    for (t = Lo; t < Hi; ++t)
        for (n=0; n<NrNetArrays; n++)
            Tree->NetAssign[Nets[n][t]] = 1;
    
    /* Permute them to the back. */
    NrEmpty = 0;
    
    for (t = 0; t < NrNets; ++t) {
        if (Tree->NetAssign[t] == 0) Tree->Pi[NrNets - (++NrEmpty)] = t;
    }
    
    Tree->Ranges[0].PiHi = NrNets - NrEmpty;

    /* Permute nonempty nets to the front. */
    NrEmpty = 0;
    
    for (t = 0; t < NrNets; ++t) {
        if (Tree->NetAssign[t] == 1) Tree->Pi[NrEmpty++] = t;
    }
    
    return TRUE;
}

int DestroyOrderTree(struct ordertree *Tree) {
    if (Tree == NULL) {
        fprintf(stderr, "DestroyOrderTree(): Null argument!\n");
        return FALSE;
    }
   
    if (Tree->NetAssign != NULL) free(Tree->NetAssign);
    if (Tree->NetToRange != NULL) free(Tree->NetToRange);
    if (Tree->Pi != NULL) free(Tree->Pi);
    if (Tree->Ranges != NULL) free(Tree->Ranges);
    if (Tree->TagList != NULL) free(Tree->TagList);
    
    memset(Tree, 0, sizeof(struct ordertree));
    
    return TRUE;
}

/* m is the remembrance struct where to store the added data,
   x is the index number to remember. */
int remembrance_add( struct remembrance * const m, const unsigned long x, const unsigned long id, const unsigned long parent ) {

    if( m == NULL ) {
        fprintf( stderr, "remembrance_add(): Null pointer supplied!\n" );
        return FALSE;
    }

    /* Adaptively grow vector. */
    if( m->size == m->alloc ) {
        m->alloc *= 2;
        m->vector = realloc( m->vector, m->alloc*sizeof( struct remembrance_item ) );
        if( m->vector == NULL ) {
            fprintf(stderr, "remembrance_add(): Error during realloc of vector\n");
            return FALSE;
        }
        if( m->size >= ULONG_MAX ) {
            fprintf(stderr, "remembrance_add(): Remembrance vector is too long\n");
            return FALSE;
        }
    }

    /* Add element and increment size field. */
    m->vector[  m->size  ].id     = id;
    m->vector[  m->size  ].parent = parent;
    m->vector[ m->size++ ].index  = x;
    return TRUE;
}

/* Remembers an empty cut, which do not get stored in the OrderTree. */
int remembrance_empty( struct remembrance * const m, const unsigned long id, const unsigned long index, const unsigned long parent ) {
    struct empty_cuts * const e = &(m->empty);

    if( m == NULL ) {
        fprintf( stderr, "remembrance_empty(): Null pointer supplied!\n" );
        return FALSE;
    }

    /* Adaptively grow vector. */
    if( e->size == e->alloc ) {
        e->alloc *= 2;
        e->vector = realloc( e->vector, e->alloc*sizeof( struct remembrance_item ) );
        if( e->vector == NULL ) {
            fprintf(stderr, "remembrance_empty(): Error during realloc of vector\n");
            return FALSE;
        }
        if( e->size >= ULONG_MAX ) {
            fprintf(stderr, "remembrance_empty(): Remembrance vector is too long\n");
            return FALSE;
        }
    }

    /* Add empty cut and increment size field. */
    e->vector[  e->size  ].id     = id;
    e->vector[  e->size  ].index  = index;
    e->vector[ e->size++ ].parent = parent;
    return TRUE;
}

/* Various functions for ordering structures used by remembrance, and hierarchy-extraction specifically.
   NOTE: SORTING THE ORDERTREERANGE INVALIDATES THE ORDERTREE DATASTRUCTURE. Use before DestroyOrderTree. */
int remembrance_ordertreerange_compare( const void *one, const void *two ) {
    const long a = ((struct ordertreerange*)one)->PiLo;
    const long b = ((struct ordertreerange*)two)->PiLo;
    if( a>b ) return  1;
    if( a<b ) return -1;
    return 0;
}
/* Combines information on remembrance and the ordertree to give
   the hierarchy together with the separator indices.
   @param empty_first 0 if ranges should have precedence over empty cuts (BBD), nonzero if otherwise. */
int remembrance_combine(struct remembrance * const m, struct ordertreerange * const r, const long r_size, const char empty_first) {
    unsigned long range_id, e_id;

    /* First sort all ranges and empty cuts on index */
    qsort( r, r_size, sizeof( struct ordertreerange ), &remembrance_ordertreerange_compare );
    qsort( m->empty.vector, m->empty.size, sizeof( struct remembrance_item ), &remembrance_item_compare );

    /* Set current indices */
    range_id = e_id = 0;

    /* Go through all elements in ordertree ranges */
    for( ; range_id < (unsigned long)r_size; range_id++ ) {
        unsigned long index;
        
        /* If this is a range inside a cut, skip */
        if( r[ range_id ].Meta > 1 ) continue;
        /* For brevity, this is the current index at a block boundary */
        index = r[ range_id ].PiLo;
        /* Separated nets should be placed first */
        /* Look if there are indices in cut nets equal to the current one, and add them */
	/* Also check if no free net precedes this index */
        while( e_id < m->empty.size && m->empty.vector[ e_id ].index <= index ) {
            remembrance_add( m, index, m->empty.vector[ e_id ].id, m->empty.vector[ e_id ].parent );
            e_id++;
        }
        /* Add the current range */
        remembrance_add( m, index, r[ range_id ].id, r[ range_id ].parent );
    }
    /* Add any remaining empty cuts; note all empty cuts should be added, in contrast to the tree ranges */
    for( ; e_id < m->empty.size; e_id++ )
       remembrance_add( m, m->empty.vector[ e_id ].index, m->empty.vector[ e_id ].id, m->empty.vector[ e_id ].parent );

    return TRUE;
}

/** Short-hand function for OrderTreeSplit with NrNets=1 (usual case). See directly below. */
int OrderTreeSplit(struct ordertree *Tree, const long Lo, const long Mid, const long Hi, long *Net, const struct opts *pOptions, struct remembrance * const m ) {
    long **Nets = &Net;
    return GeneralOrderTreeSplit(Tree, Lo, Mid, Hi, Nets, 1, pOptions, m);
}

/**
 *  Splits an order tree according to the nets given here. This is based on partitioning;
 *  partition one consists of elements Lo through Mid (exclusive), and partition two of
 *  elements Mid through Hi (exclusive). Each element is included in at least one 
 *  Connected net. A net is considered in category + if each element it contains is in
 *  partition one. A net is considered in category - if each element it contains is in
 *  partition two. Otherwise, a net is considered in category C (cut). Usually a single
 *  element maps to precisely one net (many-to-one). In the symmetric fine-grain case,
 *  a single element maps to precisely two nets. Generalising a bit, this function
 *  assumes that each element maps to <em>exactly</em> NrNets nets (of the same
 *  hypergraph).
 *
 *  @param Tree     The tree structure which keeps track of recursive bipartitioning net categories and their resulting permutation
 *  @param Lo       Lower bound on the partitioned elements, start of partition one
 *  @param Mid      Start of partition two
 *  @param Hi       Upper bound (exclusive) of the partitioned elements
 *  @param Nets     Nets corresponding to each element; element j is contained in Nets[i][j] for all 0<=i<NrNets
 *  @param NrNets   The number of nets each element is contained in (assumed to be >=1)
 *  @param pOptions Current Mondriaan options
 *  @param m        Remembrance structure which keeps track of net hierarchy and boundaries
 */
int GeneralOrderTreeSplit(struct ordertree *Tree, const long Lo, const long Mid, const long Hi, long * const * const Nets, const long NrNets, const struct opts *pOptions, struct remembrance * const m ) {
    long t, net;
    
    if (Tree == NULL || Nets == NULL || pOptions == NULL || NrNets<=0 || Nets[0] == NULL) {
        fprintf(stderr, "OrderTreeSplit(): Null arguments!\n");
        return FALSE;
    }
    
    /* printf("Start with %ld/%ld ranges, split %ld - %ld - %ld.\n", Tree->NrRanges, Tree->NrNets, Lo, Mid, Hi); */
   
    /* Clear all net assignments and mark all ranges as untagged. */
    Tree->NrTags = 0;
    
    for (t = Lo; t < Hi; ++t) {
        for (net = 0; net < NrNets; net++) {
            const long n = Nets[net][t];
            struct ordertreerange *r = &Tree->Ranges[Tree->NetToRange[n]];
        
            Tree->NetAssign[n] = 0;
        
            r->LeftCount = r->CutCount = r->RightCount = 0;
            r->LeftChild = r->CutChild = r->RightChild = -1;
            r->WasTagged = FALSE;
        }
    }
    
    /* Assign all left nets and permute them to the beginning. */
    for (t = Lo; t < Mid; ++t) {
        for (net = 0; net < NrNets; net++) {
            const long n = Nets[net][t];
            struct ordertreerange *r = &Tree->Ranges[Tree->NetToRange[n]];
        
            if (Tree->NetAssign[n] == 0) {
                Tree->NetAssign[n] = 1;
                r->LeftCount++;
            
                if (!r->WasTagged) {
                    r->WasTagged = TRUE;
                    Tree->TagList[Tree->NrTags++] = r;
                }
            }
        }
    }
    
    /* Assign all right nets and permute them to the end. */
    for (t = Mid; t < Hi; ++t) {
        for (net = 0; net < NrNets; net++) {
            const long n = Nets[net][t];
            struct ordertreerange *r = &Tree->Ranges[Tree->NetToRange[n]];
            
            if (Tree->NetAssign[n] == 0) {
                Tree->NetAssign[n] = 2;
                r->RightCount++;
                
                if (!r->WasTagged) {
                    r->WasTagged = TRUE;
                    Tree->TagList[Tree->NrTags++] = r;
                }
            }
            else if (Tree->NetAssign[n] == 1) {
                /* This net is cut: it should be permuted to the middle. If this is the case, then this range has already been tagged! */
                Tree->NetAssign[n] = 3;
                r->LeftCount--;
                r->CutCount++;
            }
        }
    }
    
    /* Update net-to-range mapping. */
    /* Check ranges and set starting offsets. */
    for (t = 0; t < Tree->NrTags; ++t) {
        struct ordertreerange *r = Tree->TagList[t];
        
        if (!r->WasTagged) {
            fprintf(stderr, "OrderTreeSplit(): Sanity check I failed!\n");
            return FALSE;
        }
        
        if (r->PiHi != r->PiLo + r->LeftCount + r->CutCount + r->RightCount) {
            fprintf(stderr, "OrderTreeSplit(): Sanity check II failed!\n");
            fprintf(stderr, "Tag %ld: Lo=%ld, Left=%ld, Cut=%ld, Right=%ld, Hi=%ld\n", t, r->PiLo, r->LeftCount, r->CutCount, r->RightCount, r->PiHi);
            return FALSE;
        }
       
        switch (pOptions->OrderPermutation) {
        case OrderPrefix:
            r->CutStart = r->PiLo;
            r->LeftStart = r->PiLo + r->CutCount;
            r->RightStart = r->PiLo + r->CutCount + r->LeftCount;
            break;
        case OrderPostfix:
            r->LeftStart = r->PiLo;
            r->RightStart = r->PiLo + r->LeftCount;
            r->CutStart = r->PiLo + r->LeftCount + r->RightCount;
            break;
        case OrderInfix:
        default:
            r->LeftStart = r->PiLo;
            r->CutStart = r->PiLo + r->LeftCount;
            r->RightStart = r->PiLo + r->LeftCount + r->CutCount;
            break;
        }
    }
    
    /* Now we determine which ranges need to be split and where their children will spawn. */
    /* Because the original range will be overwritten, the range child index may be equal to the range index. */
    /* Ranges will be overwritten with priority Left, Cut, Right (so if the left range is nonempty. this will overwrite the current range, otherwise the cut range will overwrite the current range and if both the left and cut ranges are empty, the right range will overwrite the current range. */
    for (t = Lo; t < Hi; ++t) {
        for (net = 0; net < NrNets; net++) {
            const long n = Nets[net][t], i = Tree->NetToRange[n];
            struct ordertreerange *r = &Tree->Ranges[i];
        
            /* Only recalculate the appropriate range *once* per net. */
            if (Tree->NetAssign[n] == 4) continue;
        
            switch (Tree->NetAssign[n]) {
            case 1:
                if (r->LeftChild >= 0) Tree->NetToRange[n] = r->LeftChild;
                else r->LeftChild = i;
                
                Tree->Pi[r->LeftStart++] = n;
                break;
            case 2:
                if (r->RightChild >= 0) Tree->NetToRange[n] = r->RightChild;
                else if (r->LeftCount == 0 && r->CutCount == 0) r->RightChild = i;
                else {
                    Tree->NetToRange[n] = r->RightChild = Tree->NrRanges++;
                    Tree->Ranges[r->RightChild].Meta=0;
                }
                
                Tree->Pi[r->RightStart++] = n;
                break;
            case 3:
                if (r->CutChild >= 0) Tree->NetToRange[n] = r->CutChild;
                else if (r->LeftCount == 0) r->CutChild = i;
                else {
                    Tree->NetToRange[n] = r->CutChild = Tree->NrRanges++;
                    Tree->Ranges[r->CutChild].Meta=0;
                }
            
                Tree->Pi[r->CutStart++] = n;
                break;
            default:
                fprintf(stderr, "OrderTreeSplit(): Unassigned net!\n");
                return FALSE;
                break;
            }
        
            Tree->NetAssign[n] = 4;
        }
    }
    
    /* Now actually split the ranges. */
    for (t = 0; t < Tree->NrTags; ++t) {
        /* This will be the range to split. */
        const struct ordertreerange r = *Tree->TagList[t];
        struct ordertreerange *s;
        
        /* printf("Children %ld (count %ld) - %ld (count %ld) - %ld (count %ld), permutation range %ld - %ld, tagged %d.\n", r.LeftChild, r.LeftCount, r.CutChild, r.CutCount, r.RightChild, r.RightCount, r.PiLo, r.PiHi, r.WasTagged); */
        
        if (r.LeftChild >= 0) {
            s = &Tree->Ranges[r.LeftChild];

            /* Remembrance book-keeping */
            if( r.Meta == 0 ) {
                s->Meta   = 0;
                s->id     = remembrance_getFreeID( m );
                s->parent = r.id;
            } else {
                if( r.Meta == 1 && pOptions->OrderPermutation != OrderPrefix )
                      s->Meta = 1; /* Otherwise a meta-1 cut block with this ID is lost, BBD and SBD. */
                else
                      s->Meta = 2;
                s->id     = r.id;
                s->parent = r.parent;
            }
            
            switch (pOptions->OrderPermutation) {
            case OrderPrefix:
                s->PiLo = r.PiLo + r.CutCount;
                s->PiHi = r.PiLo + r.CutCount + r.LeftCount;
                if( r.Meta == 1 ) {
                    if( r.CutCount == 0 )
                        s->Meta = 1; /* Otherwise a meta-1 cut block with this ID is lost. */
                    else
                        s->Meta = 2;
                }
                break;
            case OrderPostfix:
            case OrderInfix:
            default:
                s->PiLo = r.PiLo;
                s->PiHi = r.PiLo + r.LeftCount;
                break;
            }
        } else {
            if( r.Meta == 0 ) {
                if( !mondriaan_warn_empty_blocks ) {
                    mondriaan_warn_empty_blocks = 1;
                    fprintf( stderr, "Warning: empty non-separator blocks encountered, rowblocks/colblocks array may be incomplete!\n" );
                    fprintf( stderr, "         suggest lowering the number of parts; or otherwise, lowering epsilon.\n" );
                }
                /* Remember empty left block */
                /* Note: only works when encountered for the first time; if this empty range
                         would have been cut again, there is no way for the remembrance strategy
                         to detect that */
                switch (pOptions->OrderPermutation) {
                case OrderPrefix:
                   remembrance_empty( m, remembrance_getFreeID( m ), r.PiLo + r.CutCount ,r.id );
                   break;
                case OrderInfix:
                case OrderPostfix:
                default:
                    remembrance_empty( m, remembrance_getFreeID( m ), r.PiLo + r.LeftCount ,r.id );
                    break;
                }
            }
	}
        
        if (r.CutChild >= 0) {
            s = &Tree->Ranges[r.CutChild];
            
            /* Remembrance book-keeping */
            s->Meta = ( r.Meta==0 ? 1 : 2 );
            /* Copy parent ID and parent. */
            s->id     = r.id;
            s->parent = r.parent;

            switch (pOptions->OrderPermutation) {
            case OrderPrefix:
                s->PiLo = r.PiLo;
                s->PiHi = r.PiLo + r.CutCount;
                if( r.Meta == 1 )
                    s->Meta = 1; /* Otherwise a meta-1 cut block with this ID is lost. */
                break;
            case OrderPostfix:
                s->PiLo = r.PiLo + r.LeftCount + r.RightCount;
                s->PiHi = r.PiLo + r.LeftCount + r.RightCount + r.CutCount; 
                if( r.Meta == 1 && r.LeftCount == 0 && r.RightCount == 0 )
                    s->Meta = 1; /* Otherwise a meta-1 cut block with this ID is lost. */
                break;
            case OrderInfix:
            default:
                s->PiLo = r.PiLo + r.LeftCount;
                s->PiHi = r.PiLo + r.LeftCount + r.CutCount;
                if( r.Meta == 1 && r.LeftCount == 0 ) s->Meta = 1; /* Otherwise a meta-1 cut block with this ID is lost. */
                break;
            }
        }

        if( r.CutCount == 0 && r.Meta == 0 ) {
            /* Remembrance: empty cut! */
            switch( pOptions->OrderPermutation) {
            case OrderPostfix:
                remembrance_empty( m, r.id, r.PiHi ,r.parent );
                break;
            case OrderInfix:
                remembrance_empty( m, r.id, r.PiLo + r.LeftCount, r.parent );
                break;
            case OrderPrefix:
                remembrance_empty( m, r.id, r.PiLo ,r.parent );
                break;
            default:
                fprintf( stderr, "Undefined permutation order\n" );
                return FALSE;
            }
        }
        
        if (r.RightChild >= 0) {
            s = &Tree->Ranges[r.RightChild];

            /* Remembrance book-keeping */
            if( r.Meta == 0 ) {
                s->id     = remembrance_getFreeID( m );
                s->parent = r.id;
            } else {
                s->Meta   = 2;
                s->id     = r.id;
                s->parent = r.parent;
            }

            switch (pOptions->OrderPermutation) {
            case OrderPostfix:
                s->PiLo = r.PiLo + r.LeftCount;
                s->PiHi = r.PiLo + r.LeftCount + r.RightCount;
                if( r.Meta == 1 && r.LeftCount == 0 )
                    s->Meta = 1; /* Otherwise a meta-1 cut block with this ID is lost. */
                break;
            case OrderPrefix:
                s->PiLo = r.PiLo + r.LeftCount + r.CutCount;
                s->PiHi = r.PiLo + r.LeftCount + r.CutCount + r.RightCount;
                if( r.Meta == 1 && r.LeftCount + r.CutCount == 0 )
                    s->Meta = 1; /* Otherwise a meta-1 cut block with this ID is lost. */
                break;
            case OrderInfix:
            default:
                s->PiLo = r.PiLo + r.LeftCount + r.CutCount;
                s->PiHi = r.PiLo + r.LeftCount + r.CutCount + r.RightCount;
                if( r.Meta == 1 && r.LeftCount + r.CutCount == 0 )
                    s->Meta = 1; /* Otherwise a meta-1 cut block with this ID is lost. */
            }
        } else {
            if( r.Meta == 0 ) {
                if( !mondriaan_warn_empty_blocks ) {
                    mondriaan_warn_empty_blocks = 1;
                    fprintf( stderr, "Warning: empty non-separator blocks encountered, rowblocks/colblocks array may be incomplete!\n" );
                    fprintf( stderr, "         suggest lowering the number of parts; or otherwise, lowering epsilon.\n" );
                }
                /* remember empty right block */
                /* Note: only works when encountered for the first time; if this empty range
                         would have been cut again, there is no way for the remembrance strategy
                         to detect that */
                switch (pOptions->OrderPermutation) {
                case OrderPostfix:
                    remembrance_empty( m, remembrance_getFreeID( m ), r.PiLo + r.LeftCount ,r.id );
                    break;
                case OrderInfix:
                case OrderPrefix:
                default:
                    remembrance_empty( m, remembrance_getFreeID( m ), r.PiLo + r.LeftCount + r.CutCount, r.id );
                    break;
                }
            }
	}
    }
    
    return TRUE;
}

