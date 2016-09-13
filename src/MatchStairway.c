#include "Sort.h"
#include "MatchStairway.h"

struct stairrange {
    long InnerProduct;
    long Size;
    long IntersectCount;
	long Child;
    int Visited;
};

struct stairorderedrange {
    long InnerProduct;
    long iStart, iEnd;
};

struct stairwaydata {
	long NrIndices;
	long NrRanges;
    long NeighborhoodSize;
    long *RangeMap;
	long *VertexIndices;
	struct stairorderedrange *Ranges;
};

struct stairrange CreateStairRange(const long InnerProduct,
                                   const long Size) {
    /*
    Create a stair range with the desired parameters.
    */
    struct stairrange Range;
    
    Range.InnerProduct = InnerProduct;
    Range.Size = Size;
    Range.IntersectCount = 0;
    Range.Child = 0;
    Range.Visited = FALSE;
    
    return Range;
}

struct stairorderedrange CreateStairOrderedRange(const long InnerProduct,
                                                 const long Start,
                                                 const long End) {
    /*
    Create an ordered stair range with the desired parameters.
    */
    struct stairorderedrange Range;
    
    Range.InnerProduct = InnerProduct;
    Range.iStart = Start;
    Range.iEnd = End;
    
    return Range;
}

int MatchStairwaySetup(void **ppData,
                       struct biparthypergraph *pHG,
                       const struct opts *pOptions) {
    /*
    This function constructs a collection of ranges of vertex indices for
    a given hypergraph, such that all pairs of vertices that belong to
    the same range, have a guaranteed minimum mutual inner product.
    */
    long NrRanges;
    struct stairrange *Ranges;
    long *NrNetNonzeros;
    long *NetOrdering;
    struct stairrange **Visited;
    long t, tt;
    long Offset;
    struct stairwaydata *pStairway;
    long NeighborhoodSize = 1;
    
    pStairway = (struct stairwaydata *)malloc(sizeof(struct stairwaydata));
    
    if (!pStairway || !pHG) {
        fprintf(stderr, "MatchStairwaySetup(): Not enough memory or invalid hypergraph!\n");
        return FALSE;
    }
    
    /* Initialize the stairway data structure. */
    NeighborhoodSize = MAX(sqrt(pHG->NrVertices), (3*pHG->NrPins)/(2*pHG->NrNets));
    pStairway->NrIndices = pHG->NrVertices;
    pStairway->NrRanges = 0;
    pStairway->NeighborhoodSize = NeighborhoodSize;
    pStairway->RangeMap = (long *)malloc(pHG->NrVertices*sizeof(long));
    pStairway->VertexIndices = (long *)malloc(pStairway->NrIndices*sizeof(long));
    pStairway->Ranges = (struct stairorderedrange *)malloc(pStairway->NrIndices*sizeof(struct stairorderedrange));
    
    if (!pStairway->RangeMap || !pStairway->VertexIndices || !pStairway->Ranges) {
        fprintf(stderr, "MatchStairwaySetup(): Not enough memory!\n");
        return FALSE;
    }
    
    /* Sort adjacency lists. */
    if (!SortAdjacencyLists(pHG)) {
        fprintf(stderr, "MatchStairwaySetup(): Unable to sort adjacency lists!\n");
        return FALSE;
    }
    
    /* Create temporary arrays for the ranges. */
    Ranges = (struct stairrange *)malloc(pHG->NrVertices*sizeof(struct stairrange));
    Visited = (struct stairrange **)malloc(pHG->NrVertices*sizeof(struct stairrange *));
    NrNetNonzeros = (long *)malloc(pHG->NrNets*sizeof(long));
    
    if (!Ranges || !Visited || !NrNetNonzeros) {
        fprintf(stderr, "MatchStairwaySetup(): Not enough memory!\n");
        return FALSE;
    }
    
    /* We start with a single range containing all vertices. */
    NrRanges = 1;
    
    Ranges[0] = CreateStairRange(0, pHG->NrVertices);
    
    for (t = 0; t < pHG->NrVertices; t++) {
        pStairway->RangeMap[t] = 0;
    }
    
    /* Then, we visit all nets in the hypergraph, in order of decreasing number of nonzeros. */
    for (t = 0; t < pHG->NrNets; t++) {
        NrNetNonzeros[t] = pHG->N[t].iEnd - pHG->N[t].iStartP0;
    }
    
    NetOrdering = QSort(NrNetNonzeros, pHG->NrNets);
    
    if (!NetOrdering) {
        fprintf(stderr, "MatchStairwaySetup(): Unable to sort nets by decreasing number of nonzeros!\n");
        return FALSE;
    }
    
    /* Process each net, one-by-one. */
    for (t = 0; t < pHG->NrNets; t++) {
        const struct net Net = pHG->N[NetOrdering[t]];
        long NrVisited = 0;
        
        /* Terminate the algorithm if the nets are too small. */
        if (Net.iEnd - Net.iStartP0 <= NeighborhoodSize) {
            break;
        }
        
        /* Look at all ranges that have an intersection with the current net. */
        for (tt = Net.iStartP0; tt < Net.iEnd; tt++) {
            const long vtx = pHG->NetAdjncy[tt];
            struct stairrange *r = &Ranges[pStairway->RangeMap[vtx]];
            
            /* Do we encounter this range for the first time? */
            if (!r->Visited) {
                Visited[NrVisited++] = r;
                r->Visited = TRUE;
                r->IntersectCount = 0;
            }
            
            /* One more vertex contained both in this range and in the current net. */
            r->IntersectCount++;
        }
        
        /* Determine for each range we encountered whether or not we want to split it. */
        for (tt = 0; tt < NrVisited; tt++) {
            struct stairrange *r = Visited[tt];
            
#ifdef INFO
            /* Perform an extra sanity check. */
            if (r->IntersectCount < 1 || r->IntersectCount > r->Size || !r->Visited) {
                fprintf(stderr, "MatchStairwaySetup(): Invalid range intersection count!\n");
                return FALSE;
            }
#endif
            /* Do not split the range by default. */
            r->Child = 0;
            
            /* Is the entire range contained in this net? */
            if (r->IntersectCount == r->Size) {
                /* Move the entire range up one stair. */
                r->InnerProduct++;
            }
            /* Are we going to split this range? */
            else if (r->IntersectCount >= NeighborhoodSize && r->IntersectCount <= r->Size - NeighborhoodSize
                    && (double)rand()/(double)RAND_MAX <= (double)r->IntersectCount/(double)MIN(r->Size, Net.iEnd - Net.iStartP0)) {
                /* Create a new range. */
                Ranges[NrRanges] = CreateStairRange(r->InnerProduct + 1, r->IntersectCount);
                r->Size -= r->IntersectCount;
                r->Child = NrRanges++;
            }
            
            /* Reset auxiliary variables. */
            r->IntersectCount = 0;
            r->Visited = FALSE;
        }
        
        /* Assign vertices to the new ranges. */
        for (tt = Net.iStartP0; tt < Net.iEnd; tt++) {
            const long vtx = pHG->NetAdjncy[tt];
            const struct stairrange *r = &Ranges[pStairway->RangeMap[vtx]];
            
            if (r->Child > 0) {
                pStairway->RangeMap[vtx] = r->Child;
            }
        }
    }

#ifdef INFO
    /* Have we included all vertices in all ranges? */
    if (TRUE) {
        long Count = 0;
        
        for (t = 0; t < NrRanges; t++) {
            Count += Ranges[t].Size;
        }
        
        if (Count != pHG->NrVertices) {
            fprintf(stderr, "MatchStairwaySetup(): Not all vertices were assigned!\n");
            return FALSE;
        }
    }
    
    /* Are the range sizes corect? */
    if (TRUE) {
        for (t = 0; t < NrRanges; t++) {
            Ranges[t].IntersectCount = 0;
        }
        
        for (t = 0; t < pHG->NrVertices; t++) {
            if (pStairway->RangeMap[t] < 0 || pStairway->RangeMap[t] >= NrRanges) {
                fprintf(stderr, "MatchStairwaySetup(): Invalid range index!\n");
                return FALSE;
            }
            
            Ranges[pStairway->RangeMap[t]].IntersectCount++;
        }
        
        for (t = 0; t < NrRanges; t++) {
            if (Ranges[t].IntersectCount != Ranges[t].Size) {
                fprintf(stderr, "MatchStairwaySetup(): Invalid range size!\n");
                return FALSE;
            }
            
            Ranges[t].IntersectCount = 0;
        }
    }
#endif
    
    /* Now, we generate ranges containing the vertices in the different stairs. */
    Offset = 0;
    pStairway->NrRanges = 0;
    
    for (t = 0; t < NrRanges; t++) {
        pStairway->Ranges[pStairway->NrRanges++] = CreateStairOrderedRange(Ranges[t].InnerProduct, Offset, Offset);
        Offset += Ranges[t].Size;
    }
    
    /* And store the vertex indices via a counting sort. */
    for (t = 0; t < pHG->NrVertices; t++) {
        pStairway->VertexIndices[pStairway->Ranges[pStairway->RangeMap[t]].iEnd++] = t;
    }

#ifdef INFO
    if (TRUE) {
        /* Are the counts correct? */
        for (t = 0; t < NrRanges; t++) {
            if (Ranges[t].Size != pStairway->Ranges[t].iEnd - pStairway->Ranges[t].iStart) {
                fprintf(stderr, "MatchStairwaySetup(): Invalid ordered range size!\n");
                return FALSE;
            }
        }
        
        /* Are the ordered ranges consecutive? */
        for (t = 0; t < NrRanges - 1; t++) {
            if (pStairway->Ranges[t].iEnd != pStairway->Ranges[t + 1].iStart) {
                fprintf(stderr, "MatchStairwaySetup(): Ordered ranges are not consecutive!\n");
                return FALSE;
            }
        }
        
        if (pStairway->Ranges[NrRanges - 1].iEnd != pHG->NrVertices) {
            fprintf(stderr, "MatchStairwaySetup(): Ordered ranges are not consecutive!\n");
            return FALSE;
        }
    }
#endif
    
    /* Free data. */
    free(Ranges);
    free(NrNetNonzeros);
    free(NetOrdering);
    free(Visited);

#ifdef INFO
    printf("Successfully constructed a stairway with %ld ranges, from %ld vertices and %ld nets.\n", pStairway->NrRanges, pHG->NrVertices, pHG->NrNets);
#endif
    
    /* Return pointer to the stairway data. */
    *ppData = pStairway;
    
    return TRUE;
}

int MatchStairwayFree(void *pData) {
    /* Free arrays associated with a stairway. */
    struct stairwaydata *pStairway = (struct stairwaydata *)pData;
    
    if (pStairway)
    {
        pStairway->NrIndices = 0;
        pStairway->NrRanges = 0;
        
        if (pStairway->RangeMap) {
            free(pStairway->RangeMap);
        }
        
        pStairway->RangeMap = NULL;
        
        if (pStairway->VertexIndices) {
            free(pStairway->VertexIndices);
        }
        
        pStairway->VertexIndices = NULL;
        
        if (pStairway->Ranges) {
            free(pStairway->Ranges);
        }
        
        pStairway->Ranges = NULL;
    }
    
    return TRUE;
}

int FindNeighborStairway(long *pNeighbor, double *pNeighborIp,
                         const struct biparthypergraph *pHG, const struct contraction *pC, const struct opts *pOptions,
                         const long v, const int *Matched,
                         void *pData, long *Visited, long *Inprod, double *ScInprod)
{
    /* This function find a neighbor w of a given vertex v without changing the associated hypergraph. */
    /* This is done using a previously constructed stairway. */
    /* This method is NOT COMPATIBLE with net scaling. */
    /* The adjacency lists of vertices HAVE to be sorted by SortAdjacencyLists(). */
    /* See FindMatchInprod. */
    long t, tt;
    /* Initialise vertex weight and degree. */
    const long vtxwgt = pHG->V[v].vtxwgt;
    long w = -1;
    long maxip = -1;    /* maximum inner product found */
    long tiebreakdegree = LONG_MAX;
    long nvisited = 0;
    struct stairwaydata *pStairway = (struct stairwaydata *)pData;
    struct stairorderedrange *r;
    
    /* Initialise to an invalid neighbor. */
    *pNeighbor = -1;
    
    /* First, consider unmatched neighbors in the same range. */
    r = &pStairway->Ranges[pStairway->RangeMap[v]];
    if (r->InnerProduct > 0) {
        for (t = r->iStart; t < r->iEnd && nvisited <= pStairway->NeighborhoodSize; t++) {
            long v2 = pStairway->VertexIndices[t];
            
            /* Remove already matched vertices from ranges. */
            while (Matched[v2] != FALSE && r->iStart < r->iEnd) {
                r->iEnd--;
                v2 = (pStairway->VertexIndices[t] = pStairway->VertexIndices[r->iEnd]);
            }
            
            if (v2 != v &&
                Matched[v2] == FALSE) {
                Inprod[v2] = r->InnerProduct;
                Visited[nvisited++] = v2;
            }
        }
    }
    
    /* Compute the inner products by traversing the adjacent nets. */
    for (t = pHG->V[v].iStart; t < pHG->V[v].iEnd && nvisited <= 4*pStairway->NeighborhoodSize; t++) {
        const long n = pHG->VtxAdjncy[t];
        long Count = pStairway->NeighborhoodSize;
        
        /* Traverse net n */
        for (tt = pHG->N[n].iStartP0; tt < pHG->N[n].iStartP1 && Count > 0; tt++) {
            const long v2 = pHG->NetAdjncy[tt];
            
            if (v2 != v &&
                Matched[v2] == FALSE) {
                /* Register new visits to vertices and store them in Visited */
                if (Inprod[v2] == 0) {
                    Visited[nvisited++] = v2;
                }
                
                /* Update inner products. */
                Inprod[v2]++;
                Count--;
            }
        }
    }
    
    /* A valid neighbor w of v should have the following characteristics:
        - w != v
        - w should be unmatched
        - w should not be free
        - the vertex weight of v and w combined should not exceed the threshold.
    */
    for (t = 0; t < nvisited; t++) {
        const long v2 = Visited[t];
        
        if (!pHG->V[v2].Free &&
            vtxwgt + pHG->V[v2].vtxwgt <= pC->MaxVtxWgt) {
            const long vtxdeg2 = pHG->V[v2].iEnd - pHG->V[v2].iStart;
            const long ip = Inprod[v2];
            
            /* Is this a better candidate? */
            if (ip > maxip ||
                (ip == maxip && vtxdeg2 < tiebreakdegree)) {
                maxip = ip;
                tiebreakdegree = vtxdeg2;
                w = v2;
            }
        }
    }
    
    *pNeighbor = w;
    *pNeighborIp = maxip;
    
    /* Clear weights. */
    for (t = 0; t < nvisited; t++) {
        const long v2 = Visited[t];
        
        Inprod[v2] = 0;
    }
    
    return (maxip > 0 && w >= 0 ? TRUE : FALSE);
}

