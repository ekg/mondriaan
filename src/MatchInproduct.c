#include "MatchInproduct.h"

int MatchInproductSetup(void **ppData,
                       struct biparthypergraph *pHG,
                       const struct opts *pOptions) {
    /* We do not require any data for inner product matching. */
    *ppData = NULL;
    
    return TRUE;
}

int MatchInproductFree(void *pData) {
    /* No data should have been allocated. */
    if (pData != NULL)
    {
        fprintf(stderr, "MatchInproductFree(): No data should have been allocated!\n");
        return FALSE;
    }
    
    return TRUE;
}

int FindNeighborInproduct(long *pNeighbor, double *pNeighborIp,
                         const struct biparthypergraph *pHG, const struct contraction *pC, const struct opts *pOptions,
                         const long v, const int *Matched,
                         void *pData, long *Visited, long *Inprod, double *ScInprod)
{
    /* This function find a neighbor w of a given vertex v without changing the associated hypergraph. */
    /* See FindMatchInprod. */
    long t, tt;
    /* Initialise vertex weight and degree. */
    const long vtxwgt = pHG->V[v].vtxwgt;
    const long vtxdeg = pHG->V[v].iEnd - pHG->V[v].iStart;
    long w = -1;
    double maxip = -1.0;    /* maximum inner product found */
    long tiebreakdegree = LONG_MAX;
    long nvisited = 0;
    
    /* Initialise to an invalid neighbor. */
    *pNeighbor = -1;
    
    /* Compute the inner products by traversing the adjacent nets. */
    for (t = pHG->V[v].iStart; t < pHG->V[v].iEnd; t++) {
        const long n = pHG->VtxAdjncy[t];
        const long netsize = pHG->N[n].iEnd - pHG->N[n].iStartP0;
        /* Compute scale factor of net n. */
        const double scrowwgt = (pOptions->Coarsening_NetScaling == NetSclLinear && netsize > 0 ? 1.0/(double)netsize : 1.0);

        /* Traverse net n */
        for (tt = pHG->N[n].iStartP0; tt < pHG->N[n].iStartP1; tt++) {
            const long v2 = pHG->NetAdjncy[tt];
            
            /* Register new visits to vertices and store them in Visited */
            if (Inprod[v2] == 0) {
                Visited[nvisited++] = v2;
            }
            
            /* Update inner products. */
            Inprod[v2]++;
            ScInprod[v2] += scrowwgt;
        }
    }
    
    /* A valid neighbor w of v should have the following characteristics:
        - w != v
        - w should be unmatched
        - w should not be free
        - the vertex weight of v and w combined should not exceed the threshold.
       It is furthermore desirable (in order of appearance)
        - rescaled ip(v, w) is maximal
        - nz(w) is minimal
    */
    for (t = 0; t < nvisited; t++) {
        const long v2 = Visited[t];
        
        if (v2 != v &&
            Matched[v2] == FALSE &&
            !pHG->V[v2].Free &&
            vtxwgt + pHG->V[v2].vtxwgt <= pC->MaxVtxWgt) {
            /* Ok, so v2 is a valid candidate. */
            const long vtxdeg2 = pHG->V[v2].iEnd - pHG->V[v2].iStart;
            const long mindeg = MIN(vtxdeg, vtxdeg2);
            const long maxdeg = MAX(vtxdeg, vtxdeg2);
            double sccolwgt = 1.0;
            
            switch (pOptions->Coarsening_InprodScaling) {
                case NoIpScaling:
                    sccolwgt = 1.0;
                    break;
                case IpSclCos:
                    sccolwgt = 1.0/sqrt((double)mindeg*maxdeg);
                    break;
                case IpSclMin:
                    sccolwgt = 1.0/(double)mindeg;
                    break;
                case IpSclMax:
                    sccolwgt = 1.0/(double)maxdeg;
                    break;
                case IpSclJaccard: /* the Jaccard metric is the ratio between the
                                 number of nonzeros in the intersection of the 
                                 two columns and the number of nonzeros
                                 in the union */
                    sccolwgt = 1.0/(double)(mindeg + maxdeg - Inprod[v2]);
                    break;
            }
            
            /* Note that v2 can only occur once in the Visited array. */
            ScInprod[v2] *= sccolwgt;
            
            /* Is this a better candidate? */
            if (ScInprod[v2] > maxip ||
                (fabs(ScInprod[v2] - maxip) < 1.0e-12 && vtxdeg2 < tiebreakdegree)) {
                maxip = ScInprod[v2];
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
        ScInprod[v2] = 0.0;
    }
    
    return (maxip > 0.0 && w >= 0 ? TRUE : FALSE);
}

