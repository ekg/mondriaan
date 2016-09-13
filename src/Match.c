#include "Match.h"

int MoveVtxInNetAdjncy(struct biparthypergraph *pHG, const long v) {

    /* This function moves the vertex v to the back of the adjacency list
       in all its nets, so that it will not be eligible for matching any more.
       This saves much time during matching.

       istartP1 is abused to mark the end of the eligible candidates;
       these are stored in positions iStartP0..iStartP1-1 of each net.

       On input, v must be stored somewhere in the range iStartP0..iStartP1-1
       for all its nets. On output, iStartP1 has been decreased by 1
       for all the nets of v, and the index of v equals the new iStartP1.
        
       The function should only be used during coarsening,
       so before any partitioning.
    */
  
    long t, n, v2, iv, iv2;
    
    if (!pHG) {
        fprintf(stderr, "MoveVtxInNetAdjncy(): Null argument!\n");
        return FALSE;
    }
 
    for (t = pHG->V[v].iStart; t < pHG->V[v].iEnd; t++) {
        n = pHG->VtxAdjncy[t];
 
        /* Find vertex v in adjacency list of net n */
        iv = pHG->N[n].iStartP0; /* index of v */
        
        while (pHG->NetAdjncy[iv] != v && iv < pHG->N[n].iStartP1)
            iv++;
        
        if (iv >= pHG->N[n].iStartP1) {
            fprintf(stderr, "MoveVtxInNetAdjncy(): vertex v (%ld) not found in range!\n", v);
            /*fprintf(stderr, "Vtx2Index: %ld, Net2Index: %ld\n", pHG->Vtx2MatIndex[v], pHG->Net2MatIndex[n] );*/
            return FALSE;
        }
  
        /* Swap vertex:  v <-> Net[n].iStartP1  */
        iv2 = pHG->N[n].iStartP1 - 1;
        v2 = pHG->NetAdjncy[iv2];
        pHG->NetAdjncy[iv2] = v;
        pHG->NetAdjncy[iv] = v2;
        pHG->N[n].iStartP1--;
    }

    return TRUE;
} /* end MoveVtxInNetAdjncy */


int FindMatchArbitrary(struct biparthypergraph *pHG, struct contraction *pC, const long v, int *Matched) {
    /* This function matches the vertex v repeatedly with an arbitrary
       adjacent vertex until the weight exceeds the maximum allowed weight
       MaxVtxWgt, or the number of vertices exceeds the maximum allowed
       in a group of matched vertices, MaxNrVertices. 
       Here, the arbitrary vertex is the first one found.

       Two vertices are adjacent if they occur in a common net. */
      
    long k, t, tt, v2, n, vtxwgt;
    
    if (!pHG || !pC || !Matched) {
        fprintf(stderr, "FindMatchArbitrary(): Null arguments!\n");
        return FALSE;
    }
  
    /* Form a new group, i.e. a new match */
    k = pC->NrMatches; /* Start[0],..., Start[k] have obtained
                           their final values. */
    pC->Start[k+1] = pC->Start[k]; /* Group k will be stored in
                                       Start[k]..Start[k+1]-1 */

    pC->NrMatches++;  /* include the group of v */

    /* Add v to the group, at the first available position */
    pC->Match[pC->Start[k+1]] = v;
    pC->Start[k+1]++;
    vtxwgt = pHG->V[v].vtxwgt;
    Matched[v] = TRUE; /* if no matching vertex will be found, 
                           the vertex is still considered matched
                           and the group has size 1 */
    if (!MoveVtxInNetAdjncy(pHG, v)) return FALSE; /* move the matched vertex to the 
                                                      end (i.e. >= iStartP1) in its
                                                      net adjacency lists */
  
    /* Search in the adjacent nets for an unmatched vertex: */  
    for (t = pHG->V[v].iStart; t < pHG->V[v].iEnd; t++) {
        n = pHG->VtxAdjncy[t]; 
        
        for (tt = pHG->N[n].iStartP0; tt < pHG->N[n].iStartP1; tt++) {
            v2 = pHG->NetAdjncy[tt]; 
    
            /* Add v2 to the current match if
               - v2 is unmatched (guaranteed by construction)
               - total weight of match including v2 <= MaxVtxWgt
               - number of vertices in match excluding v2 < MaxNrVertices */
  
            if (vtxwgt + pHG->V[v2].vtxwgt <= pC->MaxVtxWgt) {
                if (pC->Start[k+1] - pC->Start[k] < pC->MaxNrVertices) {
                    pC->Match[pC->Start[k+1]] = v2;
                    pC->Start[k+1]++;
                    Matched[v2] = TRUE;
                    if (!MoveVtxInNetAdjncy(pHG, v2)) return FALSE;
                    vtxwgt += pHG->V[v2].vtxwgt;
                } else 
                    return TRUE;
            }
        }
    }

    return TRUE;
} /* end FindMatchArbitrary */


int FindMatchInprod(struct biparthypergraph *pHG, struct contraction *pC,
                     const long v, int *Matched,
                     long *Visited, long *Inprod, double *ScInprod, const struct opts *pOptions) {
  
    /* This function matches the vertex v with adjacent
       vertices based on the maximum inner product, optionally scaled,
       until the weight exceeds the maximum allowed weight MaxVtxWgt,
       or the number of vertices exceeds the maximum allowed
       in a group of matched vertices, MaxNrVertices.

       Matched is an array that marks vertices as matched.
       Visited is a working array used to store the indices of
       vertices that have been visited; these are adjacent to v and
       are guaranteed to be still unmatched.

       Inprod and ScInprod are working arrays used to store
       inner products. Initially they must be zero.
       Visited is also a working array. It need not be initialised.
       The size of these working arrays equals the number of vertices.
       On exit they will be restored to zero.
       
       The arrays Visited, Inprod, ScInprod are allocated outside this
       function to save the overhead of memory allocation for repeated calls.
       No memory is allocated inside this function.
    */
  
    long k, t, tt, v2, tmax, n, netsize, vtxwgt, vtxdeg = 0, vtxdeg2 = 0,
         nvisited, mindeg, maxdeg;
    double scrowwgt = 1.0; /* scaled row weight based on net size
                         (in this terminology, rows are identified with nets
                          and columns with vertices, irrespective of how
                          the hypergraph was constructed) */
    double sccolwgt = 1.0; /* scaled column weight based on vertex degree */
    double maxip;    /* maximum inner product found */
    
    if (!pHG || !pC || !Matched || !Visited || !Inprod || !ScInprod || !pOptions) {
        fprintf(stderr, "FindMatchInprod(): Null arguments!\n");
        return FALSE;
    }
    
    /* Form a new group, i.e. a new match */
    k = pC->NrMatches; /* Start[0],..., Start[k] have obtained
                           their final values. */
    pC->Start[k+1] = pC->Start[k]; /* Group k will be stored in
                                       Start[k]..Start[k+1]-1 */

    pC->NrMatches++;  /* include the group of v */

    /* Add v to the group, at the first available position */
    pC->Match[pC->Start[k+1]] = v;
    pC->Start[k+1]++;
    vtxwgt = pHG->V[v].vtxwgt;
    Matched[v] = TRUE; /* if no matching vertex will be found,
                           the vertex is still considered matched
                           and the group has size 1 */
    vtxdeg = pHG->V[v].iEnd - pHG->V[v].iStart; /* Calculate vertex degree. */
    
    if (!MoveVtxInNetAdjncy(pHG, v)) return FALSE; /* move the matched vertex to the
                                                      end (i.e. >= iStartP1) in its
                                                      net adjacency lists */
    
    /* Compute the inner products by traversing the adjacent nets */
    nvisited = 0;
    
    for (t = pHG->V[v].iStart; t < pHG->V[v].iEnd; t++) {
        n = pHG->VtxAdjncy[t];
        netsize = pHG->N[n].iEnd - pHG->N[n].iStartP0;

        /* Compute scale factor of net n */
        switch (pOptions->Coarsening_NetScaling) {
            case NoNetScaling:
                scrowwgt = 1.0;
                break;
            case NetSclLinear:
                if (netsize > 0) {
                    scrowwgt = 1.0/(double)netsize;
                } else {
                    scrowwgt = 1.0;
                }
                break;
            default:
                fprintf(stderr, "FindMatchInprod(): unknown net scaling!\n");
                return FALSE;
        }

        /* Traverse net n */
        for (tt = pHG->N[n].iStartP0; tt < pHG->N[n].iStartP1; tt++) {
            v2 = pHG->NetAdjncy[tt]; /* v2<>v by preceding MoveVtxInNetAdjncy
                                         of v */
            /* Register new visits to vertices and store them in Visited */
            if (Inprod[v2] == 0) {
                Visited[nvisited] = v2;
                nvisited++; /* number of visited vertices */
            }
            Inprod[v2]++;
            ScInprod[v2] += scrowwgt;
        }
    }

    /* Match identical vertices first */
    if (pOptions->Coarsening_MatchIdenticalFirst == MatchIdYes) {
        /* FIXME: I think this code can be made more efficient by escaping the for-loop as soon as pC->Start[k+1] - pC->Start[k] >= pC->MaxNrVertices. */
        for (t=0; t<nvisited; t++) {
            v2 = Visited[t];
            
            vtxdeg2 = pHG->V[v2].iEnd - pHG->V[v2].iStart;

            /* Add v2 to the current match if
               - v2 is unmatched (guaranteed by construction)
               - total weight of match including v2 <= MaxVtxWgt
               - number of vertices in match excluding v2 < MaxNrVertices
               - inner product = deg(v) = deg(v2) */

            if (vtxwgt + pHG->V[v2].vtxwgt <= pC->MaxVtxWgt &&
                 pC->Start[k+1] - pC->Start[k] < pC->MaxNrVertices &&
                 Inprod[v2] == vtxdeg &&
                 Inprod[v2] == vtxdeg2) {

                pC->Match[pC->Start[k+1]] = v2;
                pC->Start[k+1]++;
                Matched[v2] = TRUE;
                MoveVtxInNetAdjncy(pHG, v2);
                vtxwgt += pHG->V[v2].vtxwgt;
  
                /* Reset v2 and remove it from Visited */
                Inprod[v2] = 0;
                ScInprod[v2]= 0;
                Visited[t] = Visited[nvisited-1];
                nvisited--;
            }
        }
    }

    /* Match the best vertices, adding at most
       MaxNrVertices - 1 vertices to the vertex v */
    while (pC->Start[k+1] - pC->Start[k] < pC->MaxNrVertices) {
        /* Find the best remaining unmatched vertex v2 */
        maxip = -1;
        tmax = -1;
        
        /* FIXME: I think this code can be made more efficient by escaping the for-loop as soon as pC->Start[k+1] - pC->Start[k] >= pC->MaxNrVertices. */
        for (t = 0; t < nvisited; t++) {
            v2 = Visited[t];
            
            vtxdeg2 = pHG->V[v2].iEnd - pHG->V[v2].iStart;
            
            mindeg = MIN(vtxdeg, vtxdeg2);
            maxdeg = MAX(vtxdeg, vtxdeg2);
            
            switch (pOptions->Coarsening_InprodScaling) {
                case NoIpScaling:
                    sccolwgt = 1.0;
                    break;
                case IpSclCos:
                    sccolwgt = 1.0 / sqrt((double)mindeg*maxdeg);
                    break;
                case IpSclMin:
                    sccolwgt = 1.0 / (double)mindeg;
                    break;
                case IpSclMax:
                    sccolwgt = 1.0 / (double)maxdeg;
                    break;
                case IpSclJaccard: /* the Jaccard metric is the ratio between the
                                 number of nonzeros in the intersection of the 
                                 two columns and the number of nonzeros
                                 in the union */
                    sccolwgt = 1.0 / (double)(mindeg+maxdeg-Inprod[v2]);
                    break;
                default:
                    fprintf(stderr, "FindMatchInprod(): unknown inner product scaling!\n");
                    return FALSE;
            }

            ScInprod[v2] *= sccolwgt;

            if (ScInprod[v2] > maxip &&
                vtxwgt + pHG->V[v2].vtxwgt <= pC->MaxVtxWgt) {
                maxip = ScInprod[v2];
                tmax = t;
            }
        }   
  
        if (maxip == -1)
            break; /* no more candidate vertices available */
  
        /* Best match is v2 */
        v2 = Visited[tmax];

        /* Add v2 to the current match if
           - v2 is unmatched (guaranteed by construction)
           - total weight of match including v2 <= MaxVtxWgt (also guaranteed)
           - number of vertices in match excluding v2 < MaxNrVertices (guaranteed by outer loop) */
        pC->Match[pC->Start[k+1]] = v2;
        pC->Start[k+1]++;
        Matched[v2] = TRUE; 
        vtxwgt += pHG->V[v2].vtxwgt;
        
        if (!MoveVtxInNetAdjncy(pHG, v2)) return FALSE;
        
        /* Reset v2 and remove it from Visited */
        Inprod[v2] = 0;
        ScInprod[v2] = 0;
        Visited[tmax]= Visited[nvisited-1];
        nvisited--;
    }
  
    /* Reset Inprod and ScInprod for remaining vertices */
    for (t = 0; t < nvisited; t++) {
        v2 = Visited[t];
        Inprod[v2] = 0;
        ScInprod[v2] = 0;
    }
    
    return TRUE;
} /* end FindMatchInprod */

