#include "HKLFM.h"

#ifdef USE_PATOH
#include <patoh.h>
#endif

int CreateInitialBalancedPartition(struct biparthypergraph *pHG, long weightlo, long weighthi) {

    /* This function creates an initial balanced partition of
       the hypergraph, splitting it into two parts,
       such that upper bounds for the vertex weight of the smallest
       and largest parts are obeyed.

       The vertices are assigned greedily in order of decreasing weights.
       Vertices with equal weights are first randomly permuted,
       to increase the randomness of the split.

       The function also computes the resulting weights of the two parts,

       Input: weightlo = smallest upper bound for part weight, belongs to part 0
              weighthi = largest upper bound for part weight, belongs to part 1

       Output: pHG->WeightP[0] = weight of part 0 
               pHG->WeightP[1] = weight of part 1 

    */

    int P;
    long t, left, right, maxvtxwgt, totvtxwgt,
         surplus, ideallo, idealhi, *Weight, *iv;
    
    if (!pHG) {
        fprintf(stderr, "CreateInitialBalancedPartition(): Null argument!\n");
        return FALSE;
    }
  
    if (pHG->NrVertices == 0) {
        pHG->WeightP[0] = 0;
        pHG->WeightP[1] = 0;
        return TRUE;
    }

    /* Allocate memory */
    Weight = (long *) malloc(pHG->NrVertices * sizeof(long));
    
    if (Weight == NULL) {
        fprintf(stderr, "CreateInitialBalancedPartition(): Not enough memory!\n");
        return FALSE;
    }
  
    /* Copy the vertex weights into array Weight
       and compute their maximum and sum */
    /* TODO: Why calculate maxvtxwgt? */
    maxvtxwgt = 0;
    totvtxwgt = 0;
    
    for (t = 0; t < pHG->NrVertices; t++) {
        Weight[t] = pHG->V[t].vtxwgt;
        if (Weight[t] > maxvtxwgt)
            maxvtxwgt = Weight[t];
        totvtxwgt += Weight[t];
    }

    /* Compute surplus and ideal weights */
    surplus = weightlo + weighthi - totvtxwgt;
    
    if (surplus < 0) {
        /* Desired split is infeasible */
        fprintf(stderr, "CreateInitialBalancedPartition(): Desired split is infeasible!\n");
        return FALSE;
    }

    ideallo = MAX(weightlo - surplus/2, 0); /* ideal weight of smallest part */
    idealhi = MAX(weighthi - surplus/2, 0); /* ideal weight of largest part */
  
    /* Sort Weight in decreasing order
       and allocate bookkeeping array iv  */
    iv = QSort(Weight, pHG->NrVertices);
    
    if (iv == NULL) {
        fprintf(stderr, "CreateInitialBalancedPartition(): Unable to sort weights!\n");
        return FALSE;
    }
        
    /* Randomly permute each group of vertices with equal weights */
    left = right = 0;
    while (left < pHG->NrVertices) {
        right = left;
        while (right < pHG->NrVertices && Weight[right] == Weight[left])
            right++;
  
        /* Randomly permute vertices left..right-1, which have equal weights */
        RandomPermute(iv, NULL, NULL, NULL, left, right-1);

        left = right;
    }
  
    /* Assign vertices to one of the parts in a greedy way */
    pHG->WeightP[0] = 0;
    pHG->WeightP[1] = 0;
    for (t = 0; t < pHG->NrVertices; t++) {
        
        /* Compute the part P that is farthest from its ideal weight */
        if (ideallo - pHG->WeightP[0] > idealhi -  pHG->WeightP[1])
            P = 0;
        else if (ideallo - pHG->WeightP[0] < idealhi -  pHG->WeightP[1])
            P = 1;
        else P = Random1(0,1); /* random tie-breaking */
        
        /* Assign vertex iv[t] to part P */
        pHG->V[iv[t]].partition = P;
        pHG->WeightP[P] += pHG->V[iv[t]].vtxwgt;
    }
    
    /* Free memory */
    if (iv != NULL)
        free(iv);
    free(Weight);

    return TRUE;
} /* end CreateInitialBalancedPartition */
  
  
int ComputeInitialGains(struct biparthypergraph *pHG) {

    /* This function computes the initial gain of each vertex of the hypergraph.
       The gain is defined as the decrease in the number of cut nets
       (communication volume) obtained by moving the vertex from its part
       (0 or 1) to the other part. The vertex is then inserted
       in the gainbucket of the appropriate gain.

       Second, the function orders the pins (nonzeros) of each net
       such that those belonging to part 0 come first, and are followed
       by those of part 1. It sets iStartP1 of each net accordingly.
       The pins of part 0 are thus in positions iStartP0..iStartP1-1, 
       and those of part 1 are in positions iStartP1..iEnd-1.
      
       Third, this function computes the initial communication volume
       and part weights.

    */

    long v, n, t, gain, nz0, nz1;
    int P;
    long nF; /* number of vertices in the FROM part of the net */
    long nT; /* number of vertices in the TO part of the net */
    
    if (!pHG) {
        fprintf(stderr, "ComputeInitialGains(): Null argument!\n");
        return FALSE;
    }
    
    /* Sort the vertices for each net on part number and set iStartP1 */
    for (n = 0; n < pHG->NrNets; n++) {
        pHG->N[n].iStartP1 = pHG->N[n].iEnd;
        t = pHG->N[n].iStartP0;
  
        while (t < pHG->N[n].iStartP1) {
            if (pHG->V[pHG->NetAdjncy[t]].partition == 1) {
                /* Swap NetAdjncy[t] <-> NetAdjncy[iStartP1-1] */
                (pHG->N[n].iStartP1)--;
                SwapLong(pHG->NetAdjncy, t, pHG->N[n].iStartP1);
            } else /* partition = 0 */
                t++;
        }
    }
  
    /* Compute the initial gains */
    for (v = 0; v < pHG->NrVertices; v++) {
        P = pHG->V[v].partition;
  
        /* Calculate the vertex gain from the adjacent critical nets */
        gain = 0;
        for (t = pHG->V[v].iStart; t < pHG->V[v].iEnd; t++) {
            n = pHG->VtxAdjncy[t];

            /* Compute number of pins (nonzeros) in parts 0, 1 of net n */
            nz0 = pHG->N[n].iStartP1 - pHG->N[n].iStartP0;
            nz1 = pHG->N[n].iEnd - pHG->N[n].iStartP1;
  
            if (P == 0) {
                nF = nz0;
                nT = nz1;
            } else {
                nF = nz1;
                nT = nz0;
            }
      
            if (nF == 1)
                gain++; /* win 1 by moving the only vertex
                            away FROM this net part */
            if (nT == 0)
                gain--; /* lose 1 by moving a vertex
                            TO an empty net part */
        }
  
        /* Insert the gain in the appropriate gainbucket of part P */
        if (!(pHG->V[v].GBentry = BucketInsert(&(pHG->GBVtx[P]), gain, v))) {
            fprintf(stderr, "ComputeInitialGains(): Unable to insert gain!\n");
            return FALSE;
        }
    }
  
    /* Calculate the initial communication */
    pHG->CurComm = 0;
    for (n = 0; n < pHG->NrNets; n++) {
        nz0 = pHG->N[n].iStartP1 - pHG->N[n].iStartP0;
        nz1 = pHG->N[n].iEnd - pHG->N[n].iStartP1;

        if (nz0 > 0  && nz1 > 0) 
            pHG->CurComm++;
    }

    /* Calculate the initial weights */
    pHG->WeightP[0] = 0;
    pHG->WeightP[1] = 0;
    for (v = 0; v < pHG->NrVertices; v++)
        pHG->WeightP[pHG->V[v].partition] += pHG->V[v].vtxwgt;

    return TRUE;
} /* end ComputeInitialGains */


int ClearMoveLog(struct biparthypergraph *pHG) {
    
    /* This function clears the array VtxMoveLog.
       There is no need to clear the values, just to reset
       the counters pointing into the array to 0.

       CurVtxLog = the current number of vertices in the array, 
                   i.e. the number of vertices that have been moved.
                   VtxMoveLog[0..CurVtxLog-1] thus contains the
                   moved vertices.

       MinVtxLog = the current number of vertices at the time of the lowest
                   communication volume. All vertices moved after that
                   will have to be moved back.
                   VtxMoveLog[0..MinVtxLog-1] thus contains the
                   vertices that remain in their new part.
    */
    if (!pHG) {
        fprintf(stderr, "ClearMoveLog(): Null argument!\n");
        return FALSE;
    }
  
    pHG->CurVtxLog = 0;
    pHG->MinVtxLog = 0;
  
    return TRUE;
} /* end ClearMoveLog */
  
  
int MoveVertex(struct biparthypergraph *pHG, long v) {

    /* This function moves vertex v from its old part P
       to its new part 1-P, where P = 0, 1,
       i.e., it is moved from part 0 to 1 or vice versa.
       The weights of both part are updated accordingly.

       Furthermore, the vertex is locked by removing its pointer
       to the gainbucket data structure. (To keep this data structure
       consistent, the vertex must also be deleted from the structure.)
    */
    int P;
    
    if (!pHG) {
        fprintf(stderr, "MoveVertex(): Null argument!\n");
        return FALSE;
    }
    
    if (v < 0) {
        fprintf(stderr, "MoveVertex(): Invalid vertex index!\n");
        return FALSE;
    }
    
    P = pHG->V[v].partition;
  
    /* Decrease the weight for the old part */
    pHG->WeightP[P]  -= pHG->V[v].vtxwgt;
  
    /* Increase the weight for the new part */
    pHG->WeightP[1-P] += pHG->V[v].vtxwgt;
  
    /* Update the part number */
    pHG->V[v].partition = 1-P;
      
    /* Lock the vertex */
    pHG->V[v].GBentry = NULL;

    return TRUE;
} /* end MoveVertex */


int UpdateGains(struct biparthypergraph *pHG, const long v) {

    /* This function updates the gains of all vertices
       that are affected by a move of vertex v,
       immediately after the vertex has been moved,
       but before the pins have been updated to reflect the move.

       The function does this by examining all the nets adjacent to v, 
       checking whether they are critical, i.e,
       the new part has 0 or 1 pins before the move
       or the old part has 1 or 2 pins before the move
       (leaving 0 or 1 after the move). If so, the gains of the
       still movable vertices in the net are adjusted accordingly.

       Second, the function moves the pin (nonzero) corresponding
       to vertex v in each affected net, to restore the situation in which 
       the pins belonging to part 0 come first, and are followed
       by those of part 1. It sets iStartP1 of each affected net accordingly.
       The pins of part 0 are thus in positions iStartP0..iStartP1-1,
       and those of part 1 are in positions iStartP1..iEnd.

       Third, this function computes the current communication volume.
    */

    long n, nz0, nz1, vadj, newval, t, tt;
    int P, Padj;
    struct bucketentry *pE;
    long nF; /* number of vertices in the FROM (old) part of the net,
                 before the move */
    long nT; /* number of vertices in the TO (new) part of the net,
                 before the move */
    
    if (!pHG) {
        fprintf(stderr, "UpdateGains(): Null argument!\n");
        return FALSE;
    }
  
    P = pHG->V[v].partition; /* the new part of vertex v = TO part */
  
    /* Examine all the affected nets */
    for (t = pHG->V[v].iStart; t < pHG->V[v].iEnd; t++) {
        n = pHG->VtxAdjncy[t];
  
        /* Compute the number of pins of parts 0, 1 before the move */
        nz0 = pHG->N[n].iStartP1 - pHG->N[n].iStartP0;
        nz1 = pHG->N[n].iEnd - pHG->N[n].iStartP1;

        /* Quick adjustment for frequent case of singleton net */
        if (nz0 + nz1  == 1) {
            if (P == 0)
                pHG->N[n].iStartP1++;
            else
                pHG->N[n].iStartP1--;
            continue; /* with the next net */
        }

        /* Determine number of pins of TO and FROM part */
        if (P == 0) {
            /* vertex has been moved to part 0 */
            nT = nz0;
            nF = nz1;
        } else {
            /* vertex has been moved to part 1 */
            nT = nz1;
            nF = nz0;
        }
  
        /* Check the nets that are critical in the TO part */
        if (nT == 0) {
            /* Increment the gains of all unlocked vertices,
               because they can now move into the TO part
               without incurring additional cost, since v is already there.  */
            for (tt = pHG->N[n].iStartP0; tt < pHG->N[n].iEnd; tt++) {
                vadj = pHG->NetAdjncy[tt];
                pE = pHG->V[vadj].GBentry;
        
                if (pE != NULL) {
                    /* vertex is not locked */
                    newval = (pE->bucket)->value + 1;
                    Padj = pHG->V[vadj].partition;
                    
                    if (!(pHG->V[vadj].GBentry = BucketMove(&(pHG->GBVtx[Padj]), pE, newval))) {
                        fprintf(stderr, "UpdateGains(): Unable to move bucket!\n");
                        return FALSE;
                    }
                }
            }
        }

        if (nT == 1) {
            /* Decrement the gain of the single vertex in the TO part,
               provided it is unlocked, because it cannot gain any more 
               by moving away. */
            if (P == 0)
                vadj = pHG->NetAdjncy[pHG->N[n].iStartP0];
            else
                vadj = pHG->NetAdjncy[pHG->N[n].iStartP1];
      
            pE = pHG->V[vadj].GBentry;
            
            if (pE != NULL) { 
                newval = (pE->bucket)->value - 1;
                Padj = pHG->V[vadj].partition;
                
                if (!(pHG->V[vadj].GBentry = BucketMove(&(pHG->GBVtx[Padj]), pE, newval))) {
                    fprintf(stderr, "UpdateGains(): Unable to move bucket!\n");
                    return FALSE;
                }
            }
        } /* else, if nT > 1, the gains cannot change due to the TO part */
  
        /* Update the pins to reflect the move */

        if (P == 0) {
            /* Find the index of v */
            tt = pHG->N[n].iStartP1;
            
            while (pHG->NetAdjncy[tt] != v && tt < pHG->N[n].iEnd) tt++;
            
            if (tt == pHG->N[n].iEnd) {
                fprintf(stderr, "UpdateGains(): Vertex not found!\n");
                return FALSE;
            }
  
            /* Swap iStartP1 <-> tt */
            SwapLong(pHG->NetAdjncy, pHG->N[n].iStartP1, tt);
            pHG->N[n].iStartP1++;

        } else {
            /* Find the index of v */
            tt = pHG->N[n].iStartP0;
            
            while (pHG->NetAdjncy[tt] != v && tt < pHG->N[n].iStartP1) tt++;
            
            if (tt == pHG->N[n].iStartP1) {
                fprintf(stderr, "UpdateGains(): Vertex not found!\n");
                return FALSE;
            }

            /* Swap iStartP1-1 <-> tt */
            pHG->N[n].iStartP1--;
            SwapLong(pHG->NetAdjncy, pHG->N[n].iStartP1, tt);
        }
  
        /* Check the nets that are critical in the FROM part */
        if (nF == 1) {
            /* Decrement the gains of all unlocked vertices,
               because if they move into the FROM part, they will
               incur the additional cost of moving into an empty part. */
            for (tt = pHG->N[n].iStartP0; tt < pHG->N[n].iEnd; tt++) {
                vadj = pHG->NetAdjncy[tt];
                pE = pHG->V[vadj].GBentry;
  
                if (pE != NULL) {
                    newval = (pE->bucket)->value - 1;
                    Padj = pHG->V[vadj].partition;
                    
                    if (!(pHG->V[vadj].GBentry = BucketMove(&(pHG->GBVtx[Padj]), pE, newval))) {
                        fprintf(stderr, "UpdateGains(): Unable to move bucket!\n");
                        return FALSE;
                    }
                }
            }
        }
        if (nF == 2) {
            /* Increment the gain of the single remaining vertex
               in the FROM part, provided it is unlocked,
               because it can now gain by moving away. */

            if (P == 0)
                vadj = pHG->NetAdjncy[pHG->N[n].iStartP1];
            else
                vadj = pHG->NetAdjncy[pHG->N[n].iStartP0];
  
            pE = pHG->V[vadj].GBentry;
        
            if (pE != NULL) {
                newval = (pE->bucket)->value + 1;
                Padj = pHG->V[vadj].partition;
                
                if (!(pHG->V[vadj].GBentry = BucketMove(&(pHG->GBVtx[Padj]), pE, newval))) {
                    fprintf(stderr, "UpdateGains(): Unable to move bucket!\n");
                    return FALSE;
                }
            }
        } /* else, if nF > 2, the gains cannot change due to the FROMpart */

        /* Update the communication count */

        if (nT == 0) {
            /* Increment the current communication,
               caused by moving v into an empty part */
            pHG->CurComm++;
        }
        if (nF == 1) {
            /* Decrement the current communication,
               caused by moving v away and leaving an empty part */
            pHG->CurComm--;
        }
    }

    return TRUE;
} /* end UpdateGains */


int HKLFM(struct biparthypergraph *pHG, long weightlo, long weighthi,
           int MaxNrLoops, int MaxNrNoGainMoves) {

    /* This function performs one run of the
       Hypergraph Kernighan-Lin Fiduccia-Mattheyses (HKLFM) algorithm,
       which performs at most MaxNrLoops iterations.
       Each iteration tries to move all vertices to the other
       part, in order of decreasing communication gain of the move.

       After a vertex has been moved, it is locked and cannot
       move anymore in the same iteration. After a move, all
       gains are updated. 

       Input: weightlo = smallest upper bound for part weight, belongs to part 0
              weighthi = largest upper bound for part weight, belongs to part 1
              MaxNrLoops = maximum number of iterations of the main loop
              MaxNrNoGainMoves = maximum number of successive no-gain moves
                                 within one iteration

       Output: the best partition obtained so far is stored 
       in HG.OptPartVtx, together with its volume HG.OptComm .
    */

    int P;
    long iter, t, v, v0, v1, maxweight, gainP0, gainP1,
         LastComm = LONG_MAX, NrNoGainMoves;
    
    if (!pHG) {
        fprintf(stderr, "HKLFM(): Null argument!\n");
        return FALSE;
    }
    
    /* Main loop */
    iter = 0;
    while (iter < MaxNrLoops && (iter == 0 || pHG->MinComm != LastComm)) {
        /* One iteration of the main loop tries to move all vertices.
           If the lowest communication volume so far in this run equals
           the lowest volume one iteration ago, then stop the run,
           since further improvement is unlikely.  */ 

        /* Clear the initial gain data structures and unlock all vertices */
        ClearGainBucket(&(pHG->GBVtx[0]));
        ClearGainBucket(&(pHG->GBVtx[1]));
        ClearMoveLog(pHG);

        /* Compute the initial gains and CurComm,
           the current commmunication (in this iteration)  */
        if (!ComputeInitialGains(pHG)) {
            fprintf(stderr, "HKLFM(): Unable to compute initial gains!\n");
            return FALSE;
        }
  
        /* Initialise communication counters:
           OptComm = lowest communication so far in all runs of HKLFM
           MinComm = lowest commmunication so far in this run of HKLFM 
           LastComm = lowest commmunication so far in this run until the end
                      of the previous iteration */
   
        if (iter == 0) {
            /* This the first iteration of the main loop */
            if (pHG->OptComm == LONG_MAX) {
                /* This is the first run of HKLFM */
                pHG->OptComm = pHG->CurComm;
            }
            pHG->MinComm = pHG->CurComm;
        }
        LastComm = pHG->MinComm;
  
#ifdef INFO2
        /* Print current communication volume */
        printf(" %ld", pHG->CurComm);
        fflush(stdout);
#endif
        
        /* Try to move all vertices once */
        NrNoGainMoves = 0;
        while (pHG->GBVtx[0].NrBuckets + pHG->GBVtx[1].NrBuckets > 0 && 
            NrNoGainMoves < MaxNrNoGainMoves) {

            /* Try to move the unlocked vertex with the highest gain */

            v0 = GainBucketGetMaxValVertexNr(&(pHG->GBVtx[0]));
            v1 = GainBucketGetMaxValVertexNr(&(pHG->GBVtx[1]));
    
            if (v0 == LONG_MIN && v1 == LONG_MIN) {
                fprintf(stderr, "HKLFM(): No vertex to be moved!\n");
                return FALSE;
            }
  
            /* Determine part P from which to move a vertex */
            if (v0 != LONG_MIN && v1 == LONG_MIN) {
                P = 0;
            }
            else if (v0 == LONG_MIN && v1 != LONG_MIN) {
                P = 1;
            }
            else {
                /* Find part P with maximum gain */
                gainP0 = LONG_MIN;
                gainP1 = LONG_MIN;
        
                if (pHG->WeightP[1] + pHG->V[v0].vtxwgt <= weighthi) 
                    gainP0 = GainBucketGetMaxVal(&(pHG->GBVtx[0]));
                if (pHG->WeightP[0] + pHG->V[v1].vtxwgt <= weightlo) 
                    gainP1 = GainBucketGetMaxVal(&(pHG->GBVtx[1]));
        
                if (gainP0 > gainP1)
                    P = 0;
                else if (gainP1 > gainP0)
                    P = 1;
                else 
                    P = Random1(0,1); /* random tie-breaking */
            }

            /* Determine maximum weight of the other part, 1-P */
            maxweight = (P == 0 ? weighthi : weightlo); 

            /* Find vertex v with maximum gain */
            v = BucketDeleteMax(&(pHG->GBVtx[P]));
            
            if (v < 0) {
                fprintf(stderr, "HKLFM(): Unable to delete maximum vertex!\n");
                return FALSE;
            }

            /* Move v to the other part if its weight allows this */
            if (pHG->WeightP[1-P] + pHG->V[v].vtxwgt <= maxweight) {
                /* Move v from part P to part 1-P and lock it */
                MoveVertex(pHG, v);
                
                if (!UpdateGains(pHG, v)) {
                    fprintf(stderr, "HKLFM(): Unable to update gains!\n");
                    return FALSE;
                }
                
                pHG->VtxMoveLog[pHG->CurVtxLog++] = v;
            } else  {
                /* Lock v */
                pHG->V[v].GBentry = NULL;
            }
    
            /* If the current partition is the best so far, then record this */
            if (pHG->CurComm < pHG->MinComm) {
                pHG->MinComm = pHG->CurComm;
                pHG->MinVtxLog = pHG->CurVtxLog;
                NrNoGainMoves = 0; /* reset, since we had a gain move */
            } else if (pHG->CurComm == pHG->MinComm) {
                /* Choose the new partition, since this may have better
                   load balance, but do not reset */
                pHG->MinVtxLog = pHG->CurVtxLog;
                NrNoGainMoves++;
            } else
                NrNoGainMoves++;
        }
        
        /* Restore the best partition by undoing the moves done after
           the minimum communication has been achieved in this iteration */
        for (t = pHG->CurVtxLog-1; t >= pHG->MinVtxLog; t--) {
            v = pHG->VtxMoveLog[t];
            MoveVertex(pHG, v);
        }
  
        /* If this iteration has produced the lowest overall communication volume
           of all the runs, then save it. Also if the volume is equal to the best
           so far, to make sure a valid partition is saved in all cases. */
        if (pHG->MinComm <= pHG->OptComm) {
            pHG->OptComm = pHG->MinComm;
            for (t = 0; t < pHG->NrVertices; t++)
                pHG->OptPartVtx[t] = pHG->V[t].partition;
        }
        iter++;
    } /* end main loop */

    return TRUE;
} /* end HKLFM */


int RunHKLFM(struct biparthypergraph *pHG, long weightlo, long weighthi,
              int refine, const struct opts *pOptions) {

    /* This function performs several runs of the 
       Hypergraph Kernighan-Lin Fiduccia-Mattheyses (HKLFM) algorithm,
       each time starting from a different random balanced initial partition.
       If the option refine=TRUE, only one run is performed
       and the current partition is taken as initial partition.

       Input: weightlo = smallest upper bound for part weight, belongs to part 0
              weighthi = largest upper bound for part weight, belongs to part 1

    */

    long k, NrRestarts, MaxNrLoops, MaxNrNoGainMoves;
    
    if (!pHG || !pOptions) {
        fprintf(stderr, "RunHKLFM(): Null arguments!\n");
        return FALSE;
    }
    
#ifdef INFO2
    if (refine)
        printf("RefineHKLFM\n");
    else
        printf("RunHKLFM\n");
    fflush(stdout);
#endif
    
    /* OptComm = lowest communication in all runs of HKLFM
       MinComm = lowest commmunication in this run of HKLFM */
    pHG->MinComm = LONG_MAX;
    pHG->OptComm = LONG_MAX;
    
    if (refine) {
        NrRestarts = 1;
        MaxNrLoops = pOptions->KLFM_Refine_MaxNrLoops;
        MaxNrNoGainMoves = pOptions->KLFM_Refine_MaxNrNoGainMoves;
    } else {
        NrRestarts = pOptions->KLFM_InitPart_NrRestarts;
        MaxNrLoops = pOptions->KLFM_InitPart_MaxNrLoops;
        MaxNrNoGainMoves = pOptions->KLFM_InitPart_MaxNrNoGainMoves;
    }

    for (k = 0; k < NrRestarts; k++) { 
        /* Perform one run of the HKLFM algorithm, starting with
           a new random balanced initial partition */

        if (refine == FALSE) {
            if (!CreateInitialBalancedPartition(pHG, weightlo, weighthi)) {
                fprintf(stderr, "RunHKLFM(): Unable to create initial partition!\n");
                return FALSE;
            }
        }
        
        if (!HKLFM(pHG, weightlo, weighthi, MaxNrLoops, MaxNrNoGainMoves)) {
            fprintf(stderr, "RunHKLFM(): Unable to run the HKLFM algorithm!\n");
            return FALSE;
        }
          
#ifdef INFO2
        printf(" %ld (Opt: %ld)\n", pHG->MinComm, pHG->OptComm);
#endif      
    }  
    
#ifdef INFO2
    printf(" Communication: %ld\n", pHG->OptComm);
#endif

    return TRUE;
} /* end RunHKLFM */


int CoarsenGraph(struct biparthypergraph *pHG,
                 struct biparthypergraph *pCHG, struct contraction *pC,
                 long level, const struct opts *pOptions) {
                                
    /* This function matches vertices from the hypergraph HG, 
       then contracts each matched group of vertices into one vertex,
       giving the contracted hypergraph CHG.

       This is done as follows.  First the vertices to be matched
       are ordered, where the options are:
       by decreasing or increasing vertex weight;
       by decreasing or increasing vertex degree; 
       by natural order, or by random order. 

       After the ordering has been determined, each vertex to be matched
       is compared with the adjacent vertices for similarity.
       In case of the inner-product similarity measure,
       the vertex with the highest nonzero overlap with the vertex
       to be matched is chosen. If the hypergraph is the finegrain
       hypergraph corresponding to a sparse matrix, a specialised faster
       matching function can be used.

       Finally, each group of matched vertices is merged into a single
       vertex of the contracted hypergraph. This is done by adding
       the vertex weights, merging the pins (nonzeros) giving
       new vertex adjacency lists, generating corresponding net adjacency
       lists, and copying the net weights and directions.
       
       In addition to the above, free vertices are made to disappear by
       matching them with the first match group and setting their weights
       within the group to 0. Free nets are discarded as well by simply
       not adding them to the coarsened graph.
       
       Input:  level is the level of the matching, i.e. the number of times
               the original hypergraph has been contracted.
       Output: The number of vertices of CHG <= the number of vertices of HG.
               The number of nets of CHG = the number of nets of CHG.

    */

    long *VtxWeight; /* stores the weights or degrees of the vertices */
    long *iv = NULL; /* stores the order in which the vertices to be matched
                         are visited */
    long *Ip; /* working array storing inner products of vertex pairs */
    double *ScIp; /* working array storing scaled inner products */
    long *Visited; /* stores the vertices that have been visited
                       when looking for a match for a vertex v */
    int *Matched; /* Matched[v] = TRUE means vertex v has already
                      been matched */
    long *iAdjncy; /* arrays of length NrNets */
    long k, t, tt, v, n, nrpins, totalvtxwgt;
    
    long NumFreeNets, NumFreeVertices, *NetIndices; /* Counts and indices of free nets/vertices. */
    
    if (!pHG || !pCHG || !pC || !pOptions) {
        fprintf(stderr, "CoarsenGraph(): Null arguments!\n");
        return FALSE;
    }
    
    /* Allocate memory. */
    pC->Match = (long *) malloc(pHG->NrVertices * sizeof(long));
    pC->Start = (long *) malloc((pHG->NrVertices+1) * sizeof(long));
    
    if (pC->Match == NULL || pC->Start == NULL) {
        fprintf(stderr, "CoarsenGraph(): Not enough memory for Match/Start!\n");
        return FALSE;
    }
    
    /* Determine free nets and vertices. */
    NumFreeNets = 0;
    NumFreeVertices = 0;
    
    NetIndices = (long *)malloc(pHG->NrNets*sizeof(long));
    
    if (NetIndices == NULL) {
        fprintf(stderr, "CoarsenGraph(): Not enough memory!\n");
        return FALSE;
    }
    
    if (pOptions->DiscardFreeNets == FreeNetYes) {
        tt = 0;
        
        /* Determine singleton nets and count all free nets. */
        for (t = 0; t < pHG->NrNets; ++t) {
            if (pHG->N[t].Free || pHG->N[t].iStartP0 + 1 >= pHG->N[t].iEnd) {
                pHG->N[t].Free = TRUE;
                NetIndices[t] = pHG->NrNets; /* Set to an invalid index to indicate that this net is free. */
                ++NumFreeNets;
            }
            else {
                NetIndices[t] = tt++;
            }
        }
        
        /* Determine and count free vertices. */
        for (v = 0; v < pHG->NrVertices; ++v) {
            pHG->V[v].Free = TRUE;
            
            for (t = pHG->V[v].iStart; t < pHG->V[v].iEnd; ++t) {
                if (!pHG->N[pHG->VtxAdjncy[t]].Free) {
                    pHG->V[v].Free = FALSE;
                    break;
                }
            }
            
            if (pHG->V[v].Free) ++NumFreeVertices;
        }
    }
    else {
        /* Otherwise mark all vertices as non-free. */
        for (t = 0; t < pHG->NrVertices; ++t) pHG->V[t].Free = FALSE;
        for (t = 0; t < pHG->NrNets; ++t) NetIndices[t] = t;
    }
    
#ifdef INFO2
    printf("    Removed %ld/%ld free vertices and %ld/%ld free nets.\n", NumFreeVertices, pHG->NrVertices, NumFreeNets, pHG->NrNets);
#endif
    
    /* Compute the total vertex weight, without the free vertices, which will be removed from the graph. */ 
    totalvtxwgt = 0;
    
    for (t = 0; t < pHG->NrVertices; t++) 
        if (!pHG->V[t].Free) totalvtxwgt += pHG->V[t].vtxwgt;
    
    /*-------------------------------------------------------------------------
      ## Vertex ordering stage ##
      -------------------------------------------------------------------------*/
    if (level < pOptions->Coarsening_NrMatchArbitrary) {
        /* We do nothing. */
        iv = NULL;
    } else if (pOptions->Coarsening_MatchingStrategy == MatchRandom) {
    
        /* Allocate memory for iv */
        iv = (long *) malloc(pHG->NrVertices * sizeof(long));
        
        if (iv == NULL) {
            fprintf(stderr, "CoarsenGraph(): Not enough memory for iv!\n");
            return FALSE;
        }
        
        /* Fill iv with the identity permutation */
        for (t = 0; t < pHG->NrVertices; t++)
            iv[t] = t;
            
        /* Permute iv to visit vertices randomly */
        RandomPermute(iv, NULL, NULL, NULL, 0, pHG->NrVertices-1);
            
        /* Permute the vertex adjacency lists randomly to create even more randomness */
        for (t = 0; t < pHG->NrVertices; t++) {
            if (pHG->MatReValue == NULL || pHG->SplitDir == FINEGRAIN || pHG->SplitDir == SFINEGRAIN)
                /* pattern matrix or numerical values not stored;
                   or finegrain split direction, where numerical values
                   are stored with the vertices themselves, not with the adjacencies */
                RandomPermute(pHG->VtxAdjncy, NULL, NULL, NULL,
                              pHG->V[t].iStart, pHG->V[t].iEnd-1);
            else if (pHG->MatImValue == NULL)
                /* real matrix with numerical values stored */
                RandomPermute(pHG->VtxAdjncy, NULL, pHG->MatReValue, NULL, 
                              pHG->V[t].iStart, pHG->V[t].iEnd-1);
            else
                /* complex matrix with numerical values stored */
                RandomPermute(pHG->VtxAdjncy, NULL, pHG->MatReValue, pHG->MatImValue, 
                              pHG->V[t].iStart, pHG->V[t].iEnd-1);
        }
                    
    } else if (pOptions->Coarsening_MatchingStrategy == MatchInprod ||
                pOptions->Coarsening_MatchingStrategy == MatchATA) {
    
        /* Allocate memory for VtxWeight */
        VtxWeight = (long *) malloc(pHG->NrVertices * sizeof(long));
        
        if (VtxWeight == NULL) {
            fprintf(stderr, "CoarsenGraph(): Not enough memory for VtxWeight!\n");
            return FALSE;
        }
  
        if (pOptions->Coarsening_InprodMatchingOrder == DecreasingWgt ||
             pOptions->Coarsening_InprodMatchingOrder == IncreasingWgt) {
            /* Fill VtxWeight with weights */
            for (t = 0; t < pHG->NrVertices; t++) 
                VtxWeight[t] = pHG->V[t].vtxwgt;
            iv = QSort(VtxWeight, pHG->NrVertices);
            
            if (iv == NULL) {
                fprintf(stderr, "CoarsenGraph(): Unable to sort!\n");
                return FALSE;
            }
        } else if (pOptions->Coarsening_InprodMatchingOrder == DecreasingDegree ||
                    pOptions->Coarsening_InprodMatchingOrder == IncreasingDegree)  {
            /* Fill VtxWeight with degrees */
            for (t = 0; t < pHG->NrVertices; t++) 
                VtxWeight[t] = pHG->V[t].iEnd - pHG->V[t].iStart; /* degree of vertex t */
            iv = QSort(VtxWeight, pHG->NrVertices);
            
            if (iv == NULL) {
                fprintf(stderr, "CoarsenGraph(): Unable to sort!\n");
                return FALSE;
            }
        } else if (pOptions->Coarsening_InprodMatchingOrder == NaturalOrder ||
                    pOptions->Coarsening_InprodMatchingOrder == RandomOrder) {
            iv = (long *) malloc(pHG->NrVertices * sizeof(long));
            
            if (iv == NULL) {
                fprintf(stderr, "CoarsenGraph(): Not enough memory for iv!\n");
                return FALSE;
            }

            for (t = 0; t < pHG->NrVertices; t++) 
                iv[t] = t;
        } else {
            fprintf(stderr, "CoarsenGraph(): Internal error: Unknown ordering strategy!\n");
            return FALSE;
        }
        
        if (pOptions->Coarsening_InprodMatchingOrder == IncreasingWgt ||
             pOptions->Coarsening_InprodMatchingOrder == IncreasingDegree)  {
            /* Reverse the order */
            for (t = 0; t < pHG->NrVertices/2; t++) 
                SwapLong(iv, t, pHG->NrVertices-1-t);
        }
        
        if (pOptions->Coarsening_InprodMatchingOrder == RandomOrder) 
            RandomPermute(iv, NULL, NULL, NULL, 0, pHG->NrVertices-1);
      
        free(VtxWeight);
        
    } else {
        fprintf(stderr, "CoarsenGraph(): Internal error: Unknown matching strategy!\n");
        return FALSE;
    }
             
  
    /*-------------------------------------------------------------------------
      ## Matching stage ##
      -------------------------------------------------------------------------*/
  
    pC->Start[0] = 0;
    pC->NrMatches = 0;
    pC->MaxNrVertices = pOptions->Coarsening_MaxNrVtxInMatch;
    pC->MaxVtxWgt = totalvtxwgt * pOptions->Coarsening_VtxMaxFractionOfWeight;
  
    /* Set HG.iStartP1 equal to iEnd.
       The array iStartP1 of HG is abused inside the matching functions */
    for (n = 0; n < pHG->NrNets; n++)
        pHG->N[n].iStartP1 = pHG->N[n].iEnd;
    
    Matched = (int *) malloc(pHG->NrVertices * sizeof(int));
    
    if (Matched == NULL) {
        fprintf(stderr, "CoarsenGraph(): Not enough memory!\n");
        return FALSE;
    }
    
    for (t = 0; t < pHG->NrVertices; t++) {
        Matched[t] = FALSE;
    }
    
    /* Make sure that all free vertices are not eligible for matching by moving
       them to the back of their adjacency lists. */
    for (t = 0; t < pHG->NrVertices; ++t) {
        if (pHG->V[t].Free) {
            if (!MoveVtxInNetAdjncy(pHG, t)) {
                fprintf(stderr, "CoarsenGraph(): Unable to move vertex!\n");
                return FALSE;
            }
        }
    }
    
    /* Match vertices, but make sure to only include vertices that are non-free. */
    if (level < pOptions->Coarsening_NrMatchArbitrary) {
        for (v = 0; v < pHG->NrVertices; v++) {
            if (Matched[v] == FALSE && !pHG->V[v].Free) {
                if (!FindMatchArbitrary(pHG, pC, v, Matched)) {
                    fprintf(stderr, "CoarsenGraph(): Unable to find match!\n");
                    return FALSE;
                }
            }
        }
    } else if (pOptions->Coarsening_MatchingStrategy == MatchRandom ||
        (level < pOptions->Coarsening_FineSwitchLevel &&
        (pHG->SplitDir == FINEGRAIN || pHG->SplitDir == SFINEGRAIN))) {
        
        /* Sort adjacency lists if we use random matching. */
        if (!SortAdjacencyLists(pHG)) {
            fprintf(stderr, "CoarsenGraph(): Unable to sort adjacency lists!\n");
            return FALSE;
        }

        for (t = 0; t < pHG->NrVertices; t++) {
            v = iv[t];
            if (Matched[v] == FALSE && !pHG->V[v].Free) {
                if (!FindMatchArbitrary(pHG, pC, v, Matched)) {
                    fprintf(stderr, "CoarsenGraph(): Unable to find match!\n");
                    return FALSE;
                }
            }
        }
    } else if (pOptions->Coarsening_MatchingStrategy == MatchInprod) {
        Ip = (long *) malloc(pHG->NrVertices * sizeof(long));
        ScIp = (double *) malloc(pHG->NrVertices * sizeof(double));
        Visited = (long *) malloc(pHG->NrVertices * sizeof(long));

        if (Ip == NULL || ScIp == NULL || Visited == NULL) {
            fprintf(stderr, "CoarsenGraph(): Not enough memory!\n");
            return FALSE;
        }
  
        /* Initialise the working arrays that need to be initialised */
        for (t = 0; t < pHG->NrVertices; t++) {
            Ip[t] = 0;
            ScIp[t] = 0.0;
        }

        for (t = 0; t < pHG->NrVertices; t++) {
            v = iv[t];
            if (Matched[v] == FALSE && !pHG->V[v].Free) {
                if (!FindMatchInprod(pHG, pC, v, Matched, Visited, Ip, ScIp, pOptions)) {
                    fprintf(stderr, "CoarsenGraph(): Unable to find inproduct match!\n");
                    return FALSE;
                }
            }
        }
        
        free(Visited);
        free(ScIp);
        free(Ip);
    } else if (pOptions->Coarsening_MatchingStrategy == MatchATA) {
        /* Use a hybrid matching strategy which combines a graph matching algorithm with a hypergraph neighbor finding algorithm. */
        
        /* This is a pointer to the function of the used matching algorithm. */
        int (*HybridMatcher)(struct biparthypergraph *, struct contraction *,
                             const long *, int *,
                             const struct opts *,
                             int (*SetupData)(void **, struct biparthypergraph *, const struct opts *),
                             int (*FindNeighbor)(long *, double *,
                                                 const struct biparthypergraph *, const struct contraction *, const struct opts *,
                                                 const long, const int *,
                                                 void *, long *, long *, double *),
                             int (*FreeData)(void *)) = NULL;
        
        /* These are pointers to the functions of the used neighbor finding algorithm. */
        int (*SetupData)(void **, struct biparthypergraph *, const struct opts *) = NULL;
        int (*FindNeighbor)(long *, double *,
                            const struct biparthypergraph *, const struct contraction *, const struct opts *,
                            const long, const int *,
                            void *, long *, long *, double *) = NULL;
        int (*FreeData)(void *) = NULL;
        
        /* Select matcher. */
        if (pOptions->Coarsening_MatchingATAMatcher == MatchMatcherGreedy) {
            HybridMatcher = MatchUsingGreedy;
        } else if (pOptions->Coarsening_MatchingATAMatcher == MatchMatcherPGA) {
            HybridMatcher = MatchUsingPGA;
        } else {
            fprintf(stderr, "CoarsenGraph(): Selected hybrid matcher is not available!\n");
            return FALSE;
        }
        
        /* Select neighbor finder. */
        if (pOptions->Coarsening_MatchingATAFinder == MatchFinderInproduct) {
            SetupData = MatchInproductSetup;
            FindNeighbor = FindNeighborInproduct;
            FreeData = MatchInproductFree;
        } else if (pOptions->Coarsening_MatchingATAFinder == MatchFinderStairway) {
            SetupData = MatchStairwaySetup;
            FindNeighbor = FindNeighborStairway;
            FreeData = MatchStairwayFree;
        } else {
            fprintf(stderr, "CoarsenGraph(): Selected hybrid matcher is not available!\n");
            return FALSE;
        }
        
        /* Execute matching algorithm. */
        if (!HybridMatcher(pHG, pC,
                           iv, Matched,
                           pOptions,
                           SetupData,
                           FindNeighbor,
                           FreeData)) {
            fprintf(stderr, "CoarsenGraph(): Unable to create a hybrid matching!\n");
            return FALSE;
        }
    } else {
        fprintf(stderr, "CoarsenGraph(): No valid matching strategy!\n");
        return FALSE;
    }

#ifdef INFO
    /* Another sanity check. */
    for (t = 0; t < pHG->NrVertices; ++t) if (pHG->V[t].Free && Matched[t]) fprintf(stderr, "CoarsenGraph(): Matching sanity check failed!\n");

    /* Output total matching weight. */
    if (FALSE) {
        int *flags = (int *)malloc(pHG->NrNets*sizeof(int));
        long IntegerWeight = 0;
        double ScaledWeight = 0.0;
        
        if (pC->Start[0] != 0) {
            fprintf(stderr, "CoarsenGraph(): Invalid starting offset in pC!\n");
            return FALSE;
        }
        
        for (t = 0; t < pHG->NrNets; t++) flags[t] = 0;
        
        for (t = 0; t < pC->NrMatches; t++)
        {
            if (pC->Start[t + 1] == pC->Start[t] + 2) {
                const long v1 = pC->Match[pC->Start[t]];
                const long v2 = pC->Match[pC->Start[t] + 1];
                
                if (v1 < 0 || v2 < 0 || v1 >= pHG->NrVertices || v2 >= pHG->NrVertices) {
                    fprintf(stderr, "CoarsenGraph(): Invalid matched vertex indices!\n");
                    return FALSE;
                }
                
                /* Calculate inner product. */
                for (tt = pHG->V[v1].iStart; tt < pHG->V[v1].iEnd; tt++) flags[pHG->VtxAdjncy[tt]] = 1;
        
                for (tt = pHG->V[v2].iStart; tt < pHG->V[v2].iEnd; tt++)
                {
                    const long net = pHG->VtxAdjncy[tt];
                    const long netsize = pHG->N[net].iEnd - pHG->N[net].iStartP0;
                    
                    if (flags[net] == 1)
                    {
                        IntegerWeight += 1;
                        ScaledWeight += 1.0/(double)netsize;
                    }
                }
        
                for (tt = pHG->V[v1].iStart; tt < pHG->V[v1].iEnd; tt++) flags[pHG->VtxAdjncy[tt]] = 0;
            }
            else if (pC->Start[t + 1] != pC->Start[t] + 1) {
                fprintf(stderr, "CoarsenGraph(): Erroneous matching group size!\n");
                return FALSE;
            }
        }
        
        fprintf(stderr, "Total matching weight equals %ld, %f.\n", IntegerWeight, ScaledWeight);
        
        free(flags);
    }
#endif
    
    free(Matched);
    if (iv != NULL) free(iv);
  
    /* Reallocate the memory of pC->Start to the number of matches,
       possibly decreasing the size and hence saving memory space.
       This is done because C is returned by the function. */
    pC->Start = (long *) realloc(pC->Start, (pC->NrMatches+1) * sizeof(long));
    
    if (pC->Start == NULL) {
        fprintf(stderr, "CoarsenGraph(): Not enough memory for realloc of pC->Start!\n");
        return FALSE;
    }
  
    /*-------------------------------------------------------------------------
      ## Contraction stage ##
      -------------------------------------------------------------------------*/

    if (!CreateNewBiPartHyperGraph(pC->NrMatches, pHG->NrNets - NumFreeNets, pHG->NrPins, FALSE, 'P', pCHG)) {
        fprintf(stderr, "CoarsenGraph(): Unable to create new bipart hypergraph!\n");
        return FALSE;
    }
    
    pCHG->SplitDir = pHG->SplitDir;
    
    /* Initialise the working array iAdjncy */
    iAdjncy = (long *) malloc((pCHG->NrNets+1) * sizeof(long));
    
    if (iAdjncy == NULL) {
        fprintf(stderr, "CoarsenGraph(): Not enough memory for iAdjncy!\n");
        return FALSE;
    }
    
    for (t = 0; t < pCHG->NrNets; t++)
        iAdjncy[t] = DUMMY;
    
    /* Create new adjacency lists and discard all free nets. */
    nrpins = 0;
    for (k = 0; k < pC->NrMatches; k++) {
        /* Merge the vertices of group k into one vertex */

        pCHG->V[k].iStart = nrpins;
    
        for (t = pC->Start[k]; t < pC->Start[k+1]; t++) {
            /* Merge vertex v = pC->Match[t] from hypergraph HG
               into vertex k of the contracted hypergraph CHG */
            v = pC->Match[t];
            
#ifdef INFO
            /* Matched vertices should never be free! */
            if (pHG->V[v].Free) fprintf(stderr, "CoarsenGraph(): Matched vertex is free!\n");
#endif
            
            /* Add the weight of this vertex to the weight of k. */
            pCHG->V[k].vtxwgt += pHG->V[v].vtxwgt;
            
            /* Add the nets of this vertex to the vertex adjacency list of k,
               but only upon first encounter and only if the net is non-free. */
            for (tt = pHG->V[v].iStart; tt < pHG->V[v].iEnd; tt++) {
                if (!pHG->N[pHG->VtxAdjncy[tt]].Free) {
                    n = NetIndices[pHG->VtxAdjncy[tt]];
                    
                    if (iAdjncy[n] == DUMMY) { 
                        iAdjncy[n] = nrpins; /* not a dummy */
                        pCHG->VtxAdjncy[nrpins] = n;
                        nrpins++;
      
                        /* Increase the counter for the net adjacency of n.
                           iStartP1 of CHG is abused for this purpose */
                        pCHG->N[n].iStartP1++;
                    }
                }
            }
        }
  
        pCHG->V[k].iEnd = nrpins;
  
        /* Clear iAdjncy */
        for (tt = pCHG->V[k].iStart; tt < pCHG->V[k].iEnd; tt++) {
            n = pCHG->VtxAdjncy[tt];
            iAdjncy[n] = DUMMY;
        }
    }
  
    /* Reset the number of pins of the CHG to the correct value */
    pCHG->NrPins = nrpins;
    
    if (nrpins == 0) {
        /* It could be that the entire graph is free, in that case return an empty graph. */
        pCHG->NrVertices = pCHG->NrNets = pCHG->NrPins = pC->NrMatches = 0;
        free(iAdjncy);
        free(NetIndices);
        return TRUE;
    }

    /* Reallocate the memory to reflect this change */
    pCHG->VtxAdjncy = (long *) realloc(pCHG->VtxAdjncy, nrpins * sizeof(long));
    pCHG->NetAdjncy = (long *) realloc(pCHG->NetAdjncy, nrpins * sizeof(long));
    
    if (pCHG->VtxAdjncy == NULL || pCHG->NetAdjncy == NULL) {
        fprintf(stderr, "CoarsenGraph(): Not enough memory for reallocating pins!\n");
        return FALSE;
    }
  
    /* Set iStartP0, iStartP1, and iEnd for the nets of the CHG */
    nrpins = 0;
    for (t = 0; t < pCHG->NrNets; t++) {
        /* currently, CHG.N[t].istartP1 = number of adjacencies in net t */
        pCHG->N[t].iStartP0 = nrpins;
        pCHG->N[t].iEnd = nrpins;
        nrpins += pCHG->N[t].iStartP1;
    }

    for (t = 0; t < pCHG->NrNets; t++) {
        /* reset CHG.N[t].istartP1 to iStartP0, i.e.,
           the start of the adjacencies from net t */
        pCHG->N[t].iStartP1 = pCHG->N[t].iStartP0;
    }

    /* currently, iStartP0 = iStartP1 = iEnd for CHG */
  
    /* Set the net adjacency lists for the CHG and adjust iEnd */
    for (t = 0; t < pCHG->NrVertices; t++) {
        for (tt = pCHG->V[t].iStart; tt < pCHG->V[t].iEnd; tt++) {
            n = pCHG->VtxAdjncy[tt];
            pCHG->NetAdjncy[pCHG->N[n].iEnd++] = t;
        }
    }

    if (pCHG->NrNets > 0) if (pCHG->N[pCHG->NrNets - 1].iEnd != nrpins || pCHG->NrPins != nrpins) fprintf(stderr, "CoarsenGraph(): Invalid number of pins!\n");

    /* now, for all t: N[t].iEnd = end of adjacencies of net t */

    /* Copy net weights and directions to the contracted hypergraph */
    for (t = 0; t < pHG->NrNets; ++t) {
        if (!pHG->N[t].Free) {
            n = NetIndices[t];
            pCHG->N[n].netwgt = pHG->N[t].netwgt;
            pCHG->N[n].dir = pHG->N[t].dir;
        }
    }
    
    free(iAdjncy);
    free(NetIndices);
    
#ifdef INFO2
    printf("CoarsenGraph(): Coarsened graph from %ld to %ld vertices and from %ld to %ld nets in direction %d.\n", pHG->NrVertices, pCHG->NrVertices, pHG->NrNets, pCHG->NrNets, pHG->SplitDir);
#endif
 
    return TRUE;
} /* end CoarsenGraph */

int UncoarsenGraph(struct biparthypergraph *pCHG, struct contraction *pC, struct biparthypergraph *pHG, const long weightlo, const long weighthi, const struct opts *pOptions) {
    /* This function uncoarsens the given hypergraph pCHG to the hypergraph pHG,
       using the contraction information provided in pC.
       
       The partitioning information of pCHG is propagated to pHG and free vertices
       (if any) are re-assigned using a greedy algorithm to optimise balancing.
    */
    long t, tt, v;
    
    if (!pHG || !pCHG || !pC || !pOptions) {
        fprintf(stderr, "UncoarsenGraph(): Null arguments!\n");
        return FALSE;
    }
    
    if (pCHG->NrVertices != pC->NrMatches) {
        fprintf(stderr, "UncoarsenGraph(): Hypergraph and contraction are incompatible!\n");
        return FALSE;
    }
    
    /* Set all partitions to 0 to prevent problems with unassigned free vertices. */
    for (t = 0; t < pHG->NrVertices; ++t) pHG->V[t].partition = 0;
    
    /* First uncoarsen and copy partitioning information. */
    for (t = 0; t < pC->NrMatches; ++t) {
        for (tt = pC->Start[t]; tt < pC->Start[t + 1]; ++tt) {
            v = pC->Match[tt];
            pHG->V[v].partition = pCHG->OptPartVtx[t];
        }
    }
    
    /* Then recalculate part weights. */
    pHG->WeightP[0] = pHG->WeightP[1] = 0;
    
    if (pOptions->DiscardFreeNets == FreeNetYes) {
        long *Weights, *Indices, *Order;
        long NumFreeVertices = 0;
        long totwgt;
        float oldratio;
        
        /* Calculate part weights for non-free vertices. */
        for (t = 0; t < pHG->NrVertices; ++t) {
            if (pHG->V[t].Free) ++NumFreeVertices;
            else pHG->WeightP[pHG->V[t].partition] += pHG->V[t].vtxwgt;
        }
        
        if (NumFreeVertices > 0) {
        /* Assign free vertices to optimize load balancing using the playground-bully algorithm. */
        Weights = (long *)malloc(NumFreeVertices*sizeof(long));
        Indices = (long *)malloc(NumFreeVertices*sizeof(long));
        
        if (Weights == NULL || Indices == NULL) {
            fprintf(stderr, "UncoarsenGraph(): Not enough memory!\n");
            return FALSE;
        }
        
        totwgt = 0;
        tt = 0;
        
        for (t = 0; t < pHG->NrVertices; ++t) {
            if (pHG->V[t].Free) {
                totwgt += (Weights[tt] = pHG->V[t].vtxwgt);
                Indices[tt] = t;
                ++tt;
            }
        }

#ifdef INFO2
        printf("    Reassigning %ld free vertices with total weight %ld.\n", NumFreeVertices, totwgt);
#endif
        oldratio = (float)pHG->WeightP[0]/(float)(pHG->WeightP[0] + pHG->WeightP[1]);
        
        if ((Order = QSort(Weights, NumFreeVertices)) == NULL) {
            fprintf(stderr, "UncoarsenGraph(): Unable to quicksort!\n");
            return FALSE;
        }
        
        for (t = 0; t < NumFreeVertices; ++t) {
            int P;

#ifdef INFO2
            if (pHG->WeightP[0] + Weights[t] > weightlo && pHG->WeightP[1] + Weights[t] > weighthi) printf("Impossible free vertex reassignment!\n");
#endif
            
            if (pHG->WeightP[0] + Weights[t] > weightlo) P = 1;
            else if (pHG->WeightP[1] + Weights[t] > weighthi) P = 0;
            else if (pHG->WeightP[0]*weighthi < weightlo*pHG->WeightP[1]) P = 0;
            else if (pHG->WeightP[0]*weighthi > weightlo*pHG->WeightP[1]) P = 1;
            else P = Random1(0, 1);
            
            pHG->WeightP[pHG->V[Indices[Order[t]]].partition = P] += Weights[t];
        }

#ifdef INFO2        
        printf("For %ld free vertices we improved (desired %f) %f to %f.\n", NumFreeVertices, (float)weightlo/(float)(weightlo + weighthi), oldratio, (float)pHG->WeightP[0]/(float)(pHG->WeightP[0] + pHG->WeightP[1]));
#endif
        
        free(Order);
        free(Indices);
        free(Weights);
        }
    }
    else {
        /* Calculate part weights directly without reassigning free vertices. */
        for (t = 0; t < pHG->NrVertices; ++t) pHG->WeightP[pHG->V[t].partition] += pHG->V[t].vtxwgt;
    }
    
    return TRUE;
}

#ifdef USE_PATOH
int RunMLGraphPart_PaToH(struct biparthypergraph *pHG, long weightlo, long weighthi, const struct opts *pOptions) {
    /* This function uses the PaToH library to perform multilevel hypergraph partitioning. */
    /* We will use the conventions from the PaToH user manual. */
    PaToH_Parameters args;
    int _c, _n, *cwghts, *nwghts, *xpins, *pins, *partvec, cut = 0, *partweights;
    const int _nconst = 1, _k = 2; /* Number of constraints and parts to which we will split. */
    long t, t2, pinCount, weighttot, weightfree;
    float targetweights[2];
    
    if (!pHG || !pOptions) {
        fprintf(stderr, "RunMLGraphPart_PaToH(): Null arguments!\n");
        return FALSE;
    }
    
    /* Allocate data. */
    cwghts = (int *)malloc(pHG->NrVertices*sizeof(int));
    nwghts = (int *)malloc(pHG->NrNets*sizeof(int));
    xpins = (int *)malloc((pHG->NrNets + 1)*sizeof(int));
    pins = (int *)malloc(pHG->NrPins*sizeof(int));
    partvec = (int *)malloc(pHG->NrVertices*sizeof(int));
    partweights = (int *)malloc(_k*sizeof(int));
    
    if (cwghts == NULL || nwghts == NULL || xpins == NULL || pins == NULL || partvec == NULL || partweights == NULL) {
        fprintf(stderr, "RunMLGraphPart_PaToH(): Not enough memory!\n");
        return FALSE;
    }
    
    /* Set default parameters. */
    PaToH_Initialize_Parameters(&args, PATOH_CUTPART, PATOH_SUGPARAM_DEFAULT);
    
    /* Copy data from hypergraph. */
    _c = pHG->NrVertices;
    
    weighttot = 0;
    
    /* Set all vertex weights. */
    for (t = 0; t < pHG->NrVertices; ++t) {
        cwghts[t] = pHG->V[t].vtxwgt;
        weighttot += cwghts[t];
    }
    
    /* Set up all nets. Here we discard nets that are either free or have weight <= 0. */
    _n = 0;
    pinCount = 0;
    
    for (t = 0; t < pHG->NrNets; ++t) {
        if (pHG->N[t].netwgt > 0 && !pHG->N[t].Free) {
            xpins[_n] = pinCount;
            nwghts[_n] = pHG->N[t].netwgt;
            _n++;
            
            /* Add pins from this net. */
            for (t2 = pHG->N[t].iStartP0; t2 < pHG->N[t].iEnd; ++t2) {
                pins[pinCount++] = pHG->NetAdjncy[t2];
            }
        }
    }
    
    xpins[_n] = pinCount;
    
    /* Calculate args.final_imbal, args.init_imbal. */
    /* PaToH partitions into two parts with weights bounded by q*(1 + e)*Wtot and (1 - q)*(1 + e)*Wtot where e is the imbalance
       and q the desired ratios of the parts.
       We want the smaller part to have size <= weightlo and the larger part to have size <= weighthi, hence we solve
       q*(1 + e)*Wtot = weightlo, (1 - q)*(1 + e)*Wtot = weighthi,
       from which we obtain
       q = weightlo/(weightlo + weighthi), e = (weightlo + weighthi - weighttot)/weighttot.
    */
    if (pOptions->Seed > 0) args.seed = pOptions->Seed;

    args._k = _k;
    args.final_imbal = args.init_imbal = (float)(weightlo + weighthi - weighttot)/(float)weighttot;
    targetweights[0] = (float)(weightlo - 0.5)/(float)(weightlo + weighthi); /* round down */
    targetweights[1] = 1.0 - targetweights[0];
    args.MemMul_Pins = 4 + 2*args.MemMul_Pins;

    /* Allocate data and perform partitioning. */
    PaToH_Alloc(&args, _c, _n, _nconst, cwghts, nwghts, xpins, pins);
    PaToH_Part(&args, _c, _n, _nconst, 0, cwghts, nwghts, xpins, pins, targetweights, partvec, partweights, &cut);

    /* Save partitioning data. */
    /* Assign all free vertices to part 0. */
    /* TODO: Should there be unassigned vertices at all? */
    weightfree = 0;
    
    for (t = 0; t < pHG->NrVertices; ++t) {
        if (partvec[t] < 0) {
            /* This vector has not been assigned yet. */
            weightfree += cwghts[t];
            pHG->OptPartVtx[t] = 0;
        }
        else {
            pHG->OptPartVtx[t] = partvec[t];
        }
    }
    
    if (weightfree > 0) fprintf(stderr, "RunMLGraph_PaToH(): Warning: Unassigned vertices (%ld)!\n", weightfree);
    
    pHG->WeightP[0] = partweights[0] + weightfree;
    pHG->WeightP[1] = partweights[1];
    pHG->OptComm = cut; /* TODO: Is this value still correct if there are free vertices? */
    
    /* Ensure that part 0 is smallest. */
    if (pHG->WeightP[0] > pHG->WeightP[1]) {
        pHG->WeightP[1] = partweights[0] + weightfree;
        pHG->WeightP[0] = partweights[1];
        
        for (t = 0; t < pHG->NrVertices; ++t) pHG->OptPartVtx[t] = 1 - pHG->OptPartVtx[t];
    }
    
    /* Free data. */
    PaToH_Free();
    
    free(cwghts);
    free(nwghts);
    free(xpins);
    free(pins);
    free(partvec);
    free(partweights);
    
    return TRUE;
}
#endif

int RunMLGraphPart(struct biparthypergraph *pHG, long weightlo, long weighthi, const struct opts *pOptions) {
    /* This function runs the MultiLevel hypergraph partitioning
       algorithm. It first performs a number of contractions
       to coarsen the hypergraph HG, then assigns each vertex either
       to part 0 or part 1, and finally uncoarsens the hypergraph,
       copying and refining the partition information.

       Input: weightlo = smallest upper bound for part weight, belongs to part 0
              weighthi = largest upper bound for part weight, belongs to part 1

       Output: HG.OptPartVtx[v] = 0  or 1, for all vertices v of HG,
               meaning that they belong either to part 0 or to part 1.
               HG.WeightP[0] and  HG.WeightP[1], the weights of parts
               0 and 1, are set accordingly.

    */

    long v;
    long nc; /* number of contractions */
    struct biparthypergraph *CHG; /* stores the sequence of increasingly
                                      coarser hypergraphs */
    struct contraction *C = NULL; /* stores the contraction information */
    double ratio;
    
    if (!pHG || !pOptions) {
        fprintf(stderr, "RunMLGraphPart(): Null arguments!\n");
        return FALSE;
    }
    
    /* Check whether or not we are goign to use PaToH instead of Mondriaan. */
    if (pOptions->Partitioner == PartPaToH)
    {
#ifdef USE_PATOH
        return RunMLGraphPart_PaToH(pHG, weightlo, weighthi, pOptions);
#else
        fprintf(stderr, "RunMLGraphPart(): The use of PaToH is requested, but this version of Mondriaan is compiled without PaToH support!\n");
        return FALSE;
#endif
    }
    
    /* Coarsening phase */
    CHG = (struct biparthypergraph *)malloc(sizeof(struct biparthypergraph)*(pOptions->Coarsening_MaxCoarsenings + 1));
    C = (struct contraction *)malloc(sizeof(struct contraction)*pOptions->Coarsening_MaxCoarsenings);
    
    if (CHG == NULL || C == NULL) {
        fprintf(stderr, "RunMLGraphPart(): Not enough memory for coarsening!\n");
        return FALSE;
    }
  
    nc = 0;  
    CHG[0] = *pHG;
  
    while (CHG[nc].NrVertices > pOptions->Coarsening_NrVertices && nc < pOptions->Coarsening_MaxCoarsenings) {
#ifdef INFO2
        printf("Coarsen CHG[%ld] (%ld vertices) \n", nc, CHG[nc].NrVertices);
#endif
        
        /* Coarsen the hypergraph CHG[nc] into CHG[nc+1]
           and register the coarsening information in C[nc]. */
        if (!CoarsenGraph(&CHG[nc], &CHG[nc+1], &C[nc], nc, pOptions)) {
            fprintf(stderr, "RunMLGraphPart(): Unable to coarsen hypergraph!\n");
            return FALSE;
        }
  
#ifdef INFO2
        printf("CHG[%ld] (%ld vertices) obtained\n", nc+1, CHG[nc+1].NrVertices);
#endif
  
        /* Stop coarsening if the number of matched vertices becomes
           relatively small or zero, but perform at least 2 coarsenings */
        ratio = (double) (CHG[nc].NrVertices - CHG[nc+1].NrVertices) / (double) CHG[nc].NrVertices;
        nc++;
        
        if (nc > 1 && (ratio < pOptions->Coarsening_StopRatio ||
                        CHG[nc-1].NrVertices == CHG[nc].NrVertices))
            break;
    }
    
    /* Initial partitioning */
#ifdef INFO2
    printf("Compute initial partition of CHG[%ld] \n", nc);
#endif

    if (!RunHKLFM(&CHG[nc], weightlo, weighthi, FALSE, pOptions)) {
        fprintf(stderr, "RunMLGraphPart(): Unable to run HKLFM!\n");
        return FALSE;
    }
  
    while (nc > 0) {
#ifdef INFO2 
        printf("Copy partition info from CHG[%ld] to CHG[%ld]\n", nc, nc-1);
#endif
        /* Uncoarsen graph. */
        if (!UncoarsenGraph(&CHG[nc], &C[nc - 1], &CHG[nc - 1], weightlo, weighthi, pOptions)) {
            fprintf(stderr, "RunMLGraphPart(): Unable to uncoarsen graph!\n");
            return FALSE;
        }
        
        /* Delete C[nc-1] */
        free(C[nc-1].Start);
        free(C[nc-1].Match);
        
        /* Delete CHG[nc] */
        DeleteBiPartHyperGraph(&CHG[nc]);

#ifdef INFO2
       if (CHG[nc-1].WeightP[0] > weightlo || CHG[nc-1].WeightP[1] > weighthi) printf("Imbalanced partitioning (pre-refinement): WeightP[0] = %ld <= %ld, WeightP[1] = %ld <= %ld!\n", CHG[nc-1].WeightP[0], weightlo, CHG[nc-1].WeightP[1], weighthi);
       else printf("WeightP[0] = %ld WeightP[1] = %ld (pre)\n", CHG[nc-1].WeightP[0], CHG[nc-1].WeightP[1]);
#endif

#ifdef INFO2 
        printf("Refine CHG[%ld]\n", nc-1);
#endif
 
        /* Refine the copied partition by running
           a simplified HKLFM algorithm */
        if (!RunHKLFM(&CHG[nc-1], weightlo, weighthi, TRUE, pOptions)) {
            fprintf(stderr, "RunMLGraphPart(): Unable to run HKLFM!\n");
            return FALSE;
        }

        /* Calculate corresponding weights */
        CHG[nc-1].WeightP[0] = 0;
        CHG[nc-1].WeightP[1] = 0;
        for (v = 0; v < CHG[nc-1].NrVertices; v++)
            CHG[nc-1].WeightP[CHG[nc-1].V[v].partition] += CHG[nc-1].V[v].vtxwgt;

#ifdef INFO2
       if (CHG[nc-1].WeightP[0] > weightlo || CHG[nc-1].WeightP[1] > weighthi) printf("Imbalanced partitioning (post-refinement): WeightP[0] = %ld <= %ld, WeightP[1] = %ld <= %ld!\n", CHG[nc-1].WeightP[0], weightlo, CHG[nc-1].WeightP[1], weighthi);
       else printf("WeightP[0] = %ld WeightP[1] = %ld (post)\n", CHG[nc-1].WeightP[0], CHG[nc-1].WeightP[1]);
#endif

        nc--;
    }
  
    *pHG = CHG[0];
  
#ifdef INFO2
    printf("Communication: %ld \n", pHG->OptComm);
#endif
    
    if (C != NULL)
        free(C);
    if (CHG != NULL)
        free(CHG);
    
    return TRUE;
} /* end RunMLGraphPart */

