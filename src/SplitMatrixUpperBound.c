#include "SplitMatrixUpperBound.h"

struct PartWeight {
    long part;
    long weight;
};

/**
 * Comparison functions.
 */
int comparePartWeightsASC (const void *a, const void *b) {
    long diff = ((struct PartWeight*)a)->weight - ((struct PartWeight*)b)->weight;
    if(diff == 0)
        return 0;
    return (diff < 0) ? -1 : 1;
} /* end comparePartWeightsASC */

int comparePartWeightsDESC (const void *a, const void *b) {
    long diff = ((struct PartWeight*)a)->weight - ((struct PartWeight*)b)->weight;
    if(diff == 0)
        return 0;
    return (diff > 0) ? -1 : 1;
} /* end comparePartWeightsDESC */


int SplitMatrixUpperBound(struct sparsematrix *pT, int P, const struct opts *pOptions) {
    /* This function splits the sparse matrix T into P (nearly-)equal parts.
       The split is performed by distributing entire columns (or rows) over the parts,
       possibly cutting some columns (or rows) to maintain feasibility. Nonzeros may be moved.
       This function does not minimize communication volume, but does keep it below the
       upper bound (min(m,n)+1)(P-1).
       
       Input: T sparse matrix,
              P number of parts to split into,
              options the options.
       
       Output: T sparse matrix
               The nonzeros of part i are in positions pT->Pstart[i], pT->Pstart[i+1]-1,
               where 0<=i<P.
       
       NB: This function works only if the weights are based on the number of nonzeros.
       It does not work when:
        - the matrix is stored as a symmetric matrix,
        - dummy nonzeros with weight 0 are present,
        - or weights are based on matrix columns.
    */
       
    long s, t,          /* Nonzero iterators */
         j, n,          /* Column (row) index and number of columns (rows) */
         *nzInd,        /* Pointer to matrix column (row) indices */
         *w,            /* w[j] = number of nonzeros in column (row) j */
         *lastCol,      /* Last column added to a part */
         *partWeights,  /* Weight of a part */
         *partPointer,  /* Nonzero iterator per part */
         delta, end, targetWeight;
    
    int *p,      /* p[j] = part that column (row) j is assigned to */
        k,       /* General iterator over all parts */
        l,       /* First dimension iterator in ColWghtDist[] */
        q,       /* General processor number; also second dimension iterator in ColWghtDist[] */
        numMinHeap, numMaxHeap;
    
    struct PartWeight maxWeight, minWeight, *MaxHeapItems, *MinHeapItems;
    struct heap MinHeap, MaxHeap;

    if (!pT || !pOptions) {
        fprintf(stderr, "SplitMatrixUpperBound(): Null arguments!\n");
        return FALSE;
    }
    
    if ((pT->MMTypeCode[3]=='S' || pT->MMTypeCode[3]=='K' || pT->MMTypeCode[3]=='H') &&
             pOptions->SymmetricMatrix_UseSingleEntry == SingleEntYes) {
        fprintf(stderr, "SplitMatrixUpperBound(): Cannot be used with symmetric matrices!\n");
        return FALSE;
    }
    if (pT->NrDummies > 0) {
        fprintf(stderr, "SplitMatrixUpperBound(): Cannot be used with dummy nonzeros!\n");
        return FALSE;
    }
    if (pT->MMTypeCode[0] == 'W' && pT->NrColWeights > 0) {
        fprintf(stderr, "SplitMatrixUpperBound(): Cannot be used with column-weighted matrix!\n");
        return FALSE;
    }
    
    if(pT->NrProcs < P) {
        pT->NrProcs = P;
        pT->Pstart = (long *)realloc(pT->Pstart, (pT->NrProcs+1)*sizeof(long));
        if(pT->Pstart == NULL) {
            fprintf( stderr, "SplitMatrixUpperBound(): Not enough memory.\n" );
            return FALSE;
        }
    }
    
    /****
     * Determine direction
     ****/
    if(pT->m < pT->n) {
        n = pT->m;
        nzInd = pT->i;
    }
    else {
        n = pT->n;
        nzInd = pT->j;
    }
    
    /* Compute number of nonzeros in a column (row) */
    w = (long *)calloc(n, sizeof(long));
    if(w == NULL) {
        fprintf( stderr, "SplitMatrixUpperBound(): Not enough memory.\n" );
        return FALSE;
    }
    for(t=0; t<pT->NrNzElts; ++t) {
        ++w[nzInd[t]];
    }
    
    /****
     * Create min heap of part weights
     ****/
    MinHeapItems = (struct PartWeight*)malloc(P*sizeof(struct PartWeight));
    if(MinHeapItems == NULL) {
        fprintf( stderr, "SplitMatrixUpperBound(): Not enough memory.\n" );
        return FALSE;
    }
    
    HeapInit(&MinHeap, sizeof(struct PartWeight), comparePartWeightsDESC);
    
    for(k=0; k<P; ++k) {
        MinHeapItems[k].part = k;
        MinHeapItems[k].weight = 0;
    }
    Heapify(&MinHeap, MinHeapItems, P);
    
    /****
     * Assign columns (rows) to parts
     ****/
    p = (int *)malloc(n * sizeof(int));
    lastCol = (long *)malloc(P * sizeof(long));
    if(p == NULL || lastCol == NULL) {
        fprintf( stderr, "SplitMatrixUpperBound(): Not enough memory.\n" );
        return FALSE;
    }
    
    for(k=0; k<P; ++k) {
        lastCol[k] = -1;
    }
    
    for(j=0; j<n; ++j) {
        /* Assign column (row) j to part minWeight.part with currently least weight. */
        HeapPeek(&MinHeap, &minWeight);
        minWeight.weight += w[j];
        p[j] = minWeight.part;
        lastCol[minWeight.part] = j;
        HeapReplace(&MinHeap, &minWeight, NULL);
    }
    
    /****
     * Compute surpluses and shortages per part and store them in heaps, to be able to determine where we should cut columns (rows).
     * Also register weights of the columns (rows) that will be cut, to be able to keep track of how we want to cut them.
     * We reuse the variable MinHeap, which now will contain shortages, while MaxHeap will contain surpluses.
     ****/
    
    HeapInit(&MaxHeap, sizeof(struct PartWeight), comparePartWeightsASC);
    MaxHeapItems = (struct PartWeight*)malloc(P*sizeof(struct PartWeight));
    
    long **ColWghtDist = (long **)malloc(P * sizeof(long*));
    long *_ColWghtDist = (long *)calloc(P * P, sizeof(long));
    if(MaxHeapItems == NULL || ColWghtDist == NULL || _ColWghtDist == NULL) {
        fprintf( stderr, "SplitMatrixUpperBound(): Not enough memory.\n" );
        return FALSE;
    }
    
    for(l=0; l<P; ++l) {
        ColWghtDist[l] = &(_ColWghtDist[l*P]);
    }
    
    targetWeight = ceil(pT->NrNzElts / (double)P);
    numMinHeap = 0;
    numMaxHeap = 0;
    l = 0;
    for(k=0; k<P; ++k) {
        /* Here we use the fact that the Heap data structure
         * just permutes all entries, hence we can just iterate over
         * these entries as if it were an array (which it 
         * (MinHeapItems) actually still is). We could also iteratively
         * call HeapPop(), but this would be O(P log(P)) while just
         * iterating over the array is O(P).
         */
        if(MinHeapItems[k].weight > targetWeight) {
            /* This part has a surplus */
            MaxHeapItems[numMaxHeap].part = MinHeapItems[k].part;
            MaxHeapItems[numMaxHeap].weight = MinHeapItems[k].weight - targetWeight;
            ++numMaxHeap;
            
            /* Mark last column as cut and link it to ColWghtDist[l] */
            q = MinHeapItems[k].part;
            j = lastCol[q];
            p[j] = -l-1;
            ColWghtDist[l][q] = w[j];
            
            ++l;
        }
        else if(MinHeapItems[k].weight < targetWeight) {
            /* This part has a shortage */
            MinHeapItems[numMinHeap].part = MinHeapItems[k].part;
            MinHeapItems[numMinHeap].weight = MinHeapItems[k].weight - targetWeight;
            ++numMinHeap;
        }
    }
    MinHeapItems = (struct PartWeight*)realloc(MinHeapItems, numMinHeap*sizeof(struct PartWeight));
    MaxHeapItems = (struct PartWeight*)realloc(MaxHeapItems, numMaxHeap*sizeof(struct PartWeight));
    Heapify(&MinHeap, MinHeapItems, numMinHeap);
    Heapify(&MaxHeap, MaxHeapItems, numMaxHeap);
    
    free(w); w = NULL;
    
    /****
     * Cut columns (rows) where necessary, until perfect balance is achieved.
     * By cutting a column, we can move weight from a surplus part to a shortage part.
     ****/
    while(MaxHeap.numItems > 0) {
        if(MinHeap.numItems < 1) {
            fprintf( stderr, "SplitMatrixUpperBound(): Corrupt data.\n" );
            return FALSE;
        }
        HeapPeek(&MinHeap, &minWeight); /* minWeight.weight < 0 */
        HeapPeek(&MaxHeap, &maxWeight); /* maxWeight.weight > 0 */
        
        /* Determine the weight (delta) to move */
        if(maxWeight.weight < -minWeight.weight) {
            delta = maxWeight.weight;
            
            /* Remove maxWeight, reduce minWeight */
            HeapPop(&MaxHeap, NULL);
            minWeight.weight += maxWeight.weight;
            HeapReplace(&MinHeap, &minWeight, NULL);
        }
        else if(maxWeight.weight > -minWeight.weight) {
            delta = -minWeight.weight;
            
            /* Remove minWeight, reduce maxWeight */
            HeapPop(&MinHeap, NULL);
            maxWeight.weight += minWeight.weight;
            HeapReplace(&MaxHeap, &maxWeight, NULL);
        }
        else { /* max == min */
            delta = maxWeight.weight;
            
            HeapPop(&MaxHeap, NULL);
            HeapPop(&MinHeap, NULL);
        }
        
        /* Keep track of the cut */
        l = -p[lastCol[maxWeight.part]]-1;
        ColWghtDist[l][maxWeight.part] -= delta;
        ColWghtDist[l][minWeight.part] = delta;
    }
    
    free(lastCol); lastCol = NULL;
    HeapDestroy(&MaxHeap); MaxHeapItems = NULL; /* Also free()s MaxHeapItems */
    
    /* cutColPrtPtr[l] := the first (lowest) index q for which ColWghtDist[l][q] is nonzero. */
    int *cutColPrtPtr = (int *)malloc(numMaxHeap * sizeof(int));
    if(cutColPrtPtr == NULL) {
        fprintf( stderr, "SplitMatrixUpperBound(): Not enough memory.\n" );
        return FALSE;
    }
    for(l=0; l<numMaxHeap; ++l) {
        cutColPrtPtr[l] = 0;
        
        for(q=0; q<P; ++q) {
            if(ColWghtDist[l][q] > 0) {
                break;
            }
            ++cutColPrtPtr[l];
        }
    }
    
    /* At this point:
     * p[j] =  q>0    if column j is completely assigned to part q
     *         q<0    if column j cut. The distribution data is available in ColWghtDist[-q-1]
     */
    
    /****
     * Compute total part weights and build Pstart variable
     ****/
    partWeights = (long *)malloc(P * sizeof(long));
    partPointer = (long *)malloc(P * sizeof(long));
    if(partWeights == NULL || partPointer == NULL) {
        fprintf( stderr, "SplitMatrixUpperBound(): Not enough memory.\n" );
        return FALSE;
    }
    
    for(k=0; k<P; ++k) {
        partWeights[k] = targetWeight;
    }
    for(k=0; k<MinHeap.numItems; ++k) {
        /* Here again we use the fact that MinHeapItems is just
         * a valid (permuted) array to iterate over
         */
        partWeights[MinHeapItems[k].part] += MinHeapItems[k].weight;
    }
    
    HeapDestroy(&MinHeap); MinHeapItems = NULL; /* Also free()s MinHeapItems */
    
    /* Build Pstart */
    pT->Pstart[k] = partPointer[0] = 0;
    for(k=1; k<P; ++k) {
        pT->Pstart[k] = partPointer[k] = partPointer[k-1]+partWeights[k-1];
    }
    pT->Pstart[P] = partPointer[P-1]+partWeights[P-1];
    
    /* Sanity check */
    if(pT->Pstart[P] != pT->NrNzElts) {
        fprintf( stderr, "SplitMatrixUpperBound(): Corrupt sum.\n" );
        return FALSE;
    }
    
    free(partWeights); partWeights = NULL;
    
    
    /****
     * Sort nonzeros by part number in order to apply the computed partitioning
     ****/
    for(k=0; k<P; ++k) {
        end = pT->Pstart[k+1];
        
        while(partPointer[k]<end) {
            /* We move the nonzero at position s to part q */
            s = partPointer[k];
            q = p[nzInd[s]];
            
            if(q < 0) {
                /* We have a cut column (row), so determine q from the cut data */
                l = -q-1;
                
                if(ColWghtDist[l][k] > 0) {
                    /* We may still add this nonzero to the current part.
                     * Doing this, we prevent 1 swap.
                     */
                    q = k;
                }
                else {
                    /* We should add this nonzero to an other part. */
                    q = cutColPrtPtr[l];
                    
                    if(ColWghtDist[l][q] == 0) {
                        fprintf( stderr, "SplitMatrixUpperBound(): Unexpected empty column part. %d %d\n", l, q );
                        return FALSE;
                    }
                }
                
                --ColWghtDist[l][q];
                
                /* Update cutColPrtPtr to point to a part we may still assign nonzeros to */
                if(q == cutColPrtPtr[l]) {
                    while(cutColPrtPtr[l] < P && ColWghtDist[l][cutColPrtPtr[l]] == 0) {
                        ++cutColPrtPtr[l];
                    }
                }
                
            }
            /* else, column (row) nzInd[s] is completely assigned to part q */
            
            /* If nonzero is already in the right part, we do not need to swap */
            if(k == q) {
                ++partPointer[k];
                continue;
            }
            
            /* ('Else',) The nonzero needs to be swapped to the target part */
            t = partPointer[q];
            
#ifdef INFO2
            if(t > pT->Pstart[q+1]) {
                fprintf( stderr, "SplitMatrixUpperBound(): Corrupt target pointer.\n" );
                return FALSE;
            }
#endif
            
            SwapLong(pT->i, s, t);
            SwapLong(pT->j, s, t);
            if (pT->MMTypeCode[2] != 'P')
                SwapDouble(pT->ReValue, s, t);
            if (pT->MMTypeCode[2] == 'C')
                SwapDouble(pT->ImValue, s, t);
            
            /* Update partPointer to point to the next nonzero not belonging to q */
            do {
                ++partPointer[q];
            }
            while(partPointer[q] < pT->Pstart[q+1] && p[nzInd[partPointer[q]]] == q);
        }
    }
    
#ifdef INFO2
    /* Check part integrity */
    for(k=0; k<P; ++k) {
        if(partPointer[k] != pT->Pstart[k+1]) {
            fprintf( stderr, "SplitMatrixUpperBound(): Corrupt result.\n" );
            return FALSE;
        }
    }
#endif
    
    free(p);
    free(cutColPrtPtr);
    free(ColWghtDist);
    free(_ColWghtDist);
    free(partPointer);

#ifdef INFO2
    if(!CheckUpperBoundSolution(pT)) {
        return FALSE;
    }
#endif
    
    return TRUE;
} /* end SplitMatrixUpperBound */


int CheckUpperBoundSolution(struct sparsematrix *pT) {
    /* This function checks whether the distributed matrix computed by
       SplitMatrixUpperBound() is valid. I.e., it checks whether
       the computed solution has a fairly distributed work load and
       communication volume not higher than the upper bound.
       
       Input: T sparse matrix
       Output: TRUE if the solution is valid, FALSE otherwise
    */
    
    int P = pT->NrProcs;
    
    /* Calculate weights */
    long Wmax = 0, weight;
    for(int q=0; q<P; ++q) {
        weight = pT->Pstart[q+1] - pT->Pstart[q];
        if(weight > Wmax) {
            Wmax = weight;
        }
    }
    
    if(Wmax > ceil(pT->NrNzElts/(double)P)) {
#ifdef INFO2
        fprintf( stderr, "CheckUpperBoundSolution(): Invalid imbalance result.\n" );
#endif
        return FALSE;
    }
    
    /* Calculate communication volume */
    long ComVol1, ComVol2, tmp;
    CalcCom(pT, NULL, (pT->m < pT->n)?ROW:COL, &ComVol1, &tmp, &tmp, &tmp, &tmp);
    CalcCom(pT, NULL, (pT->m < pT->n)?COL:ROW, &ComVol2, &tmp, &tmp, &tmp, &tmp);
    
    long n = (pT->m < pT->n)?pT->m:pT->n;
    if(ComVol1 > n*(P-1) || ComVol2 > P-1) {
#ifdef INFO2
        fprintf( stderr, "CheckUpperBoundSolution(): Invalid communication result.\n" );
#endif
        return FALSE;
    }
    
    return TRUE;
} /* end CheckUpperBoundSolution */
