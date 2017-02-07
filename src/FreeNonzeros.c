#include "FreeNonzeros.h"

/**
 * Comparison function, taking two freeIndices a and b, and comparing the values
 * a->numProcs and b->numProcs in an ascending fashion.
 * Used in ImproveFreeNonzerosGlobal().
 * 
 * Return value:
 *      -1 if a->numProcs < b->numProcs
 *       0 if a->numProcs = b->numProcs
 *       1 if a->numProcs > b->numProcs
 */
int CompareFreeIndices (const void *a, const void *b) {
    long diff = ((struct freeIndex*)a)->numProcs - ((struct freeIndex*)b)->numProcs;
    if(diff == 0)
        return 0;
    return (diff < 0) ? -1 : 1;
} /* end CompareFreeIndices */

/**
 * Swap nonzeros s and t of matrix pM
 */
void SwapNonzero(struct sparsematrix *pM, long s, long t) {
	SwapLong(pM->i, s, t);
	SwapLong(pM->j, s, t);
	if(pM->MMTypeCode[2] != 'P')
		SwapDouble(pM->ReValue, s, t);
	if(pM->MMTypeCode[2] == 'C')
		SwapDouble(pM->ImValue, s, t);
} /* end SwapNonzero */

/**
 * Improve load balance by moving free nonzeros between all processors.
 * This function works for general and symmetric matrices.
 * It works when dummy nonzeros are present: dummies will not be moved as they do not contribute weight.
 * It does not work when weights are based on column weights.
 * 
 * Input:
 *   pOptions    Options struct
 *   P           Number of nonzeros the matrix is currently distributed over
 *   procs       abs(procs[i]) = Number of processors each part should still be distributed over
 * 
 * Input/output:
 *   pM          The matrix
 * 
 * Return value: FALSE if an error occurred, TRUE otherwise
 */
int ImproveFreeNonzerosGlobal(struct sparsematrix *pM, const struct opts *pOptions, int P, const int *procs) {
	
	long i, j, s, t, p, q, r;
	
	int symmetric = pM->m == pM->n &&
		(pM->MMTypeCode[3]=='S' || pM->MMTypeCode[3]=='K' || pM->MMTypeCode[3]=='H') &&
		pOptions->SymmetricMatrix_UseSingleEntry == SingleEntYes;
	int dummies = pM->m == pM->n && pM->NrDummies > 0;
	
	long *_rowsP = (long *)calloc(pM->m * P, sizeof(long));
	long *_colsP = (long *)calloc(pM->n * P, sizeof(long));
	
	long **rowsP = (long **)malloc(pM->m * sizeof(long *));
	long **colsP = (long **)malloc(pM->n * sizeof(long *));
	
	if(_rowsP == NULL || _colsP == NULL || rowsP == NULL || colsP == NULL) {
		fprintf(stderr, "ImproveFreeNonzerosGlobal(): Not enough memory!");
		if(_rowsP != NULL)
			free(_rowsP);
		if(_colsP != NULL)
			free(_colsP);
		if(rowsP != NULL)
			free(rowsP);
		if(colsP != NULL)
			free(colsP);
		return FALSE;
	}
	
	for(i=0; i<pM->m; ++i) {
		rowsP[i] = &(_rowsP[i*P]);
	}
	for(j=0; j<pM->n; ++j) {
		colsP[j] = &(_colsP[j*P]);
	}
	
	long *nnz = (long *)malloc(P * sizeof(long));
	
	/* Variables for freeIndices array */
	long freeIndicesSize = 16;
	long numFreeIndices = 0;
	struct freeIndex *freeIndices = (struct freeIndex *)malloc(freeIndicesSize * sizeof(struct freeIndex));
	
	if(nnz == NULL || freeIndices == NULL) {
		fprintf(stderr, "ImproveFreeNonzerosGlobal(): Not enough memory!");
		free(_rowsP);
		free(_colsP);
		free(rowsP);
		free(colsP);
		if(nnz != NULL)
			free(nnz);
		if(freeIndices != NULL)
			free(freeIndices);
		return FALSE;
	}
	
	/* Compute weight assigned to each of the processors */
	for(p=0; p<P; ++p) {
		nnz[p] = ComputeWeight(pM, pM->Pstart[p], pM->Pstart[p+1]-1, NULL, pOptions);
	}
	
	/* Determine which processors contain nonzeros in which rows/columns */
	for(p=0; p<P; ++p) {
		for(t=pM->Pstart[p]; t<pM->Pstart[p+1]; ++t) {
			
			rowsP[pM->i[t]][p]++;
			colsP[pM->j[t]][p]++;
			
			if(symmetric) {
				rowsP[pM->j[t]][p]++;
				colsP[pM->i[t]][p]++;
			}
		}
	}
	
	long numProcs; /* Free degree of a nonzero */
	
	/* Determine all free nonzeros */
	for(p=0; p<P; ++p) {
		for(t=pM->Pstart[p]; t<pM->Pstart[p+1]; ++t) {
			
			if(dummies && pM->i[t] == pM->j[t] && pM->dummy[pM->i[t]]) {
				continue;
			}
			
			numProcs = 1; /* In the loop we skip p==q, we add it here */
			for(q=0; q<P; ++q) {
				if(p == q)
					continue;
				
				if(rowsP[pM->i[t]][q] > 0 && colsP[pM->j[t]][q] > 0) {
					++numProcs;
				}
			}
			
			if(numProcs == 1) {
				continue;
			}
			
			/* Check size */
			if(numFreeIndices >= freeIndicesSize) {
				freeIndicesSize *= 2;
				struct freeIndex *newFreeIndices = (struct freeIndex *)realloc(freeIndices, freeIndicesSize*sizeof(struct freeIndex));
				if(newFreeIndices == NULL) {
					fprintf( stderr, "ImproveFreeNonzerosGlobal(): Not enough memory.\n" );
					free(_rowsP);
					free(_colsP);
					free(rowsP);
					free(colsP);
					free(nnz);
					free(freeIndices);
					return FALSE;
				}
				freeIndices = newFreeIndices;
			}
			
			/* Add free nonzero */
			freeIndices[numFreeIndices].t = t;
			freeIndices[numFreeIndices].numProcs = numProcs;
			++numFreeIndices;
			
		}
	}
	
	/* Sort the free nonzeros by their free degree, in an ascending fashion.
	 * This way, the nonzeros with the least freedom/flexibility are considered first,
	 * and the more flexibile nonzeros are considered later.
	 */
	qsort(freeIndices, numFreeIndices, sizeof(struct freeIndex), CompareFreeIndices);
	
	long minNrNzElts, minNrNzEltsP; /* The processor with the least weight */
	long nzWeight; /* Weight of the nonzero */
	
	for(s=0; s<numFreeIndices; ++s) {
		t = freeIndices[s].t;
		
		minNrNzElts = -1;
		minNrNzEltsP = -1;
		p = -1;
		
		for(q=0; q<P; ++q) {
			if(pM->Pstart[q] <= t && t < pM->Pstart[q+1]) {
				p = q;
			}
			if(rowsP[pM->i[t]][q] > 0 && colsP[pM->j[t]][q] > 0 && p != q) {
				if( nnz[q]/(double)abs(procs[q]) < minNrNzElts || minNrNzElts == -1) {
					minNrNzElts = ((double)nnz[q])/abs(procs[q]);
					minNrNzEltsP = q;
				}
			}
			
		}
		
		/* During the algorithm, the free degrees of nonzeros may decrease.
		 * Hence it may be possible we find no processor to move the nonzero to.
		 */
		if(minNrNzElts == -1) {
			continue;
		}
		
		q = minNrNzEltsP;
		nzWeight = (symmetric && pM->i[t] != pM->j[t])?2:1;
		
		if((nnz[q]+nzWeight)/(double)abs(procs[q]) > (nnz[p]-nzWeight)/(double)abs(procs[p])) {
			continue;
		}
		
		/* Move the nonzero */
		if(p < q) {
			SwapNonzero(pM, t, pM->Pstart[p+1]-1);
			
			for(r=p; r<q; ++r) {
				SwapNonzero(pM, pM->Pstart[r+1]-1, pM->Pstart[r+2]-1);
				--pM->Pstart[r+1];
			}
			
			t = pM->Pstart[q+1]-1;
		}
		else { /* p > q */
			SwapNonzero(pM, t, pM->Pstart[p]);
			
			for(r=p; r>q; --r) {
				SwapNonzero(pM, pM->Pstart[r], pM->Pstart[r-1]);
				++pM->Pstart[r];
			}
			
			t = pM->Pstart[q];
		}
		
		/* Update bookkeeping */
		rowsP[pM->i[t]][q]++;
		colsP[pM->j[t]][q]++;
		rowsP[pM->i[t]][p]--;
		colsP[pM->j[t]][p]--;
		if(symmetric) {
			rowsP[pM->j[t]][q]++;
			colsP[pM->i[t]][q]++;
			rowsP[pM->j[t]][p]--;
			colsP[pM->i[t]][p]--;
		}
		
		nnz[p] -= nzWeight;
		nnz[q] += nzWeight;
	}
	
	/* Finish */
	free(_rowsP);
	free(_colsP);
	free(rowsP);
	free(colsP);
	free(nnz);
	free(freeIndices);
	
	return TRUE;
} /* end ImproveFreeNonzerosGlobal */




/**
 * Improve load balance by moving free nonzeros between two specified processors.
 * This function works for general and symmetric matrices.
 * It works when dummy nonzeros are present: dummies will not be moved as they do not contribute weight.
 * It does not work when weights are based on column weights.
 * 
 * Input:
 *   pOptions    Options struct
 *   P           Number of nonzeros the matrix is currently distributed over
 *   procs       abs(procs[i]) = Number of processors each part should still be distributed over
 *   p1, p2      Processors to condider
 * 
 * Input/output:
 *   pM          The matrix
 * 
 * Return value: FALSE if an error occurred, TRUE otherwise
 */

int ImproveFreeNonzerosLocal(struct sparsematrix *pM, const struct opts *pOptions, int P, const int *procs, int p1, int p2) {
	
	long nnz1 = ComputeWeight(pM, pM->Pstart[p1], pM->Pstart[p1+1]-1, NULL, pOptions);
	long nnz2 = ComputeWeight(pM, pM->Pstart[p2], pM->Pstart[p2+1]-1, NULL, pOptions);
	long i, j, t, nzWeight;
	
	int symmetric = pM->m == pM->n &&
		(pM->MMTypeCode[3]=='S' || pM->MMTypeCode[3]=='K' || pM->MMTypeCode[3]=='H') &&
		pOptions->SymmetricMatrix_UseSingleEntry == SingleEntYes;
		
	int dummies = pM->m == pM->n && pM->NrDummies > 0;
	
	if(nnz1/abs(procs[p1]) > nnz2/abs(procs[p2])) {
		long tmp = p1;
		p1 = p2;
		p2 = tmp;
		tmp = nnz1;
		nnz1 = nnz2;
		nnz2 = tmp;
	}
	
	/* Now p2 is relatively larger than p1 */
	if((nnz1+1)/(double)abs(procs[p1]) >= (nnz2-1)/(double)abs(procs[p2])) {
		return TRUE;
	}
	
	char *rowsP1 = (char *)malloc(pM->m * sizeof(char));
	char *colsP1 = (char *)malloc(pM->n * sizeof(char));
	
	if(rowsP1 == NULL || colsP1 == NULL) {
		fprintf(stderr, "ImproveFreeNonzerosLocal(): Not enough memory!");
		if(rowsP1 != NULL)
			free(rowsP1);
		if(colsP1 != NULL)
			free(colsP1);
		return FALSE;
	}
	
	for(i=0; i<pM->m; ++i) {
		rowsP1[i] = 0;
	}
	for(j=0; j<pM->n; ++j) {
		colsP1[j] = 0;
	}
	
	/* In what columns/rows does P(1) have nonzeros? */
	for(t=pM->Pstart[p1]; t<pM->Pstart[p1+1]; ++t) {
		rowsP1[pM->i[t]] = 1;
		colsP1[pM->j[t]] = 1;
		
		if(symmetric) {
			rowsP1[pM->j[t]] = 1;
			colsP1[pM->i[t]] = 1;
		}
	}
	
	for(t=pM->Pstart[p2]; t<pM->Pstart[p2+1]; ++t) {
		if(rowsP1[pM->i[t]] == 1 && colsP1[pM->j[t]] == 1) {
			/* This nonzero is free. As p2 is relatively large, move it */
			
			if(dummies && pM->i[t] == pM->j[t] && pM->dummy[pM->i[t]]) {
				continue;
			}
			
			nzWeight = (symmetric && pM->i[t] != pM->j[t])?2:1;
			
			if((nnz1+nzWeight)/(double)abs(procs[p1]) > (nnz2-nzWeight)/(double)abs(procs[p2])) {
				free(rowsP1);
				free(colsP1);
				return TRUE;
			}
			
			/* Move the nonzero */
			if(p2 > p1) {
				SwapNonzero(pM, t, pM->Pstart[p2]);
				++pM->Pstart[p2];
			}
			else {
				SwapNonzero(pM, t, pM->Pstart[p2+1]-1);
				--pM->Pstart[p2+1];
			}
			
			/* Update bookkeeping */
			nnz1 += nzWeight;
			nnz2 -= nzWeight;
			
		}
		
	}
	
	/* Finish */
	free(rowsP1);
	free(colsP1);
	
	return TRUE;
} /* end ImproveFreeNonzerosLocal */
