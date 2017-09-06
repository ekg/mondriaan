#include "FreeNonzeros.h"

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
 * Improve load balance by moving free nonzeros between two specified processors.
 * This function works for general and symmetric matrices.
 * It works when dummy nonzeros are present: dummies will not be moved as they do not contribute weight.
 * It does not work when weights are based on column weights.
 * 
 * Input:
 *   pOptions    Options struct
 *   procs       abs(procs[i]) = Number of processors each part should still be distributed over
 *   p1, p2      Processors to condider
 * 
 * Input/output:
 *   pM          The matrix
 * 
 * Return value: FALSE if an error occurred, TRUE otherwise
 */

int ImproveFreeNonzeros(struct sparsematrix *pM, const struct opts *pOptions, const int *procs, int p1, int p2) {
	
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
		fprintf(stderr, "ImproveFreeNonzeros(): Not enough memory!");
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
} /* end ImproveFreeNonzeros */
