#ifndef __SplitMatrixUpperBound_h__
#define __SplitMatrixUpperBound_h__

#include "Heap.h"
#include "DistributeVecLib.h"

int SplitMatrixUpperBound(struct sparsematrix *pT, int P, const struct opts *pOptions);
int CheckUpperBoundSolution(struct sparsematrix *pT);

#endif /* __SplitMatrixUpperBound_h__ */
