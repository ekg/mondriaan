#ifndef __FreeNonzeros_h__
#define __FreeNonzeros_h__

#include "Sort.h"
#include "SparseMatrix.h"
#include "DistributeMat.h"

struct freeIndex {
	long t;
	long numProcs;
};

void SwapNonzero(struct sparsematrix *pM, long s, long t);

int ImproveFreeNonzerosGlobal(struct sparsematrix *pM, const struct opts *pOptions, int P, const int *procs);

int ImproveFreeNonzerosLocal(struct sparsematrix *pM, const struct opts *pOptions, int P, const int *procs, int p1, int p2);

#endif /* __FreeNonzeros_h__ */
