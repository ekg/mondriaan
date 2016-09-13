#ifndef __DistributeMat_h__
#define __DistributeMat_h__

#include "Graph.h"
#include "HKLFM.h"
#include "Options.h"
#include "Sort.h"
#include "SparseMatrix.h"
#include "Remembrance.h"

/* Function declarations for DistributeMat.c */
long ComputeWeight(const struct sparsematrix *pT, long lo, long hi, long *wnz, const struct opts *pOptions);
int DistributeMatrixMondriaan(struct sparsematrix *pT, int P, double eps, const struct opts *pOptions, int (*Callback)(int, int, const struct sparsematrix *));

#endif /* __DistributeMat_h__ */

