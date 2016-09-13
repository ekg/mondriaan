#ifndef __DistributeVecGreedy_h__
#define __DistributeVecGreedy_h__

#include "Sort.h"
#include "Options.h"
#include "DistributeVecLib.h"

/* Function declarations for DistributeVecGreedy.c */
long DistributeVecGreedyImprove(const struct sparsematrix *pM, long int *X, int dir, const struct opts *pOptions);

#endif /* __DistributeVecGreedy_h__ */
