#ifndef __DistributeVecOrig_h__
#define __DistributeVecOrig_h__

#include "Options.h"
#include "DistributeVecLib.h"

long DistributeVecOrig(const struct sparsematrix *pM, long int *X, int dir, const struct opts *pOptions);
int InitSums(long l, int P, long *procstart, int *procindex, long *Sums);

#endif /* __DistributeVecOrig_h__ */

