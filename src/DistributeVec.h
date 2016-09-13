#ifndef __DistributeVec_h__
#define __DistributeVec_h__

#include "Sort.h"
#include "Options.h"
#include "SparseMatrix.h"
#include "DistributeVecGreedy.h"
#include "DistributeVecLocal.h"
#include "DistributeVecOpt2.h"
#include "DistributeVecOrig.h"
#include "DistributeVecLib.h"

long DistributeVec(const struct sparsematrix *pM, long int *X, int dir, const struct opts *pOptions);

#endif /* __DistributeVec_h__ */

