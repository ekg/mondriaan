#ifndef __DistributeVecOpt2_h__
#define __DistributeVecOpt2_h__

#include "Options.h"
#include "Matalloc.h"
#include "DistributeVecLib.h"

/* Function declarations for DistributeVecOpt2.c */
int DistributeVecOpt2(const struct sparsematrix *pM, long int *X, int dir);
int FirstNzInRow (long **G, int P, int q, int start);

#endif /* __DistributeVecOpt2_h__ */
