#ifndef __DistributeVecLocal_h__
#define __DistributeVecLocal_h__

#include "Sort.h"
#include "Matalloc.h"
#include "DistributeVecLib.h"

/* Function declarations for DistributeVecLocal.c */
long DistributeVecLocal(const struct sparsematrix *pM, long int *X, int dir);
void PrintVecLocalStatistics(int P, long *Ns, long *Nr, long *Ls, long *Lr,
                             long *Jstart0, long *J, long *owner);

#endif /* __DistributeVecLocal_h__ */
