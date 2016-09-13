#ifndef __Cartesian_h__
#define __Cartesian_h__

#include "DistributeVecLib.h"
#include "SparseMatrix.h"

/* Function declarations for Cartesian.c */
int MMWriteCartesianSubmatrices(const struct sparsematrix *pM, FILE *fp);
int CRS2CCS(long m, long n, long nz, long *start, long *Index);

#endif /* __Cartesian_h__ */
