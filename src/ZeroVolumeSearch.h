#ifndef __ZeroVolumeSearch_h__
#define __ZeroVolumeSearch_h__

#include "Sort.h"
#include "SparseMatrix.h"

/* Function declarations for ZeroVolumeSearch.c */
int ZeroVolumeSearch(struct sparsematrix *pM, long weightlo, long weighthi, long *mid, const struct opts *pOptions);
int SubsetSum(long * const weights, const long N, const long weightlo, const long weighthi, long **permutation, long **select);
int SubsetSumExp(long * const weights, const long N, const long weightlo, const long weighthi, long **permutation, long **select);
int KarmarkarKarp(const long * const weights, const long N, const long weightlo, const long weighthi, long **permutation, long **select);
int DetectConnectedComponents(struct sparsematrix *pM, long maxWeight, long *numComponents, long **componentWeights, long **rowAssignments, const struct opts *pOptions);

#endif /* __ZeroVolumeSearch_h__ */

