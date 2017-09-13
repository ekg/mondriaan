#include "Sort.h"
#include "SparseMatrix.h"

void ConstructDisconnectedMatrix(struct sparsematrix *pA, int symmetric, int dummies, int colWeights,
                                 long numCompMin, long numCompMax, long *pNumComponents,
                                 long **pComponent_m, long **pComponent_n, long **pComponent_weights,
                                 long ***p_i_to_I, long ***p_j_to_J);

void DestructDisconnectedMatrix(struct sparsematrix *pA, int symmetric, int dummies, int colWeights,
                                long **pComponent_m, long **pComponent_n, long **pComponent_weights,
                                long ***p_i_to_I, long ***p_j_to_J);
