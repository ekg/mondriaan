/*
MatlabHelper.h

Created by Bas Fagginger Auer, based on the work done by Ken Stanley.
Some extentions by Albert-Jan N. Yzelman, Davide Taviani.
For copyright notifications, see the COPYING file one directory up.

This file is meant as a library that contains all the methods required for the integration between Mondriaan and Matlab.
It should be included in the actual wrapper, along with the function "void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])"
*/

#ifndef __Matlab_Helper_h
#define __Matlab_Helper_h

#include <stdlib.h>
#include <time.h>
#include <Mondriaan.h>

#ifndef USE_MATLAB
#error To be able to use the Matlab interface of Mondriaan, please make the appropriate changes in ../mondriaan.mk.
#endif

/* Matlab include files. */
#include "mex.h"
#include "matrix.h"

/* For compatibility with older versions of Matlab, we need to define mwSize and mwIndex ourselves. */
#ifndef mwSize
#define mwSize size_t
#endif
#ifndef mwIndex
#define mwIndex size_t
#endif

/* This structure hold the partitioning statistics. */
struct MondriaanStats {
	/* Maximum and minimum number of nonzeros assigned to a single processor. */
	long MaxNz, MinNz;
	double Epsilon;
	long MaxComU, MaxComV;
	long ComVolU, ComVolV;
	double Duration;
};

/* This structure will be used to translate a Mondriaan sparse matrix back to the Matlab format. */
struct MondriaanTriplet {
	long Row, Column;
	double ReValue, ImValue;
};

/* This function runs the Mondriaan algorithm on a given sparse matrix A. */
int DoMondriaan(struct sparsematrix *A, long int *u, long int *v, struct MondriaanStats *s, int NumProcessors, double Imbalance, int Permutation, int Symmetric, int SplitStrategy, int MaxIterations);

/* This function tries to extract a double precision value from the given variable. */
double ExtractDouble(const mxArray *Var);

/* This function converts a Matlab sparse matrix to a Mondriaan sparse matrix. */
struct sparsematrix *ConvertMatlabToMondriaan(const mxArray *A);

/* This function converts Mondriaan statistics to a matlab vector. */
mxArray *ConvertStatsToVector(const struct MondriaanStats *s);

/* This function returns the matlab vector [value]. */
mxArray *OneVector(const double value);

/* This function returns the matlab vector [value1;value2]. */
mxArray *TwoVector(const double value1, const double value2);

/* Thus function returns the matlab vector [1 2... Size]. */
mxArray *IncrementalVector(const int Size);

/* This function converts a vector assignment to a vector of appropriate size. */
mxArray *ConvertBoundaryIndicesToVector(const struct remembrance *m);

/* This function converts a vector assignment to a vector of appropriate size. */
mxArray *ConvertBoundaryHierarchyToVector(const struct remembrance *m);

/* This function converts a vector assignment to a vector of appropriate size. */
mxArray *ConvertLongIndicesToVector(const long *v, const size_t Size);

/* This function converts a vector assignment to a vector of appropriate size. */
mxArray *ConvertIndicesToVector(const int *v, const size_t Size);

/* This function converts a Mondriaan sparse matrix to a Matlab sparse matrix, Triplets are already created. */
mxArray *ConvertTripletsToMatlab(struct MondriaanTriplet *Triplets, const struct sparsematrix *B);

/* This function converts a Mondriaan sparse matrix to a Matlab sparse matrix. */
mxArray *ConvertMondriaanToMatlab(const struct sparsematrix *B);
	
/* This function converts a Mondriaan sparse matrix to a Matlab sparse matrix with processor indices as nonzero values. */
mxArray *ConvertMondriaanProcIndToMatlab(const struct sparsematrix *B);

#endif

