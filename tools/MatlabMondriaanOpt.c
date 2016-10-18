/*
MatlabMondriaanOpt.c

Created by Marco van Oort, based on the work done by Ken Stanley, Bas Fagginger Auer and Albert-Jan N. Yzelman.
For copyright notifications, see the COPYING file one directory up.

This is a Matlab wrapper for the Mondriaan sparse matrix library which takes as input a (sparse) matrix from Matlab and outputs the
same matrix, but with all non-zero entries replaced by the corresponding processor number to which the non-zero entry has been
assigned by the MondriaanOpt algorithm.
It requires the library file MatlabHelper.c which contains all the methods for a proper integration between Mondriaan and Matlab.

Usage:  [newVol] = MatlabMondriaan(A, epsilon, vol)
  
  Distributes the non-zero entries of A among 2 processors for use when calculating u=Av in parallel.
  The distribution resulting from the MondriaanOpt algorithm has lowest communication volume among all distributions that have imbalance at most epsilon.
  
  Parameters:
   A       : The matrix
   epsilon : Load imbalance parameter
   vol     : Initial upper bound for the communication volume

  Output:
   newVol  : The communication volume of the computed partitioning
  
  Resulting partitioning is written to MatlabMex.mtx-I2f and MatlabMex.mtx-2f.svg
  
*/

#include "MatlabHelper.h"
#include "../src/MondriaanOpt/solver.h"

int DoMondriaanOpt(struct sparsematrix *A, struct solution *sol, struct mat *a, double Imbalance, unsigned int Volume);


/* This function provides the actual Matlab interface. */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	/* These are variables we need to be able to use Mondriaan. */
	struct sparsematrix *MondriaanMatrix;
	unsigned int Volume;
	double Imbalance;
	
    /* Solution, options, matrix structure */
    struct solution sol;
    struct mat a;
	
	/* Check that the input parameters are valid. */
	if (nrhs != 3) mexErrMsgTxt("We require three input variables: the sparse matrix, the imbalance factor and the starting upper bound volume!");
	else if (nlhs < 0 || nlhs > 1) mexErrMsgTxt("We require zero or one output variables!");
	else if (!mxIsSparse(prhs[0])) mexErrMsgTxt("The matrix on which we work needs to be sparse!");
	
	/* Extract the number of processors and imbalance and verify their values. */
	Imbalance = ExtractDouble (prhs[1]);
	Volume = (int)ExtractDouble (prhs[2]);
	
	if (Imbalance <= 0.0 || Imbalance >= 1.0) mexErrMsgTxt ("The imbalance should lie between 0 and 1!");
	
	/* Convert matrix to a Mondriaan friendly format. */
	MondriaanMatrix = ConvertMatlabToMondriaan(prhs[0]);
	
	if (MondriaanMatrix == NULL) mexErrMsgTxt("Unable to convert Matlab matrix to Mondriaan!");

	/* Run the MondriaanOpt algorithm. */
	if (!DoMondriaanOpt(MondriaanMatrix, &sol, &a, Imbalance, Volume)) mexErrMsgTxt("Unable to run the MondriaanOpt algorithm!");
	
	if(nlhs > 0) {
		plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
		if (plhs[0] == NULL) mexErrMsgTxt("Unable to allocate Matlab vector!");
		
		double *wp = mxGetPr(plhs[0]);
		wp[0] = sol.maxvol;
	}
	
	/* Free the Mondriaan matrix. */
	MMDeleteSparseMatrix(MondriaanMatrix);
}

/* This function runs the MondriaanOpt algorithm on a given sparse matrix A. */
int DoMondriaanOpt(struct sparsematrix *A, struct solution *sol, struct mat *a, double Imbalance, unsigned int Volume) {
	struct options Options;
	
	if (A == NULL) return false;
	
	/* Remove double zeroes. */
	if (!SparseMatrixRemoveDuplicates(A)) return false;
	
	/* Convert sparse matrix to mat format */
	convertSparsematrixToMat(a, A);
	
	sprintf(Options.fn, "MatlabMex.mtx");
	Options.eps = Imbalance;
	Options.epsset = TRUE;
	Options.maxvol = Volume;
	
	/* Initialise solution */
	initsolution(a,sol,&Options);

	/* Start branch and bound algorithm */
	solve(a,sol);
	
	return true;
}

