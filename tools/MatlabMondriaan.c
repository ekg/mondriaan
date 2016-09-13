/*
MatlabMondriaan.c

Created by Bas Fagginger Auer, based on the work done by Ken Stanley.
Some extentions by Albert-Jan N. Yzelman.
For copyright notifications, see the COPYING file one directory up.

This is a Matlab wrapper for the Mondriaan sparse matrix library which takes as input a (sparse) matrix from Matlab and outputs the same matrix, but with all non-zero entries replaced by the corresponding processor number to which the non-zero entry has been assigned by the Mondriaan algorithm.
It requires the library file MatlabHelper.c which contains all the methods for a proper integration between Mondriaan and Matlab.

Usage:  [I, s, p, q, r, c, rh, ch, B, u, v] = MatlabMondriaan(A, p, epsilon, maxiter, permute, symm)
  [I, s, p, q, r, c, rh, ch, B] = MatlabMondriaan(A, p, epsilon, maxiter, permute, symm)
  [I, s, p, q] = MatlabMondriaan(A, p, epsilon, maxiter, permute, symm)
  [I, s] = MatlabMondriaan(A, p, epsilon, maxiter, permute, symm)
  
  Distributes the non-zero entries of A among p processors with imbalance epsilon using the Mondriaan algorithm,
  for use when calculating u=Av in parallel. Vector distribution may be skipped if Mondriaan is used to find a 
  better structure for your application other than parallel sparse matrix--vector SpMV multiplication (such as, but
  not limited to: cache-oblivious sequential SpMV, reduced fill-in sparse LU, cluster identification, et cetera);
  skipping is as easy as calling MatlabMondriaan without the final two output elements (see second suggested usage).
  The calculation of the separator indices and hierarchy is always done behind the scenes when the permute option is set;
  calling Mondriaan using the third suggested usage thus does not speed up Mondriaan (but still is a little faster since
  less return variables have to be passed to Matlab).
  Note: the fourth usage suggestion will override the user-supplied permute value and set it to none.

  -I(i,j) = #processor to which component A(i,j) is assigned to.
  -s  contains the statistics of the partitioning: distribution time, imbalance, maximum communication, and communication volume.
  -p  gives the row-permutation vector corresponding to P such that A(p,:)=PA.
  -q  gives the column-permutation vector corresponding to Q such that A(:,q)=A, and thus A(p,q)=PAQ.

  -r  corresponds to the block structure visible in PAQ; it gives the row-wise separators between the different blocks, whereas
  -c  gives the column-wise separators between the different blocks. See the Mondriaan User's guide for more detail.
  -rh gives the relative hierarchy of the blocks; the i'th rowwise block has rh(i) as its rowwise parent block.
  -ch does the same columnwise. Note that these hierarchy arrays build binary trees, which need not be complete. Separator blocks
            correspond to internal nodes in that tree, whereas the leaf nodes correspond to the uncut partition blocks.
  -B=PAQ(=A, when no permutation was requested).

  -u(i) = # of processor to which component x(i) is assigned to.
  -v(i) = # of processor to which component y(i) is assigned to.
  

  The variable maxiter contains the maximum number of iterations the algorithm should perform before exiting, useful to produce animations.
  The variable permute contains an integer which sets the desired permutation: 0 for no permutation, 1 for reverse BBD, 2 for SBD, and 3 for BBD.
  The variable symm    contains an integer value either set to 0 or 1, indicating if the matrix passed is symmetric.
*/

#include "MatlabHelper.h"

/* This function provides the actual Matlab interface. */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	/* These are variables we need to be able to use Mondriaan. */
	struct sparsematrix *MondriaanMatrix;
	struct MondriaanStats Stats;
	long int *MondriaanU, *MondriaanV;
	int NumProcessors, Permutation, Symm, SplitStrategy, MaxIterations;
	double Imbalance;
	long i;
	
	/* Check that the input parameters are valid. */
	if (nrhs != 7) mexErrMsgTxt("We require seven input variables: the sparse matrix, the number of processors, the imbalance factor, the maximum number of iterations, the permutation type, a symmetry boolean and the split strategy!");
	else if (nlhs < 2 || nlhs > 11) mexErrMsgTxt("We require at least two and at most eleven output variables!");
	else if (!mxIsSparse(prhs[0])) mexErrMsgTxt("The matrix on which we work needs to be sparse!");
	
	/* Extract the number of processors and imbalance and verify their values. */
	NumProcessors = (int)ExtractDouble (prhs[1]);
	Imbalance = ExtractDouble (prhs[2]);
	MaxIterations = (int)ExtractDouble (prhs[3]);
	Permutation = (int)ExtractDouble (prhs[4]);
	Symm = (int)ExtractDouble (prhs[5]);
	SplitStrategy = (int)ExtractDouble(prhs[6]);

	/* If not able to return any permute-related output variables, set permutation to none. */
	if( nlhs == 2 )
		Permutation = 0;
	
	if (NumProcessors <= 0) mexErrMsgTxt ("The number of processors should be a positive integer!");
	if (Imbalance <= 0.0 || Imbalance >= 1.0) mexErrMsgTxt ("The imbalance should lie between 0 and 1!");
	
	/* Convert matrix to a Mondriaan friendly format. */
	MondriaanMatrix = ConvertMatlabToMondriaan(prhs[0]);
	
	if (MondriaanMatrix == NULL) mexErrMsgTxt("Unable to convert Matlab matrix to Mondriaan!");
	
	if( Symm ) {
		MondriaanMatrix->MMTypeCode[3] = 'S';
	}

	/* Run the Mondriaan algorithm. */
	if (nlhs >= 10) {
		MondriaanU = (long int *)malloc(MondriaanMatrix->m*sizeof(long int));
		MondriaanV = (long int *)malloc(MondriaanMatrix->n*sizeof(long int));
		if (MondriaanU == NULL || MondriaanV == NULL)
			mexErrMsgTxt("Not enough memory to allocate vectors!");
	} else {
		MondriaanU = NULL;
		MondriaanV = NULL;
	}
	
	if (!DoMondriaan(MondriaanMatrix, MondriaanU, MondriaanV, &Stats, NumProcessors, Imbalance, Permutation, Symm, SplitStrategy, MaxIterations)) mexErrMsgTxt("Unable to run the Mondriaan algorithm!");

	/* And we convert the Mondriaan matrix back to Matlab. */
	/* Return processor index matrix. */
	plhs[0] = ConvertMondriaanProcIndToMatlab(MondriaanMatrix);
	/* Return stats. */
	plhs[1] = ConvertStatsToVector(&Stats);
	
	if (plhs[0] == NULL || plhs[1] == NULL) mexErrMsgTxt("Unable to convert Mondriaan matrix back to Matlab!");

	/* Return permutation vectors if desired. */
	if (nlhs >= 3) {
		if (MondriaanMatrix->row_perm != NULL)
			plhs[2] = ConvertLongIndicesToVector(MondriaanMatrix->row_perm, MondriaanMatrix->m);
		else
			plhs[2] = IncrementalVector(MondriaanMatrix->m);
	}
	if (nlhs >= 4) {
		if (MondriaanMatrix->col_perm != NULL)
			plhs[3] = ConvertLongIndicesToVector(MondriaanMatrix->col_perm, MondriaanMatrix->n);
		else
			plhs[3] = IncrementalVector(MondriaanMatrix->n);
	}
	
	/* Return rowwise block indices if desired. */
	if (nlhs >= 5) {
		if (MondriaanMatrix->rowBoundaries != NULL) {
			if (remembrance_get( MondriaanMatrix->rowBoundaries )==NULL)
				mexErrMsgTxt("Error during read-out of row boundaries!\n");
			else
				plhs[4] = ConvertBoundaryIndicesToVector( MondriaanMatrix->rowBoundaries );
		} else
			plhs[4] = TwoVector(1,MondriaanMatrix->m+1);
	}
	
	/* Return columnwise block indices if desired. */
	if (nlhs >= 6) {
		if (MondriaanMatrix->colBoundaries != NULL) {
			if (remembrance_get( MondriaanMatrix->colBoundaries )==NULL)
				mexErrMsgTxt("Error during read-out of column boundaries!\n");
			else
				plhs[5] = ConvertBoundaryIndicesToVector( MondriaanMatrix->colBoundaries );
		} else
			if (MondriaanMatrix->rowBoundaries != NULL) /* Assume symmetric finegrain */
				plhs[5] = ConvertBoundaryIndicesToVector( MondriaanMatrix->rowBoundaries );
			else
				plhs[5] = TwoVector(1,MondriaanMatrix->n+1);
	}

	/* Return rowwise block hierarchy if desired. */
	if (nlhs >= 7) {
		if (MondriaanMatrix->rowBoundaries != NULL)
			plhs[6] = ConvertBoundaryHierarchyToVector( MondriaanMatrix->rowBoundaries );
		else
			plhs[6] = OneVector(0);
	}

	/* Return columnwise block hierarchy if desired. */
	if (nlhs >= 8) {
		if (MondriaanMatrix->colBoundaries != NULL)
			plhs[7] = ConvertBoundaryHierarchyToVector( MondriaanMatrix->colBoundaries );
		else
			if (MondriaanMatrix->rowBoundaries != NULL) /* Assume symmetric finegrain */
				plhs[7] = ConvertBoundaryHierarchyToVector( MondriaanMatrix->rowBoundaries );
			else
				plhs[7] = OneVector(0);
	}

	/* Return permuted matrix, if desired. */
	if (nlhs >= 9) {
		/* Commit matrix and vector distributions to permutation. */
		if (MondriaanMatrix->row_perm_inv != NULL && MondriaanMatrix->col_perm_inv!= NULL) {
			for( i=0; i<MondriaanMatrix->NrNzElts; i++ ) {
				MondriaanMatrix->i[ i ] = MondriaanMatrix->row_perm_inv[ MondriaanMatrix->i[ i ] ];
				MondriaanMatrix->j[ i ] = MondriaanMatrix->col_perm_inv[ MondriaanMatrix->j[ i ] ];
				if( MondriaanMatrix->row_perm_inv[ MondriaanMatrix->i[ i ] ] < 0 )
					mexErrMsgTxt( "Permutation takes row index to negative values!" );
				if( MondriaanMatrix->row_perm_inv[ MondriaanMatrix->i[ i ] ] > MondriaanMatrix->m )
					mexErrMsgTxt( "Permutation takes row index over maximum!" );
				if( MondriaanMatrix->col_perm_inv[ MondriaanMatrix->j[ i ] ] < 0 )
					mexErrMsgTxt( "Permutation takes column index to negative values!" );
				if( MondriaanMatrix->col_perm_inv[ MondriaanMatrix->j[ i ] ] > MondriaanMatrix->n )
					mexErrMsgTxt( "Permutation takes column index to negative values!" );
			}
		}

		/* Write to output */
		plhs[8] = ConvertMondriaanToMatlab(MondriaanMatrix);
		if (plhs[8] == NULL) mexErrMsgTxt( "Error during conversion of Mondriaan matrix to Matlab!");
	}

	/* Return vector distributions if desired. */
	if (nlhs >= 10) {
		plhs[9] = ConvertLongIndicesToVector(MondriaanU, MondriaanMatrix->m);
		free(MondriaanU);
		if (plhs[9] == NULL) mexErrMsgTxt("Unable to convert Mondriaan vectors to Matlab!");
	}
	if (nlhs >=11) {
		plhs[10] = ConvertLongIndicesToVector(MondriaanV, MondriaanMatrix->n);
		free(MondriaanV);
		if (plhs[10] == NULL) mexErrMsgTxt("Unable to convert Mondriaan vectors to Matlab!");
	}

	/* Free the Mondriaan matrix. */
	MMDeleteSparseMatrix(MondriaanMatrix);
}
