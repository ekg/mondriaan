/*
MatlabHelper.c

Created by Bas Fagginger Auer, based on the work done by Ken Stanley.
Some extentions by Albert-Jan N. Yzelman, Davide Taviani.
For copyright notifications, see the COPYING file one directory up.

This file is meant as a library that contains all the methods required for the integration between Mondriaan and Matlab.
It should be included in the actual wrapper, along with the function "void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])"

*/

#include <stdlib.h>
#include <time.h>
#include <Mondriaan.h>

#include "MatlabHelper.h"

int CallbackMaxIterations = -1;

/* This callback ensures that we stop after MaxIterations iterations. */
int MatlabCallback(int Iteration, int SplitPart, const struct sparsematrix *A) {
	if (CallbackMaxIterations <= 0) return true;
	
	return (Iteration < CallbackMaxIterations);
}

/* This function runs the Mondriaan algorithm on a given sparse matrix A. */
int DoMondriaan(struct sparsematrix *A, long int *u, long int *v, struct MondriaanStats *s, int NumProcessors, double Imbalance, int Permutation, int Symmetric, int SplitStrategy, int MaxIterations) {
	clock_t Clock;
	struct opts Options;
	int i;
	long l;
	
	if (A == NULL || s == NULL) return false;
	
	/* Read Mondriaan options. */
	SetDefaultOptions(&Options);
	if (!SetOptionsFromFile(&Options, "Mondriaan.defaults")) {
		FILE *File;
		
		fprintf(stderr, "main(): warning, cannot set options from 'Mondriaan.defaults', using default options and creating standard 'Mondriaan.defaults'!\n");
		File = fopen("Mondriaan.defaults", "w");
		if (File != NULL) {
			ExportDefaultOptions(File);
			fclose(File);
		} else {
			fprintf(stderr, "main(): Unable to create 'Mondriaan.defaults'!\n");
		}
	}
	
	CallbackMaxIterations = MaxIterations;

	/* Check whether the Split Strategy has been explicitly provided */
	if (SplitStrategy != -1){
		if (SplitStrategy == 0) Options.SplitStrategy = Alternate;
		else if (SplitStrategy == 1) Options.SplitStrategy = LocalBest;
		else if (SplitStrategy == 2) Options.SplitStrategy = Hybrid;
		else if (SplitStrategy == 3) Options.SplitStrategy = LocalRatio;
		else if (SplitStrategy == 4) Options.SplitStrategy = OneDimRow;
		else if (SplitStrategy == 5) Options.SplitStrategy = OneDimCol;
		else if (SplitStrategy == 6) Options.SplitStrategy = FineGrain;
		else if (SplitStrategy == 7) Options.SplitStrategy = SFineGrain;
		else if (SplitStrategy == 8) Options.SplitStrategy = MediumGrain;
		else fprintf(stderr, "Invalid value for the split strategy!\n");
	}

	/* Set desired permutation ordering. */
	if (Permutation == 1) Options.OrderPermutation = OrderPrefix;
	else if (Permutation == 2) Options.OrderPermutation = OrderInfix;
	else if (Permutation == 3) Options.OrderPermutation = OrderPostfix;
	else if (Permutation == 0) Options.OrderPermutation = OrderNone;
	else mexErrMsgTxt("Invalid value for permutation!");
	
	if (Symmetric!=0 && Symmetric!=1) mexErrMsgTxt("Invalid value for symmetry!");

	Options.matrix = NULL;
	Options.P = NumProcessors;
	Options.eps = Imbalance;
	
	if (!ApplyOptions(&Options)) return false;
	
	/* Start the timer. */
	Clock = clock();
	
	/* Remove double zeroes. */
	if (!SparseMatrixRemoveDuplicates(A)) return false;
	
	/* We do not have to take weighted matrices into account, because we are working in Matlab. */
	
	/* Force MATLAB-script compatible options */
	if (Symmetric)
	{
		Options.SymmetricMatrix_UseSingleEntry = SingleEntYes;
		Options.SymmetricMatrix_SingleEntryType= ETypeLower;
		SparseMatrixSymmetricRandom2Lower(A);
	}
	
	/* Add dummies if requested. */
	if (A->m == A->n && Options.SquareMatrix_DistributeVectorsEqual == EqVecYes && Options.SquareMatrix_DistributeVectorsEqual_AddDummies == DumYes) AddDummiesToSparseMatrix(A);
	
	/* Set up parameters of A. */
	A->NrProcs = NumProcessors;
	A->Pstart = (long *)malloc((A->NrProcs + 1)*sizeof(long));
	
	if (A->Pstart == NULL) return false;
	
	/* Initialise processor array. */ 
	A->Pstart[0] = 0;
	
	for (i = 1; i <= A->NrProcs; i++) A->Pstart[i] = A->NrNzElts;
	
	/* Distribute the processors among the matrix entries. */
	if (!DistributeMatrixMondriaan(A, NumProcessors, Imbalance, &Options, CallbackMaxIterations <= 0 ? NULL : MatlabCallback)) return false;

	/* Remove dummies. */
	if (A->m == A->n && Options.SquareMatrix_DistributeVectorsEqual == EqVecYes && Options.SquareMatrix_DistributeVectorsEqual_AddDummies == DumYes) RemoveDummiesFromSparseMatrix(A);
	
	/* Restore symmetry. */
	if (Symmetric && Options.SymmetricMatrix_UseSingleEntry == SingleEntYes && Options.SymmetricMatrix_SingleEntryType == ETypeRandom) SparseMatrixSymmetricRandom2Lower(A);
	
	/* Convert to full form for vector partitioning. */
	if (Symmetric && Options.SymmetricMatrix_UseSingleEntry == SingleEntYes) SparseMatrixSymmetric2Full(A);
	
	/* Distribute vectors if requested. */
	if (u != NULL && v != NULL ) {
		if (A->m == A->n && Options.SquareMatrix_DistributeVectorsEqual == EqVecYes) {
			/* Distribute vectors equally. */
			if (Symmetric && Options.SymmetricMatrix_UseSingleEntry == SingleEntYes) {
				if ((s->MaxComV = DistributeVec(A, v, ROW, &Options)) < 0) return false;
				
				for (i = 0; i < A->m; i++) u[i] = v[i];
				
				s->MaxComU = s->MaxComV;
			}
			else {
				if ((s->MaxComU = DistributeVecOrigEq(A, u, v, &Options)) < 0) return false;
				
				s->MaxComV = 0;
			}
		}
		else {
			/* Distribute vectors independently. */
			if ((s->MaxComV = DistributeVec(A, v, ROW, &Options)) < 0 || (s->MaxComU = DistributeVec(A, u, COL, &Options)) < 0) return false;
		}
	} else {
		s->MaxComU = 0;
		s->MaxComV = -1; /* So that in the statistics output, -1 appears for MaxCom, indicating no vector was distributed. */
	}
	
	/* Calculate duration of partitioning. */
	s->Duration = (double)(clock() - Clock)/(double)CLOCKS_PER_SEC;
	
	/* Determine further statistics. */
	s->MinNz = s->MaxNz = A->Pstart[1] - A->Pstart[0];
	
	for (i = 1; i < A->NrProcs; i++) {
		l = A->Pstart[i + 1] - A->Pstart[i];
		
		if (l > s->MaxNz) s->MaxNz = l;
		if (l < s->MinNz) s->MinNz = l;
	}
	
	s->Epsilon = (double)(A->NrProcs*s->MaxNz - A->NrNzElts)/(double)A->NrNzElts;

	if (v != NULL)	
		CalcCom(A, v, ROW, &s->ComVolV, &l, &l, &l, &l);
	else
		s->ComVolV = -1;
	if (u != NULL)
		CalcCom(A, u, COL, &s->ComVolU, &l, &l, &l, &l);
	else
		s->ComVolU = 0; /* So that in the statistics output, -1 appears for MaxCom, indicating no vector was distributed. */

	return	true;
}

/* This function tries to extract a double precision value from the given variable. */
double ExtractDouble(const mxArray *Var) {
	size_t NumRows, NumCols;
	
	if (!Var) mexErrMsgTxt("Null variable pointer!");
	
	NumRows = mxGetM(Var);
	NumCols = mxGetN(Var);
	
	if (mxIsComplex(Var) || NumRows != 1 || NumCols != 1) mexErrMsgTxt("Provided variable is not a real number!");
	
	return mxGetScalar(Var);
}

/*
Let us take a look at the Matlab sparse matrix format (see mxSetJc in the Matlab help).
Denote the row indices by RI, column indices by CI, and data by D.
Then the RI array contains the indices of the rows that contain non-zero elements offset by 1, so if RI[0] = 0, then the first row of the matrix contains a non-zero element.
The D array contains the values of the non-zero entries.
Finally the CI array is the most complicated: CI[i] is the index in RI and D of the first non-zero entry in the i-th column.

So if we would have the 16x5 matrix A with non-zero entries A(2, 2) = 9, A(4, 2) = 7, A(13, 5) = 2, then
RI = 2 4 13
D = 9 7 2
CI = 0 0 2 2 2 3 (no elements in first column, two elements in second column, no elements in third and fourth column, one if fifth column).

Note that this corresponds to the CCS format.
*/
/* This function converts a Matlab sparse matrix to a Mondriaan sparse matrix. */
struct sparsematrix *ConvertMatlabToMondriaan(const mxArray *A) {
	struct sparsematrix *B;
	mwIndex *AIr, *AJc;
	double *APr, *APi;
	long i, j, k;
	
	/* Read input Matlab matrix data. */
	if (A == NULL) mexErrMsgTxt ("Null variable pointer!");
	if (mxGetNumberOfDimensions (A) != 2 || !mxIsSparse (A) || mxGetM (A) <= 1 || mxGetN (A) <= 1) mexErrMsgTxt ("Matrix that we have to retrieve is not a sparse matrix!");
	
	APr = mxGetPr(A);
	APi = (mxIsComplex(A) ? mxGetPi(A) : NULL);
	AIr = mxGetIr(A);
	AJc = mxGetJc(A);
	
	if (!APr || !AIr || !AJc) mexErrMsgTxt ("Matrix data could not be retrieved!");
	
	B = (struct sparsematrix *)malloc(sizeof (struct sparsematrix));
	
	if (B == NULL) mexErrMsgTxt("Could not allocate Mondriaan matrix!");
	
	/* First initialise the sparse matrix to its default parameters. */
	if (!MMSparseMatrixInit(B)) mexErrMsgTxt("Unable to initialise matrix!");
	
	/* Now fill in the parameters from A. */
	B->MMTypeCode[0] = 'M';
	B->MMTypeCode[1] = 'C';
	B->MMTypeCode[2] = (APi != NULL ? 'C' : 'R');
	B->MMTypeCode[3] = 'G';
	
	B->m = mxGetM(A);
	B->n = mxGetN(A);
	B->NrNzElts = AJc[B->n];
	
	/* Allocate arrays to hold the data of all matrix elements. */
	if (!MMSparseMatrixAllocateMemory(B)) mexErrMsgTxt("Unable to allocate matrix!");
	if (B->i == NULL || B->j == NULL || B->ReValue == NULL || (B->ImValue == NULL && APi != NULL)) mexErrMsgTxt ("Unable to allocate Mondriaan matrix arrays!");
	
	/* Fill the arrays with data from A. */
	i = 0;
	
	for (j = 0; j < B->n && AJc[j] < (size_t)B->NrNzElts; j++) {
		for (k = AJc[j + 1] - AJc[j]; k-- > 0; ) {
			B->i[i] = AIr[i];
			B->j[i] = j;
			B->ReValue[i] = APr[i];
			
			if (B->ImValue != NULL) B->ImValue[i] = APi[i];
			
			i++;
		}
	}
	
	return B;
}
/* Compare function necessary to order triplets in a fashion that is compatible with the Matlab sparse matrix format. */
int CompareTriplets(const void *p1, const void *p2) {
	const struct MondriaanTriplet *t1, *t2;
	
	t1 = (struct MondriaanTriplet *)p1;
	t2 = (struct MondriaanTriplet *)p2;
	
	if (t1->Column != t2->Column) return t1->Column - t2->Column;
	
	return t1->Row - t2->Row;
}

/* This function converts Mondriaan statistics to a matlab vector. */
mxArray *ConvertStatsToVector(const struct MondriaanStats *s) {
	mxArray *w;
	double *wp;
	
	if (s == NULL) mexErrMsgTxt("Null stats pointer!");
	
	w = mxCreateDoubleMatrix(4, 1, mxREAL);
	
	if (w == NULL) mexErrMsgTxt("Unable to allocate Matlab vector!");
	
	wp = mxGetPr(w);
	
	wp[0] = s->Duration;
	wp[1] = s->Epsilon;
	wp[2] = (double)(s->MaxComU + s->MaxComV);
	wp[3] = (double)(s->ComVolU + s->ComVolV);
	
	return w;
}

/* This function returns the matlab vector [value]. */
mxArray *OneVector(const double value) {
	mxArray *w = mxCreateDoubleMatrix(1, 1, mxREAL);
	if (w == NULL) mexErrMsgTxt("Unable to allocate Matlab vector!");
	mxGetPr(w)[0] = value;
	return w;
}

/* This function returns the matlab vector [value1;value2]. */
mxArray *TwoVector(const double value1, const double value2) {
	double *wp;
	mxArray *w = mxCreateDoubleMatrix(2, 1, mxREAL);
	if (w == NULL) mexErrMsgTxt("Unable to allocate Matlab vector!");
	wp		= mxGetPr(w);
	wp[0] = value1;
	wp[1] = value2;
	return w;
}

/* Thus function returns the matlab vector [1 2... Size]. */
mxArray *IncrementalVector(const int Size) {
	mxArray *w; double *wp; int i;
	w = mxCreateDoubleMatrix(Size, 1, mxREAL);
	if (w == NULL) mexErrMsgTxt("Unable to allocate Matlab vector!");
	wp = mxGetPr(w);
	for (i = 0; i < Size; i++) wp[i] = (double)(i+1);
	return w;
}

/* This function converts a vector assignment to a vector of appropriate size. */
mxArray *ConvertBoundaryIndicesToVector(const struct remembrance *m) {
	mxArray *w;
	size_t i;
	double *wp;

	if (m == NULL) mexErrMsgTxt("Null variable pointer!");

	w = mxCreateDoubleMatrix(m->size, 1, mxREAL);

	if (w == NULL) mexErrMsgTxt("Unable to allocate Matlab vector!");

	wp = mxGetPr(w);

	for (i = 0; i < m->size; i++) wp[i] = (double)(m->vector[i].index + 1);

	return w;
}

/* This function converts a vector assignment to a vector of appropriate size. */
mxArray *ConvertBoundaryHierarchyToVector(const struct remembrance *m) {
	mxArray *w;
	size_t i;
	double *wp;

	if (m == NULL) mexErrMsgTxt("Null variable pointer!");

	w = mxCreateDoubleMatrix(m->size-1, 1, mxREAL);

	if (w == NULL) mexErrMsgTxt("Unable to allocate Matlab vector!");

	wp = mxGetPr(w);

	for (i = 0; i < m->size-1; i++)
	{
		if (m->vector[i].parent == ULONG_MAX)
			wp[i] = 0;
		else
			wp[i] = (double)(m->vector[i].parent + 1);
	}

	return w;
}

/* This function converts a vector assignment to a vector of appropriate size. */
mxArray *ConvertLongIndicesToVector(const long *v, const size_t Size) {
	mxArray *w;
	size_t i;
	double *wp;

	if (v == NULL) mexErrMsgTxt("Null variable pointer!");

	w = mxCreateDoubleMatrix(Size, 1, mxREAL);

	if (w == NULL) mexErrMsgTxt("Unable to allocate Matlab vector!");

	wp = mxGetPr(w);

	for (i = 0; i < Size; i++) wp[i] = (double)(v[i] + 1);

	return w;
}

/* This function converts a vector assignment to a vector of appropriate size. */
mxArray *ConvertIndicesToVector(const int *v, const size_t Size) {
	mxArray *w;
	size_t i;
	double *wp;
	
	if (v == NULL) mexErrMsgTxt("Null variable pointer!");
	
	w = mxCreateDoubleMatrix(Size, 1, mxREAL);
	
	if (w == NULL) mexErrMsgTxt("Unable to allocate Matlab vector!");
	
	wp = mxGetPr(w);
	
	for (i = 0; i < Size; i++) wp[i] = (double)(v[i] + 1);
	
	return w;
}

/* This function converts a Mondriaan sparse matrix to a Matlab sparse matrix, Triplets are already created. */
mxArray *ConvertTripletsToMatlab(struct MondriaanTriplet *Triplets, const struct sparsematrix *B) {
	mxArray *A;
	mwIndex *AIr, *AJc;
	double *APr, *APi;
	int CurrentColumn;
	long i;
	
	if (Triplets == NULL) mexErrMsgTxt("Null variable pointer!");
	
	/* Sort them in a Matlab friendly way. */
	qsort(Triplets, B->NrNzElts, sizeof (struct MondriaanTriplet), CompareTriplets);
	
	/* Allocate Matlab matrix and obtain pointers to its real data, imaginary data, row index, and column count arrays respectively. */
	A = mxCreateSparse(B->m, B->n, B->NrNzElts, B->ImValue == NULL ? mxREAL : mxCOMPLEX);
	
	if (A == NULL) mexErrMsgTxt("Unable to allocate Matlab matrix!");
	
	APr = mxGetPr (A);
	APi = mxGetPi (A);
	AIr = mxGetIr (A);
	AJc = mxGetJc (A);
	
	/* Write triplets to the Matlab matrix. */
	CurrentColumn = 0;
	*AJc++ = 0;
	
	for (i = 0; i < B->NrNzElts; i++) {
		*AIr++ = Triplets[i].Row;
		*APr++ = Triplets[i].ReValue;
		if (APi) *APi++ = Triplets[i].ImValue;
		
		/* Check whether or not we need to advance to the next column. */
		while (CurrentColumn < Triplets[i].Column) {
			*AJc++ = i;
			CurrentColumn++;
		}
	}
	
	if (CurrentColumn > B->n) mexErrMsgTxt( "Current column advanced farther than end of matrix!" );
	if (i != B->NrNzElts)		 mexErrMsgTxt( "Not all nonzeroes are written into Matlab matrix!"	 );

	/* This ensures exactly n writes to AJc in this loop together with the
	 while loop above. Note that before the for loop above, the 0th element
	 already has been set, thus totalling n+1 elements in AJc, so the
	 correct number of elements is written. Also, the below loop fires at
	 least once, since CurrentColumn is 0-based (and exits above for loop
	 at value n-1, at most). */
	for ( ; CurrentColumn < B->n; CurrentColumn++)
		*AJc++ = B->NrNzElts;
	
	free(Triplets);
	
	return A;
}

/* This function converts a Mondriaan sparse matrix to a Matlab sparse matrix. */
mxArray *ConvertMondriaanToMatlab(const struct sparsematrix *B) {
	struct MondriaanTriplet *Triplets;
	long i;
	
	if (B == NULL) mexErrMsgTxt("Null variable pointer!");
	
	/* Extract all non-zero entries. */
	Triplets = (struct MondriaanTriplet *)malloc(B->NrNzElts*sizeof (struct MondriaanTriplet));
	
	if (Triplets == NULL) mexErrMsgTxt("Unable to allocate triplets!");
	
	for (i = 0; i < B->NrNzElts; i++) {
		Triplets[i].Row = B->i[i];
		Triplets[i].Column = B->j[i];
		Triplets[i].ReValue = B->ReValue[i];
		if (B->ImValue) Triplets[i].ImValue = B->ImValue[i];
		else Triplets[i].ImValue = 0.0;
	}
	return ConvertTripletsToMatlab( Triplets, B );
}
	
/* This function converts a Mondriaan sparse matrix to a Matlab sparse matrix with processor indices as nonzero values. */
mxArray *ConvertMondriaanProcIndToMatlab(const struct sparsematrix *B) {
	struct MondriaanTriplet *Triplets;
	long i, s, k;

	if (B == NULL) mexErrMsgTxt("Null variable pointer!");
	
	/* Extract all non-zero entries. */
	Triplets = (struct MondriaanTriplet *)malloc(B->NrNzElts*sizeof (struct MondriaanTriplet));
	
	if (Triplets == NULL) mexErrMsgTxt("Unable to allocate triplets!");
	
	i = 0;
	for (s = 0; s < B->NrProcs; s++) {
		for (k = B->Pstart[s]; k < B->Pstart[s+1]; k++) {
			Triplets[i].Row = B->i[k];
			Triplets[i].Column = B->j[k];
			Triplets[i].ReValue = (double)(s+1);
			Triplets[i].ImValue = 0.0;
			i++;
		}
	}
	return ConvertTripletsToMatlab( Triplets, B );
}

