#ifndef __SparseMatrix_h__
#define __SparseMatrix_h__

#include "Options.h"
#include "Sort.h"
#include "Remembrance.h"

#define MM_Banner           "%%MatrixMarket"
#define EMM_Banner          "%%Extended-MatrixMarket"

#define ViewTypeOriginal 0
#define ViewTypeGlobal   1
#define ViewTypeLocal    2

/* Sparse matrix data structure containing a list of triples (i,j,a)
   representing nonzeros, where i is the row index, j the column index,
   and a the numerical value */
struct sparsematrix {
  char MMTypeCode[4]; 
  char Banner[MAX_WORD_LENGTH];
  char Object[MAX_WORD_LENGTH];
  char Format[MAX_WORD_LENGTH];
  char Field[MAX_WORD_LENGTH];
  char Symmetry[MAX_WORD_LENGTH];
  char ViewType;

  /* Matrix dimensions */
  long NrNzElts;  /* Number of nonzero matrix elements */
  long m;         /* Number of matrix rows */
  long n;         /* Number of matrix columns */
  long NrDummies; /* Number of dummy diagonal nonzeros */
  int  NrProcs;   /* Number of processors for a distributed matrix */

  /* Arrays of length NrNzElts storing the triples */
  long *i;         /* List of row indices */
  long *j;         /* List of column indices */
  double *ReValue; /* List of Real part of the matrix values */
  double *ImValue; /* List of Imaginary part of the matrix values */
  
  /* Additional arrays */
  int *dummy;   /* Boolean array which says whether a diagonal element
                   is dummy or not */
  long *Pstart; /* Pointer to the start of the triples for the processors */

  char *header; /* Header comment of the matrix */
  char *comment;/* General comments included in this matrix */
  char *tail;   /* Tail comment of the matrix */

  long NrRowWeights; /* Number of row weights, = 0 or m */
  long NrColWeights; /* Number of column weights, = 0 or n */
  long *RowWeights;  /* Array of row weights */
  long *ColWeights;  /* Array of column weights */
  
  int *RowLambda;   /* Array containing the number of different parts in a certain row. */
  int *ColLambda;   /* Array containing the number of different parts in a certain column. */

  int *RowMark; /* Array for marking cut rows */
  long NumRowsMarked;
  int *ColMark;
  long NumColsMarked;
  
  long mgMid; /* For mediumgrain initial split */
  int mgDir; /* For mediumgrain iterative improvement */

  /* Keeps track of permutations. */
  long *row_perm, *col_perm;
  long *row_perm_inv, *col_perm_inv;

  /* Keeps track of permutation blocks. */
  struct remembrance *rowBoundaries, *colBoundaries;
};

/* Function declarations for Matrix Market input/output functions */
int MMSparseMatrixInit(struct sparsematrix *pM);
int MMWeightsInit(struct sparsematrix *pM);
int MMReadSparseMatrix(FILE *fp, struct sparsematrix *pM);
int MMWriteSparseMatrix(struct sparsematrix *pM, FILE *fp, const char *name, const struct opts * pOptions);
int MMDeleteSparseMatrix(struct sparsematrix *pM);
int MMSparseMatrixSetTypeCode(const char *line, struct sparsematrix *pM);
int MMSparseMatrixAppendViewType(const char view, char *line);
int MMSparseMatrixGetTypeCode(struct sparsematrix *pM, char *line, const char *name, const struct opts * pOptions);
int MMSparseMatrixAllocateMemory(struct sparsematrix *pM);
int MMSparseMatrixFreeMemory(struct sparsematrix *pM);

int MMInsertProcessorIndices(struct sparsematrix *pM);

int MMSparseMatrixReadHeader(struct sparsematrix *pM, FILE *fp);
int MMSparseMatrixReadPstart(struct sparsematrix *pM, FILE *fp);
int MMSparseMatrixReadWeights(struct sparsematrix *pM, FILE *fp);
int MMSparseMatrixReadEntries(struct sparsematrix *pM, FILE *fp);
int MMSparseMatrixReadTail(struct sparsematrix *pM, FILE *fp);

int MMSparseMatrixPrintHeader(struct sparsematrix *pM, FILE *stream, const char *name, const struct opts * pOptions);
int MMSparseMatrixPrintComment(struct sparsematrix *pM, FILE *stream);
int MMSparseMatrixPrintPstart(struct sparsematrix *pM, FILE *stream, const struct opts * pOptions);
int MMSparseMatrixPrintWeights(struct sparsematrix *pM, FILE *stream);
int MMSparseMatrixPrintEntries(struct sparsematrix *pM, FILE *stream);
int MMSparseMatrixPrintTail(struct sparsematrix *pM, FILE *stream);

/* Function declarations for input functions for other sparse matrix storage schemes */
int CRSSparseMatrixInit(struct sparsematrix *pM, const long m, const long n, const long nnz,
			const int *col, const int *row_start, const double *values,
			char base );
int CRSSparseMatrixNetWeightInit(struct sparsematrix *pM, const long m, const long n, const long nnz,
			const int *col, const int *row_start, const double *values,
			char base, double scale );

/* Utility functions */
/* Allocates and initialises the Pstart array, assumes NrNzElts and NrProcs are set. */
int PstartInit(struct sparsematrix *pM, const int P );

/* Function declarations for matrix conversion functions */
int AddDummiesToSparseMatrix(struct sparsematrix *pM);
int RemoveDummiesFromSparseMatrix(struct sparsematrix *pM);

int SparseMatrixSymmetric2Full(struct sparsematrix *pM);
int SparseMatrixSymmetricLower2Random(struct sparsematrix *pM);
int SparseMatrixSymmetricRandom2Lower(struct sparsematrix *pM);

int SparseMatrixRemoveDuplicates(struct sparsematrix *pM);
int SparseMatrixStructurallySymmetric(struct sparsematrix *pM);

int SparseMatrixMarkCutRows(struct sparsematrix *pM);
int SparseMatrixMarkCutColumns(struct sparsematrix *pM);

int SparseMatrixOriginal2Local(struct sparsematrix *pM, long int **row_perms, long int **col_perms);
int SparseMatrixLocal2Vector(struct sparsematrix *pM, long int **local2glob, long int *vec_distr,
                             long int **count, long int ***local2proc, long int ***local2index, 
                             const char dir );
int CreateInitialMediumGrainDistribution(struct sparsematrix *pM, long *mid);

#endif /* __SparseMatrix_h__ */

