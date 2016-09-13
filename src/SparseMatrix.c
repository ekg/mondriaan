#include "SparseMatrix.h"

int CRSSparseMatrixInit(struct sparsematrix *pM, const long m, const long n, const long nnz,
                        const int *col, const int *row_start, const double *values,
                        char base ) {
	return CRSSparseMatrixNetWeightInit( pM, m, n, nnz, col, row_start, values, 0, 0.0 );
}

int CRSSparseMatrixNetWeightInit(struct sparsematrix *pM, const long m, const long n, const long nnz,
                        const int *col, const int *row_start, const double *values,
			char base, double scale ) {
    int i, k;
    double *RowW = NULL, *ColW = NULL;
    long c = 0;
    
    if( !MMSparseMatrixInit(pM) ) {
        fprintf( stderr, "Error: could not initialise sparse matrix for Mondriaan!\n" );
        return FALSE;
    }

    pM->MMTypeCode[ 0 ] = 'M';
    pM->MMTypeCode[ 1 ] = 'C';
    pM->MMTypeCode[ 2 ] = 'R';
    pM->MMTypeCode[ 3 ] = 'G';
    pM->m = m; pM->n = n; pM->NrNzElts = nnz;

    if( !MMSparseMatrixAllocateMemory(pM) ) {
        fprintf( stderr, "Error: could not allocate sparse matrix (SparseMatrix.c::CRSSParseMatrixInit)!\n" );
        return FALSE;
    }

    if( !MMWeightsInit( pM ) ) {
        fprintf(stderr, "CRSSparseMatrixInit(): error during row/column weight initialisation!\n" );
        return FALSE;
    }
    
    if( scale != 0.0 ) {
        RowW = malloc( m * sizeof( double ) );
        ColW = malloc( n * sizeof( double ) );
        for( i=0; i<m; i++ ) RowW[i] = 0.0;
        for( i=0; i<n; i++ ) ColW[i] = 0.0;
    }

    for( i=0; i<m; i++ ) /* Loop over each row. */
        for( k=row_start[i]-base; k<row_start[i+1]-base; k++ ) { /* For each nonzero in that row. */
            pM->i[c] = i;
            pM->j[c] = col[k]-base;
            pM->ReValue[c++] = values[k]; /* Increment triplet count. */
            if( scale != 0.0 ) { /* Let the compiler optimise this. */
                RowW[i] += fabs( values[k] );
                ColW[col[k]-base] += fabs( values[k] );
            }
        }
    /* TODO complex? (take square norm?) */
    if( scale != 0.0 ) {
        for( i=0; i<m; i++ ) pM->RowWeights[i] = scale * ( RowW[i] + 0.5 ); /* Round up. */
        for( i=0; i<n; i++ ) pM->ColWeights[i] = scale * ( ColW[i] + 0.5 ); /* Round up. */
        free( RowW );
        free( ColW );
    }
    return TRUE;
}

int PstartInit(struct sparsematrix *pM, const int P) {
    int k;
    
    pM->NrProcs = P;
    pM->Pstart = (long *) malloc( (pM->NrProcs+1) * sizeof(long) );
    if( pM->Pstart == NULL ) {
        fprintf( stderr, "Error: out of memory!\n" );
        return FALSE;
    }
    pM->Pstart[0] = 0;
    for( k = 1; k <= pM->NrProcs; k++ )
       pM->Pstart[ k ] = pM->NrNzElts;
    return TRUE;
}

int MMWeightsInit(struct sparsematrix *pM) {
    /* This function initialises the RowWeights and
       ColWeights arrays and sets NrRowWeights/
       NrColWeights to m/n, respectively. */
    if( pM->RowWeights != NULL )
        fprintf(stderr, "Warning: MMWeightsInit(): RowWeights already initialised; is this a weighted matrix?\n");
    else {
        const long m = pM->m;
        pM->NrRowWeights = m;
        pM->RowWeights = malloc( m * sizeof( long ) );
    }
    if( pM->ColWeights != NULL )
        fprintf(stderr, "Warning: MMWeightsInit(): ColWeights already initialised; is this a weighted matrix?\n");
    else {
        const long n = pM->n;
        pM->NrColWeights = n;
        pM->ColWeights = malloc( n * sizeof( long ) );
    }
    if( pM->RowWeights == NULL || pM->ColWeights == NULL ) {
        fprintf(stderr, "MMWeightsInit(): Out of memory!\n");
        return FALSE;
    }
    return TRUE;
}

int MMSparseMatrixInit(struct sparsematrix *pM) {

    /* This function initialises a sparse matrix by setting all its
       numerical variables to 0, all character variables to a blank,
       all string variables to the empty string,
       and all pointers to NULL. It does not allocate memory. */
    if (!pM) {
        fprintf(stderr, "MMSparseMatrixInit(): No matrix provided!\n");
        return FALSE;
    }

    /* Initialise typeparams */
    pM->MMTypeCode[0] = ' ';
    pM->MMTypeCode[1] = ' ';
    pM->MMTypeCode[2] = ' ';
    pM->MMTypeCode[3] = ' ';

    /* Initialise strings */
    pM->Banner[0] = '\0';
    pM->Object[0] = '\0';
    pM->Format[0] = '\0';
    pM->Field[0] = '\0';
    pM->Symmetry[0] = '\0';
    pM->ViewType = ViewTypeOriginal;

    /* Initialise matrix dimensions */
    pM->NrNzElts = 0;  /* no nonzeros yet */
    pM->m = 0;
    pM->n = 0;
    pM->NrDummies = 0; /* no dummies added yet */
    pM->NrProcs = 0;   /* matrix has not been partitioned yet */

    /* Initialise matrix arrays of length NrNzElts */
    pM->i = NULL;
    pM->j = NULL;
    pM->ReValue = NULL;
    pM->ImValue = NULL;

    /* Initialise other matrix arrays */
    pM->dummy = NULL;
    pM->Pstart = NULL;
    pM->header = NULL;
    pM->comment= NULL;
    pM->tail = NULL;

    /* Initialise weight information */
    pM->NrRowWeights = 0;
    pM->NrColWeights = 0;
    pM->RowWeights = NULL;
    pM->ColWeights = NULL;

    /* Initialise marker arrays */
    pM->RowLambda = NULL;
    pM->ColLambda = NULL;
    pM->RowMark = NULL;
    pM->ColMark = NULL;
    
    /* Permutation data */
    pM->row_perm = NULL;
    pM->col_perm = NULL;
    pM->row_perm_inv = NULL;
    pM->col_perm_inv = NULL;

    /* Block data */
    pM->rowBoundaries = pM->colBoundaries = NULL;

    return TRUE;
} /* end MMSparseMatrixInit */
  

int MMReadSparseMatrix(FILE *fp, struct sparsematrix *pM) {

    /* This function reads the header and the nonzeros of a sparse
       matrix in Matrix Market format, and possibly a tail.
       The header is an optional banner line
           %%MatrixMarket ...
       followed by a line
           m n nz
       giving the number of rows, columns, nonzeros.

       If the matrix has been partitioned (distributed-matrix),
       the function reads distribution information.
       In that case the banner line is followed by a line
           m n nz P
       where P is the number of processors, and then by
       P+1 lines giving the start of the nonzeros of the processors,
       where the (P+1)-th line equals nz.

       If the matrix is weighted (weighted-matrix),
       the function reads m integer weights for the rows
       and then n weights for the columns. */

    if (!fp || !pM) {
        fprintf(stderr, "MMReadSparseMatrix(): Null arguments!\n");
        return FALSE;
    }
  
    if (!MMSparseMatrixInit(pM)) {
        fprintf(stderr, "MMReadSparseMatrix(): Could not initialise matrix!\n");
        return FALSE;
    }
    
    /* Open the file */
    if (ferror(fp)) {
        fprintf(stderr, "MMReadSparseMatrix(): Could not read from stream!\n");
        return FALSE;
    }
  
    /* Read the header */
    if (!MMSparseMatrixReadHeader(pM, fp)) return FALSE;
  
    /* Allocate memory for the arrays of the matrix */
    if (!MMSparseMatrixAllocateMemory(pM)) return FALSE;
  
    /* Read the distribution information */
    if (!MMSparseMatrixReadPstart(pM, fp)) return FALSE;
  
    /* Read the nonzero entries of the matrix */
    if (!MMSparseMatrixReadEntries(pM, fp)) return FALSE;
  
    /* Read the row and column weights */
    if (!MMSparseMatrixReadWeights(pM, fp)) return FALSE;
  
    /* Read the tail */
    if (!MMSparseMatrixReadTail(pM, fp)) return FALSE;
  
    return TRUE;

} /* end MMReadSparseMatrix */

int MMWriteSparseMatrix(struct sparsematrix *pM, FILE *fp, const char* name, const struct opts * pOptions) {

    /* This function prints the header and the nonzeros of a sparse
       matrix in Matrix Market format, and possibly a tail,
       into a file.
       
       The header is an optional banner line
           %%MatrixMarket ...
       followed by a line
           m n nz
       giving the number of rows, columns, nonzeros.

       If the matrix has been partitioned (distributed-matrix),
       the function prints distribution information.
       In that case the banner line is followed by a line
           m n nz P
       where P is the number of processors, and then by
       P+1 lines giving the start of the nonzeros of the processors,
       where the (P+1)-th line equals nz.

       If the matrix is weighted (weighted-matrix),
       the function prints m integer weights for the rows
       and then n weights for the columns. */ 

    if (!pM || !fp) {
        fprintf(stderr, "MMWriteSparseMatrix(): Null arguments!\n");
        return FALSE;
    }
    
    if (!pM) {
        fprintf(stderr, "MMWriteSparseMatrix(): No matrix provided!\n");
        return FALSE;
    }
  
    /* Open the file */
    if (ferror(fp)) {
        fprintf(stderr, "MMWriteSparseMatrix(): Could not write to stream!\n");
        return FALSE;
    }
  
    /* Print the header */
    if (!MMSparseMatrixPrintHeader(pM, fp, name, pOptions)) return FALSE;
 
    /* Print the distribution information */
    if (!MMSparseMatrixPrintPstart(pM, fp, pOptions)) return FALSE;
  
    /* Print the nonzero elements */
    if (!MMSparseMatrixPrintEntries(pM, fp)) return FALSE;
  
    /* Print the row and column weights */
    if (!MMSparseMatrixPrintWeights(pM, fp)) return FALSE;
  
    /* Print the tail */
    if (!MMSparseMatrixPrintTail(pM, fp)) return FALSE;

    return TRUE;

} /* end MMWriteSparseMatrix */

  
int MMDeleteSparseMatrix(struct sparsematrix *pM) {

    /* This function frees the memory for all the arrays of the matrix M,
       and reinitialises the variables so the matrix can be reused. 
       It sets the pointers to all arrays to NULL. */
    if (!pM) {
        fprintf(stderr, "MMDeleteSparseMatrix(): No matrix provided!\n");
        return FALSE;
    }
  
    MMSparseMatrixFreeMemory(pM);
    MMSparseMatrixInit(pM);

    return TRUE;
} /* end MMDeleteSparseMatrix */
  
  
int MMSparseMatrixSetTypeCode(const char *line, struct sparsematrix *pM) {

    /* This function reads a Matrix Market banner line,
       splits it into 5 strings (for banner, object, format, field, symmetry),
       reads these into the corresponding strings of the matrix M,
       and sets the type codes accordingly. 

       The first string is the banner string. The possibilities for the
       other strings are:

       Object = matrix, distributed matrix, or weighted matrix 
       TypeCode[0] = M, D, W

       Format = coordinate, array
       TypeCode[1] = C, A

       Field =  integer, real, complex, or pattern
       TypeCode[2] = I, R, C, P

       Symmetry = general, symmetric, skew-symmetric, or hermitian
       TypeCode[3] = G, S, K, H

       If there are less than 5 strings in the line,
       or the first string differs from the Matrix Market banner,
       the four type codes are set to the defaults M, C, R, G, respectively,
       and the matrix strings are initialised accordingly.  */

    char *c;
    
    if (!line || !pM) {
        fprintf(stderr, "MMSparseMatrixSetTypeCode(): No matrix or line provided!\n");
        return FALSE;
    }
  
    /* Read Matrix Market Type Codes from line */
    if (sscanf(line, "%s %s %s %s %s", pM->Banner, pM->Object, 
        pM->Format, pM->Field, pM->Symmetry) != 5 ||
        strcmp(pM->Banner, MM_Banner)) { 

        /* No correct Type Code found, assume default */
        strcpy(pM->Banner, MM_Banner);
        pM->MMTypeCode[0] = 'M'; strcpy(pM->Object, "matrix");
        pM->MMTypeCode[1] = 'C'; strcpy(pM->Format, "coordinate");
        pM->MMTypeCode[2] = 'R'; strcpy(pM->Field, "real");
        pM->MMTypeCode[3] = 'G'; strcpy(pM->Symmetry, "general");
        return TRUE;
    }
  
    /*** Set Object type in M.MMTypeCode[0] ***/

    /* convert string to lower case */
    for (c = pM->Object; *c != '\0'; c++)
        *c = tolower(*c);
  
    /* set Object type = matrix, distributed matrix, or weighted matrix */
    if (! strcmp(pM->Object, "matrix"))
        /* pM->Object = "matrix" */
        pM->MMTypeCode[0] = 'M';
    else if (! strcmp(pM->Object, "distributed-matrix") ||
              ! strcmp(pM->Object, "distributedmatrix"))
        pM->MMTypeCode[0] = 'D';      
    else if (! strcmp(pM->Object, "weighted-matrix") ||
              ! strcmp(pM->Object, "weightedmatrix"))
        pM->MMTypeCode[0] = 'W';      
    else {
        fprintf(stderr, "MMSparseMatrixSetTypeCode(): Unknown Object '%s'!\n", pM->Object);
        return FALSE;
    }
    
    /*** Set Format type in M.MMTypeCode[1] ***/

    for (c = pM->Format; *c != '\0'; c++)
        *c = tolower(*c);
  
    /* set Format type = coordinate or array */
    if (! strcmp(pM->Format, "coordinate") || 
         ! strcmp(pM->Format, "sparse"))
        pM->MMTypeCode[1] = 'C';
    else if (! strcmp(pM->Format, "array") || 
              ! strcmp(pM->Format, "dense"))
        pM->MMTypeCode[1] = 'A';
    else {
        fprintf(stderr, "MMSparseMatrixSetTypeCode(): Unknown Format '%s'!\n", pM->Format);
        return FALSE;
    }
    
    /*** Set Field type in M.MMTypeCode[2] ***/

    for (c = pM->Field; *c != '\0'; c++)
        *c = tolower(*c);
  
    /* set Field type = integer, real, complex, or pattern */
    if (! strcmp(pM->Field, "integer"))
        pM->MMTypeCode[2] = 'I';
    else if (! strcmp(pM->Field, "real") || 
              ! strcmp(pM->Field, "float"))
        pM->MMTypeCode[2] = 'R';
    else if (! strcmp(pM->Field, "complex"))
        pM->MMTypeCode[2] = 'C';
    else if (! strcmp(pM->Field, "pattern"))
        pM->MMTypeCode[2] = 'P';
    else {
        fprintf(stderr, "MMSparseMatrixSetTypeCode(): Unknown Field '%s'!\n", pM->Field);
        return FALSE;
    }
  
    /*** Set Symmetry type in M.MMTypeCode[3] ***/

    for (c = pM->Symmetry; *c != '\0'; c++)
        *c = tolower(*c);
  
    /* set Symmetry type = general, symmetric, skew-symmetric, or hermitian */
    if (! strcmp(pM->Symmetry, "general"))
        pM->MMTypeCode[3] = 'G';
    else if (! strcmp(pM->Symmetry, "symmetric"))
        pM->MMTypeCode[3] = 'S';
    else if (! strcmp(pM->Symmetry, "skew-symmetric") ||
              ! strcmp(pM->Symmetry, "skewsymmetric"))
        pM->MMTypeCode[3] = 'K';
    else if (! strcmp(pM->Symmetry, "hermitian"))
        pM->MMTypeCode[3] = 'H';
    else {
        fprintf(stderr, "MMSparseMatrixSetTypeCode(): Unknown Symmetry '%s'!\n", pM->Symmetry);
        return FALSE;
    }

    return TRUE;
} /* end MMSparseMatrixSetTypeCode */
  
int MMSparseMatrixAppendViewType(const char view, char *line) {

    /* This function appends the viewtype of a matrix to an EMM banner line. */

    if(!line) {
        fprintf(stderr, "MMSparseMatrixAppendViewType(): null argument provided!\n");
        return FALSE;
    }
    switch(view) {
    case ViewTypeOriginal:
        sprintf(line, "%s original", line);
        break;
    case ViewTypeGlobal:
        sprintf(line, "%s global", line);
        break;
    case ViewTypeLocal:
        sprintf(line, "%s local", line);
        break;
    default:
        fprintf(stderr, "MMSparseMatrixAppendViewType(): unidentified matrix ViewType!\n");
        return FALSE;
    }
    return TRUE;
} /* end MMSparseMatrixAppendViewType */

int MMSparseMatrixGetTypeCode(struct sparsematrix *pM, char *line, const char *name, const struct opts * pOptions) {

    /* This function reads the type codes from a matrix M,
       converts these into a Matrix Market banner line
       and copies this line to line.

       Furthermore, the 5 strings of the generated banner line
       (for banner, object, format, field, symmetry)
       are copied into the corresponding strings of the matrix M. 
 
       The type codes are converted into strings as follows:

       TypeCode[0] = M, D, W
       Object = matrix, distributed matrix, or weighted matrix

       TypeCode[1] = C, A
       Format = coordinate, array

       TypeCode[2] = I, R, C, P
       Field =  integer, real, complex, or pattern

       TypeCode[3] = G, S, K, H
       Symmetry = general, symmetric, skew-symmetric, or hermitian

    */
    if (!pM || !line) {
        fprintf(stderr, "MMSparseMatrixGetTypeCode(): No matrix or line provided!\n");
        return FALSE;
    }

    /* Set Banner */
    if (name != NULL)
        sprintf(pM->Banner, "%%%%%s", name);
    else {
        if (pOptions->OutputFormat == OutputDMM)
            strcpy(pM->Banner, MM_Banner);
        else if (pOptions->OutputFormat == OutputEMM)
            strcpy(pM->Banner, EMM_Banner);
        else {
            fprintf( stderr, "MMSparseMatrixGetTypeCode(): Undefined output format!\n" );
            return FALSE;
        }
    }
  
    /* Set Object type */
    if (pM->MMTypeCode[0] == 'M')
        strcpy(pM->Object, "matrix");
    else if (pM->MMTypeCode[0] == 'D')
        strcpy(pM->Object, "distributed-matrix");
    else if (pM->MMTypeCode[0] == 'W')
        strcpy(pM->Object, "weighted-matrix");
    else {
        fprintf(stderr, "MMSparseMatrixGetTypeCode(): Unknown Object Type Code '%c'!\n", pM->MMTypeCode[0]);
        return FALSE;
    }
  
    /* Set Format type */
    if (pM->MMTypeCode[1] == 'C')
        strcpy(pM->Format, "coordinate");
    else if (pM->MMTypeCode[1] == 'A')
        strcpy(pM->Format, "array");
    else {
        fprintf(stderr, "MMSparseMatrixGetTypeCode(): Unknown Format Type Code '%c'!\n", pM->MMTypeCode[1]);
        return FALSE;
    }
  
    /* Set Field type */
    if (pM->MMTypeCode[2] == 'I')
        strcpy (pM->Field, "integer");
    else if (pM->MMTypeCode[2] == 'R')
        strcpy(pM->Field, "real");
    else if (pM->MMTypeCode[2] == 'C')
        strcpy(pM->Field, "complex");
    else if (pM->MMTypeCode[2] == 'P')
        strcpy(pM->Field, "pattern");
    else {
        fprintf(stderr, "MMSparseMatrixGetTypeCode(): Unknown Field Type Code '%c'!\n", pM->MMTypeCode[2]);
        return FALSE;
    }
  
    /* Set Symmetry type */
    if (pM->MMTypeCode[3] == 'G')
        strcpy(pM->Symmetry, "general");
    else if (pM->MMTypeCode[3] == 'S')
        strcpy(pM->Symmetry, "symmetric");    
    else if (pM->MMTypeCode[3] == 'K')
        strcpy(pM->Symmetry, "skew-symmetric") ;
    else if (pM->MMTypeCode[3] == 'H')
        strcpy(pM->Symmetry, "hermitian");
    else {
        fprintf(stderr, "MMSparseMatrixGetTypeCode(): Unknown Symmetry Type Code '%c'!\n", pM->MMTypeCode[3]);
        return FALSE;
    }
  
    /* Copy the type parameters to the banner line */
    if (pOptions->OutputFormat == OutputDMM)
        sprintf(line, "%s %s %s %s %s", pM->Banner, pM->Object, pM->Format, 
                pM->Field, pM->Symmetry);
    else if (pOptions->OutputFormat == OutputEMM) {
        sprintf(line, "%s %s %s %s %s", pM->Banner, pM->Object, pM->Format, 
                pM->Field, pM->Symmetry);
        if (!MMSparseMatrixAppendViewType(pM->ViewType, line))
             return FALSE;
    } else return FALSE;

    sprintf(line, "%s\n", line);
    return TRUE;
} /* end MMSparseMatrixGetTypeCode */
  
  
int MMSparseMatrixAllocateMemory(struct sparsematrix *pM) {

    /* This function allocates memory for the following arrays of the matrix M:
           - the nonzero arrays (i, j, Val) 
           - Pstart, which represents the start of the nonzeros
             of the processors, for a distributed matrix
           - the row and column weights, for a weighted matrix.

       It does not allocate arrays for dummies or row/column marks,
       since these are only allocated where needed, by other functions. */
    
    if (!pM) {
        fprintf(stderr, "MMSparseMatrixAllocateMemory(): No matrix provided!\n");
        return FALSE;
    }

    /* Allocate memory for indices. Allocate one element extra,
       in case the matrix is empty. */
    pM->i = (long *) malloc((pM->NrNzElts+1) * sizeof(long));
    pM->j = (long *) malloc((pM->NrNzElts+1) * sizeof(long));
    
    if (pM->i == NULL || pM->j == NULL) {
        fprintf(stderr, "MMSparseMatrixAllocateMemory(): Not enough memory for matrix indices!\n");
        return FALSE;
    }
  
    /* Allocate memory for numerical values */
    if (pM->MMTypeCode[2] != 'P') { /* non-pattern matrix */
        pM->ReValue = (double *) malloc((pM->NrNzElts+1) * sizeof(double));
        if (pM->ReValue == NULL) {
            fprintf(stderr, "MMSparseMatrixAllocateMemory(): Not enough memory for matrix values!\n");
            return FALSE;
        }
    }
    if (pM->MMTypeCode[2] == 'C') { /* complex matrix */
        pM->ImValue = (double *) malloc((pM->NrNzElts+1) * sizeof(double));
        if (pM->ImValue == NULL) {
            fprintf(stderr, "MMSparseMatrixAllocateMemory(): Not enough memory for matrix values!\n");
            return FALSE;
        }
    }
  
    /* Allocate memory for Pstart, for distributed matrix */
    if (pM->MMTypeCode[0] == 'D') {
        pM->Pstart = (long *) malloc((pM->NrProcs+1) * sizeof(long));
        if (pM->Pstart == NULL) {
            fprintf(stderr, "MMSparseMatrixAllocateMemory(): Not enough memory for Pstart!\n");
            return FALSE;
        }
    }
  
    /* Allocate memory for row weights, for weighted matrix */
    if (pM->MMTypeCode[0] == 'W' && pM->NrRowWeights > 0) {
        pM->RowWeights = (long *) malloc(pM->NrRowWeights * sizeof(long));
        if (pM->RowWeights == NULL) {
            fprintf(stderr, "MMSparseMatrixAllocateMemory(): Not enough memory for RowWeights!\n");
            return FALSE;
        }
    }
  
    /* Allocate memory for column weights, for weighted matrix */
    if (pM->MMTypeCode[0] == 'W' && pM->NrColWeights > 0) {
        pM->ColWeights = (long *) malloc(pM->NrColWeights * sizeof(long));
        if (pM->ColWeights == NULL) {
            fprintf(stderr, "MMSparseMatrixAllocateMemory(): Not enough memory for ColWeights!\n");
            return FALSE;
        }
    }
    
    /* Allocate memory for row and column lambdas. */
    pM->RowLambda = (int *)malloc(pM->m * sizeof(int));
    if (pM->RowLambda == NULL) {
        fprintf(stderr, "MMSparseMatrixAllocateMemory(): Not enough memory for RowLambda!\n");
        return FALSE;
    }
    
    pM->ColLambda = (int *)malloc(pM->n * sizeof(int));
    if (pM->ColLambda == NULL) {
        fprintf(stderr, "MMSparseMatrixAllocateMemory(): Not enough memory for RowLambda!\n");
        return FALSE;
    }
    
    /* Allocate memory for row and column marks. */
    pM->RowMark = (int *) malloc(pM->m * sizeof(int));
    if (pM->RowMark == NULL) {
        fprintf(stderr, "MMSparseMatrixAllocateMemory(): Not enough memory for RowMark!\n");
        return FALSE;
    }
    
    pM->ColMark = (int *) malloc(pM->n * sizeof(int));
    if (pM->ColMark == NULL) {
        fprintf(stderr, "MMSparseMatrixAllocateMemory(): Not enough memory for ColMark!\n");
        return FALSE;
    }

    return TRUE;
} /* end MMSparseMatrixAllocateMemory */
  

int MMSparseMatrixFreeMemory(struct sparsematrix *pM) {

    /* This function frees memory for all the arrays of the matrix M:
           - the nonzero arrays (i, j, Val)
           - Pstart, which represents the start of the nonzeros
             of the processors, for a distributed matrix
           - the row and column weights, for a weighted matrix
           - dummies
           - row/column marks. */
    if (!pM) {
        fprintf(stderr, "MMSparseMatrixFreeMemory(): No matrix provided!\n");
        return FALSE;
    }

    if (pM->i != NULL)
        free(pM->i);
  
    if (pM->j != NULL)
        free(pM->j);
  
    if (pM->ReValue != NULL)
        free(pM->ReValue);
    
    if (pM->ImValue != NULL)
        free(pM->ImValue);
  
    if (pM->dummy != NULL)
        free(pM->dummy);
  
    if (pM->Pstart != NULL)
        free(pM->Pstart);
  
    if (pM->header != NULL)
        free(pM->header);
  
    if (pM->tail != NULL)
        free(pM->tail);

    if (pM->RowWeights != NULL)
        free(pM->RowWeights);
  
    if (pM->ColWeights != NULL)
        free(pM->ColWeights);
    
    if (pM->RowLambda != NULL)
        free(pM->RowLambda);

    if (pM->ColLambda != NULL)
        free(pM->ColLambda);

    if (pM->RowMark != NULL)
        free(pM->RowMark);
  
    if (pM->ColMark != NULL)
        free(pM->ColMark);
    
    if (pM->row_perm != NULL) free(pM->row_perm);
    if (pM->col_perm != NULL) free(pM->col_perm);
    if (pM->row_perm_inv != NULL) free(pM->row_perm_inv);
    if (pM->col_perm_inv != NULL) free(pM->col_perm_inv);

    if (pM->rowBoundaries != NULL) remembrance_destroy(pM->rowBoundaries);
    if (pM->colBoundaries != NULL) remembrance_destroy(pM->colBoundaries);
    
    return TRUE;
} /* end MMSparseMatrixFreeMemory */


int MMInsertProcessorIndices(struct sparsematrix *pM) {
    /* 
        This function writes the processor number to which each non-zero has been assigned into the real value of the corresponding non-zero.
    */
    int i, j;
    double *RePtr;
    
    if (!pM) {
        fprintf(stderr, "MMInsertProcessorIndices(): Null argument!\n");
        return FALSE;
    }
    
    if (pM->MMTypeCode[0] != 'D') {
        fprintf(stderr, "MMInsertProcessorIndices(): Warning, matrix has not been distributed yet!\n");
        /* return FALSE; */
    }
    
    /* Set matrix data type to integer. */
    pM->MMTypeCode[2] = 'I';
    
    if (pM->ImValue != NULL) {
        free(pM->ImValue);
        pM->ImValue = NULL;
    }
    
    if (pM->ReValue == NULL) {
        if ((pM->ReValue = (double *)malloc(pM->NrNzElts*sizeof(double))) == NULL) {
            fprintf(stderr, "MMInsertProcessorIndices(): Not enough memory!\n");
            return FALSE;
        }
    }
    
    RePtr = pM->ReValue;
    
    for (i = 0; i < pM->NrProcs; i++) {
        for (j = pM->Pstart[i]; j < pM->Pstart[i + 1]; j++) {
            *RePtr++ = (double)(i + 1);
        }
    }
    
    return TRUE;
}
 
  
int MMSparseMatrixReadHeader(struct sparsematrix *pM, FILE *fp) {

    /* This function reads the header from the input stream fp.
       The header is an optional set of lines starting with %,
       (the first is a banner line, and the remainder are comment lines),
       followed by a line with the size of the matrix,
           m n nz
       where m = the number of rows, m >=1
             n = the number of columns, n >= 1
             nz the number of nonzeros, nz >= 1.
       The size is read into variables M.m, M.n, M.NrNzElts.
       
       The banner line sets the type of the matrix.
       If no such line is present the type is set to the default. 
       The comment lines are read into the string M.header.

       For a distributed matrix, the size line has an extra parameter:
           m n nz P
       where P is the number of processors, P >= 1.

       For a weighted matrix, the size line is
           m n nz param
       where param = 0 means no row or column weights
                     1       only row weights
                     2       only column weights
                     3       both row and column weights.
       
       This function removes initial white space characters
       (blanks, tabs) from each line, and empty lines,
       before reading the actual input. */

    int first, last, SizeRead=FALSE, TypeSet=FALSE, count=0, param=0;
    char *line, linebuffer[MAX_LINE_LENGTH]; 
    long hsize;
    
    if (!pM || !fp) {
        fprintf(stderr, "MMSparseMatrixReadTail(): Null arguments!\n");
        return FALSE;
    }
  
    pM->header = NULL;
    hsize = 1;  /* header comment size (in characters) */
  
    while (SizeRead == FALSE && 
           (line = fgets(linebuffer, MAX_LINE_LENGTH, fp)) != NULL) { 
        /* a new line has been read */
        first = 0;
        last = strlen(line);
        
        /* Remove initial blank characters from the line */
        while (isspace((int) linebuffer[first])) 
            first++;
  
        line = linebuffer + first;
  
        if (first != last) {
            if (linebuffer[first] == '%') {
                /* The first part of the header, a line starting with % */
                if (TypeSet == FALSE) {
                    /* Line is a banner line. Set the type using the banner,
                       or set it to the default */
                    if (!MMSparseMatrixSetTypeCode(line, pM)) return FALSE;
                    TypeSet= TRUE;
                } else {
                    /* Add this line to the header comment */
                    hsize += last - first;
                    pM->header = (char *) realloc(pM->header, hsize*sizeof(char));
                    
                    if (pM->header == NULL) {
                        fprintf(stderr, "ReadHeader(): Not enough memory for realloc of matrix header!\n");
                        return FALSE;
                    }
                    
                    pM->header[hsize + first - last - 1] = '\0';
                    strncat(pM->header, line, last-first);
                }
            } else {
                /* The second part of the header, the size line */

                /* Read m, n, nz, and perhaps an extra parameter
                   from the line */
                count = sscanf(line, "%ld %ld %ld %d\n",
                               &(pM->m), &(pM->n), &(pM->NrNzElts), &param);
        
                if (count < 3 || pM->m < 1 || pM->n < 1 || pM->NrNzElts < 1) {
                    fprintf(stderr, "ReadHeader(): Error in matrix size!\n");
                    return FALSE;
                }
                else {
                    SizeRead = TRUE;
                }
            }
        }
    }
  
    /* If no type code has been set, set it to the default.
       This is needed if there are no lines starting with % */
    if (TypeSet == FALSE) 
        if (!MMSparseMatrixSetTypeCode("", pM)) return FALSE;
  
    if (pM->MMTypeCode[0] == 'D'){
        /* distributed matrix: param = number of processors */ 
        if (count == 4 && param >= 1) {
            pM->NrProcs = param;
        }
        else {
            fprintf(stderr, "ReadHeader(): Error in number of processors!\n");
            return FALSE;
        }
    } else if (pM->MMTypeCode[0] == 'W'){
        /* weighted matrix: param determines presence of row/column weights */

        if (count < 4 || param < 0 || param >= 4) {
            fprintf(stderr, "ReadHeader(): Error in weight param!\n");
            return FALSE;
        }

        if (param == 1 || param == 3)
            pM->NrRowWeights= pM->m;
        else
            pM->NrRowWeights= 0; 

        if (param == 2 || param == 3) 
            pM->NrColWeights= pM->n;
        else
            pM->NrColWeights= 0;
    }
    
    return TRUE;
} /* end MMSparseMatrixReadHeader */

  
int MMSparseMatrixReadPstart(struct sparsematrix *pM, FILE *fp) {

    /* This function reads the start of the nonzeros for each of the
       P processor parts, where P = M.NrProcs. The nonzeros of
       processor q start at position Pstart[q] and end at Pstart[q+1]-1,
       for 0 <= q < P.
       For convenience, Pstart[0] = 0 and Pstart[P] = nz(M). 
       This is checked, and if it does not hold the function aborts.
       Each Pstart[q] value should be on a separate line. */
       
    int t;
    
    if (!pM || !fp) {
        fprintf(stderr, "MMSparseMatrixReadPstart(): Null arguments!\n");
        return FALSE;
    }
  
    if (pM->MMTypeCode[0] == 'D'){
        for (t = 0; t <= pM->NrProcs; t++ )
            if (fscanf(fp, "%ld\n", &(pM->Pstart[t])) != 1) {
                fprintf(stderr, "MMSparseMatrixReadPstart(): Pstart: Read Error!\n");
                return FALSE;
            }

        if (pM->Pstart[0] != 0 || pM->Pstart[pM->NrProcs] != pM->NrNzElts) {
            fprintf(stderr, "MMSparseMatrixReadPstart(): Pstart: Read Error!\n");
            return FALSE;
        }
    }
    
    return TRUE;
} /* end MMSparseMatrixReadPstart */

  
int MMSparseMatrixReadWeights(struct sparsematrix *pM, FILE *fp) {

    /* This function reads the row and column weights of a weighted
       sparse matrix. The weights are integers. 
       
       The total number of row weights is M.NrRowWeights.
       The total number of column weights is M.NrColWeights. 
       Each weight should be on a separate line. */

    long t;
    
    if (!pM || !fp) {
        fprintf(stderr, "MMSparseMatrixReadWeights(): Null arguments!\n");
        return FALSE;
    }
  
    if (pM->MMTypeCode[0] == 'W'){
  
        for (t = 0; t < pM->NrRowWeights; t++ ) {
            if (fscanf(fp, "%ld\n", &(pM->RowWeights[t])) != 1) {
                fprintf(stderr, "MMSparseMatrixReadWeights(): Weights: Read Error!\n");
                return FALSE;
            }
        }
  
        for (t = 0; t < pM->NrColWeights; t++ ) {
            if (fscanf(fp, "%ld\n", &(pM->ColWeights[t])) != 1) {
                fprintf(stderr, "MMSparseMatrixReadWeights(): Weights: Read Error!\n");
                return FALSE;
            }
        }
    }

    return TRUE;
} /* end MMSparseMatrixReadWeights */
  

int MMSparseMatrixReadEntries(struct sparsematrix *pM, FILE *fp) {

    /* This function reads the nz(M) nonzero entries of a sparse matrix M
       in Matrix Market format from the input stream fp.

       An entry is:
       - a pair (i, j) for a pattern matrix,
         where i is the row index and j the column index,
         with 1 <= i <= m and 1 <= j <= n.
         Matrix Market indices are numbered starting from 1.
         The indices are converted to internal indices of the sparse
         matrix data structure, which are numbered starting from 0.
       - a triple (i, j, val) for a real or integer matrix,
         where val is the numerical value. This value is a real,
         or it is converted to a real.
       - a quadruple (i, j, ReVal, ImVal) for a complex matrix,
         where ReVal is the real part of the numerical value and
         ImVal the imaginary part. 
         
      Each entry should be on a separate line. */

    long t;
    
    if (!pM || !fp) {
        fprintf(stderr, "MMSparseMatrixReadEntries(): Null arguments!\n");
        return FALSE;
    }
  
    if (pM->MMTypeCode[2] == 'P') { /* pattern matrix */
        for (t = 0; t < pM->NrNzElts; t++) {
            if (fscanf(fp, "%ld %ld\n", &(pM->i[t]), &(pM->j[t])) != 2) {
                fprintf(stderr, "MMSparseMatrixReadEntries(): Pattern: Read Error!\n");
                return FALSE;
            }
            
            if (pM->i[t] < 1 || pM->i[t] > pM->m || pM->j[t] < 1 || pM->j[t] > pM->n) {
                fprintf(stderr, "MMSparseMatrixReadEntries(): Pattern: Read Error!\n");
                return FALSE;
            }
            
            pM->i[t]--; /* decrement by 1, since internally the
                            row and column indices of sparse matrices
                            are numbered from 0 */
            pM->j[t]--;
        }
    } else if (pM->MMTypeCode[2] == 'I') { /* integer matrix */
        for (t = 0; t < pM->NrNzElts; t++) {
            if (fscanf(fp, "%ld %ld %lg\n", &(pM->i[t]), &(pM->j[t]), &(pM->ReValue[t])) != 3) {
                fprintf(stderr, "MMSparseMatrixReadEntries(): Integer: Read Error!\n");
                return FALSE;
            }
            
            if (pM->i[t] < 1 || pM->i[t] > pM->m || pM->j[t] < 1 || pM->j[t] > pM->n) {
                fprintf(stderr, "MMSparseMatrixReadEntries(): Integer: Read Error!\n");
                return FALSE;
            }
            
            pM->i[t]--;
            pM->j[t]--;
        }
    } else if (pM->MMTypeCode[2] == 'R') { /* real matrix */
        for (t = 0; t < pM->NrNzElts; t++) {
            if (fscanf(fp, "%ld %ld %lg\n", &(pM->i[t]), &(pM->j[t]), &(pM->ReValue[t])) != 3) {
                fprintf(stderr, "MMSparseMatrixReadEntries(): Real: Read Error!\n"); 
                return FALSE;
            }
            
            if (pM->i[t] < 1 || pM->i[t] > pM->m || pM->j[t] < 1 || pM->j[t] > pM->n) {
                fprintf(stderr, "MMSparseMatrixReadEntries(): Real: Read Error!\n");
                return FALSE;
            }
            
            pM->i[t]--;
            pM->j[t]--;
        }
    } else if (pM->MMTypeCode[2] == 'C') { /* complex matrix */
        for (t = 0; t < pM->NrNzElts; t++) {
            if (fscanf(fp, "%ld %ld %lg %lg\n", &(pM->i[t]), &(pM->j[t]), &(pM->ReValue[t]), &(pM->ImValue[t])) != 4) {
                fprintf(stderr, "MMSparseMatrixReadEntries(): Complex: Read Error");
                return FALSE;
            }
            
            if (pM->i[t] < 1 || pM->i[t] > pM->m || pM->j[t] < 1 || pM->j[t] > pM->n) {
                fprintf(stderr, "MMSparseMatrixReadEntries(): Complex: Read Error");
                return FALSE;
            }
            
            pM->i[t]--;
            pM->j[t]--;
        }
    }
    
    /* Check if symmetric or hermitian matrix has entries
       only on or below the diagonal */
    if (pM->MMTypeCode[3] == 'S' || pM->MMTypeCode[3] == 'H') {
        for (t = 0; t < pM->NrNzElts; t++) {
            if (pM->i[t] < pM->j[t]) {
                fprintf(stderr, "MMSparseMatrixReadEntries(): Symmetric: Read Error!\n");
                return FALSE;
            }
        }
    }

     /* Check if skew-symmetric matrix has entries only below the diagonal */
    if (pM->MMTypeCode[3] == 'K') {
        for (t = 0; t < pM->NrNzElts; t++) {
            if (pM->i[t] <= pM->j[t]) {
                fprintf(stderr, "MMSparseMatrixReadEntries(): Skew-symmetric: Read Error");
                return FALSE;
            }
        }
    }
  
    return TRUE;
} /* end MMSparseMatrixReadEntries */
  

int MMSparseMatrixReadTail(struct sparsematrix *pM, FILE *fp) {
    
    /* This function reads the tail comment from the input stream fp
       into the string M.tail. It removes initial white space characters
       (blanks, tabs) from each line, and empty lines  */
    
    int first, last;
    char *line, linebuffer[MAX_LINE_LENGTH]; 
    long tsize;
    
    if (!pM || !fp) {
        fprintf(stderr, "MMSparseMatrixReadTail(): Null arguments!\n");
        return FALSE;
    }
  
    pM->tail = NULL;
    tsize = 1; /* tail size (in characters) */
    
    while ((line = fgets(linebuffer, MAX_LINE_LENGTH, fp)) != NULL) { 
        /* a new line has been read */
        first = 0; 
        last = strlen(line);
        
        if (linebuffer[last] != '\0') {
            fprintf(stderr, "MMSparseMatrixReadTail(): internal error!\n");
            return FALSE;
        }
  
        /* Remove initial blank characters from the line */
        while (isspace((int) linebuffer[first])) 
            first++;
  
        line = linebuffer + first;
        
        if (first != last) {
            /* Add this line to the tail */
            tsize += last - first;
            pM->tail = (char *) realloc(pM->tail, tsize * sizeof(char)); 
            
            if (pM->tail == NULL) {
                fprintf(stderr, "MMSparseMatrixReadTail(): Not enough memory for realloc of matrix tail!\n");
                return FALSE;
            }
            
            pM->tail[tsize + first - last - 1] = '\0';
            strncat(pM->tail, line, last-first);
        }
    }

    return TRUE;
} /* end MMSparseMatrixReadTail */
  
int MMSparseMatrixPrintHeader(struct sparsematrix *pM, FILE *stream, const char* name, const struct opts * pOptions) {

    /* This function prints the header of a sparse matrix M
       to the output stream.

       The header is a banner line with the types of the matrix,
       followed by the comment lines as stored in M.header,
       followed by a line with the size of the matrix,
           m n nz
       where m = the number of rows, m >=1
             n = the number of columns, n >= 1
             nz the number of nonzeros, nz >= 1.
       
       For a distributed matrix, the size line has an extra parameter:
           m n nz P
       where P is the number of processors, P >= 1.

       For a weighted matrix, the size line is
           m n nz param
       where param = 0 means no row or column weights
                     1       only row weights
                     2       only column weights
                     3       both row and column weights.
    */
       
    char line[MAX_LINE_LENGTH]; 
    int param;
    
    if (!pM || !stream) {
        fprintf(stderr, "MMSparseMatrixPrintHeader(): Null arguments!\n");
        return FALSE;
    }
  
    /* Print Matrix Market type codes to first line */
    if (!MMSparseMatrixGetTypeCode(pM, line, name, pOptions)) return FALSE;
    fprintf(stream, "%s", line);
  
    /* Print the header if name == NULL (indicating first header) */
    if (pM->header != NULL && name == NULL)
        fprintf(stream, "%s", pM->header);
  
    /* Print the comments if name == NULL (indicating first header) */
    if (name == NULL)
        if (!MMSparseMatrixPrintComment(pM, stream))
            return FALSE;
 
    /* Print m, n, nz, and perhaps an extra parameter to the first line */
    if (pM->MMTypeCode[0] == 'D') { /* distributed matrix */
        fprintf(stream, "%ld %ld %ld %d\n", pM->m, pM->n, pM->NrNzElts, pM->NrProcs);
    }
    else if (pM->MMTypeCode[0] == 'W') { /* weighted matrix */
        /* Compute parameter */
        param = 0;
        if (pM->NrRowWeights > 0)
            param |= 1;
        if (pM->NrColWeights > 0)
            param |= 2;
        fprintf(stream, "%ld %ld %ld %d\n", pM->m, pM->n, pM->NrNzElts,
                param); 
    } else {
        fprintf(stream, "%ld %ld %ld\n", pM->m, pM->n, pM->NrNzElts);
    }

    return TRUE;
} /* end MMSparseMatrixPrintHeader */
 
int MMSparseMatrixPrintComment(struct sparsematrix *pM, FILE *stream) {
    /* Simply prints the comment field of the sparse matrix to the file stream */
    if (!pM || !stream) {
        fprintf(stderr, "MMSparseMatrixPrintHeader(): Null arguments!\n");
        return FALSE;
    }
    if (pM->comment != NULL)
        fprintf(stream, "%s", pM->comment);
    return TRUE;
} /* end MMSparseMatrixPrintComment */

int MMSparseMatrixPrintPstart(struct sparsematrix *pM, FILE *stream, const struct opts * pOptions) {

    /* This function prints the start of the nonzeros for each of the
       P processor parts, where P = M.NrProcs. The nonzeros of processor q
       start at position Pstart[q] and end at Pstart[q+1]-1, for 0 <= q < P.
       For convenience, Pstart[0] = 0 and Pstart[P] = nz(M),
       but the function will print the values whatever they are.
       Each Pstart[q] value is printed on a separate line. */

    int t;
    
    if (!pM || !stream) {
        fprintf(stderr, "MMSparseMatrixPrintPstart(): Null arguments!\n");
        return FALSE;
    }
  
    if (pM->MMTypeCode[0] == 'D') 
        for (t = 0; t <= pM->NrProcs; t++ )
            if( pOptions->OutputFormat==OutputDMM )
                 fprintf(stream, "%ld\n", pM->Pstart[t]);
            else /* default to 1-based, except if explicitly in old Mondriaan mode */
                 fprintf(stream, "%ld\n", pM->Pstart[t]+1);

    return TRUE;
} /* end MMSparseMatrixPrintPstart */

  
int MMSparseMatrixPrintWeights(struct sparsematrix *pM, FILE *stream) {

    /* This function prints the row and column weights of a weighted
       sparse matrix. The weights are integers.

       The total number of row weights is M.NrRowWeights.
       The total number of column weights is M.NrColWeights.
       Each weight is printed on a separate line. */

    long t;
    
    if (!pM || !stream) {
        fprintf(stderr, "MMSparseMatrixPrintPstart(): Null arguments!\n");
        return FALSE;
    }
  
    if (pM->MMTypeCode[0] == 'W'){
        for (t = 0; t < pM->NrRowWeights; t++ )
            fprintf(stream, "%ld\n", pM->RowWeights[t]);
        for (t = 0; t < pM->NrColWeights; t++ )
            fprintf(stream, "%ld\n", pM->ColWeights[t]);
    }

    return TRUE;
} /* end MMSparseMatrixPrintWeights */
  
  
int MMSparseMatrixPrintEntries(struct sparsematrix *pM, FILE *stream) {

    /* This function prints the nz(M) nonzero entries of a sparse matrix M
       in Matrix Market format to the output stream.

       An entry is:
       - a pair (i, j) for a pattern matrix,
         where i is the row index and j the column index,
         with 1 <= i <= m and 1 <= j <= n.
         Matrix Market indices are numbered starting from 1.
         The indices are converted from internal indices of the sparse
         matrix data structure, which are numbered starting from 0.
       - a triple (i, j, val) for a real or integer matrix,
         where val is the numerical value.
       - a quadruple (i, j, ReVal, ImVal) for a complex matrix,
         where ReVal is the real part of the numerical value and
         ImVal the imaginary part.

      Each entry will be printed on a separate line. */

    long t, r;
    
    if (!pM || !stream) {
        fprintf(stderr, "MMSparseMatrixPrintEntries(): Null arguments!\n");
        return FALSE;
    }
  
    if (pM->MMTypeCode[2] == 'P') { /* pattern matrix */
        for (t = 0; t < pM->NrNzElts; t++)
            fprintf(stream, "%ld %ld\n", pM->i[t]+1, pM->j[t]+1);
    } else if (pM->MMTypeCode[2] == 'I') { /* integer matrix */
        for (t = 0; t < pM->NrNzElts; t++) {
            /* round to the nearest integer */
            if (pM->ReValue[t] >= 0)
                r = pM->ReValue[t]+0.5; 
            else 
                r = pM->ReValue[t]-0.5;
            fprintf(stream, "%ld %ld %ld\n", pM->i[t]+1, pM->j[t]+1, r);
        }
    } else if (pM->MMTypeCode[2] == 'R') { /* real matrix */
        for (t = 0; t < pM->NrNzElts; t++)
            fprintf(stream, "%ld %ld %g\n", pM->i[t]+1, pM->j[t]+1, pM->ReValue[t]);
    } else if (pM->MMTypeCode[2] == 'C') { /* complex matrix */
        for (t = 0; t < pM->NrNzElts; t++)
            fprintf(stream, "%ld %ld %g %g\n", pM->i[t]+1, pM->j[t]+1, 
                    pM->ReValue[t], pM->ImValue[t]);
    }

    return TRUE;
} /* end MMSparseMatrixPrintEntries */

  
int MMSparseMatrixPrintTail(struct sparsematrix *pM, FILE *stream) {

    /* This function prints the tail comment of a square matrix 
       to the output stream */ 
    if (!pM || !stream) {
        fprintf(stderr, "MMSparseMatrixPrintTail(): Null arguments!\n");
        return FALSE;
    }

    if (pM->tail != NULL)
        fprintf(stream, "%s", pM->tail);

    return TRUE;
} /* end MMSparseMatrixPrintTail */
  
int AddDummiesToSparseMatrix(struct sparsematrix *pM)  {
  
    /* This function adds dummies to a square matrix to make its diagonal 
       structurally nonzero. It also creates a boolean array M.dummy
       such that M.dummy[i] = TRUE if entry (i,i) is a dummy nonzero.
       
       A dummy is an entry in the sparse matrix data structure with
       numerical value 0.
       
       This function should only be used before matrix partitioning,
       because it does not update distribution information. */
  
    long t, tt;
    long i, j, n, NrNzEltsNew;
    
    if (!pM) {
        fprintf(stderr, "AddDummiesToSparseMatrix(): Null arguments!\n");
        return FALSE;
    }
  
    if (pM->m != pM->n) {
        fprintf(stderr, "AddDummiesToSparseMatrix(): Matrix is not square!\n");
        return FALSE;
    }
  
    n = pM->n;
    pM->NrDummies = n; /* maximum number of dummies needed */
    
    /* Allocate memory for dummy information */
    pM->dummy = (int *) malloc(n * sizeof(int));
    
    if (pM->dummy == NULL) {
        fprintf(stderr, "AddDummiesToSparseMatrix(): Not enough memory for dummies!\n");
        return FALSE;
    }
  
    /* Initialise dummy information */
    for (i = 0; i < n; i++)
        pM->dummy[i] = TRUE; /* dummy needed for matrix element (i,i) */
  
    /* Fill dummy and count the number of dummies */
    for (t = 0; t < pM->NrNzElts; t++) {
        i = pM->i[t];
        j = pM->j[t];
        if (i == j)
            if (pM->dummy[i]) {
                pM->dummy[i] = FALSE;
                pM->NrDummies--; /* no dummy needed because of actual
                                     diagonal nonzero */
            }
    }
  
    /* Allocate memory for dummies. Allocate one element extra
       in case the matrix becomes empty. */
    NrNzEltsNew = pM->NrNzElts + pM->NrDummies; 

    pM->i = (long *) realloc(pM->i, (NrNzEltsNew+1) * sizeof(long));
    pM->j = (long *) realloc(pM->j, (NrNzEltsNew+1) * sizeof(long));
    
    if (pM->i == NULL || pM->j == NULL) {
        fprintf(stderr, "AddDummiesToSparseMatrix(): Not enough memory for dummies!\n");
        return FALSE;
    }
    
    if (pM->MMTypeCode[2] != 'P') {
        pM->ReValue = (double *) realloc(pM->ReValue, (NrNzEltsNew+1) * sizeof(double));
        
        if (pM->ReValue == NULL) {
            fprintf(stderr, "AddDummiesToSparseMatrix(): Not enough memory for dummies!\n");
            return FALSE;
        }
    }
    
    if (pM->MMTypeCode[2] == 'C') {
        pM->ImValue = (double *) realloc(pM->ImValue, (NrNzEltsNew+1) * sizeof(double));
        
        if (pM->ImValue == NULL) {
            fprintf(stderr, "AddDummiesToSparseMatrix(): Not enough memory for dummies!\n");
            return FALSE;
        }
    }
  
    /* Add the dummies to the matrix */
    tt = pM->NrNzElts;
    for (i = 0; i < n; i++) {
        if (pM->dummy[i]) {
            pM->i[tt] = i;
            pM->j[tt] = i;
            if (pM->MMTypeCode[2] != 'P')
                pM->ReValue[tt] = 0.0;
            if (pM->MMTypeCode[2] == 'C')
                pM->ImValue[tt] = 0.0; 
            tt++;
        }
    }
  
    /* Update the number of nonempty elements of the matrix */
    pM->NrNzElts = NrNzEltsNew;
    
    return TRUE;
} /* end AddDummiesToSparseMatrix */
  
int RemoveDummiesFromSparseMatrix(struct sparsematrix *pM)  {
  
    /* This function removes all dummies from a square matrix
       and frees the memory of the boolean array M.dummy.

       A dummy is an entry in the sparse matrix data structure with
       numerical value 0.

       This function can be used before and after partitioning,
       because it updates the distribution information by adjusting Pstart. */
  
    long i, j, t, tt;
    int q;
    long NrRemoved, start;
    
    if (!pM) {
        fprintf(stderr, "RemoveDummiesFromSparseMatrix(): Null parameter!\n");
        return FALSE;
    }
    
    /* If the matrix is not square, there should be no dummies ! */
    if (pM->m != pM->n) {
        fprintf(stderr, "RemoveDummiesFromSparseMatrix(): Matrix is not square");
        return FALSE;
    }
    
    if (pM->NrProcs < 1 || pM->Pstart == NULL) { /* matrix not partitioned */
        NrRemoved = 0;
        t = 0;
        while (t < pM->NrNzElts) {
            while ((i = pM->i[t]) == pM->j[t] && NrRemoved < pM->NrDummies && pM->dummy[i]) { 
                tt = pM->NrNzElts-1;
    
                /* Swap dummy t with last entry tt */
                SwapLong(pM->i, t, tt); /* pM->i[t] <-> pM->i[tt] */
                SwapLong(pM->j, t, tt);
                if (pM->MMTypeCode[2] != 'P') 
                    SwapDouble(pM->ReValue, t, tt); 
                if (pM->MMTypeCode[2] == 'C') 
                    SwapDouble(pM->ImValue, t, tt);
      
                /* Remove last element */
                pM->NrNzElts--;
                NrRemoved++;

                /* entry t could again be a dummy, hence we need a loop */
            }
            t++;
        }
    } else { /* Update Pstart as well */
        NrRemoved = 0;
        start = pM->Pstart[0];
        tt = 0;
        for (q = 0; q < pM->NrProcs; q++) {
            for (t = start; t < pM->Pstart[q+1]; t++) {
                i = pM->i[t];
                j = pM->j[t];
                if (i == j && pM->dummy[i])
                    NrRemoved++; /* remove dummy */
                else { /* Move nondummy to front position tt */
                    pM->i[tt] = i;
                    pM->j[tt] = j;
                    if (pM->MMTypeCode[2] != 'P')
                        pM->ReValue[tt] = pM->ReValue[t];
                    if (pM->MMTypeCode[2] == 'C')
                        pM->ImValue[tt] = pM->ImValue[t];
                    tt++;
                }
            }
    
            /* Update Pstart */
            start = pM->Pstart[q+1];
            pM->Pstart[q+1] = tt;
        }
  
        /* Update NrNzElts */
        pM->NrNzElts -= NrRemoved;
    }
  
    /* Check the number of removed dummies */
    if (pM->NrDummies != NrRemoved) {
        fprintf(stderr, "RemoveDummiesFromSparseMatrix(): Internal error: NrDummies != NrRemoved!\n");
        return FALSE;
    }
  
    /* Free the memory of the dummies */
    pM->i = (long *) realloc(pM->i, (pM->NrNzElts+1) * sizeof(long));
    pM->j = (long *) realloc(pM->j, (pM->NrNzElts+1) * sizeof(long));
    
    if (pM->i == NULL || pM->j == NULL) {
        fprintf(stderr, "RemoveDummiesFromSparseMatrix(): Not enough memory!\n");
        return FALSE;
    }
    
    if (pM->MMTypeCode[2] != 'P') {
        pM->ReValue = (double *) realloc(pM->ReValue, (pM->NrNzElts+1)*sizeof(double));
        
        if (pM->ReValue == NULL) {
            fprintf(stderr, "RemoveDummiesFromSparseMatrix(): Not enough memory!\n");
            return FALSE;
        }
    }
    
    if (pM->MMTypeCode[2] == 'C') {
        pM->ImValue = (double *) realloc(pM->ImValue, (pM->NrNzElts+1)*sizeof(double));
        
        if (pM->ImValue == NULL) {
            fprintf(stderr, "RemoveDummiesFromSparseMatrix(): Not enough memory!\n");
            return FALSE;
        }
    } 
  
    if (pM->dummy != NULL) {
        free(pM->dummy);
        pM->dummy = NULL;
    }
    
    pM->NrDummies = 0;

    return TRUE;
} /* end RemoveDummiesFromSparseMatrix */

  
int SparseMatrixSymmetric2Full(struct sparsematrix *pM) {

    /* This function converts symmetric, skew-symmetric,
       and hermitian matrices to the full sparse representation.

       A symmetric matrix (a[i,j] = a[j,i]) or hermitian matrix
       (a[i,j] = conjg(a[j,i])) stores only the nonzeros on or 
       below the diagonal.

       A skew-symmetric (or anti-symmetric) matrix (a[i,j] = -a[j,i])
       stores only the nonzeros below the diagonal.

       This representation is also called the lower triangular format.
       In the full sparse representation all nonzeros are stored, 
       including those above the diagonal. 

       This function can be used before and after partitioning.
       If it is called after partitioning, it assumes
       symmetric partitioning by assigning a[i,j] and a[j,i]
       to the same processor.

    */
  
    long t, tt, Ndiag;
    long start, NrNzEltsNew;
    int q;
    
    if (!pM) {
        fprintf(stderr, "SparseMatrixSymmetric2Full(): Null parameter!\n");
        return FALSE;
    }
  
    /* Check if matrix is symmetric */
    if (pM->m != pM->n || 
         (pM->MMTypeCode[3] != 'S' && pM->MMTypeCode[3] != 'K' &&
          pM->MMTypeCode[3] != 'H')) {
        fprintf(stderr, "SparseMatrixSymmetric2Full(): matrix is not symmetric!\n");
        return FALSE;
    }
  
    /* Count the number of diagonal elements */
    Ndiag = 0;
    for (t = 0; t < pM->NrNzElts; t++)
        if (pM->i[t] == pM->j[t])
            Ndiag++;
    NrNzEltsNew = 2 * pM->NrNzElts - Ndiag;

    if (pM->MMTypeCode[3] == 'K' && Ndiag > 0) {
        fprintf(stderr, "SparseMatrixSymmetric2Full(): matrix is not skew-symmetric!\n");
        return FALSE;
    }

    /* Update Pstart */
    if (pM->NrProcs >= 1 && pM->Pstart != NULL) {
        tt = 0;
        start = pM->Pstart[0];
        for (q = 0; q < pM->NrProcs; q++) {
            for (t = start; t < pM->Pstart[q+1]; t++) {
                if (pM->i[t] == pM->j[t])
                    tt++;
                else
                    tt += 2;
            }    
           
            start = pM->Pstart[q+1];
            pM->Pstart[q+1] = tt;
        }
    }
   
    /* Allocate memory for the entries */
    pM->i = (long *) realloc(pM->i, (NrNzEltsNew+1) * sizeof(long));
    pM->j = (long *) realloc(pM->j, (NrNzEltsNew+1) * sizeof(long));
    
    if (pM->i == NULL || pM->j == NULL) {
        fprintf(stderr, "SparseMatrixSymmetric2Full(): Not enough memory to store full matrix!\n");
        return FALSE;
    }
    
    if (pM->MMTypeCode[2] != 'P') {
        pM->ReValue = (double *) realloc(pM->ReValue, (NrNzEltsNew+1)*sizeof(double));
        
        if (pM->ReValue == NULL) {
            fprintf(stderr, "SparseMatrixSymmetric2Full(): Not enough memory to store full matrix!\n");
            return FALSE;
        }
    }
    
    if (pM->MMTypeCode[2] == 'C') {
        pM->ImValue = (double *) realloc(pM->ImValue, (NrNzEltsNew+1)*sizeof(double));
        
        if (pM->ImValue == NULL) {
            fprintf(stderr, "SparseMatrixSymmetric2Full(): Not enough memory to store full matrix!\n");
            return FALSE;
        }
    }
  
    /* Copy the entries */
    tt = NrNzEltsNew -1;
    for (t = pM->NrNzElts - 1; t >= 0; t--) {
        /* Copy entry t into tt */
        if (pM->i[t] <  pM->j[t]) {
            fprintf(stderr, "SparseMatrixSymmetric2Full(): entry above diagonal!\n");
            return FALSE;
        }

        pM->i[tt] = pM->i[t];
        pM->j[tt] = pM->j[t];
        if (pM->MMTypeCode[2] != 'P')
            pM->ReValue[tt] = pM->ReValue[t];
        if (pM->MMTypeCode[2] == 'C')
            pM->ImValue[tt] = pM->ImValue[t];
        tt--;   
  
        /* Copy and transpose nondiagonal entry t */
        if (pM->i[tt+1] != pM->j[tt+1]) {
            pM->i[tt] = pM->j[tt+1]; 
            pM->j[tt] = pM->i[tt+1]; 
            if (pM->MMTypeCode[2] != 'P') {
                if (pM->MMTypeCode[3] == 'K')
                    pM->ReValue[tt] = -pM->ReValue[tt+1]; 
                else
                    pM->ReValue[tt] = pM->ReValue[tt+1]; 
            }
            if (pM->MMTypeCode[2] == 'C') { 
                if (pM->MMTypeCode[3] == 'K' || pM->MMTypeCode[3] == 'H') 
                    pM->ImValue[tt] = -pM->ImValue[tt+1];
                else 
                    pM->ImValue[tt] = pM->ImValue[tt+1];
            } 
            tt--;
        }
    }
  
    /* Update the number of nonzero entries and the matrix type code */
    pM->NrNzElts = NrNzEltsNew;
    pM->MMTypeCode[3] = 'G';
    strcpy(pM->Symmetry, "general");
   
    return TRUE;
} /* end SparseMatrixSymmetric2Full */
  

int SparseMatrixSymmetricLower2Random(struct sparsematrix *pM) {

    /* This function converts a symmetric matrix stored
       in lower triangular format (with only nonzeros a[i,j] with
       i <= j present) into a symmetric matrix stored in random format,
       where the numerical value a[i,j] with i <= j is assigned
       either to entry (i,j) or to (j,i) in the sparse matrix
       data structure, each with equal probability.

       This conversion is done by transposing the indices of every nonzero
       with a probability of 50 percent, without changing the numerical value.
       This function preserves the distribution information.
    */

    long t, tmp;
    
    if (!pM) {
        fprintf(stderr, "SparseMatrixSymmetricLower2Random(): Null parameter!\n");
        return FALSE;
    }
  
    /* Check if matrix is symmetric */
    if (pM->m != pM->n || 
         (pM->MMTypeCode[3] != 'S'  && pM->MMTypeCode[3] != 'K'  &&
          pM->MMTypeCode[3] != 'H')) {
        fprintf(stderr, "SparseMatrixSymmetricLower2Random(): matrix is not symmetric!\n");
        return FALSE;
    }
    
    for (t = 0; t < pM->NrNzElts; t++) {
        if (Random1(0,1) == 0) {
            /* Swap index i and j */
            tmp = pM->i[t];
            pM->i[t] = pM->j[t];
            pM->j[t] = tmp;
        }
    }
    
    return TRUE;
} /* end SparseMatrixSymmetricLower2Random */

int SparseMatrixSymmetricRandom2Lower(struct sparsematrix *pM) {

    /* This function converts a symmetric matrix stored
       in random format, where the numerical value a[i,j] with i <= j
       is assigned either to entry (i,j) or to (j,i) in the sparse matrix
       data structure, not necesssarily with equal probability,
       into a symmetric matrix stored in lower triangular format
       (with only nonzeros a[i,j] with i <= j present).

       This conversion is done by transposing the indices of every nonzero
       (i,j) with i < j in the data structure.

       This function preserves the distribution information.
    */


    long t, tmp;
  
    if (!pM) {
        fprintf(stderr, "SparseMatrixSymmetricRandom2Lower(): Null parameter!\n");
        return FALSE;
    }

    /* Check if matrix is symmetric */
    if (pM->m != pM->n || 
         (pM->MMTypeCode[3] != 'S'  && pM->MMTypeCode[3] != 'K'  &&
          pM->MMTypeCode[3] != 'H')) {
        fprintf(stderr, "SparseMatrixSymmetricRandom2Lower(): matrix is not symmetric!\n");
        return FALSE;
    }
   
    for (t = 0; t < pM->NrNzElts; t++) {
        if (pM->i[t] < pM->j[t]) {
            /* Swap index i and j */
            tmp = pM->i[t];
            pM->i[t] = pM->j[t];
            pM->j[t] = tmp;
        }
    }
    
    return TRUE;
} /* end SparseMatrixSymmetricRandom2Lower */

int SparseMatrixRemoveDuplicates(struct sparsematrix *pM) {

    /* This function removes the duplicate nonzeros from a sparse matrix
       by adding their numerical values, giving a single value at a location (i,j).

       This function is intended for checking the input of a sparse matrix after it has 
       been read from a file, and repairing it if necessary. 

       This function can only be used before partitioning and before adding dummies,
       as it does not preserve the distribution information nor the dummy information.
    */

    long nz, t, tt, k, kpred, *K;
    
    if (!pM) {
        fprintf(stderr, "SparseMatrixRemoveDuplicates(): Null parameter!\n");
        return FALSE;
    }
    
    nz = pM->NrNzElts;
    
    /* Allocate memory for nonzero index array */   
    K = (long *) malloc(nz * sizeof(long));
    
    if (K == NULL) {
        fprintf(stderr, "SparseMatrixRemoveDuplicates(): Not enough memory!\n");
        return FALSE;
    }
    
    for (t = 0; t < nz; t++)
        K[t] = t;
        
    /* Sort nonzero indices stably, first by column index, then by row index.
       As a result, the indices K[t] give the order of the nonzeros
       by increasing row index, and within the same row by increasing column index */

    CSort(K, pM->j, pM->n-1, 0, nz-1);
    CSort(K, pM->i, pM->m-1, 0, nz-1);

    /* Check if a nonzero is at the same location as its predecessor.
       If so, add the numerical values and mark the predecessor for removal */
    for (t = 1; t < nz; t++) {
        k = K[t];
        kpred = K[t-1];
        if (pM->i[k] == pM->i[kpred] && pM->j[k] == pM->j[kpred]) {
            /* add the numerical values */
            if (pM->MMTypeCode[2] != 'P') 
                pM->ReValue[k] +=  pM->ReValue[kpred];
            if (pM->MMTypeCode[2] == 'C')
                pM->ImValue[k] +=  pM->ImValue[kpred];
        
            /* mark the predecessor for removal */
            pM->i[kpred] = -1;
        }
    }
    
    /* Remove the marked nonzeros */
    tt = 0;
    for (t = 0; t < nz; t++) {
        if (pM->i[t] != -1){
            /* copy nonzero from position t into tt */
            pM->i[tt] = pM->i[t];
            pM->j[tt] = pM->j[t]; 
            if (pM->MMTypeCode[2] != 'P')
                pM->ReValue[tt] = pM->ReValue[t];
            if (pM->MMTypeCode[2] == 'C')
                pM->ImValue[tt] = pM->ImValue[t];
            tt++;
        }
    }
    
#ifdef INFO
    printf("  %ld duplicate nonzeros removed\n", nz-tt);
    printf("  New number of nonzeros: %ld \n", tt);
#endif

    pM->NrNzElts = tt;
    
    free(K); 
    
    return TRUE;
} /* end SparseMatrixRemoveDuplicates */


int SparseMatrixStructurallySymmetric(struct sparsematrix *pM) {

    /* This function checks whether a sparse matrix is 
       structurally symmetric, i.e, it is square and 
       A[i,j] is a nonzero (in the data structure) if and only if A[j,i] is a nonzero,
       for all i,j. If this holds, the value TRUE is returned; otherwise, FALSE. 
       
       The numerical value of the nonzero A[i,j] does not matter.
       It can even be 0.0 (an accidental zero).     

       This function does not change the matrix.
    */

    long nz, t, kr, kc, *Kr, *Kc;
    int symmetric;

    if (!pM) return FALSE;
    
    /* Check whether the matrix is square */
    if (pM->m != pM->n)
        return FALSE;

    /* Check if matrix is structurally symmetric by definition */
    if ((pM->MMTypeCode[3] == 'S' || pM->MMTypeCode[3] == 'K' ||
          pM->MMTypeCode[3] == 'H'))
              return TRUE;
        
    nz = pM->NrNzElts;
   
    /* Allocate memory for nonzero index arrays Kr and Kc */   
    Kr = (long *) malloc(nz * sizeof(long));
    Kc = (long *) malloc(nz * sizeof(long));
    
    if (Kr == NULL || Kc == NULL) return FALSE;
    
    for (t = 0; t < nz; t++) {
        Kr[t] = t;        
        Kc[t] = t;
    }

    /* Sort nonzero indices of Kr stably, first by column index, then by row index.
       As a result, the indices Kr[t] give the order of the nonzeros
       by increasing row index, and within the same row by increasing column index */
 
    CSort(Kr, pM->j, pM->n-1, 0, nz-1);
    CSort(Kr, pM->i, pM->n-1, 0, nz-1);
    
    /* Sort nonzero indices of Kc stably, first by row index, then by column index.
       As a result, the indices Kc[t] give the order
       by increasing column index, and within the same column by increasing row index */

    CSort(Kc, pM->i, pM->n-1, 0, nz-1);
    CSort(Kc, pM->j, pM->n-1, 0, nz-1);

    /* Check if each nonzero (i,j) also has a nonzero in position (j,i) */
    symmetric = TRUE;
    
    for (t = 0; t < nz; t++) {
        kr = Kr[t];
        kc = Kc[t];
        if (pM->i[kr] != pM->j[kc] || pM->j[kr] != pM->i[kc]) {
            symmetric = FALSE;
            break;
        }
    }
    
    free(Kc); 
    free(Kr); 
    
    return symmetric;
    
} /* end  SparseMatrixStructurallySymmetric */

int SparseMatrixMarkCutRows(struct sparsematrix *pM) {
    /* This function marks all rows of the given matrix that are cut
       with respect to the current processor assignments. */
    long i, j;
    
    if (!pM) {
        fprintf(stderr, "SparseMatrixMarkCutRows(): Null argument!\n");
        return FALSE;
    }
    
    if (!pM->RowMark || !pM->Pstart || !pM->i) {
        fprintf(stderr, "SparseMatrixMarkCutRows(): Faulty matrix!\n");
        return FALSE;
    }
    
    /* First clear all rows. */
    for (i = 0; i < pM->m; i++) {
        pM->RowMark[i] = -1;
    }
    
    for (i = 0; i < pM->NrProcs; i++) {
        for (j = pM->Pstart[i]; j < pM->Pstart[i + 1]; j++) {
            if (pM->RowMark[pM->i[j]] == -1) {
                /* If this row is unassigned, assign it to the current processor. */
                pM->RowMark[pM->i[j]] = i;
            }
            else if (pM->RowMark[pM->i[j]] != i) {
                /* This row was assigned to a different processor, hence it is cut. */
                pM->RowMark[pM->i[j]] = -2;
            }
        }
    }
    
    pM->NumRowsMarked = 0;
    
    for (i = 0; i < pM->m; i++) {
        if (pM->RowMark[i] == -2) {
            pM->NumRowsMarked++;
            pM->RowMark[i] = TRUE;
        }
        else {
            pM->RowMark[i] = FALSE;
        }
    }
    
#ifdef INFO2
    printf("  Marked %ld of %ld rows as cut.\n", pM->NumRowsMarked, pM->m);
#endif
    
    return TRUE;
}

int SparseMatrixMarkCutColumns(struct sparsematrix *pM) {
    /* See MMSparseMatrixMarkCutRows(). */
    long i, j;
    
    if (!pM) {
        fprintf(stderr, "SparseMatrixMarkCutColumns(): Null argument!\n");
        return FALSE;
    }
    
    if (!pM->ColMark || !pM->Pstart || !pM->j) {
        fprintf(stderr, "SparseMatrixMarkCutColumns(): Faulty matrix!\n");
        return FALSE;
    }
    
    for (i = 0; i < pM->n; i++) {
        pM->ColMark[i] = -1;
    }
    
    for (i = 0; i < pM->NrProcs; i++) {
        for (j = pM->Pstart[i]; j < pM->Pstart[i + 1]; j++) {
            if (pM->ColMark[pM->j[j]] == -1) {
                pM->ColMark[pM->j[j]] = i;
            }
            else if (pM->ColMark[pM->j[j]] != i) {
                pM->ColMark[pM->j[j]] = -2;
            }
        }
    }
    
    pM->NumColsMarked = 0;
    
    for (i = 0; i < pM->n; i++) {
        if (pM->ColMark[i] == -2) {
            pM->NumColsMarked++;
            pM->ColMark[i] = TRUE;
        }
        else {
            pM->ColMark[i] = FALSE;
        }
    }
    
#ifdef INFO2
    printf("  Marked %ld of %ld columns as cut.\n", pM->NumColsMarked, pM->n);
#endif
    
    return TRUE;
}

int SparseMatrixOriginal2Local(struct sparsematrix *pM, long int **row_perms, long int **col_perms) {
    int i, j;
    long int *row, *col;
    long int rowC, colC;
    long int minR, minC;
    if(!pM || !row_perms || !col_perms) {
        fprintf( stderr, "SparseMatrixOriginal2Local(): null pointer arguments supplied!\n" );
        return FALSE;
    }
    if(pM->MMTypeCode[0]!='D') {
        fprintf( stderr, "SparseMatrixOriginal2Local() only functions on distributed matrices!\n" );
        return FALSE;
    }
    if(pM->ViewType==ViewTypeLocal) {
        fprintf( stderr, "SparseMatrixOriginal2Local(): given matrix is already in local view!\n" );
        return FALSE;
    }
    /* allocate temporary values */
    row = malloc( pM->m * sizeof( long int ) );
    col = malloc( pM->n * sizeof( long int ) );
    /* allocate size arrays */
    if( row_perms != NULL )
        row_perms[pM->NrProcs] = malloc( pM->NrProcs * sizeof( long int ) );
    if( col_perms != NULL )
        col_perms[pM->NrProcs] = malloc( pM->NrProcs * sizeof( long int ) );
    for(i=0; i<pM->NrProcs; i++) {
        /* initialise temp vars */
        for(j=0; j<pM->m; j++) row[j]=-1; /* untouched */
        for(j=0; j<pM->n; j++) col[j]=-1;
        rowC = colC = 0;
        /* find local row/cols */
        for(j=pM->Pstart[i]; j<pM->Pstart[i+1]; j++) {
            if( row[pM->i[j]] == -1 ) {
                row[pM->i[j]] = -2; /* touched */
                rowC++;
            }
            if( col[pM->j[j]] == -1 ) {
                col[pM->j[j]] = -2;
                colC++;
            }
        }
        minR = minC = 0;
        /* create STABLE translation arrays */
        for( j = 0; j < pM->m; j++ )
            if( row[ j ] == -2 ) { /* for each filled row */
                row[ j ] = minR++; /* permute to the first available index */
            }
        for( j = 0; j < pM->n; j++ )
            if( col[ j ] == -2 ) { /* for each filled column..(see above) */
                col[ j ] = minC++;
            }
        
        /* do the translation */
        for(j=pM->Pstart[i]; j<pM->Pstart[i+1]; j++) {
            pM->i[j] = row[ pM->i[j] ];
            pM->j[j] = col[ pM->j[j] ];
        }
        /* build the local2global indices */
        if( row_perms != NULL )
            row_perms[i] = malloc( rowC * sizeof( long int ) );
        if( col_perms != NULL )
            col_perms[i] = malloc( colC * sizeof( long int ) );
        if( row_perms != NULL )
            for(j=0; j<pM->m; j++)
                if( row[j] != -1 )
                    row_perms[i][row[j]] = j;
        if( col_perms != NULL )
            for(j=0; j<pM->n; j++)
                if( col[j] != -1 ) 
                    col_perms[i][col[j]] = j;
        /* store size */
        if( row_perms != NULL )
            row_perms[pM->NrProcs][i] = rowC;
        if( col_perms != NULL )
            col_perms[pM->NrProcs][i] = colC;
    }
    /* set matrix to local view */
    pM->ViewType = ViewTypeLocal;

    /* free temporary vars */
    free(row);
    free(col);

    /* return success */
    return TRUE;
}

int SparseMatrixLocal2Vector(struct sparsematrix *pM, long int **local2glob, long int *vec_distr,
                             long int **out_count, long int ***out_local2proc, long int ***out_local2index, 
                             const char dir ) {
    /* dir = 0: row direction, dir != 0: column direction */
    int i, j;
    long int *count = malloc( pM->NrProcs * sizeof( long int ) );
    const int length = dir ? pM->n : pM->m;
    long int *inv_glob     = malloc( length * sizeof( long int ) );
    long int **local2proc  = malloc( pM->NrProcs * sizeof( long int * ) );
    long int **local2index = malloc( pM->NrProcs * sizeof( long int * ) );
    if( !count || !local2proc || !local2index || !inv_glob ) {
        fprintf(stderr, "SparseMatrixLocal2Vector: error during memory allocation\n");
        return FALSE;
    }
    for( i=0; i<pM->NrProcs; i++ ) count[i] = 0;
    for( i=0; i<length; i++ ) 
       inv_glob[i] = count[ vec_distr[i] ]++;
    for( i=0; i<pM->NrProcs; i++ ) {
        local2proc [i] = malloc( local2glob[pM->NrProcs][i] * sizeof( long int ) );
        local2index[i] = malloc( local2glob[pM->NrProcs][i] * sizeof( long int ) );
        if( !local2proc || !local2index ) {
            fprintf(stderr, "SparseMatrixLocal2Vector: error during memory allocation\n");
            return FALSE;
        }
        for( j=0; j<local2glob[pM->NrProcs][i]; j++ ) {
           local2proc [i][j] = vec_distr[ local2glob[i][j] ];
           local2index[i][j] = inv_glob [ local2glob[i][j] ];
        }
    }
    *out_count = count; *out_local2proc = local2proc; *out_local2index = local2index;
    free( inv_glob );
    return TRUE;
}

int CreateInitialMediumGrainDistribution(struct sparsematrix *pM, long *mid){

    /* This function creates a twodimensional distribution of the nonzeros
       of the sparse matrix M by splitting them into two sets, with indices
       in the ranges 0..mid-1 and mid..nz(M)-1, respectively.
 
       The first group contains nonzeros of those rows that have been assigned
       completely to a processor, and the second group nonzeros of columns
       that have been assigned completely to a processor.

       Furthermore, the first group contains the remaining
       nonzeros a(i,j) for which c(j) < r(i),
       i.e. the number of nonzeros in column j is less than the number in row i.
       The second group contains the rest. Ties are broken based on
       the dimensions of the matrix.
    */

    long *rowCnt, *colCnt;
    long nz, t, i, j, dir, left, right, nneR=0,nneC=0;
    double tmpRe = 0.0, tmpIm = 0.0 ;

    long *rowR = (long *)malloc(pM->m*sizeof(long));
    long *rowC = (long *)malloc(pM->m*sizeof(long));
    long *colR = (long *)malloc(pM->n*sizeof(long));
    long *colC = (long *)malloc(pM->n*sizeof(long));

    nz = pM->NrNzElts;
    if(nz==0){
        *mid=0;
        return TRUE;
    }

    rowCnt = (long *)malloc(pM->m*sizeof(long));
    colCnt = (long *)malloc(pM->n*sizeof(long));
    if (rowCnt == NULL || colCnt == NULL){
        fprintf(stderr, "CreateInitialMediumGrainDistribution(): Not enough memory!\n");
        return FALSE;
    }

    /* Count the nonzeros in the rows and columns */
    for(i=0; i<pM->m; i++){
        rowCnt[i] = 0;
        rowR[i] = 0;
        rowC[i] = 0;
    }
    for(j=0; j<pM->n; j++){
        colCnt[j] = 0;
        colR[j] = 0;
        colC[j] = 0;
    }

    for(t=0; t<nz; t++){
        if(rowCnt[pM->i[t]]==0) nneR++;
        rowCnt[pM->i[t]]++;
        if(colCnt[pM->j[t]]==0) nneC++;
        colCnt[pM->j[t]]++;
    }

    /* Choose uniform direction for ties between row and column count */
    if (nneC > nneR)
        dir = COL; /* columns are shorter on average */
    else if (nneC < nneR)
        dir = ROW;
    else if (Random1(0,1) == 0) /* random tie-breaking */
        dir = ROW;
    else
        dir = COL;
        
    for(t=0;t<nz;t++){
        i =   pM->i[t];
        j =   pM->j[t];
        if(rowCnt[i]==1 || (colCnt[j]!=1 && (rowCnt[i] > colCnt[j] || (rowCnt[i] == colCnt[j] && dir == COL )))){
            rowC[i]++;
            colC[j]++;
        }else{
            rowR[i]++;
            colR[j]++;
        }
    }

    /* Split nonzeros into two sets, 0..mid-1 and mid..nz-1 */ 
    left = 0;  
    right = nz - 1 ;
    while(left <= right){ /* nonzeros in range left..right still to be split */

        /* Copy left entry values into temporary space */
        i =   pM->i[left];
        j =   pM->j[left];
        if (pM->MMTypeCode[2] != 'P')
            tmpRe = pM->ReValue[left];
        if (pM->MMTypeCode[2] == 'C')
            tmpIm = pM->ImValue[left];

        if((colR[j]==1 && colCnt[j]!=1) || rowCnt[i]==1 ||  (colCnt[j] != 1 && rowC[i]!=1 &&
                            (rowCnt[i] > colCnt[j] || (rowCnt[i] == colCnt[j] &&  
                             dir == COL )))){
            left++; /* column is shorter and wins */
        } else {
            /* Copy right values into left */
            pM->i[left] =   pM->i[right];
            pM->j[left] =   pM->j[right];
            if (pM->MMTypeCode[2] != 'P')
                pM->ReValue[left] = pM->ReValue[right];
            if (pM->MMTypeCode[2] == 'C')
                pM->ImValue[left] = pM->ImValue[right];

            /* Copy temporary values into right */
            pM->i[right] = i;
            pM->j[right] = j;
            if (pM->MMTypeCode[2] != 'P')
                pM->ReValue[right] = tmpRe;
            if (pM->MMTypeCode[2] == 'C')
                pM->ImValue[right] = tmpIm; 
            right-- ;
        }
    }

    *mid = left;

    free(colCnt);
    free(rowCnt);
    free(rowR);
    free(rowC);
    free(colR);
    free(colC);

    return TRUE;

} /* end CreateInitialMediumGrainDistribution */

