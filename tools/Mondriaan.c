#include <Mondriaan.h>

int main(int argc, char **argv) {

    /* This function is the main function, called mondriaan. It partitions
       a sparse matrix, a dense input vector, and a dense output vector
       for the purpose of parallel sparse matrix-vector multiplication u = A*v.

       The function tries to minimise the total communication volume
       by suitably partitioning the matrix, while keeping the computation
       load imbalance withing a user-specified fraction epsilon.
       For the resulting communication volume, it tries to balance
       the communication work among the processors by suitably
       partitioning the vectors. The partitioning of the input and output
       vectors can either be done independently, or by imposing
       distr(u) = distr(v).

       The sparse matrix is read from file using the Matrix Market (MM) format.
       Its partitioning is written to file using a variant of the MM format,
       where all nonzeros belonging to the same processor are listed together.

       If the input file is called
           pi.mtx
       and 4 processors are required, the following output files are produced:
           pi.mtx-P4 containing the matrix distribution
           pi.mtx-v4                input vector distribution
           pi.mtx-u4                output vector distribution
           pi.mtx-C4                Cartesian submatrices corresponding
                                    to the matrix distribution
       Furthermore, useful information on the partitioning process
       and the result, such as the load imbalance and communication
       volume obtained, is written to standard output.

       The sparse matrix partitioner can also be used for partitioning
       hypergraphs. To do this, Mondriaan must be used in 1D column mode. 
       Matrix columns then represent hypergraph vertices; they can be given
       integer weights. Matrix rows then represent hyperedges; they cannot
       be given weights.

       The following additional operations are carried out by the main program:
           - timing the matrix and vector partitionings
           - removing duplicate nonzeros
           - (optionally) checking whether the matrix is structurally symmetric
           - (optionally) expanding a symmetric matrix stored in lower
             triangular format to a fully stored matrix
           - (optionally) adding diagonal dummy nonzeros to a square matrix
             to make the diagonal completely nonzero (thereby facilitating
             the vector distribution with distr(u) = distr(v)).

       All partitioning options and parameters can be set using
       the file mondriaan.defaults (provided the main function is compiled
       with the name mondriaan.) If no such file exists, it is created with
       sensible defaults.

       If the input matrix has already been partitioned, it will be
       partitioned again, from scratch.
 */

    struct opts Options; /* The Mondriaan options */
    struct sparsematrix A; /* The Matrix;-) */

    long MinNrNzElts, MaxNrNzElts; /* minimum, maximum number of nonzero elements 
                                       per processor */
    double AvgNrNzElts; /* average number of nonzero matrix elements per processor */

    long MinWeight, MaxWeight, TotWeight; /* minimum, maximum, total column weight */
    double AvgWeight;                    /* average column weight */
    long weight, i, j;
    unsigned long int ui;
    long int *temp;
    long ComU, ComV; /* communication cost for vectors u, v  */
    long q, nzq;     /* processor number and its corresponding number of nonzeros */
    int weighted, symmetric; /* boolean registering whether the input matrix
                                 was weighted, symmetric */

    long int *u_proc, *v_proc;   /* vector distribution for u, v */
    char output[MAX_WORD_LENGTH]; /* filename of the output */
    char EMMcomment[2000]; /* Comments in EMM one-file output */
    
    long int **row_perms = NULL, **col_perms = NULL; /* Will store local to global index info */

    /* Below is used for Extended Matrix-Market one-file output */
    long int *inv_count, **row_local2proc, **row_local2index, **col_local2proc, **col_local2index;

    FILE *File = NULL;
    char stdinName[] = "stdin";
    
    /* Timing variables */
    clock_t starttime, endtime;
    double cputime;
#ifdef UNIX
    struct timeval starttime1, endtime1;
#endif

    /* Get the parameters from the command line and initialise Options */
    SetDefaultOptions(&Options);
    
    if (!SetOptionsFromFile(&Options, "Mondriaan.defaults")) {
        fprintf(stderr, "main(): warning, cannot set options from 'Mondriaan.defaults', using default options and creating standard 'Mondriaan.defaults'!\n");
        
        File = fopen("Mondriaan.defaults", "w");
        
        if (File != NULL) {
            ExportDefaultOptions(File);
            SetDefaultOptions(&Options);
            fclose(File);
        }
        else {
            fprintf(stderr, "main(): Unable to create 'Mondriaan.defaults'!\n");
        }
    }
    
    if (!GetParameters(&Options, argc, argv)) {
        fprintf(stderr, "main(): invalid command line parameters!\n");
        exit(-1);
    }
    
    if (!ApplyOptions(&Options)) {
        fprintf(stderr, "main(): could not apply given options!\n");
        exit(-1);
    }
    
    /* Read matrix file from disk or standard input. */
    if (!strcmp(Options.matrix, "-") || !strcmp(Options.matrix, "stdin")) {
        if (!MMReadSparseMatrix(stdin, &A)) {
            fprintf(stderr, "main(): Could not read matrix from standard input!\n");
            exit(-1);
        }
        
        Options.matrix = stdinName;
    }
    else {
        File = fopen(Options.matrix, "r");
    
        if (!File) {
            fprintf(stderr, "main(): Unable to open '%s' for reading!\n", Options.matrix);
            exit(-1);
        }
        
        if (!MMReadSparseMatrix(File, &A)) {
            fprintf(stderr, "main(): Could not read matrix!\n");
            exit(-1);
        }
        
        fclose(File);
    }
  
#ifdef INFO
    printf("\nUsing Mondriaan version %s.\n\n", MONDRIAANVERSION);
    printf("\n**************************************************************\n");
    printf("Problem statistics:\n");
    printf("  Matrix:           : %s\n",Options.matrix);
    printf("  %s %s %s %s %s\n", A.Banner,A.Object,A.Format,A.Field,A.Symmetry);
    printf("  m = Nr rows       : %ld\n", A.m); 
    printf("  n = Nr columns    : %ld\n", A.n); 
    printf("  nz = Nr nonzeros  : %ld\n", A.NrNzElts);
  
    /* Check if matrix A is structurally symmetric */
    if (SparseMatrixStructurallySymmetric(&A))
        printf("  Matrix is structurally symmetric\n\n");
    else
        printf("  Matrix is structurally unsymmetric\n\n");
#endif

    /* Remove duplicate nonzeros by adding them */
    if (!SparseMatrixRemoveDuplicates(&A)) {
        fprintf(stderr, "main(): Unable to remove duplicates!\n");
        exit(-1);
    }

    /* Check if matrix A is already distributed */
    if (A.MMTypeCode[0] == 'D') {
        /* Matrix will be partitioned again */
        fprintf(stderr, "Warning: Matrix '%s' already distributed !\n", 
                Options.matrix);
        fprintf(stderr, "         (Ignoring current partitions)\n"); 
        
        A.NrProcs = 0;
        if (A.Pstart != NULL)
            free(A.Pstart);
        A.Pstart = NULL;
    }
   
    /* Check if matrix is weighted (thus representing a hypergraph).
       In that case, it must have n column weights (representing vertex weights), 
       0 row weights and the split strategy must be onedimcol */

    if (A.MMTypeCode[0] == 'W' && A.NrColWeights != A.n) {
        fprintf(stderr, "main(): Weighted matrix with NrColWeights != n!\n");
        exit(-1);
    }

    if (A.MMTypeCode[0] == 'W' && A.NrRowWeights != 0) {
        fprintf(stderr, "Warning: Matrix '%s' has row weights!\n",
                Options.matrix);
        fprintf(stderr, "         Row-weighted column partitioning not yet implemented\n");
        fprintf(stderr, "         (Ignoring row weights)\n");
                
        A.NrRowWeights = 0;
        if (A.RowWeights != NULL)
            free(A.RowWeights);
        A.RowWeights = NULL;
    }

    if (A.MMTypeCode[0] == 'W' && Options.SplitStrategy != OneDimCol) {
        fprintf(stderr, "Warning: Matrix '%s' is a weighted matrix!\n", Options.matrix);
        fprintf(stderr, "         must be partitioned by onedimcol strategy\n");
        fprintf(stderr, "         (Ignoring requested split strategy)\n");
        Options.SplitStrategy = OneDimCol;
    }

    /* Register whether the input matrix was weighted, since the object type code will
       be changed by the partitioning procedure to the code 'D' for a distributed matrix */
    if (A.MMTypeCode[0] == 'W')
        weighted = TRUE;
    else
        weighted = FALSE;
  
    /* Register whether the input matrix was symmetric, since the symmetry type code will
       be changed by the conversion to full, to the code 'G' for a general matrix */
    if (A.m == A.n && 
         (A.MMTypeCode[3]=='S' || A.MMTypeCode[3]=='K' || A.MMTypeCode[3]=='H')) {
        symmetric = TRUE; 
    } else {
        symmetric = FALSE;
        if (Options.SplitStrategy == SFineGrain) {
            fprintf(stderr, "Error: Symmetric finegrain can only be used on symmetric input matrices!\n");
            exit(-1);
        }
    }
    
    if (symmetric) {
        if (Options.SymmetricMatrix_UseSingleEntry == SingleEntNo)
            SparseMatrixSymmetric2Full(&A); 
        else if (Options.SplitStrategy == SFineGrain)
            SparseMatrixSymmetricRandom2Lower(&A);
        else if (Options.SymmetricMatrix_SingleEntryType == ETypeRandom)
            SparseMatrixSymmetricLower2Random(&A);
    }

    if (Options.SplitStrategy == SFineGrain && Options.SymmetricMatrix_SingleEntryType == ETypeRandom)
        printf("Warning: Symmetric finegrain requires lower-triangular format of symmetric matrix;\n         Random single entry type option is overridden.\n");
  
    /* If the matrix is square, add the dummies if requested.
       This may lead to an enhanced vector distribution in the case of
       an equal distribution of the input and output vectors.  */
    if (A.m == A.n && 
         Options.SquareMatrix_DistributeVectorsEqual == EqVecYes &&
         Options.SquareMatrix_DistributeVectorsEqual_AddDummies == DumYes)
        AddDummiesToSparseMatrix(&A);
  
    /* Set the number of processors */
    A.NrProcs = Options.P;
    
#ifdef INFO
    printf("  P = Nr processors : %d\n", A.NrProcs);
    printf("  Nr dummies added  : %ld\n", A.NrDummies);
    printf("  allowed imbalance : %g %c\n", 100*Options.eps, '%');
    printf("\n****** Mondriaan matrix distribution ******\n\n");
#endif

    /* Initialise Pstart with all nonzeros in processor 0 */
    A.Pstart = (long *) malloc((A.NrProcs+1) * sizeof(long));
    if (A.Pstart == NULL) {
        fprintf(stderr, "main(): Not enough memory for Pstart!\n");
        exit(-1);
    }
    
    A.Pstart[0] = 0;
    for (q = 1; q <= A.NrProcs; q++)
        A.Pstart[q] = A.NrNzElts;
  
    /**** Distribute the matrix (and time it) ****/
    starttime = clock();
#ifdef UNIX
    gettimeofday(&starttime1, NULL);
#endif
  
    if (!DistributeMatrixMondriaan(&A, A.NrProcs, Options.eps, &Options, 0)) {
        fprintf(stderr, "main(): Unable to distribute matrix!\n");
        exit(-1);
    }
    
    endtime = clock();
#ifdef UNIX
    gettimeofday(&endtime1, NULL);
#endif
    cputime = ((double) (endtime - starttime)) / CLOCKS_PER_SEC;
  
#ifdef TIME
    printf("  matrix distribution CPU-time    : %f seconds\n", cputime); 
#ifdef UNIX 
    printf("  matrix distribution elapsed time: %f seconds\n",
             (endtime1.tv_sec - starttime1.tv_sec) +
             (endtime1.tv_usec - starttime1.tv_usec) / 1000000.0); 
#endif 
    fflush(stdout);
#endif 
  
    /* Remove the dummies */
    if (A.m == A.n &&
         Options.SquareMatrix_DistributeVectorsEqual == EqVecYes && 
         Options.SquareMatrix_DistributeVectorsEqual_AddDummies == DumYes)
        RemoveDummiesFromSparseMatrix(&A);
  
    /* Print information about the number of matrix elements */
    AvgNrNzElts = (double) A.NrNzElts / A.NrProcs;
    MinNrNzElts = LONG_MAX;
    MaxNrNzElts = LONG_MIN;
    for (q = 0; q < A.NrProcs; q++) {
        nzq = A.Pstart[q+1] - A.Pstart[q];
        if (nzq < MinNrNzElts)
            MinNrNzElts = nzq;
        if (nzq > MaxNrNzElts)
            MaxNrNzElts = nzq;
    }

#ifdef INFO
    printf("  Nr nonzero matrix elements:\n");
    printf("    tot         : %ld \n", A.NrNzElts);
    printf("    avg = tot/P : %g  \n", AvgNrNzElts);
    printf("    max         : %ld \n", MaxNrNzElts);
    printf("    min         : %ld \n", MinNrNzElts);
    printf("    imbalance   : %g %c\n",
                100 * (MaxNrNzElts/AvgNrNzElts - 1.0), '%');
#endif

    if (weighted){
        A.MMTypeCode[0] = 'W'; /* temporarily, needed for computing weights */
        TotWeight = ComputeWeight(&A, 0, A.NrNzElts-1, NULL, &Options);
        AvgWeight = (double) TotWeight / A.NrProcs;
        MinWeight = LONG_MAX;
        MaxWeight = LONG_MIN;
        for (q = 0; q < A.NrProcs; q++) {
            weight = ComputeWeight(&A, A.Pstart[q], A.Pstart[q+1]-1, NULL, &Options);
            if (weight < MinWeight)
                MinWeight = weight;
            if (weight > MaxWeight)
                MaxWeight = weight;
        }
        A.MMTypeCode[0] = 'D'; /* back to distributed matrix */

#ifdef INFO
        printf("  Vertex (column) weight:\n");
        printf("    tot         : %ld \n", TotWeight);
        printf("    avg = tot/P : %g  \n", AvgWeight);
        printf("    max         : %ld \n", MaxWeight);
        printf("    min         : %ld \n", MinWeight);
        printf("    imbalance   : %g %c\n",
                    100 * (MaxWeight/AvgWeight - 1.0), '%');
#endif
    }

    /* Convert randomly represented symmetric matrix to standard
       lower triangular form */
    if (symmetric &&
         Options.SymmetricMatrix_UseSingleEntry == SingleEntYes &&
         Options.SymmetricMatrix_SingleEntryType == ETypeRandom)
              SparseMatrixSymmetricRandom2Lower(&A);
 
    /* Writing (distributed) matrix info */
    if( Options.OutputMode == MultipleFiles ) {
 
        /* Write the distributed matrix to file */
        sprintf(output, "%s-P%d", Options.matrix, A.NrProcs);
        File = fopen(output, "w");
        if (!File) fprintf(stderr, "main(): Unable to open '%s' for writing!\n", output);
        else {
            MMWriteSparseMatrix(&A, File, NULL, &Options);
            fclose(File);
        }

        /* Commit to permutation */
        if (A.row_perm_inv != NULL && A.col_perm_inv != NULL) {
            for( i=0; i<A.NrNzElts; i++ ) {
                A.i[ i ] = A.row_perm_inv[ A.i[ i ] ];
                A.j[ i ] = A.col_perm_inv[ A.j[ i ] ];
            }
        }
 
        /* Write permuted matrix */
        sprintf(output, "%s-reor-P%d", Options.matrix, A.NrProcs);
        File = fopen(output, "w");
        A.MMTypeCode[0] = 'M';
        if (!File) fprintf(stderr, "main(): Unable to open '%s' for writing!\n", output);
        else {
            MMWriteSparseMatrix(&A, File, NULL, &Options);
            fclose(File);
        }
        A.MMTypeCode[0] = 'D';

        /* Write out block information */
        if (A.rowBoundaries) {
            if (remembrance_get( A.rowBoundaries )==NULL) {
                fprintf(stderr, "main(): Error during read-out of row boundaries\n");
                exit( -1 );
            }
            sprintf(output, "%s-rowblocks%d", Options.matrix, A.NrProcs);
            File = fopen(output, "w");
            if (!File) fprintf(stderr, "main(): Unable to open '%s' for writing!\n", output);
            else {
                remembrance_write(A.rowBoundaries, File);
                fclose(File);
            }
        }
        if (Options.SplitStrategy != SFineGrain)
            if (A.colBoundaries) {
                if (remembrance_get( A.colBoundaries )==NULL) {
                    fprintf(stderr, "main(): Error during read-out of column boundaries\n");
                    exit( -1 );
                }
                sprintf(output, "%s-colblocks%d", Options.matrix, A.NrProcs);
                File = fopen(output, "w");
                if (!File) fprintf(stderr, "main(): Unable to open '%s' for writing!\n", output);
                else {
                   remembrance_write(A.colBoundaries, File);
                   fclose(File);
                }
            }
        

        /* Undo commit to permutation */
        if (A.row_perm_inv != NULL && A.col_perm_inv != NULL) {
            for( i=0; i<A.NrNzElts; i++ ) {
                A.i[ i ] = A.row_perm[ A.i[ i ] ];
                A.j[ i ] = A.col_perm[ A.j[ i ] ];
            }
        }

        /* Write out local matrix format */
        row_perms = malloc( (A.NrProcs+1) * sizeof( long int * ) );
        col_perms = malloc( (A.NrProcs+1) * sizeof( long int * ) );
        if( !row_perms || !col_perms || 
            !SparseMatrixOriginal2Local(&A, row_perms, col_perms) ) {
            fprintf(stderr, "main(): Unable to transform to local view!");
            exit(-1);
        }
        sprintf(output, "%s-local%d", Options.matrix, A.NrProcs);
        File = fopen(output, "w");
        if (!File) fprintf(stderr, "main(): Unable to open '%s' for writing!\n", output);
        else {
            MMWriteSparseMatrix(&A, File, NULL, &Options);
            fclose(File);
        }

        /* Transform matrix back */
        for( i=0; i<A.NrProcs; i++ )
            for(j=A.Pstart[i]; j<A.Pstart[i+1]; j++ ) {
                A.i[j] = row_perms[i][ A.i[j] ];
                A.j[j] = col_perms[i][ A.j[j] ];
            }
        A.ViewType = ViewTypeOriginal;

        /* Write out matrix the entries of which are processor indices. */
        if (!MMInsertProcessorIndices(&A)) {
            fprintf(stderr, "main(): Unable to write processor indices!\n");
            exit(-1);
        }
    
        A.MMTypeCode[0] = 'M';
        sprintf(output, "%s-I%d", Options.matrix, A.NrProcs);
    
        File = fopen(output, "w");
    
        if (!File) fprintf(stderr, "main(): Unable to open '%s' for writing!\n", output);
        else {
            MMWriteSparseMatrix(&A, File, NULL, &Options);
            fclose(File);
        }
        A.MMTypeCode[0] = 'D';
 
        /* Write out permutations. */
        if (A.col_perm != NULL) {
            sprintf(output, "%s-col%d", Options.matrix, A.NrProcs);
            File = fopen(output, "w");
            if (!File) fprintf(stderr, "main(): Unable to open '%s' for writing!\n", output);
            else {
                WriteVector(A.col_perm, 0, NULL, A.n, 0, File, &Options);
               fclose(File);
            }
        }
    
        if (A.row_perm != NULL) {
            sprintf(output, "%s-row%d", Options.matrix, A.NrProcs);
            File = fopen(output, "w");
            if (!File) fprintf(stderr, "main(): Unable to open '%s' for writing!\n", output);
            else {
                WriteVector(A.row_perm, 0, NULL, A.m, 0, File, &Options);
                fclose(File);
            }
        }

        /* Free local to global index */
        for( i=0; i<A.NrProcs+1; i++ ) {
            free( row_perms[i] );
            free( col_perms[i] );
        }
        free( row_perms ); free( col_perms );
    } else if (Options.OutputMode == OneFile) {
        /* Extended MatrixMarket output */
        if (Options.OutputFormat != OutputEMM) {
            fprintf(stderr, "main(): Single-file output requires EMM format!\n");
            fprintf(stderr, "        overriding output format to EMM.\n");
            Options.OutputFormat = OutputEMM;
        }
        sprintf(EMMcomment, "%% File generated by running Mondriaan on A=%s with p=%d and \\epsilon=%f\n", Options.matrix, A.NrProcs, Options.eps);
        sprintf(EMMcomment, "%s%% Main matrix corresponds to the distributed version of A\n", EMMcomment);
        if (A.row_perm_inv != NULL && A.col_perm_inv != NULL) {
            switch( Options.OrderPermutation ) {
                case OrderPrefix:
                    sprintf(EMMcomment, "%s%% This is followed by the (not-distributed) reverse Bordered Block Diagonal permuted matrix 'PAQ',\n", EMMcomment);
                    break;
                case OrderPostfix:
                    sprintf(EMMcomment, "%s%% This is followed by the (not-distributed) Bordered Block Diagonal permuted matrix 'PAQ',\n", EMMcomment);
                    break;
                default:
                    sprintf(EMMcomment, "%s%% This is followed by the (not-distributed) Separated Block Diagonal permuted matrix 'PAQ',\n", EMMcomment);
            }
            sprintf(EMMcomment, "%s%% followed by the row-block boundary vector 'Row-boundaries',\n", EMMcomment);
            sprintf(EMMcomment, "%s%% followed by the row-block hierarchy vector 'Row-hierarchy',\n", EMMcomment);
            sprintf(EMMcomment, "%s%% followed by the column-block boundary vector 'Column-boundaries',\n", EMMcomment);
            sprintf(EMMcomment, "%s%% followed by the column-block hierarchy vector 'Column-hierarchy',\n", EMMcomment);
        }
        sprintf(EMMcomment, "%s%% followed by a local view of the distributed form of A 'Local-A',\n", EMMcomment);
        sprintf(EMMcomment, "%s%% followed by a mapping of local row indices to global row indices 'LocalRow2Global',\n", EMMcomment);
        sprintf(EMMcomment, "%s%% followed by a mapping of local column indices to global column indices 'LocalCol2Global',\n", EMMcomment);
        if (A.row_perm != NULL)
            sprintf(EMMcomment, "%s%% followed by the row permutation vector corresponding to P in PAQ 'Row-permutation',\n", EMMcomment);
        if (A.col_perm != NULL)
            sprintf(EMMcomment, "%s%% followed by the column permutation vector corresponding to Q in PAQ 'Column-permutation',\n", EMMcomment);
        sprintf(EMMcomment, "%s%% followed by the input vector distribution 'Input-vector',\n", EMMcomment);
        sprintf(EMMcomment, "%s%% followed by the output vector distribution 'Output-vector',\n", EMMcomment);
        sprintf(EMMcomment, "%s%% followed by the local lengths of the output vectors 'OutputVectorLengths',\n", EMMcomment);
        sprintf(EMMcomment, "%s%% followed by the matrix row index to output vector processor ID mapping 'LocalRow2Processor',\n", EMMcomment);
        sprintf(EMMcomment, "%s%% followed by the matrix row index to output vector index mapping 'LocalRow2Index',\n", EMMcomment);
        sprintf(EMMcomment, "%s%% followed by the matrix column index to input vector processor ID mapping 'LocalCol2Processor',\n", EMMcomment);
        sprintf(EMMcomment, "%s%% followed by the matrix column index to input vector index mapping 'LocalCol2Index'.\n", EMMcomment);
        A.comment = EMMcomment;
        sprintf(output, "%s-P%d.emm", Options.matrix, A.NrProcs);
        File = fopen(output, "w");
        if (!File) {
            fprintf(stderr, "main(): Unable to open '%s' for writing!\n", output);
            exit( -1 );
        }
        MMWriteSparseMatrix(&A, File, NULL, &Options);

        /* Commit to permutation */
        if (A.row_perm_inv != NULL && A.col_perm_inv != NULL) {
            for( i=0; i<A.NrNzElts; i++ ) {
                A.i[ i ] = A.row_perm_inv[ A.i[ i ] ];
                A.j[ i ] = A.col_perm_inv[ A.j[ i ] ];
            }
        }
 
        /* Write permuted matrix */
        A.MMTypeCode[0] = 'M';
        MMWriteSparseMatrix(&A, File, "PAQ", &Options);
        A.MMTypeCode[0] = 'D';

        /* Write out block information */
        if (A.rowBoundaries) {
            if (remembrance_get( A.rowBoundaries )==NULL) {
                fprintf(stderr, "main(): Error during read-out of row boundaries\n");
                exit( -1 );
            }
            temp = malloc(A.rowBoundaries->size * sizeof( long int ));
            for(ui=0; ui<A.rowBoundaries->size; ui++)
                temp[ ui ] = A.rowBoundaries->vector[ ui ].index;
            WriteVector(temp, 0, "Row-boundaries", A.rowBoundaries->size, 0, File, &Options);
            for(ui=0; ui<A.rowBoundaries->size-1; ui++ )
                temp[ ui ] = A.rowBoundaries->vector[ ui ].parent == ULONG_MAX ? 0 : A.rowBoundaries->vector[ ui ].parent + 1;
            WriteVector(temp, 1, "Row-hierarchy", A.rowBoundaries->size-1, 0, File, &Options);
            free(temp);
        }
        if (Options.SplitStrategy != SFineGrain) {
            if (A.colBoundaries) {
                if (remembrance_get( A.colBoundaries )==NULL) {
                    fprintf(stderr, "main(): Error during read-out of column boundaries\n");
                    exit( -1 );
                }
                temp = malloc( A.colBoundaries->size * sizeof( long int ) );
                for(ui=0; ui<A.colBoundaries->size; ui++)
                    temp[ ui ] = A.colBoundaries->vector[ ui ].index;
                WriteVector(temp, 0, "Col-boundaries", A.colBoundaries->size, 0, File, &Options);
                for(ui=0; ui<A.colBoundaries->size-1; ui++ )
                    temp[ ui ] = A.colBoundaries->vector[ ui ].parent == ULONG_MAX ? 0 : A.colBoundaries->vector[ ui ].parent + 1;
                WriteVector(temp, 1, "Col-hierarchy", A.colBoundaries->size-1, 0, File, &Options);
                free(temp);
            }
        }
    
        /* Write out local matrix format */
        row_perms = malloc( (A.NrProcs+1) * sizeof( long int * ) );
        col_perms = malloc( (A.NrProcs+1) * sizeof( long int * ) );
        if( !row_perms || !col_perms || 
            !SparseMatrixOriginal2Local(&A, row_perms, col_perms) ) {
            fprintf(stderr, "main(): Unable to transform to local view!");
            exit(-1);
        }
        MMWriteSparseMatrix(&A, File, "Local-A", &Options);

        /* Transform matrix back from Local to global view */
        for( i=0; i<A.NrProcs; i++ )
            for(j=A.Pstart[i]; j<A.Pstart[i+1]; j++ ) {
                A.i[j] = row_perms[i][ A.i[j] ];
                A.j[j] = col_perms[i][ A.j[j] ];
            }
        A.ViewType = ViewTypeOriginal;

        /* Permutate row_perms and col_perms so that they point
         * to global indices (corresponding to A) instead of
         * permuted indices (corresponding to PAQ). 
        for( i=0; i<A.NrProcs; i++ ) {
                for( j=0; j<row_perms[A.NrProcs][i]; j++ )
                        row_perms[i][j] = row_perms[i][j]; 
                for( j=0; j<col_perms[A.NrProcs][i]; j++ )
                        col_perms[i][j] = col_perms[i][j];
        } */

        WriteVectorCollection(row_perms, "LocalRow2Global", A.NrProcs, row_perms[A.NrProcs], File);
        WriteVectorCollection(col_perms, "LocalCol2Global", A.NrProcs, col_perms[A.NrProcs], File);

        /* Write out matrix the entries of which are processor indices */
        /*if (!MMInsertProcessorIndices(&A)) {
            fprintf(stderr, "main(): Unable to write processor indices!\n");
            exit(-1);
        }
        A.MMTypeCode[0] = 'M';
        tempchar = A.MMTypeCode[2];
        A.MMTypeCode[2] = 'I';
        MMWriteSparseMatrix(&A, File, "Global-A", &Options);
        A.MMTypeCode[0] = 'D';
        A.MMTypeCode[2] = tempchar;*/
 
        /* Write out permutations */
        if (A.row_perm != NULL) {
            WriteVector(A.row_perm, 0, "Row-permutation", A.m, 0, File, &Options);
        }
        if (A.col_perm != NULL) {
            WriteVector(A.col_perm, 0, "Column-permutation", A.n, 0, File, &Options);
        }
    } else if (Options.OutputMode == DIMACS) {
        /* For DIMACS we only need the vector distribution. */
    } else {
        fprintf(stderr, "main(): Unknown output mode!\n" );
        exit( -1 );
    }

    /* Convert symmetrically partitioned, symmetric matrix to full form 
       for vector distribution */
    if (symmetric &&
         Options.SymmetricMatrix_UseSingleEntry == SingleEntYes) {
        /* Some of the above transformations can make A arbitrary symmetric
           instead of lower-traingular */
        SparseMatrixSymmetricRandom2Lower(&A);
        SparseMatrixSymmetric2Full(&A); /* now A.MMTypeCode[3]='G' */
        
        /* Print information about the number of matrix elements */
        AvgNrNzElts = (double) A.NrNzElts / A.NrProcs;
        MinNrNzElts = LONG_MAX;
        MaxNrNzElts = LONG_MIN;
        for (q = 0; q < A.NrProcs; q++) {
             nzq = A.Pstart[q+1] - A.Pstart[q];
             if (nzq < MinNrNzElts)
                 MinNrNzElts = nzq;
             if (nzq > MaxNrNzElts)
                 MaxNrNzElts = nzq;
        }

#ifdef INFO
        printf("  After conversion of distributed symmetric matrix to"
               "  distributed full matrix\n");
        printf("  Nr matrix elements:\n");
        printf("    tot         : %ld \n", A.NrNzElts);
        printf("    avg = tot/P : %g  \n", AvgNrNzElts);
        printf("    max         : %ld \n", MaxNrNzElts);
        printf("    min         : %ld \n", MinNrNzElts);

        printf("    imbalance   : %g %c\n", 
                    100 * (MaxNrNzElts / AvgNrNzElts - 1.0), '%');
        /* vertex weight information is not printed, since it is not
           meaningful here */ 
#endif   
    }
  
#ifdef INFO
    printf("\n****** Mondriaan vector distribution ******\n");
#endif
  
    /* Allocate memory for processor numbers for the vector components, 
       for u = A*v with A an m by n matrix */
    u_proc = (long int *) malloc(A.m * sizeof(long int));
    v_proc = (long int *) malloc(A.n * sizeof(long int));
    if (u_proc == NULL || v_proc == NULL) {
        fprintf(stderr, "main(): Not enough memory for vectors u_proc and v_proc!\n");
        exit(-1);
    }
      
    /**** Distribute the vector (and time it) ****/
    starttime = clock();
#ifdef UNIX
    gettimeofday(&starttime1, NULL);
#endif
   
    if (A.m == A.n && Options.SquareMatrix_DistributeVectorsEqual == EqVecYes) {
        /* Distribute the vectors equally */
        if (symmetric &&
             Options.SymmetricMatrix_UseSingleEntry == SingleEntYes) {
             
            /* Distribute v independently, and then use the same distribution for u */
            ComV = DistributeVec(&A, v_proc, ROW, &Options);
            
            if (ComV < 0) {
                fprintf(stderr, "main(): Unable to distribute vectors!\n");
                exit(-1);
            }
            
            for (i = 0; i < A.m; i++)
                u_proc[i] = v_proc[i];
            ComU = ComV;
#ifdef INFO
            printf("\nCommunication results for u are the same as for v\n");
            printf("but with sends and receives reversed\n\n");
#endif
        } else {
            /* Distribute u and v together */
            ComU = DistributeVecOrigEq(&A, u_proc, v_proc, &Options);
            
            if (ComU < 0) {
                fprintf(stderr, "main(): Unable to distribute vectors!\n");
                exit(-1);
            }
            
            ComV = 0;
        }
    } else {
        /* Distribute the vectors independently */
        ComV = DistributeVec(&A, v_proc, ROW, &Options);
        if (ComV < 0) {
            fprintf(stderr, "main(): Unable to distribute input vector!\n");
            exit(-1);
        }

        ComU = DistributeVec(&A, u_proc, COL, &Options);
        if (ComU < 0) {
            fprintf(stderr, "main(): Unable to distribute output vector!\n");
            exit(-1);
        }
    }
  
    endtime = clock();
#ifdef UNIX
    gettimeofday(&endtime1, NULL);
#endif
    cputime = ((double) (endtime - starttime)) / CLOCKS_PER_SEC;
  
#ifdef TIME
    printf("  vector distribution CPU-time      : %f seconds\n", cputime);
#ifdef UNIX  
    printf("  vector distribution elapsed time  : %f seconds\n",
             (endtime1.tv_sec - starttime1.tv_sec) +
             (endtime1.tv_usec - starttime1.tv_usec) / 1000000.0); 
#endif   
    fflush(stdout);
#endif  
      
#ifdef INFO
    printf("\nCommunication cost for u = A v: %ld (%g) \n", ComU + ComV,
           ((double) ComU / A.m  + (double) ComV / A.n) /2);  
    fflush(stdout);
#endif
 
    /* Write vector info */
    if (Options.OutputMode == MultipleFiles) {
        /* Write the vector distribution to file */
        sprintf(output, "%s-u%d", Options.matrix, A.NrProcs);
        File = fopen(output, "w");
    
        if (!File) fprintf(stderr, "main(): Unable to open '%s' for writing!\n", output);
        else {
            WriteVectorDistribution(u_proc, NULL, A.m, A.NrProcs, File, &Options);
            fclose(File);
        }
  
        sprintf(output, "%s-v%d", Options.matrix, A.NrProcs);    
        File = fopen(output, "w");
    
        if (!File) fprintf(stderr, "main(): Unable to open '%s' for writing!\n", output);
        else {
            WriteVectorDistribution(v_proc, NULL, A.n, A.NrProcs, File, &Options);
            fclose(File);
        }
    } else if (Options.OutputMode == OneFile) {
        /* Write the vector distribution to file */
        WriteVectorDistribution(v_proc, "Input-vector", A.n, A.NrProcs, File, &Options);
        WriteVectorDistribution(u_proc, "Output-vector", A.m, A.NrProcs, File, &Options);

        /* LocalRow2Processor */
        inv_count = NULL;
        row_local2proc = row_local2index = col_local2proc = col_local2index = NULL;
       if(!SparseMatrixLocal2Vector(&A, row_perms, u_proc, &inv_count, &row_local2proc, &row_local2index, 0)) {
            fprintf(stderr, "Error during derivation of local index arrays for SpMV multiplication (row direction)\n");
            exit(-1);
        }
        WriteVector(inv_count, 0, "OutputVectorLengths", A.NrProcs, 0, File, &Options);
        WriteVectorCollection(row_local2proc, "LocalRow2Processor", A.NrProcs, row_perms[A.NrProcs], File);
        WriteVectorCollection(row_local2index, "LocalRow2Index", A.NrProcs, row_perms[A.NrProcs], File);
        for( i=0; i<A.NrProcs; i++ ) {
            free( row_local2proc[i] );
            free( row_local2index[i] );
        }
        free( row_local2proc ); free( row_local2index );
        /* LocalCol2Processor */
        free(inv_count);
       if(!SparseMatrixLocal2Vector(&A, col_perms, v_proc, &inv_count, &col_local2proc, &col_local2index, 1)) {
            fprintf(stderr, "Error during derivation of local index arrays for SpMV multiplication (column direction)\n");
            exit(-1);
        }
        WriteVector(inv_count, 0, "InputVectorLengths", A.NrProcs, 0, File, &Options);
        WriteVectorCollection(col_local2proc, "LocalCol2Processor", A.NrProcs, col_perms[A.NrProcs], File);
        WriteVectorCollection(col_local2index, "LocalCol2Index", A.NrProcs, col_perms[A.NrProcs], File);
        for( i=0; i<A.NrProcs; i++ ) {
            free( col_local2proc[i] );
            free( col_local2index[i] );
        }
        free( col_local2proc ); free( col_local2index );
        free(inv_count);
        /* Also free local to global index */
        for( i=0; i<A.NrProcs+1; i++ ) {
            free( row_perms[i] );
            free( col_perms[i] );
        }
        free( row_perms ); free( col_perms );
        fclose(File);
    } else if (Options.OutputMode == DIMACS) {
        if (A.m != A.n || Options.SquareMatrix_DistributeVectorsEqual != EqVecYes) {
            fprintf(stderr, "main(): Unequal vector distributions in DIMACS mode!\n");
        }
        
        /* Only write the vector distribution to disk. */
        sprintf(output, "%s%d.part", Options.matrix, A.NrProcs);
        File = fopen(output, "w");
    
        if (!File) {
            fprintf(stderr, "main(): Unable to open '%s' for writing!\n", output);
        } else {
            for (i = 0; i < A.m; i++) fprintf(File, "%ld\n", u_proc[i]);
            fclose(File);
        }
    }
  
    if (Options.OutputMode == MultipleFiles) {
        /* Write the index sets of the Cartesian submatrices to file */
        sprintf(output, "%s-C%d", Options.matrix, A.NrProcs);
        File = fopen(output, "w");
        if (!File) fprintf(stderr, "main(): Unable to open '%s' for writing!\n", output);
        else {
            MMWriteCartesianSubmatrices(&A, File);
            fclose(File);
        }
    }

    /* Free memory */
    MMDeleteSparseMatrix(&A);  
    free(v_proc);
    free(u_proc);

    /* Exit :) */
    exit(0);
    
} /* end main */

