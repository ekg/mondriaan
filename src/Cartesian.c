#include "Cartesian.h"

int MMWriteCartesianSubmatrices(const struct sparsematrix *pM, FILE *fp){

    /* This function computes and writes the Cartesian submatrices
       I(q) x J(q) for processors q, 0 <= q < P,
       where P = M.NrProcs, for a distributed sparse matrix M. 
       The submatrix I(q) x J(q) is defined as the smallest submatrix
       which contains all nonzeros of processor q.
       Thus, the rows and columns of the submatrix are nonempty. 
       They need not be contiguous.

       The format of the output file is:
         line 1: m n P    (nr of rows, columns, processors)
         for each q :
             1 line: q m(q) n(q)
             m(q) lines with one row index i
             n(q) lines with one column index j
       All output indexing is 1-based.
       
       This function returns TRUE on success and FALSE on failure.
    */

    long P, q, i, j, t, mq, nq, totalC, totalR;

    /* Column arrays */
    int *NprocsC;     /* NprocsC[j] is the number of processors
                         that occur in column j */
    int *NprocsR;     /* same for rows i */

    /* Communication arrays. These are used to store a matrix in 
       compressed row storage (CRS) or compressed column storage (CCS).  */
       
    int *indexC;  /* indexC[startC[j]..startC[j+1]-1] contains the processor numbers
                     of the processors that occur in column j */
    int *indexR;  /* same for rows i */
    long *indexC2, *indexR2; /* copies in long */
    
    long *startC, *startR;
    
    if (!pM || !fp) {
        fprintf(stderr, "MMWriteCartesianSubmatrices(): Null arguments!\n");
        return FALSE;
    }

    P = pM->NrProcs; 

    NprocsC = (int *)malloc(pM->n*sizeof(int));
    startC = (long *)malloc( (MAX(pM->n,P)+1)*sizeof(long)); 
           /* MAX(n,P) since it will be used both for CRS and CCS storage 
              of a P x n matrix */
    NprocsR = (int *)malloc(pM->m*sizeof(int));
    startR = (long *)malloc((MAX(pM->m,P)+1)*sizeof(long));
    
    if (NprocsC == NULL || startC == NULL || NprocsR == NULL || startR == NULL || !pM) {
        fprintf(stderr, "MMWriteCartesianSubmatrices(): Not enough memory!\n");
        return FALSE;
    }
 
    /* Compute the number of processors in each row and column */
    if (!InitNprocs(pM, ROW, NprocsC)) {
                                 /* indeed ROW; confusing, but this comes from the
                                    vector distribution, where ROW corresponds 
                                    to an input vector in the row direction. 
                                    Thus, its length is the number of columns. */
        fprintf(stderr, "MMWriteCartesianSubmatrices(): Unable to initialise processors!\n");
        return FALSE;
    }
    
    if (!InitNprocs(pM, COL, NprocsR)) {
        fprintf(stderr, "MMWriteCartesianSubmatrices(): Unable to initialise processors!\n");
        return FALSE;
    }

    /* Determine the processors in each row and column. 
       This creates: P x n matrix in CCS format (for C)
                     m x P matrix in CRS format (for R) */
    totalC = 0;
    for (j=0; j<pM->n; j++)
        totalC += NprocsC[j];
    
    totalR = 0;
    
    for (i=0; i<pM->m; i++)
        totalR += NprocsR[i];

    indexC = (int *)malloc(totalC*sizeof(int));
    indexR = (int *)malloc(totalR*sizeof(int));
    
    if (indexC == NULL || indexR == NULL) {
        fprintf(stderr, "MMWriteCartesianSubmatrices(): Not enough memory!\n");
        return FALSE;
    }

    if (!InitProcindex(pM, ROW, NprocsC, startC, indexC)) {
        fprintf(stderr, "MMWriteCartesianSubmatrices(): Unable to initialise processor index!\n");
        return FALSE;
    }
    if (!InitProcindex(pM, COL, NprocsR, startR, indexR)) {
        fprintf(stderr, "MMWriteCartesianSubmatrices(): Unable to initialise processor index!\n");
        return FALSE;
    }
    
    /* Copy index arrays of ints into arrays of longs */
    indexC2 = (long *)malloc(totalC*sizeof(long));
    indexR2 = (long *)malloc(totalR*sizeof(long));
    if (indexC2 == NULL || indexR2 == NULL) {
        fprintf(stderr, "MMWriteCartesianSubmatrices(): Not enough memory!\n");
        return FALSE;
    }
        
    for (t = 0; t < totalR; t++)
        indexR2[t] = indexR[t];
    
    for (t = 0; t < totalC; t++)
        indexC2[t] = indexC[t];
        
    /* Convert m x P matrix from CRS format to CCS */
    if (!CRS2CCS(pM->m, P, totalR, startR, indexR2)) {
        fprintf(stderr, "MMWriteCartesianSubmatrices(): Unable to convert from CRS to CCS!\n");
        return FALSE;
    }

    /* Convert P x n matrix from CCS format to CRS.
       Inverse direction, hence with dimensions reversed */
    if (!CRS2CCS(pM->n, P, totalC, startC, indexC2)) {
        fprintf(stderr, "MMWriteCartesianSubmatrices(): Unable to convert from CCS to CRS!\n");
        return FALSE;
    }
    
    /*** Write to file ***/
    if (ferror(fp)) {
        fprintf(stderr, "MMWriteCartesianSubmatrices(): Could not write to stream!\n");
        return FALSE;
    }
    
    fprintf(fp, "%ld %ld %ld\n", pM->m, pM->n, P);

    for (q = 0; q < P; q++){
        /* Write dimension of submatrix */
        mq = startR[q+1]-startR[q];
        nq = startC[q+1]-startC[q];
        fprintf(fp, "%ld %ld %ld\n", q+1, mq, nq);
        
        /* Write row indices */
        for (t = startR[q]; t < startR[q+1]; t++)
            fprintf(fp, "%ld\n", indexR2[t] + 1);
            
        /* Write column indices */
        for (t = startC[q]; t < startC[q+1]; t++)
            fprintf(fp, "%ld\n", indexC2[t] + 1);
    }
   
    free(indexR2);
    free(indexC2);
    free(indexR);
    free(indexC);
    free(startR);
    free(NprocsR);
    free(startC);
    free(NprocsC);
    
    return TRUE;
} /* end MMWriteCartesianSubmatrices */


int CRS2CCS(long m, long n, long nz, long *start, long *Index) {

    /* This function transposes a sparse m x n pattern matrix
       with nz >= 0 nonzeros from compressed row storage (CRS) to
       compressed column storage (CCS) by using linked lists.
       Additionally, the nonzeros within a column in CCS are guaranteed
       to be sorted by increasing row index.
       
       CRS is defined by an array start of length m+1, 
       where start[i] gives the position of the first nonzero of row i,
       and an array index of length nz, where index[k] is the column index
       of the k-th nonzero. furthermore, start[m] = nz.

       The function can also be used to convert an m x n matrix in CCS format
       to CRS format. It is then called with the m and n parameters reversed:
           CRS2CCS(n, m, nz, start, index);       
       
       The new format (start, index) overwrites the previous one.
       
       start array has length at least  MAX(m,n)+1
       index array has length nz  */ 

    long i, j, t, tt, *First, *Last, *Next, *I;

    if (nz < 0) {
        fprintf(stderr, "CRS2CCS(): Number of nonzeros < 0!\n");
        return FALSE;
    }
    else if (nz == 0){
        for (j = 0; j <= n; j++)
            start[j] = 0;
        
        return TRUE;
    } /* else nz > 0 */

    First = (long *) malloc(n*sizeof(long));
    Last =  (long *) malloc(n*sizeof(long));
    Next =  (long *) malloc(nz*sizeof(long));
    I =     (long *) malloc(nz*sizeof(long));
    
    if (First == NULL || Last == NULL || Next == NULL || I == NULL) {
        fprintf(stderr, "CRS2CCS(): Not enough memory!\n");
        return FALSE;
    }

    /* Initialise linked lists, one for each matrix column */
    for (j = 0; j < n; j++){
        First[j] =  DUMMY;    /* First nonzero in column j */
        Last[j] =  DUMMY;     /* Last nonzero in column j */
    }
    
    for (t = 0; t < nz; t++)
        Next[t] = DUMMY;      /* Next nonzero in column */

    /* Construct linked lists. The nonzeros are read out from the rows 
       in increasing row order. As a result, the nonzeros within a column
       list are sorted by increasing row index. */

    for (i = 0; i < m; i++) {
        for (t = start[i]; t < start[i+1]; t++) {
            I[t] = i;
            j = Index[t];
            
            if (First[j] == DUMMY)
                First[j] = t;
            else 
                Next[Last[j]] = t;
            
            Last[j] = t;
        }
    }
 
    /* Use linked lists */
    tt = 0;
    for (j = 0; j < n; j++){
        /* construct column j in CCS data structure */
        start[j] = tt;
        t = First[j];    
        
        while (t != DUMMY) {
            Index[tt] = I[t];
            tt++;
            t = Next[t];
        }
    }
    
    start[n] = nz;

    free(I);
    free(Next);
    free(Last); 
    free(First);
    
    return TRUE;
} /*  end CRS2CCS */
