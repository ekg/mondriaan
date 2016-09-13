#include "SparseMatrix.h"
#include "Options.h"

int main(int argc, char **argv) {

    long t;
    int sign;
    struct sparsematrix A;

    A.MMTypeCode[2] = 'I';
    A.NrNzElts = 11;
    
    /* Allocate matrix arrays */
    A.i = (long *) malloc(A.NrNzElts* sizeof(long));
    A.j = (long *) malloc(A.NrNzElts* sizeof(long));
    A.ReValue = (double *) malloc(A.NrNzElts* sizeof(double));

    if ( A.i == NULL || A.j == NULL || A.ReValue == NULL ){
        printf("Error\n");
        exit(1);
    }

    sign= 1; /* alternating sign */
    
    /* Diagonal matrix (1,2,...,nz) */
    for ( t = 0; t < A.NrNzElts; t++ ) {
        A.i[t] = t;
        A.j[t] = t;
        A.ReValue[t] = t+1+ sign*0.001; /* small disturbance */
        sign = -sign;
    }

    MMSparseMatrixPrintEntries(&A, stdout);
    exit(0);

} /* end main */
