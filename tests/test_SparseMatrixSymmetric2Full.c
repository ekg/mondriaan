#include "SparseMatrix.h"

int main(int argc, char **argv) {

    struct sparsematrix A ;
    long n, i, j, t, q ;

    printf("Test SparseMatrixSymmetric2Full: ");
    n = 16 ; /* n by n hermitian tridiagonal matrix A, n>=3.
                 The matrix is represented by the nonzeros on or below the diagonal. */
    A.m = n;
    A.n = n;
    A.NrNzElts = 2*n-1 ; 
    A.NrProcs = n ; /* Each matrix row in a different processor */

    A.i = (long *) malloc(A.NrNzElts* sizeof(long)) ;
    A.j = (long *) malloc(A.NrNzElts* sizeof(long)) ;
    A.ReValue = (double *) malloc(A.NrNzElts* sizeof(double)) ;
    A.ImValue = (double *) malloc(A.NrNzElts* sizeof(double)) ;
    A.Pstart = (long *) malloc((A.NrProcs+1) * sizeof(long)) ;

    if ( A.i == NULL || A.j == NULL || A.ReValue  == NULL || A.ImValue  == NULL ||
         A.Pstart== NULL ){
        printf("Error\n") ;
        exit(1);
    }

    /* Initialise matrix */
    A.MMTypeCode[2]='C'; /* complex matrix */
    A.MMTypeCode[3]='H'; /* hermitian matrix */

    t= 0;
    A.Pstart[0] = 0 ;
    for (i=0; i<n; i++){
        for (j=MAX(0,i-1); j<=i; j++){
            A.i[t] = i ;
            A.j[t] = j ;
            A.ReValue[t] = i ;
            A.ImValue[t] = j ;
            t++;
        }
        A.Pstart[i+1] = t ;
    }    

    SparseMatrixSymmetric2Full(&A) ;

    /* Check result values  */
    if ( A.m != n || A.n != n || A.NrNzElts != 3*n-2 || A.NrProcs != n ||
         A.MMTypeCode[2] != 'C' || A.MMTypeCode[3] != 'G' ||
         strcmp(A.Symmetry,"general") != 0 || A.Pstart[0] != 0 ) {
        printf("Error 2\n") ;
        exit(1);
    }

    /* Check Pstart */
    for (q = 1; q <= n; q++){
        if ( A.Pstart[q] != 3*q-2) { 
            printf("Error 3\n") ;
            exit(1);
        }
    }

    /* Check nonzeros */
    for (q = 0; q < n; q++){
        /* Check nonzeros of processor q */
        for (t=A.Pstart[q]; t<A.Pstart[q+1]; t++){
            i=A.i[t] ;
            j=A.j[t] ;
            if ( ((i < j ) && (i != q-1 || j != q || A.ReValue[t] != q 
                                                  || A.ImValue[t] != -(q-1))) ||
                 ((i > j ) && (i != q || j != q-1 || A.ReValue[t] != q 
                                                  || A.ImValue[t] != q-1   )) ||
                 ((i == j) && (i != q || j != q   || A.ReValue[t] != q 
                                                  || A.ImValue[t] != q     )) ) {
                printf("Error 4\n") ;
                exit(1);
            }
        }
    }

    printf("OK\n") ;
    exit(0);

} /* end main */
