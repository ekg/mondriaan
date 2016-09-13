#include "SparseMatrix.h"

int main(int argc, char **argv) {

    struct sparsematrix A ;
    long n, i, j, t, *Val ;

    printf("Test AddDummiesToSparseMatrix: ");
    n = 180 ; /* n by n complex diagonal matrix A, n even.
                 The matrix has nonzeros in the even positions
                 on the main diagonal. */
    A.m = n;
    A.n = n;
    A.NrNzElts = n/2 ; 

    A.i = (long *) malloc(A.NrNzElts* sizeof(long)) ;
    A.j = (long *) malloc(A.NrNzElts* sizeof(long)) ;
    A.ReValue = (double *) malloc(A.NrNzElts* sizeof(double)) ;
    A.ImValue = (double *) malloc(A.NrNzElts* sizeof(double)) ;
    Val = (long *) malloc(n* sizeof(long)) ;

    if ( A.i == NULL || A.j == NULL || A.ReValue  == NULL ||
         A.ImValue == NULL || Val == NULL  ){
        printf("Error\n") ;
        exit(1);
    }

    /* Initialise matrix */
    A.MMTypeCode[2]='C'; /* complex matrix */
    A.MMTypeCode[3]='G'; /* general matrix */


    t= 0;
    for (i=0; i<n; i++){
        if (i%2 == 0 ) {
            A.i[t] = i ;
            A.j[t] = i ;
            A.ReValue[t] = i ;
            A.ImValue[t] = i ;
            t++;
        }
    }    

    AddDummiesToSparseMatrix(&A) ;

    /* Check result values  */
    if ( A.m != n || A.n != n || A.NrNzElts != n || A.NrDummies != n/2 ||
         A.MMTypeCode[2] != 'C' || A.MMTypeCode[3] != 'G') {
        printf("Error\n") ;
        exit(1);
    }

    /* Check that the dummies are in the odd positions */
    for (i=0; i<n; i++) {
        if ( (i%2 == 0 && A.dummy[i] ) ||
             (i%2 == 1 && A.dummy[i] == FALSE ) ) {
            printf("Error\n") ;
            exit(1);
        }
    }

    /* Check nonzeros */
    for (i=0; i<n; i++)
        Val[i] = 0;

    for (t = 0; t < A.NrNzElts ; t++){
        i=A.i[t] ;
        j=A.j[t] ;
        if ( i != j ||  A.ReValue[t] != A.ImValue[t]  || 
             (i%2 == 0 && A.ReValue[t] != (double)i ) ||
             (i%2 == 1 && A.ReValue[t] != 0.0       ) ) {
            printf("Error \n") ;
            exit(1);
        }
        Val[i] = i; /* mark the nonzero */
    }

    for (i=0; i<n; i++){
        if ( Val[i] != i ) {
            printf("Error \n") ;
            exit(1);
        }
    }

    printf("OK\n") ;
    exit(0);

} /* end main */
