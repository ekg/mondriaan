#include "SparseMatrix.h"

int main(int argc, char **argv) {

    struct sparsematrix A ;
    long n, P, i, q ;

    printf("Test RemoveDummiesFromSparseMatrix: ");
    P = 4 ;     /* number of processors */
    n = 73*P ;  /* n by n real diagonal matrix A, with dummies
                   on the main diagonal. */
    A.m = n;
    A.n = n;
    A.NrNzElts = n ; 
    A.NrDummies = n ; 
    A.NrProcs = P ; 

    A.i = (long *) malloc(A.NrNzElts* sizeof(long)) ;
    A.j = (long *) malloc(A.NrNzElts* sizeof(long)) ;
    A.ReValue = (double *) malloc(A.NrNzElts* sizeof(double)) ;
    A.dummy = (int *) malloc(n* sizeof(int)) ;
    A.Pstart = (long *) malloc((P+1)* sizeof(long)) ;

    if ( A.i == NULL || A.j == NULL || A.ReValue  == NULL ||
         A.dummy == NULL || A.Pstart == NULL ){
        printf("Error\n") ;
        exit(1);
    }

    /* Initialise matrix */
    A.MMTypeCode[2]='R'; /* real matrix */
    A.MMTypeCode[3]='G'; /* general matrix */

    for (q=0; q<=P; q++)
        A.Pstart[q] = q*(n/P) ; 

    for (i=0; i<n; i++){
        A.i[i] = i ;
        A.j[i] = i ;
        A.ReValue[i] = i ;
        A.dummy[i] = TRUE ;
    }    

    RemoveDummiesFromSparseMatrix(&A) ;

    /* Check result values  */
    if ( A.m != n || A.n != n || A.NrNzElts != 0 || A.NrDummies != 0 ||
         A.NrProcs != P || A.MMTypeCode[2] != 'R' || A.MMTypeCode[3] != 'G') {
        printf("Error \n") ;
        exit(1);
    }

    /* Check that the partitioning information is adjusted */
    for (q=0; q<=P ; q++){
        if ( A.Pstart[q] != 0) { 
            printf("Error\n") ;
            exit(1);
        }
    }

    printf("OK\n") ;
    exit(0);

} /* end main */
