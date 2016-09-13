#include "SparseMatrix.h"

int main(int argc, char **argv) {

    struct sparsematrix A ;
    long n, nz, i, j, t ;
    int symmetric ;

    printf("Test SparseMatrixStructurallySymmetric: ");
    n = 4 ; /* n by n pattern matrix A, checkerboard, n even */
    A.m = n;
    A.n = n;
    A.NrNzElts = n*n/2 ; 
    nz = A.NrNzElts ;

    A.i = (long *) malloc(nz* sizeof(long)) ;
    A.j = (long *) malloc(nz* sizeof(long)) ;

    if ( A.i == NULL || A.j == NULL ){
        printf("Error\n") ;
        exit(1);
    }

    /* Initialise matrix */
    A.MMTypeCode[2]='P'; /* pattern matrix */
    A.MMTypeCode[3]='G'; /* general matrix */
 
    /* Fill sparse matrix with nonzeros */
    t = 0;
    for (i=0; i<n; i++){
        for (j=0; j<n; j++){
            if ((i+j)%2 == 0) {
                A.i[t] = i ;
                A.j[t] = j ;
                t++;
            }
        }
    }    

    symmetric = SparseMatrixStructurallySymmetric(&A) ;

    if ( A.m != n || A.n != n   || A.NrNzElts != nz || 
         A.MMTypeCode[2] != 'P' || A.MMTypeCode[3] != 'G' || symmetric == FALSE) {
        printf("Error\n") ;
        exit(1);
    }

    /* Check nonzero values */
    for (t=0; t<A.NrNzElts; t++) {
        i = A.i[t] ;
        j = A.j[t] ;
        if ((i+j)%2 != 0) {
            printf("Error\n") ;
            exit(1);
        }  
    }
    
    /* Remove last two entries, to make the matrix unsymmetric */
    A.NrNzElts -= 2 ;
    
    symmetric = SparseMatrixStructurallySymmetric(&A) ;

    if ( A.m != n || A.n != n   || A.NrNzElts != nz-2 || 
         A.MMTypeCode[2] != 'P' || A.MMTypeCode[3] != 'G' || symmetric == TRUE ) {
        printf("Error\n") ;
        exit(1);
    }

    printf("OK\n") ;
    exit(0);

} /* end main */
