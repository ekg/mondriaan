#include "SparseMatrix.h"
#include "Matalloc.h"

int main(int argc, char **argv) {

    struct sparsematrix A ;
    long n, nz, i, j, t, val, **a ;

    printf("Test SparseMatrixRemoveDuplicates: ");
    n = 19 ; /* n by n real dense matrix A */
    A.m = n;
    A.n = n;
    A.NrNzElts = 2*n*n ; 
    nz = A.NrNzElts ;

    A.i = (long *) malloc(nz* sizeof(long)) ;
    A.j = (long *) malloc(nz* sizeof(long)) ;
    A.ReValue = (double *) malloc(nz* sizeof(double)) ;
    a = (long **) matallocl(n,n);
 
    if ( A.i == NULL || A.j == NULL || A.ReValue  == NULL || a == NULL ){
        printf("Error\n") ;
        exit(1);
    }

    /* Initialise matrix */
    A.MMTypeCode[2]='R'; /* real matrix */
 
    /* Initialise dense matrix */
    for (i=0; i<n; i++)
        for (j=0; j<n; j++)
            a[i][j] = 0 ;
            
    /* Fill dense and sparse matrices with nonzeros */
    t = 0;
    for (i=0; i<n; i++){
        for (j=0; j<n; j++){
            A.i[t+n*n] = A.i[t] = i ;
            A.j[t+n*n] = A.j[t] = j ;
            A.ReValue[t+n*n] = A.ReValue[t] = i+j*n ;
            a[i][j] = 2*(i+j*n) ;
            t++;
        }
    }    

    SparseMatrixRemoveDuplicates(&A) ;

    if ( A.m != n || A.n != n  || 2*A.NrNzElts != nz || 
         A.MMTypeCode[2] != 'R') {
        printf("Error\n") ;
        exit(1);
    }

    /* Check nonzero values */
    for (t=0; t<A.NrNzElts; t++) {
        i = A.i[t] ;
        j = A.j[t] ;
        val = A.ReValue[t] + 0.5  ; /* round to nearest integer value */
        if (val != 2*(i+j*n) ){
            printf("Error\n") ;
            exit(1);
        }  
        a[i][j] -= val ; /* reset dense matrix element */
    }

    /* Check whether whole dense matrix is zero again */
    for (i=0; i<n; i++){
        for (j=0; j<n; j++){
            if (a[i][j] != 0) {
                printf("Error\n") ;
                exit(1);
            }
        }    
    }

    printf("OK\n") ;
    exit(0);

} /* end main */
