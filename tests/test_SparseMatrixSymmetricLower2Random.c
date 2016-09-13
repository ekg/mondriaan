#include "SparseMatrix.h"
#include "Matalloc.h"

int main(int argc, char **argv) {

    struct sparsematrix A ;
    long n, nz, i, j, k, t, tmp, val, **a ;

    printf("Test SparseMatrixSymmetricLower2Random: ");
    n = 40 ; /* n by n real matrix A with stored nonzero diagonal
                i = j + k, where k >= 1. 
                The diagonal i = j - k is not stored. */
    k = 7 ; 
    A.m = n;
    A.n = n;
    A.NrNzElts = n-k ; 
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
    A.MMTypeCode[3]='S'; /* symmetric matrix */
 
    /* Initialise dense matrix */
    for (i=0; i<n; i++)
        for (j=0; j<n; j++)
            a[i][j] = 0 ;
            
    /* Fill dense and sparse matrices with nonzeros */
    t = 0;
    for (i=k; i<n; i++){
        j = i - k ;
        A.i[t] = i ;
        A.j[t] = j ;
        a[i][j] = i ;
        A.ReValue[t] = a[i][j] ;
        t++;
    }    

    SparseMatrixSymmetricLower2Random(&A) ;

    if ( A.m != n || A.n != n  || A.NrNzElts != nz || 
         A.MMTypeCode[2] != 'R' || A.MMTypeCode[3] != 'S' ) {
        printf("Error\n") ;
        exit(1);
    }

    /* Check nonzero values */
    for (t=0; t<nz; t++) {
        i = A.i[t] ;
        j = A.j[t] ;
        if (i < j) { 
            tmp = i ;
            i = j ;
            j = tmp ;
        }
        val = A.ReValue[t] + 0.5  ; /* round to nearest integer value */
        if (val != i){
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
