#include "SparseMatrix.h"

int main(int argc, char **argv) {

    long m, n, nz, t, sumi, sumj ;
    double sumre, sumim ;
    int P ;

    struct sparsematrix A ;

    printf("Test MMSparseMatrixAllocateMemory: ");
    m = 24 ;
    n = 17 ;
    nz = 10 ;
    P = 4 ;
    A.m = m ;
    A.n = n ;
    A.NrNzElts = nz ;
    A.NrProcs = P ;
    A.MMTypeCode[2] = 'C' ; /* complex matrix */
    A.MMTypeCode[0] = 'D' ; /* distributed matrix, we do not test
                               a weighted matrix */

    MMSparseMatrixAllocateMemory(&A) ;

    if ( A.i == NULL ||  A.j == NULL || A.ReValue == NULL ||
         A.ImValue == NULL || A.Pstart == NULL ){
        printf("Error\n") ;
        exit(1);
    }

    for (t=0; t<nz; t++){
        A.i[t] = t ;
        A.j[t] = -t ;
        A.ReValue[t] = t ;
        A.ImValue[t] = -t ;
    }
    
    /* Compute sum */
    sumi = 0;
    sumj = 0;
    sumre = 0.0;
    sumim = 0.0;
    for (t=0; t<nz; t++){
        sumi += A.i[t] ;
        sumj += A.j[t] ;
        sumre += A.ReValue[t] ;
        sumim += A.ImValue[t] ;
    }
    
    /* Check sum of array elements */
    if (sumi != (nz-1)*nz/2 ||
        sumj != -sumi ||
        sumi != (long)(sumre+0.5) || 
        sumi != (long)(-sumim+0.5) ){
        printf("Error\n") ;
        exit(1);
    }

    printf("OK\n") ;
    exit(0);

} /* end main */
