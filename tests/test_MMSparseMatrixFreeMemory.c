#include "SparseMatrix.h"

int main(int argc, char **argv) {

    long m, n, nz, t, i, j ;

    struct sparsematrix A ;

    printf("Test MMSparseMatrixFreeMemory: ");
    
    MMSparseMatrixInit(&A) ;
    m = 2400 ;
    n = 1700 ;
    nz = 1000 ;
    A.m = m ;
    A.n = n ;
    A.NrNzElts = nz ;
    A.NrRowWeights = m ;
    A.NrColWeights = n ;
    A.MMTypeCode[2] = 'C' ; /* complex matrix */
    A.MMTypeCode[0] = 'W' ; /* weighted matrix, we do not test
                               a distributed matrix */

    MMSparseMatrixAllocateMemory(&A) ;

    if ( A.i == NULL ||  A.j == NULL ||
         A.ReValue == NULL || A.ImValue == NULL || 
         A.RowWeights == NULL || A.ColWeights == NULL){
        printf("Error\n") ;
        exit(1);
    }

    for (t=0; t<nz; t++){
        A.i[t] = t ;
        A.j[t] = t ;
        A.ReValue[t] = t ;
        A.ImValue[t] = t ;
    }
    for (i=0 ; i<m ; i++)
        A.RowWeights[i] = 1;
    for (j=0 ; j<n ; j++)
        A.ColWeights[j] = 1;
    
    MMSparseMatrixFreeMemory(&A) ;

    printf("OK\n") ;
    exit(0);

} /* end main */
