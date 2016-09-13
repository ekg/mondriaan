#include "SparseMatrix.h"
#include "Options.h"

int main(int argc, char **argv) {

    long t ;
    struct sparsematrix A ;

    A.MMTypeCode[0] = 'W' ;
    A.NrRowWeights = 5 ;
    A.NrColWeights = 6 ;

    /* Allocate matrix arrays */
    A.RowWeights = (long *) malloc(A.NrRowWeights* sizeof(long)) ;
    A.ColWeights = (long *) malloc(A.NrColWeights* sizeof(long)) ;
    if ( A.RowWeights == NULL || A.ColWeights == NULL){
        printf("Error\n") ;
        exit(1);
    }

    for ( t = 0 ; t < A.NrRowWeights ; t++ ) 
        A.RowWeights[t] = t+1 ;
    for ( t = 0 ; t < A.NrColWeights ; t++ ) 
        A.ColWeights[t] = t+1 ;

    MMSparseMatrixPrintWeights(&A, stdout) ;
    exit(0) ;

} /* end main */
