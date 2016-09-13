#include "SparseMatrix.h"

int main(int argc, char **argv) {

    struct sparsematrix A ;

    printf("Test MMSparseMatrixInit: ");

    MMSparseMatrixInit(&A) ;

    /* Check result values  */
    if (A.MMTypeCode[0] != ' ' ||
        A.MMTypeCode[1] != ' ' ||
        A.MMTypeCode[2] != ' ' ||
        A.MMTypeCode[3] != ' ' ||

        strcmp(A.Banner,"") != 0 ||
        strcmp(A.Object,"") != 0 ||
        strcmp(A.Format,"") != 0 ||
        strcmp(A.Field,"") != 0 ||
        strcmp(A.Symmetry,"") != 0 ||
    
        A.NrNzElts != 0 ||
        A.m != 0 ||
        A.n != 0 ||
        A.NrDummies != 0 ||
        A.NrProcs != 0 ||
    
        A.i != NULL ||
        A.j != NULL ||
        A.ReValue != NULL ||
        A.ImValue != NULL ||
    
        A.dummy != NULL ||
        A.Pstart != NULL ||
        A.header != NULL ||
        A.tail != NULL ||
    
        A.NrRowWeights != 0 ||
        A.NrColWeights != 0 ||
        A.RowWeights != NULL ||
        A.ColWeights != NULL ||
        
        A.RowLambda != NULL ||
        A.ColLambda != NULL ||
    
        A.RowMark != NULL ||
        A.ColMark != NULL ) {

        printf("Error\n") ;
        exit(1) ;
    }
    
    printf("OK\n") ;
    exit(0) ;
    
} /* end main */
