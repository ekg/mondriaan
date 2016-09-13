#include "SparseMatrix.h"
#include "Options.h"

int main(int argc, char **argv) {

    char line[MAX_LINE_LENGTH] ;
    struct sparsematrix A ;
    struct opts Options;
    SetDefaultOptions( &Options );

    A.MMTypeCode[0] = 'W' ;
    A.MMTypeCode[1] = 'C' ;
    A.MMTypeCode[2] = 'R' ;
    A.MMTypeCode[3] = 'G' ;
    A.ViewType      = ViewTypeOriginal;
    A.m = 3 ;
    A.n = 3 ;
    A.NrNzElts = 9 ;
    A.NrRowWeights = 3 ;
    A.NrColWeights = 3 ;
    strcpy(line, "% This is a 3 x 3 matrix\n% with 9 nonzeros\n") ;
    A.header  = line ;
    A.comment = NULL;

    MMSparseMatrixPrintHeader(&A, stdout, NULL, &Options);
    exit(0) ;

} /* end main */
