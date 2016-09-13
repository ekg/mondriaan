#include "SparseMatrix.h"
#include "Options.h"

int main(int argc, char **argv) {

    int t ;
    struct sparsematrix A ;
    struct opts Options;

    SetDefaultOptions(&Options);

    A.MMTypeCode[0] = 'D' ;
    A.NrProcs = 7 ;

    /* Allocate matrix arrays */
    A.Pstart = (long *) malloc((A.NrProcs+1)* sizeof(long)) ;
    if ( A.Pstart == NULL ){
        printf("Error\n") ;
        exit(1);
    }

    for ( t = 0 ; t <= A.NrProcs ; t++ ) 
        A.Pstart[t] = 10*t ;

    MMSparseMatrixPrintPstart(&A, stdout, &Options) ;
    exit(0) ;

} /* end main */
