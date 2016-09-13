#include "SparseMatrix.h"
#include "Options.h"

int main(int argc, char **argv) {

    char line[MAX_LINE_LENGTH] ;
    struct sparsematrix A ;

    strcpy(line, "The tail wags the dog\n") ;
    A.tail = line ;

    MMSparseMatrixPrintTail(&A, stdout) ;
    exit(0); 

} /* end main */
