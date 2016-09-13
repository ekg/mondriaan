#include "SparseMatrix.h"
#include "Options.h"

int main(int argc, char **argv) {

    FILE *fp;
    char filename[MAX_WORD_LENGTH];
    struct sparsematrix A;
    long nzp;
    int q, P;

    printf("Test MMSparseMatrixReadPstart: ");

    strcpy(filename,"test_MMSparseMatrixReadPstart.inp");
    fp = fopen(filename, "r");
    
    if (!fp) {
        printf("Error\n");
        exit(1);
    }

    nzp  = 10; /* number of nonzeros per processor */
    P = 10;
    A.MMTypeCode[0] = 'D'; /* distributed matrix */
    A.NrNzElts = P * nzp;
    A.NrProcs = P;
    
    /* Allocate matrix arrays */
    A.Pstart = (long *) malloc((A.NrProcs+1)* sizeof(long));
    
    if (A.Pstart == NULL){
        printf("Error\n");
        exit(1);
    }
    
    MMSparseMatrixReadPstart(&A, fp);

    fclose(fp);

    /* Check that matrix dimensions and type are the same */
    if (A.MMTypeCode[0] != 'D' ||
        A.NrNzElts !=  P * nzp ||
        A.NrProcs != P) {

        printf("Error\n");
        exit(1);
    }
    

    /* Check that matrix entries have been read correctly */
    for (q = 0; q <= A.NrProcs; q++) {
        if (A.Pstart[q] != q * nzp) {
            printf("Error\n");
            exit(1);
        }
    }

    printf("OK\n");
    exit(0);
    
} /* end main */
