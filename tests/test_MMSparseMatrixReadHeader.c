#include "SparseMatrix.h"
#include "Options.h"

int main(int argc, char **argv) {

    FILE *fp;
    char filename[MAX_WORD_LENGTH];
    struct sparsematrix A;

    printf("Test MMSparseMatrixReadHeader: ");

    strcpy(filename,"test_MMSparseMatrixReadHeader.inp");
    fp = fopen(filename, "r");
    
    if (!fp) {
        printf("Error\n");
        exit(1);
    }

    A.m = 0;
    A.n = 0;
    A.NrNzElts = 0;
    A.NrRowWeights = 0;
    A.NrColWeights = 0;
    A.header = NULL;
    A.MMTypeCode[0] = ' ';

    MMSparseMatrixReadHeader(&A, fp);

    fclose(fp);

    /* Check that header has been read correctly */
    if (strcmp(A.header,"% This is the header comment\n") != 0 ||
         A.m != A.NrRowWeights || A.m != 5 ||
         A.n != A.NrColWeights || A.n != 7 ||
         A.NrNzElts != 13 || 
         A.MMTypeCode[0] != 'W' ||
         A.MMTypeCode[1] != 'C' ||
         A.MMTypeCode[2] != 'R' ||
         A.MMTypeCode[3] != 'G') {
        printf("Error\n");
        exit(1);
    }

    printf("OK\n");
    exit(0);
    
} /* end main */
