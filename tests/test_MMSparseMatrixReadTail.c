#include "SparseMatrix.h"
#include "Options.h"

int main(int argc, char **argv) {

    FILE *fp;
    char filename[MAX_WORD_LENGTH];
    struct sparsematrix A;

    printf("Test MMSparseMatrixReadTail: ");

    strcpy(filename,"test_MMSparseMatrixReadTail.inp");
    fp = fopen(filename, "r");
    
    if (!fp) {
        printf("Error\n");
        exit(1);
    }

    MMSparseMatrixReadTail(&A, fp);

    fclose(fp);

    /* Check that tail has been read correctly */
    if (strcmp(A.tail,"Tail 1\nTail 2\n") != 0) {
        printf("Error 1\n");
        exit(1);
    }

    printf("OK\n");
    exit(0);
    
} /* end main */
