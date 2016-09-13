#include "SparseMatrix.h"
#include "Options.h"

int main(int argc, char **argv) {

    FILE *fp;
    char filename[MAX_WORD_LENGTH];
    struct sparsematrix A;
    long m, n, i, j, weight;

    printf("Test MMSparseMatrixReadWeights: ");

    strcpy(filename,"test_MMSparseMatrixReadWeights.inp");
    fp = fopen(filename, "r");
    
    if (!fp) {
        printf("Error\n");
        exit(1);
    }

    m = 5; /* number of row weights */
    n = 4; /* number of column weights */
    A.MMTypeCode[0] = 'W'; /* weighted matrix */
    A.NrRowWeights = m;
    A.NrColWeights = n;
    
    /* Allocate matrix arrays */
    A.RowWeights = (long *) malloc(A.NrRowWeights* sizeof(long));
    A.ColWeights = (long *) malloc(A.NrColWeights* sizeof(long));
    
    if (A.RowWeights == NULL || A.ColWeights == NULL) {
        printf("Error\n");
        exit(1);
    }
    
    MMSparseMatrixReadWeights(&A, fp);

    fclose(fp);

    /* Check that matrix dimensions and type are the same */
    if (A.MMTypeCode[0] != 'W' ||
        A.NrRowWeights !=  m ||
        A.NrColWeights !=  n) {

        printf("Error\n");
        exit(1);
    }
    

    /* Check that row entries have been read correctly */
    weight =  1;
    for (i = 0; i < A.NrRowWeights; i++) {
        if (A.RowWeights[i] != weight) {
            printf("Error\n");
            exit(1);
        }
        weight *= 2;
    }

    /* Check that column entries have been read correctly */
    weight =  1;
    for (j = 0; j < A.NrColWeights; j++) {
        if (A.ColWeights[j] != weight) {
            printf("Error\n");
            exit(1);
        }
        weight *= 2;
    }

    printf("OK\n");
    exit(0);
    
} /* end main */
