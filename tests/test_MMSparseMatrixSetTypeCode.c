#include "SparseMatrix.h"

int main(int argc, char **argv) {

    char line[MAX_LINE_LENGTH];
    struct sparsematrix A;

    printf("Test MMSparseMatrixSetTypeCode: ");
    strcpy(line, "%%MatrixMarket     WeightedMatrix   Sparse     Pattern Skew-Symmetric");

    MMSparseMatrixSetTypeCode(line, &A);
    
    /* Check 5 strings and 4 type codes*/
    if (strcmp(A.Banner, MM_Banner) != 0 ||
        strcmp(A.Object, "weightedmatrix") != 0 ||
        strcmp(A.Format, "sparse") != 0 ||
        strcmp(A.Field, "pattern") != 0 ||
        strcmp(A.Symmetry, "skew-symmetric") != 0 ||
        A.MMTypeCode[0] != 'W' ||
        A.MMTypeCode[1] != 'C' ||
        A.MMTypeCode[2] != 'P' ||
        A.MMTypeCode[3] != 'K' ){
        printf("Error\n");
        exit(1);
    }

    printf("OK\n");
    exit(0);

} /* end main */
