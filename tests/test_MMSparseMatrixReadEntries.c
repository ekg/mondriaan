#include "SparseMatrix.h"
#include "Options.h"

#define TOL 0.0001 

int main(int argc, char **argv) {

    FILE *fp;
    char filename[MAX_WORD_LENGTH];
    struct sparsematrix A;

    printf("Test MMSparseMatrixReadEntries: ");

    strcpy(filename,"test_MMSparseMatrixReadEntries.inp");
    fp = fopen(filename, "r");
    
    if (!fp) {
        printf("Error\n");
        exit(1);
    }

    A.MMTypeCode[2] = 'C'; /* complex matrix */
    A.MMTypeCode[3] = 'K'; /* skew-symmetric matrix */
    A.NrNzElts = 3;
    A.m = 3;
    A.n = 3;
    
    /* Allocate matrix arrays */
    A.i = (long *) malloc(A.NrNzElts* sizeof(long));
    A.j = (long *) malloc(A.NrNzElts* sizeof(long));
    A.ReValue = (double *) malloc(A.NrNzElts* sizeof(double));
    A.ImValue = (double *) malloc(A.NrNzElts* sizeof(double));
    
    if (A.i == NULL || A.j == NULL || A.ReValue  == NULL || A.ImValue == NULL) {
        printf("Error\n");
        exit(1);
    }
    
    MMSparseMatrixReadEntries(&A, fp);

    fclose(fp);

    /* Check that matrix dimensions and type are the same */
    if (A.NrNzElts != 3 ||
        A.m != 3 ||
        A.n != 3 ||
        A.MMTypeCode[2] != 'C' ||
        A.MMTypeCode[3] != 'K') {

        printf("Error\n");
        exit(1);
    }
    

    /* Check that matrix entries have been read correctly */
    if (A.i[0] != 1 || A.j[0] != 0 || A.ReValue[0] !=  A.ImValue[0] ||
        fabs(A.ReValue[0] - 1.0) > TOL ||
        A.i[1] != 2 || A.j[1] != 0 || A.ReValue[1] !=  A.ImValue[1] ||
        fabs(A.ReValue[1] - 0.5) > TOL ||
        A.i[2] != 2 || A.j[2] != 1 || A.ReValue[2] !=  A.ImValue[2]  ||
        fabs(A.ReValue[2] - 0.025) > TOL) {

        printf("Error\n");
        exit(1);
    }

    printf("OK\n");
    exit(0);
    
} /* end main */
