#include "SparseMatrix.h"
#include "Options.h"
#include "Matalloc.h"

int main(int argc, char **argv) {

    char filename[MAX_WORD_LENGTH];
    struct sparsematrix A;
    long nz, n, i, j, k, **a;
    FILE *fp;

    printf("Test MMReadSparseMatrix: ");

    strcpy(filename,"test_MMReadSparseMatrix.inp");

    nz = 13;
    n = 5;

    /* Allocate matrix arrays */
    A.i = (long *) malloc(nz* sizeof(long));
    A.j = (long *) malloc(nz* sizeof(long));
    A.RowWeights = (long *) malloc(n* sizeof(long));
    A.ColWeights = (long *) malloc(n* sizeof(long));
    a = (long **) matallocl(n,n);

    
    if (A.i == NULL || A.j == NULL || A.RowWeights == NULL ||
         A.ColWeights == NULL || a == NULL) {
        printf("Error\n");
        exit(1);
    }
    
    fp = fopen(filename, "r");
    MMReadSparseMatrix(fp, &A);
    fclose(fp);

    /* Check matrix dimensions and types */
    if (A.NrNzElts != nz ||
        A.m != n ||
        A.n != n ||
        A.MMTypeCode[0] != 'W' ||
        A.MMTypeCode[1] != 'C' ||
        A.MMTypeCode[2] != 'P' ||
        A.MMTypeCode[3] != 'G' || 
        A.NrRowWeights != n ||
        A.NrColWeights != n ||
        strcmp(A.header,"% This is the header comment\n") != 0 /* ||
        strcmp(A.tail,"This is the tail comment\n") != 0 */) {

        printf("Error\n");
        exit(1);
    }
    

    /* Check that matrix entries have been read correctly */

    /* Initialise dense matrix */
    for (i=0; i<n; i++)
        for (j=0; j<n; j++)
            a[i][j] = 0;

    for (k=0; k<nz; k++) {
        i = A.i[k];
        j = A.j[k];
        if (i<0 || i>=n || j<0 || j>=n) {
            printf("Error\n");
            exit(1);
        }
        a[i][j]++;
    }

    for (i=0; i<n; i++)
        for (j=0; j<n; j++)
            if ((a[i][j] != 1 && (i+j)%2 == 0) ||
                 (a[i][j] != 0 && (i+j)%2 == 1)) {
                printf("Error\n");
                exit(1);
            }

    for (i=0; i<n; i++)
        if (A.RowWeights[i] != i+1 || A.ColWeights[i] != i+1) {
            printf("Error\n");
            exit(1);
        }

    printf("OK\n");
    exit(0);
    
} /* end main */
