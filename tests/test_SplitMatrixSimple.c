#include <stdlib.h>
#include <stdio.h>

#include "DistributeMat.h"

struct opts Options;

extern int SplitMatrixSimple(struct sparsematrix *pT, int k, int i,
                       long weightlo, long weighthi, const struct opts *pOptions);

int main(int argc, char **argv) {

    struct sparsematrix A;
    long n, i, j, t, P, weightlo, weighthi, weight0, weight1;

    printf("Test SplitMatrixSimple: ");
    n = 33; /* n by n dense matrix A, n odd */
    A.m = n;
    A.n = n;
    A.NrNzElts = n*n;
    P = 2; /* maximum number of parts */

    A.i = (long *) malloc(A.NrNzElts* sizeof(long));
    A.j = (long *) malloc(A.NrNzElts* sizeof(long));
    A.Pstart = (long *) malloc((P+1)* sizeof(long));
    A.RowLambda = (int *)malloc(A.m*sizeof(int));
    A.ColLambda = (int *)malloc(A.n*sizeof(int));
    A.RowMark = (int *)malloc(A.m*sizeof(int));
    A.ColMark = (int *)malloc(A.n*sizeof(int));
    
    if (!SetDefaultOptions(&Options)) {
        printf("Error\n");
        exit(1);
    }

    if (A.i == NULL || A.j  == NULL || A.Pstart  == NULL) {
        printf("Error\n");
        exit(1);
    }

    /* Fill matrix with nonzeros */
    t= 0;
    for (i=0; i<n; i++) {
        for (j=0; j<n; j++) {
            A.i[t] = i;
            A.j[t] = j;
            t++;
        }
    }

    A.MMTypeCode[3]='G';
    A.NrDummies = 0;
    A.Pstart[0] = 0;
    A.Pstart[1] = A.NrNzElts;

    weightlo = (n*n)/2;
    weighthi = (n*n)/2 + 1;
    
    if (!SplitMatrixSimple(&A, 1, 0, weightlo, weighthi, &Options)) {
        printf("Error\n");
        exit(1);
    }

    /* Check result */
    if (A.Pstart[0] != 0 || A.Pstart[1] != (n*n)/2 ||
         A.Pstart[2] != A.NrNzElts) {
        printf("Error\n");
        exit(1);
    }
    weight0 = ComputeWeight(&A, A.Pstart[0], A.Pstart[1]-1, NULL, &Options);
    weight1 = ComputeWeight(&A, A.Pstart[1], A.Pstart[2]-1, NULL, &Options);

    if (weight0 != weightlo || weight1 != weighthi || weight0 < 0 || weight1 < 0) {
        printf("Error\n");
        exit(1);
    }

    printf("OK\n");
    exit(0);

} /* end main */
