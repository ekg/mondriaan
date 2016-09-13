#include "DistributeMat.h"

struct opts Options;

int main(int argc, char **argv) {

    struct sparsematrix A;
    long n, i, j, t, weight, weight1, *wnz;

    printf("Test ComputeWeight: ");
    n = 4; /* n by n dense lower triangular matrix A, n even */
    A.m = n;
    A.n = n;
    A.NrNzElts = (n*(n+1))/2; 

    A.i = (long *) malloc(A.NrNzElts* sizeof(long));
    A.j = (long *) malloc(A.NrNzElts* sizeof(long));
    wnz = (long *) malloc(A.NrNzElts* sizeof(long));
    A.dummy = (int *) malloc(n * sizeof(int));

    if (A.i == NULL || A.j  == NULL || wnz == NULL || A.dummy == NULL) {
        printf("Error\n");
        exit(1);
    }

    /* Fill lower triangular matrix with nonzeros */
    t= 0;
    for (i=0; i<n; i++) {
        for (j=0; j<=i; j++) {
            A.i[t] = i;
            A.j[t] = j;
            t++;
        }
    }
    A.MMTypeCode[3]='S';
    Options.SymmetricMatrix_UseSingleEntry = SingleEntYes;

    /* First half of diagonal is dummy */
    for (i=0; i < n/2; i++)
        A.dummy[i] = TRUE;    
    for (i=n/2; i < n; i++)
        A.dummy[i] = FALSE;    
    A.NrDummies = n/2;

    weight = ComputeWeight(&A, 0, A.NrNzElts-1,wnz,&Options);
    
    if (weight < 0) {
        printf("Error\n");
        exit(1);
    }

    /* Check result value  */
    if (weight != n*n -n/2) {
        printf("Error\n");
        exit(1);
    }

    /* Check weights of individual nonzeros */
    t= 0;
    weight1= 0;
    for (i=0; i<n; i++) {
        for (j=0; j<=i; j++) {
            if ((i != j && wnz[t] != 2) ||
                (i == j && i<n/2 && wnz[t] != 0) ||
                (i == j && i>=n/2 && wnz[t] != 1)) {
                printf("Error\n");
                exit(1);
            }
            weight1 += wnz[t];
            t++;
        }
    }

    /* Check if individual weights add up to total weight */
    if (weight != weight1) {
        printf("Error\n"); 
        exit(1); 
    }


    printf("OK\n");
    exit(0);

} /* end main */
