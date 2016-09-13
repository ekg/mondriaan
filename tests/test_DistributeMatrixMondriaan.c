#include "DistributeMat.h"

struct opts Options;

int main(int argc, char **argv) {

    struct sparsematrix A;
    long n, neps, nz, i, j, iter, ieps, t,
         maxweight, totweight, weight;
    int P, Pmax, q;
    double eps;

    MMSparseMatrixInit( &A );

    printf("Test DistributeMatrixMondriaan: ");
    n = 1000; /* n by n symmetric tridiagonal matrix A, n even*/
    neps = 4;  /* number of different epsilons tested */
    A.m = n;
    A.n = n;
    A.NrNzElts = 2*n-1;
    Pmax = 32; /* maximum number of parts */

    A.i = (long *) malloc(A.NrNzElts* sizeof(long));
    A.j = (long *) malloc(A.NrNzElts* sizeof(long));
    A.Pstart = (long *) malloc((Pmax+1)* sizeof(long));
    A.RowLambda = (int *)malloc(A.m*sizeof(int));
    A.ColLambda = (int *)malloc(A.n*sizeof(int));
    A.RowMark = (int *)malloc(A.m*sizeof(int));
    A.ColMark = (int *)malloc(A.n*sizeof(int));
    A.dummy = (int *) malloc(n* sizeof(int));

    if (A.i == NULL || A.j == NULL || A.Pstart == NULL || A.dummy == NULL) {
        printf("Error\n");
        exit(1);
    }

    /* Initialise matrix */
    A.MMTypeCode[3]='S';
    Options.SymmetricMatrix_UseSingleEntry = SingleEntYes;
    A.NrDummies = n/2;
    for (i=0; i<n/2; i++)
        A.dummy[i] = TRUE;
    for (i=n/2; i<n; i++)
        A.dummy[i] = FALSE;
 
    /* Fill matrix with nonzeros */
    t= 0;
    for (i=0; i<n; i++) {
        for (j=MAX(0,i-1); j<=i; j++) {
            A.i[t] = i;
            A.j[t] = j;
            t++;
        }
    }    
    totweight = ComputeWeight(&A, 0, A.NrNzElts-1, NULL, &Options);

    /* Main loop for testing different options */
    Options.SplitStrategy = OneDimRow; /* does not matter, but needs a value */
    Options.SplitMethod = Simple;

    for (iter=0; iter<6; iter++) {
        if (iter%2 == 0)
            Options.LoadbalanceAdjust = AdjustYes;
        else
            Options.LoadbalanceAdjust = AdjustNo;
        if ((iter%6)/2 == 0)
            Options.LoadbalanceStrategy = Decrease;
        else if ((iter%6)/2 == 1)
            Options.LoadbalanceStrategy = Increase;
        else if ((iter%6)/2 == 2)
            Options.LoadbalanceStrategy = Constant;

        for (P=1; P<=Pmax; P++) {
            eps = 2*P*P / (double)totweight;
            for (ieps=0; ieps<neps; ieps++) {
                A.Pstart[0] = 0;
                for (q=1; q<=P; q++)
                    A.Pstart[q] = A.NrNzElts;
                
                if (!DistributeMatrixMondriaan(&A, P, eps, &Options, 0)) {
                    printf("Error\n");
                    exit(1);
                }

                /* Check that all parts have weight <= maxweight */
                maxweight = ((1 + eps) * totweight) / P; /* rounded down */
                nz= 0;
                for (q=0; q<P; q++) {
                    weight = ComputeWeight(&A, A.Pstart[q],A.Pstart[q+1]-1,NULL,&Options);
                    if (weight > maxweight || weight < 0) {
                        printf("Error: weight too large\n");
                        exit(1); 
                    }
                    nz += A.Pstart[q+1] - A.Pstart[q];
                }
                if (nz != A.NrNzElts) {
                    printf("Error\n"); 
                    exit(1);  
                }

                /* It is assumed that the simple split method does not
                   move nonzeros. Check this */
                t= 0;
                for (i=0; i<n; i++) { 
                    for (j=MAX(0,i-1); j<=i; j++) {
                        if (A.i[t] != i || A.j[t] != j) { 
                            printf("Error\n");  
                            exit(1);     
                        }  
                        t++;
                    }
                }
                eps *= 2;
            }
        }
    } 

    printf("OK\n");
    exit(0);

} /* end main */
