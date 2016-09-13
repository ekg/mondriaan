#include "DistributeVecGreedy.h"

struct opts Options;

int main(int argc, char **argv) {

    struct sparsematrix A;
    long P, k, m, n, q, r, i, j, t, nzp, 
         ComVol, MaxOut, MaxIn, MaxCompnts, TotCompnts;
    long int *X;

    printf("Test DistributeVecGreedyImprove: ");

    k = 13; /* number of processors in a column k>=2 */
    m = 2*k - 1;
    q = 100; /* number of m x m matrices in A */
    n= q*m; /* m by n matrix A  */
    P = m; /* number of  processors  */
    A.m = m;
    A.n = n;
    A.NrNzElts = k*n; 
    A.NrProcs = P;

    A.i = (long *) malloc(A.NrNzElts* sizeof(long));
    A.j = (long *) malloc(A.NrNzElts* sizeof(long));
    A.Pstart = (long *) malloc((P+1)* sizeof(long));
    X = (long int *) malloc(n* sizeof(long int));

    if ( A.i == NULL || A.j  == NULL || A.Pstart == NULL || X == NULL ){
        printf("Error\n");
        exit(1);
    }

    /* Fill matrix with nonzeros */
    t= 0;
    for (i=0; i<P; i++){
        /* nonzeros of proc i in row i */
        for (r=0; r<q; r++){
            /* nonzeros of matrix r */
            for (j=0; j<m; j++){
                if ( (0 <= i - j && i - j < k  ) ||
                     (0 >  i - j && i - j <= -k ) ){
                    A.i[t] = i;
                    A.j[t] = r*m + j;
                    t++;
                }
            }
        }
    }

    nzp = k*n/P; 
    /* Procs 0, 1, ..., P-1 have nzp nonzeros each. */
    for (i=0; i<P; i++)
        A.Pstart[i] = i*nzp;
    A.Pstart[P] = k*n;

    /* Assign all columns initially to P(0) */
    for (j=0; j<n; j++)
        X[j] = 0;
    
    /* Guarantee optimal solution. 
       Each loop iteration decreases the cost at least by 1.
       The initial cost is <= n*k. */
    Options.VectorPartition_MaxNrGreedyImproves = n*k;

    if (DistributeVecGreedyImprove(&A, X, ROW, &Options) < 0) {
        printf("Error\n");
        exit(1);
    }
    if (!CalcCom(&A, X, ROW, &ComVol, &MaxOut, &MaxIn, &MaxCompnts, &TotCompnts)) {
        printf("Error\n");
        exit(1);
    }

    /* Check result values  */
    if (ComVol != n*(k-1) || MaxOut != q*(k-1) ||
        MaxIn != MaxOut || TotCompnts != n){
        printf("Error\n");
        exit(1);
    }

    printf("OK\n");
    exit(0);

} /* end main */
