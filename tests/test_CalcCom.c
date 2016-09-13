#include "DistributeVecLib.h"

int main(int argc, char **argv) {

    struct sparsematrix A;
    long P, k, q, n, i, j, t, nzp,
         ComVol, MaxOut, MaxIn, MaxCompnts, TotCompnts;
    long int *X;

    printf("Test CalcCom: ");
    P = 8; /* number of non-empty processors */
    k= 10*P;
    q= 4;
    n= q*k; /* n by n matrix A with nonzeros in positions
                A[i,j] with i mod q = j mod q = 0 */
    A.m = n;
    A.n = n;
    A.NrNzElts = k*k; 
    A.NrProcs = P;

    A.i = (long *) malloc(A.NrNzElts* sizeof(long));
    A.j = (long *) malloc(A.NrNzElts* sizeof(long));
    A.Pstart = (long *) malloc((P+1)* sizeof(long));
    X = (long int *) malloc(n* sizeof(long int));

    if (A.i == NULL || A.j  == NULL || A.Pstart == NULL || X == NULL){
        printf("Error\n");
        exit(1);
    }

    /* Fill matrix with k*k nonzeros:
       k nonempty rows, each with k nonzeros */ 
    t= 0;
    for (i=0; i<k; i++){
        for (j=0; j<k; j++){
            A.i[t] = i*q;
            A.j[t] = j*q;
            t++;
        }
    }

    nzp = k*k/P;
    /* Procs 0, 1, ..., P-1 have nzp nonzeros each. */
    for (i=0; i<P; i++)
        A.Pstart[i] = i*nzp;
    A.Pstart[P] = k*k;
    for (j=0; j<n; j++)
        X[j] = (j/q)%P;

    if (!CalcCom(&A, X, ROW, &ComVol, &MaxOut, &MaxIn, &MaxCompnts, &TotCompnts)) {
        printf("Error\n");
        exit(1);
    }

    /* Check result values and corresponding indices */
    if (ComVol != n*(P-1)/q || MaxOut != n*(P-1)/(q*P) ||
        MaxIn != MaxOut || MaxCompnts != n/P || TotCompnts != n){
        printf("Error\n");
        exit(1);
    }

    printf("OK\n");
    exit(0);

} /* end main */
