#include "DistributeVecLocal.h"

int main(int argc, char **argv) {

    struct sparsematrix A;
    long P, n, i, j, t, nzp, maxcom,
         ComVol, MaxOut, MaxIn, MaxCompnts, TotCompnts;
    long int *X;

    printf("Test DistributeVecLocal: ");
    P = 66; /* number of  processors >= 2 */
    n = P; /* P by P dense matrix A  */
    A.m = n;
    A.n = n;
    A.NrNzElts = P*P; 
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
        for (j=0; j<P; j++){
            A.i[t] = i;
            A.j[t] = j;
            t++;
        }
    }

    nzp = P; 
    /* Procs 0, 1, ..., P-1 have nzp nonzeros each. */
    for (i=0; i<P; i++)
        A.Pstart[i] = i*nzp;
    A.Pstart[P] = P*P;

    maxcom = DistributeVecLocal(&A, X, ROW);
    if (!CalcCom(&A, X, ROW, &ComVol, &MaxOut, &MaxIn, &MaxCompnts, &TotCompnts)) {
        printf("Error\n");
        exit(1);
    }
    
    /* Check result values  */
    if (ComVol != P*(P-1) || MaxOut != P-1 ||
        MaxIn != MaxOut || MaxIn != maxcom || TotCompnts != n){
        printf("Error\n");
        exit(1);
    }

    /* Check legality of vector distribution */
    for (i=0; i<n; i++){
        if (X[i] < 0 || X[i] >= P) {
            printf("Error\n");
            exit(1);
        }
    }

    printf("OK\n");
    exit(0);

} /* end main */
