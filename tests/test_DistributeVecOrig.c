#include "DistributeVecOrig.h"

struct opts Options;

int main(int argc, char **argv) {

    struct sparsematrix A;
    long P, n, i, t, nzp, maxcom,
         ComVol, MaxOut, MaxIn, MaxCompnts, TotCompnts;
    long int *X;

    printf("Test DistributeVecOrig: ");
    P = 666; /* number of  processors >= 3 */

    /* P by 3P matrix A, consisting of three P by P blocks,
       A = [B|I|0] where B is a circulant matrix with three nonzero
       diagonals. */
    A.m = P;
    A.n = 3*P;
    n = A.n;
    A.NrNzElts = 4*P; 
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
        /* Insert the 4 nonzeros of a row */
        A.i[t] = A.i[t+1] = A.i[t+2] = A.i[t+3] = i; 
        A.j[t] = i;
        A.j[t+1] = (i+P-1) % P;
        A.j[t+2] = (i+P-2) % P;
        A.j[t+3] = i+P;
        t += 4;
    }

    nzp = 4; 
    /* Procs 0, 1, ..., P-1 have nzp nonzeros each. */
    for (i=0; i<P; i++)
        A.Pstart[i] = i*nzp;
    A.Pstart[P] = 4*P;

    Options.VectorPartition_Step3 = VecRandom;
    maxcom = DistributeVecOrig(&A, X, ROW, &Options);
    if (!maxcom < 0) {
        printf("Error\n");
        exit(1);
    }
    if (!CalcCom(&A, X, ROW, &ComVol, &MaxOut, &MaxIn, &MaxCompnts, &TotCompnts)) {
        printf("Error\n");
        exit(1);
    }

    /* Check result values. This is not a very strict check,
       because vector assignment and maxcom can vary. */
    if (ComVol != 2*P || maxcom !=  MAX(MaxIn,MaxOut) || TotCompnts != n){
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

    /* Check first block */ 
    for (i=0; i<P; i++){
        if (X[i] != i && X[i] != (i+1)%P && X[i] != (i+2)%P ) { 
            printf("Error\n"); 
            exit(1);
        } 
    }    

    /* Check second block */
    for (i=P; i<2*P; i++){
        if (X[i] != i-P) {
            printf("Error\n");
            exit(1);
        }
    }   


    printf("OK\n");
    exit(0);

} /* end main */
