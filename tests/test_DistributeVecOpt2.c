#include "DistributeVecOpt2.h"

int main(int argc, char **argv) {

    struct sparsematrix A;
    long P, k, n, q, r, i, t, nzp, maxout, 
         ComVol, MaxOut, MaxIn, MaxCompnts, TotCompnts;
    int success;
    long int *X;

    printf("Test DistributeVecOpt2: ");
    P = 6; /* number of  processors >= 2 */
    k = 7; /* number of columns each pair of processors share */
    n= P*P*k; /* 2 by n matrix A  */
    A.m = 2;
    A.n = n;
    A.NrNzElts = P*(P-1)*k; 
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
    for (q=0; q<P; q++){
        /* nonzeros of proc q in row 0 */
        for (r=q+1; r<P; r++){
            for (i=0; i<k; i++){
                A.i[t] = 0;
                A.j[t] = (q*P + r) * k + i;
                t++;
            }
        }
        /* nonzeros of proc q in row 1 */
        for (r=0; r<q; r++){
            for (i=0; i<k; i++){
                A.i[t] = 1;
                A.j[t] = (r*P + q) * k + i;
                t++;
            }
        }
    }

    nzp = (P-1)*k; 
    /* Procs 0, 1, ..., P-1 have nzp nonzeros each. */
    for (i=0; i<P; i++)
        A.Pstart[i] = i*nzp;
    A.Pstart[P] = P*(P-1)*k;

    success = DistributeVecOpt2(&A, X, ROW);
    
    if (!CalcCom(&A, X, ROW, &ComVol, &MaxOut, &MaxIn, &MaxCompnts, &TotCompnts)) {
        printf("Error\n");
        exit(1);
    }
    
    if (((P-1)*k) %2 == 0)
        maxout = ((P-1)*k)/ 2; /* ceiling */
    else 
        maxout = ((P-1)*k)/ 2 + 1;

    /* Check result values  */
    if (success != TRUE || ComVol != (P*(P-1)*k)/2 || MaxOut != maxout ||
        MaxIn != MaxOut || TotCompnts != n){
        printf("Error\n");
        exit(1);
    }

    printf("OK\n");
    exit(0);

} /* end main */
