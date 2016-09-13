#include "DistributeVecOrigEq.h"

struct opts Options ;

int main(int argc, char **argv) {

    struct sparsematrix A ;
    long P, n, k, i, maxcom;
    long int *U, *V;

    printf("Test DistributeVecOrigEq: ");
    P = 6 ; /* number of  processors */
    k= 100 ;

    /* k*P by k*P identity matrix A.
       Diagonal element a[i,i] is owned by processor i div k. */
    A.m = k*P;
    A.n = k*P;
    n = A.n ;
    A.NrNzElts = n ; 
    A.NrProcs = P ;

    A.i = (long *) malloc(A.NrNzElts* sizeof(long)) ;
    A.j = (long *) malloc(A.NrNzElts* sizeof(long)) ;
    A.Pstart = (long *) malloc((P+1)* sizeof(long)) ;
    U = (long int *) malloc(n* sizeof(long int)) ;
    V = (long int *) malloc(n* sizeof(long int)) ;

    if ( A.i == NULL || A.j  == NULL || A.Pstart == NULL ||
         U == NULL || V == NULL ){
        printf("Error\n") ;
        exit(1);
    }

    /* Fill matrix with nonzeros */
    for (i=0; i<n; i++){
        /* Insert the diagonal nonzero of row i */
        A.i[i] = i;
        A.j[i] = i ;
    }

    /* Procs 0, 1, ..., P-1 have k nonzeros each. */
    for (i=0; i<P; i++)
        A.Pstart[i] = i*k ;
    A.Pstart[P] = n ;

    Options.SplitStrategy = LocalBest ;
    maxcom = DistributeVecOrigEq(&A, U, V, &Options) ;

    /* Check result value */
    if (maxcom !=  0){
        printf("Error\n") ;
        exit(1);
    }

    /* Check vector distributions */
    for (i=0; i<n; i++){
        if (U[i] != i/k || V[i] != U[i]) {
            printf("Error\n") ;
            exit(1);
        }
    }

    printf("OK\n") ;
    exit(0);

} /* end main */
