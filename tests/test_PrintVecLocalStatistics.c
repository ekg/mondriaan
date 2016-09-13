#include "DistributeVecLocal.h"

int main(int argc, char **argv) {

    long P, n, j, t,
         *Ns, *Nr, *Ls, *Lr, *Jstart0, *J;
    int q;
    long *owner;

    printf("Test PrintVecLocalStatistics: ");
    P = 6 ; /* number of  processors >= 2 */
    n = P ; /* P by P lower triangular matrix A  */

    Ns = (long *) malloc(P* sizeof(long)) ;
    Nr = (long *) malloc(P* sizeof(long)) ;
    Ls = (long *) malloc(P* sizeof(long)) ;
    Lr = (long *) malloc(P* sizeof(long)) ;
    Jstart0 = (long *) malloc((P+1)* sizeof(long)) ;
    J = (long *) malloc((n*(n+1)/2)* sizeof(long)) ;
    owner = (long *)malloc(n*sizeof(long)) ;

    if ( Ns == NULL || Nr  == NULL || Ls == NULL || Lr  == NULL ||
         Jstart0 == NULL || J == NULL || owner == NULL ){
        printf("Error\n") ;
        exit(1);
    }

    /* Initialise statistics */
    for (q=0; q<P; q++) {
        Ns[q] = q ;
        Nr[q] = 2*q ;
        Ls[q] = 3*q ;
        Lr[q] = 4*q ;
        owner[q]= -1 ;
    }
    /* Processor q is present in columns 0, 1, ..., q.
       These columns are all unowned. */
    Jstart0[0] = 0;
    t= 0;
    for (q=0; q<P; q++) { 
        for (j=0; j<=q ; j++){
            J[t] = j ;
            t++ ;
        }
        Jstart0[q+1] = Jstart0[q] + q+1;
    }

    printf("\n");
    printf("*** Statistics for processor q should be\n");
    printf("    Ns = q, Nr = 2q, Ls = 3q, Lr = 4q, Un = q+1 \n");
    PrintVecLocalStatistics(P, Ns, Nr, Ls, Lr, Jstart0, J, owner);

    printf("OK\n") ;
    exit(0);

} /* end main */
