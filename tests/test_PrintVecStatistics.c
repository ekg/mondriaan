#include "DistributeVecLib.h"

int main(int argc, char **argv) {

    long *Ns, *Nr, *Nv ;
    int P, q ;

    printf("Test PrintVecStatistics: ");
    P = 16 ; /* number of  processors >= 1 */

    Ns = (long *) malloc(P* sizeof(long)) ;
    Nr = (long *) malloc(P* sizeof(long)) ;
    Nv = (long *) malloc(P* sizeof(long)) ;

    if ( Ns == NULL || Nr  == NULL || Nv == NULL ){
        printf("Error\n") ;
        exit(1);
    }

    /* Initialise statistics */
    for (q=0; q<P; q++) {
        Ns[q] = q ;
        Nr[q] = 2*q ;
        Nv[q] = 3*q ;
    }

    printf("\n");
    printf("*** Statistics for processor q should be\n");
    printf("    Ns = q, Nr = 2q, Nv = 3q \n");
    PrintVecStatistics(P, Ns, Nr, Nv) ;

    printf("OK\n") ;
    exit(0);

} /* end main */
