#include "DistributeVecLib.h"

int main(int argc, char **argv) {

    long *Histogram ;
    int P, q ;

    printf("Test PrintHistogram: ");
    P = 68 ; /* number of  processors >= 1 */

    Histogram = (long *) malloc((P+1)* sizeof(long)) ;

    if ( Histogram == NULL ){
        printf("Error\n") ;
        exit(1);
    }

    /* Initialise histogram*/
    for (q=0; q<P/2; q++) 
        Histogram[q] = q ;
    for (q=P/2; q<=P; q++) 
        Histogram[q] = 0 ;

    printf("\n");
    printf("*** Histogram for processor q should be q\n");
    PrintHistogram(0, P, Histogram) ;

    printf("OK\n") ;
    exit(0);

} /* end main */
