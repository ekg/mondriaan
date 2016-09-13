#include "DistributeVecLib.h"
#include "DistributeVecOrig.h"

struct opts Options;

int main(int argc, char **argv) {

    long P, P2, k, n, i, j, *procstart, *Sums;    
    int q, *procindex;

    printf("Test InitSums: ");
    P = 12; /* P is the number of processors, must be even */
    P2 = P/2;
    k= 12; /* must be even */
    n= P*k; /* P by n communication matrix C with processors
                in positions 0,1,...,P/2-1 in column j if j is even and j<n/2,
                and in positions P/2,...,P-1 in column j if j is odd and j<n/2.
             */

    procstart = (long *)malloc((n+1)*sizeof(long));
    procindex = (int *)malloc((n/2)*P2*sizeof(int));
    Sums = (long *)malloc(P*sizeof(long));

    if ( procstart == NULL || procindex == NULL || Sums ==  NULL ) {
        printf("Error\n");
        exit(1);
    }

    /* Initialise procstart */
    /* Half-full columns */
    for (j=0; j<n/2; j++) 
        procstart[j]=P2*j;
    /* Empty columns */
    for (j=n/2; j<=n; j++) 
        procstart[j]=P2*n/2;

    /* Fill procindex cyclically */
    for (i=0; i<P2*n/2; i++)
        procindex[i] = i%P;

    if (!InitSums(n, P, procstart, procindex, Sums)) {
        printf("Error\n");
        exit(1);
    }

    /* Check sums */
    for (q=0; q<P; q++) { 
        if (Sums[q] != n/4 ) {
            printf("Error\n");
            exit(1);
        }
    }
 
    printf("OK\n");
    exit(0);

} /* end main */
