#include <stdlib.h>
#include <stdio.h>

#include "DistributeMat.h"

struct opts Options ;

extern int BalanceParts(long *weight, long maxweight, int P, int k, int *procs);

int main(int argc, char **argv) {

    long maxweight, totweight, *weight ;
    int P, k, q, r, a, b, i, *procs, procs0 ;

    printf("Test BalanceParts: ");
    q = 13 ;
    k = q*3; /* number of parts */
    r = 17;
    a = 5*r; /* weight increment */
    b = 3*r; /* average weight per processor */ ;

    procs  = (int *) malloc(k * sizeof(int)) ;
    weight  = (long *) malloc(k * sizeof(long)) ;
    if ( procs == NULL || weight == NULL){
        printf("Error\n") ;
        exit(1);
    }

    for (i=0; i<k; i++)
        weight[i] = a*(i+1) ;
    totweight = a*k*(k+1)/2 ;
    P = totweight/b;
    maxweight = a-1 ; /* all parts must be split */
    
    /* Compute the number of processors for each part */
    BalanceParts (weight, maxweight, P, k, procs) ;

    /* Check result values */
    for (i=0; i<k; i++){ 
        procs0= (a*(i+1))/b ; /* number of procs rounded down */
        if ( (i%3==0 && procs[i] != procs0+1) ||
             (i%3==1 && procs[i] != procs0) ||
             (i%3==2 && procs[i] != procs0) ){
            printf("Error\n") ;
            exit(1);
        }
    }
    
    printf("OK\n") ;
    exit(0);

} /* end main */
  
