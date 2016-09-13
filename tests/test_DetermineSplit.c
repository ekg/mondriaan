#include "DistributeMat.h"

struct opts Options ;

extern int DetermineSplit(long *weight, long maxweight, int k, int *procs,
                     int *isplit, long *weightlo, long *weighthi, const struct opts *pOptions);

int main(int argc, char **argv) {

    long w, *weight, maxweight, weightlo, weighthi ;
    int k, r, a, b, i, *procs ;

    printf("Test DetermineSplit: ");

    /* k, r, a, can be varied, but b is hardwired */
    k = 20; /* number of parts >= 1 */
    r = 10 ;
    w = 5*r; /* average weight per processor */
    a = 4 ; /* size of smallest part, measured in w */
    b = 21 ; /* size of largest part, measured in w */

    procs  = (int *) malloc(k * sizeof(int)) ;
    weight  = (long *) malloc(k * sizeof(long)) ;
    if ( procs == NULL || weight == NULL){
        printf("Error\n") ;
        exit(1);
    }

    weight[0] = b*w; /* largest part */
    procs[0] = b; 
    for (i=1; i<k; i++){
        weight[i] = a*w;
        procs[i] = a; 
    }

    maxweight = 2*w;
    
    /* Compute the number of processors for each part */
    Options.LoadbalanceStrategy = Constant ;
    DetermineSplit(weight, maxweight, k, procs, &i, &weightlo, &weighthi, &Options) ;

    /* Check result values */
    if ( i != 0 || (weightlo != 60*r && weightlo != 60*r -1)
                || (weighthi != 66*r && weighthi != 66*r -1) ) {
        printf("Error\n") ;
        exit(1);
    }
    
    printf("OK\n") ;
    exit(0);

} /* end main */
  
