#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "Sort.h"

int main(int argc, char **argv) {

    /* This test function uses the Chernoff bound (from 1952) for tails of 
       the binomial distribution. It is given as theorem 4.11 on p. 205 of 
       Parallel Scientific Computation: A structured approach 
       using BSP and MPI by Rob H. Bisseling (Oxford University Press 2004).
    */

    long X, i, m, r, r1, r2 ;
    double d, eps, factor, Prob;

    printf("Test Random1: ");
    m=256 ; /* Keep this fixed. Chances of succes are about
               1 in a million.  */
    d = 0.5 ;
    eps = 0.5 ;
    factor = exp(eps) / pow(1+eps,1+eps) ;
    Prob = pow (factor, m*d) ;

#ifdef INFO2
    printf("RAND_MAX = %ld\n", (long)(RAND_MAX) );
    printf("Probability = %g\n", Prob );
#endif

    SetRandomSeed(99) ; /* Deterministic, for reproducibility of test */
    
    /* Test boundaries */
    r1 = Random1(1,0) ; /* lo > hi */
    r2 = Random1(2,2) ; /* lo = hi */
    if (r1 != LONG_MAX || r2 != 2 ) {
        printf("Error\n") ;
        exit(1);
    } 
    
    /* Perform m independent Bernoulli trials */
    X = 0 ;
    for (i=0; i<m; i++) {
        /* Flip a coin with success probability d = 0.5 */
        r = Random1(0,1) ;  
        if (r != 0 && r != 1 ) {
            printf("Error\n") ;
            exit(1);
        } 
        X += r ;
    }
    if ( (double)X > (1+eps) * m * d ) {
        printf("Error. You shouldn't be so lucky.\n") ;
        exit(1);
    } 
 
    printf("OK\n") ;
    exit(0);

} /* end main */
  
