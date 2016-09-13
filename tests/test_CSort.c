#include "Sort.h"

int main(int argc, char **argv) {

    long q, n, lo, hi, i, blocknr, blocksize, nr_in_block, *J,
         *val, maxval, P;

    printf("Test CSort: ");
    P = 13 ; /* values are from 0..P-1 */
    q = 233 ;
    n = 4*q*P; /* length of J array, n is multiple of 4P.
                  J contains the numbers 0..n-1 */

    lo = n/4 ;
    hi = (3*n)/4 -1; /* lo..hi contains n/2 numbers */
    maxval = P-1 ;

    val  = (long *) malloc(n * sizeof(long)) ;
    J  = (long *) malloc(n * sizeof(long)) ;
    if ( val == NULL || J == NULL){
        printf("Error\n") ;
        exit(1);
    }

    for (i=0; i<n; i++){
        val[i] = i%P ;
        J[i] = i;
    }
    
    /* Sort by increasing value */
    CSort(J, val, maxval, lo, hi) ;

    /* Check result values */
    for (i=0; i<lo; i++){ 
        if (J[i] != i){
            printf("Error\n") ;
            exit(1);
        }
    }
    for (i=hi+1; i<n; i++){ 
        if (J[i] != i){
            printf("Error\n") ;
            exit(1);
        }
    }
    blocksize = n/(2*P) ;
    for (i=lo; i<=hi; i++){ 
        blocknr = (i-lo)/blocksize ;
        nr_in_block = (i-lo)%blocksize ;
        if (J[i] != lo + nr_in_block*P + blocknr){
            printf("Error\n") ;
            exit(1);
        }
    }
    
    printf("OK\n") ;
    exit(0);

} /* end main */
  
