#include "Sort.h"

int main(int argc, char **argv) {

    long *X, n, i;

    printf("Test SwapLong: ");
    n=1000 ; /* must be even */

    X  = (long *) malloc(n* sizeof(long)) ;
    if ( X == NULL ){
        printf("Error\n") ;
        exit(1);
    }

    for (i=0; i<n; i++)
        X[i] = i ;
    for (i=0; i<n; i +=2)
        SwapLong(X,i,i+1) ;
    
    /* Check result values */
    for (i=0; i<n; i +=2){ 
        if (X[i] != i+1 || X[i+1] != i) {
            printf("Error\n") ;
            exit(1);
        }
    }

    printf("OK\n") ;
    exit(0);

} /* end main */
  
