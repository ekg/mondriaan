#include "Sort.h"

int main(int argc, char **argv) {

    long n, i;
    double *X, *X0;

    printf("Test SwapDouble: ");
    n=1000 ; /* must be even */

    X  = (double *) malloc(n* sizeof(double)) ;
    X0  = (double *) malloc(n* sizeof(double)) ;
    if ( X == NULL || X0 == NULL ){
        printf("Error\n") ;
        exit(1);
    }

    for (i=0; i<n; i++){
        X[i] = i ;
        X0[i] = i ; /* keep original values */
    }
    for (i=0; i<n; i +=2)
        SwapDouble(X,i,i+1) ;
    
    /* Check result values */
    for (i=0; i<n; i +=2){ 
        if (X[i] != X0[i+1] || X[i+1] != X0[i]) {
            printf("Error\n") ;
            exit(1);
        }
    }

    printf("OK\n") ;
    exit(0);

} /* end main */
  
