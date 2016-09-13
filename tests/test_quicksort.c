#include "Sort.h"

int main(int argc, char **argv) {

    long *X, *I, k, r, n, i, j;

    printf("Test quicksort: ");
    r = 100 ;
    k = 3 ;
    n = k*r;

    X  = (long *) malloc(n * sizeof(long)) ;
    I  = (long *) malloc(n * sizeof(long)) ;

    if ( X == NULL || I == NULL ){
        printf("Error\n") ;
        exit(1);
    }

    for (i=0; i<n; i++){
        X[i] = i%k ;  /* value */
        I[i] = i;     /* index */
    }
        
    /* Sort by decreasing value */
    quicksort(X, I, 0, n-1) ;

    /* Check result values and corresponding indices.
       They should be ordered as: 
       r values k-1, r values k-2, ..., r values 0. */
    for (j=0; j<k; j++){ 
        for (i=0; i<r; i++){ 
            if (X[j*r+i] != k-j-1 || I[j*r+i]%k != X[j*r+i] ||
                I[j*r+i] < 0 || I[j*r+i] >= n ) {
                printf("Error\n") ;
                exit(1);
            }
        }
    }

    /* Check whether I is a permutation of 0..n-1.
       Use X as a temporary array */ 
    for (i=0; i<n; i++)
        X[i] = -1; 

    for (i=0; i<n; i++)
        X[I[i]]= I[i] ;

    for (i=0; i<n; i++){  
        if (X[i] != i){
            printf("Error\n") ; 
            exit(1); 
        } 
    }

    printf("OK\n") ;
    exit(0);

} /* end main */
  
