#include "Sort.h"

int main(int argc, char **argv) {

    long *X, *I, LengthX, n, i;

    printf("Test QSort: ");
    n=100 ;
    LengthX= 2*n;

    X  = (long *) malloc(LengthX * sizeof(long)) ;
    if ( X == NULL ){
        printf("Error\n") ;
        exit(1);
    }

    for (i=0; i<n; i++){
        X[i] = i ;
        X[n+i] = i;
    }
    
    /* Sort by decreasing value */
    I= QSort(X, LengthX);
    
    if ( I == NULL ) {
        printf("Error\n");
        exit(1);
    }

    /* Check result values and corresponding indices */
    for (i=0; i<n; i++){ 
        if (X[2*i] != n-i-1 || X[2*i+1] != n-i-1){
            printf("Error\n") ;
            exit(1);
        }
        if (I[2*i] != n-i-1 && I[2*i] != n+ n-i-1){
            printf("Error\n") ;
            exit(1);
        }
        if (I[2*i+1] != n-i-1 && I[2*i+1] != n+ n-i-1){
            printf("Error\n") ;
            exit(1);
        }
    }

    /* Check whether I is a permutation of 0..2n-1.
       Use X as a temporary array */ 
    for (i=0; i<n; i++){ 
        X[2*i] = X[2*i+1] =-1; 
    }

    for (i=0; i<n; i++){ 
        X[I[2*i]]= I[2*i];
        X[I[2*i+1]]= I[2*i+1];
    }

    for (i=0; i<2*n; i++){  
        if (X[i] != i){
            printf("Error\n") ; 
            exit(1); 
        } 
    }

    printf("OK\n") ;
    exit(0);

} /* end main */
  
