#include "Sort.h"

int main(int argc, char **argv) {

    long *X, *Tmp, n, i;
    double *D;

    printf("Test RandomPermute: ");
    n=1000 ;

    X  = (long *) malloc(n* sizeof(long)) ;
    Tmp = (long *) malloc(n* sizeof(long)) ;
    D  = (double *) malloc(n* sizeof(double)) ;
    if ( X == NULL || Tmp  == NULL || D == NULL){
        printf("Error\n") ;
        exit(1);
    }

    for (i=0; i<n; i++){
        X[i] = i ;
        D[i] = i;
    }
#ifdef INFO2
    printf("RAND_MAX = %ld ", (long)(RAND_MAX) );
#endif
    
    /* Permute randomly */
    RandomPermute(X,NULL,D,NULL,0,n-1);

    /* Check result values and corresponding indices */
    for (i=0; i<n; i++){ 
        if (X[i] != (long)(D[i]+0.5) || X[i] < 0 || X[i] >= n){
          /* round to nearest integer */
            printf("Error\n") ;
            exit(1);
        }
    }

    /* Check whether X is a permutation of 0..n-1.
       Use Tmp as a temporary array */ 
    for (i=0; i<n; i++) 
        Tmp[i] = -1; 

    for (i=0; i<n; i++){ 
        Tmp[X[i]]= X[i];
    }

    for (i=0; i<n; i++){  
        if (Tmp[i] != i){
            printf("Error\n") ; 
            exit(1); 
        } 
    }

    printf("OK\n") ;
    exit(0);

} /* end main */
  
