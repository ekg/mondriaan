#include "Matalloc.h"

int main(int argc, char **argv) {

    long **A, n, i, j, sum;

    printf("Test matallocl: ");
    n=10 ; /* small, because we compute a value n^4 */
    
    A = matallocl(n,n) ; 

    if ( A == NULL ){
        printf("Error\n") ;
        exit(1);
    }

    for (i=0; i<n; i++){
        for (j=0; j<n; j++){
            A[i][j] = i + j*n;
        }
    }
    
    /* Compute sum */
    sum = 0;
    for (i=0; i<n; i++){
        for (j=0; j<n; j++){
            sum += A[i][j] ;
        }
    }
    
    /* Check sum of array elements */
    if (sum != n*n*(n*n-1)/2 ){
        printf("Error\n") ;
        exit(1);
    }

    printf("OK\n") ;
    exit(0);

} /* end main */
  
