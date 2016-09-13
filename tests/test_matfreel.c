#include "Matalloc.h"

int main(int argc, char **argv) {

    long **A, m, n, it, niters, i, j;

    printf("Test matfreel: ");
    m = 999 ;
    n = 101 ;
    niters = 17 ; 
    
    for (it=0; it<niters; it++){
        A = matallocl(m,n) ; 

        if ( A == NULL ){
            printf("Error\n") ;
            exit(1);
        }

        for (i=0; i<m; i++){
            for (j=0; j<n; j++){
                A[i][j] = i + j;
            }
        }
    
        matfreel(A) ;
    }
    
    printf("OK\n") ;
    exit(0);

} /* end main */
  
