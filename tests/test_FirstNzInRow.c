#include "Matalloc.h"
#include "DistributeVecOpt2.h"

int main(int argc, char **argv) {

    long **A, n, i, j, first;

    printf("Test FirstNzInRow: ");
    n=200 ; 
    
    A = matallocl(n,n) ; 

    if ( A == NULL ){
        printf("Error\n") ;
        exit(1);
    }

    /* A is an upper triangular matrix */
    for (i=0; i<n; i++){
        for (j=0; j<n; j++){
            if (i <= j)
                A[i][j] = -1;
            else
                A[i][j] = 0 ;
        }
    }
    
    /* In the lower triangular part and on the diagonal, 
       the first nonzero element in each row is the diagonal element.
       In the upper triangular part, the first nonzero element is 
       the element itself. */
 
    for (i=0; i<n; i++){
        for (j=0; j<n; j++){
            /* first element in row i with index >= j  */
            first = FirstNzInRow(A,n,i,j);
            if ( (i >= j && first != i) ||
                 (i <  j && first != j)) {
                printf("Error\n") ;
                exit(1);
            }
        }
    }
    
    matfreel(A) ;

    printf("OK\n") ;
    exit(0);

} /* end main */
  
