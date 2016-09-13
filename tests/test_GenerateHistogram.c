#include "Options.h"
#include "DistributeVecLib.h"

int main(int argc, char **argv) {

    long k, n, p, i, *ProcHistogram ;
    int *X, q;

    printf("Test GenerateHistogram: ");
    k=7;
    p=15 ;
    n=k*p ;

    X  = (int *) malloc(n * sizeof(int)) ;
    ProcHistogram  = (long *) malloc((p+1) * sizeof(long)) ;
    if ( X == NULL || ProcHistogram == NULL ){
        printf("Error\n") ;
        exit(1);
    }

    for (i=0; i<n; i++)
        X[i] = i%p ;

    GenerateHistogram(X, n, 0, p, ProcHistogram ) ;
    
    /* Check result values */
    for (q=0; q<p; q++){ 
        if ( ProcHistogram[q] != k){
            printf("Error\n") ;
            exit(1);
        }
    }
    if ( ProcHistogram[p] != 0){
        printf("Error\n") ;   
        exit(1);
    }

    printf("OK\n") ;
    exit(0);

} /* end main */
  
