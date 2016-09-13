#include "Cartesian.h"

int main(int argc, char **argv) {

    long k, m, n, nz, i, j, t, *start, *index ;

    printf("Test CRS2CCS: ");
    n = 12 ; 
    k = 7 ; /* k copies of n by n identity matrix I
                stacked above each other */
                
    m = k*n ;
    nz = m ;

    start = (long *) malloc((m+1)* sizeof(long)) ;
    index = (long *) malloc(nz* sizeof(long)) ;
 
    if ( start == NULL || index == NULL ){
        printf("Error\n") ;
        exit(1);
    }

 
    /* Initialise start, including start[m] */
    for (i=0; i<=m; i++)
        start[i] = i ;
            
    /* Initialise index */
    for (t=0; t<nz; t++)
        index[t] = t%n ;

    /* Convert from CRS to CCS */
    CRS2CCS(m, n, nz, start, index) ;

    /* Check start values for columns */
    for (j=0; j<=n; j++)
        if (start[j] != k*j){
            printf("Error\n") ;
            exit(1);
        }  
    
    /* Check index values for columns */
    for (j=0; j<n; j++)
        for (t=0; t<k; t++)
            if (index[k*j+t] != j + t*n){
                printf("Error\n") ;
                exit(1);
            }  
            
    /* Convert back to CRS */
    CRS2CCS(n, m, nz, start, index) ;
    
    /* Check start values for rows */
    for (i=0; i<=m; i++)
        if (start[i] != i){
            printf("Error\n") ;
            exit(1);
        }  

    /* Check index values for rows */
    for (t=0; t<nz; t++)
        if (index[t] != t%n ){
            printf("Error\n") ;
            exit(1);
        }  
        
    printf("OK\n") ;
    exit(0);

} /* end main */
