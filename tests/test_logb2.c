#include <stdlib.h>
#include <stdio.h>

#include "Options.h"
#include "DistributeMat.h"

struct opts Options ;

extern int logb2(int n);

int main(int argc, char **argv) {

    int k, m, n, i, res ;

    printf("Test logb2: ");
    n = 12345 ;

    m = 0;
    for (k=1; k<n; k *= 2){
        /* k = 2^m */
        for (i=0; i<k; i++) {
            res = logb2(k+i) ;
            if ((i==0 && res != m ) || (i>0 && res != m+1)) {
                printf("Error\n") ;
                exit(1);
            }
        }
        m++ ;
    }
    
    printf("OK\n") ;
    exit(0);

} /* end main */
  
