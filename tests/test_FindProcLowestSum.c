#include <stdlib.h>
#include <stdio.h>

#include "DistributeVecOrigEq.h"

struct opts Options ;

extern int FindProcLowestSum(long *procstart, int *procindex, long j, int npj,
                       long *SumsU, long *SumsV);

int main(int argc, char **argv) {

    long i, l, *procstart, *SumsU, *SumsV ;
    int P, npi, q, q0, q1, *procindex ;

    printf("Test FindProcLowestSum: ");
    P = 23 ;   /* number of non-empty processors */
    l= P*P*P ;
    npi = 2;

    procstart = (long *) malloc((l+1)* sizeof(long)) ;
    procindex = (int *) malloc(2*l* sizeof(int)) ;
    SumsU = (long *) malloc(P* sizeof(long)) ;
    SumsV = (long *) malloc(P* sizeof(long)) ;

    if ( procstart == NULL || procindex == NULL ||
         SumsU == NULL || SumsV == NULL ){
        printf("Error\n") ;
        exit(1);
    }

    /* Fill rows of communication matrix, each with two processors */
    for (i=0; i<l; i++)
        procstart[i] = 2*i;
    procstart[l] = 2*l ;
    for (i=0; i<l; i++){
        procindex[2*i] = i%P ;
        procindex[2*i+1] = (i+1)%P ;
    }

    /* Initialise sums */
    for (q=0; q<P; q++){  
        SumsU[q] = q;
        SumsV[q] = q;
    }

    /* Check result values */
    for (i=0; i<l; i++){ 
        q = FindProcLowestSum(procstart, procindex, i, npi, SumsU, SumsV ) ;
        q0 = i%P ;
        q1 = (i+1)%P ;

        if(q != MIN(q0,q1) ){
            printf("Error\n") ;
            exit(1);
        }
    }

    printf("OK\n") ;
    exit(0);

} /* end main */
