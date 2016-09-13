#include <stdlib.h>
#include <stdio.h>

#include "DistributeVecOrigEq.h"

struct opts Options;

extern int InitIntersection(long l, int P,
                     long *procstartU, int *procindexU,
                     long *procstartV, int *procindexV, int *NprocsUV);

int main(int argc, char **argv) {

    long i, b, l, *procstartU, *procstartV;
    int P, q, *procindexU, *procindexV, *NprocsUV;

    printf("Test InitIntersection: ");
    P = 17;   /* number of non-empty processors */
    b = 10*P; /* block size >= P */ 
    l= b*P;   /* P*P <= l */

    procstartU = (long *) malloc((l+1)* sizeof(long));
    procstartV = (long *) malloc((l+1)* sizeof(long));
    procindexU = (int *) malloc(2*l* sizeof(int));
    procindexV = (int *) malloc(l* sizeof(int));
    NprocsUV = (int *) malloc(l* sizeof(int));

    if (procstartU == NULL || procstartV == NULL ||
         procindexU == NULL || procindexV == NULL || NprocsUV == NULL) {
        printf("Error\n");
        exit(1);
    }

    /* Fill rows of communication matrix, each with two processors */
    for (i=0; i<l; i++)
        procstartU[i] = 2*i;
    procstartU[l] = 2*l;
    for (i=0; i<l; i++) {
        procindexU[2*i] = i%P;
        procindexU[2*i+1] = (i+1)%P;
    }
    /* Fill columns of communication matrix, each with one processor */
    for (i=0; i<l; i++)
        procstartV[i] = i;
    procstartV[l] = l;
    for (i=0; i<l; i++)
        procindexV[i] = i/b; 


    /* Initialise intersection */
    if (!InitIntersection(l, P, procstartU, procindexU, procstartV, procindexV, NprocsUV)) {
        printf("Error\n");
        exit(1);
    }

    /* Check result values */
    for (i=0; i<l; i++) { 
        q = i/b; /* the processor in column i */
        if (q == i%P) {
            /* intersection of row and column i contains 1 processor */
            if (NprocsUV[i] != 1 || procindexU[2*i] != q
                                 || procindexU[2*i+1] != (i+1)%P) {
                printf("Error\n");
                exit(1);
            }
        } else if (q == (i+1)%P) {
            if (NprocsUV[i] != 1 || procindexU[2*i] != q
                                 || procindexU[2*i+1] != i%P) {
                printf("Error\n");
                exit(1);   
            }
        } else {
            /* intersection is empty and has 0 processors */ 
            if (NprocsUV[i] != 0) {  
                printf("Error\n");
                exit(1);
            } 
        }
    }

    printf("OK\n");
    exit(0);

} /* end main */
