#include "DistributeVecLib.h"

int main(int argc, char **argv) {

    struct sparsematrix A;
    long t, LB;
    int P, q, Pact;
    
    printf("Test CalcLocalLowerBound: ");
    
    /* Hardcoded example 4.10 (p. 200) from book 
       Parallel Scientific Computation by Rob H. Bisseling */
       
    P = 4; /* number of processors */
    A.m = 4;
    A.n = 8;
    A.NrNzElts = 20; 
    A.NrProcs = P;

    A.i = (long *) malloc(A.NrNzElts* sizeof(long));
    A.j = (long *) malloc(A.NrNzElts* sizeof(long));
    A.Pstart = (long *) malloc((P+1)* sizeof(long));

    if ( A.i == NULL || A.j  == NULL || A.Pstart == NULL ){
        printf("Error\n");
        exit(1);
    }

    A.Pstart[0] = 0;
    A.Pstart[1] = 6;
    A.Pstart[2] = 12;
    A.Pstart[3] = 16;
    A.Pstart[4] = 20;

    /* Fill matrix with nonzeros (row indices) */
    for (q=0; q<P; q++)
        for (t=A.Pstart[q]; t<A.Pstart[q+1]; t++)
                A.i[t] = q;
  
    /* Fill matrix with nonzeros (column indices) */
    A.j[0] = 0; A.j[1] = 2; A.j[2] = 4; A.j[3] = 5; A.j[4] = 6;  A.j[5] = 7; 
    A.j[6] = 0; A.j[7] = 1; A.j[8] = 3; A.j[9] = 4; A.j[10] = 5; A.j[11] = 6; 
    A.j[12] = 1; A.j[13] = 5; A.j[14] = 6; A.j[15] = 7;
    A.j[16] = 2; A.j[17] = 3; A.j[18] = 4; A.j[19] = 7;

    if (!CalcLocalLowerBound(&A, ROW, &LB, &Pact)) {
        printf("Error\n");
        exit(1);
    }

    /* Check result values */
    if (LB != 4 || Pact  != P ){
        printf("Error\n");
        exit(1);
    }

    printf("OK\n");
    exit(0);

} /* end main */
