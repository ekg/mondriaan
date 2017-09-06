#include "SparseMatrix.h"
#include <math.h>

void test_SparseMatrixFull2Symmetric();

int main(int argc, char **argv) {

    printf("Test SparseMatrixFull2Symmetric: ");
    
    SetRandomSeed(111);
    /*struct timeval tv;
    gettimeofday(&tv,NULL);
    unsigned long time_in_micros = 1000000 * tv.tv_sec + tv.tv_usec;
    SetRandomSeed(time_in_micros);*/
    
    test_SparseMatrixFull2Symmetric();
    
    printf("OK\n");
    exit(0);

} /* end main */


/**
 * Test SparseMatrixFull2Symmetric()
 */
void test_SparseMatrixFull2Symmetric() {
    
    struct sparsematrix A, B;
    struct sparsematrix *pA = &A, *pB = &B;
    
    long t, p, i, j;
    long P = 5, n = 20;
    
    /* Set up matrix struct */
    MMSparseMatrixInit(pA);
    MMSparseMatrixInit(pB);
    pA->m = pB->m = n;
    pA->n = pB->n = n;
    pA->NrNzElts = pB->NrNzElts = (n*n+n)/2;
    pA->NrProcs = pB->NrProcs = P; /* maximum number of parts */

    pA->i = (long *)malloc(pA->NrNzElts*sizeof(long));
    pA->j = (long *)malloc(pA->NrNzElts*sizeof(long));
    pA->Pstart = (long *)malloc((P+1)*sizeof(long));
    
    pB->i = (long *)malloc(pA->NrNzElts*sizeof(long));
    pB->j = (long *)malloc(pA->NrNzElts*sizeof(long));
    pB->Pstart = (long *)malloc((P+1)*sizeof(long));
    
    if (pA->i == NULL || pA->j == NULL || pA->Pstart == NULL || pB->i == NULL || pB->j == NULL || pB->Pstart == NULL) {
        printf("Error\n");
        exit(1);
    }
    
    pA->MMTypeCode[0]=pB->MMTypeCode[0]='D'; /* normal matrix */
    pA->MMTypeCode[1]=pB->MMTypeCode[1]='C'; /* coordinate scheme */
    pA->MMTypeCode[2]=pB->MMTypeCode[2]='P'; /* pattern only */
    pA->MMTypeCode[3]=pB->MMTypeCode[3]='S'; /* symmetric */
    
    /* Fill matrix */
    t = 0;
    i = j = 0;
    while(TRUE) {
        j += Random1(1, 5);
        while(j > i) {
            j -= i;
            ++i;
        }
        if(i >= pA->m) {
            break;
        }
        pA->i[t] = pB->i[t] = i;
        pA->j[t] = pB->j[t] = j;
        ++t;
    }
    pA->NrNzElts = pB->NrNzElts = t;
    for(p=0; p<P; ++p) {
        pA->Pstart[p] = pB->Pstart[p] = (pA->NrNzElts/P)*p;
    }
    pA->Pstart[P] = pB->Pstart[P] = pA->NrNzElts;
    
    /* Run procedures */
    SparseMatrixSymmetric2Full(pA);
    SparseMatrixFull2Symmetric(pA, 'S');
    
    /* Check result */
    if(pA->NrNzElts != pB->NrNzElts || pA->NrProcs != P) {
        printf("Error\n");
        exit(1);
    }
    
    for(p=0; p<=P; ++p) {
        if(pA->Pstart[p] != pB->Pstart[p]) {
            printf("Error\n");
            exit(1);
        }
    }
    
    for(t=0; t<pA->NrNzElts; ++t) {
        if(pA->i[t] != pB->i[t] || pA->j[t] != pB->j[t]) {
            printf("Error\n");
            exit(1);
        }
    }
    
    MMDeleteSparseMatrix(pA);
    MMDeleteSparseMatrix(pB);
    
} /* end test_SparseMatrixFull2Symmetric */

