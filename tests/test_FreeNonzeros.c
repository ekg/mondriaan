#include "FreeNonzeros.h"
#include "DistributeMat.h"
#include "DistributeVecLib.h"
#include <math.h>

void test_FreeNonzeros(int symmetric);

int main(int argc, char **argv) {

    printf("Test FreeNonzeros: ");
    
    SetRandomSeed(111);
    /*struct timeval tv;
    gettimeofday(&tv,NULL);
    unsigned long time_in_micros = 1000000 * tv.tv_sec + tv.tv_usec;
    SetRandomSeed(time_in_micros);*/
    
    test_FreeNonzeros(FALSE);
    test_FreeNonzeros(TRUE);
    
    printf("OK\n");
    exit(0);

} /* end main */

/**
 * Compute the total communication volume of a matrix.
 */
long ComputeVolume(struct sparsematrix *pM, int symmetric) {
    if(symmetric) {
        SparseMatrixSymmetric2Full(pM);
    }
    
    long tmp, ComVolV, ComVolU;
    CalcCom(pM, NULL, ROW, &ComVolV, &tmp, &tmp, &tmp, &tmp);
    CalcCom(pM, NULL, COL, &ComVolU, &tmp, &tmp, &tmp, &tmp);
    
    if(symmetric) {
        SparseMatrixFull2Symmetric(pM, 'S');
    }
    
    return ComVolV+ComVolU;

} /* end ComputeVolume */


/**
 * Test ImproveFreeNonzeros*()
 * 
 * Input:
 * symmetric: Whether the matrix should be symmetric
 */
void test_FreeNonzeros(int symmetric) {
    
    struct sparsematrix A;
    struct sparsematrix *pA = &A;
    struct opts options;
    
    options.SymmetricMatrix_UseSingleEntry = symmetric ? SingleEntYes : SingleEntNo;
    
    long t, p, i, j;
    long P = 5, m = 20, n = 20;
    
    /* Set up matrix struct */
    MMSparseMatrixInit(pA);
    pA->m = m;
    pA->n = n;
    pA->NrNzElts = P*(m+n)+((m-P)*(n-P))/5;
    pA->NrProcs = P; /* maximum number of parts */

    pA->i = (long *)malloc(pA->NrNzElts*sizeof(long));
    pA->j = (long *)malloc(pA->NrNzElts*sizeof(long));
    pA->Pstart = (long *)malloc((P+1)*sizeof(long));
    
    if (pA->i == NULL || pA->j == NULL || pA->Pstart == NULL) {
        printf("Error\n");
        exit(1);
    }
    
    pA->MMTypeCode[0]='D'; /* normal matrix */
    pA->MMTypeCode[1]='C'; /* coordinate scheme */
    pA->MMTypeCode[2]='P'; /* pattern only */
    if(symmetric)
        pA->MMTypeCode[3]='S'; /* symmetric */
    else
        pA->MMTypeCode[3]='G'; /* general, no symmetry */
    
    /* Fill matrix */
    t = 0;
    for(p=0; p<P; ++p) {
        /* First, fill rows and columns such that we have many processors per row/column */
        pA->Pstart[p] = t;
        for(i=P; i<pA->m; ++i) {
            if(Random1(0, 2) == 0) {
                continue;
            }
            pA->i[t] = i;
            pA->j[t] = p;
            ++t;
        }
        if(!symmetric) {
            for(j=P; j<pA->n; ++j) {
                if(Random1(0, 2) == 0) {
                    continue;
                }
                pA->i[t] = p;
                pA->j[t] = j;
                ++t;
            }
        }
    }
    
    /* Add additional nonzeros to last partition */
    i = j = 0;
    while(TRUE) {
        j += Random1(5, 10);
        if(symmetric) {
            while(j > i) {
                j -= i;
                ++i;
            }
        }
        else {
            while(j >= pA->n-P) {
                j -= pA->n-P;
                ++i;
            }
        }
        if(i >= pA->m-P) {
            break;
        }
        pA->i[t] = P+i;
        pA->j[t] = P+j;
        ++t;
    }
    pA->Pstart[P] = pA->NrNzElts = t;
    
    /* We should provide a procs array */
    int *procs = (int *)malloc(P*sizeof(int));
    if (procs == NULL) {
        printf("Error\n");
        exit(1);
    }
    for(p=0; p<P; ++p) {
        procs[p] = 1;
    }
    
    /* Compute statistics before applying algorithm */
    long totalImbalanceBefore = 0, weight;
    long avgNrNzElts = pA->NrNzElts / P;
    for(p=0; p<P; ++p) {
        weight = ComputeWeight(pA, pA->Pstart[p], pA->Pstart[p+1]-1, NULL, &options);
        totalImbalanceBefore += labs(weight-avgNrNzElts);
    }
    long volumeBefore = ComputeVolume(pA, symmetric);
    
    /* Run algorithm */
    ImproveFreeNonzeros(pA, &options, procs, 3, 4);
    
    /* Compute statistics after applying algorithm */
    long totalImbalanceAfter = 0.0;
    for(p=0; p<P; ++p) {
        weight = ComputeWeight(pA, pA->Pstart[p], pA->Pstart[p+1]-1, NULL, &options);
        totalImbalanceAfter += labs(weight-avgNrNzElts);
    }
    long volumeAfter = ComputeVolume(pA, symmetric);
    
    /* Check that the total imbalance has not increased */
    if(totalImbalanceAfter > totalImbalanceBefore) {
        printf("Error1\n");
        exit(1);
    }
    
    /* Check that the volume has not increased */
    if(volumeAfter > volumeBefore) {
        printf("Error2\n");
        exit(1);
    }
    
    free(procs);
    MMDeleteSparseMatrix(pA);
    
} /* end test_FreeNonzeros */

