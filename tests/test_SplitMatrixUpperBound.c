#include "SplitMatrixUpperBound.h"
#include "testHelper_DisconnectedMatrix.h"
#include "DistributeMat.h"
#include "DistributeVecLib.h"
#include <math.h>


int main(int argc, char **argv) {

    printf("Test SplitMatrixUpperBound: ");
    
    SetRandomSeed(111);
    /*struct timeval tv;
    gettimeofday(&tv,NULL);
    unsigned long time_in_micros = 1000000 * tv.tv_sec + tv.tv_usec;
    SetRandomSeed(time_in_micros);*/
    
    
    struct sparsematrix A;
    struct opts options;
    
    long numSubmatrices, *submatrix_m, *submatrix_n, *submatrix_weights, **i_to_I, **j_to_J;
    ConstructDisconnectedMatrix(&A, FALSE, FALSE, FALSE, 1, 1, &numSubmatrices, &submatrix_m, &submatrix_n, &submatrix_weights, &i_to_I, &j_to_J);
    
    /* Determine sums */
    long totalWeight = ComputeWeight(&A, 0, A.NrNzElts-1, NULL, &options);
    
    int P = 4;
    
    if (!SplitMatrixUpperBound(&A, P, &options)) {
        printf("Error: Unable to compute upper bound solution!\n");
        exit(0);
    }
    
    /* Calculate weights */
    long Wmax = 0, weight;
    for(int q=0; q<A.NrProcs; ++q) {
        weight = A.Pstart[q+1] - A.Pstart[q];
        if(weight > Wmax) {
            Wmax = weight;
        }
        totalWeight -= weight;
    }
    
    if(totalWeight != 0) {
        printf("Error: Invalid total weight.\n");
        exit(0);
    }
    
    if(Wmax > ceil(A.NrNzElts/(double)P)) {
        printf("Error: Invalid imbalance result.\n");
        exit(0);
    }
    
    /* Calculate communication volume */
    long ComVol1, ComVol2, tmp;
    CalcCom(&A, NULL, (A.m < A.n)?ROW:COL, &ComVol1, &tmp, &tmp, &tmp, &tmp);
    CalcCom(&A, NULL, (A.m < A.n)?COL:ROW, &ComVol2, &tmp, &tmp, &tmp, &tmp);
    
    long n = (A.m < A.n)?A.m:A.n;
    if(ComVol1 > n*(P-1) || ComVol2 > P-1) {
        printf("Error: Invalid communication result.\n");
        exit(0);
    }
    
    /* Free memory */
    DestructDisconnectedMatrix(&A, FALSE, FALSE, FALSE, &submatrix_m, &submatrix_n, &submatrix_weights, &i_to_I, &j_to_J);
    
    
    printf("OK\n");
    exit(0);

} /* end main */
