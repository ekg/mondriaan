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
    
    /* Compute solution */
    int P = 4;
    if (!SplitMatrixUpperBound(&A, P, &options)) {
        printf("Error: Unable to compute upper bound solution!\n");
        exit(0);
    }
    
    /* Check solution */
    if(!CheckUpperBoundSolution(&A)) {
        printf("Error: Invalid result.\n");
        exit(0);
    }
    
    /* Free memory */
    DestructDisconnectedMatrix(&A, FALSE, FALSE, FALSE, &submatrix_m, &submatrix_n, &submatrix_weights, &i_to_I, &j_to_J);
    
    
    printf("OK\n");
    exit(0);

} /* end main */
