#include "SparseMatrix.h"
#include "testHelper_DisconnectedMatrix.h"
#include <math.h>

void test_SparseMatrixComputeVolume(int symmetric);

int main(int argc, char **argv) {

    printf("Test SparseMatrixComputeVolume: ");
    
    SetRandomSeed(111);
    /*struct timeval tv;
    gettimeofday(&tv,NULL);
    unsigned long time_in_micros = 1000000 * tv.tv_sec + tv.tv_usec;
    SetRandomSeed(time_in_micros);*/
    
    test_SparseMatrixComputeVolume(FALSE);
    test_SparseMatrixComputeVolume(TRUE);
    
    printf("OK\n");
    exit(0);

} /* end main */


/**
 * Test SparseMatrixComputeVolume()
 * 
 * Input:
 * symmetric         : Whether the matrix should be symmetric
 */
void test_SparseMatrixComputeVolume(int symmetric) {
    
    struct sparsematrix A;
    struct opts options;
    
    options.SymmetricMatrix_UseSingleEntry = symmetric ? SingleEntYes : SingleEntNo;
    
    long i, j, p, t;
    
    long numSubmatrices, *submatrix_m, *submatrix_n, *submatrix_weights, **i_to_I, **j_to_J;
    ConstructDisconnectedMatrix(&A, symmetric, FALSE, FALSE, 1, 1, &numSubmatrices, &submatrix_m, &submatrix_n, &submatrix_weights, &i_to_I, &j_to_J);
    
    A.MMTypeCode[0]='D';
    A.NrProcs = Random1(2, 8);
    A.Pstart = (long*)realloc(A.Pstart, (A.NrProcs+1)*sizeof(long));
    
    /* Set distribution */
    A.Pstart[0] = 0;
    for(p=1; p<=A.NrProcs; ++p) {
        A.Pstart[p] = (p*A.NrNzElts)/A.NrProcs;
    }
    /* Compute total volume */
    long computedVolume = SparseMatrixComputeVolume(&A, &options);
    if(computedVolume == -1) {
        printf("Error\n");
        exit(1);
    }
    
    /* Check volume */
    if(symmetric) {
        SparseMatrixSymmetric2Full(&A);
    }
    
    long *out = (long*)calloc(A.NrProcs*A.m, sizeof(long));
    long *in = (long*)calloc(A.NrProcs*A.n, sizeof(long));
    
    for(p=0; p<A.NrProcs; ++p) {
        for(t=A.Pstart[p]; t<A.Pstart[p+1]; ++t) {
            ++out[A.i[t]*A.NrProcs + p];
            ++in[A.j[t]*A.NrProcs + p];
        }
    }
    
    long checkedVolume = 0, lambda;
    for(i=0; i<A.m; ++i) {
        lambda = 0;
        
        for(p=0; p<A.NrProcs; ++p) {
            if(out[i*A.NrProcs+p] > 0) {
                ++lambda;
            }
        }
        
        if(lambda > 1)
            checkedVolume += lambda-1;
    }
    
    for(j=0; j<A.n; ++j) {
        lambda = 0;
        
        for(p=0; p<A.NrProcs; ++p) {
            if(in[j*A.NrProcs+p] > 0) {
                ++lambda;
            }
        }
        
        if(lambda > 1)
            checkedVolume += lambda-1;
    }
    
    if(computedVolume != checkedVolume) {
        printf("Error\n");
        exit(1);
    }
    
    /* Free memory */
    free(in);
    free(out);
    DestructDisconnectedMatrix(&A, symmetric, FALSE, FALSE, &submatrix_m, &submatrix_n, &submatrix_weights, &i_to_I, &j_to_J);
    
    
} /* end test_SparseMatrixComputeVolume */

