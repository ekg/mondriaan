#include "ZeroVolumeSearch.h"
#include "testHelper_DisconnectedMatrix.h"
#include "DistributeMat.h"
#include "DistributeVecLib.h"
#include <math.h>

void test_ZeroVolumeSearch(int symmetric, int dummies, int colWeights, int large, double eps);

int main(int argc, char **argv) {

    printf("Test ZeroVolumeSearch: ");
    
    SetRandomSeed(111);
    /*struct timeval tv;
    gettimeofday(&tv,NULL);
    unsigned long time_in_micros = 1000000 * tv.tv_sec + tv.tv_usec;
    SetRandomSeed(time_in_micros);*/
    
    int i;
    double eps[3] = {0.0, 0.03, 0.3};
    for(i=0; i<3; i++) {
        test_ZeroVolumeSearch(FALSE, FALSE, FALSE, FALSE, eps[i]);
        test_ZeroVolumeSearch(TRUE,  FALSE, FALSE, FALSE, eps[i]);
        test_ZeroVolumeSearch(FALSE, TRUE,  FALSE, FALSE, eps[i]);
        test_ZeroVolumeSearch(TRUE,  TRUE,  FALSE, FALSE, eps[i]);
        
        test_ZeroVolumeSearch(FALSE, FALSE, TRUE, FALSE, eps[i]);
        test_ZeroVolumeSearch(TRUE,  FALSE, TRUE, FALSE, eps[i]);
        test_ZeroVolumeSearch(FALSE, TRUE,  TRUE, FALSE, eps[i]);
        test_ZeroVolumeSearch(TRUE,  TRUE,  TRUE, FALSE, eps[i]);
        
        test_ZeroVolumeSearch(FALSE, FALSE, FALSE, TRUE, eps[i]);
        test_ZeroVolumeSearch(TRUE,  FALSE, FALSE, TRUE, eps[i]);
        test_ZeroVolumeSearch(FALSE, TRUE,  FALSE, TRUE, eps[i]);
        test_ZeroVolumeSearch(TRUE,  TRUE,  FALSE, TRUE, eps[i]);
        
        test_ZeroVolumeSearch(FALSE, FALSE, TRUE, TRUE, eps[i]);
        test_ZeroVolumeSearch(TRUE,  FALSE, TRUE, TRUE, eps[i]);
        test_ZeroVolumeSearch(FALSE, TRUE,  TRUE, TRUE, eps[i]);
        test_ZeroVolumeSearch(TRUE,  TRUE,  TRUE, TRUE, eps[i]);
    }
    
    printf("OK\n");
    exit(0);

} /* end main */


/**
 * Test DetectConnectedComponents()
 * 
 * Input:
 * symmetric         : Whether the matrix should be symmetric
 * dummies           : Whether the matrix should contain dummy diagonal nonzeros
 * colWeights        : Whether the matrix should be weighted with column weights
 * large             : Whether the matrix should be large (small), such that Karmarkar-Karp
 *                     (naive enumeration) is used for the subset-sum problem.
 * eps               : The maximum allowed load imbalance
 */
void test_ZeroVolumeSearch(int symmetric, int dummies, int colWeights, int large, double eps) {
    
    struct sparsematrix A;
    struct opts options;
    
    options.SymmetricMatrix_UseSingleEntry = symmetric ? SingleEntYes : SingleEntNo;
    
    /* In this function, we distinguish between the components we generate (calling them submatrices) and
     * the components we detect (calling them just components), to make clean what we are dealing with.
     */
    long submat, mid;
    
    long numSubmatrices, *submatrix_m, *submatrix_n, *submatrix_weights, **i_to_I, **j_to_J;
    ConstructDisconnectedMatrix(&A, symmetric, dummies, colWeights, large?20:2, large?25:5, &numSubmatrices, &submatrix_m, &submatrix_n, &submatrix_weights, &i_to_I, &j_to_J);
    
    /* Determine sums */
    long searchWeight = 0, totalWeight = 0;
    for(submat=0; submat<numSubmatrices; submat++) {
        if(submat%2 == 0) {
            searchWeight += submatrix_weights[submat];
        }
        totalWeight += submatrix_weights[submat];
    }
    
    /* Check the weights */
    if(totalWeight != ComputeWeight(&A, 0, A.NrNzElts-1, NULL, &options)) {
        printf("Error\n");
        exit(1);
    }
    
    /* Determine which part has minimum weight */
    if(searchWeight > totalWeight-searchWeight) {
        searchWeight = totalWeight-searchWeight;
    }
    
    /* Augment the search weight within what is allowed by eps */
    searchWeight = (long)ceil(((double)searchWeight)/(1+0.5*eps));
    
    /* Determine weight of small and large parts */
    long weightlo = searchWeight;
    long weighthi = totalWeight-searchWeight;
    
    /* Increase the allowed load imbalance */
    weightlo *= (1+eps);
    weighthi *= (1+eps);
    
    /* Run zero volume search */
    int foundZeroVolumePartition = ZeroVolumeSearch(&A, weightlo, weighthi, &mid, &options);
    
    /* By construction, a zero volume partition exists */
    if(foundZeroVolumePartition) {
        
        A.Pstart[1] = mid;
        
        /* Check the weights of the found partitionings */
        long weight0 = ComputeWeight(&A, 0, mid-1, NULL, &options);
        long weight1 = ComputeWeight(&A, mid, A.NrNzElts-1, NULL, &options);
        if(weight0 > weightlo || weight1 > weighthi) {
            printf("Error\n");
            exit(1);
        }
        
        /* Calculate communication volume. Change to full matrix for CalcCom() */
        if(symmetric) {
            SparseMatrixSymmetric2Full(&A);
        }
        long tmp, ComVolV, ComVolU;
        if (!CalcCom(&A, NULL, ROW, &ComVolV, &tmp, &tmp, &tmp, &tmp) ||
            !CalcCom(&A, NULL, COL, &ComVolU, &tmp, &tmp, &tmp, &tmp)) {
            printf("Error\n");
            exit(1);
        }
        
        if(ComVolV != 0 || ComVolU != 0) {
            printf("Error\n");
            exit(1);
        }
        
    }
    else if(!large) {
        printf("Error\n");
        exit(1);
        /* When using a large amount of components, the algorithm uses karmarkar-karp and
         * we cannot rely on the fact that a zero volume partitioning exists, because
         * the K-K heuristic might not find it. Hence we do not error out in that case.
         */
    }
    
    /* Free memory */
    DestructDisconnectedMatrix(&A, symmetric, dummies, colWeights, &submatrix_m, &submatrix_n, &submatrix_weights, &i_to_I, &j_to_J);
    
} /* end test_ZeroVolumeSearch */

