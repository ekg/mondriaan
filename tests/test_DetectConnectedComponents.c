#include "ZeroVolumeSearch.h"
#include "testHelper_DisconnectedMatrix.h"

void test_detectConnectedComponents(int symmetric, int dummies, int colWeights);

int main(int argc, char **argv) {

    printf("Test DetectConnectedComponents: ");
    
    SetRandomSeed(111);
    /*struct timeval tv;
    gettimeofday(&tv,NULL);
    unsigned long time_in_micros = 1000000 * tv.tv_sec + tv.tv_usec;
    SetRandomSeed(time_in_micros);*/
    
    test_detectConnectedComponents(FALSE, FALSE, FALSE);
    test_detectConnectedComponents(TRUE,  FALSE, FALSE);
    test_detectConnectedComponents(FALSE, TRUE,  FALSE);
    test_detectConnectedComponents(TRUE,  TRUE,  FALSE);
    
    test_detectConnectedComponents(FALSE, FALSE, TRUE);
    test_detectConnectedComponents(TRUE,  FALSE, TRUE);
    test_detectConnectedComponents(FALSE, TRUE,  TRUE);
    test_detectConnectedComponents(TRUE,  TRUE,  TRUE);

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
 */
void test_detectConnectedComponents(int symmetric, int dummies, int colWeights) {
    
    struct sparsematrix A;
    struct opts Options;
    
    SetDefaultOptions(&Options);
    Options.SymmetricMatrix_UseSingleEntry = symmetric ? SingleEntYes : SingleEntNo;
    
    /* In this function, we distinguish between the components we generate (calling them submatrices) and
     * the components we detect (calling them just components), to make clean what we are dealing with.
     */
    long i, I, cmpnnt, submat;
    
    long numSubmatrices, *submatrix_m, *submatrix_n, *submatrix_weights, **i_to_I, **j_to_J;
    ConstructDisconnectedMatrix(&A, symmetric, dummies, colWeights, 2, 5, &numSubmatrices, &submatrix_m, &submatrix_n, &submatrix_weights, &i_to_I, &j_to_J);
    
    long maxWeight = 0;
    for(submat=0; submat<numSubmatrices; ++submat) {
        if(submatrix_weights[submat] > maxWeight) {
            maxWeight = submatrix_weights[submat];
        }
    }
    
    long numComponentsFound = 0;
    long *componentWeightsFound = NULL;
    long *rowAssignmentsFound = NULL;
    
    /* Search for connected components (this should fail) */
    if(DetectConnectedComponents(&A, maxWeight-1, &numComponentsFound, &componentWeightsFound, &rowAssignmentsFound, &Options)) {
        printf("Error\n");
        exit(1);
    }
    
    /* Search for connected components (this should succeed) */
    if(!DetectConnectedComponents(&A, maxWeight, &numComponentsFound, &componentWeightsFound, &rowAssignmentsFound, &Options)) {
        printf("Error\n");
        exit(1);
    }
    
    /* Check number of components */
    if(numComponentsFound != numSubmatrices) {
        printf("Error\n");
        exit(1);
    }
    
    /*for(cmpnnt=0; cmpnnt<numComponentsFound; ++cmpnnt) {
        fprintf(stderr, "%ld ", submatrix_weights[cmpnnt]);
    }fprintf(stderr, "\n");
    for(cmpnnt=0; cmpnnt<numComponentsFound; ++cmpnnt) {
        fprintf(stderr, "%ld ", componentWeightsFound[cmpnnt]);
    }fprintf(stderr, "\n");*/
    
    /* Check amount of nonzeros in each component */
    long *component2submatrix = (long *)calloc(numSubmatrices, sizeof(long));
    if (component2submatrix == NULL) {
        printf("Error\n");
        exit(1);
    }
    for(submat=0; submat<numSubmatrices; ++submat) {
        int found = 0;
        for(cmpnnt=0; cmpnnt<numComponentsFound; ++cmpnnt) {
            if(component2submatrix[cmpnnt] == 0 && componentWeightsFound[cmpnnt] == submatrix_weights[submat]) {
                component2submatrix[cmpnnt] = submat;
                found = 1;
                break;
            }
        }
        if(!found) {
            printf("Error\n");
            exit(1);
        }
    }
    
    /* Check all row assigments (check rowAssignmentsFound against i_to_I) */
    long *component_m_found = (long*)calloc(numSubmatrices, sizeof(long));
    if (component_m_found == NULL) {
        printf("Error\n");
        exit(1);
    }
    long submat_nnz;
    for(I=0; I<A.m; ++I) {
        submat = component2submatrix[rowAssignmentsFound[I]];
        submat_nnz = submatrix_weights[submat];
        
        /* Search whether the current row is actually part of the submatrix */
        int found = 0;
        for(i=0; i<submat_nnz; ++i) {
            if(i_to_I[submat][i] == I) {
                found = 1;
                break;
            }
        }
        if(!found) {
            printf("Error\n");
            exit(1);
        }
        
        ++component_m_found[submat];
    }
    
    /* Check the obtained row counts */
    for(submat=0; submat<numSubmatrices; ++submat) {
        if(component_m_found[submat] != submatrix_m[submat]) {
            printf("Error\n");
            exit(1);
        }
    }
    
    free(component2submatrix);
    free(component_m_found);
    free(componentWeightsFound);
    free(rowAssignmentsFound);
    
    DestructDisconnectedMatrix(&A, symmetric, dummies, colWeights, &submatrix_m, &submatrix_n, &submatrix_weights, &i_to_I, &j_to_J);
    
} /* end test_detectConnectedComponents */

