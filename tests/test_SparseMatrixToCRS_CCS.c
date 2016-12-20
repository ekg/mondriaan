#include "SparseMatrix.h"

void test_SparseMatrixToCRS_CCS();

int main(int argc, char **argv) {

    printf("Test SparseMatrixToCRS_CCS: ");
    
    SetRandomSeed(111);
    /*struct timeval tv;
    gettimeofday(&tv,NULL);
    unsigned long time_in_micros = 1000000 * tv.tv_sec + tv.tv_usec;
    SetRandomSeed(time_in_micros);*/
    
    test_SparseMatrixToCRS_CCS();

    printf("OK\n");
    exit(0);

} /* end main */


/**
 * Test SparseMatrixToCRS_CCS()
 */
void test_SparseMatrixToCRS_CCS() {
    
    struct sparsematrix A;
    
    long i, j, k, t, start, end;
    
    MMSparseMatrixInit(&A);
    A.MMTypeCode[0]='M'; /* normal matrix */
    A.MMTypeCode[1]='C'; /* coordinate scheme */
    A.MMTypeCode[2]='P'; /* pattern only */
    A.MMTypeCode[3]='G'; /* general, no symmetry */
    A.m = Random1(10, 20);
    A.n = Random1(10, 20);
    A.NrNzElts = Random1(A.m+A.n, (A.m+A.n)*3/2);
    A.i = (long *)malloc(A.NrNzElts*sizeof(long));
    A.j = (long *)malloc(A.NrNzElts*sizeof(long));
    
    if(A.i == NULL || A.j == NULL) {
        printf("Error\n");
        exit(1);
    }
    
    /* Add random nonzeros */
    for(t=0; t<A.NrNzElts; ++t) {
        A.i[t] = Random1(0, A.m-1);
        A.j[t] = Random1(0, A.n-1);
    }
    
    /* SparseMatrixToCRS_CCS() assumes there are no duplicates */
    if(!SparseMatrixRemoveDuplicates(&A)) {
        printf("Error\n");
        exit(1);
    }
    
    /* Convert to CCS/CRS */
    struct CRCS CCS, CRS;
    if(!SparseMatrixToCRS_CCS(&A, &CCS, &CRS)) {
        printf("Error\n");
        exit(1);
    }
    
    long *visitedCCS = (long *)calloc(A.NrNzElts, sizeof(long));
    long *visitedCRS = (long *)calloc(A.NrNzElts, sizeof(long));
    
    if(visitedCCS == NULL || visitedCRS == NULL) {
        printf("Error\n");
        exit(1);
    }
    
    /* Check that all nonzeros in A are present in CCS/CRS */
    int found;
    for(t=0; t<A.NrNzElts; ++t) {
        i = A.i[t];
        j = A.j[t];
        
        /* Check CCS */
        start = CCS.starts[j];
        end = CCS.starts[j+1];
        
        found = 0;
        for(k=start; k<end; ++k) {
            if(CCS.indices[k] == i && visitedCCS[k] == 0) {
                visitedCCS[k] = 1;
                found = 1;
                break;
            }
        }
        if(!found) {
            printf("Error\n");
            exit(1);
        }
        
        /* Check CRS */
        start = CRS.starts[i];
        end = CRS.starts[i+1];
        
        found = 0;
        for(k=start; k<end; ++k) {
            if(CRS.indices[k] == j && visitedCRS[k] == 0) {
                visitedCRS[k] = 1;
                found = 1;
                break;
            }
        }
        if(!found) {
            printf("Error\n");
            exit(1);
        }
    }
    
    /* Check that all nonzeros in CCS/CRS are found */
    for(k=0; k<A.NrNzElts; ++k) {
        if(visitedCCS[k] != 1 || visitedCRS[k] != 1) {
            printf("Error\n");
            exit(1);
        }
    }
    
    free(visitedCCS);
    free(visitedCRS);
    freeCRCS(&CCS);
    freeCRCS(&CRS);
    MMSparseMatrixFreeMemory(&A);
    
    
} /* end test_SparseMatrixToCRS_CCS */

