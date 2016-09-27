#include "SparseMatrix.h"

int main(int argc, char **argv) {

    char filenameI[MAX_WORD_LENGTH], filenameA[MAX_WORD_LENGTH];
    struct sparsematrix A, I, M;
    long nz, s, t, p, P, *usedI, *usedA, matchI, matchA, found;
    FILE *fp;

    printf("Test MMReadSparseMatrixFromIndexAndValueMatrixFiles: ");
    
    /* Check sorting in general */
    
    strcpy(filenameA,"test_MMReadSparseMatrixFromIndexAndValueMatrixFilesA.inp");
    fp = fopen(filenameA, "r");
    MMReadSparseMatrix(fp, &A);
    fclose(fp);
    
    strcpy(filenameI,"test_MMReadSparseMatrixFromIndexAndValueMatrixFilesI.inp");
    fp = fopen(filenameI, "r");
    MMReadSparseMatrix(fp, &I);
    fclose(fp);

    nz = 12;
    P = 2;
    
    
    /* Use values as processor numbers */
    if (!MMReadSparseMatrixFromIndexAndValueMatrixFiles(filenameA, filenameI, &M)) {
        printf("Error!\n");
        exit(1);
    }
    
    
    usedI = (long *) malloc(nz * sizeof(long));
    usedA = (long *) malloc(nz * sizeof(long));
    for(s=0; s<nz; s++) {
        usedI[s] = 0;
        usedA[s] = 0;
    }
    
    matchI = matchA = 0;
    for(p=1; p<=P; p++) {
        for(t=M.Pstart[p-1]; t<M.Pstart[p]; t++) {
            
            found = 0;
        
            for(s=0; s<nz; s++) {
                if(usedI[s] == 0 && M.i[t] == I.i[s] && M.j[t] == I.j[s] && (long)I.ReValue[s]+1 == p) {
                    matchI++;
                    usedI[s] = 1;
                    found++;
                    break;
                }
            }
            
            for(s=0; s<nz; s++) {
                if(usedA[s] == 0 && M.i[t] == A.i[s] && M.j[t] == A.j[s] && M.ReValue[t] == A.ReValue[s]) {
                    matchA++;
                    usedA[s] = 1;
                    found++;
                    break;
                }
            }
            
            if(found != 2) {
                printf("Error\n");
                exit(1);
            }
            
        }
    }
    
    if(matchI != nz || matchA != nz) {
        printf("Error\n");
        exit(1);
    }
    
    MMDeleteSparseMatrix(&A);
    MMDeleteSparseMatrix(&I);
    MMDeleteSparseMatrix(&M);
    free(usedI);
    free(usedA);
    

    printf("OK\n");
    exit(0);
    
} /* end main */
