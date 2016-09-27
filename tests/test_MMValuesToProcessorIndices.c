#include "SparseMatrix.h"

int main(int argc, char **argv) {

    char filename[MAX_WORD_LENGTH];
    struct sparsematrix A_orig, A_dist;
    long nz, P, s, t, p, *used, match, found;
    FILE *fp;

    printf("Test MMValuesToProcessorIndices: ");

    strcpy(filename,"test_MMValuesToProcessorIndices.inp");

    nz = 13;
    P = 3;
    
    fp = fopen(filename, "r");
    MMReadSparseMatrix(fp, &A_orig);
    rewind(fp);
    MMReadSparseMatrix(fp, &A_dist);
    fclose(fp);
    
    /* Use values as processor numbers */
    if (!MMValuesToProcessorIndices(&A_dist)) {
        printf("Error!\n");
        exit(1);
    }
    
    used = (long *) malloc(nz * sizeof(long));
    for(s=0; s<nz; s++) {
        used[s] = 0;
    }
    
    match = 0;
    found = 0;
    for(t=0; t<nz; t++) {
        
        found = 0;
        
        for(p=1; p<=P; p++) {
            for(s=A_dist.Pstart[p-1]; s<A_dist.Pstart[p]; s++) {
                if(used[s] == 0 && A_orig.i[t] == A_dist.i[s] && A_orig.j[t] == A_dist.j[s] && (long)A_orig.ReValue[t]+1 == p) {
                    match++;
                    used[s] = 1;
                    found = 1;
                    break;
                }
            }
            if(found) {
                break;
            }
        }
    }
    
    if(match != nz) {
        printf("Error!\n");
        exit(1);
    }

    MMDeleteSparseMatrix(&A_orig);
    MMDeleteSparseMatrix(&A_dist);
    free(used);

    printf("OK\n");
    exit(0);
    
} /* end main */
