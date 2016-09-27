#include "SparseMatrix.h"

int main(int argc, char **argv) {

    char filename[MAX_WORD_LENGTH];
    struct sparsematrix A_orig, A_sort;
    long nz, s, t, p, P, *used, match;
    FILE *fp;

    printf("Test SpMatSortNonzeros: ");
    
    /* Check sorting in general */
    
    strcpy(filename,"test_SpMatSortNonzeros1.inp");

    nz = 12;
    
    fp = fopen(filename, "r");
    MMReadSparseMatrix(fp, &A_orig);
    rewind(fp);
    MMReadSparseMatrix(fp, &A_sort);
    fclose(fp);
    
    /* Use values as processor numbers */
    if (!SpMatSortNonzeros(&A_sort, 0)) {
        printf("Error!\n");
        exit(1);
    }
    
    used = (long *) malloc(nz * sizeof(long));
    for(s=0; s<nz; s++) {
        used[s] = 0;
    }
    
    /* Check completeness */
    match = 0;
    for(t=0; t<nz; t++) {
        
        for(s=0; s<nz; s++) {
            if(used[s] == 0 && A_orig.i[t] == A_sort.i[s] && A_orig.j[t] == A_sort.j[s] && A_orig.ReValue[t] == A_sort.ReValue[s]) {
                match++;
                used[s] = 1;
                break;
            }
        }
    }
    
    if(match != nz) {
        printf("Error!\n");
        exit(1);
    }
    
    /* Check order */
    for(s=1; s<nz; s++) {
        if(A_sort.i[s] < A_sort.i[s-1]) {
            printf("Error!\n");
            exit(1);
        }
        
        if(A_sort.i[s] == A_sort.i[s-1] && A_sort.j[s] < A_sort.j[s-1]) {
            printf("Error!\n");
            exit(1);
        }
    }
    
    MMDeleteSparseMatrix(&A_orig);
    MMDeleteSparseMatrix(&A_sort);
    free(used);
    
    /* Check sorting per part/processor */
    
    strcpy(filename,"test_SpMatSortNonzeros2.inp");

    nz = 13;
    P = 3;
    
    fp = fopen(filename, "r");
    MMReadSparseMatrix(fp, &A_orig);
    rewind(fp);
    MMReadSparseMatrix(fp, &A_sort);
    fclose(fp);
    
    /* Use values as processor numbers */
    if (!SpMatSortNonzeros(&A_sort, 1)) {
        printf("Error!\n");
        exit(1);
    }
    
    used = (long *) malloc(nz * sizeof(long));
    for(s=0; s<nz; s++) {
        used[s] = 0;
    }
    
    /* Check completeness */
    match = 0;
    for(p=1; p<=P; p++) {
        for(t=A_orig.Pstart[p-1]; t<A_orig.Pstart[p]; t++) {
        
            for(s=A_orig.Pstart[p-1]; s<A_orig.Pstart[p]; s++) {
                if(used[s] == 0 && A_orig.i[t] == A_sort.i[s] && A_orig.j[t] == A_sort.j[s] && A_orig.ReValue[t] == A_sort.ReValue[s]) {
                    match++;
                    used[s] = 1;
                    break;
                }
            }
        }
    }
    
    if(match != nz) {
        printf("Error!\n");
        exit(1);
    }
    
    /* Check order */
    for(p=1; p<=P; p++) {
        for(s=A_sort.Pstart[p-1]+1; s<A_sort.Pstart[p]; s++) {
            if(A_sort.i[s] < A_sort.i[s-1]) {
                printf("Error!\n");
                exit(1);
            }
            
            if(A_sort.i[s] == A_sort.i[s-1] && A_sort.j[s] < A_sort.j[s-1]) {
                printf("Error!\n");
                exit(1);
            }
        }
    }
    
    MMDeleteSparseMatrix(&A_orig);
    MMDeleteSparseMatrix(&A_sort);
    free(used);


    printf("OK\n");
    exit(0);
    
} /* end main */
