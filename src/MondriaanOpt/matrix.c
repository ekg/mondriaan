#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include "../include/Mondriaan.h"
#include "matrix.h"

int compare (const void *a, const void *b) {

    /* This functions compares the values a and b, and returns a-b.
       It is needed for quicksort. */

    return ( *(unsigned int*)a - *(unsigned int*)b );

} /* end compare */


char readmatrixfromfile(struct mat *in, const char *fn){

    /* This function reads a matrix from file 'fn' using Mondriaan I/O functions,
       and puts information in the matrix structure 'in', in CRS+CCS format
       (i.e., compressed row storage + compressed column storage).
       The function returns FALSE if there is an error, and TRUE otherwise */

    /* File pointer */
    FILE *fp;

    /* Sparse matrix in Mondriaan format */
    struct sparsematrix A;

    /* Open file and read matrix */
    fp = fopen(fn,"r");
    if(fp==NULL) {
        fprintf(stderr, "Could not read matrix!\n");
        exit(EXIT_FAILURE);
    }

    if (!MMReadSparseMatrix(fp, &A)) {
        fprintf(stderr, "Could not read matrix!\n");
        exit(EXIT_FAILURE);
    }
    fclose(fp);
    
    return convertSparsematrixToMat(in, &A);
} /* end readmatrixfromfile */

char convertSparsematrixToMat(struct mat *in, struct sparsematrix *A) {

    int i;

    /* Position of a matrix nonzero, 0 <= x < m and 0 <= y < n,
       where A is an m by n matrix */
    unsigned int x, y;
    
    /* Remove duplicate nonzeros */
    if (!SparseMatrixRemoveDuplicates(A)){
        fprintf(stderr, "Could not remove duplicates!\n");
        exit(EXIT_FAILURE);
    } 

    /* Check symmetry and expand */
    if (A->m == A->n && (A->MMTypeCode[3] == 'S' || A->MMTypeCode[3] == 'K' || A->MMTypeCode[3] == 'H')){
        SparseMatrixSymmetric2Full(A);
    }

    in->indices[ROW] = malloc(sizeof(unsigned int)*A->NrNzElts);
    in->starts[ROW] = malloc(sizeof(unsigned int)*(A->m+1));
    in->size[ROW] = malloc(sizeof(unsigned int)*A->m);
    in->indices[COL] = malloc(sizeof(unsigned int)*A->NrNzElts);
    in->starts[COL] = malloc(sizeof(unsigned int)*(A->n+1));
    in->size[COL] = malloc(sizeof(unsigned int)*A->n);
    if (in->indices[ROW] == NULL || in->starts[ROW] == NULL || in->size[ROW] == NULL ||
        in->indices[COL] == NULL || in->starts[COL] == NULL || in->size[COL]  == NULL){
        fprintf(stderr, "Not enough memory!\n");
        exit(EXIT_FAILURE);
    }

    /* Fill starts array */
    for(i=0;i<=A->m;i++)
        in->starts[ROW][i] = 0;
    for(i=0;i<=A->n;i++)
        in->starts[COL][i] = 0;

    for(i=0;i<A->NrNzElts;i++){
        x = A->i[i];
        y = A->j[i];
        in->starts[ROW][x+1]++;
        in->starts[COL][y+1]++;
    }

    for(i=2;i<=A->m;i++)
        in->starts[ROW][i] += in->starts[ROW][i-1];
    for(i=2;i<=A->n;i++)
        in->starts[COL][i] += in->starts[COL][i-1];

    /* Fill index arrays */
    unsigned int *nInRow = malloc(sizeof(unsigned int)*A->m);
    unsigned int *nInCol = malloc(sizeof(unsigned int)*A->n);
    if (nInRow == NULL || nInCol == NULL){
        fprintf(stderr, "Not enough memory!\n");
        exit(EXIT_FAILURE);
    }

    for(i=0;i<A->m;i++)
        nInRow[i] = 0;
    for(i=0;i<A->n;i++)
        nInCol[i] = 0;


    for(i=0;i<A->NrNzElts;i++){
        x = A->i[i];
        y = A->j[i];
        /* Fill first free position in row x with column index y */
        in->indices[ROW][in->starts[ROW][x]+nInRow[x]]=y;
        /* Fill first free position in column y with row index x */
        in->indices[COL][in->starts[COL][y]+nInCol[y]]=x;
        nInRow[x]++;
        nInCol[y]++;
    }

    /* Sort the indices within the rows and columns and
       determine cmax, the maximum number of nonzeros in a row or column */
    in->cmax=0;
    for(i=0;i<A->m;i++){
        in->size[ROW][i] = in->starts[ROW][i+1]-in->starts[ROW][i];
        if(in->size[ROW][i]>in->cmax)
            in->cmax = in->size[ROW][i];
        qsort(in->indices[ROW]+in->starts[ROW][i], in->size[ROW][i], sizeof(unsigned int), compare);
    }
    for(i=0;i<A->n;i++){
        in->size[COL][i] = in->starts[COL][i+1]-in->starts[COL][i];
        if(in->size[COL][i]>in->cmax)
            in->cmax = in->size[COL][i];
        qsort(in->indices[COL]+in->starts[COL][i], in->size[COL][i], sizeof(unsigned int), compare);
    }

    in->nnz = A->NrNzElts;
    in->m = A->m;
    in->n = A->n;

    free(nInRow);
    free(nInCol);

    return TRUE;

} /* end convertSparsematrixToMat */
