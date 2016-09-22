#ifndef MATRIX_H_DONE

#define MATRIX_H_DONE

#include "../include/Mondriaan.h"
#include <stdlib.h>
#include <stdio.h>

#define ROW 0
#define COL 1

struct mat {
    /* CRS+CCS structure holding input matrix information.
       This should be set during initialisation and not changed afterwards */

    /* number of matrix nonzeros */
    unsigned int nnz;
    
    /* maximum number of nonzeros in a row or column */
    unsigned int cmax;

    /* number of matrix rows and columns */
    unsigned int m, n;

    unsigned int* indices[2];
    unsigned int* starts[2];
    unsigned int* size[2];

};

char readmatrixfromfile(struct mat *in, const char *fn);
char convertSparsematrixToMat(struct mat *in, struct sparsematrix *A);

#endif
