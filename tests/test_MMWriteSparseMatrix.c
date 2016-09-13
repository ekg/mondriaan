#include "SparseMatrix.h"
#include "Options.h"

int main(int argc, char **argv) {

    long i, j, q, t, nzp;
    char filename[MAX_WORD_LENGTH];
    char line[MAX_LINE_LENGTH], line2[MAX_LINE_LENGTH];
    struct sparsematrix A;
    struct opts Options;
    FILE *fp;

    SetDefaultOptions(&Options);

    strcpy(filename,"outMMWriteSparseMatrix");
    A.MMTypeCode[0] = 'D';
    A.MMTypeCode[1] = 'C';
    A.MMTypeCode[2] = 'R';
    A.MMTypeCode[3] = 'G';
    A.ViewType      = ViewTypeOriginal;
    A.m = 5;
    A.n = 6;
    A.NrNzElts = 10;
    A.NrProcs = 2;
    nzp = A.NrNzElts / A.NrProcs; /* should be integer */
    strcpy(line, "% This is a 5 x 6 matrix\n% with 10 nonzeros\n");
    A.header  = line;
    A.comment = NULL;

    /* Allocate matrix arrays */
    A.i = (long *) malloc(A.NrNzElts* sizeof(long));
    A.j = (long *) malloc(A.NrNzElts* sizeof(long));
    A.ReValue = (double *) malloc(A.NrNzElts* sizeof(double));
    A.Pstart = (long *) malloc((A.NrProcs+1)* sizeof(long));

    if (A.i == NULL || A.j == NULL || A.ReValue == NULL || A.Pstart == NULL) {
        printf("Error\n");
        exit(1);
    }

    /* Matrix with nonzeros in first and last columns */
    t = 0; 
    for (i = 0; i< A.m; i++) {
        for (j = 0; j< A.n; j += A.n-1) {
            A.i[t] = i;
            A.j[t] = j;
            A.ReValue[t] = t;
            t++;
        }
    }

    for (q = 0; q <= A.NrProcs; q++) 
        A.Pstart[q] = nzp*q;

    strcpy(line2, "and with 2 processors\n");
    A.tail = line2;
    
    fp = fopen(filename, "w");
    MMWriteSparseMatrix(&A, fp, NULL, &Options);
    fclose(fp);
    exit(0);

} /* end main */
