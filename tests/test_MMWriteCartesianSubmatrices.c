#include "Cartesian.h"

int main(int argc, char **argv) {

    char filename[MAX_WORD_LENGTH];
    struct sparsematrix A;
    long n, nz, i, j, t, q, P;
    FILE *fp;

    strcpy(filename,"outMMWriteCartesianSubmatrices");

    n = 40; /* n by n pattern matrix A, checkerboard, n even */
    A.m = n;
    A.n = n;
    A.NrNzElts = n*n/2; 
    nz = A.NrNzElts;
    P = 4;

    A.i = (long *) malloc(nz* sizeof(long));
    A.j = (long *) malloc(nz* sizeof(long));

    if (A.i == NULL || A.j == NULL) {
        printf("Error\n");
        exit(1);
    }

    /* Fill sparse matrix with nonzeros at locations (i,j) 
       with i, j both even */
    t = 0;
    for (i=0; i<n; i+=2) {
        for (j=0; j<n; j+=2) {
            A.i[t] = i;
            A.j[t] = j;
            t++;
        }
    }    

    /* Fill sparse matrix with nonzeros at locations (i,j)
       with i, j both odd */
    for (i=1; i<n; i+=2) {
        for (j=1; j<n; j+=2) {
            A.i[t] = i;
            A.j[t] = j;
            t++;
        }
    }

    /* Create distribution, where the last P-2 parts are empty */
    A.NrProcs = P;
    A.Pstart = (long *) malloc((P+1) * sizeof(long));
    A.Pstart[0] = 0;
    A.Pstart[1] = nz/2;
    for (q=2; q<=P; q++)
        A.Pstart[q] = nz;
    
    fp = fopen(filename, "w");
    MMWriteCartesianSubmatrices(&A, fp);
    fclose(fp);

    if (A.m != n || A.n != n   || A.NrNzElts != nz || A.NrProcs != P) {
        printf("Error\n");
        exit(1);
    }

    exit(0);

} /* end main */
