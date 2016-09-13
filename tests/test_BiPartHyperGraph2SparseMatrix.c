#include "Graph.h"
#include "Matalloc.h"

struct opts Options;

int main(int argc, char **argv) {

    struct sparsematrix A;
    struct biparthypergraph HG;
    long n, nz, i, j, t, mid, val, **a;
    
    SetDefaultOptions(&Options);

    printf("Test BiPartHyperGraph2SparseMatrix: ");
    n = 4; /* n by n real tridiagonal matrix A, n even */
    A.m = n;
    A.n = n;
    A.NrNzElts = 3*n-2; 
    nz = A.NrNzElts; /* even */

    A.i = (long *) malloc(nz* sizeof(long));
    A.j = (long *) malloc(nz* sizeof(long));
    A.RowLambda = (int *)malloc(A.m*sizeof(int));
    A.ColLambda = (int *)malloc(A.n*sizeof(int));
    A.RowMark = (int *)malloc(A.m*sizeof(int));
    A.ColMark = (int *)malloc(A.n*sizeof(int));
    A.ReValue = (double *) malloc(nz* sizeof(double));
    A.dummy= (int *) malloc(n* sizeof(int));
    a = (long **) matallocl(n,n);
 
    if ( A.i == NULL || A.j == NULL || A.ReValue  == NULL || A.dummy == NULL ||
         a == NULL ){
        printf("Error\n");
        exit(1);
    }

    /* Initialise matrix */
    A.MMTypeCode[2]='R'; /* real matrix */
    A.MMTypeCode[3]='G'; /* general matrix */
    A.NrDummies = 0;
 
    /* Initialise dense matrix */
    for (i=0; i<n; i++)
        for (j=0; j<n; j++)
            a[i][j] = 0;
            
    /* Fill dense and sparse matrices with nonzeros */
    t= 0;
    for (i=0; i<n; i++){
        for (j = MAX(0,i-1); j<=MIN(n-1, i+1); j++) {
            A.i[t] = i;
            A.j[t] = j;
            a[i][j] = i*n+j;
            A.ReValue[t] = a[i][j];
            t++;
        }
    }    

    Options.SplitStrategy = OneDimRow; /* does not matter, but needs a value
                                           unequal to LocalRatio */

    if (!SparseMatrix2BiPartHyperGraph(&A, FINEGRAIN, &Options, &HG)) {
        printf("Error\n");
        exit(1);
    }
    
    /* Assign half the vertices (nonzeros) to part 1 */ 
    for ( t = 0; t < HG.NrVertices; t++ ) 
        HG.OptPartVtx[t] = t%2;

    if (!BiPartHyperGraph2SparseMatrix(&HG, 0, nz-1, &mid, &A)) {
        printf("Error\n");
        exit(1);
    }
  
    /* Check result values  */
    if ( A.m != n || A.n != n  || A.NrNzElts != nz || 
         A.MMTypeCode[2] !='R' || A.MMTypeCode[3] !='G' ||
         A.NrDummies != 0 || mid != nz/2) {
        
        printf("Error\n");
        exit(1);
    }

    /* Check nonzero values */
    for (t=0; t<nz; t++) {
        i = A.i[t];
        j = A.j[t];
        val = A.ReValue[t] + 0.5 ; /* round to nearest integer value */
        if (val != i*n+j){
            printf("Error\n");
            exit(1);
        }  
        a[i][j] -= val; /* reset dense matrix element */
    }

    /* Check whether whole dense matrix is zero again */
    for (i=0; i<n; i++){
        for (j=0; j<n; j++){
            if (a[i][j] != 0) {
                printf("Error\n");
                exit(1);
            }
        }    
    }

    printf("OK\n");
    exit(0);

} /* end main */
