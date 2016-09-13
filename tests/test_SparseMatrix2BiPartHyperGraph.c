#include "Graph.h"

struct opts Options;

int main(int argc, char **argv) {

    /* This test function is geared towards the current SparseMatrix2BiPartHyperGraph
       function, and will not work if this function changes the order of the vertices,
       nets, or pins. */
       
    struct sparsematrix A;
    struct biparthypergraph HG;
    long n, i, j, t;
    
    SetDefaultOptions(&Options);

    printf("Test SparseMatrix2BiPartHyperGraph: ");
    n = 100; /* n by n symmetric tridiagonal matrix A, n even.
                 The matrix is represented by the nonzeros on or below the diagonal. */
    A.m = n;
    A.n = n;
    A.NrNzElts = 2*n-1; 

    A.i = (long *) malloc(A.NrNzElts* sizeof(long));
    A.j = (long *) malloc(A.NrNzElts* sizeof(long));
    A.RowLambda = (int *)malloc(A.m*sizeof(int));
    A.ColLambda = (int *)malloc(A.n*sizeof(int));
    A.RowMark = (int *)malloc(A.m*sizeof(int));
    A.ColMark = (int *)malloc(A.n*sizeof(int));
    A.ReValue = (double *) malloc(A.NrNzElts* sizeof(double));
    A.dummy= (int *) malloc(n* sizeof(int));

    if ( A.i == NULL || A.j == NULL || A.ReValue  == NULL || A.dummy == NULL ){
        printf("Error\n");
        exit(1);
    }

    /* Initialise matrix */
    A.MMTypeCode[2]='R'; /* real matrix */
    A.MMTypeCode[3]='S'; /* symmetric matrix */

    Options.SymmetricMatrix_UseSingleEntry = SingleEntYes;
    A.NrDummies = n/2;
    for (i=0; i<n/2; i++)
        A.dummy[i] = TRUE;
    for (i=n/2; i<n; i++)
        A.dummy[i] = FALSE;
 
    /* Fill matrix with nonzeros */
    t= 0;
    for (i=0; i<n; i++){
        for (j=MAX(0,i-1); j<=i; j++){
            A.i[t] = i;
            A.j[t] = j;
            if (i==j)
                A.ReValue[t] = -2.0;
            else
                A.ReValue[t] = 1.0;
            t++;
        }
    }    

    Options.SplitStrategy = OneDimRow; /* does not matter, but needs a value
                                           unequal to LocalRatio */

    if (!SparseMatrix2BiPartHyperGraph(&A, ROW, &Options, &HG)) {
        printf("Error\n");
        exit(1);
    }
    

    
    /* Check result values  */
    if (HG.NrVertices != n ||
        HG.NrNets != n ||       
        HG.NrPins != 2*n-1 ) {
        
        printf("Error\n");
        exit(1);
    }

    /* Check vertex values */
    for (t=0; t<n; t++){
        if ((t==0 && HG.V[t].vtxwgt != 0) ||
            (t>0 && t<n/2 && HG.V[t].vtxwgt != 2) ||
            (t>=n/2 && HG.V[t].vtxwgt != 3) ||
            (t==0 && HG.V[t].iStart != 0) ||
            (t>0 && HG.V[t].iStart != 2*t-1) ||
            HG.V[t].iEnd != 2*t+1 || 
            HG.Vtx2MatIndex[t] != t ) {
            
            printf("Error\n");
            exit(1);
        }
    }

    /* Check net values */
    for (t=0; t<n; t++){
        if (HG.N[t].netwgt != 1 ||   /* uniform netweights assumed */
            HG.N[t].iStartP0 != 2*t ||
            HG.N[t].iStartP1 != HG.N[t].iStartP0 ||
            (t<n-1 && HG.N[t].iEnd != 2*t+2) ||
            (t==n-1 && HG.N[t].iEnd != 2*n-1) || 
            HG.N[t].dir != COL || 
            HG.Net2MatIndex[t] != t ) {
 
            printf("Error\n");
            exit(1); 
        }
    }    
    
    /* Check pin values */ 
    for (t=0; t<2*n-1; t++){
        if (HG.VtxAdjncy[t] != t/2 ||               
            (t%2==0 && HG.NetAdjncy[t] != t/2 ) ||
            (t%2==1 && HG.NetAdjncy[t] != t/2 + 1) ||
            (t%2==0 && HG.MatReValue[t] != -2.0) ||
            (t%2==1 && HG.MatReValue[t] != 1.0) ) {
 
            printf("Error\n");
            exit(1);
        }
    }    

    printf("OK\n");
    exit(0);

} /* end main */
