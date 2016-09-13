#include "DistributeVec.h"

struct opts Options;

int main(int argc, char **argv) {

    struct sparsematrix A;
    long P, k, n, i, j, t, nzp, maxcom,
         ComVol, MaxOut, MaxIn, MaxCompnts, TotCompnts;
    int proc;
    long int *X;

    printf("Test DistributeVec: ");
    k = 14; /* k even */
    P = 3*k; /* number of  processors >= 3 */

    /* P by P matrix A, containing three k by k diagonal blocks.
       The first block is empty; the second block is defined by b[i,j] = i,
       for 0 <= i, j < k; and the third block by b[i,j] = k + 2*i + j/(k/2) */

    A.m = 3*k;
    A.n = 3*k;
    n = A.m;
    A.NrNzElts = 2*k*k; 
    A.NrProcs = P;

    A.i = (long *) malloc(A.NrNzElts* sizeof(long));
    A.j = (long *) malloc(A.NrNzElts* sizeof(long));
    A.Pstart = (long *) malloc((P+1)* sizeof(long));
    X = (long int *) malloc(n* sizeof(long int));

    if ( A.i == NULL || A.j  == NULL || A.Pstart == NULL || X == NULL ){
        printf("Error\n");
        exit(1);
    }

    /* Fill second block of matrix with nonzeros */
    t= 0;
    for (i=k; i<2*k; i++){
        /* Insert the k nonzeros of a row */
        for (j=k; j<2*k; j++){
            A.i[t] = i;
            A.j[t] = j;
            t++;
        }
    }

    /* Fill third block of matrix with nonzeros */
    for (i=2*k; i<3*k; i++){ 
        /* Insert the k nonzeros of a row */
        for (j=2*k; j<3*k; j++){
            A.i[t] = i;
            A.j[t] = j;
            t++;
        } 
    }   

    nzp = k; 
    /* Procs 0, 1, ..., k-1 have k nonzeros each. */
    for (i=0; i<k; i++)
        A.Pstart[i] = i*nzp;
    nzp = k/2; 
    /* Procs k, k+1, ..., 3*k-1 have k/2 nonzeros each. */
    for (i=k; i<3*k; i++)  
        A.Pstart[i] = k*k + (i-k)*nzp;

    A.Pstart[P] = 2*k*k;

    Options.VectorPartition_Step3 = VecRandom;
    Options.VectorPartition_MaxNrLoops = 10;
    Options.VectorPartition_MaxNrGreedyImproves = 10;
    maxcom = DistributeVec(&A, X, ROW, &Options);
    
    if (maxcom < 0) {
        printf("Error\n");
        exit(1);
    }
    
    if (!CalcCom(&A, X, ROW, &ComVol, &MaxOut, &MaxIn, &MaxCompnts, &TotCompnts)) {
        printf("Error\n");
        exit(1);
    }

    /* Check result values. This is not a very strict check,
       because vector assignment and maxcom can vary. */
    if (ComVol != 2*(k-1)*k || maxcom !=  MAX(MaxIn,MaxOut) || TotCompnts != n){
        printf("Error\n");
        exit(1);
    }

    /* Check legality of vector distribution */
    for (i=0; i<n; i++){
        if (X[i] < 0 || X[i] >= P) {
            printf("Error\n");
            exit(1);
        }
    }

    /* Check second block */ 
    for (i=k; i<2*k; i++){
        if (X[i] >= k ) { 
            printf("Error\n"); 
            exit(1);
        } 
    }    

    /* Check third block */
    for (i=2*k; i<2*k +k/2; i++){
        if (X[i] < k || X[i]%2 == 1 ) {
            printf("Error\n");
            exit(1);
        }
    }   
    for (i=2*k +k/2; i<3*k; i++){
        if (X[i] < k || X[i]%2 == 0 ) {
            printf("Error\n");
            exit(1);
        }
    }    

    maxcom = DistributeVec(&A, X, COL, &Options);
    
    if (maxcom < 0) {
        printf("Error\n");
        exit(1);
    }
    
    if (!CalcCom(&A, X, COL, &ComVol, &MaxOut, &MaxIn, &MaxCompnts, &TotCompnts)) {
        printf("Error\n");
        exit(1);
    }
 
    /* Check result values. This is not a very strict check,
       because vector assignment and maxcom can vary. */
    if (ComVol != k || maxcom !=  MAX(MaxIn,MaxOut) || TotCompnts != n){
        printf("Error\n");
        exit(1); 
    }
 
    /* Check legality of vector distribution */
    for (i=0; i<n; i++){  
        if (X[i] < 0 || X[i] >= P) {
            printf("Error\n");
            exit(1);
        }
    }   
 
    /* Check second block */ 
    for (i=k; i<2*k; i++){
        if (X[i] != i-k ) {
            printf("Error\n");  
            exit(1);  
        }
    }    
 
    /* Check third block */
    for (i=2*k; i<3*k; i++){
        proc = 2*i -3*k;
        if (X[i] < proc || X[i] > proc + 1 ) {
            printf("Error\n");
            exit(1);
        }
    }    

    printf("OK\n");
    exit(0);

} /* end main */
