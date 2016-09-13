/*
This is a test of the l*(l - 1)-metric partitioning done by Mondriaan + PaToH.

The matrix

test_lambdalambda1.mtx

1 1 0 0 0 0 0
1 1 0 0 0 0 0
0 1 1 1 1 1 1

partitioned into 4 parts (ONEDIMCOL) with imbalance 1/7, has optimal partitionings
0 0 1 2 2 3 3 (volumes 3 and 12 for (l - 1) and l*(l - 1)) and
0 1 1 2 2 3 3 (volumes 4 and 10).

test_lambdalambda2.mtx

1 0 0 0 1 0 1 0 0 1
0 0 1 0 1 0 1 0 0 1
1 0 0 1 1 1 0 0 0 1
0 1 1 0 0 1 1 1 1 1
0 1 0 0 1 0 0 1 0 1
0 0 0 1 1 1 0 0 0 0
1 0 0 0 1 0 1 0 0 0
0 1 0 0 1 0 1 0 1 1
1 0 1 0 1 1 0 0 1 0
1 0 1 0 1 0 0 0 1 1
1 0 0 1 1 0 0 0 0 1
0 0 0 0 0 1 1 0 0 0
1 0 1 0 1 0 1 1 1 1
1 0 1 1 0 0 0 1 1 0
0 1 1 0 1 0 0 1 0 1
1 0 0 1 0 1 1 0 1 1
1 0 0 1 0 1 1 0 0 0
1 1 0 0 0 1 0 0 1 1
0 1 1 0 0 0 0 0 0 1
1 0 0 1 1 0 0 0 1 1
   l - 1 : 43 (154) -- 4 2 1 4 0 3 3 2 1 0
l*(l - 1): 150 (44) -- 1 3 2 4 0 4 3 2 1 0 

test_lambdalambda3.mtx

0 0 1 1 1 1 1 0 
1 0 1 1 1 0 1 0 
1 0 1 1 0 1 1 1 
1 0 1 1 1 0 0 1 
0 0 0 0 0 1 0 1 
0 1 0 1 0 0 1 1 
1 1 0 0 0 0 0 0 
1 1 0 0 0 0 0 1 
   l - 1 : 13 (44) -- 3 3 2 1 2 0 1 0 
l*(l - 1): 42 (14) -- 3 0 3 2 2 1 1 0 

Mondriaan should find these optima.
*/
#include "Mondriaan.h"

struct opts Options;

#define TRYMAT 2

void ReadMatrix(struct sparsematrix *A) {
    long i;

#if TRYMAT==0
    FILE *file = fopen("test_lambdalambda1.mtx", "r");
    
    Options.P = 4;
    Options.eps = 1.0/6.9;
#elif TRYMAT==1
    FILE *file = fopen("test_lambdalambda2.mtx", "r");
    
    Options.P = 5;
    Options.eps = 0.1;
#elif TRYMAT==2
    FILE *file = fopen("test_lambdalambda3.mtx", "r");
    
    Options.P = 4;
    Options.eps = 0.1;
#endif
    
    /* Read file from disk. */
    if (file == NULL || !MMReadSparseMatrix(file, A)) {
        printf("Error\n");
        exit(1);
    }
    
    fclose(file);
    
    /* Set weights. */
    A->MMTypeCode[0] = 'W';
    
    A->ColWeights = (long *)malloc(A->n*sizeof(long));
    
    if (A->ColWeights == NULL) {
        printf("Error\n");
        exit(1);
    }
    
    A->NrColWeights = A->n;
    
    for (i = 0; i < A->n; i++) A->ColWeights[i] = 1;
    
    /* Create processor array. */
    A->NrProcs = Options.P;
    A->Pstart = (long *)malloc((A->NrProcs + 1)*sizeof(long));
    
    if (A->Pstart == NULL) {
        printf("Error\n");
        exit(1);
    }
    
    A->Pstart[0] = 0;
    
    for (i = 1; i <= A->NrProcs; i++) A->Pstart[i] = A->NrNzElts;
}

void CheckVol(long *_Vol1, long *_Vol2, const long *Partition, const struct sparsematrix *A) {
    long *PartMark;
    long LargestPart, Vol1, Vol2;
    long i;
    double epsilon;
    
    PartMark = (long *)malloc(A->NrProcs*sizeof(long));
    
    if (PartMark == NULL) {
        printf("Error\n");
        exit(1);
    }
    
    /* Check imbalance. */
    for (i = 0; i < A->NrProcs; i++) PartMark[i] = 0;
    for (i = 0; i < A->n; i++) PartMark[Partition[i]]++;
    
    LargestPart = 0;
    for (i = 0; i < A->NrProcs; i++) LargestPart = (LargestPart < PartMark[i] ? PartMark[i] : LargestPart);
    
    epsilon = (double)(A->NrProcs*LargestPart - A->n)/(double)A->n;
    
    if (epsilon > Options.eps) {
        printf("Error\n");
        exit(1);
    }
    
    /* Calculate volumes. */
    Vol1 = 0;
    Vol2 = 0;
    
    for (i = 0; i < A->m; i++) {
        Vol1 += A->RowLambda[i] - 1;
        Vol2 += A->RowLambda[i]*(A->RowLambda[i] - 1);
    }
    
    *_Vol1 = Vol1;
    *_Vol2 = Vol2;
    
    /* Free memory. */
    free(PartMark);
}

int main(int argc, char **argv) {

    struct sparsematrix A;
    long ComVol = -1;
    long *Partition;
    long Vol11 = -1, Vol12 = -1, Vol21 = -1, Vol22 = -1;
    long i;

    printf("Test LambdaLambdaMinusOneMetric: ");

#ifndef USE_PATOH
    printf("Untested (no PaToH)\n");
    exit(0);
#endif
    
    /* Initialise relevant options */
    if (!SetDefaultOptions(&Options)) {
        printf("Error\n");
        exit(1);
    }
    
    Options.P = -1; /* set by ReadMatrix() */
    Options.eps = -1.0;
    Options.Seed = 12345;
    Options.SplitStrategy = OneDimCol;
    Options.Partitioner = PartPaToH;
    Options.SymmetricMatrix_UseSingleEntry = SingleEntNo;
    Options.SquareMatrix_DistributeVectorsEqual = EqVecNo;
    Options.Metric = MetricLambda;
    
    if (!ApplyOptions(&Options)) {
        printf("Error\n");
        exit(1);
    }
    
    /* Try the (l - 1)-metric. */
    Options.Metric = MetricLambda;
    ReadMatrix(&A);
    
    /* Allocate partitioning. */
    Partition = (long *)malloc(A.n*sizeof(long));
    
    if (Partition == NULL) {
        printf("Error\n");
        exit(1);
    }
    
    for (i = 0; i < A.n; i++) Partition[i] = 0;
    
    if (!DistributeMatrixMondriaan(&A, Options.P, Options.eps, &Options, NULL)) {
        printf("Error\n");
        exit(1);
    }
    
    ComVol = DistributeVec(&A, Partition, ROW, &Options);
    
    if (ComVol < 0) {
        printf("Error\n");
        exit(1);
    }
    
    CheckVol(&Vol11, &Vol12, Partition, &A);
    
    /* Free memory. */
    MMDeleteSparseMatrix(&A);
    
    /* Try the l*(l - 1)-metric. */
    Options.Metric = MetricLambdaLambdaMinusOne;
    ReadMatrix(&A);
    
    for (i = 0; i < A.n; i++) Partition[i] = 0;
    
    if (!DistributeMatrixMondriaan(&A, Options.P, Options.eps, &Options, NULL)) {
        printf("Error\n");
        exit(1);
    }
    
    ComVol = DistributeVec(&A, Partition, ROW, &Options);
    
    if (ComVol < 0) {
        printf("Error\n");
        exit(1);
    }
    
    CheckVol(&Vol21, &Vol22, Partition, &A);
    
    /*
    printf("Volumes   |   (l - 1) | l*(l - 1) |\n");
    printf("----------+-----------+-----------+\n");
    printf("(l - 1)   | % 9ld | % 9ld |\n", Vol11, Vol12);
    printf("l*(l - 1) | % 9ld | % 9ld |\n", Vol21, Vol22);
    printf("----------+-----------+-----------+\n");
    */
    
    /* Free memory. */
    free(Partition);
    MMDeleteSparseMatrix(&A);
    
    if (Vol11 >= Vol21 || Vol22 >= Vol12 || Vol11 < 0 || Vol12 < 0 || Vol21 < 0 || Vol22 < 0) {
        printf("Error\n");
        exit(1);
    }
    
    printf("OK\n");
    exit(0);

}
