#include "DistributeMat.h"
#include "DistributeVecLib.h"

struct opts Options;

extern int SplitMatrixKLFM(struct sparsematrix *pT, int k, int i, int dir, 
                     long weightlo, long weighthi, const struct opts *pOptions);

int main(int argc, char **argv) {

    struct sparsematrix A;
    long n, i, j, t, P, weightlo, weighthi, weight0, weight1,
         ComVol, MaxOut, MaxIn, MaxCompnts, TotCompnts;

    printf("Test SplitMatrixKLFM: ");
    n = 33; /* n by n dense matrix A, n odd */
    A.m = n;
    A.n = n;
    A.NrNzElts = n*n;
    P = 2; /* maximum number of parts */

    A.i = (long *)malloc(A.NrNzElts*sizeof(long));
    A.j = (long *)malloc(A.NrNzElts*sizeof(long));
    A.Pstart = (long *)malloc((P+1)*sizeof(long));
    A.RowLambda = (int *)malloc(A.m*sizeof(int));
    A.ColLambda = (int *)malloc(A.n*sizeof(int));
    A.RowMark = (int *)malloc(A.m*sizeof(int));
    A.ColMark = (int *)malloc(A.n*sizeof(int));
    
    if (!SetDefaultOptions(&Options)) {
        printf("Error\n");
        exit(1);
    }

    if (A.i == NULL || A.j  == NULL || A.Pstart  == NULL) {
        printf("Error\n");
        exit(1);
    }

    /* Fill matrix with nonzeros */
    t= 0;
    for (i=0; i<n; i++) {
        for (j=0; j<n; j++) {
            A.i[t] = i;
            A.j[t] = j;
            t++;
        }
    }

    A.MMTypeCode[0]='M'; /* matrix */
    A.MMTypeCode[1]='C'; /* coordinate scheme */
    A.MMTypeCode[2]='P'; /* pattern only */
    A.MMTypeCode[3]='G'; /* general, no symmetry */
    A.NrDummies = 0;
    A.dummy = NULL;
    A.Pstart[0] = 0;
    A.Pstart[1] = A.NrNzElts;

    weightlo = n*(n-1)/2;
    weighthi = n*(n+1)/2;

    /* Initialise relevant options */
    Options.SplitStrategy = OneDimRow;
    Options.Coarsening_NrVertices = 3;
    Options.Coarsening_MaxNrVtxInMatch = 4;
    Options.Coarsening_StopRatio = 0.1;
    Options.Coarsening_VtxMaxFractionOfWeight = 0.1;
    Options.Coarsening_MatchingStrategy = MatchRandom;
    Options.Coarsening_InprodMatchingOrder = IncreasingDegree;
    Options.Coarsening_FineSwitchLevel = 2;
    Options.Coarsening_NetScaling = NoNetScaling;
    Options.Coarsening_InprodScaling = IpSclMin;
    Options.Coarsening_MatchIdenticalFirst = MatchIdNo;
    Options.KLFM_InitPart_NrRestarts = 10;
    Options.KLFM_InitPart_MaxNrLoops = 10;
    Options.KLFM_InitPart_MaxNrNoGainMoves = 50;
    Options.KLFM_Refine_MaxNrLoops = 10;
    Options.KLFM_Refine_MaxNrNoGainMoves = 50;

    if (!SplitMatrixKLFM(&A, 1, 0, ROW, weightlo, weighthi, &Options)) {
        printf("Error\n");
        exit(1);
    }

    /* Check nonzeros */
    if (A.Pstart[0] != 0 || A.Pstart[1] != n*(n-1)/2 ||
         A.Pstart[2] != A.NrNzElts) {
        printf("Error\n");
        exit(1);
    }

    /* Check part weights */
    weight0 = ComputeWeight(&A, A.Pstart[0], A.Pstart[1]-1, NULL, &Options);
    weight1 = ComputeWeight(&A, A.Pstart[1], A.Pstart[2]-1, NULL, &Options);

    if (weight0 != weightlo || weight1 != weighthi || weight0 < 0 || weight1 < 0) {
        printf("Error\n");
        exit(1);
    }

    /* Check communication */
    if (!CalcCom(&A, NULL, ROW, &ComVol, &MaxOut, &MaxIn, &MaxCompnts, &TotCompnts)) {
        printf("Error\n");
        exit(1);
    }
    if (ComVol !=  n) {
        printf("Error\n");
        exit(1);
    }

    printf("OK\n");
    exit(0);

} /* end main */
