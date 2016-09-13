#include "Graph.h"
#include "Match.h"

struct opts Options;

int main(int argc, char **argv) {

    struct biparthypergraph HG;
    struct contraction C;

    long i, j, n, t, v, *Visited, *Inprod;
    int *Matched;
    double *ScInprod;

    printf("Test FindMatchInprod: ");
    n= 40; 

    /* Hypergraph, corresponding to n by n upper triangular matrix */
    HG.NrVertices = n;
    HG.NrNets = n;
    HG.NrPins = n*(n+1)/2;

    HG.V = (struct vertex *) malloc(HG.NrVertices * sizeof(struct vertex));
    HG.N = (struct net *) malloc(HG.NrNets * sizeof(struct net));
    HG.VtxAdjncy = (long *) malloc(HG.NrPins * sizeof(long));
    HG.NetAdjncy = (long *) malloc(HG.NrPins * sizeof(long));
    Matched = (int *) malloc(HG.NrVertices * sizeof(int));
    Inprod =  (long *) malloc(HG.NrVertices * sizeof(long));
    Visited =  (long *) malloc(HG.NrVertices * sizeof(long));
    ScInprod =  (double *) malloc(HG.NrVertices * sizeof(double));

    if (HG.V == NULL || HG.N == NULL || 
         HG.VtxAdjncy == NULL || HG.NetAdjncy == NULL || Matched == NULL ||
         Inprod == NULL || Visited == NULL || ScInprod == NULL) {
        fprintf(stderr, "test_FindMatchInprod(): Not enough memory!\n");
        printf("Error\n");
        exit(1);
    }

    /* Contraction, at initial stage with no vertices matched yet */
    C.NrMatches = 0;
    C.MaxNrVertices = 2; /* pairwise matching */
    C.MaxVtxWgt = n*n;     /* no limit on weight */
    C.Match = (long *) malloc(HG.NrVertices * sizeof(long));
    C.Start = (long *) malloc((HG.NrVertices+1) * sizeof(long));
    if (C.Match == NULL ||  C.Start == NULL) {
        fprintf(stderr, "test_FindMatchInprod(): Not enough memory!\n");
        printf("Error\n");
        exit(1);
    }
    C.Start[0] = 0;

    /* Initialise vertices */
    for (t=0; t<HG.NrVertices; t++) {
        HG.V[t].vtxwgt =  t+1;
        HG.V[t].iStart = t*(t+1)/2; /* 1+2+3+...+t */
        HG.V[t].iEnd = (t+1)*(t+2)/2;
        Matched[t] =  FALSE;
    }

    /* Initialise nets */
    for (t=0; t<HG.NrNets; t++) {
        HG.N[t].iStartP0 = t*n - (t-1)*t/2; /* n+(n-1)+...+(n-t+1) */
        HG.N[t].iStartP1 = (t+1)*n - t*(t+1)/2; /* n+(n-1)+...+(n-t) */
        HG.N[t].iEnd = (t+1)*n - t*(t+1)/2;
    }    

    /* Initialise each vertex adjacency list j to 0, 1, 2, 3, ..., j */
    for (j= 0; j<n; j++) {
        for (i=0; i<=j; i++) {
            t = j*(j+1)/2 + i;
            HG.VtxAdjncy[t] = i;
        }
    }    
    /* Initialise each net adjacency list i to i, i+1, ..., n-1 */
    for (i=0; i<n; i++) {
        for (j=i; j<n; j++) {
            t = i*n - (i-1)*i/2 + j-i;
            HG.NetAdjncy[t] = j;
        }
    }    

    /* Initialise Inprod and ScInprod */
    for (j= 0; j<n; j++) {
        Inprod[j] = 0;
        ScInprod[j] = 0;
    }

    Options.Coarsening_NetScaling = NoNetScaling;
    Options.Coarsening_MatchIdenticalFirst = MatchIdYes;
    Options.Coarsening_InprodScaling = IpSclCos;

    v = n-2;
    FindMatchInprod(&HG, &C, v, Matched, Visited, Inprod, ScInprod, &Options);

    /* Check hypergraph dimensions */
    if (HG.NrVertices != n ||
        HG.NrNets != n ||
        HG.NrPins != n*(n+1)/2) {

        printf("Error\n");
        exit(1);
    }

    /* Check Matched, Inprod, ScInprod arrays */
    for (j=0; j<HG.NrVertices; j++) {
        if (((j >= n-2) && (Matched[j] == FALSE)) ||
             (j < n-2 && Matched[j]) ||
             Inprod[j] != 0 || ScInprod[j] != 0) {
            printf("Error\n");
            exit(1);
        }
    }

    /* Check number of matches */
    if (C.NrMatches != 1 ||
         C.MaxNrVertices != 2 ||
         C.MaxVtxWgt != n*n ||
         C.Start[0] != 0 ||
         C.Start[1] != 2) {
        printf("Error\n");
        exit(1);
    }

    /* Check matches in contraction */
    for (t= C.Start[0]; t < C.Start[1]; t++) {
        j = C.Match[t];
        if (j < n-2 || Matched[j] == FALSE) {
            printf("Error\n");
            exit(1);
        }
    }

    /* Reset all matched vertices */
    for (t= C.Start[0]; t < C.Start[1]; t++) {
        j = C.Match[t];
        Matched[j] = FALSE; 
    }

    /* Check reset Matched array */
    for (j=0; j<HG.NrVertices; j++) {
        if (Matched[j])  {
            printf("Error\n");
            exit(1);
        }
    }

    printf("OK\n");
    exit(0);

} /* end main */
