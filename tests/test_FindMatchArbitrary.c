#include "Graph.h"
#include "Match.h"

struct opts Options;

int main(int argc, char **argv) {

    struct biparthypergraph HG;
    struct contraction C;

    long i, j, n, t, v;
    int *Matched;

    printf("Test FindMatchArbitrary: ");
    n= 24; /* must be multiple of 4 */

    /* Hypergraph, corresponding to n by n checkerboard matrix
       with a(i,j) nonzero if i+j is even */
    HG.NrVertices = n;
    HG.NrNets = n;
    HG.NrPins = n*n/2;

    HG.V = (struct vertex *) malloc(HG.NrVertices * sizeof(struct vertex));
    HG.N = (struct net *) malloc(HG.NrNets * sizeof(struct net));
    HG.VtxAdjncy = (long *) malloc(HG.NrPins * sizeof(long));
    HG.NetAdjncy = (long *) malloc(HG.NrPins * sizeof(long));
    Matched = (int *) malloc(HG.NrVertices * sizeof(int));
    if (HG.V == NULL ||  HG.N == NULL || 
         HG.VtxAdjncy == NULL || HG.NetAdjncy == NULL || Matched == NULL) {
        fprintf(stderr, "test_FindMatchArbitrary(): Not enough memory!\n");
        printf("Error\n");
        exit(1);
    }

    /* Contraction, at initial stage with no vertices matched yet */
    C.NrMatches = 0;
    C.MaxNrVertices = n; /* all vertices can end up in the same group */
    C.MaxVtxWgt = n/2; 
    C.Match = (long *) malloc(HG.NrVertices * sizeof(long));
    C.Start = (long *) malloc((HG.NrVertices+1) * sizeof(long));
    if (C.Match == NULL ||  C.Start == NULL) {
        fprintf(stderr, "test_FindMatchArbitrary(): Not enough memory!\n");
        printf("Error\n");
        exit(1);
    }
    C.Start[0] = 0;

    /* Initialise vertices */
    for (t=0; t<HG.NrVertices; t++) {
        HG.V[t].vtxwgt =  1;
        HG.V[t].iStart = t*n/2; /* each column has n/2 nonzeros */
        HG.V[t].iEnd = (t+1)*n/2;
        Matched[t] =  FALSE;
    }
    HG.V[0].vtxwgt =  n; /* vertex 0 is very heavy, cannot match */

    /* Initialise nets */
    for (t=0; t<HG.NrNets; t++) {
        HG.N[t].iStartP0 = t*n/2; 
        HG.N[t].iStartP1 = (t+1)*n/2;
        HG.N[t].iEnd = (t+1)*n/2;
    }    

    /* Initialise each adjacency list to 0,2,4,...,n-2 (even rows/columns)
       or 1,3,5,...,n-1 (odd rows/columns) */ 
    for (j= 0; j<n; j++) {
        for (i=j%2; i<n; i +=2) {
            t = j*n/2 + i/2;
            HG.VtxAdjncy[t] = i;
        }
    }    
    for (i=0; i<n; i++) {
        for (j= i%2; j<n; j +=2) {
            t = i*n/2 + j/2;
            HG.NetAdjncy[t] = j;
        }
    }    

    v = n/2;
    FindMatchArbitrary(&HG, &C, v, Matched);

    /* Check hypergraph dimensions */
    if (HG.NrVertices != n ||
        HG.NrNets != n ||
        HG.NrPins != n*n/2) {

        printf("Error\n");
        exit(1);
    }

    /* Check Matched array */
    for (j=0; j<HG.NrVertices; j++) {
        if (((j%2==0) && j!=0 && (Matched[j] == FALSE)) ||
             (j == 0 && Matched[j]) ||
             ((j%2==1) && Matched[j])   ) {
            printf("Error\n");
            exit(1);
        }
    }

    /* Check number of matches */
    if (C.NrMatches != 1 ||
         C.MaxNrVertices != n ||
         C.MaxVtxWgt != n/2 ||
         C.Start[0] != 0 ||
         C.Start[1] != n/2-1) {
        printf("Error\n");
        exit(1);
    }

    /* Check matches in contraction */
    for (t= C.Start[0]; t < C.Start[1]; t++) {
        j = C.Match[t];
        if (j == 0 || (j%2==1) || Matched[j] == FALSE) {
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
