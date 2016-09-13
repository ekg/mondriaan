#include "Graph.h"
#include "Match.h"

struct opts Options;

int main(int argc, char **argv) {

    struct biparthypergraph HG;

    long i, j, n, t, v;

    printf("Test MoveVtxInNetAdjncy: ");
    n= 13;

    /* dense hypergraph, corresponding to n by n dense matrix */
    HG.NrVertices = n;
    HG.NrNets = n;
    HG.NrPins = n*n;

    HG.V = (struct vertex *) malloc(HG.NrVertices * sizeof(struct vertex));
    HG.N = (struct net *) malloc(HG.NrNets * sizeof(struct net));
    HG.VtxAdjncy = (long *) malloc(HG.NrPins * sizeof(long));
    HG.NetAdjncy = (long *) malloc(HG.NrPins * sizeof(long));
    if (HG.V == NULL ||  HG.N == NULL || 
         HG.VtxAdjncy == NULL || HG.NetAdjncy == NULL) {
        fprintf(stderr, "test_MoveVtxInNetAdjncy(): Not enough memory!\n");
        printf("Error\n");
        exit(1);
    }

    /* Initialise vertex start/end */
    for (t=0; t<HG.NrVertices; t++) {
        HG.V[t].iStart = t*n;
        HG.V[t].iEnd = (t+1)*n;
    }

    /* Initialise net start/end */
    for (t=0; t<HG.NrNets; t++) {
        HG.N[t].iStartP0 = t*n; 
        HG.N[t].iStartP1 = (t+1)*n;
        HG.N[t].iEnd = (t+1)*n;
    }    

    /* Initialise each adjacency list to 0,1,2, ...,n-1 */ 
    for (i=0; i<n; i++) {
        for (j=0; j<n; j++) {
            t = i*n + j;
            HG.VtxAdjncy[t] = j;
            HG.NetAdjncy[t] = j;
        }
    }    

    v = 0;
    MoveVtxInNetAdjncy(&HG, v); 

    /* Check hypergraph dimensions */
    if (HG.NrVertices != n ||
        HG.NrNets != n ||
        HG.NrPins != n*n) {

        printf("Error\n");
        exit(1);
    }

    /* Check vertex values */
    for (t=0; t<HG.NrVertices; t++) {
        if (HG.V[t].iStart != t*n ||
            HG.V[t].iEnd != (t+1)*n ) {

            printf("Error\n");
            exit(1);
        }
    }

    /* Check net values */
    for (t=0; t<HG.NrNets; t++) {
        if (HG.N[t].iStartP0 != t*n ||
            HG.N[t].iStartP1 != t*n + n-1 ||
            HG.N[t].iEnd != (t+1)*n) {

            printf("Error\n");
            exit(1);
        }
    }

    /* Check pin values */
    for (i=0; i<n; i++) {
        for (j=0; j<n; j++) {
            t = i*n + j;
            if (HG.VtxAdjncy[t] != j) {
                printf("Error\n");
                exit(1);
            }
            /* 0 must be in last position */
            if ((j<n-1 && HG.NetAdjncy[t] == 0) ||
                (j==n-1 && HG.NetAdjncy[t] != 0)) {
                printf("Error\n");
                exit(1);
            }
        }
    }

    printf("OK\n");
    exit(0);

} /* end main */
