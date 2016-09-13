#include "Graph.h"
#include "HKLFM.h"
#include "Match.h"

struct opts Options;

extern int ComputeInitialGains(struct biparthypergraph *pHG);

int main(int argc, char **argv) {

    struct biparthypergraph HG;

    long n, t, i, P;

    printf("Test ComputeInitialGains: ");
    n= 60; /* must be a multiple of 6 */

    HG.NrNets = n;
    HG.NrVertices = 2;
    HG.NrPins = n;

    HG.V = (struct vertex *) malloc(HG.NrVertices * sizeof(struct vertex));
    HG.N = (struct net *) malloc(HG.NrNets * sizeof(struct net));
    HG.VtxAdjncy = (long *) malloc(HG.NrPins * sizeof(long));
    HG.NetAdjncy = (long *) malloc(HG.NrPins * sizeof(long));
    if (HG.V == NULL ||  HG.N == NULL ||
         HG.VtxAdjncy == NULL || HG.NetAdjncy == NULL) {
        fprintf(stderr, "test_ComputeInitialGains(): Not enough memory!\n");
        printf("Error\n");
        exit(1);
    }

    /* Initialise partition of vertices */
    HG.V[0].partition = 0;
    HG.V[1].partition = 1;

    /* Initialise start and end of vertices */
    HG.V[0].iStart = 0;
    HG.V[0].iEnd = n/2;
    HG.V[1].iStart = n/2;
    HG.V[1].iEnd = n;

    /* Initialise adjacency lists of vertices */
    /* column 0 */
    for (t=0; t<n/6; t++) {
        i = 3*t;
        HG.VtxAdjncy[i] = 6*t + 2;
        HG.VtxAdjncy[i+1] = 6*t + 4;
        HG.VtxAdjncy[i+2] = 6*t + 5;
    }
    /* column 1 */
    for (t=0; t<n/6; t++) {
        i = 3*t;
        HG.VtxAdjncy[n/2 + i] = 6*t + 1;
        HG.VtxAdjncy[n/2 + i+1] = 6*t + 2;
        HG.VtxAdjncy[n/2 + i+2] = 6*t + 5;
    }

    /* Initialise start and end of nets */
    for (t=0; t<n/3; t++) {
        /* no nonzeros in row i = 3t */
        i = 3*t;
        HG.N[i].iStartP0 = i;
        HG.N[i].iEnd = i;
        /* 1 nonzero in row i+1 */
        HG.N[i+1].iStartP0 = i;
        HG.N[i+1].iEnd = i+1;
        /* 2 nonzeros in row i+2 */
        HG.N[i+2].iStartP0 = i+1;
        HG.N[i+2].iEnd = i+3;
    }

    /* Initialise adjacency lists of nets */
    for (t=0; t<n; t++) 
        HG.NetAdjncy[t] = (t+1)%2;

    /* Initialise gainbucket data structure */
    for (P=0; P<2; P++) {
        HG.GBVtx[P].NrBuckets  = 0;
        HG.GBVtx[P].Root = NULL;
    }

    ComputeInitialGains(&HG);

    /* Check hypergraph dimensions */
    if (HG.NrVertices != 2 || HG.NrNets != n || HG.NrPins != n) {
        printf("Error\n");
        exit(1);
    }

    /* Check gainbucket data structure */
    for (P=0; P<2; P++) {
        if (HG.GBVtx[P].NrBuckets  != 1 ||
             (HG.GBVtx[P].Root)->value != n/3 ||
             ((HG.GBVtx[P].Root)->entry)->vtxnr != P) {
            printf("Error\n");
            exit(1);
        }
    }

    /* Check current communication */
    if (HG.CurComm != n/3) {
        printf("Error\n");
        exit(1);
    }

    printf("OK\n");
    exit(0);

} /* end main */
