#include "Graph.h"
#include "HKLFM.h"
#include "Match.h"

struct opts Options;

int main(int argc, char **argv) {

    struct biparthypergraph HG;

    long i, j, n, nz, t, P, weightlo, weighthi, n0, n1, v0, v1;

    printf("Test RunHKLFM: ");
    n = 12; /* parameter can be changed, must be even */
    nz = n*(n-1) + 2;

    HG.NrNets = n;
    HG.NrVertices = n;
    HG.NrPins = nz;

    HG.V = (struct vertex *) malloc(HG.NrVertices * sizeof(struct vertex));
    HG.OptPartVtx = (int *) malloc(HG.NrVertices * sizeof(int));
    HG.VtxMoveLog = (long *) malloc(HG.NrVertices * sizeof(long));
    HG.N = (struct net *) malloc(HG.NrNets * sizeof(struct net));
    HG.VtxAdjncy = (long *) malloc(HG.NrPins * sizeof(long));
    HG.NetAdjncy = (long *) malloc(HG.NrPins * sizeof(long));
    if (HG.V == NULL || HG.OptPartVtx == NULL || HG.VtxMoveLog == NULL ||
         HG.N == NULL || HG.VtxAdjncy == NULL || HG.NetAdjncy == NULL) {
        fprintf(stderr, "test_RunHKLFM(): Not enough memory!\n");
        printf("Error\n");
        exit(1);
    }

    /* Initialise weight of vertices */
    for (j=0; j<n; j++) 
        HG.V[j].vtxwgt = 1;

    /* Initialise start and end of vertices */
    for (j=0; j<2; j++) {
        HG.V[j].iStart =  j*n;
        HG.V[j].iEnd =  (j+1)*n;
    }
    for (j=2; j<n; j++) {
        HG.V[j].iStart =  j*(n-1) + 2;
        HG.V[j].iEnd =  (j+1)*(n-1) + 2;
    }

    /* Initialise adjacency lists of vertices */
    for (j=0; j<n; j++)
        for (t=HG.V[j].iStart; t<HG.V[j].iEnd; t++) 
            HG.VtxAdjncy[t] = t - HG.V[j].iStart; 

    /* Initialise start and end of nets */
    for (i=0; i<n-1; i++) {
        HG.N[i].iStartP0 = i*n;
        HG.N[i].iEnd = (i+1)*n;
    }
    HG.N[n-1].iStartP0 = (n-1)*n;
    HG.N[n-1].iEnd = (n-1)*n + 2;

    /* Initialise adjacency lists of nets */        
    for (i=0; i<n; i++)
        for (t=HG.N[i].iStartP0; t<HG.N[i].iEnd; t++) 
            HG.NetAdjncy[t] = t - HG.N[i].iStartP0; 

    /* Initialise gainbucket data structure */
    for (P=0; P<2; P++) {
        HG.GBVtx[P].NrBuckets  = 0;
        HG.GBVtx[P].Root = NULL;
    }    

    /* Initialise partition bounds */
    weightlo = n/2 + 1; /* hence always at least one vertex can move */
    weighthi = n/2 + 1; 

    /* Perform several HKLFM runs, each with several iterations */
    Options.KLFM_InitPart_NrRestarts = 10;
    Options.KLFM_InitPart_MaxNrLoops = 10;
    Options.KLFM_InitPart_MaxNrNoGainMoves = n;

    RunHKLFM(&HG, weightlo, weighthi, FALSE, &Options);

    /* Check hypergraph dimensions */
    if (HG.NrNets != n || HG.NrVertices != n || HG.NrPins != nz) {
        printf("Error\n");
        exit(1);
    }
        
    /* Check start and end of first n-1 nets and their net parts */
    for (i=0; i<n-1; i++) {
         if(HG.N[i].iStartP0 != i*n ||
             HG.N[i].iEnd != (i+1)*n ||
             HG.N[i].iStartP1 < HG.N[i].iStartP0 + n/2 - 1 ||
             HG.N[i].iStartP1 > HG.N[i].iStartP0 + n/2 + 1) {

             printf("Error\n");
             exit(1);
         }
    }

    /* Check range of net adjacencies */
    for (i=0; i<n; i++)
        for (t=HG.N[i].iStartP0; t<HG.N[i].iEnd; t++) 
            if (HG.NetAdjncy[t] < 0 || HG.NetAdjncy[t] >= n) {
                printf("Error\n");
                exit(1);
            }
            
    /* Check gainbucket data structures. They should be empty */
    v0 = GainBucketGetMaxValVertexNr(&(HG.GBVtx[0])); 
    v1 = GainBucketGetMaxValVertexNr(&(HG.GBVtx[1]));
    for (P=0; P<2; P++) 
        if (HG.GBVtx[P].NrBuckets != 0 || HG.GBVtx[P].Root != NULL) {
            printf("Error\n");
            exit(1);
        }    

    if (v0 != LONG_MIN || v1 != LONG_MIN) {
        printf("Error\n");
        exit(1);
    }

    /* Check communication and weights */
    if (HG.OptComm != n-1 || HG.MinComm != n-1 ||
         HG.WeightP[0] > n/2 + 1 || HG.WeightP[1] > n/2 + 1 ||
         HG.WeightP[0] +  HG.WeightP[1] != n) {

        printf("Error\n");
        exit(1);
    }
    
    /* Check saved partition */
    n0 = 0;
    n1 = 0;
    for (j=0; j<n; j++) {
        if (HG.OptPartVtx[j] == 0)
            n0++;
        else if (HG.OptPartVtx[j] == 1) 
            n1++; 
        else {
            printf("Error\n");
            exit(1);
        }
    }

    if (n0 != HG.WeightP[0] || n1 != HG.WeightP[1]) {
        printf("Error\n");
        exit(1);
    }
        
    printf("OK\n");
    exit(0);

} /* end main */
