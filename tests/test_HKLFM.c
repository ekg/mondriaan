#include <stdlib.h>
#include <stdio.h>

#include "Graph.h"
#include "HKLFM.h"
#include "Match.h"

struct opts Options;

extern int HKLFM(struct biparthypergraph *pHG, long weightlo, long weighthi,
           int MaxNrLoops, int MaxNrNoGainMoves);

int main(int argc, char **argv) {

    struct biparthypergraph HG;

    long i, j, n, nz, t, P, weightlo, weighthi, v0, v1;

    printf("Test HKLFM: ");
    n = 75; /* parameter can be changed */
    nz = n*(n+1)/2;

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
        fprintf(stderr, "test_HKLFM(): Not enough memory!\n");
        printf("Error\n");
        exit(1);
    }

    /* Initialise partition of vertices */
    for (j=0; j<n-1; j++) 
        HG.V[j].partition = 0;
    HG.V[n-1].partition = 1;
    
    /* Initialise weight of vertices */
    for (j=0; j<n-1; j++) 
        HG.V[j].vtxwgt = 1;
    HG.V[n-1].vtxwgt = n;

    /* Initialise start and end of vertices */
    for (j=0; j<n; j++) {
        HG.V[j].iStart =  j*(j+1)/ 2;
        HG.V[j].iEnd =  (j+1)*(j+2)/ 2;
    }

    /* Initialise adjacency lists of vertices */
    for (j=0; j<n; j++)
        for (t=HG.V[j].iStart; t<HG.V[j].iEnd; t++) 
            HG.VtxAdjncy[t] = t - HG.V[j].iStart; 

    /* Initialise start and end of nets */
    for (i=0; i<n; i++) {
        HG.N[i].iStartP0 = n*(n+1)/2 - (n-i)*(n-i+1)/2;
        HG.N[i].iEnd = n*(n+1)/2 - (n-i-1)*(n-i)/2;
        HG.N[i].iStartP1 = HG.N[i].iEnd - 1;
    }

    /* Initialise adjacency lists of nets */        
    for (i=0; i<n; i++)
        for (t=HG.N[i].iStartP0; t<HG.N[i].iEnd; t++) 
            HG.NetAdjncy[t] = t - HG.N[i].iStartP0 + i; 

    /* Initialise gainbucket data structure */
    for (P=0; P<2; P++) {
        HG.GBVtx[P].NrBuckets  = 0;
        HG.GBVtx[P].Root = NULL;
    }    

    HG.OptComm = LONG_MAX; /* this is the first HKLFM run */
    HG.MinComm = LONG_MAX;
    
    /* Initialise partition weights and bounds */
    HG.WeightP[0] = n-1; /* the first n-1 vertices have weight 1 */
    HG.WeightP[1] = n;   /* the last vertex has weight n */
    weightlo = n-1; /* hence vertex n-1 cannot move into part 0 */
    weighthi = 2*n-1; /* hence all vertices can move to part 1 */

    /* Perform an HKLFM run with 1 pass through the vertices,
       and with at most 1 gainless move allowed in the pass */
    if (!HKLFM(&HG, weightlo, weighthi, 1, 1)) {
        printf("Error\n");
        exit(1);
    }

    /* Check hypergraph dimensions */
    if (HG.NrNets != n || HG.NrVertices != n || HG.NrPins != nz) {
        printf("Error\n");
        exit(1);
    }
        
    /* Check start and end of nets and net parts */
    for (i=0; i<n; i++) {
         if(HG.N[i].iStartP0 != n*(n+1)/2 - (n-i)*(n-i+1)/2 ||
             HG.N[i].iEnd != n*(n+1)/2 - (n-i-1)*(n-i)/2 ||
             HG.N[i].iStartP1 != HG.N[i].iStartP0) { /* all pins are
                                                         in part 1 */
             printf("Error\n");
             exit(1);
         }
    }

    /* Check range of net parts */
    for (i=0; i<n; i++)
        for (t=HG.N[i].iStartP0; t<HG.N[i].iEnd; t++) 
            if (HG.NetAdjncy[t] < i || HG.NetAdjncy[t] >= n) {
                printf("Error\n");
                exit(1);
            }
            
    /* Check partitions and locks of vertices. The vertices should all
       have moved to part 1 and they all should have been locked */
    for (j=0; j<n; j++) 
        if (HG.V[j].partition != 1 || HG.V[j].GBentry != NULL) {
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

    /* Check communication, weights, number of moves */
    if (HG.OptComm != 0 || HG.MinComm != 0 ||
         HG.WeightP[0] != 0 || HG.WeightP[1] != 2*n-1 ||
         HG.CurVtxLog != n-1 || HG.MinVtxLog != n-1) {

        printf("Error\n");
        exit(1);
    }
    
    /* Check saved partition */
    for (j=0; j<n; j++) 
        if (HG.OptPartVtx[j] != 1) {
            printf("Error\n");
            exit(1);
        }
        
    /* Check movelog. It should contain moves of vertices n-2, n-3, ..., 1, 0.
       Vertex n-1 should have been locked first, but not moved */
    for (j=0; j<n-1; j++) 
        if (HG.VtxMoveLog[j] != n-j-2) {
            printf("Error\n");
            exit(1);
        }

    printf("OK\n");
    exit(0);

} /* end main */
