#include "Graph.h"
#include "HKLFM.h"
#include "Match.h"

struct opts Options;

int main(int argc, char **argv) {

    struct biparthypergraph HG;

    long i, j, m, n, nz, t, P, weightlo, weighthi, n0, n1, v0, v1;

    printf("Test RunMLGraphPart: ");
    m = 62; /* parameter can be changed, must be even */
    n = 92; /* parameter can be changed, must be even */
    nz = m*n/2;
    
    if (!SetDefaultOptions(&Options)) {
        printf("Error\n");
        exit(1);
    }

    if (!CreateNewBiPartHyperGraph(n, m, nz, FALSE, 'P', &HG)) {
        printf("Error\n");
        exit(1);
    }

    /* Test hypergraph represents an m by n checkerboard matrix,
       with nonzeros in positions (i,j) with i+j even */

    /* Initialise weight of vertices */
    for (j=0; j<n; j++) 
        HG.V[j].vtxwgt = m/2;

    /* Initialise start and end of vertices */
    for (j=0; j<n; j++) {
        HG.V[j].iStart =  j*(m/2);
        HG.V[j].iEnd =  (j+1)*(m/2);
    }

    /* Initialise adjacency lists of vertices */
    for (j=0; j<n; j++)
        for (t=HG.V[j].iStart; t<HG.V[j].iEnd; t++) 
            HG.VtxAdjncy[t] = 2*(t - HG.V[j].iStart) + j%2; 

    /* Initialise start and end of nets */
    for (i=0; i<m; i++) {
        HG.N[i].iStartP0 = i*(n/2);
        HG.N[i].iEnd = (i+1)*(n/2);
    }

    /* Initialise adjacency lists of nets */        
    for (i=0; i<m; i++)
        for (t=HG.N[i].iStartP0; t<HG.N[i].iEnd; t++) 
            HG.NetAdjncy[t] = 2* (t - HG.N[i].iStartP0) + i%2; 

    /* Initialise partition bounds */
    weightlo = nz/2 + m/2; /* hence always at least one vertex can move */
    weighthi = nz/2 + m/2; 

    /* Initialise relevant options */ 
    Options.Coarsening_NetScaling = NoNetScaling;
    Options.Coarsening_MatchIdenticalFirst = MatchIdNo;
    Options.Coarsening_InprodScaling = IpSclMin;
    Options.Coarsening_MatchingStrategy = MatchRandom;
    Options.Coarsening_InprodMatchingOrder = DecreasingDegree;
    Options.Coarsening_MaxNrVtxInMatch = 2;
    Options.Coarsening_VtxMaxFractionOfWeight = 0.2;
    Options.Coarsening_FineSwitchLevel = 2;
    Options.Coarsening_NrVertices =  10;
    Options.Coarsening_StopRatio = 0.05;

    Options.KLFM_InitPart_NrRestarts = 10;
    Options.KLFM_InitPart_MaxNrLoops = 10;
    Options.KLFM_InitPart_MaxNrNoGainMoves = n;
    Options.KLFM_Refine_MaxNrLoops = 10;
    Options.KLFM_Refine_MaxNrNoGainMoves = n;

    if (!RunMLGraphPart(&HG, weightlo, weighthi, &Options)) {
        printf("Error\n");
        exit(1);
    }

    /* Check hypergraph dimensions */
    if (HG.NrNets != m || HG.NrVertices != n || HG.NrPins != nz ) {
        printf("Error\n");
        exit(1);
    }
        
    /* Check vertex values of HG */
    for (j=0; j<n; j++){
        if (HG.V[j].vtxwgt != m/2 ||
            HG.V[j].iStart != j*(m/2) ||
            HG.V[j].iEnd != (j+1)*(m/2) ) {

            printf("Error\n");
            exit(1);
        }
    }

    /* Check start and end of nets and their net parts.
       iStartP1-iStartP0 must be the same for all even rows, and for all odd rows. */
    for (i=0; i<m; i++) {
         if( HG.N[i].iStartP0 != i*(n/2) ||
             HG.N[i].iEnd != (i+1)*(n/2) ||
             HG.N[i].iStartP1 - HG.N[i].iStartP0 != HG.N[i%2].iStartP1 - HG.N[i%2].iStartP0) {

             printf("Error\n");
             exit(1);
         }
    }

    /* Check range of net adjacencies */
    for (i=0; i<m; i++)
        for (t=HG.N[i].iStartP0; t<HG.N[i].iEnd; t++) 
            if (HG.NetAdjncy[t] < 0 || HG.NetAdjncy[t] >= n ) {
                printf("Error\n");
                exit(1);
            }
            
    /* Check gainbucket data structures. They should be empty */
    v0 = GainBucketGetMaxValVertexNr(&(HG.GBVtx[0])); 
    v1 = GainBucketGetMaxValVertexNr(&(HG.GBVtx[1]));
    for (P=0; P<2; P++) 
        if ( HG.GBVtx[P].NrBuckets != 0 || HG.GBVtx[P].Root != NULL ) {
            printf("Error\n");
            exit(1);
        }    

    if ( v0 != LONG_MIN || v1 != LONG_MIN ) {
        printf("Error\n");
        exit(1);
    }

    /* Check weights */
    if ( HG.WeightP[0] > weightlo || HG.WeightP[1] > weighthi ||
         HG.WeightP[0] +  HG.WeightP[1] != m*n/2 ) {

        printf("Error\n");
        exit(1);
    }
    
    /* Check saved partition */
    n0 = 0;
    n1 = 0;
    for (j=0; j<n; j++) {
        if ( HG.OptPartVtx[j] == 0 )
            n0++;
        else if ( HG.OptPartVtx[j] == 1 ) 
            n1++; 
        else {
            printf("Error\n");
            exit(1);
        }
    }

    if (HG.WeightP[0] != n0*m/2 || HG.WeightP[1] != n1*m/2 ){
        printf("Error\n");
        exit(1);
    }
        
    printf("OK\n");
    exit(0);

}
