#include "Graph.h"
#include "HKLFM.h"
#include "Match.h"

struct opts Options;

int main(int argc, char **argv) {

    struct biparthypergraph HG, CHG;
    struct contraction C;

    long i, j, t, n, NrVertices, NrNets, NrPins, max, min;
    int *Mark;

    printf("Test CoarsenGraph: ");
    n = 10; /* must be even */
    NrVertices = n;
    NrNets = n;
    NrPins = 2*n;

    if (!CreateNewBiPartHyperGraph(NrVertices, NrNets, NrPins, FALSE, 'P', &HG)) {
        printf("Error\n");
        exit(1);
    }

    /* Initialise vertex values */
    for (t=0; t<NrVertices; t++) {
        HG.V[t].vtxwgt = 2;
        HG.V[t].iStart = 2*t;
        HG.V[t].iEnd = 2*(t+1);
    }

    /* Initialise net values */
    for (t=0; t<NrNets; t++) {
        HG.N[t].iStartP0 = 2*t;
        HG.N[t].iEnd = 2*(t+1);
        HG.N[t].dir = ROW; 
    }    

    /* Initialise pin values */ 
    for (t=0; t<NrPins; t++) {
        HG.VtxAdjncy[t] = 2*(t/4) + t%2;
        HG.NetAdjncy[t] = HG.VtxAdjncy[t];
    }    

    /* Initialise relevant options */ 
    Options.Coarsening_NetScaling = NetSclLinear;
    Options.Coarsening_MatchIdenticalFirst = MatchIdYes;
    Options.Coarsening_InprodScaling = IpSclMax;
    Options.Coarsening_MatchingStrategy = MatchInprod;
    Options.Coarsening_InprodMatchingOrder = RandomOrder;
    Options.Coarsening_MaxNrVtxInMatch = 2;
    Options.Coarsening_VtxMaxFractionOfWeight = 0.2;
    Options.Coarsening_FineSwitchLevel = 0;
    Options.Coarsening_NrVertices = 2;
    Options.Coarsening_StopRatio = 0.05;
    
    if (!CoarsenGraph(&HG, &CHG, &C, 0, &Options)) {
        printf("Error\n");
        exit(1);
    }
    
    /* Check original hypergraph HG */
    if (HG.NrVertices != NrVertices ||
        HG.NrNets != NrNets ||       
        HG.NrPins != NrPins ) {

        printf("Error\n");
        exit(1);
    }

    /* Check vertex values of HG */
    for (t=0; t<NrVertices; t++) {
        if (HG.V[t].vtxwgt != 2 ||
            HG.V[t].iStart != 2*t ||
            HG.V[t].iEnd != 2*(t+1) ) {

            printf("Error\n");
            exit(1);
        }
    }
    
    /* Check net values of HG */
    for (t=0; t<NrNets; t++) {
        if (HG.N[t].iStartP0 != 2*t ||
            HG.N[t].iEnd != 2*(t+1) || 
            HG.N[t].dir != ROW ) {
 
            printf("Error\n");
            exit(1); 
        }
    }    

        
    /* Check contracted hypergraph HG */
    if (CHG.NrVertices != NrVertices/2 ||
        CHG.NrNets != NrNets ||       
        CHG.NrPins != NrPins/2 ||
        CHG.CurComm != 0 ||
        CHG.MinComm != 0 ||
        CHG.OptComm != 0 ||
        CHG.SplitDir != ROW ||
        CHG.CurVtxLog != 0 ||
        CHG.MinVtxLog != 0 ||
        CHG.WeightP[0] != 0 ||
        CHG.WeightP[1] != 0 ||
        CHG.GBVtx[0].NrBuckets != 0 ||
        CHG.GBVtx[1].NrBuckets != 0 ||
        CHG.GBVtx[0].Root != NULL ||
        CHG.GBVtx[1].Root != NULL ) {

        printf("Error\n");
        exit(1);
    }

    /* Allocate boolean marker array for checking unordered values */
    Mark = (int *) malloc(n * sizeof(int));
    if ( Mark == NULL ) {
        printf("Error\n");
        fprintf(stderr, "test_CoarsenGraph(): Not enough memory!\n");
        exit(1);
    }
    
    for (t=0; t<n; t++)
        Mark[t] = FALSE;

    /* Check vertex values of CHG */
    for (t=0; t<NrVertices/2; t++) {
        if (CHG.V[t].vtxwgt != 4 ||
            CHG.V[t].iEnd - CHG.V[t].iStart != 2 ||
            CHG.V[t].partition != 0 ||
            CHG.V[t].GBentry != NULL ) {

            printf("Error\n");
            exit(1);
        }
    }

    /* Check net values of CHG */
    for (t=0; t<NrNets; t++) {
        if (CHG.N[t].iEnd - CHG.N[t].iStartP0 != 1 ||
            CHG.N[t].dir != ROW ) {
 
            printf("Error\n");
            exit(1); 
        }
    }

    /* Check pin values of CHG */
    for (j=0; j<NrVertices/2; j++)
        for (t=CHG.V[j].iStart; t<CHG.V[j].iEnd; t++)
            Mark[CHG.VtxAdjncy[t]] = TRUE;

    for (t=0; t<n; t++)
        if (Mark[t] == FALSE ) {
            printf("Error\n");
            exit(1);
        }

    /* Now, Mark = TRUE */

    for (i=0; i<NrNets; i++)
        for (t=CHG.N[i].iStartP0; t<CHG.N[i].iEnd; t++)
            Mark[CHG.NetAdjncy[t]] = FALSE;

    /* Check only the first half of the marker array */
    for (t=0; t<n/2; t++)
        if (Mark[t] == TRUE ) {
            printf("Error\n");
            exit(1);
        }

    /* Reset the second half as well */
    for (t=n/2; t<n; t++)
        Mark[t] = FALSE;

    /* Now, Mark = FALSE */

    /* Check contraction C */
    if (C.NrMatches != NrVertices/2 ||
        C.MaxNrVertices != 2 ||
        C.MaxVtxWgt != 0.2 * (2*n) ) {

        printf("Error\n");
        exit(1);
    }

    /* Check C.Match */
    for (t=0; t<C.NrMatches; t++) {
        max  = MAX(C.Match[2*t], C.Match[2*t+1]);
        min  = MIN(C.Match[2*t], C.Match[2*t+1]);
        Mark[max]= TRUE;
        Mark[min]= TRUE;
        if (max-min != 1) {
            printf("Error\n");
            exit(1);
        }
    }

    for (t=0; t<n; t++)
        if (Mark[t] == FALSE) {
            printf("Error\n");
            exit(1);
        }

    /* Now, Mark = TRUE */

    /* Check C.start */
    for (t=0; t<C.NrMatches; t++) {
        if (C.Start[t]%2 == 1 || C.Start[t+1] - C.Start[t] != 2) {
            printf("Error\n");
            exit(1);
        }
        Mark[C.Match[C.Start[t]]] = FALSE;
        Mark[C.Match[C.Start[t]+1]] = FALSE;
    }

    for (t=0; t<n; t++)
        if (Mark[t] == TRUE) {            
            printf("Error\n");
            exit(1);
        }

    printf("OK\n");
    exit(0);

} /* end main */
