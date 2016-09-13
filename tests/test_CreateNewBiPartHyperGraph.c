#include "Graph.h"

struct opts Options;

int main(int argc, char **argv) {

    struct biparthypergraph HG;

    long t, NrVertices, NrNets, NrPins;

    printf("Test CreateNewBiPartHyperGraph: ");
    NrVertices = 17;
    NrNets = 95;
    NrPins = 200;

    if (!CreateNewBiPartHyperGraph(NrVertices, NrNets, NrPins, TRUE , 'C', &HG)) {
        printf("Error\n");
        exit(1);
    }

    /* Check result value  */
    if (HG.NrVertices != NrVertices ||
        HG.NrNets != NrNets ||       
        HG.NrPins != NrPins ||
        HG.CurComm != 0 ||
        HG.MinComm != 0 ||
        HG.OptComm != 0 ||
        HG.SplitDir != ROW ||
        HG.CurVtxLog != 0 ||
        HG.MinVtxLog != 0 ||
        HG.WeightP[0] != 0 ||
        HG.WeightP[1] != 0 ||
        HG.GBVtx[0].NrBuckets != 0 ||
        HG.GBVtx[1].NrBuckets != 0 ||
        HG.GBVtx[0].Root != NULL ||
        HG.GBVtx[1].Root != NULL ) {

        printf("Error\n");
        exit(1);
    }

    /* Check vertex values */
    for (t=0; t<NrVertices; t++){
        if (HG.V[t].vtxwgt != 0 ||
            HG.V[t].iStart != 0 ||
            HG.V[t].iEnd != 0 ||
            HG.V[t].partition != 0 ||
            HG.V[t].GBentry != NULL ||
            HG.Vtx2MatIndex[t] != 0 ||
            HG.OptPartVtx[t] != 0) {

            printf("Error\n");
            exit(1);
        }
    }

    /* Check net values */
    for (t=0; t<NrNets; t++){
        if (HG.N[t].netwgt != 0 ||
            HG.N[t].iStartP0 != 0 ||
            HG.N[t].iStartP1 != 0 ||
            HG.N[t].iEnd != 0 || 
            HG.N[t].dir != ROW || 
            HG.Net2MatIndex[t] != 0) {
 
            printf("Error\n");
            exit(1); 
        }
    }    

    /* Check pin values */ 
    for (t=0; t<NrPins; t++){
        if (HG.VtxAdjncy[t] != 0 || 
            HG.NetAdjncy[t] != 0 || 
            HG.MatReValue[t] != 0.0 || 
            HG.MatImValue[t] != 0.0) {
 
            printf("Error\n");
            exit(1);
        }
    }    

    printf("OK\n");
    exit(0);

} /* end main */
