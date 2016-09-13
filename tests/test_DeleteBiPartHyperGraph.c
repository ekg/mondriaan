#include "Graph.h"

struct opts Options;

int main(int argc, char **argv) {

    struct biparthypergraph HG;

    long NrVertices, NrNets, NrPins;

    printf("Test DeleteBiPartHyperGraph: ");
    NrVertices = 12;
    NrNets = 5;
    NrPins = 101;

    if (!CreateNewBiPartHyperGraph(NrVertices, NrNets, NrPins, TRUE , 'C', &HG)) {
        printf("Error\n");
        exit(1);
    }
    
    if (!DeleteBiPartHyperGraph(&HG)) {
        printf("Error\n");
        exit(1);
    }

    /* Check result value  */
    if (HG.NrVertices != 0 ||
        HG.NrNets != 0 ||       
        HG.NrPins != 0 ||
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

    printf("OK\n");
    exit(0);

} /* end main */
