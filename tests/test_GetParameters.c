#include "Options.h"

int main(int argc, char **argv) {

    /* This test function should be called by 
           test_GetParameters pi.mtx 4 0.25 -SplitStrategy=onedimrow -SplitMethod simple
       to overrule the defaults for SplitStrategy and SplitMethod.
    */
    
    struct opts Options ;
    
    printf("Test GetParameters: ") ;
    
    SetDefaultOptions(&Options);
    SetOptionsFromFile(&Options, "test_GetParameters.defaults");
    GetParameters(&Options, argc, argv);
    
    /* Check option values */
    if ( Options.P !=  4 || 
         strcmp(Options.matrix, "pi.mtx")  ||
         Options.eps !=  0.25 ||
         Options.LoadbalanceStrategy != Constant ||
         Options.LoadbalanceAdjust != AdjustYes ||
         Options.SplitStrategy != OneDimRow ||
         Options.Alternate_FirstDirection != FirstDirRatio ||
         Options.SplitMethod != Simple ||
         Options.SquareMatrix_DistributeVectorsEqual != EqVecNo ||
         Options.SquareMatrix_DistributeVectorsEqual_AddDummies != DumYes ||
         Options.SymmetricMatrix_UseSingleEntry != SingleEntNo ||
         Options.SymmetricMatrix_SingleEntryType != ETypeLower ||
         Options.Coarsening_NrVertices != 200 ||
         Options.Coarsening_MaxNrVtxInMatch != 2 ||
         Options.Coarsening_StopRatio != 0.05 ||
         Options.Coarsening_VtxMaxFractionOfWeight != 0.2 ||
         Options.Coarsening_MatchingStrategy != MatchInprod ||
         Options.Coarsening_InprodMatchingOrder != DecreasingWgt ||
         Options.Coarsening_FineSwitchLevel != 2 ||
         Options.Coarsening_NetScaling != NetSclLinear ||
         Options.Coarsening_InprodScaling != IpSclMin ||
         Options.Coarsening_MatchIdenticalFirst != MatchIdYes ||
         Options.KLFM_InitPart_NrRestarts != 8 ||
         Options.KLFM_InitPart_MaxNrLoops != 25 ||
         Options.KLFM_InitPart_MaxNrNoGainMoves != 200 ||
         Options.KLFM_Refine_MaxNrLoops != 25 ||
         Options.KLFM_Refine_MaxNrNoGainMoves != 200 ||
         Options.VectorPartition_Step3 != VecRandom ||
         Options.VectorPartition_MaxNrLoops != 10 ||
         Options.VectorPartition_MaxNrGreedyImproves  != 10 ||
         Options.Seed != 99 ) {
        printf("Error\n") ;
        exit(1);
    }

    printf("OK\n") ;
    exit(0);

} /* end main */
