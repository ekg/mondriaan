#include "Options.h"

int main(int argc, char **argv) {

    struct opts Options ;
    
    Options.argc = argc ;
    Options.argv = argv ;
    printf("Test SetDefaults: ") ;
    
    SetDefaults(&Options) ;
    
    /* Check option values */
    if ( Options.LoadbalanceStrategy != Constant ||
         Options.LoadbalanceAdjust != AdjustNo ||
         Options.SplitStrategy != Alternate ||
         Options.Alternate_FirstDirection != Row ||
         Options.SplitMethod != Simple ||
         Options.SquareMatrix_DistributeVectorsEqual != EqVecNo ||
         Options.SquareMatrix_DistributeVectorsEqual_AddDummies != DumNo ||
         Options.SymmetricMatrix_UseSingleEntry != SingleEntNo ||
         Options.SymmetricMatrix_SingleEntryType != ETypeLower ||
         Options.Coarsening_NrVertices != 1 ||
         Options.Coarsening_MaxNrVtxInMatch != 2 ||
         Options.Coarsening_StopRatio != 0.25 ||
         Options.Coarsening_VtxMaxFractionOfWeight != 0.25 ||
         Options.Coarsening_MatchingStrategy != Random ||
         Options.Coarsening_InprodMatchingOrder != RandomOrder ||
         Options.Coarsening_FineSwitchLevel != 2 ||
         Options.KLFM_InitPart_NrRestarts != 10 ||
         Options.KLFM_InitPart_MaxNrLoops != 10 ||
         Options.KLFM_InitPart_MaxNrNoGainMoves != 10 ||
         Options.KLFM_Refine_MaxNrLoops != 10 ||
         Options.KLFM_Refine_MaxNrNoGainMoves != 10 ||
         Options.VectorPartition_Step3 != VecIncrease ||
         Options.VectorPartition_MaxNrLoops != 10 ||
         Options.VectorPartition_MaxNrGreedyImproves  != 10 ||
         Options.Seed != -1 ) {
        printf("Error\n") ;
        exit(1);
    }

    printf("OK\n") ;
    exit(0);

} /* end main */
