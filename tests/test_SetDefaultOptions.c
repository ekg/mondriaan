#include "Options.h"

void CheckOptions(struct opts Options) {
    if ( Options.LoadbalanceStrategy != Constant ||
         Options.LoadbalanceAdjust != AdjustYes ||
         Options.SplitStrategy != MediumGrain ||
         Options.Alternate_FirstDirection != FirstDirRatio ||
         Options.SplitMethod != KLFM ||
         Options.Partitioner != PartMondriaan ||
         Options.Metric != MetricLambda ||
         Options.DiscardFreeNets != FreeNetYes ||
         Options.SquareMatrix_DistributeVectorsEqual != EqVecNo ||
         Options.SquareMatrix_DistributeVectorsEqual_AddDummies != DumYes ||
         Options.SymmetricMatrix_UseSingleEntry != SingleEntNo ||
         Options.SymmetricMatrix_SingleEntryType != ETypeLower ||
         Options.Coarsening_NrVertices != 200 ||
         Options.Coarsening_MaxCoarsenings != 128 ||
         Options.Coarsening_NrMatchArbitrary != 0 ||
         Options.Coarsening_MaxNrVtxInMatch != 2 ||
         Options.Coarsening_StopRatio != 0.05 ||
         Options.Coarsening_VtxMaxFractionOfWeight != 0.2 ||
         Options.Coarsening_MatchingStrategy != MatchATA ||
         Options.Coarsening_MatchingATAMatcher != MatchMatcherPGA ||
         Options.Coarsening_MatchingATAFinder != MatchFinderInproduct ||
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
         Options.Iterative_Refinement != IR_After_MG ||
         Options.VectorPartition_Step3 != VecRandom ||
         Options.VectorPartition_MaxNrLoops != 10 ||
         Options.VectorPartition_MaxNrGreedyImproves  != 10 ||
         Options.OutputFormat != OutputDMM ||
         Options.OutputMode != MultipleFiles ||
         Options.Seed != 99 ||
         Options.OrderPermutation != OrderNone ||
         Options.SymmetricPermute != FALSE) {
        printf("Error\n") ;
        exit(1);
    }
}

int main(int argc, char **argv) {
    struct opts Options;
    FILE *File;
    
    printf("Test SetDefaultOptions: ");
    
    SetDefaultOptions(&Options);
    /* SetOptionsFromFile(&Options, "test_SetDefaultOptions.defaults"); */
    CheckOptions(Options);
    
    File = fopen("test_SetDefaultOptions.check", "w");
    if (File != NULL) {
        ExportOptions(File, &Options);
        fclose(File);
        
        memset(&Options, 231, sizeof(struct opts));
        SetOptionsFromFile(&Options, "test_SetDefaultOptions.check");
        
        File = fopen("test_SetDefaultOptions.check2", "w");
        ExportOptions(File, &Options);
        fclose(File);
        
        CheckOptions(Options);
    }
    else {
        printf("Error\n");
        exit(1);
    }
 
    /* Check option values */

    printf("OK\n") ;
    exit(0);

} /* end main */
