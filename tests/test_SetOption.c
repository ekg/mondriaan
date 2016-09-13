#include "Options.h"

int main(int argc, char **argv) {

    struct opts Options;

    printf("Test SetOption: ");
    SetOption(&Options, "LoadbalanceStrategy", "constant");
    SetOption(&Options, "LoadbalanceAdjust", "no");
    SetOption(&Options, "SplitStrategy", "alternate");
    SetOption(&Options, "Alternate_FirstDirection", "row");
    SetOption(&Options, "SplitMethod", "simple");
    SetOption(&Options, "SquareMatrix_DistributeVectorsEqual", "no");
    SetOption(&Options, "SquareMatrix_DistributeVectorsEqual_AddDummies", "no");
    SetOption(&Options, "SymmetricMatrix_UseSingleEntry", "no");
    SetOption(&Options, "SymmetricMatrix_SingleEntryType", "lower");
    SetOption(&Options, "Coarsening_NrVertices", "1");
    SetOption(&Options, "Coarsening_MaxCoarsenings", "128");
    SetOption(&Options, "Coarsening_NrMatchArbitrary", "0");
    SetOption(&Options, "Coarsening_MaxNrVtxInMatch", "2");
    SetOption(&Options, "Coarsening_StopRatio", "0.25");
    SetOption(&Options, "Coarsening_VtxMaxFractionOfWeight", "0.25");
    SetOption(&Options, "Coarsening_MatchingStrategy", "random");
    SetOption(&Options, "Coarsening_MatchingATAMatcher", "greedy");
    SetOption(&Options, "Coarsening_MatchingATAFinder", "stairway");
    SetOption(&Options, "Coarsening_InprodMatchingOrder", "random");
    SetOption(&Options, "Coarsening_FineSwitchLevel", "2");
    SetOption(&Options, "Coarsening_NetScaling", "linear");
    SetOption(&Options, "Coarsening_InprodScaling", "jaccard");
    SetOption(&Options, "Coarsening_MatchIdenticalFirst", "no");
    SetOption(&Options, "KLFM_InitPart_NrRestarts", "10");
    SetOption(&Options, "KLFM_InitPart_MaxNrLoops", "10");
    SetOption(&Options, "KLFM_InitPart_MaxNrNoGainMoves", "10");
    SetOption(&Options, "KLFM_Refine_MaxNrLoops", "10");
    SetOption(&Options, "KLFM_Refine_MaxNrNoGainMoves", "10");
    SetOption(&Options, "VectorPartition_Step3", "increase");
    SetOption(&Options, "VectorPartition_MaxNrLoops", "10");
    SetOption(&Options, "VectorPartition_MaxNrGreedyImproves", "10");
    SetOption(&Options, "Seed", "random");
    
    /* Check option values */
    if ( Options.LoadbalanceStrategy != Constant ||
         Options.LoadbalanceAdjust != AdjustNo ||
         Options.SplitStrategy != Alternate ||
         Options.Alternate_FirstDirection != FirstDirRow ||
         Options.SplitMethod != Simple ||
         Options.SquareMatrix_DistributeVectorsEqual != EqVecNo ||
         Options.SquareMatrix_DistributeVectorsEqual_AddDummies != DumNo ||
         Options.SymmetricMatrix_UseSingleEntry != SingleEntNo ||
         Options.SymmetricMatrix_SingleEntryType != ETypeLower ||
         Options.Coarsening_NrVertices != 1 ||
         Options.Coarsening_MaxCoarsenings != 128 ||
         Options.Coarsening_NrMatchArbitrary != 0 ||
         Options.Coarsening_MaxNrVtxInMatch != 2 ||
         Options.Coarsening_StopRatio != 0.25 ||
         Options.Coarsening_VtxMaxFractionOfWeight != 0.25 ||
         Options.Coarsening_MatchingStrategy != MatchRandom ||
         Options.Coarsening_MatchingATAMatcher != MatchMatcherGreedy ||
         Options.Coarsening_MatchingATAFinder != MatchFinderStairway ||
         Options.Coarsening_InprodMatchingOrder != RandomOrder ||
         Options.Coarsening_FineSwitchLevel != 2 ||
         Options.Coarsening_NetScaling != NetSclLinear ||
         Options.Coarsening_InprodScaling != IpSclJaccard ||
         Options.Coarsening_MatchIdenticalFirst != MatchIdNo ||
         Options.KLFM_InitPart_NrRestarts != 10 ||
         Options.KLFM_InitPart_MaxNrLoops != 10 ||
         Options.KLFM_InitPart_MaxNrNoGainMoves != 10 ||
         Options.KLFM_Refine_MaxNrLoops != 10 ||
         Options.KLFM_Refine_MaxNrNoGainMoves != 10 ||
         Options.VectorPartition_Step3 != VecIncrease ||
         Options.VectorPartition_MaxNrLoops != 10 ||
         Options.VectorPartition_MaxNrGreedyImproves  != 10 ||
         Options.Seed != -1 ) {
        printf("Error\n");
        exit(1);
    }

    printf("OK\n");
    exit(0);

} /* end main */
