#include "Options.h"
#include "Sort.h"

int GetParameters(struct opts *pOptions, int argc, char **argv) {

    /* This function gets the parameters from the command line
       and enters these into the Options data structure.
       The command line has the syntax
           program matrix P eps -option0=value0 -option1=value1 etc.
       For example:
           mondriaan pi.mtx 4 0.1 -SplitStrategy=onedimrow
       
       It returns TRUE on success and FALSE on failure.
       This function does NOT initialise the entire pOptions structure,
       but only changes to parameters specified in the command line
       options!
    */
    int i;
    
    /* Check whether or not the given options are empty. */
    if (argc <= 0 || !argv || !pOptions) {
        fprintf(stderr, "GetParameters(): null arguments!\n");
        return FALSE;
    }
   
    /* Print Help if no option is given or a help option */
    if (argc == 1 ||  ! strcmp(argv[1], "-h") || !strcmp(argv[1], "-help") || ! strcmp(argv[1], "--help")) {
        PrintHelp(argc, argv);
        return FALSE;
    }
  
    /* Check the number of arguments */
    if (argc < 4) {
        fprintf(stderr, "GetParameters(): invalid number of arguments!\n");
        PrintHelp(argc, argv);
        return FALSE;
    }
  
    /* Read the required arguments */
    pOptions->matrix = argv[1];
    
    if (atoi(argv[2]) < 0) fprintf(stderr, "GetParameters(): invalid number of processors!\n");
    else pOptions->P = atoi(argv[2]);
    
    if (atof(argv[3]) < 0.0) fprintf(stderr, "GetParameters(): imbalance out of range!\n");
    else pOptions->eps = atof(argv[3]);
  
    /* Read the optional command line arguments, e.g.
           -SplitStrategy=onedimrow
       to overrule the corresponding defaults */
    i = 3;
    
    while (++i < argc) {
        if (argv[i][0] != '-') {
            fprintf(stderr, "GetParameters(): warning, invalid argument '%s'!\n", argv[i]);
        }
        else {
            const char *Option = argv[i] + 1, *Value = 0;

            /* Determine the form of the option: -option=value or -option value. */
            if (strchr(Option, '=')) {
                Value = strchr(Option, '=') + 1;
        *(strchr(Option, '=')) = 0;
            }
            else if (++i < argc) {
                Value = argv[i];
            }
            else {
                fprintf(stderr, "GetParameters(): warning, missing value for '%s'!\n", Option);
            }

            /* Set the read option. */
            SetOption(pOptions, Option, Value);
        }
    }
    
    return TRUE;
} /* end GetParameters */
  
int SetOptions(struct opts *pOptions, const char *Text) {

    /* This function sets the options to their default values
       as specified in a given string which may for example be
       extracted from a file on the harddisk.
       
       Each line of the defaults string contains an option and its desired value.
       Empty lines and lines starting with a # are ignored. 
       
       This function returns TRUE on success and FALSE on failure.
       */

    char Option[MAX_WORD_LENGTH], Value[MAX_WORD_LENGTH];
    
    /* First verify that we have valid input. */
    if (!pOptions || !Text) {
        fprintf(stderr, "SetOptions(): null arguments!\n");
        return FALSE;
    }

    /* Now loop through all the options and set them individually. */
    while (*Text) {
        /* Go to the first non-whitespace character. */
        while (*Text && isspace(*Text)) Text++;

        /* Discard this line if it is a comment, otherwise read it. */
        if (*Text && *Text != '#') {
            if (sscanf(Text, "%s %s", Option, Value) != 2) {
                fprintf(stderr, "SetOptions(): warning, unable to read the option from this line!\n");
            }
            else {
                SetOption(pOptions, Option, Value);
            }
        }

        /* Advance to the next line. */
        while (*Text && *Text != '\n') Text++;
    }

    return TRUE;
}

char* GetDefaultOptionText() {
    /* This function sets the default Mondriaan options. */
    return "# Mondriaan default settings: \n"
"LoadbalanceStrategy                            constant \n"
"LoadbalanceAdjust                              yes \n"
"SplitStrategy                                  mediumgrain \n"
"Alternate_FirstDirection                       ratio \n"
"SplitMethod                                    KLFM \n"
"Partitioner                                    mondriaan \n"
"Metric                                         lambda1 \n"
"Discard_Free_Nets                              yes \n"
"SquareMatrix_DistributeVectorsEqual            no \n"
"SquareMatrix_DistributeVectorsEqual_AddDummies yes \n"
"SymmetricMatrix_UseSingleEntry                 no \n"
"SymmetricMatrix_SingleEntryType                lower \n"
"Coarsening_NrVertices                          200 \n"
"Coarsening_MaxCoarsenings                      128 \n"
"Coarsening_NrMatchArbitrary                    0 \n"
"Coarsening_MaxNrVtxInMatch                     2 \n"
"Coarsening_StopRatio                           0.05 \n"
"Coarsening_VtxMaxFractionOfWeight              0.2 \n"
"Coarsening_MatchingStrategy                    ata \n"
"Coarsening_MatchingATAMatcher                  pga \n"
"Coarsening_MatchingATAFinder                   inproduct \n"
"Coarsening_InprodMatchingOrder                 decrwgt \n"
"Coarsening_FineSwitchLevel                     2 \n"
"Coarsening_NetScaling                          linear \n"
"Coarsening_InprodScaling                       min \n"
"Coarsening_MatchIdenticalFirst                 yes \n"
"KLFM_InitPart_NrRestarts                       8 \n"
"KLFM_InitPart_MaxNrLoops                       25 \n"
"KLFM_InitPart_MaxNrNoGainMoves                 200 \n"
"KLFM_Refine_MaxNrLoops                         25 \n"
"KLFM_Refine_MaxNrNoGainMoves                   200 \n"
"Iterative_Refinement                           aftermg \n"
"VectorPartition_Step3                          random \n"
"VectorPartition_MaxNrLoops                     10 \n"
"VectorPartition_MaxNrGreedyImproves            10 \n"
"OutputFormat                                   original \n"
"OutputMode                                     original \n"
"Seed                                           99 \n"
"Permute                                        none \n"
"EnforceSymmetricPermutation                    no ";
}

int SetDefaultOptions(struct opts *pOptions) {
    /* This function sets the default Mondriaan options. */
    const char* DefaultText = GetDefaultOptionText();
   
    pOptions->matrix = 0;
    
    return SetOptions(pOptions, DefaultText);
}

int ExportDefaultOptions(FILE *Out) {
    if (fprintf(Out, "%s", GetDefaultOptionText()) < 0)
        return FALSE;
    return TRUE;
}

/* The following function writes the Mondriaan options to a
filestream (i.e. stdout or an opened file on disk). */
int ExportOptions(FILE *Out, const struct opts *Opts) {
    if (Out == NULL || Opts == NULL) return FALSE;
    
    fprintf(Out, "# Mondriaan %s options:\n", MONDRIAANVERSION);
    
    fprintf(Out, "LoadbalanceStrategy ");
    if (Opts->LoadbalanceStrategy == Constant) fprintf(Out, "constant");
    else if (Opts->LoadbalanceStrategy == Increase) fprintf(Out, "increase");
    else if (Opts->LoadbalanceStrategy == Decrease) fprintf(Out, "decrease");
    else return FALSE;
    fprintf(Out, "\n");
    
    fprintf(Out, "LoadbalanceAdjust ");
    if (Opts->LoadbalanceAdjust == AdjustNo) fprintf(Out, "no");
    else if (Opts->LoadbalanceAdjust == AdjustYes) fprintf(Out, "yes");
    else return FALSE;
    fprintf(Out, "\n");
    
    fprintf(Out, "SplitStrategy ");
    if (Opts->SplitStrategy == Alternate) fprintf(Out, "alternate");
    else if (Opts->SplitStrategy == LocalBest) fprintf(Out, "localbest");
    else if (Opts->SplitStrategy == Hybrid) fprintf(Out, "hybrid");
    else if (Opts->SplitStrategy == LocalRatio) fprintf(Out, "localratio");
    else if (Opts->SplitStrategy == OneDimRow) fprintf(Out, "onedimrow");
    else if (Opts->SplitStrategy == OneDimCol) fprintf(Out, "onedimcol");
    else if (Opts->SplitStrategy == FineGrain) fprintf(Out, "finegrain");
    else if (Opts->SplitStrategy == SFineGrain) fprintf(Out, "symfinegrain");
    else if (Opts->SplitStrategy == MediumGrain) fprintf(Out, "mediumgrain");
    else return FALSE;
    fprintf(Out, "\n");

    fprintf(Out, "Alternate_FirstDirection ");
    if (Opts->Alternate_FirstDirection == FirstDirRow) fprintf(Out, "row");
    else if (Opts->Alternate_FirstDirection == FirstDirCol) fprintf(Out, "col");
    else if (Opts->Alternate_FirstDirection == FirstDirRatio) fprintf(Out, "ratio");
    else return FALSE;
    fprintf(Out, "\n");
    
    fprintf(Out, "SplitMethod ");
    if (Opts->SplitMethod == Simple) fprintf(Out, "simple");
    else if (Opts->SplitMethod == KLFM) fprintf(Out, "KLFM");
    else return FALSE;
    fprintf(Out, "\n");
    
    fprintf(Out, "Partitioner ");
    if (Opts->Partitioner == PartMondriaan) fprintf(Out, "mondriaan");
    else if (Opts->Partitioner == PartPaToH) fprintf(Out, "patoh");
    else if (Opts->Partitioner == FullPaToH) fprintf(Out, "fullpatoh");
    else return FALSE;
    fprintf(Out, "\n");
    
    fprintf(Out, "Metric ");
    if (Opts->Metric == MetricLambda) fprintf(Out, "lambda1");
    else if (Opts->Metric == MetricCut) fprintf(Out, "cutnet");
    else if (Opts->Metric == MetricLambdaLambdaMinusOne) fprintf(Out, "lambdalambda1");
    else return FALSE;
    fprintf(Out, "\n");
    
    fprintf(Out, "Discard_Free_Nets ");
    if (Opts->DiscardFreeNets == FreeNetYes) fprintf(Out, "yes");
    else if (Opts->DiscardFreeNets == FreeNetNo) fprintf(Out, "no");
    else return FALSE;
    fprintf(Out, "\n");
    
    fprintf(Out, "SquareMatrix_DistributeVectorsEqual ");
    if (Opts->SquareMatrix_DistributeVectorsEqual == EqVecNo) fprintf(Out, "no");
    else if (Opts->SquareMatrix_DistributeVectorsEqual == EqVecYes) fprintf(Out, "yes");
    else return FALSE;
    fprintf(Out, "\n");
    
    fprintf(Out, "SquareMatrix_DistributeVectorsEqual_AddDummies ");
    if (Opts->SquareMatrix_DistributeVectorsEqual_AddDummies == DumNo) fprintf(Out, "no");
    else if (Opts->SquareMatrix_DistributeVectorsEqual_AddDummies == DumYes) fprintf(Out, "yes");
    else return FALSE;
    fprintf(Out, "\n");
    
    fprintf(Out, "SymmetricMatrix_UseSingleEntry ");
    if (Opts->SymmetricMatrix_UseSingleEntry == SingleEntNo) fprintf(Out, "no");
    else if (Opts->SymmetricMatrix_UseSingleEntry == SingleEntYes) fprintf(Out, "yes");
    else return FALSE;
    fprintf(Out, "\n");
    
    fprintf(Out, "SymmetricMatrix_SingleEntryType ");
    if (Opts->SymmetricMatrix_SingleEntryType == ETypeLower) fprintf(Out, "lower");
    else if (Opts->SymmetricMatrix_SingleEntryType == ETypeRandom) fprintf(Out, "random");
    else return FALSE;
    fprintf(Out, "\n");
    
    fprintf(Out, "Coarsening_NrVertices %ld \n", Opts->Coarsening_NrVertices);
    fprintf(Out, "Coarsening_MaxCoarsenings %ld \n", Opts->Coarsening_MaxCoarsenings);
    fprintf(Out, "Coarsening_NrMatchArbitrary %ld \n", Opts->Coarsening_NrMatchArbitrary);
    fprintf(Out, "Coarsening_MaxNrVtxInMatch %ld \n", Opts->Coarsening_MaxNrVtxInMatch);
    fprintf(Out, "Coarsening_StopRatio %f \n", Opts->Coarsening_StopRatio);
    fprintf(Out, "Coarsening_VtxMaxFractionOfWeight %f \n", Opts->Coarsening_VtxMaxFractionOfWeight);
    
    fprintf(Out, "Coarsening_MatchingStrategy ");
    if (Opts->Coarsening_MatchingStrategy == MatchRandom) fprintf(Out, "random");
    else if (Opts->Coarsening_MatchingStrategy == MatchInprod) fprintf(Out, "inproduct");
    else if (Opts->Coarsening_MatchingStrategy == MatchATA) fprintf(Out, "ata");
    else return FALSE;
    fprintf(Out, "\n");
    
    fprintf(Out, "Coarsening_MatchingATAMatcher ");
    if (Opts->Coarsening_MatchingATAMatcher == MatchMatcherGreedy) fprintf(Out, "greedy");
    else if (Opts->Coarsening_MatchingATAMatcher == MatchMatcherPGA) fprintf(Out, "pga");
    else return FALSE;
    fprintf(Out, "\n");
    
    fprintf(Out, "Coarsening_MatchingATAFinder ");
    if (Opts->Coarsening_MatchingATAFinder == MatchFinderInproduct) fprintf(Out, "inproduct");
    else if (Opts->Coarsening_MatchingATAFinder == MatchFinderStairway) fprintf(Out, "stairway");
    else return FALSE;
    fprintf(Out, "\n");
    
    fprintf(Out, "Coarsening_InprodMatchingOrder ");
    if (Opts->Coarsening_InprodMatchingOrder == DecreasingWgt) fprintf(Out, "decrwgt");
    else if (Opts->Coarsening_InprodMatchingOrder == IncreasingWgt) fprintf(Out, "incrwgt");
    else if (Opts->Coarsening_InprodMatchingOrder == DecreasingDegree) fprintf(Out, "decrdegree");
    else if (Opts->Coarsening_InprodMatchingOrder == IncreasingDegree) fprintf(Out, "incrdegree");
    else if (Opts->Coarsening_InprodMatchingOrder == NaturalOrder) fprintf(Out, "natural");
    else if (Opts->Coarsening_InprodMatchingOrder == RandomOrder) fprintf(Out, "random");
    else return FALSE;
    fprintf(Out, "\n");
    
    fprintf(Out, "Coarsening_FineSwitchLevel %ld \n", Opts->Coarsening_FineSwitchLevel);
    
    fprintf(Out, "Coarsening_NetScaling ");
    if (Opts->Coarsening_NetScaling == NoNetScaling) fprintf(Out, "no");
    else if (Opts->Coarsening_NetScaling == NetSclLinear) fprintf(Out, "linear");
    else return FALSE;
    fprintf(Out, "\n");
    
    fprintf(Out, "Coarsening_InprodScaling ");
    if (Opts->Coarsening_InprodScaling == NoIpScaling) fprintf(Out, "no");
    else if (Opts->Coarsening_InprodScaling == IpSclCos) fprintf(Out, "cos");
    else if (Opts->Coarsening_InprodScaling == IpSclMin) fprintf(Out, "min");
    else if (Opts->Coarsening_InprodScaling == IpSclMax) fprintf(Out, "max");
    else if (Opts->Coarsening_InprodScaling == IpSclJaccard) fprintf(Out, "jaccard");
    else return FALSE;
    fprintf(Out, "\n");
    
    fprintf(Out, "Coarsening_MatchIdenticalFirst ");
    if (Opts->Coarsening_MatchIdenticalFirst == MatchIdNo) fprintf(Out, "no");
    else if (Opts->Coarsening_MatchIdenticalFirst == MatchIdYes) fprintf(Out, "yes");
    else return FALSE;
    fprintf(Out, "\n");
    
    fprintf(Out, "KLFM_InitPart_NrRestarts %ld \n", Opts->KLFM_InitPart_NrRestarts);
    fprintf(Out, "KLFM_InitPart_MaxNrLoops %ld \n", Opts->KLFM_InitPart_MaxNrLoops);
    fprintf(Out, "KLFM_InitPart_MaxNrNoGainMoves %ld \n", Opts->KLFM_InitPart_MaxNrNoGainMoves);
    fprintf(Out, "KLFM_Refine_MaxNrLoops %ld \n", Opts->KLFM_Refine_MaxNrLoops);
    fprintf(Out, "KLFM_Refine_MaxNrNoGainMoves %ld \n", Opts->KLFM_Refine_MaxNrNoGainMoves);
    
    fprintf(Out, "Iterative_Refinement " );
    if (Opts->Iterative_Refinement == IR_Never)
        fprintf(Out, "never" );
    else if (Opts->Iterative_Refinement == IR_After_MG)
        fprintf(Out, "aftermg" );
    else if (Opts->Iterative_Refinement == IR_Always)
        fprintf(Out, "always" );
    else return FALSE;
    fprintf(Out, "\n");
    
    fprintf(Out, "VectorPartition_Step3 ");
    if (Opts->VectorPartition_Step3 == VecIncrease) fprintf(Out, "increase");
    else if (Opts->VectorPartition_Step3 == VecDecrease) fprintf(Out, "decrease");
    else if (Opts->VectorPartition_Step3 == VecRandom) fprintf(Out, "random");
    else return FALSE;
    fprintf(Out, "\n");
    
    fprintf(Out, "VectorPartition_MaxNrLoops %ld \n", Opts->VectorPartition_MaxNrLoops);
    fprintf(Out, "VectorPartition_MaxNrGreedyImproves %ld \n", Opts->VectorPartition_MaxNrGreedyImproves);
    
    fprintf(Out, "OutputFormat " );
    if (Opts->OutputFormat == OutputDMM)
        fprintf(Out, "original" );
    else if (Opts->OutputFormat == OutputEMM)
        fprintf(Out, "emm" );
    else return FALSE;
    fprintf(Out, "\n");

    fprintf(Out, "OutputMode ");
    if (Opts->OutputMode == OneFile)
        fprintf(Out, "onefile");
    else if (Opts->OutputMode == MultipleFiles)
        fprintf(Out, "original");
    else if (Opts->OutputMode == DIMACS)
        fprintf(Out, "DIMACS");
    else return FALSE;
    fprintf(Out, "\n");
    
    fprintf(Out, "Seed %ld \n", Opts->Seed);
    
    fprintf(Out, "Permute ");
    if (Opts->OrderPermutation == OrderNone) fprintf(Out, "none");
    else if (Opts->OrderPermutation == OrderPrefix) fprintf(Out, "reverseBBD");
    else if (Opts->OrderPermutation == OrderInfix) fprintf(Out, "SBD");
    else if (Opts->OrderPermutation == OrderPostfix) fprintf(Out, "BBD");
    else return FALSE;
    fprintf(Out, "\n");

    fprintf(Out, "EnforceSymmetricPermutation ");
    if (Opts->SymmetricPermute) fprintf(Out, "yes");
    else fprintf(Out, "no");
    fprintf(Out, "\n");
    
    return TRUE;
}

int ExportOptionsToLaTeX(FILE *Out, const struct opts *Opts) {
    /* This function prints all set options to a LaTeX table. */
    if (Out == NULL || Opts == NULL) return FALSE;
    
    fprintf(Out, "\\begin{tabular}{l|l}\nOption & Value \\\\\n\\hline\n");
    
    fprintf(Out, "P & %d \\\\\n", Opts->P);
    fprintf(Out, "eps & %e \\\\\n", Opts->eps);
    
    fprintf(Out, "LoadbalanceStrategy & ");
    if (Opts->LoadbalanceStrategy == Constant) fprintf(Out, "constant");
    else if (Opts->LoadbalanceStrategy == Increase) fprintf(Out, "increase");
    else if (Opts->LoadbalanceStrategy == Decrease) fprintf(Out, "decrease");
    else fprintf(Out, "?");
    fprintf(Out, " \\\\\n");
    
    fprintf(Out, "LoadbalanceAdjust & ");
    if (Opts->LoadbalanceAdjust == AdjustNo) fprintf(Out, "no");
    else if (Opts->LoadbalanceAdjust == AdjustYes) fprintf(Out, "yes");
    else fprintf(Out, "?");
    fprintf(Out, " \\\\\n");
    
    fprintf(Out, "SplitStrategy & ");
    if (Opts->SplitStrategy == Alternate) fprintf(Out, "alternate");
    else if (Opts->SplitStrategy == LocalBest) fprintf(Out, "localbest");
    else if (Opts->SplitStrategy == Hybrid) fprintf(Out, "hybrid");
    else if (Opts->SplitStrategy == LocalRatio) fprintf(Out, "localratio");
    else if (Opts->SplitStrategy == OneDimRow) fprintf(Out, "onedimrow");
    else if (Opts->SplitStrategy == OneDimCol) fprintf(Out, "onedimcol");
    else if (Opts->SplitStrategy == FineGrain) fprintf(Out, "finegrain");
    else if (Opts->SplitStrategy == SFineGrain) fprintf(Out, "symfinegrain");
    else if (Opts->SplitStrategy == MediumGrain) fprintf(Out, "mediumgrain");
    else fprintf(Out, "?");
    fprintf(Out, " \\\\\n");
    
    fprintf(Out, "Alternate-FirstDirection & ");
    if (Opts->Alternate_FirstDirection == FirstDirRow) fprintf(Out, "row");
    else if (Opts->Alternate_FirstDirection == FirstDirCol) fprintf(Out, "col");
    else if (Opts->Alternate_FirstDirection == FirstDirRatio) fprintf(Out, "ratio");
    else fprintf(Out, "?");
    fprintf(Out, " \\\\\n");
    
    fprintf(Out, "SplitMethod & ");
    if (Opts->SplitMethod == Simple) fprintf(Out, "simple");
    else if (Opts->SplitMethod == KLFM) fprintf(Out, "KLFM");
    else fprintf(Out, "?");
    fprintf(Out, " \\\\\n");
    
    fprintf(Out, "Partitioner & ");
    if (Opts->Partitioner == PartMondriaan) fprintf(Out, "mondriaan");
    else if (Opts->Partitioner == PartPaToH) fprintf(Out, "patoh");
    else if (Opts->Partitioner == FullPaToH) fprintf(Out, "fullpatoh");
    else fprintf(Out, "?");
    fprintf(Out, " \\\\\n");
    
    fprintf(Out, "Metric & ");
    if (Opts->Metric == MetricLambda) fprintf(Out, "lambda1");
    else if (Opts->Metric == MetricCut) fprintf(Out, "cutnet");
    else if (Opts->Metric == MetricLambdaLambdaMinusOne) fprintf(Out, "lambdalambda1");
    else fprintf(Out, "?");
    fprintf(Out, " \\\\\n");
    
    fprintf(Out, "Discard-Free-Nets & ");
    if (Opts->DiscardFreeNets == FreeNetYes) fprintf(Out, "yes");
    else if (Opts->DiscardFreeNets == FreeNetNo) fprintf(Out, "no");
    else fprintf(Out, "?");
    fprintf(Out, " \\\\\n");
    
    fprintf(Out, "DistributeVectorsEqual & ");
    if (Opts->SquareMatrix_DistributeVectorsEqual == EqVecNo) fprintf(Out, "no");
    else if (Opts->SquareMatrix_DistributeVectorsEqual == EqVecYes) fprintf(Out, "yes");
    else fprintf(Out, "?");
    fprintf(Out, " \\\\\n");
    
    fprintf(Out, "DistributeVectorsEqual-AddDummies & ");
    if (Opts->SquareMatrix_DistributeVectorsEqual_AddDummies == DumNo) fprintf(Out, "no");
    else if (Opts->SquareMatrix_DistributeVectorsEqual_AddDummies == DumYes) fprintf(Out, "yes");
    else fprintf(Out, "?");
    fprintf(Out, " \\\\\n");
    
    fprintf(Out, "SymmetricMatrix-UseSingleEntry & ");
    if (Opts->SymmetricMatrix_UseSingleEntry == SingleEntNo) fprintf(Out, "no");
    else if (Opts->SymmetricMatrix_UseSingleEntry == SingleEntYes) fprintf(Out, "yes");
    else fprintf(Out, "?");
    fprintf(Out, " \\\\\n");
    
    fprintf(Out, "SymmetricMatrix-SingleEntryType & ");
    if (Opts->SymmetricMatrix_SingleEntryType == ETypeLower) fprintf(Out, "lower");
    else if (Opts->SymmetricMatrix_SingleEntryType == ETypeRandom) fprintf(Out, "random");
    else fprintf(Out, "?");
    fprintf(Out, " \\\\\n");
    
    fprintf(Out, "Coarsening-NrVertices & %ld \\\\\n", Opts->Coarsening_NrVertices);
    fprintf(Out, "Coarsening-MaxCoarsenings & %ld \\\\\n", Opts->Coarsening_MaxCoarsenings);
    fprintf(Out, "Coarsening-MaxNrMatchArbitrary & %ld \\\\\n", Opts->Coarsening_NrMatchArbitrary);
    fprintf(Out, "Coarsening-MaxNrVtxInMatch & %ld \\\\\n", Opts->Coarsening_MaxNrVtxInMatch);
    fprintf(Out, "Coarsening-StopRatio & %f \\\\\n", Opts->Coarsening_StopRatio);
    fprintf(Out, "Coarsening-VtxMaxFractionOfWeight & %f \\\\\n", Opts->Coarsening_VtxMaxFractionOfWeight);
    
    fprintf(Out, "Coarsening-MatchingStrategy & ");
    if (Opts->Coarsening_MatchingStrategy == MatchRandom) fprintf(Out, "random");
    else if (Opts->Coarsening_MatchingStrategy == MatchInprod) fprintf(Out, "inproduct");
    else if (Opts->Coarsening_MatchingStrategy == MatchATA) fprintf(Out, "ata");
    else fprintf(Out, "?");
    fprintf(Out, " \\\\\n");
    
    fprintf(Out, "Coarsening-MatchingATAMatcher & ");
    if (Opts->Coarsening_MatchingATAMatcher == MatchMatcherGreedy) fprintf(Out, "greedy");
    else if (Opts->Coarsening_MatchingATAMatcher == MatchMatcherPGA) fprintf(Out, "pga");
    else fprintf(Out, "?");
    fprintf(Out, " \\\\\n");
    
    fprintf(Out, "Coarsening-MatchingATAFinder & ");
    if (Opts->Coarsening_MatchingATAFinder == MatchFinderInproduct) fprintf(Out, "inproduct");
    else if (Opts->Coarsening_MatchingATAFinder == MatchFinderStairway) fprintf(Out, "stairway");
    else fprintf(Out, "?");
    fprintf(Out, " \\\\\n");
    
    fprintf(Out, "Coarsening-InprodMatchingOrder & ");
    if (Opts->Coarsening_InprodMatchingOrder == DecreasingWgt) fprintf(Out, "decrwgt");
    else if (Opts->Coarsening_InprodMatchingOrder == IncreasingWgt) fprintf(Out, "incrwgt");
    else if (Opts->Coarsening_InprodMatchingOrder == DecreasingDegree) fprintf(Out, "decrdegree");
    else if (Opts->Coarsening_InprodMatchingOrder == IncreasingDegree) fprintf(Out, "incrdegree");
    else if (Opts->Coarsening_InprodMatchingOrder == NaturalOrder) fprintf(Out, "natural");
    else if (Opts->Coarsening_InprodMatchingOrder == RandomOrder) fprintf(Out, "random");
    else fprintf(Out, "?");
    fprintf(Out, " \\\\\n");
    
    fprintf(Out, "Coarsening-FineSwitchLevel & %ld \\\\\n", Opts->Coarsening_FineSwitchLevel);
    
    fprintf(Out, "Coarsening-NetScaling & ");
    if (Opts->Coarsening_NetScaling == NoNetScaling) fprintf(Out, "no");
    else if (Opts->Coarsening_NetScaling == NetSclLinear) fprintf(Out, "linear");
    else fprintf(Out, "?");
    fprintf(Out, " \\\\\n");
    
    fprintf(Out, "Coarsening-InprodScaling & ");
    if (Opts->Coarsening_InprodScaling == NoIpScaling) fprintf(Out, "no");
    else if (Opts->Coarsening_InprodScaling == IpSclCos) fprintf(Out, "cos");
    else if (Opts->Coarsening_InprodScaling == IpSclMin) fprintf(Out, "min");
    else if (Opts->Coarsening_InprodScaling == IpSclMax) fprintf(Out, "max");
    else if (Opts->Coarsening_InprodScaling == IpSclJaccard) fprintf(Out, "jaccard");
    else fprintf(Out, "?");
    fprintf(Out, " \\\\\n");
    
    fprintf(Out, "Coarsening-MatchIdenticalFirst & ");
    if (Opts->Coarsening_MatchIdenticalFirst == MatchIdNo) fprintf(Out, "no");
    else if (Opts->Coarsening_MatchIdenticalFirst == MatchIdYes) fprintf(Out, "yes");
    else fprintf(Out, "?");
    fprintf(Out, " \\\\\n");
    
    fprintf(Out, "KLFM-InitPart-NrRestarts & %ld \\\\\n", Opts->KLFM_InitPart_NrRestarts);
    fprintf(Out, "KLFM-InitPart-MaxNrLoops & %ld \\\\\n", Opts->KLFM_InitPart_MaxNrLoops);
    fprintf(Out, "KLFM-InitPart-MaxNrNoGainMoves & %ld \\\\\n", Opts->KLFM_InitPart_MaxNrNoGainMoves);
    fprintf(Out, "KLFM-Refine-MaxNrLoops & %ld \\\\\n", Opts->KLFM_Refine_MaxNrLoops);
    fprintf(Out, "KLFM-Refine-MaxNrNoGainMoves & %ld \\\\\n", Opts->KLFM_Refine_MaxNrNoGainMoves);
    
    fprintf(Out, "VectorPartition-Step3 & ");
    if (Opts->VectorPartition_Step3 == VecIncrease) fprintf(Out, "increase");
    else if (Opts->VectorPartition_Step3 == VecDecrease) fprintf(Out, "decrease");
    else if (Opts->VectorPartition_Step3 == VecRandom) fprintf(Out, "random");
    else fprintf(Out, "?");
    fprintf(Out, " \\\\\n");
    
    fprintf(Out, "VectorPartition-MaxNrLoops & %ld \\\\\n", Opts->VectorPartition_MaxNrLoops);
    fprintf(Out, "VectorPartition-MaxNrGreedyImproves & %ld \\\\\n", Opts->VectorPartition_MaxNrGreedyImproves);
    
    fprintf(Out, "Seed & %ld \\\\\n", Opts->Seed);
    
    fprintf(Out, "Permute & ");
    if (Opts->OrderPermutation == OrderNone) fprintf(Out, "none");
    else if (Opts->OrderPermutation == OrderPrefix) fprintf(Out, "reverseBBD");
    else if (Opts->OrderPermutation == OrderInfix) fprintf(Out, "SBD");
    else if (Opts->OrderPermutation == OrderPostfix) fprintf(Out, "BBD");
    else fprintf(Out, "?");
    fprintf(Out, "\\\\\n");

    fprintf(Out, "EnforceSymmetricPermutation & ");
    if (Opts->SymmetricPermute) fprintf(Out, "yes");
    else fprintf(Out, "no");
    
    fprintf(Out, "\n\\\\\\hline\n");
    
    fprintf(Out, "\\end{tabular}\n");
    
    return TRUE;
}

/* The following function reads an options file from disk to memory
and passes it to SetOptions().

The function returns TRUE on success and FALSE on failure.*/
int SetOptionsFromFile(struct opts *pOptions, const char *File) {
    FILE *pFile;
    char *Text;
    size_t Size;
    int Failure = FALSE;
    
    if (!pOptions || !File) {
        fprintf(stderr, "SetOptionsFromFile(): null arguments!\n");
        return FALSE;
    }
        
    if ((pFile = fopen(File, "rb"))) {
        fseek(pFile, 0, SEEK_END);
        Size = ftell(pFile);
        fseek(pFile, 0, SEEK_SET);

        if ((Text = (char *)malloc(Size + 1))) {
            if (fread(Text, Size, 1, pFile) == 1) {
                Text[Size] = 0;
                SetOptions(pOptions, Text);
            }
            else {
                fprintf(stderr, "SetOptionsFromFile(): could not read data from options file!\n");
                Failure = TRUE;
            }

            free(Text);
        }
        else {
            fprintf(stderr, "SetOptionsFromFile(): could not allocate text buffer for options file!\n");
            Failure = TRUE;
        }

        fclose(pFile);
    }
    else {
        fprintf(stderr, "SetOptionsFromFile(): could not open options file '%s'!\n", File);
        Failure = TRUE;
    }

    return (!Failure);
}

int ApplyOptions(const struct opts *pOptions) {
    /* This function checks the option values in the given options structure and sets the random seed.
       Returns TRUE on success and FALSE on failure.
    */
    if (!pOptions) {
        fprintf(stderr, "ApplyOptions(): null arguments!\n");
        return FALSE;
    }
    
    /* Check whether the numerical option values are within the proper range */
    if (pOptions->Coarsening_NrVertices < 1) {
        fprintf(stderr, "ApplyOptions(): Coarsening_NrVertices out of range!\n");
        return FALSE;
    }
    
    if (pOptions->Coarsening_MaxCoarsenings < 1) {
        fprintf(stderr, "ApplyOptions(): Coarsening_MaxCoarsenings out of range!\n");
        return FALSE;
    }
  
    if (pOptions->Coarsening_MaxNrVtxInMatch < 2) {
        fprintf(stderr, "ApplyOptions(): Coarsening_MaxNrVtxInMatch out of range!\n");
        return FALSE;
    }
  
    if (pOptions->Coarsening_StopRatio < 0 || 
        pOptions->Coarsening_StopRatio > 1) {
        fprintf(stderr, "ApplyOptions(): Coarsening_StopRatio out of range!\n");
        return FALSE;
    }
    
    if (pOptions->Coarsening_VtxMaxFractionOfWeight <= 0 || 
        pOptions->Coarsening_VtxMaxFractionOfWeight > 1) {
        fprintf(stderr, "ApplyOptions(): Coarsening_VtxMaxFractionOfWeight out of range!\n");
        return FALSE;
    }
  
    if (pOptions->KLFM_InitPart_NrRestarts < 1) {
        fprintf(stderr, "ApplyOptions(): KLFM_InitPart_NrRestarts out of range!\n");
        return FALSE;
    }
    
    if (pOptions->KLFM_InitPart_MaxNrLoops < 1) {
        fprintf(stderr, "ApplyOptions(): KLFM_InitPart_MaxNrLoops out of range!\n");
        return FALSE;
    }
    
    if (pOptions->KLFM_InitPart_MaxNrNoGainMoves < 0) {
        fprintf(stderr, "ApplyOptions(): KLFM_InitPart_MaxNrNoGainMoves out of range!\n");
        return FALSE;
    }
    
    if (pOptions->KLFM_Refine_MaxNrLoops < 1) {
        fprintf(stderr, "ApplyOptions(): KLFM_Refine_MaxNrLoops out of range!\n");
        return FALSE;
    }
  
    if (pOptions->KLFM_Refine_MaxNrNoGainMoves < 0) {
        fprintf(stderr, "ApplyOptions(): KLFM_Refine_MaxNrNoGainMoves out of range!\n");
        return FALSE;
    }
  
    if (pOptions->VectorPartition_MaxNrGreedyImproves < 0) {
        fprintf(stderr, "ApplyOptions(): VectorPartition_MaxNrGreedyImproves out of range!\n");
        return FALSE;
    }
  
    if (pOptions->VectorPartition_MaxNrLoops < 1) {
        fprintf(stderr, "ApplyOptions(): VectorPartition_MaxNrLoops out of range!\n");
        return FALSE;
    }
  
    /* Set the random number seed */
    if (pOptions->Seed < -1) {
        fprintf(stderr, "ApplyOptions(): Random number seed out of range!\n");
        return FALSE;
    }
    
    SetRandomSeed(pOptions->Seed);
    
    if (pOptions->Metric == MetricCut && pOptions->DiscardFreeNets == FreeNetNo) {
        fprintf(stderr, "ApplyOptions(): For the cut-net metric it is mandatory that free nets are discarded!\n");
        return FALSE;
    }
    
    if (pOptions->Metric == MetricLambdaLambdaMinusOne && pOptions->Partitioner != PartPaToH) {
        fprintf(stderr, "ApplyOptions(): PaToH is required for the use of the lambda*(lambda - 1) metric!\n");
        return FALSE;
    }
    
    if (pOptions->Coarsening_MatchingStrategy == MatchATA && pOptions->Coarsening_MaxNrVtxInMatch != 2) {
        fprintf(stderr, "ApplyOptions(): Hybrid matching is only supported for matching groups of two vertices!\n");
        return FALSE;
    }

#ifndef USE_PATOH
    if (pOptions->Partitioner == PartPaToH || pOptions->Partitioner == FullPaToH) {
        fprintf(stderr, "ApplyOptions(): The use of PaToH is requested, but this version of Mondriaan is compiled without PaToH support!\n");
        return FALSE;
    }
#endif
    
    return TRUE;
}

int SetOption(struct opts *pOptions, const char *option, const char *value) {

    /* This function sets the given option to its desired value
       in the Options data structure. The function does this
       by checking all the possible options to see whether
       the given option matches, and then checking all the
       possible values. 
       It returns TRUE when the given option has succesfully
       been set to the desired value and FALSE otherwise.*/
    if (!pOptions || !option || !value) return FALSE;
    if (!strcmp(option, "")) return FALSE;
    
    if (!strcmp(option, "LoadbalanceStrategy")) {
        if (!strcmp(value, "constant"))
            pOptions->LoadbalanceStrategy = Constant;
        else if (!strcmp(value, "increase"))
            pOptions->LoadbalanceStrategy = Increase;
        else if (!strcmp(value, "decrease"))
            pOptions->LoadbalanceStrategy = Decrease;
        else {
            fprintf(stderr, "SetOption(): unknown %s '%s'!\n", option, value);
            return FALSE;
        }
    } else if (!strcmp(option, "LoadbalanceAdjust")) {
        if (!strcmp(value, "no") || !strcmp(value, "0"))
          pOptions->LoadbalanceAdjust = AdjustNo;
        else if (!strcmp(value, "yes") || !strcmp(value, "1"))
          pOptions->LoadbalanceAdjust = AdjustYes;
        else {
            fprintf(stderr, "SetOption(): unknown %s '%s'!\n", option, value);
            return FALSE;
        }
    } else if (!strcmp(option, "SplitStrategy")) {
        if (!strcmp(value, "alternate"))
            pOptions->SplitStrategy = Alternate;
        else if (!strcmp(value, "localbest"))
            pOptions->SplitStrategy = LocalBest;
        else if (!strcmp(value, "hybrid"))
            pOptions->SplitStrategy = Hybrid;
        else if (!strcmp(value, "localratio"))
            pOptions->SplitStrategy = LocalRatio;
        else if (!strcmp(value, "onedimrow"))
            pOptions->SplitStrategy = OneDimRow;
        else if (!strcmp(value, "onedimcol"))
            pOptions->SplitStrategy = OneDimCol;
        else if (!strcmp(value, "finegrain"))
            pOptions->SplitStrategy = FineGrain;
        else if (!strcmp(value, "symfinegrain"))
            pOptions->SplitStrategy = SFineGrain;
        else if (!strcmp(value, "mediumgrain"))
            pOptions->SplitStrategy = MediumGrain;
        else {
            fprintf(stderr, "SetOption(): unknown %s '%s'!\n", option, value);
            return FALSE;
        }
    } else if (!strcmp(option, "Alternate_FirstDirection")) {
        if (!strcmp(value, "row"))
            pOptions->Alternate_FirstDirection = FirstDirRow;
        else if (!strcmp(value, "col"))
            pOptions->Alternate_FirstDirection = FirstDirCol;
        else if (!strcmp(value, "ratio"))
            pOptions->Alternate_FirstDirection = FirstDirRatio;
        else {
            fprintf(stderr, "SetOption(): unknown %s '%s'!\n", option, value);
            return FALSE;
        }
    } else if (!strcmp(option, "SplitMethod")) {
        if (!strcmp(value, "simple"))
          pOptions->SplitMethod = Simple;
        else if (!strcmp(value, "KLFM") || !strcmp(value, "klfm"))
          pOptions->SplitMethod = KLFM;
        else {
            fprintf(stderr, "SetOption(): unknown %s '%s'!\n", option, value);
            return FALSE;
        }
    } else if (!strcmp(option, "Partitioner")) {
        if (!strcmp(value, "mondriaan"))
          pOptions->Partitioner = PartMondriaan;
        else if (!strcmp(value, "patoh"))
          pOptions->Partitioner = PartPaToH;
        else if (!strcmp(value, "fullpatoh"))
          pOptions->Partitioner = FullPaToH;
        else {
            fprintf(stderr, "SetOption(): unknown %s '%s'!\n", option, value);
            return FALSE;
        }
    } else if (!strcmp(option, "Metric")) {
        if (!strcmp(value, "lambda1"))
          pOptions->Metric = MetricLambda;
        else if (!strcmp(value, "cutnet"))
          pOptions->Metric = MetricCut;
        else if (!strcmp(value, "lambdalambda1"))
          pOptions->Metric = MetricLambdaLambdaMinusOne;
        else {
            fprintf(stderr, "SetOption(): unknown %s '%s'!\n", option, value);
            return FALSE;
        }
    } else if (!strcmp(option, "Discard_Free_Nets")) {
        if (!strcmp(value, "yes")) {
            pOptions->DiscardFreeNets = FreeNetYes;
        } else if (!strcmp(value, "no")) {
            pOptions->DiscardFreeNets = FreeNetNo;
        } else {
            fprintf(stderr, "SetOptions(): unknown %s '%s'!\n", option, value);
            return FALSE;
        }
    } else if (!strcmp(option, "SquareMatrix_DistributeVectorsEqual")) {
        if (!strcmp(value, "no") || !strcmp(value, "0"))
            pOptions->SquareMatrix_DistributeVectorsEqual = EqVecNo;
        else if (!strcmp(value, "yes") || ! strcmp(value, "1"))
            pOptions->SquareMatrix_DistributeVectorsEqual = EqVecYes;
        else {
            fprintf(stderr, "SetOption(): unknown %s '%s'!\n", option, value);
            return FALSE;
        }
    } else if (!strcmp(option, "SquareMatrix_DistributeVectorsEqual_AddDummies")) {
        if (!strcmp(value, "no") || !strcmp(value, "0"))
            pOptions->SquareMatrix_DistributeVectorsEqual_AddDummies = DumNo;
        else if (!strcmp(value, "yes") || !strcmp(value, "1"))
            pOptions->SquareMatrix_DistributeVectorsEqual_AddDummies = DumYes;
        else {
            fprintf(stderr, "SetOption(): unknown %s '%s'!\n", option, value);
            return FALSE;
        }
    } else if (!strcmp(option, "SymmetricMatrix_UseSingleEntry")) {
        if (!strcmp(value, "no") || ! strcmp(value, "0"))
            pOptions->SymmetricMatrix_UseSingleEntry = SingleEntNo;
        else if (!strcmp(value, "yes") || ! strcmp(value, "1"))
            pOptions->SymmetricMatrix_UseSingleEntry = SingleEntYes;
        else {
            fprintf(stderr, "SetOption(): unknown %s '%s'!\n", option, value);
            return FALSE;
        }
    } else if (!strcmp(option, "SymmetricMatrix_SingleEntryType")) {
        if (!strcmp(value, "lower"))
            pOptions->SymmetricMatrix_SingleEntryType = ETypeLower;
        else if (!strcmp(value, "random"))
            pOptions->SymmetricMatrix_SingleEntryType = ETypeRandom;
        else {
            fprintf(stderr, "SetOption(): unknown %s '%s'!\n", option, value);
            return FALSE;
        }
    } else if (!strcmp(option, "Coarsening_NrVertices")) {
        pOptions->Coarsening_NrVertices = atol(value);
    } else if (!strcmp(option, "Coarsening_MaxCoarsenings")) {
        pOptions->Coarsening_MaxCoarsenings = atol(value);
    } else if (!strcmp(option, "Coarsening_NrMatchArbitrary")) {
        pOptions->Coarsening_NrMatchArbitrary = atol(value);
    } else if (!strcmp(option, "Coarsening_MaxNrVtxInMatch")) {
        pOptions->Coarsening_MaxNrVtxInMatch = atol(value);
    } else if (!strcmp(option, "Coarsening_StopRatio")) {
        pOptions->Coarsening_StopRatio = atof(value);
    } else if (!strcmp(option, "Coarsening_VtxMaxFractionOfWeight")) {
        pOptions->Coarsening_VtxMaxFractionOfWeight = atof(value);
    } else if (!strcmp(option, "Coarsening_MatchingStrategy")) {
        if (!strcmp(value, "random"))
            pOptions->Coarsening_MatchingStrategy = MatchRandom;
        else if (!strcmp(value, "inproduct"))
            pOptions->Coarsening_MatchingStrategy = MatchInprod;
        else if (!strcmp(value, "ata"))
            pOptions->Coarsening_MatchingStrategy = MatchATA;
        else {
            fprintf(stderr, "SetOption(): unknown %s '%s'!\n", option, value);
            return FALSE;
        }
    } else if (!strcmp(option, "Coarsening_MatchingATAMatcher")) {
        if (!strcmp(value, "greedy"))
            pOptions->Coarsening_MatchingATAMatcher = MatchMatcherGreedy;
        else if (!strcmp(value, "pga"))
            pOptions->Coarsening_MatchingATAMatcher = MatchMatcherPGA;
        else {
            fprintf(stderr, "SetOption(): unknown %s '%s'!\n", option, value);
            return FALSE;
        }
    } else if (!strcmp(option, "Coarsening_MatchingATAFinder")) {
        if (!strcmp(value, "inproduct"))
            pOptions->Coarsening_MatchingATAFinder = MatchFinderInproduct;
        else if (!strcmp(value, "stairway"))
            pOptions->Coarsening_MatchingATAFinder = MatchFinderStairway;
        else {
            fprintf(stderr, "SetOption(): unknown %s '%s'!\n", option, value);
            return FALSE;
        }
    } else if (!strcmp(option, "Coarsening_FineSwitchLevel")) {
        pOptions->Coarsening_FineSwitchLevel = atol(value);
    } else if (!strcmp(option, "Coarsening_InprodMatchingOrder")) {
        if (!strcmp(value, "decrwgt"))
            pOptions->Coarsening_InprodMatchingOrder = DecreasingWgt;
        else if (!strcmp(value, "incrwgt"))
            pOptions->Coarsening_InprodMatchingOrder = IncreasingWgt;
        else if (!strcmp(value, "decrdeg"))
            pOptions->Coarsening_InprodMatchingOrder = DecreasingDegree;
        else if (!strcmp(value, "incrdeg"))
            pOptions->Coarsening_InprodMatchingOrder = IncreasingDegree;
        else if (!strcmp(value, "natural"))
            pOptions->Coarsening_InprodMatchingOrder = NaturalOrder;
        else if (!strcmp(value, "random"))
            pOptions->Coarsening_InprodMatchingOrder = RandomOrder;
        else {
            fprintf(stderr, "SetOption(): unknown %s '%s'!\n", option, value);
            return FALSE;
        }
    } else if (!strcmp(option, "Coarsening_NetScaling")) {
        if (!strcmp(value, "nonetscaling") || ! strcmp(value, "no"))
            pOptions->Coarsening_NetScaling = NoNetScaling;
        else if (!strcmp(value, "linear"))
            pOptions->Coarsening_NetScaling = NetSclLinear;
        else {
            fprintf(stderr, "SetOption(): unknown %s '%s'!\n", option, value);
            return FALSE;
        }
    } else if (!strcmp(option, "Coarsening_InprodScaling")) {
        if (!strcmp(value, "noipscaling") || ! strcmp(value, "no"))
            pOptions->Coarsening_InprodScaling = NoIpScaling;
        else if (!strcmp(value, "cos"))
            pOptions->Coarsening_InprodScaling = IpSclCos;
        else if (!strcmp(value, "min"))
            pOptions->Coarsening_InprodScaling = IpSclMin;
        else if (!strcmp(value, "max"))
            pOptions->Coarsening_InprodScaling = IpSclMax;
        else if (!strcmp(value, "jaccard"))
            pOptions->Coarsening_InprodScaling = IpSclJaccard;
        else {
            fprintf(stderr, "SetOption(): unknown %s '%s'!\n", option, value);
            return FALSE;
        }
    } else if (!strcmp(option, "Coarsening_MatchIdenticalFirst")) {
        if (!strcmp(value, "no") || ! strcmp(value, "0"))
            pOptions->Coarsening_MatchIdenticalFirst = MatchIdNo;
        else if (!strcmp(value, "yes") || ! strcmp(value, "1"))
            pOptions->Coarsening_MatchIdenticalFirst = MatchIdYes;
        else {
            fprintf(stderr, "SetOption(): unknown %s '%s'!\n", option, value);
            return FALSE;
        }
    } else if (!strcmp(option, "KLFM_InitPart_NrRestarts")) {
        pOptions->KLFM_InitPart_NrRestarts = atol(value);
    } else if (!strcmp(option, "KLFM_InitPart_MaxNrLoops")) {
        pOptions->KLFM_InitPart_MaxNrLoops = atol(value);
    } else if (!strcmp(option, "KLFM_InitPart_MaxNrNoGainMoves")) {
        pOptions->KLFM_InitPart_MaxNrNoGainMoves = atol(value);
    } else if (!strcmp(option, "KLFM_Refine_MaxNrLoops")) {
        pOptions->KLFM_Refine_MaxNrLoops = atol(value);
    } else if (!strcmp(option, "KLFM_Refine_MaxNrNoGainMoves")) {
        pOptions->KLFM_Refine_MaxNrNoGainMoves = atol(value);
    } else if (!strcmp(option, "VectorPartition_Step3")) {
        if (!strcmp(value, "increase"))
            pOptions->VectorPartition_Step3 = VecIncrease;
        else if (!strcmp(value, "decrease"))
            pOptions->VectorPartition_Step3 = VecDecrease;
        else if (!strcmp(value, "random"))
            pOptions->VectorPartition_Step3 = VecRandom;
        else {
            fprintf(stderr, "SetOption(): unknown %s '%s'!\n", option, value);
            return FALSE;
        }
    } else if (!strcmp(option, "VectorPartition_MaxNrLoops")) {
        pOptions->VectorPartition_MaxNrLoops = atol(value);
    } else if (!strcmp(option, "VectorPartition_MaxNrGreedyImproves")) {
        pOptions->VectorPartition_MaxNrGreedyImproves  = atol(value);
    } else if (!strcmp(option, "OutputFormat")) {
        if (!strcmp(value, "original"))
            pOptions->OutputFormat = OutputDMM;
        else if(!strcmp(value, "EMM") || !strcmp(value, "emm") )
            pOptions->OutputFormat = OutputEMM;
        else {
            fprintf(stderr, "SetOption(): unkown %s '%s'!\n", option, value);
            return FALSE;
        }
    } else if (!strcmp(option, "Iterative_Refinement")) {
        if (!strcmp(value, "never"))
            pOptions->Iterative_Refinement = IR_Never;
        else if(!strcmp(value, "aftermg"))
            pOptions->Iterative_Refinement = IR_After_MG;
        else if(!strcmp(value, "always"))
            pOptions->Iterative_Refinement = IR_Always;
        else {
            fprintf(stderr, "SetOption(): unkown %s '%s'!\n", option, value);
            return FALSE;
        }
    } else if (!strcmp(option, "OutputMode")) {
        if (!strcmp(value, "onefile"))
            pOptions->OutputMode = OneFile;
        else if(!strcmp(value, "original"))
            pOptions->OutputMode = MultipleFiles;
        else if(!strcmp(value, "DIMACS"))
            pOptions->OutputMode = DIMACS;
        else {
            fprintf(stderr, "SetOption(): unknown %s '%s'!\n", option, value);
            return FALSE;
        }
    } else if (!strcmp(option, "Seed")) {
        if (!strcmp(value, "random") || atol(value) < 0)
            pOptions->Seed = -1;
        else
            pOptions->Seed = labs(atol(value)); 
    } else if (!strcmp(option, "Permute")) {
        if (!strcmp(value, "none")) {
            pOptions->OrderPermutation = OrderNone;
        } else if (!strcmp(value, "reverseBBD") || !strcmp(value, "prefix")) {
            pOptions->OrderPermutation = OrderPrefix;
        } else if (!strcmp(value, "SBD") || !strcmp(value, "infix")) {
            pOptions->OrderPermutation = OrderInfix;
        } else if (!strcmp(value, "BBD") || !strcmp(value, "postfix")) {
            pOptions->OrderPermutation = OrderPostfix;
        } else {
            fprintf(stderr, "SetOption(): unknown %s '%s'!\n", option, value);
            return FALSE;
        }
    } else if (!strcmp(option, "EnforceSymmetricPermutation")) {
        if (!strcmp(value, "yes")) {
            pOptions->SymmetricPermute = TRUE;
        } else if (!strcmp(value, "no")) {
            pOptions->SymmetricPermute = FALSE;
        } else {
            fprintf(stderr, "SetOption(): unknown %s '%s'!\n", option, value);
            return FALSE;
        }
    } else {
        fprintf(stderr, "SetOption(): unknown option '%s'!\n", option);
        return FALSE;
    }
    
    return TRUE;
} /* end SetOption */
  
  
void PrintHelp(int argc, char **argv) {

    /* This function prints a help message on how to use Mondriaan
       from the command line. This function can for instance be called
       in case there are too few arguments. */

    printf("\nMondriaan version %s.\n", MONDRIAANVERSION);
    if (argv && argv[0]) printf("Usage: %s [matrix] [P] [eps] <options>\n\n", argv[0]);
    else printf("Usage: ./Mondriaan [matrix] [P] [eps] <options>\n\n");
    printf(" [matrix] - the matrix to partition\n");
    printf(" [P]      - the number of processors\n");
    printf(" [eps]    - the maximum allowed load imbalance\n\n");
    printf(" <options>:");
    printf(" for details, see the User's Guide at\n");
    printf("            http://www.math.uu.nl/people/bisseling");
    printf("/Mondriaan/users_guide.html\n");
    printf("\n");
    fflush(stdout);

} /* end PrintHelp */

