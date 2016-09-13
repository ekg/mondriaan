#ifndef __Options_h__
#define __Options_h__

/* This header file includes standard libraries,
   defines universal constants ROW, COL, FALSE, TRUE,
   and maximum word and line lengths and macros MIN, MAX. */

/* Standard libraries */
#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <math.h>
#include <limits.h>
#include <string.h>
#include <time.h>

#ifdef UNIX
#include <sys/time.h>
#endif

/* Directions in matrix */
#define ROW   0
#define COL   1 
#define FINEGRAIN 2
#define SFINEGRAIN 3
#define MEDIUMGRAIN 4

/* Logical constants */
#ifndef FALSE
#define FALSE 0
#endif
#ifndef TRUE
#define TRUE  1
#endif

/* Dummy negative value, cannot be an index.
Not to be confused with dummy diagonal nonzero */
#define DUMMY -1

/* Maximum lengths */
#define MAX_WORD_LENGTH 1024
#define MAX_LINE_LENGTH 16384

struct opts {
    /* General options */
    char *matrix;  /* name of the matrix */
    int P;         /* number of processors */
    double eps;    /* allowed load imbalance */
  
    /* Split strategy options */
    enum {Constant, Increase, Decrease} LoadbalanceStrategy; 
    enum {AdjustNo, AdjustYes} LoadbalanceAdjust; 
    enum {Alternate, LocalBest, Hybrid, LocalRatio,
          OneDimRow, OneDimCol, FineGrain, SFineGrain, MediumGrain} SplitStrategy;
    enum {FirstDirRow, FirstDirCol, FirstDirRatio} Alternate_FirstDirection; 
    enum {Simple, KLFM} SplitMethod;
    enum {PartMondriaan, PartPaToH, FullPaToH} Partitioner;
    enum {MetricLambda, MetricCut, MetricLambdaLambdaMinusOne} Metric;
    enum {FreeNetYes, FreeNetNo} DiscardFreeNets;
  
    /* Matrix options */
    enum {EqVecNo, EqVecYes} SquareMatrix_DistributeVectorsEqual;
    enum {DumNo, DumYes} SquareMatrix_DistributeVectorsEqual_AddDummies;
    enum {SingleEntNo, SingleEntYes} SymmetricMatrix_UseSingleEntry;
    enum {ETypeLower, ETypeRandom} SymmetricMatrix_SingleEntryType;
  
    /* Output options */
    enum {OutputEMM, OutputDMM} OutputFormat;
    enum {OneFile, MultipleFiles, DIMACS} OutputMode; /* Extended, (original) distributed MatrixMarket, or DIMACS output. */
  
    /* Coarsening options */
    long Coarsening_NrVertices, Coarsening_MaxCoarsenings;
    long Coarsening_NrMatchArbitrary;
    long Coarsening_MaxNrVtxInMatch;
    double Coarsening_StopRatio;
    double Coarsening_VtxMaxFractionOfWeight; 
    enum {MatchRandom, MatchInprod, MatchATA} Coarsening_MatchingStrategy;
    enum {MatchMatcherGreedy, MatchMatcherPGA} Coarsening_MatchingATAMatcher;
    enum {MatchFinderInproduct, MatchFinderStairway} Coarsening_MatchingATAFinder;
    enum {DecreasingWgt, IncreasingWgt, DecreasingDegree, IncreasingDegree,
          NaturalOrder, RandomOrder} Coarsening_InprodMatchingOrder;
    long Coarsening_FineSwitchLevel;
    enum {NoNetScaling, NetSclLinear} Coarsening_NetScaling;
    enum {NoIpScaling, IpSclCos, IpSclMin, IpSclMax, IpSclJaccard} Coarsening_InprodScaling;
    enum {MatchIdNo, MatchIdYes} Coarsening_MatchIdenticalFirst;
  
    /* Initial partitioning options */
    long KLFM_InitPart_NrRestarts;
    long KLFM_InitPart_MaxNrLoops;
    long KLFM_InitPart_MaxNrNoGainMoves;
  
    /* Refinement options */
    long KLFM_Refine_MaxNrLoops;
    long KLFM_Refine_MaxNrNoGainMoves;
    
    /* Iterative refinement options */
    enum {IR_Never, IR_After_MG, IR_Always} Iterative_Refinement;
  
    /* Vector partitioning options */
    enum {VecIncrease, VecDecrease, VecRandom} VectorPartition_Step3;
    long VectorPartition_MaxNrLoops;
    long VectorPartition_MaxNrGreedyImproves;
    
    long Seed;
    
    /* Permutation options */
    enum {OrderNone, OrderPrefix, OrderInfix, OrderPostfix} OrderPermutation;
    char SymmetricPermute;
};
  
/* Function declarations for Options.c */
int   GetParameters(struct opts *pOptions, int argc, char **argv);
char* GetDefaultOptionText();
int   SetOptions(struct opts *pOptions, const char *Text);
int   SetDefaultOptions(struct opts *pOptions);
int   SetOptionsFromFile(struct opts *pOptions, const char *File);
int   ExportOptions(FILE *Out, const struct opts *pOptions);
int   ExportDefaultOptions(FILE *Out);
int   ExportOptionsToLaTeX(FILE *Out, const struct opts *Opts);
int   SetOption(struct opts *pOptions, const char *option, const char *value);
int   ApplyOptions(const struct opts *pOptions);
void  PrintHelp(int argc, char **argv);

#endif /* __Options_h__ */
