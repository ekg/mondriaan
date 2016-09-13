#ifndef __HKLFM_h__
#define __HKLFM_h__

#include "Options.h"
#include "Graph.h"
#include "GainBucket.h"
#include "Match.h"
#include "MatchMatchers.h"
#include "MatchInproduct.h"
#include "MatchStairway.h"

/*** function declarations for HKLFM.c ***/
/*
int CreateInitialBalancedPartition(struct biparthypergraph *pHG, long weightlo, long weighthi);
int ComputeInitialGains(struct biparthypergraph *pHG);
int ClearMoveLog(struct biparthypergraph *pHG);
int MoveVertex(struct biparthypergraph *pHG, long v);
int UpdateGains(struct biparthypergraph *pHG, long v);
int HKLFM(struct biparthypergraph *pHG, long weightlo, long weighthi, int MaxNrLoops, int MaxNrNoGainMoves);
*/

int RunHKLFM(struct biparthypergraph *pHG, long weightlo, long weighthi, int refine, const struct opts *pOptions);
int CoarsenGraph(struct biparthypergraph *pHG, struct biparthypergraph *pCHG, struct contraction *pC, long level, const struct opts *pOptions);
int UncoarsenGraph(struct biparthypergraph *pCHG, struct contraction *pC, struct biparthypergraph *pHG, const long weightlo, const long weighthi, const struct opts *pOptions);
int RunMLGraphPart(struct biparthypergraph *pHG, long weightlo, long weighthi, const struct opts *pOptions);

#endif /* __HKLFM_h__ */
