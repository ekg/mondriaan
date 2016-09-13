#ifndef __Match_Stairway_h__
#define __Match_Stairway_h__

#include "Graph.h"
#include "Options.h"
#include "Match.h"

int MatchStairwaySetup(void **ppData,
                       struct biparthypergraph *pHG,
                       const struct opts *pOptions);

int FindNeighborStairway(long *pNeighbor, double *pNeighborIp,
                         const struct biparthypergraph *pHG, const struct contraction *pC, const struct opts *pOptions,
                         const long v, const int *Matched,
                         void *pData, long *Visited, long *Inprod, double *ScInprod);

int MatchStairwayFree(void *pData);

#endif /* __Match_Stairway_h__ */

