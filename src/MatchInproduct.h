#ifndef __Match_Inproduct_h__
#define __Match_Inproduct_h__

#include "Graph.h"
#include "Options.h"
#include "Match.h"

int MatchInproductSetup(void **ppData,
                       struct biparthypergraph *pHG,
                       const struct opts *pOptions);

int FindNeighborInproduct(long *pNeighbor, double *pNeighborIp,
                         const struct biparthypergraph *pHG, const struct contraction *pC, const struct opts *pOptions,
                         const long v, const int *Matched,
                         void *pData, long *Visited, long *Inprod, double *ScInprod);

int MatchInproductFree(void *pData);

#endif /* __Match_Inproduct_h__ */

