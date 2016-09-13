#ifndef __Match_Matchers_h__
#define __Match_Matchers_h__

#include "Graph.h"
#include "Options.h"
#include "Match.h"


int MatchUsingGreedy(struct biparthypergraph *pHG, struct contraction *pC,
                     const long *iv, int *Matched,
                     const struct opts *pOptions,
                     int (*SetupData)(void **, struct biparthypergraph *, const struct opts *),
                     int (*FindNeighbor)(long *, double *,
                                         const struct biparthypergraph *, const struct contraction *, const struct opts *,
                                         const long, const int *,
                                         void *, long *, long *, double *),
                     int (*FreeData)(void *));


long FindOptimalPathMatching(long *matchings[3], const double *PathWeights, const long PathLength);

int MatchUsingPGA(struct biparthypergraph *pHG, struct contraction *pC,
                  const long *iv, int *Matched,
                  const struct opts *pOptions,
                  int (*SetupData)(void **, struct biparthypergraph *, const struct opts *),
                  int (*FindNeighbor)(long *, double *,
                                      const struct biparthypergraph *, const struct contraction *, const struct opts *,
                                      const long, const int *,
                                      void *, long *, long *, double *),
                  int (*FreeData)(void *));

#endif /* __Match_Matchers_h__ */

