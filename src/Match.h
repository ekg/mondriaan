#ifndef __Match_h__
#define __Match_h__

#include "Options.h"
#include "Graph.h"

/*** Contraction data structure ***/
struct contraction {
  long *Match;        /* Array containing the vertex numbers of
                          matched vertices, with the vertices of each
                          group stored contiguously */
  long *Start;        /* Start[k] = the start in array Match of the
                          vertices of group k */
  long NrMatches;     /* Nr of groups of matched vertices */
  long MaxNrVertices; /* Maximum number of vertices allowed in a group.
                          For pairwise matching, MaxNrVertices = 2 */
  long MaxVtxWgt;     /* Maximum total vertex weight allowed in a group */
};

/*** Function declarations for Match.c ***/
int MoveVtxInNetAdjncy(struct biparthypergraph *pHG, const long v);
int FindMatchArbitrary(struct biparthypergraph *pHG, struct contraction *pC,
                        const long v, int *Matched);
int FindMatchInprod(struct biparthypergraph *pHG, struct contraction *pC,
                     const long v, int *Matched,
                     long *Visited, long *Inprod, double *ScInprod, const struct opts *pOptions);

#endif /* __Match_h__ */
