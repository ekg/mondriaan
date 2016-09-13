/* This header file defines the structures vertex, net,
   and bipartitioned hypergraph.  */

#ifndef __Graph_h__
#define __Graph_h__

#include "Sort.h"
#include "Options.h"
#include "GainBucket.h"
#include "SparseMatrix.h"

/*** Hypergraph data structures ***/

/* A vertex contains a weight, a list of adjacent nets 
   (i.e. those nets that contain the vertex), the partition number
   in a bipartitioning, and a pointer to its bucketentry
   in the gainbucket data structure. */

struct vertex {
    long vtxwgt;   /* weight of the vertex */
    long iStart;   /* index into VtxAdjncy where the list begins */
    long iEnd  ;   /* index where it ends, such that iStart <= i < iEnd
                       gives the complete list */
    int partition; /* partition of the vertex (either 0 or 1) */
    int Free; /* whether or not this is a free vertex. */
    struct bucketentry *GBentry; /* pointer to the bucketentry in the
                                     gainbucket data structure */
};
  

/* A net contains a weight, a list of vertices, with those
   of partition 0 coming first, and then those of partition 1,
   and the direction at its creation from a matrix row or column. */

struct net {
    long netwgt;   /* weight of the net */
    long iStartP0; /* index into NetAdjncy where the list for partition 0 begins */
    long iStartP1; /* index where the list for partition 1 begins */
    long iEnd;     /* index where the list ends */
    int dir;       /* net direction (ROW for a row-net, or
                                      COL for a column-net) */
    int Free;      /* is this a free net (that is, should it be discarded by HKLFM?) */
};


/* A bipartitioned hypergraph contains vertices, nets, and 
   adjacency lists corresponding to the pins of the hypergraph.

   Furthermore, it can store information on the current situation
   in the bipartitioning process by the Kernighan-Lin/Fiduccia-Mattheyses
   algorithm. 

   Finally, it can store the numerical values of the nonzeros
   corresponding to the pins. This is usually only done for the hypergraph
   created from the matrix to be bipartitioned. Hypergraphs created by coarsening
   other hypergraphs usually do not store numerical information. */

struct biparthypergraph {
    /* General hypergraph variables */
    long NrVertices; /* number of vertices of the hypergraph */
    long NrNets;     /* number of nets */
    long NrPins;      /* number of pins (nonzeros in the corresponding matrix) */
  
    struct vertex *V; /* stores the information about the vertices */
    long *VtxAdjncy;  /* stores the adjacency lists of the vertices */
  
    struct net *N;    /* stores the information about the nets */
    long *NetAdjncy;  /* stores the adjacency lists of the nets */

    long WeightP[2];  /* stores the total weight for the vertex partitions */
  
    /* Variables of the Kernighan-Lin/Fiduccia-Mattheyses algorithm */
    struct gainbucket GBVtx[2]; /* stores the gainbuckets for the vertex
                                    partitions */
  
    long *VtxMoveLog; /* stores a log of the moved vertices */
    long CurVtxLog, MinVtxLog; /* indices for VtxMoveLog */
  
    long CurComm; /* stores current communication volume */
    long MinComm; /* stores minimum volume of one iteration */
    long OptComm; /* stores optimal volume (minimum over all iterations) */
    int *OptPartVtx; /* stores the optimal partition for the vertices */
  
    int SplitDir; /* stores the split direction of the whole hypergraph
                      (ROW, COL, or FINEGRAIN) created from a matrix */
  
    /* Matrix values and indices (not used for a coarser hypergraph) */
    double *MatReValue; /* stores the real part of the values of the
                            pins (nonzero matrix entries) in the same way
                            as the adjacency lists of the vertices,
                            for a hypergraph created in the ROW or COL
                            direction. 

                            stores the values in the same way as the vertices,
                            for a hypergraph created in the FINEGRAIN direction
                            (where each vertex represents a nonzero).  */

    double *MatImValue; /* stores the imaginary part of the values */
    long *Vtx2MatIndex; /* stores the index of the corresponding
                            matrix row/column for the vertices,
                            in case of the ROW/COL direction, respectively.
                            Not used for FINEGRAIN.  */
    long *Net2MatIndex; /* stores the index of the corresponding
                            matrix column/row for the nets, 
                            in case of the ROW/COL direction, respectively.
                            For FINEGRAIN, each net has its own
                            direction (ROW or COL). */
};
  

/*** Function declarations for Graph.c ***/
int CreateNewBiPartHyperGraph(long NrVertices, long NrNets, long NrPins, int StoreMat, char MType, struct biparthypergraph *pHG);
int DeleteBiPartHyperGraph(struct biparthypergraph *pHG);
int SparseMatrix2BiPartHyperGraph(const struct sparsematrix *pM, int dir, const struct opts *pOptions, struct biparthypergraph *pHG);
int BiPartHyperGraph2SparseMatrix(const struct biparthypergraph *pHG, long lo, long hi, long *mid, struct sparsematrix *pM);
int SortAdjacencyLists(struct biparthypergraph *pHG);
  
#endif /* __Graph_h__ */
