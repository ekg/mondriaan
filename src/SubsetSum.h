#ifndef __SubsetSum_h__
#define __SubsetSum_h__

#include "Sort.h"

/*** Graph data structures ***/
struct Node {
    long key;
    long val;
    long deg;
    long adjSize;
    long *adjList;
};

struct Graph {
    long numNodes;
    long numEdges;
    struct Node *nodes;
};

/* Weight used in KarmarkarKarp() */
struct KKWeight {
    long key;
    long weight;
};

/* Function declarations for SubsetSum.c */
int SubsetSum(long * const weights, const long N, const long weightlo, const long weighthi, long **permutation, long **select);
int SubsetSumExp(long * const weights, const long N, const long weightlo, const long weighthi, long **permutation, long **select);
int KarmarkarKarp(const long * const weights, const long N, const long weightlo, const long weighthi, long **permutation, long **select);

/* Function declarations for Graph, used in KarmarkarKarp() */
void GraphInit(struct Graph* pGraph, long numNodes);
void GraphDestroy(struct Graph* pGraph);
void GraphCheckAdjSize(struct Graph* pGraph, long i);
void GraphAddEdge(struct Graph* pGraph, long u, long v);
void TwoColorTree(struct Graph* pGraph);

#endif /* __SubsetSum_h__ */

