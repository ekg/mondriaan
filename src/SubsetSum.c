
#include "SubsetSum.h"
#include "Heap.h"

int _SubsetSumExp(const long * const weights, long * const select, const long N, const long pos, const long weightSelected, const long weightNotSelected, const long totWeight, const long minWeight, const long maxWeight);

/**
 * Perform a subset sum calculation.
 * This is a variant of subset sum where the desired subset should have sum between minWeight and maxWeight,
 * as opposed to more common variants where the sum should be equal to a fixed weight, or where it should be
 * simply maximised (in non-decision form).
 * 
 * Input:
 * N                 : Number of items/weights
 * weightlo          : The weight the smaller subset should at most have
 * weighthi          : The weight the larger subset should at most have
 * 
 * Input/Output:
 * weights           : The weights of the items (may be permuted)
 * 
 * Output:
 * permutation       : The permutation of the items used by the algorithm (Only if return value == TRUE)
 * select            : Array defining which items are included in the small subset (1) or in the large subset (0) (Only if return value == TRUE)
 * 
 * Returned memory allocations:
 *     permutation:    O(N) (Only if return value == TRUE)
 *     select:         O(N) (Only if return value == TRUE)
 * 
 * Return value: TRUE if a feasible subset has been found, FALSE if this has not been found or an error occurred.
 */
int SubsetSum(long * const weights, const long N, const long weightlo, const long weighthi, long **permutation, long **select) {
    
    if(N > 16) {
        return KarmarkarKarp(weights, N, weightlo, weighthi, permutation, select);
    }
    
    return SubsetSumExp(weights, N, weightlo, weighthi, permutation, select);
    
} /* end SubsetSum */

/**
 * Perform a subset sum calculation, using an exponential time algorithm. This should not be used for large N!
 * 
 * See SubsetSum()
 */
int SubsetSumExp(long * const weights, const long N, const long weightlo, const long weighthi, long **permutation, long **select) {
    
    long totalWeight = 0, i;
    for(i=0; i<N; i++) {
        totalWeight += weights[i];
    }
    
    /* Compute min and max weight of smallest desired partition */
    long maxWeight = weightlo;
    long minWeight = totalWeight - weighthi;
    
    /* Sort the weights in ascending order */
    *permutation = QSort(weights, N);
    
    if(*permutation == 0) {
        return FALSE;
    }
    
    for (long t = 0; t < N/2; ++t) {
        SwapLong(weights, t, N-1-t);
        SwapLong(*permutation, t, N-1-t);
    }
    
    /* Allocate and zero-initialize the decision array */
    *select = (long*)calloc(N, sizeof(long));
    
    if( *select == NULL ) {
        free(*permutation);
        fprintf( stderr, "SubsetSumExp(): Not enough memory.\n" );
        return FALSE;
    }
    
    /* Run subset sum */
    int result = _SubsetSumExp(weights, *select, N, 0, 0, 0, totalWeight, minWeight, maxWeight);
    if(!result) {
        free(*permutation);
        free(*select);
        return FALSE;
    }
    
    return TRUE;
    
} /* end SubsetSumExp */

/**
 * Perform a subset sum recursion. See also SubsetSumExp().
 * 
 * Input:
 * weights           : The weights of the items
 * N                 : Number of items/weights
 * pos               : Current item under consideration
 * weightSelected    : Weight selected from elements [0,pos-1]
 * weightNotSelected : Weight not selected from elements [0,pos-1]
 * totWeight         : Total weight of all N items
 * minWeight         : The weight the desired subset should at least have
 * maxWeight         : The weight the desired subset should at most have
 * 
 * Input/Output:
 * select            : Array defining which items are included in the subset (1) or not (0)
 * 
 * Return value: TRUE if a feasible subset has been found in this recursion or further down
 *               in the current recursion tree, FALSE otherwise.
 */
int _SubsetSumExp(const long * const weights, long * const select, const long N, const long pos, const long weightSelected, const long weightNotSelected, const long totWeight, const long minWeight, const long maxWeight) {
    
    /* (A)
     * By (C), we know that weightSelected <= maxWeight. If also
     * minWeight <= weightSelected, we have found a feasible subset.
     */
    if(weightSelected >= minWeight) {
        return TRUE;
    }
    
    /* (B)
     * End of tree
     */
    if(pos >= N) {
        return FALSE;
    }
    
    /* (C)
     * All items are sorted by increasing weight. Hence, if adding
     * the new weight will exceed the maximum weight, any subsequent
     * weight will exceed the maximum weight too. So we may return FALSE.
     */
    if(weightSelected+weights[pos] > maxWeight) {
        return FALSE;
    }
    
    /* (D)
     * Check whether there is enough weight left to reach minWeight.
     */
    if(minWeight > totWeight-weightNotSelected) {
        return FALSE;
    }
    
    /* Branch: Select the current item */
    select[pos] = 1;
    if(_SubsetSumExp(weights, select, N, pos+1, weightSelected+weights[pos], weightNotSelected, totWeight, minWeight, maxWeight)) {
        return TRUE;
    }
    
    /* Branch: Don't select the current item */
    select[pos] = 0;
    if(_SubsetSumExp(weights, select, N, pos+1, weightSelected, weightNotSelected+weights[pos], totWeight, minWeight, maxWeight)) {
        return TRUE;
    }
    
    /* If in both branches no feasible subset is found, return FALSE */
    return FALSE;
    
} /* end _SubsetSumExp */


/**
 * Comparison function, taking two KKWeights, and comparing the values
 * a->weight and b->weight in an ascending fashion.
 * Used in KarmarkarKarp().
 * 
 * Input:
 * a           : The left value
 * b           : The right value
 * 
 * Return value:
 *      -1 if a->weight < b->weight
 *       0 if a->weight = b->weight
 *       1 if a->weight > b->weight
 */
int compareKKWeights (const void *a, const void *b) {
    long diff = ((struct KKWeight*)a)->weight - ((struct KKWeight*)b)->weight;
    if(diff == 0)
        return 0;
    return (diff < 0) ? -1 : 1;
} /* end compareKKWeights */

/**
 * Perform a subset sum using the Karmarkar-Karp differencing algorithm.
 * This is a variant of the Karmarkar Karp algorithm that can be applied to the variant on Subset Sum as
 * explained in SubsetSum(). This heuristic is faster than SubsetSum(), but may not find a feasible subset
 * even if one exists.
 * A dummy weight is used to convert the subset sum problem to the number partitioning problem, on which
 * Karmarkar-Karp can be applied.
 * 
 * See SubsetSum()
 */
int KarmarkarKarp(const long * const weights, const long N, const long weightlo, const long weighthi, long **permutation, long **select) {
    
    long totalWeight = 0, i;
    for(i=0; i<N; i++) {
        totalWeight += weights[i];
    }
    
    /* Compute min and max weight of smallest desired partition */
    long maxWeight = weightlo;
    long minWeight = totalWeight - weighthi;
    
    /* If the two desired subsets have unequal size, we should add a dummy weight to make them equal */
    int useDummyWeight = (weightlo != weighthi)?1:0;
    
    /* Initialization */
    struct KKWeight *heapItems = (struct KKWeight*)malloc((N+useDummyWeight)*sizeof(struct KKWeight));
    if(heapItems == NULL) {
        fprintf( stderr, "KarmarkarKarp(): Not enough memory.\n" );
        return FALSE;
    }
    
    struct Graph Tree;
    GraphInit(&Tree, N+useDummyWeight);
    
    /* Register all weights */
    for(i=0; i<N; ++i) {
        heapItems[i].key = i;
        heapItems[i].weight = weights[i];
    }
    /* If necessary, also register the dummy weight */
    if(useDummyWeight) {
        heapItems[N].key = N;
        heapItems[N].weight = totalWeight - (minWeight + maxWeight);
    }
    
    /* Create the heap */
    struct heap Heap;
    HeapInit(&Heap, sizeof(struct KKWeight), compareKKWeights);
    Heapify(&Heap, heapItems, N+useDummyWeight);
    
    /* A heap item exists of two long values, the first holding an identification number (key),
     * the second holding the corresponding weight.
     */
    struct KKWeight max1, max2; /* Key and weight of items with highest and second highest weight */
    
    /* The Karmarkar Karp differencing loop.
     * In each iteration, the two largest weights are taken, of which the smaller weight is
     * subtracted from both, leaving the larger weight in the heap (with reduced weight),
     * and removing the smaller weight from the heap.
     * Then between the two selected weights an edge is created, to signify that we put them in
     * different sets.
     */
    while(Heap.numItems > 1) {
        /* Highest weight */
        HeapPop(&Heap, &max1);
        
        /* Second highest weight */
        HeapPeek(&Heap, &max2);
        
        /* Compute difference and replace max2 with it */
        max1.weight -= max2.weight;
        HeapReplace(&Heap, &max1, NULL);
        
        /* Connect the two weights */
        GraphAddEdge(&Tree, max1.key, max2.key);
    }
    
    /* We are left with one differencing value. The weight of it
     * equals the difference in total weight of the two subsets.
     * We check that it is sufficiently small.
     */
    HeapPeek(&Heap, &max1);
    
    if(max1.weight > (maxWeight-minWeight)/2) {
        HeapDestroy(&Heap); /* Also free()s heapItems */
        GraphDestroy(&Tree);
        return FALSE;
    }
    
    HeapDestroy(&Heap); /* Also free()s heapItems */
    
    /* Assign colours (0/1) to the tree such that no two neighbours have equal colours */
    TwoColorTree(&Tree);
    
    /* select[i] should be 1 if weight[i] belongs to the small subset, and 0 if it belongs
     * to the large subset. Register the colour of the small subset.
     */
    long smallSubsetColour = (useDummyWeight) ? Tree.nodes[N].val : 1;
    
    /* Generate the output */
    *permutation = (long*)malloc(N*sizeof(long));
    *select = (long*)malloc(N*sizeof(long));
    
    for(i=0; i<N; ++i) {
        (*permutation)[i] = i;
        (*select)[i] = (smallSubsetColour==0) ? 1 - Tree.nodes[i].val : Tree.nodes[i].val;
    }
    
    GraphDestroy(&Tree);
    
    return TRUE;
} /* end KarmarkarKarp */


/**
 * Initialize a graph struct
 * This assumes that the struct is not already initialized.
 * 
 * Input:
 * numNodes          : The number of nodes to allocate memory for
 * 
 * Output:
 * pGraph            : The graph struct
 * 
 * Returned memory allocations:
 *     pGraph:       O(numNodes)    (multiple mallocs)
 */
void GraphInit(struct Graph* pGraph, long numNodes) {
    pGraph->numNodes = numNodes;
    pGraph->numEdges = 0;
    pGraph->nodes = (struct Node*)malloc(numNodes*sizeof(struct Node));
    
    long i;
    for(i=0; i<pGraph->numNodes; ++i) {
        pGraph->nodes[i].key = i;
        pGraph->nodes[i].val = 0;
        pGraph->nodes[i].deg = 0;
        pGraph->nodes[i].adjSize = 4;
        pGraph->nodes[i].adjList = (long*)malloc(pGraph->nodes[i].adjSize*sizeof(long));
    }
}

/**
 * Free a graph struct
 * 
 * Input/Output:
 * pGraph            : The graph struct
 */
void GraphDestroy(struct Graph* pGraph) {
    long i;
    for(i=0; i<pGraph->numNodes; ++i) {
        free(pGraph->nodes[i].adjList);
    }
    free(pGraph->nodes);
    pGraph->numNodes = 0;
}

/**
 * Check the size of the adjacency list of a node. If the adjacency
 * list is full, resize it.
 * 
 * Input:
 * i                 : The node number
 * 
 * Input/Output:
 * pGraph            : The graph struct
 */
void GraphCheckAdjSize(struct Graph* pGraph, long i) {
    if(pGraph->nodes[i].deg == pGraph->nodes[i].adjSize) {
        pGraph->nodes[i].adjSize *= 2;
        pGraph->nodes[i].adjList = (long*)realloc(pGraph->nodes[i].adjList, pGraph->nodes[i].adjSize*sizeof(long));
        if(pGraph->nodes[i].adjList == NULL) {
            fprintf( stderr, "GraphCheckAdjSize(): Not enough memory.\n" );
            exit(1);
        }
    }
}

/**
 * Add an edge to the graph, by adding the nodes of the
 * edge to both their adjacency lists.
 * 
 * Input:
 * u, v              : The node numbers of the nodes in the edge
 * 
 * Input/Output:
 * pGraph            : The graph struct
 */
void GraphAddEdge(struct Graph* pGraph, long u, long v) {
    GraphCheckAdjSize(pGraph, u);
    GraphCheckAdjSize(pGraph, v);
    
    pGraph->nodes[u].adjList[pGraph->nodes[u].deg++] = v;
    pGraph->nodes[v].adjList[pGraph->nodes[v].deg++] = u;
}

/**
 * Assign values 0 and 1 to nodes, in such a way that every two
 * adjacent nodes have different values (colours). Node 0 is
 * assigned colour 0, with which the entire colouring is
 * unambiguously defined.
 * This function assumes the inputted graph is a tree. Else,
 * the obtained colouring may not be valid.
 * 
 * Input/Output:
 * pGraph            : The graph struct
 */
void TwoColorTree(struct Graph* pGraph) {
    long color = 0;
    long *nodeList = (long*)malloc(pGraph->numNodes * sizeof(long));
    long *currListStart, *currListEnd, *nextListEnd;
    nodeList[0] = 0;
    currListStart = nodeList;
    currListEnd = &nodeList[1];
    nextListEnd = &nodeList[1];
    long nodeKey, i;
    
    for(i=0; i<pGraph->numNodes; ++i) {
        pGraph->nodes[i].val = -1;
    }
    
    /* nodeList is of size such that each node may appear exactly once in it.
     * It is constructed as follows:
     * +------------------+--------------+-------------+------------+
     * |     Visited      |   Visiting   |  To visit   | Free space |
     * +------------------+--------------+-------------+------------+
     * ^                  ^              ^             ^
     * nodeList           currListStart  currListEnd   nextListEnd
     * 
     * `Visiting` contains all nodes that we are currently applying `color` to.
     * `To visit` contains all nodes that we will be applying the other color to in the next iteration.
     */
    
    while(currListStart < currListEnd) {
        
        /* Visit all nodes we should be `Visiting` */
        while(currListStart < currListEnd) {
            nodeKey = *(currListStart++);
            pGraph->nodes[nodeKey].val = color;
            
            /* We are required `To visit` the adjacent nodes in the next iteration */
            for(i=0; i<pGraph->nodes[nodeKey].deg; ++i) {
                if(pGraph->nodes[pGraph->nodes[nodeKey].adjList[i]].val == -1) {
                    *(nextListEnd++) = pGraph->nodes[nodeKey].adjList[i];
                }
            }
            
        }
        
        /* Next iteration: change `To visit` to `Visiting` */
        color = (color+1)%2;
        currListEnd = nextListEnd;
    }
    
    free(nodeList);
}
