#include "MatchMatchers.h"

/* Greedy matching algorithm. Visits all vertices in the given order and matches them to their heaviest unmatched neighbor. */
int MatchUsingGreedy(struct biparthypergraph *pHG, struct contraction *pC,
                     const long *iv, int *Matched,
                     const struct opts *pOptions,
                     int (*SetupData)(void **, struct biparthypergraph *, const struct opts *),
                     int (*FindNeighbor)(long *, double *,
                                         const struct biparthypergraph *, const struct contraction *, const struct opts *,
                                         const long, const int *,
                                         void *, long *, long *, double *),
                     int (*FreeData)(void *)) {
    /* This function creates a simple greedy matching.
        
        pHG :       The hypergraph in which we create the matching.
        pC :        The structure containing the groups of matched pairs of vertices.
        iv :        Order in which we visit pHG's vertices.
        Matched :   Array indicating which vertices are matched.
        pOptions :  Selected options for matching.
        SetupData : Function to setup the data for the neighbor finding algorithm (e.g. MatchStairwaySetup).
        FindNeighbor : Function to find an unmatched neighbor of a given vertex (e.g. FindNeighborStairway).
        FreeData :  Function to free the allocated data for the neighbor finding algorithm (e.g. MatchStairwayFree).
    */
    /* Working array to contain inner products. */
    long *Ip = (long *)malloc(pHG->NrVertices*sizeof(long));
    /* Working array to contain scaled inner products. */
    double *ScIp = (double *)malloc(pHG->NrVertices*sizeof(double));
    /* Working array to contain visited vertices. */
    long *Visited = (long *)malloc(pHG->NrVertices*sizeof(long));
    /* Data for neighbor finding. */
    void *pData = NULL;
    
    long t;
    
    if (!pHG || !pC || !iv || !Matched || !pOptions ||
        !Ip || !ScIp || !Visited) {
        fprintf(stderr, "MatchUsingGreedy(): Null arguments or not enough memory!\n");
        return FALSE;
    }
    
    if (pOptions->Coarsening_MaxNrVtxInMatch != 2) {
        fprintf(stderr, "MatchUsingGreedy(): This only works when matching groups of two vertices!\n");
        return FALSE;
    }
    
    if (!SetupData(&pData, pHG, pOptions)) {
        fprintf(stderr, "MatchUsingGreedy(): Unable to create neighbor finding data!\n");
        return FALSE;
    }
    
    /* Start with all vertices being unmatched. We will abuse the Matched array to indicate which vertices are part of the path. */
    for (t = 0; t < pHG->NrVertices; t++) {
        Matched[t] = FALSE;
        Ip[t] = 0;
        ScIp[t] = 0.0;
    }
    
    /* Greedily create a maximal matching. */
    for (t = 0; t < pHG->NrVertices; t++) {
        const long v = iv[t];
        
        if (Matched[v] == FALSE && !pHG->V[v].Free) {
            long w = -1;
            double ip = 0.0;
            const long k = pC->NrMatches;
            
            /* Remove v from the hypergraph. */
            pC->NrMatches++;
            pC->Start[k + 1] = pC->Start[k];
            pC->Match[pC->Start[k + 1]++] = v;
            Matched[v] = TRUE;
            
            if (!MoveVtxInNetAdjncy(pHG, v)) return FALSE;
            
            if (FindNeighbor(&w, &ip,
                             pHG, pC, pOptions,
                             v, Matched,
                             pData, Visited, Ip, ScIp)) {
                if (w < 0 || v == w) {
                    fprintf(stderr, "MatchUsingGreedy(): Unable to find a neighbor of an unmatched vertex!\n");
                    return FALSE;
                }
                
                /* Remove w from the hypergraph. */
                pC->Match[pC->Start[k + 1]++] = w;
                Matched[w] = TRUE;
                
                if (!MoveVtxInNetAdjncy(pHG, w)) return FALSE;
            }
        }
    }

    /* Free memory. */
    free(Ip);
    free(ScIp);
    free(Visited);
    
    FreeData(pData);
    
    return TRUE;
} /* end MatchUsingGreedy */

long FindOptimalPathMatching(char **Build, double *Weight[2], char *MatchingOpt,
                             const double *PathWeights, const long NrVtxInPath, const long MaxNrVtxInMatch)
{
    /* This function uses dynamic programming to create an optimal (generalized) matching on the given path.
     * Build           : Working array for constructing the maximum matching.
     * Weight          : Working array for registering accumulated weights during the process.
     * MatchingOpt     : Working array to return the maximum matching in (output).
     * PathWeights     : Weights of the edges
     * NrVtxInPath     : Number of vertices in the given path (= number of edges/weights plus 1)
     * MaxNrVtxInMatch : Maximum number of vertices in a generalized match (2 for the standard maximum matching problem)
     */
    
    long t, k, prev, curr;
    
    /* Change from number of vertices to number of edges */
    const long PathLength = NrVtxInPath-1;
    const long MaxCnsqEdges = MaxNrVtxInMatch-1;
    
    /* Initialize */
    Weight[0][0] = 0;
    for(k=1; k<=MaxCnsqEdges; ++k) {
        Weight[0][k] = PathWeights[0];
        Build[0][k] = 1;
    }
    
    /* Walk through path */
    for(t=1; t<PathLength; ++t) {
        curr = t%2; prev = (curr+1)%2;
        
        Weight[curr][0] = Weight[prev][MaxCnsqEdges];
        
        for(k=1; k<=MaxCnsqEdges; ++k) {
            if(Weight[prev][k-1]+PathWeights[t] > Weight[curr][0]) {
                Weight[curr][k] = Weight[prev][k-1]+PathWeights[t];
                Build[t][k] = 1;
            }
            else {
                Weight[curr][k] = Weight[curr][0];
                Build[t][k] = 0;
            }
        }
    }
    
    /* Build[PathLength-1][MaxCnsqEdges] now contains the maximum matching,
     * with total value equal to Weight[(PathLength-1)%2][MaxCnsqEdges].
     * Backtrack from there.
     */
    
    long nrMatchedEdges = 0;
    k = MaxCnsqEdges;
    for(t=PathLength-1; t>=0; --t) {
        if(k > 0 && Build[t][k]) {
            /* We chose to include edge t, and include at most k-1 previous edges */
            ++nrMatchedEdges;
            MatchingOpt[t] = 1;
            --k;
        }
        else {
            /* We chose to not include edge t, so we may include MaxCnsqEdges previous edges */
            MatchingOpt[t] = 0;
            k = MaxCnsqEdges;
        }
    }
    
    return nrMatchedEdges;
}

int ApplyPathMatching(struct contraction *pC, int *Matched,
                      const char *PathMatching, const long *PathIndices, const long _PathLength)
{
    /* Applies a (generalized) path matching generated by FindOptimalPathMatching, and adds this matching to the given list of contractions. */
    
    long t;
    
    /* Change from number of vertices to number of edges */
    const long PathLength = _PathLength-1;
    
    for (t = 0; t < PathLength; ++t) {
        if(!PathMatching[t]) {
            continue;
        }
        
        const long k = pC->NrMatches;
        
        pC->NrMatches++;
        pC->Start[k+1] = pC->Start[k];
        
        /* For each included edge, add the first indicent vertex to the (generalized) match */
        while(t < PathLength && PathMatching[t]) {
            const long v = PathIndices[t];
            pC->Match[pC->Start[k+1]++] = v;
            Matched[v] = TRUE;
            ++t;
        }
        
        /* The current edge is not included, so add first indicent vertex to
         * the (generalized) match, and continue to the next match */
        const long v = PathIndices[t];
        pC->Match[pC->Start[k+1]++] = v;
        Matched[v] = TRUE;
        /* Do not increase t, this is done by the for loop */
        
    }
    
    return TRUE;
}

int MatchUsingPGA(struct biparthypergraph *pHG, struct contraction *pC,
                  const long *iv, int *Matched,
                  const struct opts *pOptions,
                  int (*SetupData)(void **, struct biparthypergraph *, const struct opts *),
                  int (*FindNeighbor)(long *, double *,
                                      const struct biparthypergraph *, const struct contraction *, const struct opts *,
                                      const long, const int *,
                                      void *, long *, long *, double *),
                  int (*FreeData)(void *)) {
    /* This function creates a (1/2)-optimal matching using the PGA' algorithm by Drake and Hougardy (2003).
        
        pHG :       The hypergraph in which we create the matching.
        pC :        The structure containing the groups of matched pairs of vertices.
        iv :        Order in which we visit pHG's vertices.
        Matched :   Array indicating which vertices are matched.
        pOptions :  Selected options for matching.
        SetupData : Function to setup the data for the neighbor finding algorithm (e.g. MatchStairwaySetup).
        FindNeighbor : Function to find an unmatched neighbor of a given vertex (e.g. FindNeighborStairway).
        FreeData :  Function to free the allocated data for the neighbor finding algorithm (e.g. MatchStairwayFree).
    */
    /* Working array to contain inner products. */
    long *Ip = (long *)malloc(pHG->NrVertices*sizeof(long));
    /* Working array to contain scaled inner products. */
    double *ScIp = (double *)malloc(pHG->NrVertices*sizeof(double));
    /* Working array to contain visited vertices. */
    long *Visited = (long *)malloc(pHG->NrVertices*sizeof(long));
    
    /* Working arrays for creating optimal path matchings. */
    char *_PathMatchings = NULL, **PathMatchings = NULL, *PathMatchingOpt = NULL;
    double *_PathMatchingWeights = NULL, *PathMatchingWeights[2];
    /* Array containing the matching weights along the constructed path. */
    double *PathWeights = (double *)malloc(pHG->NrVertices*sizeof(double));
    /* Array containing the indices of the vertices along this path. */
    long *PathIndices = (long *)malloc(pHG->NrVertices*sizeof(long));
    long PathLength = 0;
    
    /* Data for neighbor finding. */
    void *pData = NULL;
    
    long t, tt;

    /* Check arguments */
    if (!pHG || !pC || !iv || !Matched || !pOptions) {
        fprintf(stderr, "MatchUsingPGA(): Null arguments!\n");
        return FALSE;
    }
    
    if (!Ip || !ScIp || !Visited || !PathWeights || !PathIndices) {
        fprintf(stderr, "MatchUsingPGA(): Not enough memory!\n");
        return FALSE;
    }
    
    long MaxNrVtxInMatch = pOptions->Coarsening_MaxNrVtxInMatch;
    if(MaxNrVtxInMatch == -1) {
        MaxNrVtxInMatch = log2(pHG->NrVertices / pOptions->Coarsening_NrVertices);
        if(MaxNrVtxInMatch < 2) {
            MaxNrVtxInMatch = 2;
        }
        if(MaxNrVtxInMatch > 4) {
            MaxNrVtxInMatch = 4;
        }
    }
    
    if (MaxNrVtxInMatch < 2) {
        fprintf(stderr, "MatchUsingPGA(): Invalid number of vertices per match!\n");
        return FALSE;
    }
    
    /* Allocation for FindOptimalPathMatching()
     * PathMatchings[t][0] for 'no edge selected at end'
     * PathMatchings[t][k] for 'max k edges selected at end', where k=1..MaxCnsqEdges
     * PathMatchingWeights[0/1] for weights in previous and current iteration
     * PathMatchingOpt for maximum matching
     */
    const long MaxCnsqEdges = MaxNrVtxInMatch - 1;
    
    /* Set up PathMatchings */
    _PathMatchings = (char *)malloc((MaxCnsqEdges+1) * pHG->NrVertices * sizeof(char));
    PathMatchings = (char **)malloc(pHG->NrVertices * sizeof(char *));
    if (!_PathMatchings || !PathMatchings) {
        fprintf(stderr, "MatchUsingPGA(): Not enough memory!\n");
        return FALSE;
    }
    for(t=0; t<pHG->NrVertices; ++t) {
        PathMatchings[t] = &(_PathMatchings[t*(MaxCnsqEdges+1)]);
    }
    
    /* Set up PathMatchingWeights */
    _PathMatchingWeights = (double *)malloc(2*(MaxCnsqEdges+1) * sizeof(double));
    if (!_PathMatchingWeights) {
        fprintf(stderr, "MatchUsingPGA(): Not enough memory!\n");
        return FALSE;
    }
    PathMatchingWeights[0] = &(_PathMatchingWeights[0]);
    PathMatchingWeights[1] = &(_PathMatchingWeights[MaxCnsqEdges+1]);
    
    /* Set up PathMatchingOpt */
    PathMatchingOpt = (char *)malloc(pHG->NrVertices * sizeof(char));
    if (!PathMatchingOpt) {
        fprintf(stderr, "MatchUsingPGA(): Not enough memory!\n");
        return FALSE;
    }
    
    if (!SetupData(&pData, pHG, pOptions)) {
        fprintf(stderr, "MatchUsingPGA(): Unable to create neighbor finding data!\n");
        return FALSE;
    }
    
    /* Start with all vertices being unmatched. We will abuse the Matched array to indicate which vertices are part of the path. */
    for (t = 0; t < pHG->NrVertices; t++) {
        Matched[t] = FALSE;
        Ip[t] = 0;
        ScIp[t] = 0.0;
    }
    
    /* Visit the vertices in the prescribed order and start a path at each unmatched vertex. */
    for (t = 0; t < pHG->NrVertices; t++) {
        long PathMatchingLength = -1;
        long CurrentVertex = iv[t];
        long NextVertex = -1;
        double NextWeight = 0.0;
        
        /* If this vertex has already been matched, we do nothing. */
        if (Matched[CurrentVertex] != FALSE || pHG->V[CurrentVertex].Free) {
            continue;
        }
        
        /* Start path. */
        Matched[CurrentVertex] = TRUE;
        PathIndices[0] = CurrentVertex;
        PathLength = 1;
        
        if (!MoveVtxInNetAdjncy(pHG, CurrentVertex)) return FALSE;

        /* Extend the path along the unmatched neighbor that has the highest inner product. */
        while (FindNeighbor(&NextVertex, &NextWeight,
                            pHG, pC, pOptions,
                            CurrentVertex, Matched,
                            pData, Visited, Ip, ScIp)) {
            if (NextVertex < 0 || NextVertex == CurrentVertex || NextWeight <= 0.0) {
                fprintf(stderr, "MatchUsingPGA(): An error has occurred while finding the next neighbor!\n");
                return FALSE;
            }
            
            if (Matched[NextVertex] != FALSE) {
                fprintf(stderr, "MatchUsingPGA(): Extended path using a matched vertex!\n");
                return FALSE;
            }
            
            /* Extend path. */
            Matched[NextVertex] = TRUE;
            PathWeights[PathLength - 1] = NextWeight;
            PathIndices[PathLength] = NextVertex;
            PathLength++;
            
            if (!MoveVtxInNetAdjncy(pHG, NextVertex)) return FALSE;

            CurrentVertex = NextVertex;
            NextVertex = -1;
        }
        
        /* Reset matched flags. */
        for (tt = 0; tt < PathLength; tt++) {
            Matched[PathIndices[tt]] = FALSE;
        }
        
        /* Only match if the path contains at least two vertices. */
        if (PathLength <= 1) {
            for (tt = 0; tt < PathLength; tt++) {
                if (!MoveVtxBackInNetAdjncy(pHG, PathIndices[tt])) return FALSE;
            }
            continue;
        }
        
        /* Find optimal weight (generalized) matching along the path using dynamic programming. */
        PathMatchingLength = FindOptimalPathMatching(PathMatchings, PathMatchingWeights, PathMatchingOpt,
                                                     PathWeights, PathLength, MaxNrVtxInMatch);
        
        if (PathMatchingLength < 0) {
            fprintf(stderr, "MatchUsingPGA(): Unable to find a matching along the path!\n");
            return FALSE;
        }
        
        if (!ApplyPathMatching(pC, Matched, PathMatchingOpt, PathIndices, PathLength)) {
            fprintf(stderr, "MatchUsingPGA(): Unable to apply path matching!\n");
            return FALSE;
        }
        
        for (tt = 0; tt < PathLength; tt++) {
            if(!Matched[PathIndices[tt]]) {
                if (!MoveVtxBackInNetAdjncy(pHG, PathIndices[tt])) return FALSE;
            }
        }
    }
    
    /* Greedily extend matching to a maximal matching. */
    for (t = 0; t < pHG->NrVertices; t++) {
        const long v = iv[t];
        
        if (Matched[v] == FALSE && !pHG->V[v].Free) {
            long w = -1;
            double ip = 0.0;
            const long k = pC->NrMatches;
            
            if (FindNeighbor(&w, &ip,
                             pHG, pC, pOptions,
                             v, Matched,
                             pData, Visited, Ip, ScIp)) {
                if (w < 0 || v == w) {
                    fprintf(stderr, "MatchUsingPGA(): Unable to find a neighbor of an unmatched vertex!\n");
                    return FALSE;
                }
                
                /* Create a new matching group containing v and w. */
                pC->NrMatches++;
                pC->Start[k + 1] = pC->Start[k];
                pC->Match[pC->Start[k + 1]++] = v;
                pC->Match[pC->Start[k + 1]++] = w;
                Matched[v] = TRUE;
                Matched[w] = TRUE;
                
                /* Remove v and w from the hypergraph. */
                if (!MoveVtxInNetAdjncy(pHG, v)) return FALSE;
                if (!MoveVtxInNetAdjncy(pHG, w)) return FALSE;
            }
            else {
                /* This vertex has no available neighbor. */
                pC->NrMatches++;
                pC->Start[k + 1] = pC->Start[k];
                pC->Match[pC->Start[k + 1]++] = v;
                Matched[v] = TRUE;
                
                if (!MoveVtxInNetAdjncy(pHG, v)) return FALSE;
            }
        }
    }

    /* Free memory. */
    free(Ip);
    free(ScIp);
    free(Visited);
    free(PathWeights);
    free(PathIndices);
    free(_PathMatchings);
    free(PathMatchings);
    free(_PathMatchingWeights);
    free(PathMatchingOpt);
    
    FreeData(pData);
    
    return TRUE;
} /* end MatchUsingPGA */


