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


long FindOptimalPathMatching(long *matchings[3], const double *PathWeights, const long PathLength)
{
    /* This function uses dynamic programming to create an optimal matching on the given path and adds this matching to the given list of contractions. */
    /* Based in Fig. 7 from Maue and Sanders (2007). */
    double weights[2];
    long sizes[2];
    long t;
    int cur = 1;
    long offset = 0;
    
    /* M[0] = empty. */
    weights[0] = 0.0;
    sizes[0] = 0;
    
    /* M[1] = first edge. */
    weights[1] = PathWeights[0];
    sizes[1] = 1;
    matchings[1][0] = 0;
    
    for (t = 1; t < PathLength - 1; t++) {
        if (PathWeights[t] + weights[1 - cur] > weights[cur]) {
            weights[1 - cur] += PathWeights[t];
            matchings[1 - cur][sizes[1 - cur]++] = t;
        } else {
            weights[1 - cur] = weights[cur];
            sizes[1 - cur] = sizes[cur];
            
            for ( ; offset < sizes[cur]; offset++) {
                matchings[2][offset] = matchings[cur][offset];
            }
        }
        
        cur = 1 - cur;
    }
    
    for ( ; offset < sizes[cur]; offset++) {
        matchings[2][offset] = matchings[cur][offset];
    }

#ifdef INFO
    /* Sanity check. */
    if (TRUE) {
        double oddeven[2] = {0.0, 0.0};
        double opt = 0.0;
        
        for (t = 0; t < PathLength - 1; t++) {
            oddeven[t & 1] += PathWeights[t];
        }
        
        for (t = 0; t < offset; t++) {
            opt += PathWeights[matchings[2][t]];
        }
        
        /* Matching should be at least as good as all even or all odd edges. */
        if (oddeven[0] > opt || oddeven[1] > opt || weights[1 - cur] > opt || weights[cur] > opt) {
            fprintf(stderr, "FindOptimalPathMatching(): Weight sanity check failed!\n");
            return -1;
        }
    }
#endif
    
    return offset;
}

int ApplyPathMatching(struct biparthypergraph *pHG, struct contraction *pC, int *Matched,
                      const long *PathMatching, const long *PathIndices, const long PathMatchingLength)
{
    /* Applies a path matching generated by FindOptimalPathMatching. */
    
    /* Create groups of two vertices based on this matching. */
    long t;
    
    for (t = 0; t < PathMatchingLength; t++) {
        const long v = PathIndices[PathMatching[t]];
        const long w = PathIndices[PathMatching[t] + 1];
        const long k = pC->NrMatches;
        
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
    long *PathMatchings[3];
    /* Array containing the matching weights along the constructed path. */
    double *PathWeights = (double *)malloc(pHG->NrVertices*sizeof(double));
    /* Array containing the indices of the vertices along this path. */
    long *PathIndices = (long *)malloc(pHG->NrVertices*sizeof(long));
    long PathLength = 0;
    
    /* Data for neighbor finding. */
    void *pData = NULL;
    
    long t, tt;
    
    PathMatchings[0] = (long *)malloc(pHG->NrVertices*sizeof(long));
    PathMatchings[1] = (long *)malloc(pHG->NrVertices*sizeof(long));
    PathMatchings[2] = (long *)malloc(pHG->NrVertices*sizeof(long));
    
    if (!pHG || !pC || !iv || !Matched || !pOptions ||
        !Ip || !ScIp || !Visited ||
        !PathMatchings[0] || !PathMatchings[1] || !PathMatchings[2] ||
        !PathWeights || !PathIndices) {
        fprintf(stderr, "MatchUsingPGA(): Null arguments or not enough memory!\n");
        return FALSE;
    }
    
    if (pOptions->Coarsening_MaxNrVtxInMatch != 2) {
        fprintf(stderr, "MatchUsingPGA(): This only works when matching groups of two vertices!\n");
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
            
            CurrentVertex = NextVertex;
            NextVertex = -1;
        }
        
        /* Reset matched flags. */
        for (tt = 0; tt < PathLength; tt++) {
            Matched[PathIndices[tt]] = FALSE;
        }
        
        /* Only match if the path contains at least two vertices. */
        if (PathLength <= 1) {
            continue;
        }
        
        /* Find optimal weight matching along the path using dynamic programming. */
        PathMatchingLength = FindOptimalPathMatching(PathMatchings, PathWeights, PathLength);
        
        if (PathMatchingLength < 0) {
            fprintf(stderr, "MatchUsingPGA(): Unable to find a matching along the path!\n");
            return FALSE;
        }
        
        if (!ApplyPathMatching(pHG, pC, Matched, PathMatchings[2], PathIndices, PathMatchingLength)) {
            fprintf(stderr, "MatchUsingPGA(): Unable to apply path matching!\n");
            return FALSE;
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
    free(PathMatchings[0]);
    free(PathMatchings[1]);
    free(PathMatchings[2]);
    
    FreeData(pData);
    
    return TRUE;
} /* end MatchUsingPGA */


