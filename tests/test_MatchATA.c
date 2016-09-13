#include "Options.h"
#include "Sort.h"
#include "SparseMatrix.h"
#include "Graph.h"
#include "MatchMatchers.h"
#include "MatchInproduct.h"
#include "MatchStairway.h"

long testATAMatcher(
        const struct opts *pOptions,
        const struct biparthypergraph *pHGOrig,
        int (*HybridMatcher)(struct biparthypergraph *, struct contraction *,
                             const long *, int *,
                             const struct opts *,
                             int (*SetupData)(void **, struct biparthypergraph *, const struct opts *),
                             int (*FindNeighbor)(long *, double *,
                                                 const struct biparthypergraph *, const struct contraction *, const struct opts *,
                                                 const long, const int *,
                                                 void *, long *, long *, double *),
                             int (*FreeData)(void *)),
        int (*SetupData)(void **, struct biparthypergraph *, const struct opts *),
        int (*FindNeighbor)(long *, double *,
                            const struct biparthypergraph *, const struct contraction *, const struct opts *,
                            const long, const int *,
                            void *, long *, long *, double *),
        int (*FreeData)(void *))
{
    struct biparthypergraph G;
    struct contraction C;
    long *VtxDegree;
    long *iv;
    int *Matched;
    int *Flags;
    long MatchingWeight;
    long t, tt;
    
    /* Create a copy of the supplied hypergraph (SetupData() may very well change it). */
    if (!CreateNewBiPartHyperGraph(pHGOrig->NrVertices,
                                   pHGOrig->NrNets,
                                   pHGOrig->NrPins,
                                   FALSE,
                                   'P',
                                   &G)) {
        return -1;
    }
    
    G.SplitDir = pHGOrig->SplitDir;
    memcpy(G.V, pHGOrig->V, G.NrVertices*sizeof(struct vertex));
    memcpy(G.VtxAdjncy, pHGOrig->VtxAdjncy, G.NrPins*sizeof(long));
    memcpy(G.N, pHGOrig->N, G.NrNets*sizeof(struct net));
    memcpy(G.NetAdjncy, pHGOrig->NetAdjncy, G.NrPins*sizeof(long));
    
    for (t = 0; t < G.NrNets; t++)
         G.N[t].iStartP1 = G.N[t].iEnd; 
    
    /* Order vertices by decreasing number of nonzeros. */
    VtxDegree = (long *)malloc(G.NrVertices*sizeof(long));
    
    for (t = 0; t < G.NrVertices; t++) 
        VtxDegree[t] = G.V[t].iEnd - G.V[t].iStart;
    
    iv = QSort(VtxDegree, G.NrVertices);
    free(VtxDegree);
    
    /* Create Matched vertex array. */
    Matched = (int *)malloc(G.NrVertices*sizeof(int));
    
    for (t = 0; t < G.NrVertices; t++)
        Matched[t] = FALSE;
    
    /* Create contraction structure. */
    C.Match = (long *)malloc(G.NrVertices*sizeof(long));
    C.Start = (long *)malloc((G.NrVertices + 1)*sizeof(long));
    
    C.Start[0] = 0;
    C.NrMatches = 0;
    C.MaxNrVertices = 2;
    C.MaxVtxWgt = 2*G.NrPins;
    
    /* Execute matching algorithm. */
    if (!HybridMatcher(&G, &C,
                       iv, Matched,
                       pOptions,
                       SetupData,
                       FindNeighbor,
                       FreeData)) {
        fprintf(stderr, "Unable to create a hybrid matching!\n");
        return -1;
    }
    
    /* Verify matching and calculate its weight. */
    Flags = (int *)malloc(G.NrNets*sizeof(int));
    MatchingWeight = 0;
        
    if (C.Start[0] != 0) {
        fprintf(stderr, "Invalid starting offset in C!\n");
        return -1;
    }
    
    for (t = 0; t < G.NrNets; t++)
        Flags[t] = 0;
    
    for (t = 0; t < C.NrMatches; t++)
    {
        if (C.Start[t + 1] == C.Start[t] + 2) {
            const long v1 = C.Match[C.Start[t]];
            const long v2 = C.Match[C.Start[t] + 1];
            
            if (v1 < 0 || v2 < 0 || v1 >= G.NrVertices || v2 >= G.NrVertices) {
                fprintf(stderr, "Invalid matched vertex indices!\n");
                return -1;
            }
            
            /* Calculate inner product. */
            for (tt = G.V[v1].iStart; tt < G.V[v1].iEnd; tt++)
                Flags[G.VtxAdjncy[tt]] = 1;
    
            for (tt = G.V[v2].iStart; tt < G.V[v2].iEnd; tt++)
                if (Flags[G.VtxAdjncy[tt]] == 1) MatchingWeight++;
    
            for (tt = G.V[v1].iStart; tt < G.V[v1].iEnd; tt++)
                Flags[G.VtxAdjncy[tt]] = 0;
        }
        else if (C.Start[t + 1] != C.Start[t] + 1) {
            fprintf(stderr, "Erroneous matching group size!\n");
            return -1;
        }
    }
    
#ifdef INFO
    fprintf(stderr, "Total matching weight equals %ld.\n", MatchingWeight);
#endif
    
    /* Free data. */
    DeleteBiPartHyperGraph(&G);
    free(C.Match);
    free(C.Start);
    free(iv);
    free(Matched);
    free(Flags);
    
    return MatchingWeight;
}

int main(int argc, char **argv) {
    struct opts Options;
    struct sparsematrix A;
    struct biparthypergraph G;
    long volume = 0;
    FILE *file;
    
    printf("Test MatchATA: ");
    
    /* Read matrix and convert it to a hypergraph. */
    if (!SetDefaultOptions(&Options)) {
        printf("Error\n");
        exit(1);
    }
    
    Options.SplitStrategy = OneDimCol;
    Options.Coarsening_NetScaling = NoNetScaling;
    Options.Coarsening_InprodScaling = NoIpScaling;
    
    if (!ApplyOptions(&Options)) {
        printf("Error\n");
        exit(1);
    }
    
    srand(123456);
    
    if (!(file = fopen(argc != 2 ? "n4c6-b2.mtx" : argv[1], "r"))) {
        printf("Error\n");
        exit(1);
    }
    
    if (!MMReadSparseMatrix(file, &A)) {
        printf("Error\n");
        exit(1);
    }
    
    fclose(file);
    
    if (A.m == A.n && (A.MMTypeCode[3] == 'S' || A.MMTypeCode[3] == 'K' || A.MMTypeCode[3] == 'H')) {
        if (!SparseMatrixSymmetric2Full(&A)) {
            printf("Error\n");
            exit(1);
        }
    }
    
    if (!SparseMatrix2BiPartHyperGraph(&A, COL, &Options, &G)) {
        printf("Error\n");
        exit(1);
    }
    
    /* Test the different hybrid matchers. */
    /* For n4c6-b2, we should have (on average):
        Greedy-Stairway     100
        Greedy-Inproduct    100
        PGA-Stairway        105
        PGA-Inproduct       105
       For dfl001.mtx, we should have (on average):
        Greedy-Stairway     7753
        Greedy-Inproduct    7824
        PGA-Stairway        8011
        PGA-Inproduct       8050
     */
    volume = testATAMatcher(&Options, &G, MatchUsingGreedy, MatchStairwaySetup, FindNeighborStairway, MatchStairwayFree);
    
    if (volume < 0 || (argc != 2 && volume < 95)) {
        printf("Error\n");
        exit(1);
    }
    
    volume = testATAMatcher(&Options, &G, MatchUsingGreedy, MatchInproductSetup, FindNeighborInproduct, MatchInproductFree);
    
    if (volume < 0 || (argc != 2 && volume < 95)) {
        printf("Error\n");
        exit(1);
    }
    
    volume = testATAMatcher(&Options, &G, MatchUsingPGA, MatchStairwaySetup, FindNeighborStairway, MatchStairwayFree);
    
    if (volume < 0 || (argc != 2 && volume != 105)) {
        printf("Error\n");
        exit(1);
    }
    
    volume = testATAMatcher(&Options, &G, MatchUsingPGA, MatchInproductSetup, FindNeighborInproduct, MatchInproductFree);
    
    if (volume < 0 || (argc != 2 && volume != 105)) {
        printf("Error\n");
        exit(1);
    }
    
    printf("OK\n");
    exit(0);
} /* end main */

