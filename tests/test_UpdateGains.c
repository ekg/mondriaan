#include <stdlib.h>
#include <stdio.h>

#include "Graph.h"
#include "HKLFM.h"
#include "Match.h"

struct opts Options;

extern int UpdateGains(struct biparthypergraph *pHG, long v);

int main(int argc, char **argv) {

    struct biparthypergraph HG;

    long m, n, nz, t, P, v, v0, v1, v2, val0, val1;

    long vtxadj[27] = { 1, 2, 4, 5, 7, 8,   /* adjacency list of vertex 0 */
                        2, 5, 8,            /*                   vertex 1 */
                        0, 1, 2, 3, 4, 5, 6, 7, 8,
                        3, 4, 5, 6, 7, 8,  
                        6, 7, 8
                       };
    long netadj[27] = { 2,     /* adjacency list of net 0 */
                        0, 2,  /*                   net 1 */
                        0, 1, 2,
                        2, 3,
                        0, 2, 3,
                        0, 1, 2, 3, 
                        2, 3, 4,
                        0, 2, 3, 4,
                        0, 1, 2, 3, 4
                       };
    long gain[5] = { -2, -1, 0, 3, 0};

    printf("Test UpdateGains: ");
    m = 9; /* all parameters are fixed */
    n = 5;
    nz = 27;

    HG.NrNets = m;
    HG.NrVertices = n;
    HG.NrPins = nz;

    HG.V = (struct vertex *) malloc(HG.NrVertices * sizeof(struct vertex));
    HG.N = (struct net *) malloc(HG.NrNets * sizeof(struct net));
    HG.VtxAdjncy = (long *) malloc(HG.NrPins * sizeof(long));
    HG.NetAdjncy = (long *) malloc(HG.NrPins * sizeof(long));
    if (HG.V == NULL ||  HG.N == NULL ||
         HG.VtxAdjncy == NULL || HG.NetAdjncy == NULL) {
        fprintf(stderr, "test_UpdateGains(): Not enough memory!\n");
        printf("Error\n");
        exit(1);
    }

    /* Initialise partition of vertices */
    HG.V[0].partition = 0;
    HG.V[1].partition = 0;
    HG.V[2].partition = 1; /* has just been moved to part 1 */
    HG.V[3].partition = 1;
    HG.V[4].partition = 1;

    /* Initialise start and end of vertices */
    HG.V[0].iStart =  0; HG.V[0].iEnd =  6;
    HG.V[1].iStart =  6; HG.V[1].iEnd =  9;
    HG.V[2].iStart =  9; HG.V[2].iEnd = 18;
    HG.V[3].iStart = 18; HG.V[3].iEnd = 24;
    HG.V[4].iStart = 24; HG.V[4].iEnd = 27;

    /* Initialise adjacency lists of vertices */
    for (t=0; t<nz; t++) 
        HG.VtxAdjncy[t] = vtxadj[t]; 

    /* Initialise start and end of nets */
    HG.N[0].iStartP0 =  0; HG.N[0].iStartP1 =  1; HG.N[0].iEnd =  1;
    HG.N[1].iStartP0 =  1; HG.N[1].iStartP1 =  3; HG.N[1].iEnd =  3;
    HG.N[2].iStartP0 =  3; HG.N[2].iStartP1 =  6; HG.N[2].iEnd =  6;
    HG.N[3].iStartP0 =  6; HG.N[3].iStartP1 =  7; HG.N[3].iEnd =  8;
    HG.N[4].iStartP0 =  8; HG.N[4].iStartP1 = 10; HG.N[4].iEnd = 11;
    HG.N[5].iStartP0 = 11; HG.N[5].iStartP1 = 14; HG.N[5].iEnd = 15;
    HG.N[6].iStartP0 = 15; HG.N[6].iStartP1 = 16; HG.N[6].iEnd = 18;
    HG.N[7].iStartP0 = 18; HG.N[7].iStartP1 = 20; HG.N[7].iEnd = 22;
    HG.N[8].iStartP0 = 22; HG.N[8].iStartP1 = 25; HG.N[8].iEnd = 27;

    /* Initialise adjacency lists of nets */
    for (t=0; t<nz; t++) 
        HG.NetAdjncy[t] = netadj[t]; 

    /* Initialise gainbucket data structure */
    for (P=0; P<2; P++) {
        HG.GBVtx[P].NrBuckets  = 0;
        HG.GBVtx[P].Root = NULL;
    }
    
    /* Insert vertices 0, 1, 3, 4 into gainbucket data structure */
    for (v = 0; v < 2; v++)
        HG.V[v].GBentry = BucketInsert(&(HG.GBVtx[0]), gain[v], v);
    for (v = 3; v < 5; v++)
        HG.V[v].GBentry = BucketInsert(&(HG.GBVtx[1]), gain[v], v);
    HG.V[2].GBentry = NULL; /* vertex 2 is locked, and is excluded
                                  from the gainbucket data structure */
    HG.CurComm = 6;
    
    UpdateGains(&HG, 2);

    /* Check hypergraph dimensions */
    if (HG.NrNets != m || HG.NrVertices != n || HG.NrPins != nz) {
        printf("Error\n");
        exit(1);
    }
    
    /* Check start and end of nets and net parts */
    if(HG.N[0].iStartP0 !=  0 || HG.N[0].iStartP1 !=  0 || HG.N[0].iEnd !=  1 ||
        HG.N[1].iStartP0 !=  1 || HG.N[1].iStartP1 !=  2 || HG.N[1].iEnd !=  3 ||
        HG.N[2].iStartP0 !=  3 || HG.N[2].iStartP1 !=  5 || HG.N[2].iEnd !=  6 ||
        HG.N[3].iStartP0 !=  6 || HG.N[3].iStartP1 !=  6 || HG.N[3].iEnd !=  8 ||
        HG.N[4].iStartP0 !=  8 || HG.N[4].iStartP1 !=  9 || HG.N[4].iEnd != 11 ||
        HG.N[5].iStartP0 != 11 || HG.N[5].iStartP1 != 13 || HG.N[5].iEnd != 15 ||
        HG.N[6].iStartP0 != 15 || HG.N[6].iStartP1 != 15 || HG.N[6].iEnd != 18 ||
        HG.N[7].iStartP0 != 18 || HG.N[7].iStartP1 != 19 || HG.N[7].iEnd != 22 ||
        HG.N[8].iStartP0 != 22 || HG.N[8].iStartP1 != 24 || HG.N[8].iEnd != 27) {

        printf("Error\n");
        exit(1);
    }

    /* Check net parts with one adjacency (these are easy to check) */
    if (HG.NetAdjncy[HG.N[0].iStartP1] != 2 ||
        HG.NetAdjncy[HG.N[1].iStartP0] != 0 ||
        HG.NetAdjncy[HG.N[1].iStartP1] != 2 ||
        HG.NetAdjncy[HG.N[2].iStartP1] != 2 ||
        HG.NetAdjncy[HG.N[4].iStartP0] != 0 ||
        HG.NetAdjncy[HG.N[7].iStartP0] != 0) {

        printf("Error\n");
        exit(1);
    }

    /* Check gainbucket data structure for part 0 */
    v0 = GainBucketGetMaxValVertexNr(&(HG.GBVtx[0]));
    val0 = GainBucketGetMaxVal(&(HG.GBVtx[0]));
    BucketDeleteMax(&(HG.GBVtx[0]));
    v1 = GainBucketGetMaxValVertexNr(&(HG.GBVtx[0]));
    val1 = GainBucketGetMaxVal(&(HG.GBVtx[0]));
    BucketDeleteMax(&(HG.GBVtx[0]));
    v2 = GainBucketGetMaxValVertexNr(&(HG.GBVtx[0])); /* bucket should be empty */
    if (v0 != 0 || val0 != 3 ||  v1 != 1 || val1 != 0 || v2 != LONG_MIN) {
            printf("Error\n");
            exit(1);
    }

    /* Check gainbucket data structure for part 1 */
    v0 = GainBucketGetMaxValVertexNr(&(HG.GBVtx[1]));
    val0 = GainBucketGetMaxVal(&(HG.GBVtx[1]));
    BucketDeleteMax(&(HG.GBVtx[1]));
    v1 = GainBucketGetMaxValVertexNr(&(HG.GBVtx[1]));
    val1 = GainBucketGetMaxVal(&(HG.GBVtx[1]));
    BucketDeleteMax(&(HG.GBVtx[1]));
    v2 = GainBucketGetMaxValVertexNr(&(HG.GBVtx[1])); /* bucket should be empty */
    if (v0 != 4 || val0 != -1 ||  v1 != 3 || val1 != -2 || v2 != LONG_MIN) {
            printf("Error\n");
            exit(1);
    }

    /* Check current communication */
    if (HG.CurComm != 6) {
        printf("Error\n");
        exit(1);
    }

    printf("OK\n");
    exit(0);

} /* end main */
