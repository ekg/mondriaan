#include <stdlib.h>
#include <stdio.h>

#include "Graph.h"
#include "HKLFM.h"
#include "Match.h"

struct opts Options;

extern int MoveVertex(struct biparthypergraph *pHG, long v);

int main(int argc, char **argv) {

    long n, v;
    struct biparthypergraph HG;
    struct bucketentry be;

    printf("Test MoveVertex: ");

    n = 16; /* must be even */
    HG.NrVertices = n;
    HG.WeightP[0] = 5*n;
    HG.WeightP[1] = n;
    HG.V = (struct vertex *) malloc(HG.NrVertices * sizeof(struct vertex));
    if (HG.V == NULL) {
        fprintf(stderr, "test_MoveVertex(): Not enough memory!\n");
        printf("Error\n");
        exit(1);
    }

    /* Initialise even vertices */
    for (v = 0; v < HG.NrVertices; v +=2) {
        HG.V[v].partition = 0;
        HG.V[v].vtxwgt = 10;
        HG.V[v].GBentry = &be;
    }

    /* Initialise odd vertices */
    for (v = 1; v < HG.NrVertices; v +=2) {
        HG.V[v].partition = 1;
        HG.V[v].vtxwgt = 2;
        HG.V[v].GBentry = &be;
    }

    /* Move all vertices */ 
    for (v = 0; v < HG.NrVertices; v++)
        MoveVertex(&HG, v);

    /* Check total weights */
    if (HG.NrVertices != n || HG.WeightP[0] != n || HG.WeightP[1] != 5*n) {
        printf("Error\n");
        exit(1);
    }

    /* Check even vertices */
    for (v = 0; v < HG.NrVertices; v +=2) {
        if (HG.V[v].partition != 1 || HG.V[v].vtxwgt != 10 ||
            HG.V[v].GBentry != NULL) { 
            printf("Error\n");
            exit(1);
        }
    }
    /* Check odd vertices */
    for (v = 1; v < HG.NrVertices; v +=2) {
        if (HG.V[v].partition != 0 || HG.V[v].vtxwgt != 2 ||
            HG.V[v].GBentry != NULL) {
            printf("Error\n");
            exit(1);
        }
    }

    printf("OK\n");
    exit(0);

} /* end main */
