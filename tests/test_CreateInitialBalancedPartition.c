#include <stdlib.h>
#include <stdio.h>

#include "Graph.h"
#include "HKLFM.h"
#include "Match.h"

struct opts Options;

extern int CreateInitialBalancedPartition(struct biparthypergraph *pHG, long weightlo, long weighthi);

int main(int argc, char **argv) {

    struct biparthypergraph HG;

    long n, t, totwgt, weightlo, weighthi, sum0, sum1;

    printf("Test CreateInitialBalancedPartition: ");
    n= 40; /* must be a multiple of 4 */

    HG.NrVertices = 2*n;

    HG.V = (struct vertex *) malloc(HG.NrVertices * sizeof(struct vertex));
    if (HG.V == NULL) {
        fprintf(stderr, "test_CreateInitialBalancedPartition(): Not enough memory!\n");
        printf("Error\n");
        exit(1);
    }

    /* Initialise vertex weights: 1,1,2,2,3,3,..., n,n */
    for (t=0; t<HG.NrVertices; t++)
        HG.V[t].vtxwgt = t/2 + 1;
    totwgt = n*(n+1);
    weightlo = totwgt/4 + n; 
    weighthi = 3*totwgt/4 + n; 

    if (!CreateInitialBalancedPartition(&HG, weightlo, weighthi)) {
        printf("Error\n");
        exit(1);
    }

    /* Check hypergraph dimensions */
    if (HG.NrVertices != 2*n) {
        printf("Error\n");
        exit(1);
    }

    /* Check vertex weights and partitions */
    sum0 = 0;
    sum1 = 0;
    for (t=0; t<HG.NrVertices; t++) {
        if (HG.V[t].vtxwgt != t/2 + 1 ||
            HG.V[t].partition < 0 ||
            HG.V[t].partition > 1) {

            printf("Error\n");
            exit(1);
        }
       if (HG.V[t].partition == 0)
           sum0 += HG.V[t].vtxwgt; 
       else if (HG.V[t].partition == 1)
           sum1 += HG.V[t].vtxwgt; 
    }

    /* Check part weights */
    if (sum0 > weightlo || sum1 > weighthi) {
        printf("Error\n");
        exit(1);
    }

    printf("OK\n");
    exit(0);

} /* end main */
