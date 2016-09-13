#include "Options.h"
#include "MatchMatchers.h"

#define LENGTH 4

int main(int argc, char **argv) {
    long *PathMatchings[3];
    double *PathWeights;
    long Size;
    long t;
    
    for (t = 0; t < 3; t++) PathMatchings[t] = (long *)malloc(LENGTH*sizeof(long));
    PathWeights = (double *)malloc(LENGTH*sizeof(double));

    printf("Test FindOptimalPathMatching: ");
    
    /* Solution *=*-*=*-* */
    for (t = 0; t < LENGTH; t++) PathMatchings[2][t] = -1;
    PathWeights[0] = 2.0;
    PathWeights[1] = 1.9;
    PathWeights[2] = 0.5;
    PathWeights[3] = 0.4;
    
    Size = FindOptimalPathMatching(PathMatchings, PathWeights, LENGTH + 1);
    
    if (Size != 2 || PathMatchings[2][0] != 0 || PathMatchings[2][1] != 2) {
        printf("Error\n");
        exit(1);
    }

    /* Solution *-*=*-*=* */
    for (t = 0; t < LENGTH; t++) PathMatchings[2][t] = -1;
    PathWeights[0] = 1.9;
    PathWeights[1] = 2.0;
    PathWeights[2] = 0.5;
    PathWeights[3] = 0.5;
    
    Size = FindOptimalPathMatching(PathMatchings, PathWeights, LENGTH + 1);
    
    if (Size != 2 || PathMatchings[2][0] != 1 || PathMatchings[2][1] != 3) {
        printf("Error\n");
        exit(1);
    }

    /* Solution *=*-*-*=* */
    for (t = 0; t < LENGTH; t++) PathMatchings[2][t] = -1;
    PathWeights[0] = 2.0;
    PathWeights[1] = 1.9;
    PathWeights[2] = 0.5;
    PathWeights[3] = 0.6;
    
    Size = FindOptimalPathMatching(PathMatchings, PathWeights, LENGTH + 1);
    
    if (Size != 2 || PathMatchings[2][0] != 0 || PathMatchings[2][1] != 3) {
        printf("Error\n");
        exit(1);
    }

    printf("OK\n");
    exit(0);

} /* end main */

