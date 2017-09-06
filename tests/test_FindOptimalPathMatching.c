#include "Options.h"
#include "MatchMatchers.h"

#define LENGTH 4


int main(int argc, char **argv) {
    long Size;
    long t;
    
    /* Working arrays for creating optimal path matchings. */
    char *_PathMatchings = NULL, **PathMatchings = NULL, *PathMatchingOpt = NULL;
    double *_PathMatchingWeights = NULL, *PathMatchingWeights[2];
    /* Array containing the matching weights along the constructed path. */
    double *PathWeights = (double *)malloc(LENGTH*sizeof(double));
    /* Array containing the indices of the vertices along this path. */
    long *PathIndices = (long *)malloc(LENGTH*sizeof(long));
    long MaxNrVtxInMatch = 2;
    
    /* Set up PathMatchings */
    _PathMatchings = (char *)malloc(MaxNrVtxInMatch * LENGTH * sizeof(char));
    PathMatchings = (char **)malloc(LENGTH * sizeof(char *));
    for(t=0; t<LENGTH; ++t) {
        PathMatchings[t] = &(_PathMatchings[t*MaxNrVtxInMatch]);
    }
    /* Set up PathMatchingWeights */
    _PathMatchingWeights = (double *)malloc(2*MaxNrVtxInMatch * sizeof(double));
    PathMatchingWeights[0] = &(_PathMatchingWeights[0]);
    PathMatchingWeights[1] = &(_PathMatchingWeights[MaxNrVtxInMatch]);
    /* Set up PathMatchingOpt */
    PathMatchingOpt = (char *)malloc(LENGTH * sizeof(char));
    
    
    printf("Test FindOptimalPathMatching: ");
    
    /* Solution *=*-*=*-* */
    for (t = 0; t < LENGTH; t++) PathMatchingOpt[t] = PathMatchings[t][0] = PathMatchings[t][1] = -1;
    PathWeights[0] = 2.0;
    PathWeights[1] = 1.9;
    PathWeights[2] = 0.5;
    PathWeights[3] = 0.4;
    
    Size = FindOptimalPathMatching(PathMatchings, PathMatchingWeights, PathMatchingOpt, PathWeights, LENGTH + 1, MaxNrVtxInMatch);
    
    if (Size != 2 || PathMatchingOpt[0] != 1 || PathMatchingOpt[1] != 0 || PathMatchingOpt[2] != 1 || PathMatchingOpt[3] != 0) {
        printf("Error\n");
        exit(1);
    }

    /* Solution *-*=*-*=* */
    for (t = 0; t < LENGTH; t++) PathMatchingOpt[t] = PathMatchings[t][0] = PathMatchings[t][1] = -1;
    PathWeights[0] = 1.9;
    PathWeights[1] = 2.0;
    PathWeights[2] = 0.5;
    PathWeights[3] = 0.5;
    
    Size = FindOptimalPathMatching(PathMatchings, PathMatchingWeights, PathMatchingOpt, PathWeights, LENGTH + 1, MaxNrVtxInMatch);
    
    if (Size != 2 || PathMatchingOpt[0] != 0 || PathMatchingOpt[1] != 1 || PathMatchingOpt[2] != 0 || PathMatchingOpt[3] != 1) {
        printf("Error\n");
        exit(1);
    }

    /* Solution *=*-*-*=* */
    for (t = 0; t < LENGTH; t++) PathMatchingOpt[t] = PathMatchings[t][0] = PathMatchings[t][1] = -1;
    PathWeights[0] = 2.0;
    PathWeights[1] = 1.9;
    PathWeights[2] = 0.5;
    PathWeights[3] = 0.6;
    
    Size = FindOptimalPathMatching(PathMatchings, PathMatchingWeights, PathMatchingOpt, PathWeights, LENGTH + 1, MaxNrVtxInMatch);
    
    if (Size != 2 || PathMatchingOpt[0] != 1 || PathMatchingOpt[1] != 0 || PathMatchingOpt[2] != 0 || PathMatchingOpt[3] != 1) {
        printf("Error\n");
        exit(1);
    }
    
    
    free(PathWeights);
    free(PathIndices);
    free(_PathMatchings);
    free(PathMatchings);
    free(_PathMatchingWeights);
    free(PathMatchingOpt);

    printf("OK\n");
    exit(0);

} /* end main */

