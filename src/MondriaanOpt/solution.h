#ifndef SOLUTION_H_DONE

#define SOLUTION_H_DONE

#include "options.h"
#include "matrix.h"

#include <math.h>
#include <time.h>

/* Define values (categories) of rows and columns */
#define TO0 0  /* Assigned to part 0 */
#define TO1 1  /* Assigned to part 1 */
#define CUT 2  /* Explicitly or implicitly cut */
#define PAR0 3 /* Partially assigned to part 0 */
#define PAR1 4 /* Partially assigned to part 1 */
#define UNAS 5 /* Unassigned */
#define ASGN 6 /* Assigned */

#define TRUE 1
#define FALSE 0

struct indexedint{
    unsigned int val;
    unsigned int key;
};

struct priority{
    unsigned int idx;
    struct priority *next;
};

/* Bipartite graph for matching */
struct match_graph {
    int *match; /* match[i] = j and match[j]=i if vertex i is matches to j,
                   otherwise match[i]=i */
    int nmatch; /* number of matches */
    int *adj;   /* stores the adjacency lists of the vertices */
    int *start; /* adjacency list of vertex i is stored in adj[start[i]..start[i]+length[i]-1] */
    int *length;
    int *visited; /* visited[i] = TRUE if vertex i has been visited in the breadth-first search (BFS) */
    int *vertex1; /* visited vertices in subset V1 */
    int *pred;    /* predecessors of vertices in V1 */
    int *vertex2; /* visited vertices in subset V2 */
    int *dstart;  /* vertices at distance 2*d edge steps from the initial vertex are stored in
                     vertex1[dstart[d]..dstart[d+1]-1] */
    unsigned char *inGraph; /* inGraph[ii] = TRUE if row/column ii is in the bipartite graph,
                               with 0 <= ii < m+n */
};


struct solution {
    /* Structure that holds a solution of the partitioning problem */

    /* Current row and column values */
    unsigned char* vals[2];

    unsigned int *n[2][2];

    unsigned int *startUnAssigned[2];

    unsigned char *partialstate[2];

    unsigned int *justpartial;
    unsigned int njustpartial;
    unsigned int *justpartialoff;
    unsigned int njustpartialoff;

    unsigned int *partial[2];
    unsigned int npartial[2];

    unsigned int *prior;
    unsigned int nprior;

    struct priority *nexts[2];
    struct indexedint *idxi[2];

    /* Bipartite graph for conflict nonzeros */
    struct match_graph mGraph;

    /* tot[0] = total number of nonzeros assigned to part 0 */
    unsigned int tot[2];

    /* Current depth in the tree */
    unsigned int depth;

    /* Number of tree branches searched */
    unsigned int nbranches;

    /* Previous direction */
    unsigned char prevdir;

    /* Boolean stating that a row/column index
       has already been assigned to part 0 */
    unsigned char FirstIdxAssignedTo0;

    /* Maximum allowed number of nonzeros for a part  */
    unsigned int max1;

    /* Current communication volume */
    unsigned int vol;

    /* Maximum allowed volume (= maximum allowed k) */
    unsigned int maxvol;


    /* Matrix name */
    char *matname;

    /* Start time of solution search */
    time_t starttime;

    /* Execution is stopped if maximum runtime (in s) has been exceeded */
    double maxruntime;

};

void insert_priority(struct solution *s, unsigned char dir, unsigned int n, unsigned int idx);
void delete_priority(struct solution *s, unsigned char dir, unsigned int n, unsigned int idx);

/* Initialise solution and allocate arrays for use with matrix */
void initsolution(const struct mat *a, struct solution *s, struct options *o);

/* Reset values of solution to zero */
void resetsolution(const struct mat *a, struct solution *s);

/* Function for printing a picture in Scalable Vector Graphics (SVG) format */
void printToSVG(const char *fn, struct solution *s, const struct mat *a);

/* Function for output of a solution as an integer sparse matrix in Matrix Market (MM) format */
void printToMM(const char *fn, struct solution *s, const struct mat *a);

void printConverted(struct options *o);
void fillFreeNonzeros(struct sparsematrix *pA);

#endif
