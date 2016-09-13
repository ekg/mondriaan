#include "matching.h"
#include "matrix.h"

void update_graph_after_setting(struct mat *a, struct solution *sol, int i,
                                unsigned char dir, unsigned char part){

    /* This function updates the graph after setting row i to part 
       if dir = ROW, and setting column i if dir = COL.
       Comments below are for the case dir = ROW. */

    int k, ki;

    remove_vertex_from_graph(a,sol,i,dir);

    for(ki=0; ki<sol->njustpartial; ki++){
        /* add column k that just became partially assigned 
           as a column vertex to the graph */
        k = sol->justpartial[ki];
        add_vertex_to_graph(a,sol,k,1-dir,part);
    }

    for(ki=0; ki<sol->njustpartialoff; ki++){
        /* remove column k that just became not partially assigned */
        k = sol->justpartialoff[ki];
        remove_vertex_from_graph(a,sol,k,1-dir);
    }

} /* end update_graph_after_setting */


void update_graph_after_unsetting(struct mat *a, struct solution *sol, int i,
                                  unsigned char dir, unsigned char part){

    /* This function updates the graph after unsetting row i from part
       if dir = ROW, and unsetting column i if dir = COL.
       Comments below are for the case dir = ROW. */

    int k, ki;

    for(ki=0; ki<sol->njustpartial; ki++){
        /* add column k that just became partially assigned
           to the other part (1-part) as a column vertex to the graph */
        k = sol->justpartial[ki];
        add_vertex_to_graph(a,sol,k,1-dir,1-part);
    }

    for(ki=0; ki<sol->njustpartialoff; ki++){
        /* remove column k that just became not partially assigned */
        k = sol->justpartialoff[ki];
        remove_vertex_from_graph(a,sol,k,1-dir);
    }

    /* n[part][dir][i] = the number of nonzeros in row i and in a column partially
                         assigned to part */
    if(sol->n[TO0][dir][i]>0 && sol->n[TO1][dir][i]==0){
        add_vertex_to_graph(a,sol,i,dir,TO0);
    } else if(sol->n[TO1][dir][i]>0 && sol->n[TO0][dir][i]==0){
        add_vertex_to_graph(a,sol,i,dir,TO1);
    }

} /* end update_graph_after_unsetting */


void add_vertex_to_graph(struct mat *a, struct solution *sol, int i,
                         unsigned char dir, unsigned char part){

    /* This function adds the vertex corresponding to a partially assigned
       row i to the graph if dir = ROW (and column i if dir = COL).
       The row is partially assigned to part (either 0 or 1).
       Comments below are for the case dir = ROW. */


    int j, k, ii, ki;

    /* Translate (i,dir) to a single vertex number ii, 0 <= ii < m+n */
    ii = i+dir*a->m; 

    /* Add vertex ii, even if it will have no edges */
    sol->mGraph.inGraph[ii] = TRUE;

    for(j=a->starts[dir][i]; j<a->starts[dir][i+1]; j++){
        k = a->indices[dir][j];
        if(sol->partialstate[1-dir][k] == PAR0+(1-part)){
            /* column k is partially assigned to the other part, so we have a conflict */
            ki = k+(1-dir)*a->m;
            if(sol->mGraph.inGraph[ki]){
                /* both row i and column k are in the graph, so we add edge (ii,ki)
                   to the adjacency lists of vertices ii and ki */
                sol->mGraph.adj[sol->mGraph.start[ii]+sol->mGraph.length[ii]] = ki;
                sol->mGraph.length[ii]++;
                sol->mGraph.adj[sol->mGraph.start[ki]+sol->mGraph.length[ki]] = ii;
                sol->mGraph.length[ki]++;
            }
        }
    }

    /* Restore an optimal matching by searching for an augmenting path starting in vertex ii */
    AugmentPath(ii, sol->mGraph.match, &(sol->mGraph.nmatch),
                sol->mGraph.adj, sol->mGraph.start, sol->mGraph.length,
                sol->mGraph.visited, sol->mGraph.vertex1, sol->mGraph.vertex2,
                sol->mGraph.pred, sol->mGraph.dstart);

} /* end add_vertex_to_graph */


void remove_vertex_from_graph(struct mat *a, struct solution *sol, int i,
                              unsigned char dir){

    /* This function removes the vertex corresponding to row i and all its edges
       from the graph if dir = ROW (and column i if dir = COL).
       Comments below are for the case dir = ROW. */


    int j, l, ii, ki;

    ii = i+dir*a->m;

    sol->mGraph.inGraph[ii] = FALSE;

    for(j=sol->mGraph.start[ii]; j<sol->mGraph.start[ii]+sol->mGraph.length[ii]; j++){
        ki = sol->mGraph.adj[j];
        /* Search for row i in column k, i.e. index ii in the adjacency list of ki */
        for(l=sol->mGraph.start[ki]; l<sol->mGraph.start[ki]+sol->mGraph.length[ki]; l++){
            if(sol->mGraph.adj[l]==ii){
                /* Found! Delete edge (ii,ki) from the adjacency list
                   by copying the last adjacency into its position */
                sol->mGraph.length[ki]--;
                sol->mGraph.adj[l] = sol->mGraph.adj[sol->mGraph.start[ki]+sol->mGraph.length[ki]];
                break;
            }
        }
    }
    /* Empty the adjacency list of ii */
    sol->mGraph.length[ii] = 0;

    if(sol->mGraph.match[ii]!=ii){
        /* removed vertex ii was matched to vertex ki */
        ki = sol->mGraph.match[ii];

        /* Break the match */
        sol->mGraph.match[ki]=ki;
        sol->mGraph.match[ii]=ii;
        sol->mGraph.nmatch--;

        /* Restore an optimal matching by searching for an augmenting path starting in 
           the former match ki of vertex ii */
        AugmentPath(ki, sol->mGraph.match, &(sol->mGraph.nmatch),
                    sol->mGraph.adj, sol->mGraph.start, sol->mGraph.length,
                    sol->mGraph.visited, sol->mGraph.vertex1, sol->mGraph.vertex2,
                    sol->mGraph.pred, sol->mGraph.dstart);
    }

} /* end remove_vertex_from_graph */


void AugmentPath (int v, int *match, int *nmatch,
                  int *adj, int *start, int *length,
                  int *visited, int *vertex1, int *vertex2,
                  int *pred, int *dstart) {

/*  This function finds an augmenting path in the current matching M of a bipartite graph
    starting from vertex v, provided such a path exists, and then modifies the matching M,
    increasing the number of matches by 1. If no augmenting path exists,
    the matching remains unchanged. This is done by Breadth-First Search (BFS).

    v =     a vertex, corresponding to a matrix row (if 0 <= v < m)
            or a matrix column (if m <= v < m+n)
            v is in V1, where the bipartite graph has vertices in V1 and V2.
    match = array of length m+n storing the matches:
            match[i]=j and match[j]=i if i and j are matched, i.e. (i,j) in M
            match[i]=i if vertex i is unmatched
    nmatch = number of matches

    adj = array of size 2*nnz storing the adjacencies of the vertices in CRS format;
          those for vertex i are stored in adj[start[i].. start[i]+length[i]-1]

    visited = boolean array of lengthm m+n, FALSE on input and output,
            used to see quickly whether a vertex in V2
            has already been visited by the BFS

    WORKING SPACE (allocated once outside this function):
    vertex1 = array of length m+n, to be used as working space for storing
             the vertices in V1 encountered. The vertices at distance d from v
             are stored in vertex1[dstart[d]..dstart[d+1]-1].
    pred = predecessors of vertices in V1 stored in the same way as vertex1.
    dstart = array of length m+n+2 (to be safe).
    vertex2 = array of length m+n, to be used as working space for storing
             the vertices in V2 encountered. */


    int nvisited = 0 ; /* number of vertices in V2 visited so far in the BFS */
    int done = FALSE ; /* boolean denoting whether the search is finished */
    int dist ;         /* distance from v in steps of size 2 (back and forth) */
    int u ;            /* current vertex in V1 */
    int w ;            /* a neighbour of u */
    int mu ;           /* the match of u */
    int mw ;           /* the match of w */
    int j, k ;

    /* Initialise vertex v with distance 0 */
    dist = 0;
    dstart[0] = 0 ;
    dstart[1] = 1 ;
    vertex1[0] = v ;
    pred[0] = 0;

    while ((!done) && dstart[dist] != dstart[dist+1]){

       /* Handle all vertices at distance dist */

       dstart[dist+2] = dstart[dist+1] ;
       for (k = dstart[dist]; k < dstart[dist+1]; k++){

           u = vertex1[k] ; /* in V1 */

           for (j = start[u]; j < start[u] + length[u]; j++){
               w = adj[j] ; /* in V2 */

               if (visited[w])
                   continue ;
               visited[w] = TRUE ;
               vertex2[nvisited] = w ;
               nvisited++ ;

               if (match[w] == w){
                  /* w is unmatched so we have found an augmenting path from v to w */

                  /* Follow the links back to v and flip the matches */
                  match[w] = u ;
                  while (u != v) {
                      mu = match[u] ;
                      match[u] = w ;
                      k = pred[k] ;
                      u = vertex1[k] ;
                      w = mu ;
                      match[w] = u ;
                  }
                  match[v] = w ;
                  (*nmatch)++ ;
                  done = TRUE ;
                  break ;
               } else {
                   mw = match[w] ; /* in V1 */
                   vertex1[dstart[dist+2]] = mw ;
                   pred[dstart[dist+2]] = k ; /* predecessor of mw in V1 is u */
                   dstart[dist+2]++ ;
               }
            }

            if (done)
                break ;
        }
        dist++ ;
    }

    /* Clean up */
    for (k = 0; k < nvisited  ; k++)
        visited[vertex2[k]] = FALSE ;

    return ;

} /* end AugmentPath */
