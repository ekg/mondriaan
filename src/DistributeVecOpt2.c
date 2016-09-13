#include "DistributeVecOpt2.h"

int DistributeVecOpt2(const struct sparsematrix *pM, long int *X, int dir) {
    /* This function computes a vector distribution X with provably minimal
       cost for the distributed matrix A in the special case
       where all matrix columns have at most two processors. 
       The cost is defined as the maximum over all processors
       of max(sends,recvs).
       
       Input:  A distributed matrix.
               dir = ROW (for distribution of v) or dir = COL (for u).
       Output: X vector distribution, with 0 <= X[j] < P, for 0 <= j < pM->n,
                   where P= pM->NrProcs. X[j] is the owner of vector component j.
                   
       If the number of processors in each column of A is <= 2, 
       the function finds an optimal solution X and returns TRUE;
       otherwise the function does not do anything and returns FALSE.
            
       This function uses 3 arrays of length P*P. */
 
    long j, l=0, total;
    int P, q, r, tmp, count;

    /* Arrays of length pM->n */
    int *Nprocs;    /* Nprocs[j] is the number of processors that occur
                       in column j */
                        
    /* Communication matrix */
    int *procindex; /* array of length at most V+pM->n, where V is the
                       communication volume for the vector v.
                       procindex[procstart[j]..procstart[j+1]-1]
                       contains the processor numbers of the processors
                       that occur in column j */
    long *procstart; /* array of length pM->n+1 containing the starts of the
                       columns of the communication matrix */
                       
    /* Arrays of length P */
    long *Nv;        /* Nv[q] = current number of vector components (columns)
                                owned by processor q */
    int *degree;     /* degree[q] = current number of unowned columns 
                                    with Nproc=2 in which processor q occurs
                                  = current number of edges of vertex q in graph
                        0 <= degree[q] < P. */ 
    int *first_adj;  /* first_adj[q]  = adjacent vertex of vertex q with
                                        smallest vertex number,
                                        or P if no adjacent vertex exists.
                        0 <= first_adj[q] <= P. */

    /* P x P matrices */
    long **G;   /* Current adjacency matrix representing the connectivity graph:
                   G[q][r] = number of remaining edges connecting processors
                             q and r (corresponding to the unowned columns
                             with Nprocs = 2 in which processors q and r occur).
                   G is a symmetric matrix with zeros on the main diagonal. 
                   G changes because matrix elements decrease in value 
                   when edges are removed from the graph. */
    long **G2;  /* Reduced adjacency matrix at the start of phase 1
                   (when all columns with Nprocs = 2 are still unowned):
                        G2[q][r] = G[q][r] mod 2
                   G2 is a symmetric matrix with elements 0 and 1.
                   The diagonal is all zero. G2 does not change. */
    long **Col;  /* Matrix containing the column numbers of the unowned columns
                    after phase 1. 
                    Col[q][r] = the column number j of the remaining unowned
                                column shared by processors q and r,
                                if such a column exists, and -1 otherwise. 
                    Col does not change. */
    if (!pM || !X) {
        fprintf(stderr, "DistributeVecOpt2(): Null arguments!\n");
        return FALSE;
    }
   
    if (dir == COL) {
        l = pM->m;
    }
    else if (dir == ROW) {
        l = pM->n;
    }
    else {
        fprintf(stderr, "DistributeVecOpt2: Direction unknown!\n");
        return FALSE;
    }
    
    P= pM->NrProcs;

    /* ## Initialise columns ## */
    Nprocs = (int *)malloc(l*sizeof(int));
    if (Nprocs == NULL) {
        fprintf(stderr, "DistributeVecOpt2(): Not enough memory!\n");
        return FALSE;
    }

    if (!InitNprocs(pM, dir, Nprocs)) {
        fprintf(stderr, "DistributeVecOpt2(): Unable to initialise processor array!\n");
        return FALSE;
    }

    total = 0;
    for (j=0; j<l; j++) {
        if (Nprocs[j] > 2) {
             free(Nprocs);
             return (int)(FALSE);
        }
        total += Nprocs[j];
    }
    
    procstart = (long *)malloc((l+1)*sizeof(long));
    procindex = (int *)malloc(total*sizeof(int));
    if (procstart == NULL || procindex == NULL) {
        fprintf(stderr, "DistributeVecOpt2(): Not enough memory!\n");
        return FALSE;
    }

    if (!InitProcindex(pM, dir, Nprocs, procstart, procindex)) {
        fprintf(stderr, "DistributeVecOpt2(): Unable to initialise processor index!\n");
        return FALSE;
    }
    	
    for (j=0; j<l; j++)
        X[j] = -1;
        
    /********* Phase 1 *********/  
 
    /* If processors q and r share k columns j with Nprocs[j] = 2,
       then we let q and r both own k div 2 columns */
    
    /* The P x P matrix G stores the current number of columns that
       processors q and r share, for all pairs (q, r) */
    G = matallocl(P,P);
    
    if (!G) {
        fprintf(stderr, "DistributeVecOpt2(): Not enough memory!\n");
        return FALSE;
    }
    
    for (q=0; q<P; q++)
        for (r=0; r<P; r++)
            G[q][r] = 0;

    for (j=0; j<l; j++) {
        if (Nprocs[j] == 2) {
            /* find processors q, r sharing column j */
            q= procindex[procstart[j]];
            r= procindex[procstart[j]+1];
            if (q == r) {
                fprintf(stderr, "DistributeVecOpt2(): internal error q = r!\n");
                return FALSE;
            }

            G[q][r]++;
            G[r][q]++;
        }
    }

    /* Compute G2 = initial G mod 2 */
    G2 = matallocl(P,P);
    
    if (!G2) {
        fprintf(stderr, "DistributeVecOpt2(): Not enough memory!\n");
        return FALSE;
    }
    
    for (q=0; q<P; q++)
        for (r=0; r<P; r++)
            G2[q][r] = G[q][r] % 2;
    
    /* Initialise number of vector components Nv */
    Nv = (long *)malloc(P*sizeof(long));
    if (Nv == NULL) {
        fprintf(stderr, "DistributeVecOpt2(): Not enough memory!\n");
        return FALSE;
    }
    for (q=0; q<P; q++)
        Nv[q] = 0;

    /* Assign the shared columns such that processors q and r 
       both own G[q][r] div 2 columns */
    for (j=0; j<l; j++) {
        if (Nprocs[j] == 2) {
            q= procindex[procstart[j]];
            r= procindex[procstart[j]+1];
            if (q>r) {
                /* swap q and r, to make q the smallest processor number */
                tmp= q;
                q= r;
                r= tmp;
            }
            /* G[q][r] >= 1 */
            if (G[q][r] > 1 || G2[q][r] == 0) {
                /* if G[q,r] = 1, we can perform one more assignment
                   provided the initial G[q,r] was even */
               
                /* Alternate the assignments between the processors
                   with smallest (q) and largest (r) processor number */
                if (G[q][r] % 2 == 0) {
                    if (!AssignColumnToProc(X, NULL,NULL,NULL,NULL, Nv, NULL, j,q)) {
                        fprintf(stderr, "DistributeVecOpt2(): Unable to assign column!\n");
                        return FALSE;
                    }
                }
                else {
                    if (!AssignColumnToProc(X, NULL,NULL,NULL,NULL, Nv, NULL, j,r)) {
                        fprintf(stderr, "DistributeVecOpt2(): Unable to assign column!\n");
                        return FALSE;
                    }
                }
                G[q][r]--;
                G[r][q]--;
            }
        }
    }

    /* Now G = G2 must hold. Check this. */
    for (q=0; q<P; q++) {
        for (r=0; r<P; r++) {
            if (G[q][r] != G2[q][r]) {
                fprintf(stderr, "DistributeVecOpt2(): G != G2!\n");
                return FALSE;
            }
        }
    }

    /********* Phases 2, 3 *********/
    
    /* Assign all columns with Nprocs= 2 that are still unowned
       by using a graph algorithm.
       The vertices of the undirected graph are the processors 0,...,P-1.
       The edges of the graph are the pairs (q,r) where q,r share
       an unowned column.
       The data structure used is the P x P adjacency matrix G of the graph.
       This matrix is symmetric and has zero diagonal. */

    degree= (int *)malloc(P*sizeof(int));
    first_adj= (int *)malloc(P*sizeof(int));
    if (degree == NULL || first_adj == NULL) {
        fprintf(stderr, "DistributeVecOpt2(): Not enough memory!\n");
        return FALSE;
    }

    /* Initialise degrees and store column info in Col */
    

    for (q=0; q<P; q++)
        degree[q] = 0;
   
    Col = matallocl(P,P);
    
    if (!Col) {
        fprintf(stderr, "DistributeVecOpt2(): Not enough memory!\n");
        return FALSE;
    }
    
    for (q=0; q<P; q++)
        for (r=0; r<P; r++)
            Col[q][r] = -1;
   
    for (j=0; j<l; j++) {
        if (Nprocs[j] == 2 && X[j] == -1) {
            q= procindex[procstart[j]];
            r= procindex[procstart[j]+1];
            Col[q][r] = j;
            Col[r][q] = j;
            degree[q]++;
            degree[r]++;
        }
    }

    /* Initialise first adjacent vertex. The array first_adj takes care 
       that each row in G has to be traversed only once. */
    for (q=0; q<P; q++)
        first_adj[q]= FirstNzInRow(G,P,q,0);

    /* P loop iterations for odd-degree starting vertices (phase 2),
       followed by P iterations for the remaining starting vertices,
       which can be shown to have even degree (phase 3) */ 
   
    for (count=0; count<2*P; count++) {
        q= count%P; /* starting vertex */
        if (degree[q] > 0 && (degree[q]%2 == 1 || count >= P)) {
            r= first_adj[q];   /* r < P */

            while (r < P) {
                /* Handle edge (q,r) */

                /* Assign column */
                j= Col[q][r];
                if (j == -1) {
                    fprintf(stderr, "DistributeVecOpt2(): Error, column j = -1!\n");
                    return FALSE;
                }
                if (!AssignColumnToProc(X, NULL,NULL,NULL,NULL, Nv, NULL, j,q)) {
                    fprintf(stderr, "DistributeVecOpt2(): Unable to assign column!\n");
                    return FALSE;
                }
                G[q][r]--;
                G[r][q]--;
          
                /* Update degree and first_adj */
                degree[q]--;
                degree[r]--;
                first_adj[q]= FirstNzInRow(G,P,q,first_adj[q]);
                first_adj[r]= FirstNzInRow(G,P,r,first_adj[r]);

                /* Move on */
                q= r;
                r= first_adj[q];
            }
        }
    }

    /* Now G = 0 must hold. Check this. */
    for (q=0; q<P; q++) {
        for (r=0; r<P; r++) {
            if (G[q][r] != 0) {
                fprintf(stderr, "DistributeVecOpt2(): G != 0!\n");
                return FALSE;
            }
        }
    }

    /********* End of phase 3 *********/  

    /* Assign all columns with Nprocs = 1 */
    for (j=0; j<l; j++)
        if (Nprocs[j] == 1)
            X[j] = procindex[procstart[j]]; 
            
    /* Assign all remaining columns */
    if (!AssignRemainingColumns(l,P,X,Nv)) {
        fprintf(stderr, "DistributeVecOpt2(): Unable to assign remaining columns!\n");
        return FALSE;
    }
    
    /* Free all allocated memory */
    matfreel(Col);
    free(first_adj);
    free(degree);
    free(Nv);
    matfreel(G2);
    matfreel(G);
    free(procindex);
    free(procstart);
    free(Nprocs);
     
    return TRUE;
} /* end DistributeVecOpt2 */


int FirstNzInRow (long **G, int P, int q, int start) {
    
    /* This function finds the index of the first nonzero in row q
       of the P x P matrix G, with index >= start. If all elements 
       with index >= start are zero, the function returns the value P.  */

    int Index, i;
    
    if (!G) {
        fprintf(stderr, "FirstNzInRow(): Null argument!\n");
        return 0;
    }

    Index= P;
    for (i=start; i<P; i++)
        if (G[q][i] != 0) {
            Index= i;
            break;
        }

    return Index;

} /* end FirstNzInRow  */

