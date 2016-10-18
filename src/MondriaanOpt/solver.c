#include "solver.h"
#include "matching.h"

inline char isAssigned(unsigned char val){
    /* This function checks whether a value is an assigned value */ 
    return (val==TO0 || val==TO1 || val==CUT) ;
}

unsigned int getNextToAssign(struct mat *a, struct solution *sol){

    /* This function returns the key (0 <= key < m+n) of a row or column
       that is the next to be assigned */

    int i, c;

    /* First check the priority queue for implicitly assigned rows or columns,
       as these can immediately be assigned to CUT */
    while(sol->nprior>0){
        /* Priority queue is not empty. Remove a key. */
        sol->nprior--;
        unsigned int key = sol->prior[sol->nprior];

        /* Determine corresponding row or column index */
        i = key;
        unsigned char dir = ROW;
        if(key>=a->m){
          dir=COL;
          i-=a->m;
        }

        /* Find first implicitly cut row/column, with nonzeros assigned to both parts */
        if(sol->vals[dir][i]==UNAS && sol->n[TO0][dir][i]>0 && sol->n[TO1][dir][i]>0)
            return key;
    }

    /* Then check the main in decreasing order of row/column length */

    /* Search in the other direction than the previous assignment,
       as we prefer alternating assignment of rows and column */
    for(c=a->cmax; c>=0; c--){
      if(sol->nexts[1-sol->prevdir][c].next != NULL)
          return sol->nexts[1-sol->prevdir][c].next->idx;
    }
    /* Search in the direction of the previous assignment */
    for(c=a->cmax; c>=0; c--){
      if(sol->nexts[sol->prevdir][c].next != NULL)
          return sol->nexts[sol->prevdir][c].next->idx;
    }

    /* Otherwise, simply return the first non-empty unassigned row/col */
    for(i=0; i<a->m; i++){
        if(sol->vals[ROW][i]!=UNAS || a->starts[ROW][i]==a->starts[ROW][i+1])
            continue;
        return i;
    }
    for(i=0; i<a->n; i++){
        if(sol->vals[COL][i]!=UNAS || a->starts[COL][i]==a->starts[COL][i+1])
            continue;
        return i+a->m;
    }

    /* Otherwise, we are done with assigning */
    return a->m+a->n;

} /* getNextToAssign */


unsigned int getLB3(struct mat *a, struct solution *sol){

    /* This function computes the lower bound L_3, reducing load imbalance
       to within the allowed bounds, by cutting the largest rows and columns */
  
    int c, part;

    /* Initialise number of nonzeros in rows/columns partially assigned to 0 and 1 */
    int n[2];
    n[0] = sol->npartial[TO0];
    n[1] = sol->npartial[TO1];
  
    /* Initialise lower bound L_3 */
    unsigned int lb=0;

    /* Loop over both parts */
    for(part=0; part<2; part++){ 
        if(sol->tot[part]+n[part] > sol->max1){
            /* Total number of nonzeros already assigned to this part
                  +  number of nonzeros in rows/columns partially assigned
               is larger than the allowed number max1 */
    
            /* Loop over sorted array, cutting rows/cols until smaller than bound */
            c = a->cmax-1; /* length of largest row/column */

            /* Initialise current number of partially assigned rows/columns of length c */
            int curC = sol->partial[part][c];
            while(sol->tot[part]+n[part] > sol->max1){
                /* Find first row/column to cut */
                while(curC==0){
                  c--;
                  curC = sol->partial[part][c];
                }
  
                /* Cut, reducing the number of nonzeros by c */
                lb++;
                curC--;
                n[part] -= c;
            }
        }
    }

    return lb;

} /* end getLB3 */


unsigned char isBranchFeasible(struct mat *a, struct solution *sol, unsigned char val){

    /* This function checks whether a branch just set to a value is feasible */
  
    /* Check the current load balance if a row or column has just been set */
    if(val==TO0){
        if(sol->tot[TO0] > sol->max1)
            return FALSE;
    }else if(val==TO1){
        if(sol->tot[TO1] > sol->max1)
            return FALSE;
    }
  
    /* Check the lower bounds on the communication volume */
    unsigned int lb1 = sol->vol; /* Number of explicitly cut rows/columns.
                                    There are no implictly cut rows, 
                                    as these will be cut immediately */
    unsigned int lb3 = getLB3(a,sol); /* Number to be cut because of load imbalance */
    unsigned int lb4 = sol->mGraph.nmatch; /* Number to be cut because of conflicting nonzeros */
 
    /* Take maximum of bounds L_3 and L_4 */
    unsigned int lb34 ; 
    lb34 = ( lb4>lb3 ? lb4 : lb3 );
  
    /* Check if current lower bound >= current maximum volume */
    unsigned int curLB = lb1 + lb34;
    if(curLB >= sol->maxvol)
        /* no better solution can be found in this branch of the search tree */
        return FALSE;

    return TRUE;

} /* end isBranchFeasible */


void setPartial(struct mat *a, struct solution *sol, unsigned int idx,
                unsigned char dir, unsigned char val, unsigned char set){

    /* For set=TRUE, this function updates the partial counters (for LB3) based on
       the currently assigned row/col idx where 0 <= idx < m for dir=ROW,
       and 0 <= idx < n for dir=COL. The assignment is either to part 0,
       to part1 or cut. For set=FALSE, the updates are undone.
       Comments are for dir=ROW and for set=TRUE. */

    int j, k, c, incr;
    unsigned char part;

    /* Determine increment for set and unset */
    if(set)
        incr = 1;
    else
        incr = -1;

    for(part=0; part<2; part++){
        if( sol->n[part][dir][idx]>0 && sol->n[1-part][dir][idx]==0 ){
            /* row idx was partially assigned to part */
            
            /* Determine the number of unassigned nonzeros in row idx */
            c = a->size[dir][idx] - sol->n[part][dir][idx] ;
            
            /* There is one row less of size c */
            sol->partial[part][c] -= incr;
                
            /* Update the total number of unassigned nonzeros in rows or columns
               partially assigned to part. First, consider all nonzeros of 
               row idx to be assigned now */
            sol->npartial[part] -= incr*c; 
            
            if(val==CUT){
                /* Put nonzeros back in, if their column is partially assigned */
                for (k = a->starts[dir][idx]; k < a->starts[dir][idx+1]; k++){
                    j = a->indices[dir][k];
                    if(sol->vals[1-dir][j]==UNAS && 
                       sol->n[part][1-dir][j]>0 && sol->n[1-part][1-dir][j]==0)
                           sol->npartial[part] += incr;
                }
            }
        }
    }

} /* end setPartial */


void set01(struct mat *a, struct solution *sol, unsigned int idx,
           unsigned char dir, unsigned char vl){

    /* This function sets the current row/column idx to part value vl (0 or 1),
       where 0 <= idx < m for dir=ROW, and 0 <= idx < n for dir=COL.
       Comments are for dir=ROW. */


    int j, k, i, l, c;

    setPartial(a,sol,idx,dir,vl,TRUE);

    /* Loop over all nonzeros in row idx */
    for( j = a->starts[dir][idx]; j < a->starts[dir][idx+1]; j++){

        k =  a->indices[dir][j] ; /* nonzero a(idx,k) */
        sol->n[vl][1-dir][k]++; /* column k has one more nonzero assigned to vl */

        if(sol->vals[1-dir][k]!=UNAS)
            continue; 

        /* column k is unassigned */

        if(sol->n[1-vl][1-dir][k]==0){
            /* Column k is partially assigned to vl */

            /* Determine number of unassigned nonzeros in column k */
            c = a->size[1-dir][k] - sol->n[vl][1-dir][k];

            /* one more partially assigned column with c unassigned nonzeros */
            sol->partial[vl][c]++; 

            if(sol->n[vl][1-dir][k] > 1){
                /* column k was already partially assigned to vl */
                sol->partial[vl][c+1]--; /* it now has one unassigned nonzero less */ 

                if(sol->n[vl][dir][idx]==0)
                    sol->npartial[vl]--;
            } else {
                /* column k just became partially assigned, so update the counters */
                sol->partialstate[1-dir][k] = PAR0+vl; /* PAR0 or PAR1 */
                sol->justpartial[sol->njustpartial] = k;
                sol->njustpartial++;

                /* Loop over all nonzeros in column k */
                for( i = a->starts[1-dir][k]; i < a->starts[1-dir][k+1]; i++ ){
                    l =  a->indices[1-dir][i];
                    if(l==idx)
                        continue; /* row l is the assigned row idx */

                    if(sol->vals[dir][l]==UNAS &&
                       sol->n[vl][dir][l]>0 && sol->n[1-vl][dir][l]==0)
                        /* row l is partially asssigned to part vl */
                        continue; /* this nonzero has already been counted */
                    sol->npartial[vl]++;
                }
            }
        }

        if(sol->partialstate[1-dir][k]==PAR0+1-vl){
            /* column k was partially assigned to 1-vl, so now implicitly cut */
            sol->partialstate[1-dir][k] = CUT;
            sol->justpartialoff[sol->njustpartialoff] = k;
            sol->njustpartialoff++;

            /* Add to priority queue */
            sol->prior[sol->nprior] = k+(1-dir)*a->m;
            sol->nprior++;

            /* Update partial counters */
            c = a->size[1-dir][k] - sol->n[1-vl][1-dir][k];
            sol->partial[1-vl][c]--;
            sol->npartial[1-vl] -= c ;

            /* Loop over all nonzeros in column k */
            for( i = a->starts[1-dir][k]; i < a->starts[1-dir][k+1]; i++ ){ 
                l =  a->indices[1-dir][i] ;
                if(l==idx)
                    continue;

                if(sol->vals[dir][l]==UNAS && sol->n[1-vl][dir][l]>0 && sol->n[vl][dir][l]==0)
                    /* count removed nonzero again, because its row l is partially assigned */ 
                    sol->npartial[1-vl]++;
            }
        }
    }
  
    sol->tot[vl] += a->size[dir][idx] - sol->n[vl][dir][idx];
    sol->vals[dir][idx] = vl;

} /* end set01 */


void unset01(struct mat *a, struct solution *sol, unsigned int idx,
             unsigned char dir, unsigned char oldval, unsigned char vl){

    /* This function resets the current row/column idx from part value vl (0 or 1)
       to the old value oldval, where 0 <= idx < m for dir=ROW,
       and 0 <= idx < n for dir=COL.  It undoes the working of function set01.
       Comments are for dir=ROW. */

    int j, k, i, l, c;
  
    /* Loop over all nonzeros in row idx */
    for (j = a->starts[dir][idx]; j < a->starts[dir][idx+1]; j++ ){ 

        k = a->indices[dir][j]; /* nonzero a(idx,k) */
        sol->n[vl][1-dir][k]--; /* column k has one nonzero less assigned to vl */

        if(sol->vals[1-dir][k]!=UNAS)
            continue; /* column was assigned so no action was taken by set01 */
 
        if(sol->n[1-vl][1-dir][k]==0){
            /* Column k was partially assigned to vl */

            /* Determine number of unassigned nonzeros in column k */
            c = a->size[1-dir][k] - sol->n[vl][1-dir][k] ;

            /* one less partially assigned row with c-1 unassigned nonzeros */
            sol->partial[vl][c-1]--;

            if(sol->n[vl][1-dir][k] > 0){
                /* column k was already partially assigned to vl before calling set01 */
                sol->partial[vl][c]++; /* it now has one unassigned nonzero more */

                if(sol->n[vl][dir][idx]==0)
                    sol->npartial[vl]++;
            } else {
                /* column k is not partially assigned anymore */
                sol->justpartialoff[sol->njustpartialoff] = k;
                sol->njustpartialoff++;
                sol->partialstate[1-dir][k] = UNAS;

                /* Loop over all nonzeros in column k */
                for( i = a->starts[1-dir][k]; i < a->starts[1-dir][k+1]; i++ ){ 
                    l = a->indices[1-dir][i];
                    if(l==idx)
                        continue; /* row l is the row that was assigned */
  
                    if(sol->vals[dir][l]==UNAS &&
                       sol->n[vl][dir][l]>0 && sol->n[1-vl][dir][l]==0)
                        /* row l is partially assigned to part vl */
                        continue; /* nothing was done by set01 */
                    sol->npartial[vl]--;
                }
            }
        }

        if(sol->n[vl][1-dir][k]==0 && sol->n[1-vl][1-dir][k]>0){
            /* column k was partially assigned to 1-vl before calling set01,
               became implicitly cut and now reverts */
            sol->justpartial[sol->njustpartial] = k;
            sol->njustpartial++;
            sol->partialstate[1-dir][k] = PAR0+1-vl;

            /* No need to remove column k from the piority queue,
               as this has already been done */

            /* Bring partial counters back to the previous state */
            c  = a->size[1-dir][k] - sol->n[1-vl][1-dir][k] ;
            sol->partial[1-vl][c]++;
            sol->npartial[1-vl] += c;

            /* Loop over all nonzeros in column k */
            for(i = a->starts[1-dir][k]; i < a->starts[1-dir][k+1]; i++ ){ 
                l =  a->indices[1-dir][i];
                if(l==idx)
                    continue;

                if(sol->vals[dir][l]==UNAS &&
                   sol->n[1-vl][dir][l]>0 && sol->n[vl][dir][l]==0)
                    /* count nonzero again, because its row l is partially assigned */
                    sol->npartial[1-vl]--;
            }
        }
    }

    sol->tot[vl] -= a->size[dir][idx] - sol->n[vl][dir][idx];
    sol->vals[dir][idx] = oldval;
    setPartial(a,sol,idx,dir,vl,FALSE);

} /* end unset01 */


void setCut(struct mat *a, struct solution *sol, unsigned int idx, unsigned char dir,
            unsigned char oldval, unsigned char set){

    /* For set=TRUE, this function sets row/column idx to CUT,
       increments the communication volume, and updates the partial counters
       where 0 <= idx < m for dir=ROW, and 0 <= idx < n for dir=COL.
       For set=FALSE, this function undoes the setting and
       restores the solution value of row/column to oldval */

    /* The following three operations can be done in any order */ 
    setPartial(a,sol,idx,dir,CUT,set);
    if(set){
        sol->vals[dir][idx] = CUT;
        sol->vol++;
    } else {
       sol->vals[dir][idx] = oldval;
       sol->vol--;
    }

} /* end setCut */


void setNextToAssign(struct mat *a, struct solution *sol, unsigned int idx,
                     unsigned char dir, unsigned char set){
    /* This function updates the priority list of the next row or column to pick. 
       If set=TRUE, it removes the row/column that was just assigned,
       otherwise it adds that row/column back. */

    int i, c;

    i = idx;
    if(dir==COL)
        i += a->m;
    c = sol->idxi[dir][idx].key ; /* length of row/column idx */

    if(set)
        delete_priority(sol, dir, c, i);
    else 
        insert_priority(sol, dir, c, i);

} /* setNextToAssign */


void assign(struct mat *a, struct solution *sol, unsigned int idx, unsigned char val){
  
    /* This recursive function assigns row/column idx to part val (0, 1, or CUT)
       and then calls itself for the remaining rows and columns until a solution is found */

    /* Convert row/column index idx, 0 <= idx < m+n, to idx for row or column separately,
       and compute the corresponding direction */
    unsigned char dir = ROW;
    if(idx >= a->m){
        idx -= a->m;
        dir = COL;
    }
  
    unsigned char oldval = sol->vals[dir][idx];
  
    /* Partially assigned rows/columns have only 1 choice */
    if(val==TO0 && sol->n[TO1][dir][idx]>0)
        return;
    
    if(val==TO1 && sol->n[TO0][dir][idx]>0)
        return;

    /* A row to be cut cannot have 1 nonzero, or have all nonzeros
       partially assigned to the same part */
    if(val==CUT && (a->size[dir][idx]==1 || sol->n[TO0][dir][idx]==a->size[dir][idx] ||
                                            sol->n[TO1][dir][idx]==a->size[dir][idx]))
        return;
  
    /* Register current depth and total number of tree nodes (branches) encountered */
    sol->depth++;
    sol->nbranches++;
    if(sol->nbranches%PERIOD == 0){
        /* Print current depth and check whether maximum run time has been exceeded */
        fprintf(stderr,"%d/%d \n", sol->depth, a->m+a->n);
        if(sol->maxruntime > 1){
            if(difftime(time(NULL),sol->starttime) > sol->maxruntime){
                printf("(STOPPED) %d\n",sol->maxvol);
                exit(EXIT_SUCCESS);
            }
        }
    }
  
    /* Set the curent row/column to the assigned value */
    sol->njustpartial = 0;
    sol->njustpartialoff = 0;
    if(val==TO0){
        set01(a,sol,idx,dir,TO0);
    }else if(val==TO1){
        set01(a,sol,idx,dir,TO1);
    }else{
        setCut(a,sol,idx,dir,CUT,TRUE);
    }
  
    /* Update the graph for the matching bound L4 */
    update_graph_after_setting(a,sol,idx,dir,val);
  
    if(isBranchFeasible(a,sol,val)){
        /* Remove idx from the priority list */
        setNextToAssign(a,sol,idx,dir,TRUE);

        /* Keep previous direction in case of assignment to a part */
        unsigned char olddir = sol->prevdir;
        if(val==TO0 || val==TO1)
            sol->prevdir = dir;

        /* Get next row or column to assign */
        unsigned int next = getNextToAssign(a,sol);

        if(next == a->m+a->n){
            /* all rows and columns have been assigned and a better solution was found */
            unsigned int newvol = sol->vol;
            fprintf(stderr,"NEW SOLUTION FOUND! (vol: %d, load: %d, %d, %d)\n", newvol, sol->tot[0],sol->tot[1],a->nnz - sol->tot[0] - sol->tot[1]);
            sol->maxvol = newvol;
            char fn[MAXFNSIZE];

            /* Output solution as SVG picture */
            sprintf(fn,"%s-2f.svg",sol->matname);
            printToSVG(fn,sol,a);

            /* Output solution as sparse integer matrix in Matrix Market format */
            sprintf(fn,"%s-I2f",sol->matname);
            printToMM(fn,sol,a);
        }else{
            /* Preserve old state of next */
            unsigned char oldstate;
            if(next >= a->m){
                oldstate = sol->partialstate[COL][next-a->m];
                sol->partialstate[COL][next-a->m] = ASGN;
            }else{
                oldstate = sol->partialstate[ROW][next];
                sol->partialstate[ROW][next] = ASGN;
            }

            /* Exploit symmetry by assigning first non-cut row/column only to part 0 */
            if(val==TO0 && sol->FirstIdxAssignedTo0==FALSE){
                sol->FirstIdxAssignedTo0 = TRUE;
                assign(a,sol,next,TO1);
                assign(a,sol,next,TO0);
                assign(a,sol,next,CUT);
                sol->FirstIdxAssignedTo0 = FALSE;
            }else{
                /* Assign first to the part with the least nonzeros */
                if(sol->tot[TO0] <= sol->tot[TO1]){
                    assign(a,sol,next,TO0);
                    if(sol->FirstIdxAssignedTo0==TRUE)
                        assign(a,sol,next,TO1);
                }else{
                    if(sol->FirstIdxAssignedTo0==TRUE)
                        assign(a,sol,next,TO1);
                    assign(a,sol,next,TO0);
                }
                assign(a,sol,next,CUT);
            }

            /* Restore old state of next */
            if(next >= a->m)
                sol->partialstate[COL][next-a->m] = oldstate;
            else
                sol->partialstate[ROW][next] = oldstate;
            
        }
        sol->prevdir = olddir;

        /* Add idx back to the priority list */
        setNextToAssign(a,sol,idx,dir,FALSE);
    }
  
    /* Unset the current row/column from the assigned value to the old value */
    sol->njustpartial = 0;
    sol->njustpartialoff = 0;
    if(val==TO0){
        unset01(a,sol,idx,dir,oldval,TO0);
    }else if(val==TO1){
        unset01(a,sol,idx,dir,oldval,TO1);
    }else{
        setCut(a,sol,idx,dir,oldval,FALSE); /* FALSE means unsetCut */
    }

    /* Update the graph for the matching bound L4 */
    update_graph_after_unsetting(a,sol,idx,dir,val);

    sol->depth--;

} /* end assign */


void solve(struct mat *a, struct solution *sol){
    
    /* This function solves the bipartitioning problem for the sparse matrix a
     by a branch-and-bound algorithm and stores the solution in sol */
    
    /* Initialise the solution, setting all values to 0 */
    resetsolution(a,sol);
    
    /* Get the index of the first row/column to assign, 0 <= first < m+n */
    unsigned int first = getNextToAssign(a,sol);
    
    /* Try assigning the row/column to part 0 and CUT.
       No need to try part 1, because of symmetry */
    assign(a,sol,first,TO0);
    assign(a,sol,first,CUT);
    
} /* end solve */

