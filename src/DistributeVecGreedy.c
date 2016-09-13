#include "DistributeVecGreedy.h"

/* All comments only consider the distribution of the vector v,
   where all vector-length arrays are of length pM->n.
   This corresponds to the direction dir = ROW,
   meaning that the vector is in the row direction. 
   In case of the vector u =  Av, i.e., dir = COL,
   one should read pM->m instead. */

long DistributeVecGreedyImprove(const struct sparsematrix *pM, long int *X, int dir, const struct opts *pOptions) {
    /*  This function improves a given vector distribution X by randomly choosing
        a vector component and reassigning it to the best processor, reducing
        the communication cost or keeping it equal. This is repeated until no
        further improvement is possible or a given maximum number of tries
        is reached.
        Input:  initialised A,
                array X of length pM->n with values between 0 and P-1, 
                where P= pM->NrProcs.
        Output: improved vector distribution X.
                The function returns maxcom = max(sends,recvs)
                over all processors. */

 
    long j, k, l=0, s, total, ncols, ncols0,
         before_q, after_q, before_r, after_r,
         min, count, maxcom;
    int P, q, r, rmin, q_occurs_in_col;

    /* Arrays of length pM->n */
    int *Nprocs;    /* Nprocs[j] is the number of processors that occur
                       in column j */
    long *cols;     /* contains the column indices j=0,..,pM->n-1 in permuted
                       order. The first ncols columns are candidates for
                       improvement */

    /* Communication matrix */
    int *procindex; /* array of length at most V+pM->n, where V is the
                       communication volume for the vector v.
                       procindex[procstart[j]..procstart[j+1]-1]
                       contains the processor numbers of the processors
                       that occur in column j */
    long *procstart; /* array of length pM->n+1 containing the starts of the
                       columns of the communication matrix */

    /* Arrays of length P */
    long *Ns;       /* Ns[q] = the current number of vector components
                       sent by processor q to other processors */
    long *Nr;       /* Nr[q] = the current number of vector components
                       received by processor q from other processors */
    long *Nv;       /* Nv[q] = the current number of vector components
                       owned by processor q */
    if (!pM || !X || !pOptions) {
        fprintf(stderr, "DistributeVecGreedyImprove(): Null parameters!\n");
        return -1;
    }
    
    if (dir == COL)
        l = pM->m;
    else if (dir == ROW)
        l = pM->n;
    else {
        fprintf(stderr, "DistributeVecGreedyImprove(): Unknown direction!\n");
        return -1;
    }
    P = pM->NrProcs; 

    cols = (long *)malloc(l*sizeof(long));
    Nprocs = (int *)malloc(l*sizeof(int));
    procstart = (long *)malloc((l+1)*sizeof(long));

    Ns = (long *)malloc(P*sizeof(long));
    Nr = (long *)malloc(P*sizeof(long));
    Nv = (long *)malloc(P*sizeof(long));

    if (cols == NULL || Nprocs == NULL || procstart == NULL ||
        Ns == NULL || Nr == NULL || Nv == NULL) {
        fprintf(stderr, "DistributeVecGreedyImprove(): Not enough memory!\n");
        return -1;
    }
 
    if (!InitNprocs(pM, dir, Nprocs)) {
        fprintf(stderr, "DistributeVecGreedyImprove(): Unable to initialise processor array!\n");
        return -1;
    }

    total = 0;
    for (j=0; j<l; j++)
        total += Nprocs[j];
    procindex = (int *)malloc(total*sizeof(int));
    if (procindex == NULL) {
        fprintf(stderr, "DistributeVecGreedyImprove(): Not enough memory!\n");
        return -1;
    }

    if (!InitProcindex(pM, dir, Nprocs, procstart, procindex)) {
        fprintf(stderr, "DistributeVecGreedyImprove(): Unable to initialise processor index!\n");
        return -1;
    }

    /* Initialise */
    for (q=0; q<P; q++) {
        Ns[q] = 0;
        Nr[q] = 0;
        Nv[q] = 0;
    }
    
    ncols= 0; 

    /* Set Ns, Nr, Nv for initial assignement by X */
    for (j=0; j<l; j++) {
        if (X[j] < 0 || X[j] >= P) {
            fprintf(stderr, "DistributeVecGreedyImprove(): X[j] not in 0..P-1!\n");
            return -1;
        }
        q = X[j];
        X[j] = -1; /* temporarily */

        /* Check whether q occurs in column j */
        q_occurs_in_col = FALSE;
        for (s=procstart[j]; s<procstart[j+1]; s++) {
            if (procindex[s] == q) {
                q_occurs_in_col = TRUE;
                break;
            }
        }

        /* Assign column j to a processor in the column, if possible.
           Motivation: we do not want to increase the communication volume. */

        if (Nprocs[j] == 0  || q_occurs_in_col) {
            if (!AssignColumnToProc(X, procstart, procindex, Ns, Nr, Nv, NULL, j, q)) {
                fprintf(stderr, "DistributeVecGreedyImprove(): Unable to assign column!\n");
                return -1;
            }
        }
        else {
            /* Choose a random processor r in the column */
            s = Random1(0,Nprocs[j]-1);
            r = procindex[procstart[j]+s];
            if (!AssignColumnToProc(X, procstart, procindex, Ns, Nr, Nv, NULL, j, r)) {
                fprintf(stderr, "DistributeVecGreedyImprove(): Unable to assign column!\n");
                return -1;
            }
        }

        /* Only columns with > 1 processors can still be improved
           without increasing the communication volume */
        if (Nprocs[j] > 1) { 
            cols[ncols] = j;
            ncols++;
        }
    }

    ncols0 = ncols; /* number of columns with Nprocs >1 */
    count=0;

#ifdef INFO2
    printf("Improve greedily\n");
    printf(" Ns, Nr = number of sends/recvs\n");
    printf(" Nv = number of vector components\n");
    printf(" After %ld columns tried:\n", count);
    PrintVecStatistics(P, Ns, Nr, Nv);
#endif

    while (ncols>0 && count< pOptions->VectorPartition_MaxNrGreedyImproves*l) {
        /* Choose a column */
        k = Random1(0,ncols-1); 
        j = cols[k];
        q = X[j];

        before_q = MAX(Ns[q], Nr[q]); /* cost of proc q before moving column j 
                                         from q to another proc */

        /* q occurs in column j and Nprocs[j] > 1 */
        after_q = MAX(Ns[q]- (Nprocs[j]-1), Nr[q]+1);

        /* Choose best processor r to own column j (based on Ns, Nr, Nv).
           Initialisation: best move is no move */
        rmin = q;  
        min = LONG_MAX;
 
        for (s = procstart[j]; s<procstart[j+1]; s++) {
            r = procindex[s]; 
            if (r==q)
                continue;
            before_r = MAX(Ns[r], Nr[r]); /* cost of proc r before moving
                                             column j from q to r */
            after_r = MAX(Ns[r]+Nprocs[j]-1, Nr[r]-1);

            /* Primary criterion: to achieve a reduction in cost.
               The size of the reduction does not matter much for
               choosing the receiver, since the reduction is either
               determined by the sender (and it is >= 1),
               or by the receivers (=1).
               Secondary criterion: minimal number of sends. */
            if (MAX(after_q, after_r) < MAX(before_q, before_r) && Ns[r] < min) {
                /* A move to r is the current best move */ 
                rmin = r;
                min = Ns[r];
            }
        } 

        if (rmin!=q) {
            /* Improvement in communication balance achieved */
            if (!RemoveColumnFromProc(X, procstart, procindex, Ns, Nr, Nv, NULL, j, q)) {
                fprintf(stderr, "DistributeVecGreedyImprove(): Unable to remove column!\n");
                return -1;
            }
            if (!AssignColumnToProc(X, procstart, procindex, Ns, Nr, Nv, NULL, j, rmin)) {
                fprintf(stderr, "DistributeVecGreedyImprove(): Unable to assign column!\n");
                return -1;
            }
            /* Reset. All columns with Nprocs>1 are candidates again */
            ncols = ncols0;
        } else {
            /* No move. Try another column */
            SwapLong(cols, k, ncols-1);
            ncols--;
        }

        count++; /* counts number of columns tried */
        
#ifdef INFO2
        if (l>0 && (count%l==0)) {
            printf(" After %ld columns tried:\n", count);
            PrintVecStatistics(P, Ns, Nr, Nv);
        }
#endif
    }
    

    /* Assign all columns with Nprocs = 0 */
    for (j=0; j<l; j++) {
        if (Nprocs[j] == 0) {
            q = X[j];
            if (!RemoveColumnFromProc(X, procstart, procindex, Ns, Nr, Nv, NULL, j, q)) {
                fprintf(stderr, "DistributeVecGreedyImprove(): Unable to remove column!\n");
                return -1;
            }
        }
    }
    
    if (!AssignRemainingColumns(l,P,X,Nv)) {
        fprintf(stderr, "DistributeVecGreedyImprove(): Unable to assign column!\n");
        return -1;
    }

    /* Check result */
    for (j=0; j<l; j++) {
        if (X[j] < 0 || X[j] >= P) {
            fprintf(stderr, "DistributeVecGreedyImprove(): column still unowned!\n");
            return -1;
        }
    }
#ifdef INFO2
    printf(" After %ld columns tried:\n", count);
    PrintVecStatistics(P, Ns, Nr, Nv);
#endif

    maxcom = 0;
    for (q=0; q<P; q++) {
        if (Ns[q] > maxcom)
            maxcom = Ns[q];
        if (Nr[q] > maxcom)
            maxcom = Nr[q];
    }


    free(procindex);
    free(Nv);
    free(Nr);
    free(Ns);
    free(procstart);
    free(Nprocs);
    free(cols);

    return maxcom;
} /* end DistributeVecGreedyImprove */

