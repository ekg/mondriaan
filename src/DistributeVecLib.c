#include "DistributeVecLib.h"

/* All comments only consider the distribution of the vector v,
   where all vector-length arrays are of length A.n.
   This corresponds to the direction dir = ROW,
   meaning that the vector is in the row direction. 
   In case of the vector u =  Av, one should read A.m instead. */

int AssignColumnToProc(long int *X, long *procstart, int *procindex,
                        long *Ns, long *Nr, long *Nv, long *Sums,
                        long j, int q) {

    /* This function assigns column j to processor q, updating the
       owner array X and the number of sends Ns, receives Nr,
       vector elements Nv of all processors. 
       If Sums <> NULL, it also updates the sum of the number of sends
       and receives, assuming that inevitable communications have already
       been included (1 communication for every processor in a matrix column). 
       Column j must not yet been owned.
       Input: 0 <= j < n, 0 <= q < P.
       This function can also be used to update only X and Nv.
       In that case, the parameters procstart, procindex, Ns, Nr, Sums
       must be NULL upon input. */

    long r;
    int q_occurs_in_col, npj;
    
    if (!X || !Nv) {
        fprintf(stderr, "AssignColumnToProc(): Null arguments!\n");
        return FALSE;
    }

    if (X[j] != -1) {
        fprintf(stderr, "AssignColumnToProc(): Column already owned!\n");
        return FALSE;
    }

    /* Assign column to processor */
    X[j] = q;
    Nv[q]++;

    /* Return if sends and receives need not be updated */
    if (procstart == NULL && procindex == NULL &&
        Ns == NULL && Nr == NULL && Sums == NULL)
        return TRUE;

    /* Check whether q occurs in column j and
       add receives to all receiving processors */
    q_occurs_in_col = FALSE;
    for (r=procstart[j]; r<procstart[j+1]; r++) {
        if (procindex[r] == q)
            q_occurs_in_col = TRUE;
        else
            Nr[procindex[r]]++;
    }

    /* Add sends to sending processor and update sums */
    npj = procstart[j+1] - procstart[j]; /* number of processors in column j */
    if (q_occurs_in_col) {
        Ns[q] += npj - 1; /* This is >=0 */
        if (Sums != NULL && npj >= 3)
            Sums[q] += npj - 2;
    } else {
        Ns[q] += npj;
        if (Sums != NULL) 
            Sums[q] += npj;  
    }

    return TRUE;
} /* end AssignColumnToProc */


int RemoveColumnFromProc(long int *X, long *procstart, int *procindex,
                          long *Ns, long *Nr, long *Nv, long *Sums,
                          long j, int q) {
    /* This function removes column j from processor q, updating the 
       owner array X and the number of sends Ns, receives Nr,
       vector elements Nv of all processors.
       If Sums <> NULL, it also updates the sum of the number of sends
       and receives, assuming that inevitable communications have already
       been included (1 communication for every processor in a matrix column). 
       Column j must be owned by q. 
       Input: 0 <= j < n, 0 <= q < P. 
       This function can also be used to update only X and Nv.
       In that case, the parameters procstart, procindex, Ns, Nr, Sums
       must be NULL upon input. */  

    long r;
    int q_occurs_in_col, npj;
    
    if (!X || !Nv) {
        fprintf(stderr, "RemoveColumnFromProc(): Null arguments!\n");
        return FALSE;
    }

    if (X[j] != q) {
        fprintf(stderr, "RemoveColumnFromProc(): Column removed from nonowner!\n");
        return FALSE;
    }

    /* Remove column from processor */
    X[j] = -1;
    Nv[q]--;

    /* Return if sends and receives need not be updated */
    if (procstart == NULL && procindex == NULL &&
        Ns == NULL && Nr == NULL && Sums == NULL)
        return TRUE;

    /* Check whether q occurs in column j and 
       subtract receives from all receiving processors */
    q_occurs_in_col = FALSE;
    for (r=procstart[j]; r<procstart[j+1]; r++) {
        if (procindex[r] == q)
            q_occurs_in_col = TRUE;
        else
            Nr[procindex[r]]--;
    }

    /* Subtract sends from sending processor and update sums */
    npj = procstart[j+1] - procstart[j]; /* number of processors in column j */
    if (q_occurs_in_col) {
        Ns[q] -= npj - 1; /* This is >=0 */
        if (Sums != NULL && npj >= 3)
            Sums[q] -= npj - 2;
    } else {
        Ns[q] -= npj;
        if (Sums != NULL)
            Sums[q] -= npj;  
    }
 
    return TRUE;
} /* end RemoveColumnFromProc */


int AssignRemainingColumns(long l, int P, long int *X, long *Nv) { 
    /* This function assigns all the remaining unowned vector components
       (i.e. columns) trying to balance the number of vector components
       among the processors. 
       Input: l is the number of vector components.
              P is the number of processors.
              X is an array of length l. 
                    If X[j] = -1, then component j is still unowned.
                    Otherwise, 0 <= X[j] < P and X[j] will not be changed.
              Nv is an array of length P. It need not be initialised;
              it will be overwritten.  
       Output: 0 <= X[j] < P for all j. 
               Nv[q] is the number of vector components j with X[j] = q,
                        for 0 <= q < P. */
    
    long max, j;
    int q;
    
    if (!X || !Nv) {
        fprintf(stderr, "AssignRemainingColumns(): Null arguments!\n");
        return FALSE;
    }
    
    /* Count the number of vector components Nv[q] owned by processor q */
    for (q=0; q < P; q++)
        Nv[q] = 0;

    for (j=0; j < l; j++) {
        if (X[j] != -1) {
            if (X[j] < 0 || X[j] > P) {
                fprintf(stderr, "AssignRemainingColumns(): Error in owner!\n");
                return FALSE;
            }
            Nv[X[j]]++;
        }
    }        

    /* Determine current maximum number of vector components per processor */
    max= 0;
    for (q=0; q < P; q++)
        if (Nv[q] > max)
            max = Nv[q];
    
    /* Assign columns until all processors have reached max */
    q=0;
    for (j=0; j < l; j++) {
        if (X[j] == -1) {
            /* find first available processor */
            while(q < P && Nv[q] == max)
                q++;
            
            if (q==P) {
                /* All processors have reached max */
                break;
            }
            else {
                /* Nv[q] < max */
                if (!AssignColumnToProc(X, NULL,NULL,NULL,NULL, Nv, NULL,j,q)) {
                    fprintf(stderr, "AssignRemainingColumns(): Unable to assign column!\n");
                    return FALSE;
                }
            }
        }
    }

    /* Assign remaining columns cyclically */
    q=0;
    for (j=0; j < l; j++) {
        if (X[j] == -1) {
            if (!AssignColumnToProc(X, NULL,NULL,NULL,NULL, Nv, NULL, j,q)) {
                fprintf(stderr, "AssignRemainingColumns(): Unable to assign column!\n");
                return FALSE;
            }
            q++;
            if (q == P)
                q = 0;
        }
    }
    
    return TRUE;
} /* end AssignRemainingColumns */


int AssignRemainingNonemptyColumns(long l, int P, long int *X, 
                                    long *procstart, int *procindex,
                                    long *Ns, long *Nr, long *Nv) {

    /* This function greedily assigns all unowned columns with Nprocs >= 1 
       to the processor q with minimal cost, updating the 
       owner array X and the number of sends Ns, receives Nr,
       vector elements Nv of all processors.
       Input: l is the number of vector components.
              P is the number of processors.
              X is an array of length l. 
                If X[j] = -1, then component j is still unowned.
                Otherwise, 0 <= X[j] < P and X[j] will not be changed.
 
       Output: 0 <= X[j] < P for all j, except if Nprocs[j]=0. */
    
    long j, s, min;
    int q, r, npj;
    
    if (!X || !Nv) {
        fprintf(stderr, "AssignRemainingNonemptyColumns(): Null arguments!\n");
        return FALSE;
    }
       
    for (j=0; j<l; j++) {
        if (X[j] != -1)
            continue;
        npj = procstart[j+1] - procstart[j];
        if (npj == 1) {
            q = procindex[procstart[j]];
            if (!AssignColumnToProc(X, procstart, procindex, Ns, Nr, Nv, NULL, j, q)) {
                fprintf(stderr, "AssignRemainingNonemptyColumns(): Unable to assign column!\n");
                return FALSE;
            }
        } else if (npj > 1) {
            /* Add receives temporarily to all processors in the column,
               so we need to consider only one processor at a time */
            for (s=procstart[j]; s<procstart[j+1]; s++) {
                r = procindex[s];
                Nr[r]++;
            }

            /* Find processor q with minimum cost */
            q = procindex[procstart[j]];
            min = MAX(Ns[q]+npj-1, Nr[q]-1);

            for (s=procstart[j]+1; s<procstart[j+1]; s++) {
                r = procindex[s];
                if (MAX(Ns[r]+npj-1, Nr[r]-1) < min) {
                    q = r;
                    min = MAX(Ns[r]+npj-1,Nr[r]-1);
                }
            }

            /* Remove temporary receives */
            for (s=procstart[j]; s<procstart[j+1]; s++) {
                r = procindex[s];
                Nr[r]--;
            }  

            if (!AssignColumnToProc(X, procstart, procindex, Ns, Nr, Nv, NULL, j, q)) {
                fprintf(stderr, "AssignRemainingNonemptyColumns(): Unable to assign column!\n");
                return FALSE;
            }
        }
    }
    
    return TRUE;
} /* end AssignRemainingNonemptyColumns */


int InitNprocs(const struct sparsematrix *pM, int dir, int *Nprocs) {
    /* This function initialises the array Nprocs such that 
       Nprocs[j] is the number of processors that own
       nonzeros in matrix column j (in case dir = ROW),
       or matrix row j (in case dir = COL).

       The nonzeros of the matrix have been sorted by increasing
       processor number, as given by the array pM->Pstart[0..P],
       where P=pM->NrProcs and pM->Pstart[P] = pM->NrNzElts. */


    int *occur;     /* occur[j]= q means processor q occurs in column j
                       and is its current highest-numbered processor.
                       occur[j]= -1 means no processor occurs yet in column j */
    int q;          /* q is the current processor number */

    long j, t, l=0, *Index=NULL;
    
    if (!pM || !Nprocs) {
        fprintf(stderr, "InitNprocs(): Null arguments!\n");
        return FALSE;
    }

    if (dir == COL) {
        l = pM->m;
        Index = pM->i;
    } else if (dir == ROW) {
        l = pM->n;
        Index = pM->j;
    } else {
        fprintf(stderr, "InitNprocs(): Unknown direction!\n");
        return FALSE;
    }

    occur = (int *)malloc(l*sizeof(int));
    if (occur == NULL) {
        fprintf(stderr, "InitNprocs(): Not enough memory!\n");
        return FALSE;
    }

    for (j=0; j<l; j++) {
        Nprocs[j] = 0;
        occur[j] = -1;
    }

    /* Initialise Nprocs to the number of processors in the column */
    q=0;
    for (t=0; t<pM->NrNzElts; t++) {
        while(t == pM->Pstart[q+1])
            q++;  /* this terminates because t != pM->Pstart[P] = pM->NrNzElts */

        j = Index[t];
        if (occur[j] != q) {
            occur[j] = q;
            Nprocs[j]++;
        }
    }

    free(occur);
    
    return TRUE;
} /* end InitNprocs */


int InitProcindex(const struct sparsematrix *pM, int dir, int *Nprocs,
                   long *procstart, int *procindex) {

    /* This function initialises procstart and procindex using Nprocs.
       Input:  Nprocs, where Nprocs[j] is the number of processors that own
                       nonzeros in matrix column j (in case dir = ROW).
       Output: procindex, contains the processor numbers occurring in 
                       the matrix columns. The processor numbers for one
                       column j are stored consecutively.
               procstart, where procstart[j] is the position in array procindex
                       of the first processor number of column j, 0 <= j < pM->n.
                       procstart[pM->n] = total number of processor numbers
                       of all the columns together.
       Together, InitNprocs and InitProcindex compute the P x n sparse
       communication matrix C, stored in compressed column storage (CCS)
       format, and corresponding to the m x n distributed sparse matrix pM->  */

    int *occur;     /* occur[j]= q means processor q occurs in column j
                       and is its current highest-numbered processor.
                       occur[j]= -1 means no processor occurs yet in column j */
    int q;          /* q is the current processor number */

    long j, t, l=0, total, *Index=NULL;
    
    if (!pM || !Nprocs || !procstart || !procindex) {
        fprintf(stderr, "InitProcindex(): Null arguments!\n");
        return FALSE;
    }

    if (dir == COL) {
        l = pM->m;
        Index = pM->i;
    } else if (dir == ROW) {
        l = pM->n;
        Index = pM->j;
    } else {
        fprintf(stderr, "InitProcindex(): Unknown direction!\n");
        return FALSE;
    }

    /* Initialise procstart */
    total= 0;
    for (j=0; j<l; j++) {
        procstart[j] = total;
        total += Nprocs[j];
    }
    procstart[l] = total;

    occur = (int *)malloc(l*sizeof(int));
    if (occur == NULL) {
        fprintf(stderr, "InitProcindex(): Not enough memory!\n");
    }
    
    for (j=0; j<l; j++)
        occur[j] = -1;

    /* Initialise procindex, modifying procstart */
    q=0;
    for (t=0; t<pM->NrNzElts; t++) {
        while(t == pM->Pstart[q+1])
            q++;  /* this terminates because t != pM->Pstart[P] = pM->NrNzElts */

        j = Index[t];
        if (occur[j] != q) {
            occur[j] = q;
            procindex[procstart[j]] = q;
            procstart[j]++;
        }
    }

    /* Reinitialise procstart */
    total= 0;
    for (j=0; j<l; j++) {
        procstart[j] = total;
        total += Nprocs[j];
    }
    procstart[l] = total;

    free(occur);
    
    return TRUE;
} /* end InitProcindex */


int GenerateHistogram(int *X, long l, int lo, int hi, long *Histogram) {

    /* This function computes a histogram for the array X of length l,
       for function values X[j] in the interval lo..hi.
       Output: Histogram[q-lo] = |{j: 0 <= j < l and X[j] = q }|
                                 for q=lo,..,hi.
       The array Histogram must be at least of length hi-lo+1. */

    long j;
    int q;
    
    if (!X || !Histogram) {
        fprintf(stderr, "GenerateHistogram(): Null arguments!\n");
        return FALSE;
    }
  
    /* Initialise histogram */
    for (q = lo; q <= hi; q++)
        Histogram[q-lo]= 0;
    
    /* Fill histogram */
    for (j = 0; j < l; j++)
        if (lo <= X[j] && X[j] <= hi)
            Histogram[X[j]-lo]++;
    
    return TRUE;
} /* end GenerateHistogram */


int PrintHistogram(int lo, int hi, long *Histogram) {

    /* This function prints a histogram for values lo..hi
       but only up to the maximum value that occurs.
       The function  prints q and Histogram[q-lo] for lo <= q <= qmax.
       Here, lo <= qmax <= hi.
       The array Histogram must be at least of length hi-lo+1. */

    int q, qmax;
    
    if (!Histogram) {
        fprintf(stderr, "PrintHistogram(): Null argument!\n");
        return FALSE;
    }
  
    /* Determine maximum occurring value qmax */
    qmax = lo; /* qmax >= lo, to print at least one value */
    for (q = hi; q >= lo; q--)
        if (Histogram[q-lo]>0) {
            qmax = q;
            break;
        }
             
    for (q = lo; q <= qmax; q++)
        printf("  p=%d: %ld \n", q, Histogram[q-lo]);
 
    return TRUE;
} /* end PrintHistogram */


int CalcCom(const struct sparsematrix *pM, long int *X, int dir, long *ComVol,
             long *MaxOut, long *MaxIn, long *MaxCompnts, long *TotCompnts) {

    /* This function calculates the total and maximum communication volume
       for a given distributed sparse matrix A and a distributed vector X.
       Note that these numbers are for only one direction (vector).
       If X=NULL, the volume is only based on the matrix 
                  and no vector-based statistics are computed,
                  i.e. no MaxOut, MaxIn, MaxCompnts, TotCompnts.
       The total volume may be higher than the matrix-based volume.

       This function can also be used if the vector is only partially
       distributed. Unowned components (with X[j] = -1) are ignored.

       Input: A distributed matrix,
              X distribution of vector, with 0 <= X[j] < P, for 0 <= j < pM->n.
                  X[j] is the owner of vector component j.
              dir = ROW (for distribution of v) or dir = COL (for u).
       Output: ComVol is total communication volume for direction dir,
               MaxOut is maximum number of data words sent by a processor,
               MaxIn is maximum number of data words received by a processor,
               MaxCompnts is maximum number of vector components of a processor.
               TotCompnts is total number of owned vector components.
    */
    
    int P, q, *Nprocs, *procindex;
    long total, maxs, maxr, maxv, tots, totr, totv,
         l=0, j, *procstart, *Ns, *Nr, *Nv;
    
    if (!pM || !ComVol || !MaxOut || !MaxIn || !MaxCompnts || !TotCompnts) {
        fprintf(stderr, "CalcCom(): Null arguments!\n");
        return FALSE;
    }
    
    if (dir == COL) {
        l = pM->m;
    } else if (dir == ROW) {
        l = pM->n;
    } else {
        fprintf(stderr, "CalcCom(): Unknown direction!\n");
        return FALSE;
    }
    
    P = pM->NrProcs;
    
    Nprocs = (int *)malloc(l*sizeof(int));
    if (Nprocs == NULL) {
        fprintf(stderr, "CalcCom(): Not enough memory!\n");
        return FALSE;
    }
    
    if (!InitNprocs(pM, dir, Nprocs)) {
        fprintf(stderr, "CalcCom(): Unable to initialise processor array!\n");
        return FALSE;
    }

    if (X == NULL) { 
        /* Compute matrix-based communication volume */
        tots= 0;
        for (j=0; j<l; j++)
            tots += MAX(Nprocs[j] - 1,0); 
        *ComVol = tots;
        *MaxOut = 0;
        *MaxIn = 0;
        *MaxCompnts = 0;
        *TotCompnts = 0;
        free(Nprocs);
        return TRUE;
    }
    
    total = 0;
    for (j=0; j<l; j++)
        total += Nprocs[j];
    procindex = (int *)malloc(total*sizeof(int));
    procstart = (long *)malloc((l+1)*sizeof(long));
    if (procindex == NULL || procstart == NULL) {
        fprintf(stderr, "CalcCom(): Not enough memory!\n");
        return FALSE;
    }
    if (!InitProcindex(pM, dir, Nprocs, procstart, procindex)) {
        fprintf(stderr, "CalcCom(): Cannot initialise processor index!\n");
        return FALSE;
    }

    Ns = (long *)malloc(P*sizeof(long));
    Nr = (long *)malloc(P*sizeof(long));
    Nv = (long *)malloc(P*sizeof(long));
    if (Ns == NULL || Nr == NULL || Nv == NULL) {
        fprintf(stderr, "CalcCom(): Not enough memory!\n");
        return FALSE;
    }

    /* Compute the number of components sent, received, owned 
       by the P processors */
    for (q=0; q<P; q++) {
        Ns[q] = 0;
        Nr[q] = 0;
        Nv[q] = 0;
    } 
    
    for (j=0; j<l; j++) {
        q = X[j];
        if (q >= 0 && q < P) {
            /* Disassign for a moment and reassign to update the statistics */
            X[j] = -1;
            if (!AssignColumnToProc(X, procstart, procindex, Ns, Nr, Nv, NULL, j, q)) {
                fprintf(stderr, "CalcCom(): Unable to assign column!\n");
                return FALSE;
            }
        }
    }

    /* Compute the maximum and total number of components
       sent, received, owned */
    maxs= 0; maxr= 0; maxv= 0;
    tots= 0; totr= 0; totv=0;
    for (q=0; q<P; q++) {
        if (Ns[q] > maxs)
            maxs = Ns[q];
        if (Nr[q] > maxr)
            maxr= Nr[q];
        if (Nv[q] > maxv)
            maxv= Nv[q];
        tots += Ns[q];
        totr += Nr[q];
        totv += Nv[q];
    } 
    *MaxOut= maxs;
    *MaxIn= maxr;
    *MaxCompnts= maxv;
    if (tots != totr) {
        fprintf(stderr, "CalcCom(): Total sent != total received!\n");
        return FALSE;
    }
    *ComVol= totr;
    *TotCompnts= totv;
    
    free(Nv); free(Nr); free(Ns); 
    free(procstart); free(procindex);
    free(Nprocs);
    
    return TRUE;
} /* end CalcCom */


int PrintCom(int P, long l, int dir, long ComVol, long MaxOut, long MaxIn, 
              long MaxCompnts, long TotCompnts) {

    /* This function prints the total and maximum communication volume
       and total and maximum number of vector components
       for a given distributed sparse matrix A and a distributed vector X
       as computed by CalcCom. Note that these numbers are for only
       one direction (vector). 

       The function also prints scaled numbers (between brackets),
       which have been normalised by the length l of the vector.

       Input: P is the number of processors,
              l is the length of the vector,
              dir = ROW (for distribution of v) or dir = COL (for u)
              ComVol is total communication volume for direction dir,
              MaxOut is maximum number of data words sent by a processor,
              MaxIn is maximum number of data words received by a processor,
              MaxCompnts is maximum number of vector components of a processor.
              TotCompnts is total number of owned vector components
                            (can be < l).
    */

    double ComVolSc, MaxOutSc, MaxInSc;

    /* Scaled values, normalised by length of vector */
    ComVolSc = (double) ComVol  / (double) l;
    MaxOutSc = (double) MaxOut / (double) l;
    MaxInSc  = (double) MaxIn  / (double) l;
  
    printf("Communication for vector"
           " %c (after vector distribution):\n", dir == ROW ? 'v' : 'u');
    printf("  vol           : %ld (%g)\n", ComVol, ComVolSc);
    printf("  avg = vol/P   : %g (%g)\n", (double)ComVol / (double)P,
                                         ComVolSc / (double)P);
    if (dir == ROW) {
        printf("  max sent      : %ld (%g)\n", MaxOut, MaxOutSc);
        printf("  max recd      : %ld (%g)\n", MaxIn, MaxInSc);
    } else {
        /* roles of sends and receives are interchanged */
        printf("  max sent      : %ld (%g)\n", MaxIn, MaxInSc);
        printf("  max recd      : %ld (%g)\n", MaxOut, MaxOutSc);
    }
    printf("Nr components for vector "
           "%c (after vector distribution):\n", dir == ROW ? 'v' : 'u');
    printf("  avg           : %g (%g) \n", (double) TotCompnts / (double) P,
           (double) TotCompnts / (double) (P*l));
    printf("  max           : %ld (%g)\n", MaxCompnts,
           (double) MaxCompnts /(double) l);
    
    return TRUE;
} /* end PrintCom */


int CalcLocalLowerBound(const struct sparsematrix *pM, int dir, long *LB, int *Pact) {

    /* This function computes the local lower bound on the maximum communication
       volume per processor for a distributed sparse matrix pM-> 
       This bound for a vector distribution is solely based on the matrix distribution.
       Input:  distributed matrix A,
               dir = ROW (for distribution of v) or dir = COL (for u).
       Output: LB = local lower bound (max over all processors),
               Pact = number of active processors, i.e., processors
               that will have to communicate. 0 <= Pact <= pM->NrNprocs.
    */
   
    int P, q, r, Pactive;
    long l=0, t, i, j, local, Ls, Lr, *Index = NULL;

    /* Arrays of length pM->n */
    int *Nprocs;     /* Nprocs[j] is the number of processors that occur in column j */    
                       
    int *occur;      /* occur[j]= q means processor q occurs in column j
                       and is its current highest-numbered processor.
                       occur[j]= -1 means no processor occurs yet in column j */
    
    /* Array of length P+1 */
    long *Available; /* Available[r] = the number of columns with r processors
                        that contain the current processor, 0 <= r <= P. */
    
    if (!pM || !LB || !Pact) {
        fprintf(stderr, "CalcLocalLowerBound(): Null arguments!\n");
        return FALSE;
    }
    
    if (dir == COL) {
        l = pM->m;
        Index = pM->i;        
    } else if (dir == ROW) {
        l = pM->n;
        Index = pM->j;
    } else {
        fprintf(stderr, "CalcLocalLowerBound(): Unknown direction!\n");
        return FALSE;
    }
    P = pM->NrProcs;


    Nprocs = (int *)malloc(l*sizeof(int));
    occur = (int *)malloc(l*sizeof(int));
    Available = (long *)malloc((P+1)*sizeof(long));
    if (Nprocs == NULL || occur == NULL || Available == NULL) {
        fprintf(stderr, "CalcLocalLowerBound(): Not enough memory!\n");
        return FALSE;
    }

    if (!InitNprocs(pM, dir, Nprocs)) {
        fprintf(stderr, "CalcLocalLowerBound(): Unable to initialise processor array!\n");
        return FALSE;
    }

    /* Initialise occur */
    for (j=0; j<l; j++)
        occur[j] = -1;

#ifdef INFO2
    printf("Local bounds for vector %c:\n", dir == ROW ? 'v' : 'u');
#endif
    local = 0;
    Pactive = 0;
    for (q=0; q<P; q++) {
        /*### Determine bound for processor q ###*/
    
        /* Initialise Available[r] at number of columns with r processors
           that contain processor q */
        for (r=0; r<=P; r++)
            Available[r] = 0;

        for (t=pM->Pstart[q]; t<pM->Pstart[q+1]; t++) {
            j= Index[t];
            if (occur[j] != q) {
                occur[j] = q;
                Available[Nprocs[j]]++;
            }
        }

        /* Compute egoistic lower bound for processor q */
        Ls= 0; /* number of sends for processor q */
        Lr= 0; /* number of receives */
        for (r=2; r<=P; r++)
            Lr += Available[r]; /* Initialise Lr at total number of columns
                                   with Nprocs >= 2 that contain processor q */
             
        if (Lr > 0)
            Pactive++;
      
        for (r=2; r<=P; r++) {        
            if (Ls + (r-1)*Available[r] <=  Lr - Available[r]) {
                /* Grab all available columns with r processors */
                Ls += (r-1)*Available[r];
                Lr -= Available[r];
            } else {
                /* Grab i columns with i as large as possible, 
                   such that new Ls <= Lr */
                i = (Lr-Ls)/r; /* 0 <= i <= Available[r] */
                Ls += (r-1)*i;
                Lr -= i;
                break;
            }
        }

       if (Ls > Lr) {
            fprintf(stderr, "CalcLocalLowerBound(): Ls > Lr!\n");
            return FALSE;
       }
#ifdef INFO2
        printf(" Proc[%d]: %ld\n", q, Lr); /* lower bound = Lr */
#endif
 
        if (Lr > local)
            local = Lr;
    }
    
    *LB = local;
    *Pact = Pactive;
    
    free(Available);
    free(occur);
    free(Nprocs);

    return TRUE;
} /* end CalcLocalLowerBound */


void PrintVecStatistics(int P, long *Ns, long *Nr, long *Nv) {

    /* This function prints for each processor the send, receive values
       and the number of components.
       Ns, Nr = number of sends/recvs;
       Nv = number of components. */

    int q;
    
    if (!Ns || !Nr || !Nv) return;

    for (q=0; q<P; q++)
        printf("  Proc[%d]: Ns = %ld, Nr = %ld, Nv = %ld\n",
               q, Ns[q], Nr[q], Nv[q]);

} /* end PrintVecStatistics */


int WriteVector(const long int *X, const char base, const char *name, long l, int P, FILE *fp, const struct opts *pOptions) {

    /* Base vector-writing function. Automatically adapts base of array to write.
       Assumes that we are writing a distribution vector when P>0. P=0 enables
       simple writing of a vector.
       l is the length of the vector to be written (X).
       base is the base value of X (usually 0).
       name is used only when EMM output is enabled.
       name is NULL if vector is the main object of the EMM file.
    */

    long t;

    if (!X || !fp) {
        fprintf(stderr, "WriteVectorDistribution(): Null arguments!\n");
        return FALSE;
    }

    if (ferror(fp)) {
        fprintf(stderr, "WriteVectorDistribution(): Unable to write to stream!\n");
        return FALSE;
    }

    if (pOptions->OutputFormat == OutputEMM) {
        if( P == 0 ) {
            if (name==NULL)
                fprintf(fp, "%%%%Extended-MatrixMarket vector array integer general original\n");
            else
                fprintf(fp, "%%%%%s vector array integer general original\n", name);
        } else {
            if (name==NULL)
                fprintf(fp, "%%%%Extended-MatrixMarket distributed-vector array integer general global\n");
            else
                fprintf(fp, "%%%%%s distributed-vector array integer general global\n", name);
        }
    }

    if( P > 0 ) {
        fprintf(fp, "%ld %d\n", l, P);
        for (t = 0; t < l; t++)
            fprintf(fp, "%ld %ld\n", t+1, X[t]+(1-base));
    } else {
        if( pOptions->OutputFormat == OutputEMM )
            fprintf(fp, "%ld\n", l);
        for (t = 0; t < l; t++)
            fprintf(fp, "%ld\n", X[t]+(1-base));
    }
    return TRUE;
} /* end WriteVector */

int WriteVectorDistribution(const long int *X, const char *name, long l, int P, FILE *fp, const struct opts *pOptions) {

    /* This function writes the owners of the vector components to a file.
       Input: l is the number of vector components.
              P is the number of processors.
              X is an array of length l, containing the owners 
                   of the vector components. 
              The input numbering is 0-based: 0 <= X[i] < P, for 0 <= i < l.
       The output file format is:
              1 line with l , P
              1 line per vector component, giving the index i and the owner X[i].
              The output numbering is 1-based: 1 <= X[i] <= P, for 1 <= i <= l.
              This numbering is consistent with the Matrix Market format. */

    return WriteVector( X, 0, name, l, P, fp, pOptions );
} /* end WriteVectorDistribution */

int WriteVectorCollection(long int **X, const char *name, const long i, const long *j, FILE *fp) {
        int l, k;
        if(name == NULL)
            fprintf(fp, "%%%%Extended-MatrixMarket vector-collection array integer general original\n");
        else
            fprintf(fp, "%%%%%s vector-collection array integer general original\n", name);
        fprintf(fp, "%ld\n", i);
        for(l=0; l<i; l++) {
            fprintf(fp, "%ld\n", j[l]);
            for(k=0; k<j[l]; k++)
                fprintf(fp, "%ld\n", (X[l][k]+1)); /* base 1 */
        }
        return TRUE;
} /* end WriteVectorCollection */

