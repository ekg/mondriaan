#include "DistributeVecOrig.h"

long DistributeVecOrig(const struct sparsematrix *pM, long int *X, int dir, const struct opts *pOptions) {
    /*  This function determines a vector distribution X by the
        vector partitioning algorithm, Algorithm 2,
        in the original SIAM Review paper [1],
        which tries to minimise the communication cost, defined as
        the maximum over all processors of max(sends,recvs).

        [1] B. Vastenhouw and R. H. Bisseling,
            A two-dimensional data distribution method
            for parallel sparse matrix-vector multiplication",
            SIAM Review, Vol. 47, No. 1 (2005) pp. 67-95.

        Input:  initialised A,
                uninitialised array X of length pM->n,
                dir is the direction:
                dir = ROW means that the vector is the input vector v,
                      COL means that the vector is the output vector u = Av.

        Output: vector distribution X, which is an array of length pM->n
                with values between 0 and P-1, where P= pM->NrProcs.
                The function returns maxcom = max(sends,recvs). */

    long j, k, s, t, total, l=0, lo, hi, min, Nstotal, Nrtotal, maxcom,
         *ProcHistogram,  *iv = NULL;
    int P, npj, q, r, proc;

    /* Arrays of length pM->n */
    int *Nprocs;    /* Nprocs[j] is the number of processors
                       that occur in column j */
    long *Nprocs2;  /* Same array, but sorted in decreasing order.
                       Type is long, enabling use of our quicksort for longs */
 
    /* Communication matrix, stored in compressed column storage (CCS).
       Includes all columns of the original matrix. */
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
    long *Sums;     /* Sums[q] = the current number of communications
                       (sends or receives) of processor q, taking inevitable
                       communications into account from the start */
    if (!pM || !X || !pOptions) {
        fprintf(stderr, "DistributeVecOrig(): Null parameters!\n");
        return -1;
    }
    
    if (dir == COL) {
        l = pM->m;
    } else if (dir == ROW) {
        l = pM->n;
    } else {
        fprintf(stderr, "DistributeVecOrig(): Unknown direction!\n");
        return -1;
    }
    P = pM->NrProcs;
    
    /* Initialise processor arrays: */
    ProcHistogram = (long *) malloc((P+1) * sizeof(long));
    Nv = (long *) malloc(P * sizeof(long));
    Sums = (long *) malloc(P * sizeof(long));
    Ns = (long *) malloc(P * sizeof(long));
    Nr = (long *) malloc(P * sizeof(long));
    if (ProcHistogram == NULL || Nv == NULL || 
         Sums == NULL || Ns == NULL || Nr == NULL) {
        fprintf(stderr, "DistributeVecOrig(): Not enough memory!\n");
        return -1;
    }
  
    for (q = 0; q < P; q++) {
        Ns[q] = 0;
        Nr[q] = 0;
        Nv[q] = 0;
    }
    
    /* All columns are unowned at the start */
    for (j=0; j<l; j++)
        X[j] = -1;

    /*## STEP 1: 
      ##  Generate communication matrix and histogram, and initialise sums: ##*/

    /* Initialise CCS data structure for communication matrix  */
    Nprocs = (int *)malloc(l*sizeof(int));
    Nprocs2 = (long *)malloc(l*sizeof(long));
    procstart = (long *)malloc((l+1)*sizeof(long));
    if (Nprocs == NULL || Nprocs2 == NULL || procstart == NULL) {
        fprintf(stderr, "DistributeVecOrig(): Not enough memory!\n");
        return -1;
    }

    if (!InitNprocs(pM, dir, Nprocs)) {
        fprintf(stderr, "DistributeVecOrig(): Unable to initialise processor array!\n");
        return -1;
    }

    total = 0;
    for (j=0; j<l; j++)
        total += Nprocs[j];
    procindex = (int *)malloc(total*sizeof(int));
    if (procindex == NULL) {
        fprintf(stderr, "DistributeVecOrig(): Not enough memory!\n");
        return -1;
    }

    if (!InitProcindex(pM, dir, Nprocs, procstart, procindex)) {
        fprintf(stderr, "DistributeVecOrig(): Unable to initialise processor index!\n");
        return -1;
    }

    if (!GenerateHistogram(Nprocs, l, 0, P, ProcHistogram)) {
        fprintf(stderr, "DistributeVecOrig(): Unable to create histogram!\n");
        return -1;
    }
    if (!InitSums(l, P, procstart, procindex, Sums)) {
        fprintf(stderr, "DistributeVecOrig(): Unable to initialise sums!\n");
        return -1;
    }
  
    /* Copy Nprocs into Nprocs2 and sort Nprocs2 in decreasing order */
    for (j=0; j<l; j++)
        Nprocs2[j] = Nprocs[j];
    iv = QSort(Nprocs2, l); /* iv stores the original indices */
    
    if (iv == NULL) {
        fprintf(stderr, "DistributeVecOrig(): Sorting failed!\n");
        return -1;
    }
  
    /*## STEP 2: 
      ##  Assign columns with one processor :##*/
    lo = l - ProcHistogram[0] - ProcHistogram[1];
    hi = l - ProcHistogram[0];
    for (t = lo; t < hi; t++) {      
        j = iv[t]; 
  
        /* Check number of processors */
        if (Nprocs2[t] != 1) {
            fprintf(stderr, "DistributeVecOrig(): Internal error step 2: Nprocs != 1!\n");
            return -1;
        }
  
        /* Assign column j to unique owner */
        q = procindex[procstart[j]];
        X[j] = q;
        Nv[q]++;
    }
  
  
    /*## STEP 3:
      ##  Assign columns with more than two processors: ##*/
    if (P > 2) {
        lo = 0;
        hi = l - ProcHistogram[0] - ProcHistogram[1] - ProcHistogram[2];
        if (pOptions->VectorPartition_Step3 == VecIncrease) {
            /* Reverse the order in arrays iv and Nprocs2 */
            for (k = 0; k < (hi-lo) / 2; k++) {
                SwapLong(iv, lo+k, hi-k-1);
                SwapLong(Nprocs2, lo+k, hi-k-1);
            }
        } else if (pOptions->VectorPartition_Step3 == VecRandom) {
            RandomPermute(iv, Nprocs2, NULL, NULL, lo, hi-1);
        } else if (pOptions->VectorPartition_Step3 == VecDecrease) {
           ; /* do nothing, arrays are already in decreasing order */
        } else  {
            fprintf(stderr, "DistributeVecOrig(): Unknown order in step 3!\n");
            return -1;
        }
        
        for (t = lo; t < hi; t++) {
            j = iv[t]; 
            npj = Nprocs2[t];
     
            /* Check number of processors */
            if (npj < 3 || npj > P) {
                fprintf(stderr, "DistributeVecOrig(): Internal error step 3: Nprocs < 3 or Nprocs > P!\n");
                return -1;
            }
          
            /* Find a processor with the minimum sum */
            proc = procindex[procstart[j]];
            min = Sums[proc];
            for (s=procstart[j]+1; s<procstart[j+1]; s++) {
                q = procindex[s];
                if (Sums[q] < min) {
                    proc = q;
                    min = Sums[q];      
                }
            }

            /* If there is more than one processor with the minimum sum
               we could look at a secondary objective like best balance.
               However, we do not use secondary objectives here. */
        
            /* Assign column j to this processor: */
            if (!AssignColumnToProc(X, procstart, procindex, Ns, Nr, Nv, Sums, j, proc)) {
                fprintf(stderr, "DistributeVecOrig(): Unable to assign column!\n"); 
                return -1;
            }
        }
    }          
  
  
    /*## STEP 4:
      ##  Assign columns with two processors: ##*/
    if (P > 1) {
        lo = l - ProcHistogram[0] - ProcHistogram[1] - ProcHistogram[2];
        hi = l - ProcHistogram[0] - ProcHistogram[1];
  
        for (t = lo; t < hi; t++) {
            j = iv[t]; 
            
            /* Check number of processors */
            if (Nprocs2[t] != 2) {
                fprintf(stderr, "DistributeVecOrig(): Internal error step 4: Nprocs != 2!\n");
                return -1;
            }
            
            /* Find the direction q->r or r->q with the minimum
               value Ns(sender) + Nr(receiver).
               This is the least loaded send-receive direction. */
            q = procindex[procstart[j]];
            r = procindex[procstart[j]+1];

            if (Ns[q] + Nr[r] < Ns[r] + Nr[q]) {
                proc = q;
            } else if (Ns[q] + Nr[r] > Ns[r] + Nr[q]) {
                proc = r;
            } else if (Random1(0,1) == 0) { /* random tie-breaking for
                                              Ns[q]+Nr[r] = Ns[r]+Nr[q] */
                proc = q;
            } else
                proc = r;
      
            /* Assign column j to this processor */
            if (!AssignColumnToProc(X, procstart, procindex, Ns, Nr, Nv, Sums, j, proc)) {
                fprintf(stderr, "DistributeVecOrig(): Unable to assign column!\n"); 
                return -1;
            }
        }
    }
    
  
    /*## STEP 5:
      ##  Assign columns with no processors: ##*/
    if (!AssignRemainingColumns(l, P, X, Nv)) {
        fprintf(stderr, "DistributeVecOrig(): Unable to assign remaining columns!\n");
        return -1;
    }

    /* Compute the communication volume and cost */
    Nstotal = 0;
    Nrtotal = 0;
    maxcom = 0;
    for (q=0; q<P; q++) {
        Nstotal += Ns[q];
        Nrtotal += Nr[q];
        if (Ns[q] > maxcom)
            maxcom = Ns[q];
        if (Nr[q] > maxcom)
            maxcom = Nr[q];
    }
    if (Nstotal != Nrtotal) {
        fprintf(stderr, "DistributeVecOrig(): Total sends != total recvs!\n");
        return -1;
    }
 
    /* Clean up memory */
    
    free(iv);
    free(procindex); /* no check needed */
    free(procstart);
    free(Nprocs2);
    free(Nprocs);
    free(Sums);
    free(Nv);
    free(Nr);
    free(Ns);
    free(ProcHistogram);
  
    return maxcom;
} /* end DistributeVecOrig */


int InitSums(long l, int P, long *procstart, int *procindex, long *Sums) {
    /* This function counts the inevitable communication operations
       (sends or receives) that are known at the start of the
       vector distribution. If a column is owned by more than one processor,
       each processor has to send or receive something. 
       Input: l is the number of vector components.
              P is the number of processors.
              procindex, contains the processor numbers occurring in
                         the matrix columns. The processor numbers for one
                         column j are stored consecutively.
              procstart, where procstart[j] is the position in array procindex
                         of the first processor number of column j, 0 <= j < pM->n.
                         procstart[pM->n] = total number of processor numbers
                         of all the columns together.
       Output: Sums[q] = the number of inevitable communications of processor q.
    */

    long j, q, r;
    
    if (!procstart || !procindex || !Sums) {
        fprintf(stderr, "InitSums(): Null arguments!\n");
        return FALSE;
    }

    for (q = 0; q < P; q++)
        Sums[q] = 0;
  
    for (j = 0; j < l; j++) {
        if (procstart[j+1] - procstart[j] >= 2)
            for (r = procstart[j]; r < procstart[j+1]; r++)
                Sums[procindex[r]]++;
    }

    return TRUE;
} /* end InitSums */

