#include "DistributeVecOrigEq.h"

int InitIntersection(long l, int P,
                     long *procstartU, int *procindexU,
                     long *procstartV, int *procindexV, int *NprocsUV);

int FindProcLowestSum(long *procstart, int *procindex, long j, int npj,
                       long *SumsU, long *SumsV);

long DistributeVecOrigEq(const struct sparsematrix *pM, long int *U, long int *V, const struct opts *pOptions) {
    /*  This function determines a vector distribution U = V by the
        vector partitioning algorithm for equal input and output distribution
        in the original SIAM Review paper [1],
        which tries to minimise the communication cost, defined as
        the maximum over all processors of max(sends,recvs).

        The function has not been fully optimised, since in general
        there will not be much choice in assigning the vector components.

        [1] B. Vastenhouw and R. H. Bisseling,
            A two-dimensional data distribution method
            for parallel sparse matrix-vector multiplication",
            SIAM Review, Vol. 47, No. 1 (2005) pp. 67-95.

        Input:  initialised A,
                uninitialised arrays U, V of length pM->m = pM->n,
                representing output and input distributions for u = Av.

        Output: vector distributions U, V, with U=V, 
                which are arrays of length pM->m
                with values between 0 and P-1, where P= pM->NrProcs.
                The function returns maxcom = max(sends,recvs) and a value < 0 on failure. */
 
    int P, Pact, q, r, proc;
    long l, j, total, LB, actbound, bestbound, 
         ComVol, MaxOut, MaxIn, MaxCompnts, TotCompnts, maxcom,
         *ProcHistogramU, *ProcHistogramV , *ProcHistogramUV;
    
    /* Arrays of length pM->n for vector v */
    int *NprocsV;   /* NprocsV[j] is the number of processors
                       that occur in column j */
 
    /* Communication matrix for vector v, stored in compressed column storage (CCS).
       Includes all columns of the original matrix. */
    int *procindexV; /* array of length at most V+pM->n, where V is the
                       communication volume for the vector v.
                       procindexV[procstartV[j]..procstartV[j+1]-1]
                       contains the processor numbers of the processors
                       that occur in column j */
    long *procstartV; /* array of length pM->n+1 containing the starts of the
                       columns of the communication matrix */

    /* Arrays of length P for vector v */
    long *NsV;       /* NsV[q] = the current number of vector components
                       sent by processor q to other processors */
    long *NrV;       /* NrV[q] = the current number of vector components
                       received by processor q from other processors */
    long *NvV;       /* NvV[q] = the current number of vector components
                       owned by processor q */
    long *SumsV;     /* SumsV[q] = the current number of communications
                       (sends or receives) of processor q, taking inevitable
                       communications into account from the start */
    
    /* The same for vector u */
                       
    int *NprocsU;   /* NprocsU[j] is the number of processors
                       that occur in row j */
    int *procindexU; /* procindexU[procstartU[j]..procstartU[j+1]-1]
                       contains the processor numbers of the processors
                       that occur in row j */
    long *procstartU; 
    long *NsU, *NrU, *NvU, *SumsU;     
        
    int *NprocsUV;  /* NprocsUV[j] is the number of processors that occur 
                       both in row j and in column j */
    
    if (!pM || !U || !V || !pOptions) {
        fprintf(stderr, "DistributeVecOrigEq(): Null parameters!\n");
        return -1;
    }
    
    /*## Check that the matrix is square: ##*/
    if (pM->m != pM->n) {
        fprintf(stderr, "DistributeVecOrigEq(): U and V cannot be distributed the same way because the matrix is not square!\n");
        return -1;
    }
  
    l = pM->n;
    P = pM->NrProcs;
 
    /*## Allocate processor arrays: ##*/
    ProcHistogramU = (long *) malloc((P+1) * sizeof(long));
    ProcHistogramV = (long *) malloc((P+1) * sizeof(long));
    ProcHistogramUV = (long *) malloc((P+1) * sizeof(long));
    NsU = (long *) malloc(P * sizeof(long));
    NsV = (long *) malloc(P * sizeof(long));
    NrU = (long *) malloc(P * sizeof(long));
    NrV = (long *) malloc(P * sizeof(long));
    NvU = (long *) malloc(P * sizeof(long));
    NvV = (long *) malloc(P * sizeof(long));
    SumsU = (long *) malloc(P * sizeof(long));
    SumsV = (long *) malloc(P * sizeof(long));
    if (ProcHistogramU == NULL || ProcHistogramV == NULL || 
         ProcHistogramUV == NULL ||  
         NsU == NULL || NrU == NULL || NvU == NULL || SumsU == NULL ||
         NsV == NULL || NrV == NULL || NvV == NULL || SumsV == NULL) {
        fprintf(stderr, "DistributeVecOrigEq(): Not enough memory for processor arrays!\n");
        return -1;
    }
    
    /*## Clear memory: ##*/
    for (q = 0; q < P; q++) {
        NsU[q] = NrU[q] = NvU[q] = 0;
        NsV[q] = NrV[q] = NvV[q] = 0;
    }
  
    /*## STEP 1: 
      ##  Generate communication matrices, histograms and
          initialise intersections and sums: ##*/
      
    /* Initialise CRS data structure for communication matrix of U */
    NprocsU = (int *)malloc(l*sizeof(int));
    procstartU = (long *)malloc((l+1)*sizeof(long));
    if (NprocsU == NULL || procstartU == NULL) {
        fprintf(stderr, "DistributeVecOrigEq(): Not enough memory!\n");
        return -1;
    }

    if (!InitNprocs(pM, COL, NprocsU)) {
        fprintf(stderr, "DistributeVecOrigEq(): Unable to initialise processor array!\n");
        return -1;
    }

    total = 0;
    for (j=0; j<l; j++)
        total += NprocsU[j];
    procindexU = (int *)malloc(total*sizeof(int));
    if (procindexU == NULL) {
        fprintf(stderr, "DistributeVecOrigEq(): Not enough memory!\n");
        return -1;
    }

    if (!InitProcindex(pM, COL, NprocsU, procstartU, procindexU)) {
        fprintf(stderr, "DistributeVecOrigEq(): Unable to initialise processor index!\n");
        return -1;
    }
    
    /* Initialise CCS data structure for communication matrix of V */
    NprocsV = (int *)malloc(l*sizeof(int));
    procstartV = (long *)malloc((l+1)*sizeof(long));
    if (NprocsV == NULL || procstartV == NULL) {
        fprintf(stderr, "DistributeVecOrigEq(): Not enough memory!\n");
        return -1;
    }

    if (!InitNprocs(pM, ROW, NprocsV)) {
        fprintf(stderr, "DistributeVecOrigEq(): Unable to initialise processor array!\n");
        return -1;
    }

    total = 0;
    for (j=0; j<l; j++)
        total += NprocsV[j];
    procindexV = (int *)malloc(total*sizeof(int));
    if (procindexV == NULL) {
        fprintf(stderr, "DistributeVecOrigEq(): Not enough memory!\n");
        return -1;
    }

    if (!InitProcindex(pM, ROW, NprocsV, procstartV, procindexV)) {
        fprintf(stderr, "DistributeVecOrigEq(): Unable to initialise processor index!\n"); 
        return -1;
    }

    NprocsUV = (int *) malloc(l * sizeof(int));
    if (NprocsUV == NULL) {
        fprintf(stderr, "DistributeVecOrigEq(): Not enough memory!\n");
        return -1;
    }
    if (!InitIntersection(l, P, procstartU, procindexU, procstartV, procindexV, NprocsUV)) {
        fprintf(stderr, "DistributeVecOrigEq(): Unable to initialise intersection!\n");
        return -1;
    }
    
    if (!GenerateHistogram(NprocsU, l, 0, P, ProcHistogramU)) {
        fprintf(stderr, "DistributeVecOrigEq(): Unable to create histogram!\n");
        return -1;
    }
    if (!GenerateHistogram(NprocsV, l, 0, P, ProcHistogramV)) {
        fprintf(stderr, "DistributeVecOrigEq(): Unable to create histogram!\n");
        return -1;
    }
    if (!GenerateHistogram(NprocsUV, l, 0, P, ProcHistogramUV)) {
        fprintf(stderr, "DistributeVecOrigEq(): Unable to create histogram!\n");
        return -1;
    }
  
    /* Calculate and print communication statistics for vector v */ 
    if (!CalcCom(pM, NULL, ROW, &ComVol, &MaxOut, &MaxIn, &MaxCompnts, &TotCompnts)) {
        fprintf(stderr, "DistributeVecOrigEq(): Unable to calculate communication!\n");
        return -1;
    }
    if (!CalcLocalLowerBound(pM, ROW, &LB, &Pact)) {
        fprintf(stderr, "DistributeVecOrigEq(): Unable to calculate communication lower bound!\n");
        return -1;
    }
    
    if (Pact > 0)
        actbound= (ComVol%Pact == 0 ? ComVol/Pact : ComVol/Pact +1);
    else
        actbound= 0;
    
    bestbound= MAX(LB,actbound);

#ifdef INFO
    printf("\nNr matrix columns owned by p processors:\n");
    PrintHistogram(0, P, ProcHistogramV);
    printf("Communication for vector v (after matrix distribution):\n");
    printf("  vol           : %ld \n", ComVol);
    printf("  active procs  : %d \n", Pact);
    printf("  best bound    : %ld \n", bestbound);
#endif
#ifdef INFO2
    printf("    active bound: %ld \n", actbound);
    printf("    local bound : %ld \n", LB);
#endif 
  
    /* Calculate and print communication statistics for vector u */  
    if (!CalcCom(pM, NULL, COL, &ComVol, &MaxOut, &MaxIn, &MaxCompnts, &TotCompnts)) {
        fprintf(stderr, "DistributeVecOrigEq(): Unable to calculate communication!\n");
        return -1;
    }
    if (!CalcLocalLowerBound(pM, COL, &LB, &Pact)) {
        fprintf(stderr, "DistributeVecOrigEq(): Unable to calculate communication lower bound!\n");
        return -1;
    }
    
    if (Pact > 0) 
        actbound= (ComVol%Pact == 0 ? ComVol/Pact : ComVol/Pact +1); 
    else 
        actbound= 0;
    
    bestbound= MAX(LB,actbound); 

#ifdef INFO
    printf("\nNr matrix rows owned by p processors:\n");
    PrintHistogram(0, P, ProcHistogramU); 
    printf("Communication for vector u (after matrix distribution):\n");
    printf("  vol           : %ld \n", ComVol); 
    printf("  active procs  : %d \n", Pact); 
    printf("  best bound    : %ld \n", bestbound);
#endif 
#ifdef INFO2 
    printf("    active bound: %ld \n", actbound); 
    printf("    local bound : %ld \n", LB); 
#endif 

/* Print statistics for intersection of vectors u and v */
#ifdef INFO 
    printf("\nNr rows/columns with p processors in the intersection\n");
    printf("                           (after matrix distribution):\n");
    PrintHistogram(0, P, ProcHistogramUV); 
    printf("\n");
#endif
  
    /* Initialise sums and vector distributions */
    if (!InitSums(l, P, procstartU, procindexU, SumsU)) {
        fprintf(stderr, "DistributeVecOrigEq(): Unable to initialise sums!\n");
        return -1;
    }
    if (!InitSums(l, P, procstartV, procindexV, SumsV)) {
        fprintf(stderr, "DistributeVecOrigEq(): Unable to initialise sums!\n");
        return -1;
    }
    
    /* All columns are unowned at the start */
    for (j=0; j<l; j++) {
        U[j] = -1;
        V[j] = -1;
    }
  
    /*## STEP 2: 
      ##  Assign vector components with one processor in the intersection : ##*/ 
    for (j=0; j<l; j++) {
        if (NprocsUV[j] == 1) {
            q = procindexU[procstartU[j]];
  
            /* Assign U[j] and V[j] to processor q: */
            if (!AssignColumnToProc(U, procstartU, procindexU, NsU, NrU, NvU, SumsU, j, q)) {
                fprintf(stderr, "DistributeVecOrigEq(): Unable to assign column!\n");
                return -1;
            }
            if (!AssignColumnToProc(V, procstartV, procindexV, NsV, NrV, NvV, SumsV, j, q)) {
                fprintf(stderr, "DistributeVecOrigEq(): Unable to assign column!\n");
                return -1;
            }
        }
    }
  
    /*## STEP 3: 
      ##  Assign vector components with more than one processor in the 
      ##  intersection :##*/ 
    for (j=0; j<l; j++) {
        if (NprocsUV[j] > 1) {
      
            /* Find a processor with the lowest total sum: */
            q = FindProcLowestSum(procstartU, procindexU, j, NprocsUV[j], SumsU, SumsV);
 
            /* If there is more than one processor having the minimum sum
               we could look at a secondary objective. 
               However, we do not use secondary objectives here. */
    
            /* Assign U[j] and V[j] to this processor: */
            if (!AssignColumnToProc(U, procstartU, procindexU, NsU, NrU, NvU, SumsU, j, q)) {
                fprintf(stderr, "DistributeVecOrigEq(): Unable to assign column!\n");
                return -1;
            }
            if (!AssignColumnToProc(V, procstartV, procindexV, NsV, NrV, NvV, SumsV, j, q)) {
                fprintf(stderr, "DistributeVecOrigEq(): Unable to assign column!\n");
                return -1;
            }
        }
    }   
  
  
    /*## STEP 4: 
      ##  Assign vector components with no processors in the
          intersection of the processor lists of u and v: ##*/ 
  
    /*## Step 4A:
      ##  Assign vector components j with NprocsU[j]=0 XOR NprocsV[j]=0: ##*/
    for (j=0; j<l; j++) {
        if (NprocsUV[j] == 0 && (NprocsU[j] == 0 || NprocsV[j] == 0) && (NprocsU[j] != 0 || NprocsV[j] != 0) ) {

            /* Find a processor with the lowest total sum: */
            if (NprocsU[j] != 0 && NprocsV[j] == 0) 
                q = FindProcLowestSum(procstartU, procindexU, j,NprocsU[j], SumsU,SumsV);
            else if (NprocsU[j] == 0 && NprocsV[j] != 0) 
                q = FindProcLowestSum(procstartV, procindexV, j, NprocsV[j], SumsU,SumsV);
      
            /* Assign U[j] and V[j] to this processor: */
            if (!AssignColumnToProc(U, procstartU, procindexU, NsU, NrU, NvU, SumsU, j, q)) {
                fprintf(stderr, "DistributeVecOrigEq(): Unable to assign column!\n");
                return -1;
            }
            if (!AssignColumnToProc(V, procstartV, procindexV, NsV, NrV, NvV, SumsV, j, q)) {
                fprintf(stderr, "DistributeVecOrigEq(): Unable to assign column!\n");
                return -1;
            }
        }
    }
  
    
    /*## Step 4B:
      ##  Assign vector components j with NprocsU[j]!=0 and NprocsV[j]!=0: ##*/
    for (j=0; j<l; j++) {
        if (NprocsUV[j] == 0 && NprocsU[j] != 0 && NprocsV[j] != 0) {

            /* Find a processor with the lowest total sum: */
            q = FindProcLowestSum (procstartU, procindexU, j, NprocsU[j], SumsU,SumsV);
            r = FindProcLowestSum (procstartV, procindexV, j, NprocsV[j], SumsU,SumsV);   
    
            /* If the matrix is partitioned one-dimensionally,
               then the vector is partitioned conformally.
               Otherwise, the processor with the lowest total sum is chosen.
               Ties are decided by flipping a coin, to avoid favouring 
               the input or output vector. */
            if (pOptions->SplitStrategy == OneDimRow)
               proc = q;
            else if (pOptions->SplitStrategy == OneDimCol)
               proc = r;
            else if (SumsU[q] + SumsV[q] < SumsU[r] + SumsV[r])
               proc = q;
            else if (SumsU[q] + SumsV[q] > SumsU[r] + SumsV[r])
               proc = r;
            else if (Random1(0,1) == 0)
               proc = q;
            else
               proc = r;
  
            /* Assign U[j] and V[j] to this processor: */
            if (!AssignColumnToProc(U, procstartU, procindexU, NsU, NrU, NvU, SumsU, j, proc)) {
                fprintf(stderr, "DistributeVecOrigEq(): Unable to assign column!\n");
                return -1;
            }
            if (!AssignColumnToProc(V, procstartV, procindexV, NsV, NrV, NvV, SumsV, j, proc)) {
                fprintf(stderr, "DistributeVecOrigEq(): Unable to assign column!\n");
                return -1;
            }
        }
    }
    
    
    /*## Step 4C:
      ##  Assign vector components j with NprocsU[j]=0 and NprocsV[j]=0: ##*/
    if (!AssignRemainingColumns(l, P, U, NvU)) {
        fprintf(stderr, "DistributeVecOrigEq(): Unable to assign remaining columns!\n");
        return -1;
    }
    
    for (j=0; j<l; j++) {
        if (V[j] == -1) {
            V[j] = U[j];
        }
        else if (V[j] != U[j]) {
            fprintf(stderr, "DistributeVecOrigEq(): Internal error step 4C!\n");
            return -1;
        }
    }
  
    if (!CalcCom(pM, V, ROW, &ComVol, &MaxOut, &MaxIn, &MaxCompnts, &TotCompnts)) {
        fprintf(stderr, "DistributeVecOrigEq(): Unable to calculate communication!\n");
        return -1;
    }
    maxcom =  MAX(MaxOut,MaxIn);
    
    /* Check whether total number of components equals length: */
    if (TotCompnts != l) {
        fprintf(stderr, "DistributeVecOrigEq(): Internal error: TotCompnts != length!\n");
        return -1;
    }
  
#ifdef INFO
    PrintCom(P, l, ROW, ComVol, MaxOut, MaxIn, MaxCompnts, TotCompnts);
#endif
    
    if (!CalcCom(pM, U, COL, &ComVol, &MaxOut, &MaxIn, &MaxCompnts, &TotCompnts)) {
        fprintf(stderr, "DistributeVecOrigEq(): Unable to calculate communication!\n");
        return -1;
    }
    maxcom +=  MAX(MaxOut,MaxIn);

    /* Check whether total number of components equals length: */
    if (TotCompnts != l) {
        fprintf(stderr, "DistributeVecOrigEq(): Internal error: TotCompnts != length!\n");
        return -1;
    }

#ifdef INFO
    /* Interchange roles of sends and receives for u  */
    PrintCom(P, l, COL, ComVol, MaxIn, MaxOut, MaxCompnts, TotCompnts);
#endif
  
    /*## Cleanup processor arrays: ##*/
    free(NprocsUV); /* no check needed */
    free(procindexV);
    free(procstartV);
    free(NprocsV);
    free(procindexU); 
    free(procstartU);
    free(NprocsU);
    free(SumsV);
    free(SumsU);
    free(NvV);
    free(NvU);
    free(NrV);
    free(NrU);
    free(NsV);
    free(NsU);
    free(ProcHistogramUV);
    free(ProcHistogramV);
    free(ProcHistogramU);
  
    return maxcom;
} /* end DistributeVecOrigEq */
  
  
int InitIntersection(long l, int P,
                      long *procstartU, int *procindexU,
                      long *procstartV, int *procindexV, int *NprocsUV) {

    /* This function computes the intersection of the processor lists for
       row and column j, for each j, 0 <= j < l, and brings the processors
       of the intersection to the front in the list of row j.

       Input:  procindexU, contains the processor numbers occurring in
                       the matrix rows. The processor numbers for one
                       row j are stored consecutively.
               procstartU, where procstartU[j] is the position in array
                       procindexU of the first processor number of row j,
                       0 <= j < l.
                       procstartU[l] = total number of processor numbers
                       of all the rows together.

               procstartV, procindexV, the same for columns.

       Output: NprocsUV, where NprocsUV[j] is the number of processors that own 
                       nonzeros both in matrix row j and column j.

               procindexU is modified such that
                       procindexU[procstartU[j]..procstartU[j]+NprocsUV[j]-1]
                       contains the processors in the intersection for row/col j.
                       The remaining processors of row j are in
                       procindexU[procstartU[j]+NprocsUV[j]..procstartU[j+1]-1] 

               procindexV is not modified.  */

    int q, t, *Tmp;
    long j, s, *Present;
    
    if (!procstartU || !procindexU || !procstartV || !procindexV || !NprocsUV) {
        fprintf(stderr, "InitIntersection(): Null arguments!\n");
        return FALSE;
    }

    Present = (long *) malloc(P * sizeof(long));
    Tmp = (int *) malloc(P * sizeof(int));
    if (Present == NULL || Tmp == NULL) {
        fprintf(stderr, "InitIntersection(): Not enough memory!\n");
        return FALSE;
    }

    /* Initialise boolean array */
    for (q=0; q<P; q++)
        /* processor q is not present in current row */
        Present[q] = -1;

    for (j=0; j<l; j++) {
        /* Read processors of current row j */
        for (s=procstartU[j]; s<procstartU[j+1]; s++) {
            q = procindexU[s];
            /* mark processor q present in row j, using the time stamp j */
            Present[q] = j; 
        }

        /* Read processors of column j and copy intersection into Tmp */
        t = 0;
        NprocsUV[j] = 0;
        for (s=procstartV[j]; s<procstartV[j+1]; s++) { 
            q = procindexV[s]; 
            if (Present[q] == j) {
                /* q is in intersection */
                Tmp[t] = q;
                t++;
                NprocsUV[j]++;
                Present[q] = -1; /* reset time stamp */
            }
        }

        /* Copy remaining processors of row j into Tmp */
        for (s=procstartU[j]; s<procstartU[j+1]; s++) { 
            q = procindexU[s]; 
            if (Present[q] == j) {  
                /* q is in row j but not in the intersection */ 
                Tmp[t] = q;
                t++;
            }
            /* no need to reset */
        }

        /* Copy processors from Tmp back into procindexU, in the new order */
        for (s=procstartU[j], t=0; s<procstartU[j+1]; s++, t++)
            procindexU[s] = Tmp[t]; 
    }

    free(Tmp);
    free(Present);
    
    return TRUE;
} /* end InitIntersection */


int FindProcLowestSum(long *procstart, int *procindex, long j, int npj,
                      long *SumsU, long *SumsV) {
    /* This function returns the processor with the lowest total sum
       in row/column j, searching the first npj processors 
       procindex[procstart[j]..procstart[j]+npj-1]. */

    int proc, q;
    long s, sum, min;
    
    if (!procstart || !procindex || !SumsU || !SumsV) {
        fprintf(stderr, "FindProcLowestSum(): Null arguments!\n");
        return INT_MIN;
    }
    
    if (npj < 1)
        return INT_MIN;

    proc = procindex[procstart[j]];
    if (npj == 1)
        return proc;
    min = SumsU[proc] + SumsV[proc];

    for (s = procstart[j]+1; s < procstart[j] + npj; s++) {
        q = procindex[s];
        sum = SumsU[q] + SumsV[q];
        if (sum < min) {
            proc = q;
            min = sum;
        }
    }
    
    return proc;
} /* end FindProcLowestSum */
