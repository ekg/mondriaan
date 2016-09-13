#include "DistributeVec.h"

long DistributeVec(const struct sparsematrix *pM, long int *X, int dir, const struct opts *pOptions) {
    /* This function distributes the vector components for a vector
       that can be freely distributed (hence with distr(u) <> distr(v)).

       If possible, the function uses the optimal graph algorithm Opt2
       for the special case Nprocs <=2. Otherwise it tries both the
       original Mondriaan method and the local-bound method,
       taking the best result. Each method is enhanced by
       greedy improvement. This is repeated several times,
       where the best result is kept.
   
       The local lower bound is used to stop the search for a better solution, 
       in case the bound is achieved, since this means that an optimal solution
       has been found.

       Input:  initialised A,
               uninitialised array X of length pM->n,
               dir is the direction:
               dir = ROW means that the vector is the input vector v,
                     COL means that the vector is the output vector u = Av.

       Output: vector distribution X, which is an array of length pM->n
               with values between 0 and P-1, where P = pM->NrProcs.
               The function returns maxcom = max(sends,recvs)
               over all processors. */

    int P, Pact, q,  useOpt2= FALSE, useMond, *Xbest, *Nprocs;
    long int *XC;
    long total, LB, actbound, bestbound, 
         ComVol, MaxOut, MaxIn, MaxCompnts, TotCompnts,
         maxcom, bestcom = LONG_MAX,
         l=0, i, j, k, s, t, tt, ncols, lC, nzC,
         *J, *Jinv, *JA, *Nv, *ProcHistogram,
         *Tmpi, *Tmpj, *Nz;

    /* Communication matrix, stored in compressed column storage (CCS).
       Includes all columns of the original matrix pM-> */
    int *procindex; /* array of length at most V+pM->n, where V is the
                       communication volume for the vector v.
                       procindex[procstart[j]..procstart[j+1]-1]
                       contains the processor numbers of the processors
                       that occur in column j */
    long *procstart; /* array of length pM->n+1 containing the starts of the
                       columns of the communication matrix */

    struct sparsematrix C; /* Smaller communication matrix in triple
                              format which only contains the columns from A
                              with >= 2 processors and with one nonzero
                              per processor in a column. This speeds up
                              the vector distribution. */
    if (!pM || !X || !pOptions) {
        fprintf(stderr, "DistributeVec(): Null parameters!\n");
        return -1;
    }
    
    if (dir == COL) {
        l = pM->m;
    } else if (dir == ROW) {
        l = pM->n;
    } else {
        fprintf(stderr, "DistributeVec(): Unknown direction!\n");
        return -1;
    }
    P = pM->NrProcs;

    /* Nprocs[j] is the number of processors that occur in column j */
    Nprocs = (int *)malloc(l*sizeof(int));
    ProcHistogram = (long *)malloc((P+1)*sizeof(long));
    JA = (long *)malloc(l*sizeof(long));
    Nv = (long *)malloc(P*sizeof(long));
    if (Nprocs == NULL || ProcHistogram == NULL ||
         JA == NULL || Nv == NULL) {
        fprintf(stderr, "DistributeVec(): Not enough memory!\n");
        return -1;
    }

    /* Initialise CCS data structure for communication matrix of A  */
    procstart = (long *)malloc((l+1)*sizeof(long));
    if (procstart == NULL) {
        fprintf(stderr, "DistributeVec(): Not enough memory!\n");
        return -1;
    }
    
    if (!InitNprocs(pM, dir, Nprocs)) {
        fprintf(stderr, "DistributeVec(): Unable to initialise processor array!\n");
        return -1;
    }

    total = 0;
    for (j=0; j<l; j++)
        total += Nprocs[j];
    procindex = (int *)malloc(total*sizeof(int));
    if (procindex == NULL) {
        fprintf(stderr, "DistributeVec(): Not enough memory!\n");
        return -1;
    }
    if (!InitProcindex(pM, dir, Nprocs, procstart, procindex)) {
        fprintf(stderr, "DistributeVec(): Unable to initialise processor index!\n");
        return -1;
    }

    if (!GenerateHistogram(Nprocs, l, 0, P, ProcHistogram)) {
        fprintf(stderr, "DistributeVec(): Unable to generate histogram!\n");
        return -1;
    }

#ifdef INFO
    if (dir == ROW)
        printf("\nNr matrix columns owned by p processors:\n");
    else
        printf("\nNr matrix rows owned by p processors:\n");
    PrintHistogram(0, P, ProcHistogram);
#endif
  
    /* Calculate the matrix-based communication volume V,
       the local lower bound, and a lower bound based
       on the number of active processors */

    if (!CalcCom(pM, NULL, dir, &ComVol, &MaxOut, &MaxIn, &MaxCompnts, &TotCompnts)) {
        fprintf(stderr, "DistributeVec(): Unable to calculate communication!\n");
        return -1;
    }
    if (!CalcLocalLowerBound(pM, dir, &LB, &Pact)) {
        fprintf(stderr, "DistributeVec(): Unable to calculate lower bound!\n");
        return -1;
    }
    if (Pact > 0)
        actbound= (ComVol%Pact == 0 ? ComVol/Pact : ComVol/Pact +1);
    else
        actbound= 0;
    bestbound= MAX(LB,actbound);
 
#ifdef INFO
    printf("Communication for vector"
           " %c (after matrix distribution):\n", dir == ROW ? 'v' : 'u');
    printf("  vol           : %ld \n", ComVol);
    printf("  active procs  : %d \n", Pact);
    printf("  active bound  : %ld \n", actbound);
    printf("  local bound   : %ld \n", LB);
    printf("  best bound    : %ld \n", bestbound);
#endif 

    /* All columns are unowned at the start */
    for (j=0; j<l; j++)
        X[j] = -1;

    /* Assign columns with 1 processor, remove empty columns,
       and determine columns of reduced P x lC matrix C */
    ncols = 0;
    nzC = 0;
    for (j=0; j<l; j++) {
        if (Nprocs[j] == 1) {
            /* Assign column j to unique owner */
            q = procindex[procstart[j]];
            X[j] = q;
        } else if (Nprocs[j] > 1) { 
            /* Column will be part of smaller communication matrix C */

            /* Store column number from original matrix */
            JA[ncols] = j;
            ncols++;

            /* Count the number of nonzeros in C */
            nzC += Nprocs[j];


        } /* else Nprocs[j] = 0, so column is empty
                  and will be assigned at the end */
    }
    lC = ncols;

    /* Create sparse matrix C */
    if (!MMSparseMatrixInit(&C)) {
        fprintf(stderr, "DistributeVec(): Unable to initialise sparse matrix!\n");
        return -1;
    }
    
    C.m = P;
    C.n = lC;
    C.NrNzElts = nzC; 
    C.NrProcs = P;
    C.MMTypeCode[0] = 'D'; /* distributed matrix */
    C.MMTypeCode[2] = 'P'; /* pattern only */

    /* Create smaller matrix C and distribute the corresponding vector */
    if (nzC > 0) {
        /* Matrix C is not empty */
        if (!MMSparseMatrixAllocateMemory(&C)) {
            fprintf(stderr, "DistributeVec(): Unable to allocate sparse matrix!\n");
            return -1;
        }
        
        Tmpi = (long *) malloc(nzC* sizeof(long));
        Tmpj = (long *) malloc(nzC* sizeof(long));
        Nz = (long *) malloc(P* sizeof(long));
        if (Tmpi == NULL || Tmpj == NULL || Nz == NULL) {
            fprintf(stderr, "DistributeVec(): Not enough memory!\n");
            return -1;
        }

        /* Create nonzeros of C and insert them into a temporary array */
        ncols = 0;
        t = 0;
        for (q=0; q<P; q++)
            Nz[q] = 0; /* Number of nonzeros of processor q */

        for (j=0; j<l; j++) {
            if (Nprocs[j] > 1) {
                for (s=procstart[j]; s<procstart[j+1]; s++) {
                    q = procindex[s];
                    Nz[q]++;
                    Tmpi[t] = q;
                    Tmpj[t] = ncols;
                    t++;
                }
                ncols++;
            }
        }
        
        /* Initialise C.Pstart */
        C.Pstart[0] = 0;
        for (q=0; q<P; q++)
            C.Pstart[q+1] = C.Pstart[q] + Nz[q];
            
        /* Copy nonzeros into C in order of increasing processor number */
        for (q=0; q<P; q++)
            Nz[q] = 0;
        for (t=0; t<nzC; t++) {
            q = Tmpi[t]; /* processor number */
            tt = C.Pstart[q] + Nz[q];
            C.i[tt] = q;
            C.j[tt] = Tmpj[t];
            Nz[q]++;
        }
        free(Nz);
        free(Tmpj);
        free(Tmpi);

#ifdef INFO2
        printf("  Reduced communication matrix has size %d x %ld \n", P, lC);
#endif

        /* Allocate and initialise distribution vector for 
           communication matrix C */
        XC= (long int *)malloc(lC*sizeof(long int));
        if (XC == NULL) {
            fprintf(stderr, "DistributeVec(): Not enough memory!\n");
            return -1;
        }

        for (j=0; j<lC; j++)
            XC[j] = -1;
 
        /**** Distribute the vector XC ****/

        /* Determine if algorithm Opt2 can be used, and if so, run it.
           ROW direction because C is a P x lC matrix */
        useOpt2 = DistributeVecOpt2(&C, XC, ROW);

        /* Communication cost of best solution found so far is infinity */

        if (useOpt2) {
#ifdef INFO2      
            printf("  **** Opt2 algorithm used: Nprocs[j] <= 2, all j.\n");
#endif
        } else {
            Xbest = (int *)malloc(lC*sizeof(int));
            J = (long *)malloc(lC*sizeof(long));
            Jinv = (long *)malloc(lC*sizeof(long));
            if (Xbest == NULL || J == NULL || Jinv== NULL) {
                fprintf(stderr, "DistributeVec(): Not enough memory!\n");
                return -1;
            }

            /* Initialise J to the identity permutation */
            for (j=0; j < lC; j++) 
                J[j] = j;

            bestcom = LONG_MAX;


            /* useMond registers whether best solution is obtained by
               original Mondriaan method */
            useMond = FALSE; 

            for (k=0; k<pOptions->VectorPartition_MaxNrLoops &&
                       bestbound < bestcom; k++) {

                /* For k=0, use the natural ordering of the matrix columns,
                   since this might have an advantage.
                   For k>0, randomly permute the columns. */

                if (k > 0) {
                    /* Determine a random permutation J */
                    RandomPermute(J, NULL, NULL, NULL, 0, lC-1);
 
                    /* Permute the columns of the matrix C by J */
                    for (t=0; t < C.NrNzElts; t++)
                        C.j[t] = J[C.j[t]];
                }

                /** Use the original Mondriaan algorithm **/
#ifdef INFO2      
                printf("  **** Try original Mondriaan algorithm");
                printf(" for vector distribution\n");
#endif
                maxcom = DistributeVecOrig(&C, XC, ROW, pOptions);
                
                if (maxcom < 0) {
                    fprintf(stderr, "DistributeVec(): Unable to use DistributeVecOrig!\n");
                    return -1;
                }
   
                if (maxcom<bestcom) {
                    /* Save best vector distribution in original order */
                    for (j=0; j<lC; j++)
                        Xbest[j] = XC[J[j]];
                    bestcom = maxcom;
                    useMond = TRUE;
                }

                /* Improve greedily */
                for (i=0; i<pOptions->VectorPartition_MaxNrGreedyImproves &&
                          bestbound < bestcom; i++) {
                    maxcom = DistributeVecGreedyImprove(&C, XC, ROW, pOptions);
                    
                    if (maxcom < 0) {
                        fprintf(stderr, "DistributeVec(): Unable to use DistributeVecGreedyImprove!\n");
                        return -1;
                    }
                    
                    if (maxcom<bestcom) {
                        for (j=0; j<lC; j++)
                            Xbest[j] = XC[J[j]];
                        bestcom = maxcom;
                        useMond = TRUE; 
                    }
                }
                if (bestcom == bestbound)
                    break; /* from k-loop. 
                              No need to try the local lower bound algorithm, or
                              to permute C back to the original column order */
        
                /** Use the local lower bound algorithm **/
#ifdef INFO2      
                printf("  **** Try local lower bound algorithm");
                printf("       for vector distribution\n");
#endif
        
                maxcom = DistributeVecLocal(&C, XC, ROW);
                
                if (maxcom < 0) {
                    fprintf(stderr, "DistributeVec(): Unable to use DistributeVecLocal!\n");
                    return -1;
                }
                
                if (maxcom<bestcom) {
                    for (j=0; j<lC; j++)
                        Xbest[j] = XC[J[j]];
                    bestcom = maxcom;
                    useMond = FALSE; 
                }
            
                /* Improve greedily */
                for (i=0; i<pOptions->VectorPartition_MaxNrGreedyImproves &&
                          bestbound < bestcom; i++) {
                    maxcom = DistributeVecGreedyImprove(&C, XC, ROW, pOptions);
                    
                    if (maxcom < 0) {
                        fprintf(stderr, "DistributeVec(): Unable to use DistributeVecGreedyImprove!\n");
                        return -1;
                    }
                    
                    if (maxcom<bestcom) {
                        for (j=0; j<lC; j++)
                           Xbest[j] = XC[J[j]];
                        bestcom = maxcom;
                        useMond = FALSE; 
                    }
                }   
                if (k>0) {
                    /* Determine the inverse of J */
                    for (j=0; j < lC; j++)
                        Jinv[J[j]] = j;

                    /* Permute the columns of the matrix C back by Jinv */
                    for (t=0; t < C.NrNzElts; t++)
                        C.j[t] = Jinv[C.j[t]];
                }

            } /* end k-loop */

            /* Copy best solution into XC; it is already in the right order */
            for (j=0; j<lC; j++)
                 XC[j] = Xbest[j];
            free(Jinv);
            free(J);
            free(Xbest);

#ifdef INFO2      
            if (useMond)
                printf("  **** Original Mondriaan"
                       " vector distribution chosen\n");
            else
                printf("  **** Local lower bound"
                       " vector distribution chosen\n"); 
#endif

        }

        /* Transfer solution for C into solution for A */
        for (j=0; j<lC; j++) 
            X[JA[j]] =  XC[j];

        free(XC);
        MMDeleteSparseMatrix(&C);

        /* end if nzC > 0 */
    } else
        bestcom = 0;

    if (!AssignRemainingColumns(l, P, X, Nv)) {
        fprintf(stderr, "DistributeVec(): Unable to assign remaining columns!\n");
        return -1;
    }

    /* Compute and print results */
    if (!CalcCom(pM, X, dir, &ComVol, &MaxOut, &MaxIn, &MaxCompnts, &TotCompnts)) {
        fprintf(stderr, "DistributeVec(): Unable to calculate communication!\n");
        return -1;
    }
    
    maxcom = MAX(MaxOut,MaxIn);
    
    if (useOpt2 == FALSE && maxcom != bestcom) {
        fprintf(stderr, "DistributeVec(): Error in maxcom computation!\n");
        return -1;
    }

#ifdef INFO
    if (maxcom == bestbound)
        printf("  Vector distribution is optimal\n");
    PrintCom(P, l, dir, ComVol, MaxOut, MaxIn, MaxCompnts, TotCompnts);
#endif

    free(procindex);
    free(procstart);
    free(Nv);
    free(JA);
    free(ProcHistogram);
    free(Nprocs);

    return maxcom;
} /* end DistributeVec */
