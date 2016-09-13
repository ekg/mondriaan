#include "DistributeVecLocal.h"

/* All comments only consider the distribution of the vector v,
   where all vector-length arrays are of length pM->n.
   This corresponds to the direction dir = ROW,
   meaning that the vector is in the row direction. 
   In case of the vector u =  Av, i.e., dir = COL,
   one should read pM->m instead. */

long DistributeVecLocal(const struct sparsematrix *pM, long int *X, int dir) {
    /*  This function determines a vector distribution X by first running 
        the local-bound algorithm. Remaining components are assigned greedily.
        
        Input:  initialised A, 
                uninitialised array X of length pM->n,
                dir is the direction: 
                dir = ROW means that the vector is the input vector v,
                      COL means that the vector is the output vector u = Av. 

        Output: vector distribution X, which is an array of length pM->n
                with values between 0 and P-1, where P= pM->NrProcs.
                The function returns maxcom = max(sends,recvs)
                over all processors. */



    long total, tmp, maxcom, maxval,
         s, t, j, l=0, Nstotal, Nrtotal, localmax;
    int P, qmax, q, q1, r, r1;

    /* Arrays of length pM->n */
    int *Nprocs;    /* Nprocs[j] is the number of processors
                       that occur in column j */
    long *Nprocs2;  /* a copy in long */

    /* Communication matrix, stored in compressed column storage (CCS).
       Includes all columns of the original matrix. */
    int *procindex; /* array of length at most V+pM->n, where V is the
                       communication volume for the vector v.
                       procindex[procstart[j]..procstart[j+1]-1]
                       contains the processor numbers of the processors
                       that occur in column j */
    long *procstart; /* array of length pM->n+1 containing the starts of the
                       columns of the communication matrix */
 
    /* Communication matrix, stored in compressed row storage (CRS).
       Includes only columns with Nprocs[j] > 1. 
       Columns are removed from the matrix after becoming owned. */
    long *J;        /* array of length at most V+pM->n. 
                       J[Jstart[q]..Jstart[q+1]-1] contains the still unowned
                       columns j with Nprocs[j]>1 in which processor q occurs */
    long *Jstart;   /* array of length P+1 containing the current starts of the
                       rows of the communication matrix */
    long *Jstart0;  /* array of length P+1. A copy of the initial Jstart,
                       needed to find the end of the rows. */
    
    /* Arrays of length P */
    long *Ns;       /* Ns[q] = the current number of vector components
                       sent by processor q to other processors */
    long *Nr;       /* Nr[q] = the current number of vector components
                       received by processor q from other processors */
    long *Nv;       /* Nv[q] = the current number of vector components
                       owned by processor q */
    long *Ls;       /* Ls[q] = the number of vector components
                       that processor q should send to other processors
                       in order to obtain the generalised local lower bound */
    long *Lr;       /* Lr[q] = the number of vector components
                       that processor q should receive from other processors
                       in order to obtain the generalised local lower bound */

    /* Matrices of size P x P+1 */
    long **Needed;    /* Needed[q,r]= the number of columns with 
                         r processors that processor q needs to obtain 
                         the generalised local lower bound */
    long **Available; /* Available[q,r]= the number of columns with 
                         r processors and with processor q occurring 
                         in them, that are owned by q or are still unowned */
  
    /* Array of length P */
    int *Prev;     /* Prev[q] is the most recently checked column size r
                      for which Needed[q,r] < Available[q,r].
                      The use of Prev means that we only have to go once
                      through the rows of Needed and Available */
    if (!pM || !X) {
        fprintf(stderr, "DistributeVecLocal(): Null arguments!\n");
        return -1;
    }

    if (dir == COL) {
        l = pM->m;
    } else if (dir == ROW) {
        l = pM->n;
    } else {
        fprintf(stderr, "DistributeVecLocal(): Unknown direction!\n");
        return -1;
    }
    P = pM->NrProcs; 

    /* ## Initialise CCS data structure for communication matrix ## */
    Nprocs = (int *)malloc(l*sizeof(int));
    Nprocs2 = (long *)malloc(l*sizeof(long));
    procstart = (long *)malloc((l+1)*sizeof(long));
    if (Nprocs == NULL || Nprocs2 == NULL || procstart == NULL) {
        fprintf(stderr, "DistributeVecLocal(): Not enough memory!\n");
        return -1;
    }

    if (!InitNprocs(pM, dir, Nprocs)) {
        fprintf(stderr, "DistributeVecLocal(): Unable to initialise processor array!\n");
        return -1;
    }
    
    /* make a copy with long integers instead of integers (to enable use of CSort) */
    for (j=0; j<l; j++)
        Nprocs2[j] = Nprocs[j];

    total = 0;
    for (j=0; j<l; j++)
        total += Nprocs[j];
    procindex = (int *)malloc(total*sizeof(int));
    if (procindex == NULL) {
        fprintf(stderr, "DistributeVecLocal(): Not enough memory!\n");
        return -1;
    }

    if (!InitProcindex(pM, dir, Nprocs, procstart, procindex)) {
        fprintf(stderr, "DistributeVecLocal(): Unable to initialise processor index!\n");
        return -1;
    }

    /* ## Initialise processor arrays ## */
    Ns = (long *)malloc(P*sizeof(long));
    Nr = (long *)malloc(P*sizeof(long));
    Nv = (long *)malloc(P*sizeof(long));
    J = (long *)malloc(total*sizeof(long));
    Jstart = (long *)malloc((P+1)*sizeof(long));
    Jstart0 = (long *)malloc((P+1)*sizeof(long));
    Ls = (long *)malloc(P*sizeof(long));
    Lr = (long *)malloc(P*sizeof(long));
    Available = matallocl(P,P+1);
    Needed = matallocl(P,P+1);
    Prev = (int *)malloc(P*sizeof(int));

    if (Ns == NULL || Nr == NULL || Nv == NULL || J == NULL ||
        Jstart == NULL || Jstart0 == NULL || Ls == NULL || Lr == NULL ||
        Available == NULL || Needed == NULL || Prev == NULL) {
        fprintf(stderr, "DistributeVecLocal(): Not enough memory!\n");
        return -1;
    }
 
    /* Initialise arrays */
    for (q=0; q<P; q++) {
        Ns[q] = 0;
        Nr[q] = 0;
        Nv[q] = 0;
        Jstart[q] = 0;
        Ls[q] = 0;
        Lr[q] = 0;
        Prev[q] = 0;
    }
    Jstart[P] = 0;
    
    /* All columns are unowned at the start */
    for (j=0; j<l; j++)
        X[j] = -1;

    for (q=0; q<P; q++) {
        for (r=0; r<=P; r++) {
            Available[q][r] = 0;
            Needed[q][r] = 0;  
        }
    }

    /* Initialise Jstart, Lr, Available */
    for (j=0; j<l; j++) {
        if (Nprocs[j] > 1) {
            for (s=procstart[j]; s<procstart[j+1]; s++) {
                q = procindex[s];
                
                /* Jstart[q] is initialised to the number of columns 
                   with Nprocs > 1 in which q occurs */
                Jstart[q]++;

                /* Lr is initialised assuming that each unowned column
                   with Nprocs > 1 in which q occurs causes a receive for q */
                Lr[q]++;
                
                /* Available is initialised assuming that every column is available
                   at the start */
                Available[q][Nprocs[j]]++;
            }
        }
    }
    
    /* Determine Needed, decrease Lr and increase Ls
       while maintaining Ls <= Lr */
    for (q=0; q<P; q++) {
        if (Ls[q] < Lr[q]) {
            for (r=2; r<=P; r++) {           
                if (Ls[q] + (r-1)*Available[q][r] <= Lr[q] - Available[q][r]) {
                    /* All available columns with r processors are needed to achieve
                       the local lower bound of processor q */
                    Needed[q][r] = Available[q][r];
                    Ls[q] += (r-1)*Needed[q][r];
                    Lr[q] -= Needed[q][r];
                } else {
                    /* Only part of the columns is needed.
                       More columns will actually hurt. */
                    Needed[q][r] = (Lr[q]-Ls[q])/r;
                    Ls[q] += (r-1)*Needed[q][r];
                    Lr[q] -= Needed[q][r];
                    break;
                }
            }
            /* Ls[q] <= Lr[q] */
        }
    }
    
    /* ## Initialise CRS data structure for communication matrix ## */

    /* Make Jstart cumulative and make a copy */
    total = 0;
    for (q=0; q<P; q++) {
        tmp = total;
        total += Jstart[q];
        Jstart[q] = tmp;
        Jstart0[q] = Jstart[q];
    }
    Jstart0[P] = Jstart[P] = total;

    /* Initialise J, modifying Jstart */
    for (j=0; j<l; j++) {
        if (Nprocs[j] > 1) {
            for (s=procstart[j]; s<procstart[j+1]; s++) {
                q = procindex[s];
                J[Jstart[q]] = j;
                Jstart[q]++;
            }
        }
    }

    /* Repair Jstart */
    for (q=0; q<=P; q++)
        Jstart[q] = Jstart0[q];

    /* Sort each processor's part of the J array by increasing Nprocs
       using a counting sort */
    maxval = P;
    for (q=0; q<P; q++)
        CSort(J, Nprocs2, maxval, Jstart[q], Jstart[q+1]-1);

#ifdef INFO2
    if (dir == ROW) {
        printf("\n Ns, Nr = number of sends/recvs for owned columns\n");
        printf(" Ls, Lr = number of sends/recvs with lower bound\n");
        printf(" Un = number of unowned columns with Nprocs>=2\n");
    } else {
        printf("\n Ns, Nr = number of recvs/sends for owned rows\n");
        printf(" Ls, Lr = number of recvs/sends with lower bound\n");
        printf(" Un = number of unowned rows with Nprocs>=2\n");
    }
    printf("Initial local bounds:\n");
    PrintVecLocalStatistics(P, Ns, Nr, Ls, Lr, Jstart0, J, X);
#endif

    /* First phase: assign as many columns as possible with Nprocs > 1
       in accordance with the sorted J array */

    q = -1;
    while(TRUE) {   
        /* Find a new processsor q, if we are in the first iteration or 
           if q has already attained its local lower bound. */
        if (q == -1 || Ns[q] == Ls[q]) {
            localmax = -1;

            for (r=0; r<P; r++) {
                if (Ns[r] < Ls[r]) {
                    /* Processor r is only taken as a candidate if it
                       has not yet attained its local lower bound */
                    if (Lr[r] > localmax) {
                        q = r;
                        localmax = Lr[r]; /* Lr[r] >= Ls[r],
                                             hence Lr[r] is the bound */
                    }
                }
            }

            if (localmax == -1)
                /* All processors have reached their local lower bound */
                break;  /* from the while loop */
        }
        
        /* Find the first unowned column j in which q occurs.
           The begin point of the search moves to the right, 
           but the end point remains fixed. */
        j = -1;
        for (t=Jstart[q]; t< Jstart0[q+1]; t++) {
            j = J[t];
            Jstart[q]++; /* next time, start one place beyond t
                            since j = J[t] will be owned then */
            if (X[j] == -1)
                break; /* unowned column found */
        }
        
        if (j == -1 || X[j] != -1) {
            /* Did not find an unowned column. This should not happen! */
            fprintf(stderr, "DistributeVecLocal(): Did not find unowned column!\n");
            return -1;
        }

        /* Assign column j to processor q */
        if (!AssignColumnToProc(X, procstart, procindex, Ns, Nr, Nv, NULL, j, q)) {
            fprintf(stderr, "DistributeVecLocal(): Unable to assign column!\n");
            return -1;
        }

        /* Update the Ls, Lr, Available and Needed values of the processors
           receiving in column j and find the new processor qmax with the
           maximum lower bound */
        qmax = q;
        r = Nprocs[j];
        for (s=procstart[j]; s<procstart[j+1]; s++) {
            q1 = procindex[s];

            if (q1 != q) {
                Available[q1][r]--;

                /* If processor q1 needs more columns with r processors than 
                   available, then adjust its local values */
                if (Needed[q1][r] > Available[q1][r]) {
                    Needed[q1][r]--; 
                    /* restored Needed[q1][r] <= Available[q1][r] */
                    Ls[q1] -= r - 1;
                    Lr[q1]++;

                    /* Search for a column size r1 for which columns are still 
                       available and not needed, and use these to increase Ls
                       and decrease Lr, if possible */
                    for (r1 = Prev[q1]; r1 <= P && 
                         Needed[q1][r1] == Available[q1][r1]; r1++) ;              /* equality holds for r1 = 0,1 */
                    
                    if (r1 <= P) {
                        /* Needed[q1][r1] < Available[q1][r1] */
                        if (Ls[q1] + r1-1 <= Lr[q1]-1) {
                            Ls[q1] += r1-1;
                            Lr[q1]--;
                            Needed[q1][r1]++;
                        }
                        /* Ls[q1] <= Lr[q1] */
                    }
                    Prev[q1] = r1;

                    /* Processor q1 may have the largest local bound */
                    if (Lr[q1] > Lr[qmax])
                        qmax = q1;
                }
            }
        }
        q = qmax; /* new processor with largest lower bound */
    } /* end while loop */

#ifdef INFO2
    printf("Local bounds after first phase of local bound algorithm:\n");
    PrintVecLocalStatistics(P, Ns, Nr, Ls, Lr, Jstart0, J, X);
#endif

    /* Second phase: greedily assign all columns left over from phase 1 */
 
    /* Assign all remaining columns with Nprocs > 0 */
    if (!AssignRemainingNonemptyColumns(l,P,X,procstart,procindex,Ns,Nr,Nv)) {
        fprintf(stderr, "DistributeVecLocal(): Unable to assign remaining nonempty columns!\n");
        return -1;
    }

    /* Assign all columns with Nprocs = 0 */
    if (!AssignRemainingColumns(l,P,X,Nv)) {
        fprintf(stderr, "DistributeVecLocal(): Unable to assign remaining columns!\n");
        return -1;
    }

    /* Check if any columns are still unowned */
    for (j=0; j<l; j++) {
        if (X[j]  < 0 || X[j] >= P) {
            fprintf(stderr, "DistributeVecLocal(): Columns left unowned!\n");
            return -1;
        }
    }

#ifdef INFO2
    for (q=0; q<P; q++) { 
        Ls[q] = Ns[q];
        Lr[q] = Nr[q];
    }
    
    printf("Final local bounds:\n");
    PrintVecLocalStatistics(P, Ns, Nr, Ls, Lr, Jstart0, J, X);
#endif

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
        fprintf(stderr, "DistributeVecLocal(): Total sends != total recvs!\n");
        return -1;
    }

    /* wrapping up */
    free(Prev);
    matfreel(Needed);
    matfreel(Available);
    free(Lr);
    free(Ls);
    free(Jstart0);
    free(Jstart);
    free(J);
    free(Nv);
    free(Nr);
    free(Ns);
    free(procindex);
    free(procstart);
    free(Nprocs2);
    free(Nprocs);

    return maxcom;
} /* end DistributeVecLocal */


void PrintVecLocalStatistics(int P, long *Ns, long *Nr, long *Ls, long *Lr,
                             long *Jstart0, long *J, long *owner) {
                          
    /* This function prints for each processor the send, receive,
       and local-bound values and the number of unowned columns:
       Ns, Nr = number of sends/recvs for owned columns;
       Ls, Lr = number of sends/recvs with lower bound;
       Un = number of unowned columns with Nprocs>=2. */


    long t, unowned;
    int q;
    
    if (!Ns || !Nr || !Ls || !Lr || !Jstart0 || !J || !owner) return;

    for (q=0; q<P; q++) {
        /* Count the number of unowned columns in which processor q occurs */
        unowned= 0;
        for (t=Jstart0[q];  t<Jstart0[q+1]; t++) {
            if (owner[J[t]] == -1)
                unowned++;
        }

        printf("  Proc[%d]: Ns = %ld, Nr = %ld, Ls = %ld, Lr = %ld, "
               "Un = %ld\n",q,Ns[q],Nr[q],Ls[q],Lr[q], unowned);
    }
    
}
/* end PrintVecLocalStatistics */
