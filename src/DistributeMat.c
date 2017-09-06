#include "DistributeMat.h"
#include "Permute.h"
#include "SplitMatrixUpperBound.h"

#ifdef USE_PATOH
#include <patoh.h>
#endif

int BalanceParts(long *weight, long maxweight, int P, int k, int *procs);
int DetermineSplit(long *weight, long maxweight, int k, int *procs,
                     int *isplit, long *weightlo, long *weighthi, const struct opts *pOptions);

int logb2(int n);

int SplitMatrixKLFM(struct sparsematrix *pT, int k, int i, int dir, 
                     long weightlo, long weighthi, const struct opts *pOptions);
int SplitMatrixSimple(struct sparsematrix *pT, int k, int i,
                       long weightlo, long weighthi, const struct opts *pOptions);

int SplitMatrixZeroVolume(struct sparsematrix *pT, int k, int i,
                     long weightlo, long weighthi, const struct opts *pOptions);

#ifdef USE_PATOH
struct patohnz {
    int P;
    long i, j;
    double Re, Im;
};

int PaToHNzCompare(const void *_a, const void *_b) {
    const struct patohnz *a = (const struct patohnz *)_a, *b = (const struct patohnz *)_b;
    
    return a->P - b->P;
}

int DistributeMatrixPaToH(struct sparsematrix *pT, int P, double eps, const struct opts *pOptions) {
    /* This function partitions the nonzeros of the matrix T
       into P parts, allowing a load imbalance of a fraction eps,
       using the PaToH hypergraph partitioner.
    */
    
    int _c, _n, *cwghts, *nwghts, *xpins, *pins, *partvec, cut, *partweights;
    float *targetweights;
    const int _nconst = 1, _k = P; /* Number of constraints and parts to which we will split. */
    PaToH_Parameters args;
    
    struct biparthypergraph HG, *pHG;
    struct patohnz *Ordering;
    long t, tt, n;
    int Result;
    
    if (pOptions->OrderPermutation != OrderNone) {
        fprintf(stderr, "DistributeMatrixPaToH(): PaToH partitioning does not support permutations at the moment!\n");
        return FALSE;
    }
    
    switch (pOptions->SplitStrategy) {
    case OneDimRow:
        Result = SparseMatrix2BiPartHyperGraph(pT, ROW, pOptions, &HG);
        break;
    case OneDimCol:
        Result = SparseMatrix2BiPartHyperGraph(pT, COL, pOptions, &HG);
        break;
    case SFineGrain:
        /* FIXME */
        fprintf(stderr, "*WARNING* DistributeMatrixPaToH(): Symmetric finegrain in PaToH is untested!\n"); 
        Result = SparseMatrix2BiPartHyperGraph(pT, SFINEGRAIN, pOptions, &HG);
        break;
    case FineGrain:
        Result = SparseMatrix2BiPartHyperGraph(pT, FINEGRAIN, pOptions, &HG);
        break;
    default:
        fprintf(stderr, "DistributeMatrixPaToH(): Invalid splitting strategy: valid choices are onedimrow, onedimcol, and finegrain!\n");
        return FALSE;
    }
    
    if (!Result) {
        fprintf(stderr, "DistributeMatrixPaToH(): Unable to convert matrix to hypergraph!\n");
        return FALSE;
    }
    
    pHG = &HG;
    
    cwghts = (int *)malloc(pHG->NrVertices*sizeof(int));
    nwghts = (int *)malloc(pHG->NrNets*sizeof(int));
    xpins = (int *)malloc((pHG->NrNets + 1)*sizeof(int));
    pins = (int *)malloc(pHG->NrPins*sizeof(int));
    partvec = (int *)malloc(pHG->NrVertices*sizeof(int));
    partweights = (int *)malloc(_k*sizeof(int));
    targetweights = (float *)malloc(_k*sizeof(float));
    
    if (cwghts == NULL || nwghts == NULL || xpins == NULL || pins == NULL || partvec == NULL || partweights == NULL || targetweights == NULL) {
        fprintf(stderr, "DistributeMatrixPaToH(): Not enough memory!\n");
        return FALSE;
    }
    
    /* Set default parameters. */
    PaToH_Initialize_Parameters(&args, pOptions->Metric == MetricCut ? PATOH_CUTPART : PATOH_CONPART, PATOH_SUGPARAM_DEFAULT);
    
    /* Copy data from hypergraph. */
    _c = pHG->NrVertices;
    _n = pHG->NrNets;
    
    for (t = 0; t < pHG->NrVertices; ++t) {
        partvec[t] = 0;
        cwghts[t] = pHG->V[t].vtxwgt;
    }
    
    for (t = 0; t < pHG->NrNets; ++t) {
        xpins[t] = pHG->N[t].iStartP0;
        nwghts[t] = pHG->N[t].netwgt;
    }
    
    xpins[pHG->NrNets] = pHG->NrPins;
    
    for (t = 0; t < pHG->NrPins; ++t) {
        pins[t] = pHG->NetAdjncy[t];
    }
    
    args._k = _k;
    args.seed = pOptions->Seed;
    args.final_imbal = args.init_imbal = eps;
    
    for (t = 0; t < _k; ++t) {
        targetweights[t] = 1.0f;
        partweights[t] = 0;
    }
    
    /* Allocate data and perform partitioning. */
    PaToH_Alloc(&args, _c, _n, _nconst, cwghts, nwghts, xpins, pins);
    PaToH_Part(&args, _c, _n, _nconst, 0, cwghts, nwghts, xpins, pins, targetweights, partvec, partweights, &cut);
    
    /* Free unused data. */
    PaToH_Free();
    
    free(cwghts);
    free(nwghts);
    free(xpins);
    free(pins);
    free(partweights);
    free(targetweights);
    
    /* Order nonzeros by processor index. */
    Ordering = (struct patohnz *)malloc(pT->NrNzElts*sizeof(struct patohnz));
    
    if (Ordering == NULL) {
        fprintf(stderr, "DistributeMatrixPaToH(): Not enough memory!\n");
        return FALSE;
    }
    
    switch (pOptions->SplitStrategy) {
    case OneDimRow:
        n = 0;
        
        for (t = 0; t < pHG->NrVertices; ++t) {
            for (tt = pHG->V[t].iStart; tt < pHG->V[t].iEnd; ++tt) {
                Ordering[n].i = pHG->Vtx2MatIndex[t];
                Ordering[n].j = pHG->Net2MatIndex[pHG->VtxAdjncy[tt]];
                Ordering[n].P = partvec[t];
                if (pHG->MatReValue != NULL) Ordering[n].Re = pHG->MatReValue[tt];
                if (pHG->MatImValue != NULL) Ordering[n].Im = pHG->MatImValue[tt];
                ++n;
            }
        }
        
        break;
    case OneDimCol:
        n = 0;
        
        for (t = 0; t < pHG->NrVertices; ++t) {
            for (tt = pHG->V[t].iStart; tt < pHG->V[t].iEnd; ++tt) {
                Ordering[n].i = pHG->Net2MatIndex[pHG->VtxAdjncy[tt]];
                Ordering[n].j = pHG->Vtx2MatIndex[t];
                Ordering[n].P = partvec[t];
                if (pHG->MatReValue != NULL) Ordering[n].Re = pHG->MatReValue[tt];
                if (pHG->MatImValue != NULL) Ordering[n].Im = pHG->MatImValue[tt];
                ++n;
            }
        }
        
        break;
    case SFineGrain:
        fprintf(stderr, "*WARNING* DistributeMatrixPaToH(): Symmetric finegrain in PaToH translates to normal (unsymmetric) finegrain!\n");
        /* FIXME */
    case FineGrain:
        for (t = 0; t < pHG->NrVertices; ++t) {
            for (tt = pHG->V[t].iStart; tt < pHG->V[t].iEnd; ++tt) {
                n = pHG->VtxAdjncy[tt];
                
                if (pHG->N[n].dir == ROW) Ordering[t].i = pHG->Net2MatIndex[n];
                else Ordering[t].j = pHG->Net2MatIndex[n];
            }
            
            Ordering[t].P = partvec[t];
            if (pHG->MatReValue != NULL) Ordering[t].Re = pHG->MatReValue[t];
            if (pHG->MatImValue != NULL) Ordering[t].Im = pHG->MatImValue[t];
        }
        break;
    default:
        break;
    }
    
    qsort(Ordering, pT->NrNzElts, sizeof(struct patohnz), PaToHNzCompare);
    
    /* Free data. */
    free(partvec);
    DeleteBiPartHyperGraph(&HG);
    
    /* Find processor indices and reorder nonzeros. */
    if (pT->Pstart == NULL) {
        fprintf(stderr, "Warning (DistributeMatrixPaToH): Pstart array was not initialised. Doing this now.\n");
        if (!PstartInit(pT, P)) {
             fprintf(stderr, "DistributeMatrixPaToH(): error during initialisation of Pstart!\n");
             return FALSE;
        }
    }
    
    pT->Pstart[0] = pT->Pstart[1] = 0;
    tt = 0;
    
    for (t = 0; t < pT->NrNzElts; ++t) {
        /* printf("P[%ld] = %d\n", t, Ordering[t].P); */
        
        pT->i[t] = Ordering[t].i;
        pT->j[t] = Ordering[t].j;
        if (pT->ReValue != NULL) pT->ReValue[t] = Ordering[t].Re;
        if (pT->ImValue != NULL) pT->ImValue[t] = Ordering[t].Im;
        
        if (Ordering[t].P < 0 || Ordering[t].P >= P) {
            fprintf(stderr, "DistributeMatrixPaToH(): Wrongly assigned nonzero (0 <= %d < %d)!\n", Ordering[t].P, P);
            return FALSE;
        }
        
        if (Ordering[t].P == tt) {
            pT->Pstart[tt + 1]++;
        }
        else {
            /* printf("PSTART[%ld] = %ld, PSTART[%ld] = %ld.\n", tt, pT->Pstart[tt], tt + 1, pT->Pstart[tt + 1]); */
            ++tt;
            pT->Pstart[tt + 1] = pT->Pstart[tt] + 1;
        }
    }
    
    if (pT->Pstart[P] != pT->NrNzElts) fprintf(stderr, "DistributeMatrixPaToH(): Sanity check failed!\n");
    
    /* Free data. */
    free(Ordering);
    
    return TRUE;
}
#else
int DistributeMatrixPaToH(struct sparsematrix *pT, int P, double eps, const struct opts *pOptions) {
    fprintf(stderr, "DistributeMatrixPaToH(): This version of Mondriaan was compiled without PaToH support! Please use the Mondriaan partitioner.\n");
    return FALSE;
}
#endif

void VerifyLambdas(const int *lambdas, const long n, const int P) {
    /* Verification function of the calculated lambdas. */
    long *hist = (long *)malloc(P*sizeof(long));
    long t, c;
    
    if (hist == NULL) return;
    
    for (t = 0; t < P; t++) hist[t] = 0;
    
    c = 0;
    
    for (t = 0; t < n; t++) {
        long l = lambdas[t];
        
         c += l;
         hist[l]++;
    }
    
    for (t = 0; t < P; t++) printf("% 4ld %ld\n", t, hist[t]);
    printf("Sum: %ld (%ld)\n", c, c - n);
    
    free(hist);
}

int DistributeMatrixMondriaan(struct sparsematrix *pT, int P, double eps, const struct opts *pOptions, int (*Callback)(int, int, const struct sparsematrix *)) {
    /* This function partitions the nonzeros of the matrix T
       into P parts, allowing a load imbalance of a fraction eps.

       Input: T sparse matrix,
              P number of parts, P >= 1,
              eps load imbalance fraction, eps >= 0.
       Output: T distributed sparse matrix.
               The nonzeros of part i are in positions
                   pT->Pstart[i], pT->Pstart[i+1]-1.

       The function picks the largest part, splits it into two parts
       of nearly-equal weight, and repeats this until all k parts
       are small enough to satisfy the load balance criterion,
           weight[i] <= (1 + eps) * totweight / P for i=0, ..., k-1,
       or when P parts have been reached.

    */

    int *procs;     /* procs[i] = number of processors assigned to part i */
    int procslo;    /* number of processors for smaller part after split */
    int procshi;    /*    same for largest part */

    long totweight; /* total weight of whole matrix */
    long maxweight; /* maximum allowed weight of a matrix part */
    long *weight;   /* weight[i] = weight of part i */
    long weightlo;  /* smallest upper bound for part weight */
    long weighthi;  /* largest upper bound */
    
    struct ordertree RowTree, ColTree; /* For use with permutation, hierarchy */
    int Symmetric = FALSE; /* Indicates symmetric permutation */

    long **Nets;    /* For use with symmetric finegrain and permutation */

    long t;
    int i, j, k, dir=ROW, done;
   
    Nets = NULL;

    if (pT == NULL || pOptions == NULL) {
        fprintf(stderr, "DistributeMatrixMondriaan(): Null arguments!\n");
        return FALSE;
    }
 
    /* Use full PaToH partitioning if so desired. */
    if (pOptions->Partitioner == FullPaToH) return DistributeMatrixPaToH(pT, P, eps, pOptions);
    
    if( pT->Pstart == NULL ) {
        fprintf(stderr, "Warning (DistributeMatrixMondriaan): Pstart array was not initialised. Doing this now.\n");
        if( !PstartInit(pT, P) ) {
             fprintf(stderr, "DistributeMatrixMondriaan(): error during initialisation of Pstart!\n");
             return FALSE;
        }
    }
    
    if ((totweight = ComputeWeight(pT, 0, pT->NrNzElts-1, NULL, pOptions)) < 0) {
        fprintf(stderr, "DistributeMatrixMondriaan(): Unable to compute weight!\n");
        return FALSE;
    }
    
    /*
       Determine whether or not we are using symmetric partitioning.
       This is true when the matrix is symmetric and the UseSingleEntry option is set,
       *or* when the symmetric permute option is set (regardless of input matrix symmetry).
    */
    if (pOptions->SymmetricPermute || (
        (pT->MMTypeCode[3] == 'S' || pT->MMTypeCode[3] == 'K' || pT->MMTypeCode[3] == 'H') && 
        pOptions->SymmetricMatrix_UseSingleEntry == SingleEntYes)
       ) {
#ifdef INFO
        if (!pOptions->SymmetricPermute)
            fprintf(stderr, "Info: resulting permutation will be symmetric.\n");
#endif
        Symmetric = TRUE;
    } else
        Symmetric = FALSE;
    
    /* Initially all rows and columns contain nonzeros from a single part. */
    for (t = 0; t < pT->m; ++t) pT->RowLambda[t] = 1;
    for (t = 0; t < pT->n; ++t) pT->ColLambda[t] = 1;
    
    /* Mark all rows and columns as not being cut. */
    for (t = 0; t < pT->m; ++t) pT->RowMark[t] = 0;
    for (t = 0; t < pT->n; ++t) pT->ColMark[t] = 0;
    
    /* Allocate memory for storing the sorting depth and order of the splitting procedure. */
    if (pOptions->OrderPermutation != OrderNone) {
        if (pT->row_perm == NULL) {
            pT->row_perm = (long *)malloc(pT->m*sizeof(long));
            for (t = 0; t < pT->m; ++t) pT->row_perm[t] = t;
        }
        if (pT->row_perm_inv == NULL) {
            pT->row_perm_inv = (long *)malloc(pT->m*sizeof(long));
            for (t = 0; t < pT->m; ++t) pT->row_perm_inv[t] = t;
        }
        if (pT->col_perm == NULL) {
            pT->col_perm = (long *)malloc(pT->n*sizeof(long));
            for (t = 0; t < pT->n; ++t) pT->col_perm[t] = t;
        }
        if (pT->col_perm_inv == NULL) {
            pT->col_perm_inv = (long *)malloc(pT->n*sizeof(long));
            for (t = 0; t < pT->n; ++t) pT->col_perm_inv[t] = t;
        }
        
        pT->rowBoundaries = remembrance_init();
        
        if (Symmetric) {
            Nets = (long**)malloc(2*sizeof(long*)); /* Each vertex corresponds to two nets of the same type */
            Nets[0] = pT->i;
            Nets[1] = pT->j; /* That is, both rows and columns map back to the same nets */
            if (!GeneralCreateOrderTree(&RowTree, pT->m, 0, pT->NrNzElts, Nets, 2)) {
                fprintf(stderr, "DistributeMatrixMondriaan(): Unable to create row permutation tree!\n");
                return FALSE;
            }
        } else {
            /* Only initialise colBoundaries when not in symmetric mode */
            pT->colBoundaries = remembrance_init();
            /* Create regular order trees */
            if (!CreateOrderTree(&RowTree, pT->m, 0, pT->NrNzElts, pT->i)) {
                fprintf(stderr, "DistributeMatrixMondriaan(): Unable to create row permutation tree!\n");
                return FALSE;
            }
            if (!CreateOrderTree(&ColTree, pT->n, 0, pT->NrNzElts, pT->j)) {
                fprintf(stderr, "DistributeMatrixMondriaan(): Unable to create column permutation tree!\n");
                return FALSE;
            }
        }
    }
    
    /* Setup Mondriaan options. */
    maxweight = ((1 + eps) * totweight) / P;  /* rounded down */
    
    if (pOptions->SplitStrategy == OneDimRow)
        dir = ROW;
  
    if (pOptions->SplitStrategy == OneDimCol)
        dir = COL;
  
    if (pOptions->SplitStrategy == FineGrain)
        dir = FINEGRAIN;

    if (pOptions->SplitStrategy == SFineGrain)
        dir = SFINEGRAIN;
  
    if (pOptions->SplitStrategy == Alternate &&
        pOptions->Alternate_FirstDirection == FirstDirRow)
        dir = ROW;

    if (pOptions->SplitStrategy == Alternate &&
        pOptions->Alternate_FirstDirection == FirstDirCol)
        dir = COL;
  
    if (pOptions->SplitStrategy == Alternate &&
        pOptions->Alternate_FirstDirection == FirstDirRatio) {  
        if (pT->m > pT->n)
            dir = ROW;
        else if (pT->n > pT->m)
            dir = COL;
        else if (Random1(0,1) == 0) /* random tie-breaking */
            dir = ROW;
        else
            dir = COL;
    }

    procs = (int *) malloc(P * sizeof(int));
    weight = (long *) malloc(P * sizeof(long));
    
    if (procs == NULL || weight == NULL) {
        fprintf(stderr, "DistributeMatrixMondriaan(): Not enough memory!\n");
        return FALSE;
    }
  
    /**** Split until all parts are small enough ****/  
    done = FALSE;
    k = 1; /* k = current number of parts */ 
    procs[0] = P;
    weight[0] = totweight;
    
    while (done == FALSE && k < P) {
        /* If desired, the number of processors assigned to each
           part is adjusted to reflect the actual outcome of
           previous splittings. This may allow for a larger
           imbalance parameter in future splits, aiming
           at lower communication. There is a trade-off:
               Load imbalance <--> communication volume. */
        if (pOptions->LoadbalanceAdjust == AdjustYes) {
            if (!BalanceParts(weight, maxweight, P, k, procs)) {
                fprintf(stderr, "DistributeMatrixMondriaan(): Unable to balance parts!\n");
                return FALSE;
            }
        }

        if (!DetermineSplit(weight, maxweight, k, procs, &i, &weightlo, &weighthi, pOptions)) {
            fprintf(stderr, "DistributeMatrixMondriaan(): Unable to determine split!\n");
            return FALSE;
        }

        if (i == -1) { 
            /* no more need to split  */
            done = TRUE;
            break;
        }
        
        /* Part i will be split into i, i+1 and
           the old parts i+1,..,k-1 will be shifted into i+2,..,k.
           Pstart will be  adjusted accordingly (within SplitMatrix) */
#ifdef INFO2
        printf("  ******** Split part %d from %d parts******** \n", i, k);
#endif
        if (pOptions->ZeroVolumeSearch == ZeroVolYes && SplitMatrixZeroVolume(pT, k, i, weightlo, weighthi, pOptions)) {
#ifdef INFO
            printf("Found zero volume partition!\n");
#endif
        }
        else if (pOptions->SplitMethod == Simple) {
            /* Simple split of the matrix only based on load balance,
               not minimising communication volume. 
               Useful for testing and debugging */
            if (!SplitMatrixSimple(pT, k, i, weightlo, weighthi, pOptions)) {
                fprintf(stderr, "DistributeMatrixMondriaan(): Unable to split using simple splitting!\n");
                return FALSE;
            }
        }
        else if (pOptions->SplitMethod == KLFM) {
            if (!SplitMatrixKLFM(pT, k, i, dir, weightlo, weighthi, pOptions)) {
                fprintf(stderr, "DistributeMatrixMondriaan(): Unable to split using KLFM!\n");
                return FALSE;
            }
        }
        else {
            fprintf(stderr, "DistributeMatrixMondriaan(): Unknown SplitMethod!\n");
            return FALSE;
        }
#ifdef INFO2
        printf("  Pstart[%d] = %ld ", i,  pT->Pstart[i]);
        printf("Pstart[%d] = %ld ", i+1,  pT->Pstart[i+1]);
        printf("Pstart[%d] = %ld ", i+2,  pT->Pstart[i+2]);
        printf("\n");
#endif
        /* Update lambdas for all rows and columns. */
        if (TRUE) {
            long countrows = 0, countcols = 0;
            
            /* Mark rows/columns in first part. */
            for (t = pT->Pstart[i]; t < pT->Pstart[i + 1]; ++t) {
                pT->RowMark[pT->i[t]] = 2;
                pT->ColMark[pT->j[t]] = 2;
            }
            
            /* Check if these rows/columns also occur in the second part and if so, increase their lambda. */
            for (t = pT->Pstart[i + 1]; t < pT->Pstart[i + 2]; ++t) {
                /* Make sure we only count each row/column once. */
                if (pT->RowMark[pT->i[t]] == 2) {
                    pT->RowMark[pT->i[t]] = 0;
                    pT->RowLambda[pT->i[t]]++;
                    countrows++;
                }
                if (pT->ColMark[pT->j[t]] == 2) {
                    pT->ColMark[pT->j[t]] = 0;
                    pT->ColLambda[pT->j[t]]++;
                    countcols++;
                }
            }
            
            /* Clear remaining flags. */
            for (t = pT->Pstart[i]; t < pT->Pstart[i + 1]; ++t) {
                pT->RowMark[pT->i[t]] = 0;
                pT->ColMark[pT->j[t]] = 0;
            }

#ifdef INFO2
            printf("DistributeMatrixMondriaan(): Updated lambdas of %ld rows and %ld columns.\n", countrows, countcols);
#endif
        }
        
        /* Keep track of the way the collection of nonzeros is being split and permute matrix. */
        if (pOptions->OrderPermutation != OrderNone) {
            if (!Symmetric) {
                OrderTreeSplit(&RowTree, pT->Pstart[i], pT->Pstart[i + 1], pT->Pstart[i + 2], pT->i, pOptions, pT->rowBoundaries);
                OrderTreeSplit(&ColTree, pT->Pstart[i], pT->Pstart[i + 1], pT->Pstart[i + 2], pT->j, pOptions, pT->colBoundaries);
            } else {
                Nets[0] = pT->i;
                Nets[1] = pT->j;
                GeneralOrderTreeSplit(&RowTree, pT->Pstart[i], pT->Pstart[i + 1], pT->Pstart[i + 2], Nets, 2l, pOptions, pT->rowBoundaries);
            }
            
#ifdef INFO2
            /* Check whether we actually generate permutations. */
            for (t = 0; t < pT->m; ++t) pT->row_perm[t] = 0;
            for (t = 0; t < pT->m; ++t) pT->row_perm[RowTree.Pi[t]] = 1;
            for (t = 0; t < pT->m; ++t) if (pT->row_perm[t] != 1) fprintf(stderr, "DistributeMat(): Faulty row permutation!\n");
            if (!Symmetric) {
                for (t = 0; t < pT->n; ++t) pT->col_perm[t] = 0;
                for (t = 0; t < pT->n; ++t) pT->col_perm[ColTree.Pi[t]] = 1;
                for (t = 0; t < pT->n; ++t) if (pT->col_perm[t] != 1) fprintf(stderr, "DistributeMat(): Faulty column permutation!\n");
            }
#endif
            
            for (t = 0; t < pT->m; ++t) {
                pT->row_perm[t] = RowTree.Pi[t];
                pT->row_perm_inv[RowTree.Pi[t]] = t;
            }
            if (!Symmetric)
                for (t = 0; t < pT->n; ++t) {
                    pT->col_perm[t] = ColTree.Pi[t];
                    pT->col_perm_inv[ColTree.Pi[t]] = t;
                }
            else { /* Copy from row */
                for (t = 0; t < pT->n; ++t) {
                    pT->col_perm[t]     = pT->row_perm[t];
                    pT->col_perm_inv[t] = pT->row_perm_inv[t];
                }
            }
        }
        
        /* Shift weight and procs */
        for (j = k; j > i+1; j--) {
            weight[j] = weight[j-1];
            procs[j] = procs[j-1];
        }

        k++; /* new number of parts */

        weight[i] = ComputeWeight(pT, pT->Pstart[i], pT->Pstart[i+1]-1, NULL, pOptions);
        weight[i+1] = ComputeWeight(pT, pT->Pstart[i+1], pT->Pstart[i+2]-1, NULL, pOptions);
        
        if (weight[i] < 0 || weight[i + 1] < 0) {
            fprintf(stderr, "DistributeMatrixMondriaan(): Unable to compute weights!\n");
            return FALSE;
        }
        
        procslo = procs[i]/2;
        procshi = (procs[i]%2==0 ? procslo : procslo+1);
        
        if (weight[i] <= weight[i+1]) {
            procs[i] = procslo;
            procs[i+1] = procshi;
        } else { 
            procs[i] = procshi;
            procs[i+1] = procslo;
        }


        /* Check if there is a part that is too large */
        done = TRUE;
        for (j = 0; j < k; j++) {
            if (weight[j] > maxweight) {
                done = FALSE;
                break;
            }
        }
  
        /* Alternate between split directions, each time after
           the number of parts reaches a power of two.
           This strategy resembles the strategy in Mondriaan version 1 
           for powers of two, but it is not identical.
           (The part to be cut is the largest, and not the one
           determined by the numbering as in version 1.) */
        if (pOptions->SplitStrategy == Alternate && logb2(k+1) != logb2(k))
            dir = ! dir;
	
        /* Issue callback. */
        if (Callback != NULL) {
            if (!Callback(k, i, pT)) break;
        }
    }
    
#ifdef INFO2
    /* Check if each part is small enough */
    for (j = 0; j < k; j++)
        if (weight[j] > maxweight) 
            printf("WARNING: part %d too large \n",j);
#endif

    /* Set Pstart for empty parts  */
    for (j = k+1; j <= P; j++)
        pT->Pstart[j] = pT->Pstart[k];

#ifdef INFO2
    printf("  Number of parts = %d \n", k);
    printf("  Pstart = ");
    for (j = 0; j <= P; j++)
        printf("%ld ", pT->Pstart[j]);
    printf("\n\n");
#endif      

#ifdef INFO2
    /* Print all lambdas. */
    printf("  Row lambda histogram:\n");
    VerifyLambdas(pT->RowLambda, pT->m, P);
    printf("  Column lambda histogram:\n");
    VerifyLambdas(pT->ColLambda, pT->n, P);
#endif
    
    if(pOptions->CheckUpperBound == CheckUpperBoundYes) {
        /* Compute volume. It should be at most (min(m,n)+1)(P-1) */
        long ComVol1, ComVol2, tmp;
        CalcCom(pT, NULL, (pT->m < pT->n)?ROW:COL, &ComVol1, &tmp, &tmp, &tmp, &tmp);
        CalcCom(pT, NULL, (pT->m < pT->n)?COL:ROW, &ComVol2, &tmp, &tmp, &tmp, &tmp);
        long upperBound = (((pT->m < pT->n)?pT->m:pT->n)+1)*(P-1);
        
        if(ComVol1+ComVol2 > upperBound) {
#ifdef INFO
            printf("Info: Achieved volume %ld is larger than upper bound %ld. Attempting to generate upper bound solution.\n", ComVol1+ComVol2, upperBound);
#endif
            if (!SplitMatrixUpperBound(pT, P, pOptions)) {
                fprintf(stderr, "DistributeMatrixMondriaan(): Unable to compute upper bound solution!\n");
            }
            k = P;
        }
    }
    

    /* Set matrix type code to distributed */
    pT->MMTypeCode[0] = 'D';


    /* Free local arrays */
    if (pOptions->OrderPermutation != OrderNone) {
        /* First the row permutations */
        /* Save hierarchy info into remembrance struct first */
        if( !remembrance_combine( pT->rowBoundaries, RowTree.Ranges, RowTree.NrRanges, 1 ) ) {
            fprintf(stderr, "DistributeMatrixMondriaan: Error deriving block hierarchy (row direction)!\n");
            return FALSE;
        }
        /* Add last bounds */
        remembrance_add( pT->rowBoundaries, pT->m, ULONG_MAX, ULONG_MAX );
        pT->rowBoundaries->vector[pT->rowBoundaries->size-1].id = ULONG_MAX;
        pT->rowBoundaries->vector[pT->rowBoundaries->size-1].parent = ULONG_MAX;
        /* Then delete */
        DestroyOrderTree(&RowTree);

        /* Now column permutations, if applicable */
        if (!Symmetric) {
            /* Save hierarchy info into remembrance struct first */
            if( !remembrance_combine( pT->colBoundaries, ColTree.Ranges, ColTree.NrRanges, 1 ) ) {
                fprintf(stderr, "DistributeMatrixMondriaan: Error deriving block hierarchy (column direction)!\n");
                return FALSE;
            }
            /* Add last bounds */
            remembrance_add( pT->colBoundaries, pT->n, ULONG_MAX, ULONG_MAX );
            pT->colBoundaries->vector[pT->colBoundaries->size-1].id = ULONG_MAX;
            pT->colBoundaries->vector[pT->colBoundaries->size-1].parent = ULONG_MAX;
            /* Then delete */
            DestroyOrderTree(&ColTree);
        } else
            free(Nets);
    }
    
    free(weight);
    free(procs);

    return TRUE;
} /* end DistributeMatrixMondriaan */  


long ComputeWeight(const struct sparsematrix *pT, long lo, long hi, long *wnz, const struct opts *pOptions) {
    /* This function computes the weight of the matrix part containing
       the nonzeros in positions lo..hi. Each nonzero has weight 1,
       except the dummies which have weight 0. 

       In case of a symmetric, skew-symmetric, or hermitian matrix
       and a single entry representing a pair of nonzeros,
       the off-diagonal nonzeros have weight 2.

       If wnz <> NULL, the weights of the individual nonzeros are also
       returned: wnz[t] = weight of nonzero lo+t, for 0 <= t <= hi-lo. 

       In case of a column-weighted matrix, the returned weight is the total
       weight of the columns that have a nonzero in a position lo..hi.
       No individual weights are returned.

   */

    int *marked;
    long t, j, weight, Ndum, Ndiag;
    
    if (!pT) {
        fprintf(stderr, "ComputeWeight(): Null argument!\n");
        return -1;
    }
  
    if (lo > hi) return 0;
    
    /*
     DBG
    fprintf(stdout, "NrColWeights: %ld\n", pT->NrColWeights);
    fprintf(stdout, "MMTypeCode: %s\n", pT->MMTypeCode[0]);
    fflush(stdout);
    */
   
    if (pT->MMTypeCode[0] == 'W' && pT->NrColWeights > 0) {
        /* The weight is determined by the weights 
           of columns containing nonzeros and not by the nonzeros */
        
        /* Mark processed columns to avoid counting their weight twice */
        marked = (int *) malloc(pT->NrColWeights * sizeof(int));
        
        if (marked == NULL) {
            fprintf(stderr, "ComputeWeight(): Not enough memory!\n");
            return -1;
        }
        
        for (j = 0; j < pT->NrColWeights; j++)
            marked[j] = FALSE;
    
        weight = 0;
        for (t = lo; t <= hi; t++) {
            j = pT->j[t];
            if (marked[j]==FALSE) {
                weight +=  pT->ColWeights[j]; 
                marked[j] = TRUE;
            }
        }
        
        free(marked);
        
    } else {
                        
        weight = hi - lo + 1;
    
        /* Initialise individual weights to standard value */
        if (wnz != NULL)
            for (t = lo; t <= hi; t++)
                wnz[t-lo] = 1;

        /* Count the dummies and set their individual weights to 0 */
        Ndum = 0;
        if (pT->m == pT->n && pT->NrDummies > 0) { 
            for (t = lo; t <= hi; t++) {
                if (pT->i[t] == pT->j[t] && pT->dummy[pT->i[t]]) {
                      /* nonzero t is a dummy */
                      Ndum++;
                      if (wnz != NULL) 
                          wnz[t-lo] = 0; 
                }
            }
            
            weight -= Ndum; /* weight = number of nondummy nonzeros */
        }
  
        /* Count the number of diagonal entries and set the individual
           weights of off-diagonal entries to 2 */
        if (pT->m == pT->n &&
             (pT->MMTypeCode[3]=='S' || pT->MMTypeCode[3]=='K' ||
              pT->MMTypeCode[3]=='H') &&
             pOptions->SymmetricMatrix_UseSingleEntry == SingleEntYes) {
            Ndiag = 0;
            for (t = lo; t <= hi; t++) {
                if (pT->i[t] == pT->j[t]) 
                    Ndiag++; /* individual weight is already OK */
                else if (wnz != NULL)
                    wnz[t-lo] = 2;   
            }

            /* Double the total weight, except for the Ndiag - Ndum
               nondummy nonzeros on the diagonal */
            weight = 2*weight - (Ndiag - Ndum); 
        }
    }

    return weight;
} /* end ComputeWeight */


int BalanceParts (long *weight, long maxweight, int P, int k, int *procs) {
    /* This function balances the parts of a partitioning by assigning
       a number of processors to each part in proportion to its weight. 
       If a part has weight <= maxweight, its number of processors will be set
       to -1 meaning that it need not be split anymore.

       Input: weight[i] = weight of part i, weight[i] >= 0, 0 <= i < k,
              maxweight = maximum allowed weight per processor
                          after partitioning, maxweight >= 0,
              P = total number of processors available,
              k = number of parts, 1 <= k <= P. 
       Output: procs[i] = number of processors assigned to part i. */
 
    long totweight, *surplus, *I;
    int i, j, nfinished, nextra, p, q;
    
    if (!weight || !procs) {
        fprintf(stderr, "Null arguments!\n");
        return FALSE;
    }
  
    if (k <= 0)
        return TRUE; /* no need for balancing */
    else if (k == 1) {
        if (weight[0] <= maxweight)
            procs[0] = -1;
        else
            procs[0] = P;
    } else {
        /* Assign finished parts and compute total remaining weight */
        totweight = 0;
        nfinished = 0; /* number of finished parts */
        
        for (i = 0; i < k; i++) {
            if (weight[i] <= maxweight) { 
                procs[i] = -1;
                nfinished++;
            } else {     /* weight[i] > 0 */
                procs[i] = 0; /* to be set later */
                totweight += weight[i];
            }
        }
        p = P - nfinished; /* number of processors still available,
                               i.e., not assigned to finished parts */
                
        surplus = (long *)malloc(k*sizeof(long));
        
        if (surplus == NULL) {
            fprintf(stderr, "BalanceParts(): Not enough memory!\n");
            return FALSE;
        }
        
        q = 0;
        
        for (i=0; i < k; i++) {
            if (procs[i] == 0 && totweight > 0) {
                /* Use long longs for intermediate results to prevent overflow. */
                procs[i] = ((long long)weight[i] * (long long)p) / totweight;
                           /* = weight[i] / avgweight */

                surplus[i] = (long long)weight[i] * (long long)p -
                             (long long)procs[i] * (long long)totweight;
                             /* = (weight[i] - procs[i]*avgweight) * p.
                             The scaling by a factor p is to obtain integers. */

                q += procs[i]; /* total number of processors assigned so far
                                   for unfinished parts */
            } else {
                surplus[i] = 0;
            }
        }
        
        I = QSort(surplus, (long)k);
        
        if (I == NULL) {
            fprintf(stderr, "BalanceParts(): Unable to sort surplus!\n");
            return FALSE;
        }
        
        nextra = 0;
        
        for (j=0; j < k; j++) {
            i = I[j];
            if (procs[i] >= 0 && nextra < p - q) {
                /* Assign an extra processor to this part because of a large surplus.
                   The number of extra processors available is p-q.  */
                procs[i]++;
                nextra++;
            }
        }
       
        free(I);
        free(surplus);

    }
    
    return TRUE;
} /* end BalanceParts */


int DetermineSplit(long *weight, long maxweight, int k, int *procs, int *isplit, long *weightlo, long *weighthi, const struct opts *pOptions) {
    /* This function determines the part to split in the next bipartitioning
       and computes upper bounds for the weight of the smallest and largest
       parts resulting from the split. 
       
       The weight bounds are based on a split into equal parts if procs[isplit]
       is even, and nearly equal parts if it is odd. All splits of part isplit
       that remain to be done allow the weight imbalance to grow with a constant
       factor, or, alternatively, by an increasing or decreasing factor.
    .  The factor is chosen such that all parts resulting from isplit
       at the end of the partitioning have weight <= maxweight.
       
       Input: weight[i] = weight of part i, weight[i] >= 0, 0 <= i < k,
              maxweight = maximum allowed weight per processor
                          after partitioning, maxweight >= 0,
              k = number of parts,
              procs[i] = number of processors assigned to part i,
                         or -1 if part need not be split anymore.
              
       Output: isplit = index of part to be split,
                        or -1 if no need to split any more.
               weightlo = smallest upper bound for part weight,
               weighthi = largest upper bound for part weight.
               procs is modified: some finished parts are marked by -1
                     and the rest must be reinitialised afterwards. */
       
    int plo, phi, qlo, qhi, i, ilarge;
    long wlarge;
    double wlo, whi, eps, deltalo, deltahi;
       
    if (!weight || !procs || !isplit || !weightlo || !weighthi || !pOptions) {
        fprintf(stderr, "DetermineSplit(): Null arguments!\n");
        return FALSE;
    }
    
    if (k <= 0) {
        *isplit = -1;
        return TRUE; /* no need for splitting */
    }
       
    /* Find part of largest weight */
    wlarge = -1;
    ilarge = -1;
    for (i = 0; i < k; i++) {
        if (weight[i] > wlarge && procs[i] >= 0) {
             wlarge = weight[i];
             ilarge = i;
        }
    }
    i = ilarge;
    if (wlarge == -1)
        *isplit = -1;
    else {
        *isplit = i;
        if (weight[i] > procs[i]*maxweight) {
            /* Desired split is infeasible, so grab an extra processor,
               which is almost always available. */
            procs[i]++;
#ifdef INFO2
            printf("    Warning: extra processor grabbed for part %d\n",i);
#endif
        }
        if (weight[i] > procs[i]*maxweight) {
            /* Desired split is still infeasible */
            fprintf(stderr, "DetermineSplit(): No feasible split possible!\n");
            return FALSE;
        }

        if (procs[i]==1) {
            /* No split needed. Part is feasible, so mark it as finished. */
            procs[i] = -1;
            return TRUE;
        }

        plo = procs[i]/2;
        phi = (procs[i]%2==0 ? plo : plo+1);
        
        /* maximum number of splits (including current split)
           needed to reach leaf of splitting tree */
        qlo = logb2(plo) + 1; 
        qhi = logb2(phi) + 1;
      
        eps =  (maxweight*procs[i])/(double)weight[i] - 1.0; /* eps >= 0 */

        deltalo = eps; /* default value for plo = 1 */
        deltahi = eps;
           
        /* Splitting imbalance delta(0) for part i is determined by inequality
                (1+delta(0))...(1+delta(q-1)) weight[i]/procs[i] <= maxweight,
           i.e. (1+delta(0))...(1+delta(q-1)) <= 1 + eps. */
    
        if (pOptions->LoadbalanceStrategy == Decrease) {
            /* take delta(0) large */
            if (plo > 1)
                deltalo = eps/2;
            if (phi > 1)
                deltahi = eps/2;
        } else if (pOptions->LoadbalanceStrategy == Increase) {
            /* take delta(0) small */ 
            if (plo > 1)
                deltalo = eps/(2*plo);
            if (phi > 1)
                deltahi = eps/(2*phi);
        } else if (pOptions->LoadbalanceStrategy == Constant) {
            /* take delta(0) assuming constant growth */
            if (plo > 1)
                deltalo = eps/qlo;
            if (phi > 1)
                deltahi = eps/qhi;
        } else {
            fprintf(stderr, "DetermineSplit(): Unknown LoadbalanceStrategy!\n");
            return FALSE;
        }

        wlo = (weight[i] / (double)procs[i]) * plo;
        whi = (weight[i] / (double)procs[i]) * phi;
        *weightlo = (long) ((1.0 + deltalo) * wlo);
        *weighthi = (long) ((1.0 + deltahi) * whi);
    }
 
    return TRUE;
 } /* end DetermineSplit */


int logb2(int n) {
    /* This function computes the base-2 log of n, n >= 1,
       rounded up to the nearest integer. */
    int k, i;
    char roundup;

    k = 0;
    i = n; /* n = (b(m-1) ...b(0)) in binary notation
                                   with b(0) the least significant bit
                                   and b(m-1) = 1 the most significant bit */
    roundup = FALSE;

    while (i >= 2) {
        /* i = (b(m-1) ...b(k)) */
        if (i%2 == 1)
            roundup = TRUE;
        i /= 2;
        k++;
    }
    /* i = 1 = (b(m-1)),
       k = number of divisions performed = m-1*/

    if (roundup)
        k++;

    return k;
} /* end logb2 */


int SplitMatrixKLFM(struct sparsematrix *pT, int k, int i, int dir,
                     long weightlo, long weighthi, const struct opts *pOptions) {
    
    /* This function splits part i of the sparse matrix T into two parts,
       the first with weight <= weightlo and the second with weight <= weighthi.
       The split is performed by translating the matrix into a hypergraph
       and bipartitioning this hypergraph, trying to keep vertices
       from the same hyperedge together as much as possible.
       The bipartitioning employs a multilevel hypergraph partitioner
       based on the Kernighan-Lin/Fiduccia-Mattheyses algorithm.

       Input: T sparse matrix,
              k current number of parts, 1 <= k < P,
              i number of part to be split, 0 <= i < k,
              dir = direction of split. The input value is used,
                    except for the LocalBest, Hybrid strategies,
                    where the value is overruled by this function,
                    and for LocalRatio, where it is set elsewhere.
              weightlo = smallest upper bound for part weight, belongs to part 0
              weighthi = largest upper bound for part weight, belongs to part 1.
              
       Output: T sparse matrix.
               The nonzeros of the new part i (the first part) are in positions
                   pT->Pstart[i], pT->Pstart[i+1]-1.
               The nonzeros of the new part i+1 (the second) are in positions
                   pT->Pstart[i+1], pT->Pstart[i+2]-1.
               All parts > i+1 have been shifted.
       
    */

    long lo, hi, mid, nz, weight,reverseDir;
    int j, bestdir=ROW, bestcomm;
    struct sparsematrix A;
    struct biparthypergraph HG, HGRow, HGCol, HGFine;
    
    if (!pT || !pOptions) {
        fprintf(stderr, "SplitMatrixKLFM(): Null arguments!\n");
        return FALSE;
    }
    
    lo = pT->Pstart[i];
    hi = pT->Pstart[i+1]-1;
    
    /* Shift Pstart for parts > i */
    for (j = k; j > i; j--)
        pT->Pstart[j+1] = pT->Pstart[j];
        
    nz = hi-lo+1;
    
    if (nz == 0) 
        return TRUE; /* Pstart[i] = Pstart[i+1] */
    
    weight = ComputeWeight(pT, lo, hi, NULL, pOptions);
    
    if (weight > weightlo + weighthi || weight < 0) {
        /* Desired split is infeasible */
        fprintf(stderr, "SplitMatrixKLFM(): desired hypergraph split is infeasible!\n");
        return FALSE;
    }
     
    /* Copy info from T to A */ 
    A = *pT;            /* A has same size as T, and same other parameters, */
    A.i = &(pT->i[lo]); /* but only a subset of the nonzeros */
    A.j = &(pT->j[lo]);
    if (A.MMTypeCode[2] != 'P')
        A.ReValue = &(pT->ReValue[lo]);
    if (A.MMTypeCode[2] == 'C')
        A.ImValue = &(pT->ImValue[lo]);
    A.NrNzElts = nz;

    /*
    DBG
    fprintf(stdout, "KLFM: %ld==%ld\n", A.NrColWeights, pT->NrColWeights );
    */
  
    if (pOptions->SplitStrategy == LocalBest ||
         pOptions->SplitStrategy == Hybrid) {
         long bestvalidcomm = LONG_MAX;
         int bestvaliddir = -1;
         
#ifdef INFO2
        printf("**** Try split in row direction \n");
#endif    
        if (!SparseMatrix2BiPartHyperGraph(&A, ROW, pOptions, &HGRow)) {
            fprintf(stderr, "SplitMatrixKLFM(): Unable to convert matrix to hypergraph!\n");
            return FALSE;
        }
        
        if (!RunMLGraphPart(&HGRow, weightlo, weighthi, pOptions)) {
            fprintf(stderr, "SplitMatrixKLFM(): Unable to run KLFM (in row-direction)!\n");
            return FALSE;
        }
  
#ifdef INFO2
        printf("**** Try split in col direction \n");
#endif   
        if (!SparseMatrix2BiPartHyperGraph(&A, COL, pOptions, &HGCol)) {
            fprintf(stderr, "SplitMatrixKLFM(): Unable to convert matrix to hypergraph!\n");
            return FALSE;
        }
        
        if (!RunMLGraphPart(&HGCol, weightlo, weighthi, pOptions)) {
            fprintf(stderr, "SplitMatrixKLFM(): Unable to run KLFM (in column-direction)!\n");
            return FALSE;
        }
  
#ifdef INFO2
        printf("Communication: row=%ld col=%ld\n", HGRow.OptComm, HGCol.OptComm);
#endif
        
        if (HGRow.OptComm < HGCol.OptComm)
            bestdir = ROW;
        else if (HGCol.OptComm < HGRow.OptComm)
            bestdir = COL;
        else if (Random1(0,1) == 0) /* no preference, random tie-breaking */
            bestdir = ROW;
        else
            bestdir = COL;
        bestcomm = MIN(HGRow.OptComm, HGCol.OptComm);
        
        /* Select the best direction which still satisfies the balance constraints. */
        if (HGRow.WeightP[0] <= weightlo && HGRow.WeightP[1] <= weighthi && HGRow.OptComm < bestvalidcomm) {
            bestvalidcomm = HGRow.OptComm;
            bestvaliddir = ROW;
            
            if (HGCol.WeightP[0] <= weightlo && HGCol.WeightP[1] <= weighthi) {
                if (HGCol.OptComm < bestvalidcomm) {
                    bestvalidcomm = HGCol.OptComm;
                    bestvaliddir = COL;
                }
                else if (HGCol.OptComm == HGRow.OptComm) {
                    /* Random tie-breaking. */
                    bestvaliddir = (Random1(0, 1) == 0 ? ROW : COL);
                }
            }
        }
        else if (HGCol.WeightP[0] <= weightlo && HGCol.WeightP[1] <= weighthi && HGCol.OptComm < bestvalidcomm) {
            bestvalidcomm = HGCol.OptComm;
            bestvaliddir = COL;
        }

        if (pOptions->SplitStrategy == Hybrid) {

#ifdef INFO2
            printf("**** Try split in finegrain direction \n");
#endif
            if (!SparseMatrix2BiPartHyperGraph(&A, FINEGRAIN, pOptions, &HGFine)) {
                fprintf(stderr, "SplitMatrixKLFM(): Unable to convert matrix to hypergraph!\n");
                return FALSE;
            }
        
            if (!RunMLGraphPart(&HGFine, weightlo, weighthi, pOptions)) {
                fprintf(stderr, "SplitMatrixKLFM(): Unable to run KLFM (in finegrain mode)!\n");
                return FALSE;
            }

#ifdef INFO2
            printf("Communication: finegrain=%ld\n", HGFine.OptComm);
#endif
            if (HGFine.OptComm < bestcomm) {
                bestdir = FINEGRAIN;
                bestcomm = HGFine.OptComm;
            }
            /* in case of equality, we prefer the row or column direction */
            
            if (HGFine.WeightP[0] <= weightlo && HGFine.WeightP[1] <= weighthi && HGFine.OptComm < bestvalidcomm) {
                bestvalidcomm = HGFine.OptComm;
                bestvaliddir = FINEGRAIN;
            }
        }
        
        /* Override the overall best choice by the overall best valid choice (in case the balance constraint could not be satisfied for all splitting methods). */
        if (bestvaliddir != -1) {
#ifdef INFO2
            if (bestvalidcomm > bestcomm) printf("Warning, the imbalance constraints were not met in all split directions!\n");
#endif
            bestdir = bestvaliddir;
            bestcomm = bestvalidcomm;
        }
#ifdef INFO
        else {
            printf("Warning, no splitting direction meets the imbalance constraints!\n");
        }
#endif
        
#ifdef INFO2
        if (bestdir == ROW)
            printf("**** Row direction chosen\n");
        else if (bestdir == COL)
            printf("**** Col direction chosen\n");
        else
            printf("**** Finegrain direction chosen\n");
#endif   
  
        if (bestdir == ROW) {
            if (!BiPartHyperGraph2SparseMatrix(&HGRow, lo, hi, &mid, pT)) {
                fprintf(stderr, "SplitMatrixKLFM(): Unable to convert hypergraph to matrix!\n");
                return FALSE;
            }
        }
        else if (bestdir == COL) {
            if (!BiPartHyperGraph2SparseMatrix(&HGCol, lo, hi, &mid, pT)) {
                fprintf(stderr, "SplitMatrixKLFM(): Unable to convert hypergraph to matrix!\n");
                return FALSE;
            }
        }
        else {
            if (!BiPartHyperGraph2SparseMatrix(&HGFine, lo, hi, &mid, pT)) {
                fprintf(stderr, "SplitMatrixKLFM(): Unable to convert hypergraph to matrix!\n");
                return FALSE;
            }
        }
  
        if (pOptions->SplitStrategy == Hybrid) 
            DeleteBiPartHyperGraph(&HGFine);
        DeleteBiPartHyperGraph(&HGCol);
        DeleteBiPartHyperGraph(&HGRow);
     } else if(pOptions->SplitStrategy == MediumGrain){

#ifdef INFO2
        printf("**** Split extended mediumgrain matrix in column direction \n");
#endif

        A.mgMid = hi+1;
        A.mgDir=0;
        if(!CreateInitialMediumGrainDistribution(&A, &(A.mgMid))){
            fprintf(stderr, "SplitMatrixKLFM(): Unable to perform initial mediumgrain split!\n");
            return FALSE;
        }

        /* Translate to a hypergraph by 1D columns, bipartition, and translate back */
        if (!SparseMatrix2BiPartHyperGraph(&A, MEDIUMGRAIN, pOptions, &HG)) {
            fprintf(stderr, "SplitMatrixKLFM(): Unable to convert mediumgrain matrix to hypergraph!\n");
            return FALSE;
        }

        if (!RunMLGraphPart(&HG, weightlo, weighthi, pOptions)) {
            fprintf(stderr, "SplitMatrixKLFM(): Unable to run HKLFM!\n");
            return FALSE;
        }
        if (!BiPartHyperGraph2SparseMatrix(&HG, lo, hi, &mid, pT)) {
            fprintf(stderr, "SplitMatrixKLFM(): Unable to convert hypergraph to matrix!\n");
            return FALSE;
        }

        bestcomm = HG.OptComm;
        DeleteBiPartHyperGraph(&HG);
        
    } else {

#ifdef INFO2
        if (pOptions->SplitStrategy != LocalRatio) {
            if (dir == ROW)
                printf("**** Split in row direction \n");
            else
                printf("**** Split in col direction \n");
        } /* For LocalRatio, dir will be set inside
             SparseMatrix2BiPartHyperGraph */
#endif
  
        if (!SparseMatrix2BiPartHyperGraph(&A, dir, pOptions, &HG)) {
            fprintf(stderr, "SplitMatrixKLFM(): Unable to convert matrix to hypergraph!\n");
            return FALSE;
        }
        
        if (!RunMLGraphPart(&HG, weightlo, weighthi, pOptions)) {
            fprintf(stderr, "SplitMatrixKLFM(): Unable to run HKLFM!\n");
            return FALSE;
        }
        
        if (!BiPartHyperGraph2SparseMatrix(&HG, lo, hi, &mid, pT)) {
            fprintf(stderr, "SplitMatrixKLFM(): Unable to convert hypergraph to matrix!\n");
            return FALSE;
        }
        bestcomm = HG.OptComm;
        DeleteBiPartHyperGraph(&HG);
    }
    
    /* Iterative refinement */
    if(pOptions->Iterative_Refinement==IR_Always || (pOptions->Iterative_Refinement==IR_After_MG && pOptions->SplitStrategy == MediumGrain)){
        reverseDir=0;
        while(1){
            A.mgMid = mid-lo;
            if(reverseDir==1){
                A.mgDir=(A.mgDir+1)%2;
            }
            
            if (!SparseMatrix2BiPartHyperGraph(&A, MEDIUMGRAIN, pOptions, &HG)) {
                fprintf(stderr, "SplitMatrixKLFM(): Unable to convert mediumgrain matrix to hypergraph!\n");
                return FALSE;
            }
            
            /* Set columns of large matrix to old partitioning */
            HG.WeightP[0]=0;
            HG.WeightP[1]=0;
            for(j=0;j<HG.NrVertices;j++){
                if(HG.Vtx2MatIndex[j]<A.n){
                    HG.V[j].partition=0;
                    HG.WeightP[0]+=HG.V[j].vtxwgt;
                }else{
                    HG.V[j].partition=1;
                    HG.WeightP[1]+=HG.V[j].vtxwgt;
                }
            }
            
            /* Run single KLFM improvement */
            if (!RunHKLFM(&HG,weightlo,weighthi,1,pOptions)) {
                fprintf(stderr, "SplitMatrixKLFM(): Unable to run HKLFM!\n");
                return FALSE;
            }
            if (!BiPartHyperGraph2SparseMatrix(&HG, lo, hi, &mid, pT)) {
                fprintf(stderr, "SplitMatrixKLFM(): Unable to convert hypergraph to matrix!\n");
                return FALSE;
            }
            if(HG.OptComm==bestcomm){
                if(reverseDir==1){
                    DeleteBiPartHyperGraph(&HG);
                    break;
                }
                reverseDir++;
            }else{
                reverseDir=0;
            }
            bestcomm = HG.OptComm;
            DeleteBiPartHyperGraph(&HG);
        }
    }

    /* register new splitting point */ 
    pT->Pstart[i+1] = mid;

    return TRUE;
} /* end SplitMatrixKLFM */


int SplitMatrixSimple(struct sparsematrix *pT, int k, int i,
                       long weightlo, long weighthi, const struct opts *pOptions) {
    /* This function splits part i of the sparse matrix T into two parts,
       the first with weight <= weightlo and the second with weight <= weighthi.
       The split is performed by simply cutting somewhere in the middle
       of the array of nonzeros. Nonzeros are not moved.
       Input: T sparse matrix,
              k current number of parts, 1 <= k < P,
              i number of part to be split, 0 <= i < k,
              weightlo = smallest upper bound for part weight, belongs to part 0,
              weighthi = largest upper bound for part weight, belongs to part 1.
              
       Output: T sparse matrix.
               The nonzeros of the new part i (the first part) are in positions
                   pT->Pstart[i], pT->Pstart[i+1]-1.
               The nonzeros of the new part i+1 (the second) are in positions
                   pT->Pstart[i+1], pT->Pstart[i+2]-1.
               All parts > i+1 have been shifted.
               
       NB: This function works only if the weights are based on the number
       of nonzeros.  Dummy nonzeros with  weight 0 are allowed, and so are
       nonzeros with weight 2 caused by a symmetric matrix with single-entry
       storage. The function does not work, however, for weights based on
       matrix columns.
       
    */
       
    long lo, hi, nz, weight, surplus, ideal, wtot, t, *wnz;
    int j;

    if (!pT || !pOptions) {
        fprintf(stderr, "SplitMatrixSimple(): Null arguments!\n");
        return FALSE;
    }
    
    if (pT->MMTypeCode[0] == 'W' && pT->NrColWeights > 0) {
        fprintf(stderr, "SplitMatrixSimple(): Cannot be used with column-weighted matrix!\n");
        return FALSE;
    }

    lo = pT->Pstart[i];
    hi = pT->Pstart[i+1]-1;
    
    /* Shift Pstart for parts > i */
    for (j = k; j > i; j--)
        pT->Pstart[j+1] = pT->Pstart[j];
        
    nz = hi-lo+1;
    if (nz == 0) 
        return TRUE; /* Pstart[i] = Pstart[i+1] */
        
    /* Allocate memory for storing weights 0, 1, 2 of individual nonzeros */
    wnz = (long *) malloc(nz * sizeof(long));
    
    if (wnz == NULL) {
        fprintf(stderr, "SplitMatrixSimple(): Not enough memory!\n");
        return FALSE;
    }
 
    weight = ComputeWeight(pT, lo, hi, wnz, pOptions);
       
    if (weight > weightlo + weighthi || weight < 0) {
        /* Desired split is infeasible */
        fprintf(stderr, "SplitMatrixSimple(): desired simple split is infeasible!\n");
        return FALSE;
    }
    
    /* The surplus is equally divided over the two parts */
    surplus = weightlo + weighthi - weight;
    ideal = MAX(weightlo - surplus/2,0); /* ideal weight of new part i */
    
    /* Find first legal split */
    wtot = 0;
    t = 0;
    
    while (t < nz && wtot + weighthi < weight) {
        wtot += wnz[t];
        t++; /* t = number of nonzeros included in smallest part */
    }
    
    if (wtot > weightlo) {
        /* Desired simple split is infeasible, because a nonzero must be split.
           This is only likely to occur if bounds are too tight, and 
           then only if weights of 2 are used.*/
        fprintf(stderr, "SplitMatrixSimple(): simple split causes split of a nonzero!\n");
        return FALSE;
    }
         
    /* If possible, improve until ideal split */
    while (t < nz && wtot + wnz[t] <= ideal) {
        wtot += wnz[t];
        t++;
    }  
     
    pT->Pstart[i+1] = pT->Pstart[i] + t;
    
    free(wnz);
    
    return TRUE;
} /* end SplitMatrixSimple */


int SplitMatrixZeroVolume(struct sparsematrix *pT, int k, int i,
                     long weightlo, long weighthi, const struct opts *pOptions) {
    
    /* This function splits part i of the sparse matrix T into two parts,
       the first with weight <= weightlo and the second with weight <= weighthi.
       The split is performed by searching for a split with zero communication
       volume. If such a split is found, it is applied to pT. Otherwise, pT is
       left untouched.

       Input: T sparse matrix,
              k current number of parts, 1 <= k < P,
              i number of part to be split, 0 <= i < k,
              weightlo = smallest upper bound for part weight, belongs to part 0
              weighthi = largest upper bound for part weight, belongs to part 1.
              
       Output: T sparse matrix.
               The following applies if a zero volume split is found:
               The nonzeros of the new part i (the first part) are in positions
                   pT->Pstart[i], pT->Pstart[i+1]-1.
               The nonzeros of the new part i+1 (the second) are in positions
                   pT->Pstart[i+1], pT->Pstart[i+2]-1.
               All parts > i+1 have been shifted.
       
    */

    long lo, hi, mid = 0, nz, weight;
    int j;
    struct sparsematrix A;
    
    if (!pT || !pOptions) {
        fprintf(stderr, "SplitMatrixZeroVolume(): Null arguments!\n");
        return FALSE;
    }
    
    lo = pT->Pstart[i];
    hi = pT->Pstart[i+1]-1;
    
    nz = hi-lo+1;
    
    if (nz > 0) {
        weight = ComputeWeight(pT, lo, hi, NULL, pOptions);
        
        if (weight > weightlo + weighthi || weight < 0) {
            /* Desired split is infeasible */
            fprintf(stderr, "SplitMatrixZeroVolume(): desired split is infeasible!\n");
            return FALSE;
        }
        
        /* Copy info from T to A */
        A = *pT;            /* A has same size as T, and same other parameters, */
        A.i = &(pT->i[lo]); /* but only a subset of the nonzeros */
        A.j = &(pT->j[lo]);
        if (A.MMTypeCode[2] != 'P')
            A.ReValue = &(pT->ReValue[lo]);
        if (A.MMTypeCode[2] == 'C')
            A.ImValue = &(pT->ImValue[lo]);
        A.NrNzElts = nz;
        
#ifdef TIME
        clock_t starttime, endtime;
        double cputime;
        starttime = clock();
#ifdef UNIX
        struct timeval starttime1, endtime1;
        gettimeofday(&starttime1, NULL);
#endif
#endif
        /* Run zero volume search. If a zero volume split is found, ZeroVolumeSearch()
         * will apply this split directly; we then only need to update Pstart.
         */
        int foundZeroVolumePartition = ZeroVolumeSearch(&A, weightlo, weighthi, &mid, pOptions);
        
#ifdef TIME
        endtime = clock();
        cputime = ((double) (endtime - starttime)) / CLOCKS_PER_SEC;
        printf("  matrix distribution zeroVolumeSearch CPU-time    : %f seconds\n", cputime);
#ifdef UNIX
        gettimeofday(&endtime1, NULL);
        printf("  matrix distribution zeroVolumeSearch elapsed time: %f seconds\n",
                (endtime1.tv_sec - starttime1.tv_sec) +
                (endtime1.tv_usec - starttime1.tv_usec) / 1000000.0);
#endif
        fflush(stdout);
#endif

        if(!foundZeroVolumePartition) {
            return FALSE;
        }
    }
    else {
        mid = 0; /* Pstart[i] = Pstart[i+1] */
    }
    
    /* Shift Pstart for parts > i */
    for (j = k; j > i; j--)
        pT->Pstart[j+1] = pT->Pstart[j];
    
    /* Register new splitting point.
     * mid lies in [0,hi-lo], translate it to [lo,hi] */
    pT->Pstart[i+1] = lo + mid;

    return TRUE;
} /* end SplitMatrixZeroVolume */
