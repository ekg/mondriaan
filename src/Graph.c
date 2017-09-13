#include "Graph.h"

int CreateNewBiPartHyperGraph(long NrVertices, long NrNets,
                                                  long NrPins,
                                                  int StoreMat, char MType, struct biparthypergraph *pHG) { 

    /* This function creates a new bipartitioned hypergraph
       with NrVertices vertices, NrNets nets (hyperedges),
       and NrPins pins (nonzeros, if the hypergraph represents a sparse matrix).

       If the boolean StoreMat is true, the numerical values of the pins
       (nonzeros) are stored as well:
           - a real part, if the matrix is not pattern-only, i.e. MType <> 'P';
           - an imaginary part if the matrix is complex, i.e. MType = 'C'.

       The function initialises all variables of the hypergraph to 0,
       sets all pointers to NULL, and it allocates memory to store the nets,
       vertices, and pins.
    */
    register long t;
    
    if (!pHG) {
        fprintf(stderr, "CreateNewBiPartHyperGraph(): Null argument provided!\n");
        return FALSE;
    }

    pHG->NrVertices = NrVertices;
    pHG->NrNets = NrNets;
    pHG->NrPins = NrPins;

    /* Initialise array pointers to NULL */
    pHG->V = NULL;
    pHG->VtxMoveLog = NULL;
    pHG->OptPartVtx = NULL;
    pHG->Vtx2MatIndex = NULL;
    pHG->N = NULL;
    pHG->Net2MatIndex = NULL;
    pHG->VtxAdjncy = NULL;
    pHG->NetAdjncy = NULL;
    pHG->MatReValue = NULL;
    pHG->MatImValue = NULL;

    /* Allocate memory for the vertices */
    if (pHG->NrVertices > 0) {
        pHG->V = (struct vertex *) malloc(pHG->NrVertices * sizeof(struct vertex));
        pHG->VtxMoveLog = (long *) malloc(pHG->NrVertices * sizeof(long));
        pHG->OptPartVtx = (int *) malloc(pHG->NrVertices * sizeof(int));

        if (pHG->V == NULL || pHG->VtxMoveLog == NULL || pHG->OptPartVtx == NULL) {
            fprintf(stderr, "CreateNewBiPartHyperGraph(): Not enough memory!\n");
            return FALSE;
        }

        if (StoreMat) {
            pHG->Vtx2MatIndex = (long *) malloc(pHG->NrVertices * sizeof(long));
    
            if (pHG->Vtx2MatIndex == NULL) {
                fprintf(stderr, "CreateNewBiPartHyperGraph(): Not enough memory!\n");
                return FALSE;
            }
        }
    }

    /* Allocate memory for the nets */
    if (pHG->NrNets > 0) {
        pHG->N = (struct net *) malloc(pHG->NrNets * sizeof(struct net));

        if (pHG->N == NULL) {
            fprintf(stderr, "CreateNewBiPartHyperGraph(): Not enough memory!\n");
            return FALSE;
        }

        if (StoreMat) {
            pHG->Net2MatIndex = (long *) malloc(pHG->NrNets * sizeof(long));
  
            if (pHG->Net2MatIndex == NULL) {
                fprintf(stderr, "CreateNewBiPartHyperGraph(): Not enough memory!\n");
                return FALSE;
            }
        }
    }

    /* Allocate memory for the pins */

    if (pHG->NrPins > 0) {
        pHG->VtxAdjncy = (long *) malloc(pHG->NrPins * sizeof(long));
        pHG->NetAdjncy = (long *) malloc(pHG->NrPins * sizeof(long));
        
        if (pHG->VtxAdjncy == NULL || pHG->NetAdjncy == NULL) {
            fprintf(stderr, "CreateNewBiPartHyperGraph(): Not enough memory!\n");
            return FALSE;
        }

        if (StoreMat) {
            if (MType != 'P') {
                pHG->MatReValue = (double *) malloc(pHG->NrPins * sizeof(double));
  
                if (pHG->MatReValue == NULL) {
                    fprintf(stderr, "CreateNewBiPartHyperGraph(): Not enough memory!\n");
                    return FALSE;
                }
            }
            if (MType == 'C') {
                pHG->MatImValue = (double *) malloc(pHG->NrPins * sizeof(double));
  
                if (pHG->MatImValue == NULL) {
                    fprintf(stderr, "CreateNewBiPartHyperGraph(): Not enough memory!\n");
                    return FALSE;
                }
            }
        }
    }

    /* Clear all values for the vertices */
    for (t = 0; t < pHG->NrVertices; t++) {
        pHG->V[t].vtxwgt = 0;
        pHG->V[t].iStart = 0;
        pHG->V[t].iEnd = 0;
        pHG->V[t].partition = 0;
        pHG->V[t].Free = FALSE;
        pHG->V[t].GBentry = NULL;
        pHG->OptPartVtx[t] = 0;
        if (StoreMat)
            pHG->Vtx2MatIndex[t] = 0;
    }

    /* Clear all values for the nets */
    for (t = 0; t < pHG->NrNets; t++) {
        pHG->N[t].netwgt = 0;
        pHG->N[t].iStartP0 = 0;
        pHG->N[t].iStartP1 = 0;
        pHG->N[t].iEnd = 0;
        pHG->N[t].dir = ROW;
        pHG->N[t].Free = FALSE;
        
        if (StoreMat)
            pHG->Net2MatIndex[t] = 0;
    }

    /* Clear all values for the pins */
    for (t = 0; t < pHG->NrPins; t++) {
        pHG->VtxAdjncy[t] = 0;
        pHG->NetAdjncy[t] = 0;
        if (StoreMat) {
            if (MType != 'P')
                pHG->MatReValue[t] = 0.0;
            if (MType == 'C')
                pHG->MatImValue[t] = 0.0;
        }
    }

    /* Initialize the gain buckets */
    pHG->GBVtx[0].NrBuckets = pHG->GBVtx[1].NrBuckets = 0;
    pHG->GBVtx[0].Root = pHG->GBVtx[1].Root = NULL;

    /* Clear weights */
    pHG->WeightP[0] = pHG->WeightP[1] = 0;

    /* Clear comm, log, dir */
    pHG->CurComm = pHG->MinComm = pHG->OptComm = 0;
    pHG->CurVtxLog = pHG->MinVtxLog = 0;
    pHG->SplitDir = ROW;

    return TRUE;

} /* end CreateNewBiPartHyperGraph */


int DeleteBiPartHyperGraph(struct biparthypergraph *pHG) {

    /* This function deletes a bipartitioned hypergraph
       by freeing all its memory arrays and resetting all
       its variables to 0.

       The order of deallocation is the reverse of the order
       in which memory was allocated, to improve the chance
       of freeing large contiguous memory chunks.
    */
    if (!pHG) {
        fprintf(stderr, "DeleteBiPartHyperGraph(): Null argument provided!\n");
        return FALSE;
    }

    ClearGainBucket(&(pHG->GBVtx[0]));
    ClearGainBucket(&(pHG->GBVtx[1]));
  
    if (pHG->MatImValue != NULL)
        free(pHG->MatImValue);
    if (pHG->MatReValue != NULL)
        free(pHG->MatReValue);
    if (pHG->NetAdjncy != NULL)
        free(pHG->NetAdjncy);
    if (pHG->VtxAdjncy != NULL)
        free(pHG->VtxAdjncy);
    if (pHG->Net2MatIndex != NULL)
        free(pHG->Net2MatIndex);
    if (pHG->N != NULL)
        free(pHG->N);
    if (pHG->Vtx2MatIndex != NULL)
        free(pHG->Vtx2MatIndex);
    if (pHG->OptPartVtx != NULL)
        free(pHG->OptPartVtx);
    if (pHG->VtxMoveLog != NULL)
        free(pHG->VtxMoveLog);
    if (pHG->V != NULL)
        free(pHG->V);      
  
    pHG->NrVertices = pHG->NrNets = pHG->NrPins = 0;
  
    /* Clear weights */
    pHG->WeightP[0] = pHG->WeightP[1] = 0;
  
    /* Clear comm, log, dir */
    pHG->CurComm = pHG->MinComm = pHG->OptComm = 0;
    pHG->CurVtxLog = pHG->MinVtxLog = 0;
    pHG->SplitDir = ROW;

    return TRUE;
} /* end DeleteBiPartHyperGraph */


int SparseMatrix2BiPartHyperGraph(const struct sparsematrix *pM, 
                                                      int dir, const struct opts *pOptions, struct biparthypergraph *pHG) {

    /* This function translates a sparse matrix A into a hypergraph,
       which is returned as output.

       If the split direction dir = ROW, matrix rows become vertices
       of the hypergraph, and columns become nets.
       If the split direction dir = COL, matrix columns become vertices,
       and rows become nets.
       If the split direction dir = FINEGRAIN, matrix nonzeros become vertices,
       and rows as well as columns become nets.
 
       If the split strategy is local ratio, the direction given by
       the input parameter is ignored, and the direction is determined
       based on the aspect ratio between nonempty rows and nonempty columns.
       Only nonempty rows and columns are included in the hypergraph.

       The numerical values of the matrix  nonzeros are stored as well
       (with the pins of the hypergraph):
           - a real part, if the matrix is not pattern-only, i.e.
             pM->MMTypeCode[2] <> 'P';
           - an imaginary part if the matrix is complex, i.e.
             pM->MMTypeCode[2] = 'C'.

    */
    register long i, j, t;

    long NrNeRows, NrNeCols, *nnzR, *nnzC, *Row2HGIndex, *Col2HGIndex,
         pin_index, vtxwgt, v=-1, n=-1, nC=-1, Acstart, Acend, Arstart,
         Arend;
    long *mgNeR = NULL, *mgNeC = NULL, nDummies = 0, nVertices = 0, nNets = 0;
    int *Diag;

    if (!pM || !pHG || !pOptions) {
        fprintf(stderr, "SparseMatrix2BiPartHyperGraph(): Null arguments provided!\n");
        return FALSE;
    }

    /* Initialise relevant array pointers to NULL,
       to prevent compilation warnings. */
    pHG->V = NULL;
    pHG->Vtx2MatIndex = NULL;
    pHG->N = NULL;
    pHG->Net2MatIndex = NULL;
    pHG->VtxAdjncy = NULL;
    pHG->NetAdjncy = NULL;
    pHG->MatReValue = NULL;
    pHG->MatImValue = NULL;

    /* Allocate memory */
    nnzR = (long *) malloc(pM->m * sizeof(long));
    nnzC = (long *) malloc(pM->n * sizeof(long));

    if (dir!=MEDIUMGRAIN){
        Row2HGIndex = (long *) malloc(pM->m * sizeof(long));
        Col2HGIndex = (long *) malloc(pM->n * sizeof(long));
    }else{
        Col2HGIndex = (long *) malloc((pM->m+pM->n) * sizeof(long));	
        Row2HGIndex = (long *) malloc((pM->m+pM->n) * sizeof(long));	
        for (i = 0; i < pM->m+pM->n; i++) {
            Row2HGIndex[i] = -1;
            Col2HGIndex[i] = -1;
        }
    }

    Diag =  (int *) malloc(pM->m * sizeof(int));

    if (nnzR == NULL || nnzC == NULL || Row2HGIndex == NULL ||
         Col2HGIndex == NULL || Diag == NULL) {
        fprintf(stderr, "SparseMatrix2BiPartHyperGraph(): Not enough memory");
        return FALSE;
    }

    /* Initialise */
    for (i = 0; i < pM->m; i++) {
        nnzR[i] = 0;
        Row2HGIndex[i] = -1;
        Diag[i] = FALSE;
    }

    for (j = 0; j < pM->n; j++) {
        nnzC[j] = 0;
        Col2HGIndex[j] = -1;
    }

    /* Count the number of elements in each row/column of the matrix pM-> 
       Check whether each row contains a dummy or a regular diagonal element */

    NrNeRows = 0;
    NrNeCols = 0;

    for (t = 0; t < pM->NrNzElts; t++) {
        i = pM->i[t];
        j = pM->j[t];
        
        if (nnzR[i] == 0)
            NrNeRows++; /* increment number of nonempty rows */
        
        nnzR[i]++;      /* increment number of nonzeros of row i */

        if (nnzC[j] == 0)
            NrNeCols++;
        
        nnzC[j]++;

        if (i == j &&  pM->m == pM->n &&  
             (pM->MMTypeCode[3]=='S' || pM->MMTypeCode[3]=='K' || pM->MMTypeCode[3]=='H')) 
            Diag[i] = TRUE; /* row i has a nonzero diagonal element */
    }

    /* We need the overlap between the nonempty rows and nonempty columns
       to count the amount of nets required in symmetric finegrain. This is
       easiest accomplished by checking for all i if nnzC[i]>0 while nnzR 
       is zero, and incrementing the number of non-empty rows if so.

       Note that in symmetric finegrain, most row-net variables are abused
       to represent the symmetric nets instead. */
    if (pOptions->SplitStrategy == SFineGrain)
        for (i=0; i<pM->m; i++)
            if (nnzR[i]==0 && nnzC[i]>0)
                NrNeRows++; /* Take this symmetric net into account */

    /* If the SplitStrategy is mediumgrain, we need to calculate the number
       of nonempty columns and which dummies have to be added.  */
    if (dir == MEDIUMGRAIN) {
        nDummies=0;
        mgNeR = (long *)malloc((pM->m+pM->n)*sizeof(long));
        mgNeC = (long *)malloc((pM->m+pM->n)*sizeof(long));
        nVertices=0;
        nNets=0;
        for(i=0;i<pM->m+pM->n;i++){
            mgNeR[i]=0;
            mgNeC[i]=0;
        }
        if(pM->mgDir==0){
            Acstart=0;
            Acend = pM->mgMid;
            Arstart = pM->mgMid;
            Arend = pM->NrNzElts;
        }else{
            Arstart=0;
            Arend = pM->mgMid;
            Acstart = pM->mgMid;
            Acend = pM->NrNzElts;
        }
        for(t=Acstart;t<Acend;t++){ /* in A^c */
            i = pM->i[t];
            j = pM->j[t];
            if(mgNeC[j]==0) nVertices++;
            mgNeC[j]++;
            if(mgNeR[i+pM->n]==0) nNets++;
            mgNeR[i+pM->n]++;
        }
        for(t=Arstart;t<Arend;t++){ /* in A^r */
            i = pM->i[t];
            j = pM->j[t];
            if(mgNeC[i+pM->n]==0) nVertices++;
            mgNeC[i+pM->n]++;
            if(mgNeR[j]==0) nNets++;
            mgNeR[j]++;
        }
        for(i=0;i<pM->m+pM->n;i++){
            if(mgNeR[i]!=0&&mgNeC[i]!=0) nDummies++;
        }
    }

    /* If the SplitStrategy is localratio, set the direction,
       overriding the direction given on input  */
    if (pOptions->SplitStrategy == LocalRatio) {
        if (NrNeRows > NrNeCols)   
            dir = ROW;
        else if (NrNeCols > NrNeRows)   
            dir = COL;
        else if (Random1(0,1) == 0) /* random tie-breaking */
            dir = ROW;
        else
            dir = COL;
#ifdef INFO2
        if (dir == ROW)
            printf("**** Split in row direction \n");
        else
            printf("**** Split in col direction \n");
#endif
    }

    /* Create new bipartitioned hypergraph */
    if (dir == ROW) {
        if (!CreateNewBiPartHyperGraph(NrNeRows, NrNeCols, pM->NrNzElts,
                                       TRUE, pM->MMTypeCode[2], pHG)) {
            fprintf(stderr, "SparseMatrix2BiPartHyperGraph(): Could not create hypergraph!\n");
            return FALSE;
        }
    }
    else if (dir == COL) {
        if (!CreateNewBiPartHyperGraph(NrNeCols, NrNeRows, pM->NrNzElts,
                                       TRUE, pM->MMTypeCode[2], pHG)) {
            fprintf(stderr, "SparseMatrix2BiPartHyperGraph(): Could not create hypergraph!\n");
            return FALSE;
        }
    }
    else if (dir == FINEGRAIN) {
        if (!CreateNewBiPartHyperGraph(pM->NrNzElts, NrNeRows + NrNeCols, 2*pM->NrNzElts,
                                       TRUE, pM->MMTypeCode[2], pHG)) {
            fprintf(stderr, "SparseMatrix2BiPartHyperGraph(): Could not create hypergraph!\n");
            return FALSE;
        }
    } else if (dir == SFINEGRAIN) { 
        if (!CreateNewBiPartHyperGraph(pM->NrNzElts, NrNeRows, 2*pM->NrNzElts,
                                       TRUE, pM->MMTypeCode[2], pHG)) {
            fprintf(stderr, "SparseMatrix2BiPartHyperGraph(): Could not create hypergraph!\n");
            return FALSE;
        }
    } else if (dir == MEDIUMGRAIN) { 
        if (!CreateNewBiPartHyperGraph(nVertices, nNets, pM->NrNzElts+nDummies,
                                       TRUE, pM->MMTypeCode[2], pHG)) {
            fprintf(stderr, "SparseMatrix2BiPartHyperGraph(): Could not create hypergraph!\n");
            return FALSE;
        }
    } else {
        fprintf(stderr,"SparseMatrix2BiPartHypergraph(): Undefined splitting strategy\n");
        return FALSE;
    }

    /* Store the split direction */
    pHG->SplitDir = dir;

    if (pHG->SplitDir == ROW) {
        /* Insert the nonempty rows of the matrix as vertices in the hypergraph and
           register the corresponding row number in Vtx2MatIndex */
        t = 0;
        pin_index = 0;
        
        for (i = 0; i < pM->m; i++) {
            if (nnzR[i] > 0) {
                pHG->V[t].iStart = pHG->V[t].iEnd = pin_index;
                vtxwgt = nnzR[i];

                if (pM->dummy !=  NULL && pM->dummy[i])
                    vtxwgt--;

                if (pM->m == pM->n && 
                    (pM->MMTypeCode[3]=='S' || pM->MMTypeCode[3]=='K' ||
                     pM->MMTypeCode[3]=='H') &&
                     pOptions->SymmetricMatrix_UseSingleEntry == SingleEntYes) {
                    /* matrix is symmetric and each off-diagonal nonzero
                       has weight 2 */
                    if (Diag[i])
                        vtxwgt += (nnzR[i] -1);
                    else
                        vtxwgt += nnzR[i]; 
                }

                pHG->V[t].vtxwgt = vtxwgt;
                pin_index += nnzR[i];
                pHG->Vtx2MatIndex[t] = i;
                Row2HGIndex[i] = t; /* temporarily register inverse mapping */ 
                t++;
            }
        }

        /* Insert the nonempty columns of the matrix as nets in the hypergraph and
           register the corresponding column number in Net2MatIndex */
        t = 0;
        pin_index = 0;
        
        for (j = 0; j < pM->n; j++) {
            if (nnzC[j] > 0) {
                if (pOptions->Metric == MetricCut) pHG->N[t].netwgt = (pM->ColLambda[j] > 1 ? 0 : 1);
                else if (pOptions->Metric == MetricLambdaLambdaMinusOne) pHG->N[t].netwgt = 2*pM->ColLambda[j];
                else pHG->N[t].netwgt = 1; /* uniform netweights for the time being */
                
                pHG->N[t].iStartP0 = pHG->N[t].iStartP1 = pHG->N[t].iEnd = pin_index;
                pHG->N[t].dir = COL;
                if (pOptions->Metric == MetricCut) pHG->N[t].Free = (pM->ColLambda[j] > 1 ? TRUE : FALSE);
                pin_index += nnzC[j];
                pHG->Net2MatIndex[t] = j;
                Col2HGIndex[j] = t;
                t++;
            }
        }
    } else if (pHG->SplitDir == COL) {
        /* Insert the nonempty columns of the matrix as vertices in the hypergraph and
           register the corresponding column number in Vtx2MatIndex */
        t = 0;
        pin_index = 0;
        
        for (j = 0; j < pM->n; j++) {
            if (nnzC[j] > 0) {
                pHG->V[t].iStart = pHG->V[t].iEnd = pin_index;
                
                if (pM->MMTypeCode[0] == 'W' && pM->NrColWeights > 0) {
                    pHG->V[t].vtxwgt = pM->ColWeights[j];
                }
                else {
                    vtxwgt = nnzC[j];

                    if (pM->dummy != NULL && pM->dummy[j])
                        vtxwgt--;

                    if (pM->m == pM->n && 
                        (pM->MMTypeCode[3]=='S' || pM->MMTypeCode[3]=='K' ||
                         pM->MMTypeCode[3]=='H') &&
                         pOptions->SymmetricMatrix_UseSingleEntry == SingleEntYes){
                        if (Diag[j])
                            vtxwgt += (nnzC[j] -1);
                        else
                            vtxwgt += nnzC[j];
                    }
                    pHG->V[t].vtxwgt = vtxwgt;
                }  

                pin_index += nnzC[j];
                pHG->Vtx2MatIndex[t] = j;
                Col2HGIndex[j] = t; /* temporarily register inverse mapping */ 
                t++;
            }
        }

        /* Insert the nonempty rows of the matrix as nets in the hypergraph and
           register the corresponding row number in Net2MatIndex */
        t = 0;
        pin_index = 0;
        
        for (i = 0; i < pM->m; i++) {
            if (nnzR[i] > 0) {
                if (pOptions->Metric == MetricCut) pHG->N[t].netwgt = (pM->RowLambda[i] > 1 ? 0 : 1);
                else if (pOptions->Metric == MetricLambdaLambdaMinusOne) pHG->N[t].netwgt = 2*pM->RowLambda[i];
                else pHG->N[t].netwgt = 1; /* uniform netweights for the time being */
                
                pHG->N[t].iStartP0 = pHG->N[t].iStartP1 = pHG->N[t].iEnd = pin_index;
                pHG->N[t].dir = ROW;
                if (pOptions->Metric == MetricCut) pHG->N[t].Free = (pM->RowLambda[i] > 1 ? TRUE : FALSE);
                pin_index += nnzR[i];
                pHG->Net2MatIndex[t] = i;
                Row2HGIndex[i] = t;
                t++;
            }
        }
    } else if (pHG->SplitDir == MEDIUMGRAIN) {
        /* Insert the nonempty rows and columns of the matrix as vertices in the hypergraph and
           register the corresponding index in Vtx2MatIndex */
        t = 0;
        pin_index = 0;
        
        for (j = 0; j < pM->m+pM->n; j++) {
            if (mgNeC[j] > 0) {
                pHG->V[t].iStart = pHG->V[t].iEnd = pin_index;
                pHG->V[t].vtxwgt = mgNeC[j];
                if (mgNeR[j]!=0){
                    pin_index += mgNeC[j]+1; /* Added dummy */
                }else{
                    pin_index += mgNeC[j];
                }
                pHG->Vtx2MatIndex[t] = j;
                Col2HGIndex[j] = t; /* temporarily register inverse mapping */ 
                t++;
            }
        }

        t = 0;
        pin_index = 0;
        /* Insert the nonempty columns of the matrix as nets in the hypergraph and
           register the corresponding index in Net2MatIndex */       
        for (i = 0; i < pM->n; i++) {
            if (mgNeR[i] > 0) {
                if (pOptions->Metric == MetricCut) pHG->N[t].netwgt = (pM->ColLambda[i] > 1 ? 0 : 1);
                else if (pOptions->Metric == MetricLambdaLambdaMinusOne) pHG->N[t].netwgt = 2*pM->ColLambda[i];
                else pHG->N[t].netwgt = 1; /* uniform netweights for the time being */
                pHG->N[t].iStartP0 = pHG->N[t].iStartP1 = pHG->N[t].iEnd = pin_index;
                pHG->N[t].dir = ROW;
                if (pOptions->Metric == MetricCut) pHG->N[t].Free = (pM->ColLambda[i] > 1 ? TRUE : FALSE);
                if (mgNeC[i]!=0){
                    pin_index += mgNeR[i]+1; /* Added dummy */
                }else{
                    pin_index += mgNeR[i];
                }
                pHG->Net2MatIndex[t] = i;
                Row2HGIndex[i] = t;
                t++;
            }
        }
        /* Insert the nonempty rows of the matrix as nets in the hypergraph and
           register the corresponding index in Net2MatIndex */
        for (i = pM->n; i < pM->m+pM->n; i++) {
            if (mgNeR[i] > 0) {
                if (pOptions->Metric == MetricCut) pHG->N[t].netwgt = (pM->RowLambda[i-pM->n] > 1 ? 0 : 1);
                else if (pOptions->Metric == MetricLambdaLambdaMinusOne) pHG->N[t].netwgt = 2*pM->RowLambda[i-pM->n];
                else pHG->N[t].netwgt = 1; /* uniform netweights for the time being */
                pHG->N[t].iStartP0 = pHG->N[t].iStartP1 = pHG->N[t].iEnd = pin_index;
                pHG->N[t].dir = ROW;
                if (pOptions->Metric == MetricCut) pHG->N[t].Free = (pM->RowLambda[i-pM->n] > 1 ? TRUE : FALSE);
                if (mgNeC[i]!=0){
                    pin_index += mgNeR[i]+1; /* Added dummy */
                }else{
                    pin_index += mgNeR[i];
                }
                pHG->Net2MatIndex[t] = i;
                Row2HGIndex[i] = t;
                t++;
            }
        }
    } else if (pHG->SplitDir == FINEGRAIN) {
        /* Insert the nonempty rows of the matrix as nets in the hypergraph and
           register the corresponding row number in Net2MatIndex */
        t = 0;
        pin_index = 0;
        
        for (i = 0; i < pM->m; i++) {
            if (nnzR[i] > 0) {
                if (pOptions->Metric == MetricCut) pHG->N[t].netwgt = (pM->RowLambda[i] > 1 ? 0 : 1);
                else if (pOptions->Metric == MetricLambdaLambdaMinusOne) pHG->N[t].netwgt = 2*pM->RowLambda[i];
                else pHG->N[t].netwgt = 1; /* uniform netweights for the time being */
                
                pHG->N[t].iStartP0 = pHG->N[t].iStartP1 = pHG->N[t].iEnd = pin_index;
                pHG->N[t].dir = ROW;
                if (pOptions->Metric == MetricCut) pHG->N[t].Free = (pM->RowLambda[i] > 1 ? TRUE : FALSE);
                pin_index += nnzR[i];
                pHG->Net2MatIndex[t] = i;
                Row2HGIndex[i] = t;
                t++;
            }
        }
        /* Insert the nonempty columns of the matrix as nets in the hypergraph and
           register the corresponding column number in Net2MatIndex */
        for (j = 0; j < pM->n; j++){
            if (nnzC[j] > 0) {
                if (pOptions->Metric == MetricCut) pHG->N[t].netwgt = (pM->ColLambda[j] > 1 ? 0 : 1);
                else if (pOptions->Metric == MetricLambdaLambdaMinusOne) pHG->N[t].netwgt = 2*pM->ColLambda[j];
                else pHG->N[t].netwgt = 1; /* uniform netweights for the time being */
                
                pHG->N[t].iStartP0 = pHG->N[t].iStartP1 = pHG->N[t].iEnd = pin_index;
                pHG->N[t].dir = COL;
                if (pOptions->Metric == MetricCut) pHG->N[t].Free = (pM->ColLambda[j] > 1 ? TRUE : FALSE);
                pin_index += nnzC[j];
                pHG->Net2MatIndex[t] = j;
                Col2HGIndex[j] = t;
                t++;
            }
        }

        /* Insert the nonzeros of the matrix as vertices in the hypergraph */
        pin_index = 0;
        
        for (t = 0; t < pM->NrNzElts; t++) {
            pHG->V[t].iStart = pHG->V[t].iEnd = pin_index;

            i = pM->i[t];
            j = pM->j[t];
            vtxwgt = 1; /* the default value */

            if (i == j && pM->dummy != NULL &&  pM->dummy[i])
                vtxwgt = 0;
            if (i != j && pM->m == pM->n &&
                 (pM->MMTypeCode[3]=='S' || pM->MMTypeCode[3]=='K' ||
                  pM->MMTypeCode[3]=='H') &&
                 pOptions->SymmetricMatrix_UseSingleEntry == SingleEntYes)
                    /* matrix is symmetric and each off-diagonal nonzero has weight 2 */
                    vtxwgt = 2;

            pHG->V[t].vtxwgt = vtxwgt;
            pin_index += 2;
        }
    } else if (pHG->SplitDir == SFINEGRAIN) {
        if (pM->MMTypeCode[3]!='S' && pM->MMTypeCode[3]!='K' &&
                  pM->MMTypeCode[3]!='H') {
            fprintf(stderr, "SparseMatrix2BiPartHypergraph(): Symmetric finegrain used on a non-symmetric matrix (did you enable the SymmetricMatrix_UseSingleEntry option?)\n");
            return FALSE;
        }

        /* Count number of nonzeros (pins) in each net */
        for (i = 0; i < pM->n; i++)
            Col2HGIndex[i] = 0;
        for (t = 0; t < pM->NrNzElts; t++) {
            if( pM->i[t] == pM->j[t] )
                Col2HGIndex[pM->i[t]]++; /* One diagonal element */
            else {
                Col2HGIndex[pM->i[t]]++; /* Two off-diagonal elements */
                Col2HGIndex[pM->j[t]]++;
            }
        }

        /* Insert the nonempty rows of the matrix as nets in the hypergraph and
           register the corresponding row number in Net2MatIndex */
        t = 0;
        pin_index = 0;
        
        /* Store symmetric nets */
        for (i = 0; i < pM->m; i++) {
            if (nnzR[i] > 0 || nnzC[i]>0) {
                if (pOptions->Metric == MetricCut) pHG->N[t].netwgt = (pM->RowLambda[i] > 1 ? 0 : 1);
                else if (pOptions->Metric == MetricLambdaLambdaMinusOne) pHG->N[t].netwgt = 2*pM->RowLambda[i];
                else pHG->N[t].netwgt = 2; /* uniform, 2 since each row net also represents a column */
                
                pHG->N[t].iStartP0 = pHG->N[t].iStartP1 = pHG->N[t].iEnd = pin_index;
                pHG->N[t].dir = ROW;
                if (pOptions->Metric == MetricCut) pHG->N[t].Free = (pM->RowLambda[i] > 1 ? TRUE : FALSE);
                pin_index += Col2HGIndex[i]; /* There will be that many pins in this net (see above) */
                pHG->Net2MatIndex[t] = i;
                Row2HGIndex[i] = t;
                t++;
            }
        }

        /* Insert the nonzeros of the matrix as vertices in the hypergraph */
        pin_index = 0;
        
        for (t = 0; t < pM->NrNzElts; t++) {
            pHG->V[t].iStart = pHG->V[t].iEnd = pin_index;

            i = pM->i[t];
            j = pM->j[t];
            vtxwgt = 1; /* the default value */

            if (i == j && pM->dummy != NULL &&  pM->dummy[i])
                vtxwgt = 0;
            if (i != j && pM->m == pM->n &&
                 (pM->MMTypeCode[3]=='S' || pM->MMTypeCode[3]=='K' ||
                  pM->MMTypeCode[3]=='H') &&
                 pOptions->SymmetricMatrix_UseSingleEntry == SingleEntYes)
                    /* matrix is symmetric and each off-diagonal nonzero has weight 2 */
                    vtxwgt = 2;

            pHG->V[t].vtxwgt = vtxwgt;
            /* Also remember column location of nonzero, for back translation */
            pHG->Vtx2MatIndex[t] = j;
            pin_index += 2;
        }
    }


    /* Insert the nonzero elements of the matrix as pins in the 
       hypergraph and fill the adjacency lists */
    for (t = 0; t < pM->NrNzElts; t++) {
        if (pHG->SplitDir == ROW) {
            v =  Row2HGIndex[pM->i[t]];
            n =  Col2HGIndex[pM->j[t]];
        } else if (pHG->SplitDir == COL) {
            v =  Col2HGIndex[pM->j[t]];
            n =  Row2HGIndex[pM->i[t]];
        } else if (pHG->SplitDir == MEDIUMGRAIN) {
            if(pM->mgDir==0){
                if(t<pM->mgMid){
                    v =  Col2HGIndex[pM->j[t]];
                    n =  Row2HGIndex[pM->i[t]+pM->n];
                }else{
                    v =  Col2HGIndex[pM->i[t]+pM->n];
                    n =  Row2HGIndex[pM->j[t]];
                }
            }else{
                if(t<pM->mgMid){
                    v =  Col2HGIndex[pM->i[t]+pM->n];
                    n =  Row2HGIndex[pM->j[t]];
                }else{
                    v =  Col2HGIndex[pM->j[t]];
                    n =  Row2HGIndex[pM->i[t]+pM->n];
                }
            }
        } else if (pHG->SplitDir == FINEGRAIN) {
            v = t; /* vertex = nonzero */
            n = Row2HGIndex[pM->i[t]];  /* row-net */
            nC = Col2HGIndex[pM->j[t]]; /* column-net */
        } else if (pHG->SplitDir == SFINEGRAIN) {
            v = t; /* vertex = nonzero */
            n = Row2HGIndex[pM->i[t]]; /* symm. finegrain */
            nC= Row2HGIndex[pM->j[t]]; /* symm. finegrain */
        }
        if ( n == -1 || ((pHG->SplitDir == FINEGRAIN || pHG->SplitDir == SFINEGRAIN) && nC == -1) ) {
            fprintf(stderr, "SparseMatrix2BiPartHypergraph(): Invalid net number (net corresponding to current vertex is not in hypergraph)\n");
            return FALSE;
        }

        /* Store a numerical value with the vertex */
        if (pHG->SplitDir == ROW || pHG->SplitDir == COL || pHG->SplitDir == MEDIUMGRAIN) {
            /* in the order of the adjacency lists */
            if (pM->MMTypeCode[2] != 'P')
                pHG->MatReValue[pHG->V[v].iEnd] = pM->ReValue[t]; 
            if (pM->MMTypeCode[2] == 'C')
                pHG->MatImValue[pHG->V[v].iEnd] = pM->ImValue[t]; 
         } else if (pHG->SplitDir == FINEGRAIN || pHG->SplitDir == SFINEGRAIN) {
            /* in the order of the vertices */
            if (pM->MMTypeCode[2] != 'P')
                pHG->MatReValue[v] = pM->ReValue[t]; 
            if (pM->MMTypeCode[2] == 'C')
                pHG->MatImValue[v] = pM->ImValue[t]; 
         }
        
        /* Store the net with the vertex.
           This is the row-net in case of the finegrain direction */
        pHG->VtxAdjncy[pHG->V[v].iEnd++] = n;

        /* Store the vertex with the net */
        pHG->NetAdjncy[pHG->N[n].iEnd++] = v;
        if ((n+1) < pHG->NrNets && pHG->N[n].iEnd > pHG->N[n+1].iStartP0) {
            fprintf(stderr, "SparseMatrix2BiPartHypergraph(): not enough space reserved (%ld) for net adjacency list! (row direction)\n", (pHG->N[n].iEnd-pHG->N[n].iStartP0));
            return FALSE;
        }

        if (pHG->SplitDir == FINEGRAIN || (pHG->SplitDir == SFINEGRAIN && n != nC)) {
            /* Store the column-net with the vertex */
            pHG->VtxAdjncy[pHG->V[v].iEnd++] = nC;
            /* Store the vertex with the column-net */
            pHG->NetAdjncy[pHG->N[nC].iEnd++] = v;
            if ((nC+1) < pHG->NrNets && pHG->N[nC].iEnd > pHG->N[nC+1].iStartP0) {
                fprintf(stderr, "SparseMatrix2BiPartHypergraph(): not enough space reserved (%ld) for net adjacency list! (column direction)\n", (pHG->N[n].iEnd-pHG->N[n].iStartP0));
                return FALSE;
            }
        }
    }

    if(dir==MEDIUMGRAIN){ /* Add dummies */
        for(i=0;i<pM->m+pM->n;i++){
            if(mgNeC[i]!=0 && mgNeR[i]!=0){
                v = Col2HGIndex[i];
                n = Row2HGIndex[i];
                if (pM->MMTypeCode[2] != 'P')
                    pHG->MatReValue[pHG->V[v].iEnd] = 0; 
                if (pM->MMTypeCode[2] == 'C')
                    pHG->MatImValue[pHG->V[v].iEnd] = 0; 
                /* Store the net with the vertex.
                    This is the row-net in case of the finegrain direction */
                pHG->VtxAdjncy[pHG->V[v].iEnd++] = n;

                /* Store the vertex with the net */
                pHG->NetAdjncy[pHG->N[n].iEnd++] = v;
                if ((n+1) < pHG->NrNets && pHG->N[n].iEnd > pHG->N[n+1].iStartP0) {
                    fprintf(stderr, "SparseMatrix2BiPartHypergraph(): not enough space reserved (%ld) for net adjacency list! (row direction)\n", (pHG->N[n].iEnd-pHG->N[n].iStartP0));
                    return FALSE;
                }
            }
        }
    }

    /* Deallocate memory */
    free(Diag);
    free(Col2HGIndex);
    free(Row2HGIndex);
    free(nnzC);
    free(nnzR);
    if(dir==MEDIUMGRAIN){
        free(mgNeC);
        free(mgNeR);
    }

    return TRUE;

} /* end SparseMatrix2BiPartHyperGraph */


int BiPartHyperGraph2SparseMatrix(const struct biparthypergraph *pHG, long lo, long hi, long *mid, struct sparsematrix *pM) {

    /* This function translates a hypergraph HG into a subset of the nonzeros
       of the sparse matrix pM The subset is given by a range lo..hi.
       The nonzeros of this range are given new values for the row index i,
       column index j, and the associated numerical values.
       The remainder of pM remains unchanged.

       The nonzeros of the range are split into two parts:
       lo..mid-1 belonging to part 0 (representing one processor)
       mid..hi belonging to part 1 (representing another processor).

       The nonzero arrays are filled starting from both sides simultaneously.

       If the split direction pHG->SplitDir = ROW or COL, nonzeros are obtained
       from the net adjacency lists of the vertices.
       If the split direction pHG->SplitDir = FINEGRAIN, each vertex becomes a 
       matrix nonzero.

       If present, the numerical values of the matrix nonzeros are retrieved
       from the hypergraph as well:
           - a real part, if the matrix is not pattern-only, i.e.
             pM->MMTypeCode[2] <> 'P';
           - an imaginary part if the matrix is complex, i.e.
             pM->MMTypeCode[2] = 'C'.
       This function works both when numerical values are stored and
       when they are not.

    */

    int P;
    long new, n,
         left = lo,   /* lowest index of nonzero in range */ 
         right = hi; /* highest index */
    register long t, tt;
    
    if (!pM || !pHG || !mid) {
        fprintf(stderr, "BiPartHyperGraph2SparseMatrix(): Null arguments provided!\n");
        return FALSE;
    }
  
    for (t = 0; t < pHG->NrVertices; t++) {
        P = pHG->OptPartVtx[t]; /* the part of vertex t in
                                  the optimal partitioning */
  
        if (pHG->SplitDir == ROW || pHG->SplitDir == COL) {
            /* Consider all nets adjacent to vertex t */
            for (tt = pHG->V[t].iStart; tt < pHG->V[t].iEnd; tt++) {
                /* create a new nonzero for the current net tt */
                if (P == 0) 
                    new = left;
                else  /* P = 1 */
                    new = right;
                if (new < lo || new > hi) {
                    fprintf(stderr, "BiPartHyperGraph2SparseMatrix(): nonzero outside range!\n");
                    return FALSE;
                }

                if (pHG->SplitDir == ROW) { 
                    pM->i[new] = pHG->Vtx2MatIndex[t];
                    pM->j[new] = pHG->Net2MatIndex[pHG->VtxAdjncy[tt]];
                } else if (pHG->SplitDir == COL) {
                    pM->j[new] = pHG->Vtx2MatIndex[t];
                    pM->i[new] = pHG->Net2MatIndex[pHG->VtxAdjncy[tt]];
                }

                if (pHG->MatReValue != NULL)
                    pM->ReValue[new] = pHG->MatReValue[tt];
                if (pHG->MatImValue != NULL)
                    pM->ImValue[new] = pHG->MatImValue[tt];

                if (P == 0)
                    left++;
                else
                    right--;
            }
        } else if (pHG->SplitDir == MEDIUMGRAIN) {
            /* Consider all nets adjacent to vertex t */
            for (tt = pHG->V[t].iStart; tt < pHG->V[t].iEnd; tt++) {
                
                if(pHG->Vtx2MatIndex[t]==pHG->Net2MatIndex[pHG->VtxAdjncy[tt]]) continue;
                
                /* create a new nonzero for the current net tt */
                if (P == 0) 
                    new = left;
                else  /* P = 1 */
                    new = right;
                if (new < lo || new > hi) {
                    fprintf(stderr, "BiPartHyperGraph2SparseMatrix(): nonzero outside range!\n");
                    return FALSE;
                }
                
                

                if (pHG->Vtx2MatIndex[t] <pM->n) { 
                    pM->j[new] = pHG->Vtx2MatIndex[t];
                    pM->i[new] = pHG->Net2MatIndex[pHG->VtxAdjncy[tt]]-pM->n;
                } else {
                    pM->i[new] = pHG->Vtx2MatIndex[t]-pM->n;
                    pM->j[new] = pHG->Net2MatIndex[pHG->VtxAdjncy[tt]];
                }
                

                if (pHG->MatReValue != NULL)
                    pM->ReValue[new] = pHG->MatReValue[tt];
                if (pHG->MatImValue != NULL)
                    pM->ImValue[new] = pHG->MatImValue[tt];

                if (P == 0)
                    left++;
                else
                    right--;
            }
        } else if (pHG->SplitDir == FINEGRAIN) {
            if (P == 0)
                new = left;
            else  /* P = 1 */
                new = right;

            /* Create one new nonzero for the two nets adjacent to vertex t */
            if (pHG->V[t].iEnd - pHG->V[t].iStart != 2) {
                fprintf(stderr, "BiPartHyperGraph2SparseMatrix(): corrupted finegrain hypergraph!\n");
                return FALSE;
            }

            for (tt = pHG->V[t].iStart; tt < pHG->V[t].iEnd; tt++) {
                n = pHG->VtxAdjncy[tt]; 
                if (pHG->N[n].dir == ROW)
                    /* net is row-net */
                    pM->i[new] = pHG->Net2MatIndex[n];
                else /* column-net */
                    pM->j[new] = pHG->Net2MatIndex[n];
            }

            if (pHG->MatReValue != NULL)
                pM->ReValue[new] = pHG->MatReValue[t];
            if (pHG->MatImValue != NULL)
                pM->ImValue[new] = pHG->MatImValue[t];

            if (P == 0)
                left++;
            else
                right--;
        } else if (pHG->SplitDir == SFINEGRAIN) {
            if (!(pM->MMTypeCode[3]=='S' || pM->MMTypeCode[3]=='K' ||
                     pM->MMTypeCode[3]=='H')) {
                fprintf(stderr, "BiPartHyperGraph2SparseMatrix(): supplied sparse matrix pointer is not symmetric while in symmetric finegrain mode!\n");
                return FALSE;
            }
            if (P == 0)
                new = left;
            else  /* P = 1 */
                new = right;

            /* Create one or two new nonzeros for the first net adjacent to vertex t */
            if (pHG->V[t].iEnd - pHG->V[t].iStart <= 0 || pHG->V[t].iEnd - pHG->V[t].iStart > 2) {
                fprintf(stderr, "BiPartHyperGraph2SparseMatrix(): corrupted symmetric finegrain hypergraph!\n");
                return FALSE;
            }
            tt = pHG->V[t].iStart;
            n = pHG->VtxAdjncy[tt];
            if (pHG->V[t].iEnd - pHG->V[t].iStart == 2 && pHG->Net2MatIndex[n] <= pHG->Vtx2MatIndex[t]) { /* The next net might be the correct one */
                if (pHG->Net2MatIndex[pHG->VtxAdjncy[tt+1]] >= pHG->Vtx2MatIndex[t]) { /* The next one indeed is lower triangular or diagonal */
                    tt++;
                    n = pHG->VtxAdjncy[tt];
                }
            }
            if (pHG->Net2MatIndex[n] < pHG->Vtx2MatIndex[t]) { /* Still not lower triangular */
                fprintf(stderr, "BiPartHyperGraph2SparseMatrix(): corrupted symmetric finegrain hypergraph\n");
                return FALSE;
            }
            if (pHG->N[n].dir != ROW) {
                fprintf(stderr, "BiPartHyperGraph2SparseMatrix(): corrupted symmetric finegrain hypergraph\n");
                return FALSE;
            }
            /* net is row-net, vertex stores column */
            pM->i[new] = pHG->Net2MatIndex[n];
            pM->j[new] = pHG->Vtx2MatIndex[t];

            if (pHG->MatReValue != NULL)
                pM->ReValue[new] = pHG->MatReValue[t];
            if (pHG->MatImValue != NULL)
                pM->ImValue[new] = pHG->MatImValue[t];

            if (P == 0)
                left++;
            else
                right--;
        }
    }
    
    /* output splitting point */
    *mid = left;
 
    return TRUE;
} /* end BiPartHyperGraph2SparseMatrix */


int SortAdjacencyLists(struct biparthypergraph *pHG) {
    /* This function rearranges the vertex adjacency lists such that the nets are sorted in increasing netsize order.
       Hopefully, this will improve random vertex matching. */
    long *Temp, *NetSizes;
    long t, tt;
    
    Temp = (long *)malloc(pHG->NrNets*sizeof(long));
    NetSizes = (long *)malloc(pHG->NrNets*sizeof(long));
    
    if (Temp == NULL || NetSizes == NULL) {
        fprintf(stderr, "SortAdjacencyLists(): Not enough memory!\n");
        return FALSE;
    }
    
    for (t = 0; t < pHG->NrVertices; ++t) {
        const long NrAdjncy = pHG->V[t].iEnd - pHG->V[t].iStart, Start = pHG->V[t].iStart;
        long *Order;
        
        /* Copy net sizes. */
        for (tt = 0; tt < NrAdjncy; ++tt) {
            const long n = pHG->VtxAdjncy[Start + tt];
            
            NetSizes[tt] = pHG->N[n].iStartP0 - pHG->N[n].iEnd; /* We want to sort from small to large. */
        }
        
        /* Sort them. */
        if ((Order = QSort(NetSizes, NrAdjncy)) == NULL) {
            fprintf(stderr, "SortAdjacencyLists(): Unable to QSort()!\n");
            return FALSE;
        }
        
        /* Permute indices. */
        for (tt = 0; tt < NrAdjncy; ++tt) Temp[tt] = pHG->VtxAdjncy[Start + Order[tt]];
        
        free(Order);
        
        /* Copy permuted adjacency list. */
        memcpy(&pHG->VtxAdjncy[Start], Temp, NrAdjncy*sizeof(long));
    }
    
    /* Free memory. */
    free(NetSizes);
    free(Temp);
            
    return TRUE;
}

