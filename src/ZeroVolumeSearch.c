
#include "ZeroVolumeSearch.h"
#include "SubsetSum.h"

/**
 * Search for a zero-volume split of a (sub)matrix.
 * If such a split is found, the array values in pM->i, pM->j,
 * pM->ReValue and pM->ImValue are rearranged such that the
 * nonzeros in [0,mid-1] form the smaller partition, and the
 * nonzeros in [mid,nnz-1] form the larger partition.
 * 
 * Input:
 * weightlo          : The maximum weight of the smallest partition
 * weighthi          : The maximum weight of the largest partition
 * pOptions          : Options struct
 * 
 * Input/Output:
 * pM                : The matrix (m-by-n)
 * 
 * Output:
 * mid               : The index of the first nonzero in the second partition. This value 
 *                     is only set if return value == TRUE.
 * 
 * Return value: TRUE if a feasible zero-volume split has been found, FALSE otherwise.
 *               If TRUE, the found zero-volume split has already been applied to pM, i.e.
 *               the nonzeros are reordered such that nonzeros [0,mid-1] belong to partition 0,
 *               and nonzeros [mid,nnz-1] belong to partition 1.
 */
int ZeroVolumeSearch(struct sparsematrix *pM, long weightlo, long weighthi, long *mid, const struct opts *pOptions) {
    
    long numComponents = 0;
    long *componentWeights = NULL;
    long *rowToComponent = NULL;
    
    /* Search for connected components */
    if(!DetectConnectedComponents(pM, weighthi, &numComponents, &componentWeights, &rowToComponent, pOptions)) {
        return FALSE;
    }
    
    if(numComponents < 2) {
        free(componentWeights);
        free(rowToComponent);
        if(numComponents == 1) {
            /* As DetectConnectedComponents() returned TRUE, the single component the matrix consists of has
             * a weight lower than weighthi, which implies that we may put all nonzeros in one single partition.
             */
            *mid = 0;
            return TRUE;
        }
        return FALSE;
    }
    
    long *componentIds = NULL;
    long *componentsSelected = NULL;
    
    /* Determine whether the components can be distributed over the partitions */
    int found = SubsetSum(componentWeights, numComponents, weightlo, weighthi, &componentIds, &componentsSelected);
    free(componentWeights);
    
    if(found) {
        
        /* Form a component-to-partition mapping */
        long *componentToPartition = (long*)malloc(numComponents*sizeof(long));
        if( componentToPartition == NULL ) {
            fprintf( stderr, "ZeroVolumeSearch(): Not enough memory.\n" );
            free(rowToComponent);
            free(componentIds);
            free(componentsSelected);
            return FALSE;
        }
        
        long i;
        for(i=0; i<numComponents; ++i) {
            /* In componentsSelected, those entries with value 1 belong to the small partition.
             * For conceptual ease, we place the small partition before the large partition,
             * and make 0 the small partition and 1 the large partition.
             */
            componentToPartition[componentIds[i]] = 1-componentsSelected[i];
        }
        
        /* 'Collapse the walls'; a single iteration of quicksort
         * Swap all zeros to the left, and all ones to the right.
         */
        long lo = 0, hi = pM->NrNzElts-1;
        do {
            
            while(componentToPartition[rowToComponent[pM->i[lo]]] == 0 && lo <= hi) {
                ++lo;
            }
            while(componentToPartition[rowToComponent[pM->i[hi]]] == 1 && lo <= hi) {
                --hi;
            }
            
            if(lo < hi) {
                /* Swap entries */
                SwapLong(pM->i, lo, hi);
                SwapLong(pM->j, lo, hi);
                if(pM->ReValue != NULL)
                    SwapDouble(pM->ReValue, lo, hi);
                if(pM->ImValue != NULL)
                    SwapDouble(pM->ImValue, lo, hi);
                ++lo;
                --hi;
            }
            
        }
        while(lo <= hi);
        /* At this point, hi+1 = lo; and part0 = [0,hi], part1 = [lo,nnz-1] */
        
        /* Register where new partition starts */
        *mid = lo;
        
        free(componentToPartition);
        free(componentIds);
        free(componentsSelected);
    }
    
    free(rowToComponent);
    
    return (found == TRUE);
    
} /* end ZeroVolumeSearch */

/**
 * Detect connected components in a (sub)matrix.
 * When viewing the (sub)matrix as a (bipartite) graph, with the columns and rows as vertices
 * and the nonzeros as edges, this amounts to finding connected components in this graph.
 * This is done using a breadth-first search, starting at a column and then switching between
 * columns and rows in a breadth-first fashion.
 * In case pM contains colWeights, this function assumes these are all positive.
 * 
 * Input:
 * pM                : The matrix (m-by-n)
 * maxWeight         : Maximum weight of a component. If a component is detected of weight
 *                     larger than this value, it will stop and return FALSE.
 * pOptions          : Options struct
 * 
 * Output:
 * numComponents     : The number of components found (Only if return value == TRUE)
 * componentWeights  : The number of nonzeros in each found component (Only if return value == TRUE)
 * rowAssignments    : The components to which each row is assigned (Only if return value == TRUE)
 * 
 * Returned memory allocations:
 *     componentWeights:  O(numComponents) (Only if return value == TRUE)
 *     rowAssignments:    O(m)             (Only if return value == TRUE)
 * 
 * Return value: FALSE if an error occurs, TRUE otherwise
 */

int DetectConnectedComponents(struct sparsematrix *pM, long maxWeight, long *numComponents, long **componentWeights, long **rowAssignments, const struct opts *pOptions) {
    
    struct CRCS CRS, CCS;
    
    if(!SparseMatrixToCRS_CCS(pM, &CCS, &CRS)) {
        return FALSE;
    }
    
    /* Array of numbers of nonzeros in each component */
    long compArrLen = 8;
    long *component_weights = (long*)malloc(sizeof(long)*compArrLen);
    long component = -1; /* Current component */
    
    /* Arrays containing the components the rows/columns are assigned to. */
    long *rowAssign = (long*)malloc(sizeof(long)*pM->m);
    long *colAssign = (long*)malloc(sizeof(long)*pM->n);
    
    /* Arrays containing the rows/columns we still have to process. */
    long *rowStack = (long*)malloc(sizeof(long)*pM->m);
    long *colStack = (long*)malloc(sizeof(long)*pM->n);
    long *rowStackHead = rowStack;
    long *colStackHead = colStack;
    
    /* Check memory */
    if( component_weights == NULL || rowAssign == NULL || colAssign == NULL ||
        rowStack == NULL || colStack == NULL ) {
        if(component_weights != NULL)
            free(component_weights);
        if(rowAssign != NULL)
            free(rowAssign);
        if(colAssign != NULL)
            free(colAssign);
        if(rowStack != NULL)
            free(rowStack);
        if(colStack != NULL)
            free(colStack);
        fprintf( stderr, "DetectConnectedComponents(): Not enough memory.\n" );
        return FALSE;
    }
    
    long i, k, col, row, start, end; /* Indices and iterators */
    
    /* Boolean variables */
    int symmetric = (pM->m == pM->n && (pM->MMTypeCode[3] == 'S' || pM->MMTypeCode[3] == 'K' || pM->MMTypeCode[3] == 'H') && pOptions->SymmetricMatrix_UseSingleEntry == SingleEntYes);
    int dummies = (pM->m == pM->n && pM->NrDummies > 0);
    int colWeights = (pM->MMTypeCode[0] == 'W' && pM->NrColWeights > 0);
    /* If colWeights is false, we count the nonzeros as the weight of a submatrix.
     * Every nonzero is counted exactly twice (once in each direction), hence we divide by 2 in the end.
     * If colWeights is true, the weight of a submatrix is based on the column weights.
     */
    
    long numAssignedRows = 0;
    int broken = FALSE;
    
    /* Initialize */
    for(i=0; i<pM->m; ++i)
        rowAssign[i] = -1;
    for(i=0; i<pM->n; ++i)
        colAssign[i] = -1;
    
    /* Walk through all rows. As soon as all rows have been processed,
     * all nonzeros have been assigned to a component and we are finished.
     */
    for(i=0; i<pM->m; ) {
        /* If we have no new rows to process, we have finished a component.
         * In that case, find the next unassigned row to assign it to the next component.
         */
        if(rowStackHead-rowStack == 0) {
            if(rowAssign[i] != -1) {
                /* Row i is already assigned to a component */
                ++i;
                continue;
            }
            
            /* Found an unassigned row */
            if(component == -1 || component_weights[component] > 0) {
                /* We have an unassigned row unrelated to any previous components, and
                 * the last component is non-empty, hence create a new one.
                 * (Or we have no other components at all; in which case we als should create one.)
                 */
                ++component;
                
                /* Check component array length */
                if(component >= compArrLen) {
                    /* Array is too small; enlarge it */
                    compArrLen *= 2;
                    long *newArr = (long*)realloc(component_weights, sizeof(long)*compArrLen);
                    if (newArr == NULL) {
                        fprintf(stderr, "DetectConnectedComponents(): Not enough memory.\n");
                        broken = TRUE;
                        break;
                    }
                    else {
                        component_weights = newArr;
                    }
                }
                
                component_weights[component] = 0;
                numAssignedRows = 0;
            }
            
            *(rowStackHead++) = i;
            rowAssign[i] = component;
        }
        
        /* Process all found rows, and register new columns */
        while(rowStackHead-rowStack>0) {
            row = *(--rowStackHead);
            start = CRS.starts[row];
            end = CRS.starts[row+1];
            if(start == end) {
                continue;
            }
            if(!colWeights) {
                component_weights[component] += end-start;
            }
            ++numAssignedRows;
            
            if(symmetric) {
                /* When a row is added, the column should also be added */
                if(!colWeights) {
                    component_weights[component] += end-start;
                }
                if(colAssign[row] == -1) {
                    *(colStackHead++) = row;
                    colAssign[row] = component;
                }
            }
            
            /* Walk through all nonzeros in this row */
            for(k=start; k<end; ++k) {
                col = CRS.indices[k];
                if(colAssign[col] == -1) {
                    *(colStackHead++) = col;
                    colAssign[col] = component;
                }
                if(col == row && !colWeights) {
                    if(symmetric) {
                        /* We only count diagonal elements once */
                        --component_weights[component];
                    }
                    if(dummies && pM->dummy[row]) {
                        /* Dummy element */
                        --component_weights[component];
                    }
                }
            }
            
        } /* end while */
        
        /* Process all found columns, and register new rows */
        while(colStackHead-colStack>0) {
            col = *(--colStackHead);
            start = CCS.starts[col];
            end = CCS.starts[col+1];
            if(start == end) {
                continue;
            }
            if(!colWeights) {
                component_weights[component] += end-start;
            }
            else {
                component_weights[component] += pM->ColWeights[col];
            }
            
            if(symmetric) {
                /* When a column is added, the row should also be added */
                if(!colWeights) {
                    component_weights[component] += end-start;
                }
                if(rowAssign[col] == -1) {
                    *(rowStackHead++) = col;
                    rowAssign[col] = component;
                }
            }
            
            /* Walk through all nonzeros in this column */
            for(k=start; k<end; ++k) {
                row = CCS.indices[k];
                if(rowAssign[row] == -1) {
                    *(rowStackHead++) = row;
                    rowAssign[row] = component;
                }
                if(row == col && !colWeights) {
                    if(symmetric) {
                        /* We only count diagonal elements once */
                        --component_weights[component];
                    }
                    if(dummies && pM->dummy[row]) {
                        /* Dummy element */
                        --component_weights[component];
                    }
                }
            }
        } /* end while */
        
        if(component_weights[component] > maxWeight*(colWeights?1:2)) {
            broken = TRUE;
            break;
        }
        
    } /* end for */
    
    free(rowStack);
    free(colStack);
    free(colAssign);
    freeCRCS(&CCS);
    freeCRCS(&CRS);
    
    if(broken == TRUE) {
        free(component_weights);
        free(rowAssign);
        return FALSE;
    }
    
    /* If the last component is empty, discard it. */
    if(component_weights[component] == 0) {
        --component;
        if(numAssignedRows > 0) {
            /* Rows were added, but the weight is 0. As we assume all weights are positive, we may assign the rows to another component. */
            for(i=0; i<pM->m; ++i) {
                if(rowAssign[i] > component) {
                    rowAssign[i] = 0;
                }
            }
        }
    }
    
    if(!colWeights) {
        /* Every nonzero is counted twice, so divide by two */
        for(i=0; i<component+1; ++i) {
            component_weights[i] /= 2;
        }
    }
    
    /* Assign values to return variables */
    *componentWeights = component_weights;
    *numComponents = component + 1;
    *rowAssignments = rowAssign;
    
    return TRUE;
    
} /* end DetectConnectedComponents */
