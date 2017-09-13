#include "testHelper_DisconnectedMatrix.h"

#define min(x,y) ( ( (x)<(y) ) ? (x) : (y) )
#define max(x,y) ( ( (x)>(y) ) ? (x) : (y) )


/**
 * Construct a random (disconnected) matrix.
 * 
 * Input:
 * symmetric         : Whether the matrix should be symmetric
 * dummies           : Whether the matrix should contain dummy diagonal nonzeros
 * colWeights        : Whether the matrix should be weighted with column weights
 * numCompMin        : Minimum number of connected components
 * numCompMax        : Maximum number of connected components
 * 
 * Output:
 * pA                : The generated matrix
 * pNumComponents    : The number of components generated
 * pComponent_m      : Height of each connected component
 * pComponent_n      : Width of each connected component
 * pComponent_weights: Weight of each connected component (effective number of nonzeros or sum of column weights)
 * p_i_to_I          : I = p_i_to_I[cmpnnt][i] is the row index I in the large matrix which corresponds to the row index i in component cmpnnt
 * p_j_to_J          : J = p_j_to_J[cmpnnt][j] is the column index J in the large matrix which corresponds to the column index j in component cmpnnt
 * 
 * Returned memory allocations (to be freed with DestructDisconnectedMatrix()):
 *     pA:
 *         pA->i             O(pA->m)
 *         pA->j             O(pA->n)
 *         pA->Pstart        O(2+1)
 *         pA->dummy         O(pA->m)            (only if dummy == TRUE)
 *         pA->ColWeights    O(pA->n)            (only if colWeights == TRUE)
 *     pComponent_m:         O(numComponents)
 *     pComponent_n:         O(numComponents)    (only if symmetric == FALSE)
 *     pComponent_weights:   O(numComponents)
 *     p_i_to_I:             O(numComponents)
 *     p_i_to_I[0]:          O(pA->m)
 *     p_j_to_J:             O(numComponents)    (only if symmetric == FALSE)
 *     p_j_to_J[0]:          O(pA->n)            (only if symmetric == FALSE)
 */
void ConstructDisconnectedMatrix(struct sparsematrix *pA, int symmetric, int dummies, int colWeights,
                                 long numCompMin, long numCompMax, long *pNumComponents,
                                 long **pComponent_m, long **pComponent_n, long **pComponent_weights,
                                 long ***p_i_to_I, long ***p_j_to_J) {
    
    long i, j, I, J, m, n;
    long cmpnnt;
    
    /* Determine number of components (submatrices), the dimensions of each submatrix, and their numbers of nonzeros */
    long numComponents = Random1(numCompMin, numCompMax);
    long *component_dummies = (long *)malloc(numComponents*sizeof(long));
    long *component_m = (long *)malloc(numComponents*sizeof(long));
    long *component_n;
    if(symmetric || dummies) {
        component_n = component_m;
    }
    else {
        component_n = (long *)malloc(numComponents*sizeof(long));
    }
    long *component_nnz = (long *)malloc(numComponents*sizeof(long));
    
    if (component_dummies == NULL || component_m  == NULL || component_n == NULL || component_nnz  == NULL) {
        printf("Error\n");
        exit(1);
    }
    
    long M=0, N=0, NNZ=0, totalDummies = 0;
    for(cmpnnt=0; cmpnnt<numComponents; ++cmpnnt) {
        
        /* component_nnz[] should uniquely identify a component; i.e. all components should have different nnz.
         * This helps when identifying the components in test functions.
         */
        int present;
        long nonDummyDiagNonZeros[2];
        long full_nnz[2];
        
        do {
            /* The dimensions should at least be above 10, such that generating the
             * matrices won't take too long (or, when below 5, may become impossible).
             * Additionally, we should have enough freedom to ensure each component has
             * a unique number of nonzeros
             */
            component_m[cmpnnt] = Random1(14, 20) + Random1(numComponents, numCompMax);
            if(!symmetric && !dummies) {
                component_n[cmpnnt] = Random1(14, 20) + Random1(numComponents, numCompMax);
            }
            long maxdim = (component_m[cmpnnt]>component_n[cmpnnt]) ? component_m[cmpnnt] : component_n[cmpnnt];
            
            /* Determine number of nonzeros */
            if(symmetric) {
                component_nnz[cmpnnt] = 2*maxdim-1 + Random1(maxdim/2, maxdim);
            }
            else {
                component_nnz[cmpnnt] = 2*maxdim-1 + Random1(maxdim, 2*maxdim);
            }
            
            component_dummies[cmpnnt] = dummies ? Random1(maxdim/4, 3*maxdim/4) : 0;
            component_nnz[cmpnnt] -= component_dummies[cmpnnt];
            
            /* Check that this number of non-dummy nonzeros is unique */
            present = 0;
            for(i=0; i<cmpnnt; ++i) {
                if(symmetric) {
                    nonDummyDiagNonZeros[0] = component_m[cmpnnt] - component_dummies[cmpnnt];
                    full_nnz[0] = (component_nnz[cmpnnt]-nonDummyDiagNonZeros[0])*2 + nonDummyDiagNonZeros[0];
                    
                    nonDummyDiagNonZeros[1] = component_m[i] - component_dummies[i];
                    full_nnz[1] = (component_nnz[i]-nonDummyDiagNonZeros[1])*2 + nonDummyDiagNonZeros[1];
                    
                    if(full_nnz[0] == full_nnz[1]) {
                        present = 1;
                        break;
                    }
                }
                else {
                    
                    if(component_nnz[i] == component_nnz[cmpnnt]) {
                        present = 1;
                        break;
                    }
                }
            }
        }
        while(present);
        
        /* Update totals */
        M += component_m[cmpnnt];
        N += component_n[cmpnnt];
        NNZ += component_nnz[cmpnnt];
        totalDummies += component_dummies[cmpnnt];
    }
    
    /* Allocate space for a mapping from submatrix indices (i,j) to global indices (I,J) */
    long *i_to_I_alloc, *j_to_J_alloc, **i_to_I, **j_to_J;
    i_to_I_alloc = (long *)malloc(M*sizeof(long));
    i_to_I = (long **)malloc(numComponents*sizeof(long *));
    if(symmetric || dummies) {
        j_to_J_alloc = i_to_I_alloc;
        j_to_J = i_to_I;
    }
    else {
        j_to_J_alloc = (long *)malloc(N*sizeof(long));
        j_to_J = (long **)malloc(numComponents*sizeof(long *));
    }
    
    if (i_to_I_alloc == NULL || i_to_I  == NULL || j_to_J_alloc == NULL || j_to_J  == NULL) {
        printf("Error\n");
        exit(1);
    }
    
    i_to_I[0] = i_to_I_alloc;
    if(!symmetric && !dummies)
        j_to_J[0] = j_to_J_alloc;
    for(cmpnnt=1; cmpnnt<numComponents; ++cmpnnt) {
        i_to_I[cmpnnt] = &i_to_I[cmpnnt-1][component_m[cmpnnt-1]];
        if(!symmetric && !dummies)
            j_to_J[cmpnnt] = &j_to_J[cmpnnt-1][component_n[cmpnnt-1]];
    }
    
    /* Create a permutation of the identity */
    for(I=0; I<M; ++I) {
        i_to_I_alloc[I] = I;
    }
    if(!symmetric && !dummies) {
        for(J=0; J<N; ++J) {
            j_to_J_alloc[J] = J;
        }
    }
    
    RandomPermute(i_to_I_alloc, NULL, NULL, NULL, 0, M-1);
    if(!symmetric && !dummies)
        RandomPermute(j_to_J_alloc, NULL, NULL, NULL, 0, N-1);
    
    /* Set up matrix struct */
    MMSparseMatrixInit(pA);
    pA->m = M;
    pA->n = N;
    pA->NrNzElts = NNZ + totalDummies;
    pA->NrProcs = 2; /* maximum number of parts */

    pA->i = (long *)malloc(pA->NrNzElts*sizeof(long));
    pA->j = (long *)malloc(pA->NrNzElts*sizeof(long));
    pA->Pstart = (long *)malloc((pA->NrProcs+1)*sizeof(long));
    pA->RowLambda = (int *)malloc(pA->m*sizeof(int));
    pA->ColLambda = (int *)malloc(pA->n*sizeof(int));
    
    if (pA->i == NULL || pA->j == NULL || pA->Pstart == NULL || pA->RowLambda == NULL || pA->ColLambda == NULL) {
        printf("Error\n");
        exit(1);
    }
    
    if(colWeights)
        pA->MMTypeCode[0]='W'; /* weighted matrix */
    else
        pA->MMTypeCode[0]='M'; /* normal matrix */
    pA->MMTypeCode[1]='C'; /* coordinate scheme */
    pA->MMTypeCode[2]='P'; /* pattern only */
    if(symmetric)
        pA->MMTypeCode[3]='S'; /* symmetric */
    else
        pA->MMTypeCode[3]='G'; /* general, no symmetry */
    pA->NrDummies = 0;
    pA->dummy = NULL;
    pA->Pstart[0] = 0;
    pA->Pstart[1] = pA->NrNzElts;
    pA->Pstart[2] = pA->NrNzElts;
    pA->NrDummies = totalDummies;
    
    if(dummies) {
        pA->dummy = (int *)malloc(pA->m*sizeof(int));
        if (pA->dummy == NULL) {
            printf("Error\n");
            exit(1);
        }
    }
    
    /* Generate each component */
    long maxdim, t = 0, t2, T, nnzFilled, dummiesUsed;
    
    for(cmpnnt=0; cmpnnt<numComponents; ++cmpnnt) {
        m = component_m[cmpnnt];
        n = component_n[cmpnnt];
        maxdim = (m>n) ? m : n;
        nnzFilled = 0;
        dummiesUsed = 0;
        T = t;
        
        /* Make sure we create connected components */
        if(Random1(0,1) == 0) {
            /* Here, we generate a submatrix with nonzero diagonal
             * and one nonzero subdiagonal. This way, each row is
             * connected with the next, and hence all are connected.
             */
            for(i=0; i<maxdim; ++i) {
                I = i_to_I[cmpnnt][i%m];
                J = j_to_J[cmpnnt][i%n];
                pA->i[t] = symmetric ? max(I,J) : I;
                pA->j[t] = symmetric ? min(I,J) : J;
                ++t;
                ++nnzFilled;
                
                if(dummies) {
                    pA->dummy[I] = (dummiesUsed < component_dummies[cmpnnt]) ? 1 : 0;
                    if(pA->dummy[I] == 1) {
                        ++dummiesUsed;
                        --nnzFilled;
                    }
                }
                
                if(i < maxdim-1) {
                    I = i_to_I[cmpnnt][(i+1)%m];
                    J = j_to_J[cmpnnt][i%n];
                    pA->i[t] = symmetric ? max(I,J) : I;
                    pA->j[t] = symmetric ? min(I,J) : J;
                    ++t;
                    ++nnzFilled;
                }
            }
        }
        else {
            /* Here, we generate a submatrix with nonzero diagonal
             * and a dense first column. This way, all rows are
             * connected with the first column, hence all rows and
             * columns are connected.
             */
            for(i=0; i<maxdim; ++i) {
                I = i_to_I[cmpnnt][i%m];
                J = j_to_J[cmpnnt][i%n];
                pA->i[t] = symmetric ? max(I,J) : I;
                pA->j[t] = symmetric ? min(I,J) : J;
                ++t;
                ++nnzFilled;
                
                if(dummies) {
                    pA->dummy[I] = (dummiesUsed < component_dummies[cmpnnt]) ? 1 : 0;
                    if(pA->dummy[I] == 1) {
                        ++dummiesUsed;
                        --nnzFilled;
                    }
                }
                
                if(i > 0) {
                    I = i_to_I[cmpnnt][i%m];
                    J = j_to_J[cmpnnt][0];
                    pA->i[t] = symmetric ? max(I,J) : I;
                    pA->j[t] = symmetric ? min(I,J) : J;
                    ++t;
                    ++nnzFilled;
                }
            }
        }
        
        /* Fill the submatrix until the desired number of nonzeros */
        while(nnzFilled < component_nnz[cmpnnt]) {
            i = Random1(0, m-1);
            j = Random1(0, n-1);
            I = i_to_I[cmpnnt][i%m];
            J = j_to_J[cmpnnt][j%n];
            
            I = symmetric ? max(I,J) : I;
            J = symmetric ? min(I,J) : J;
            
            /* Make sure this element is not already nonzero.
             * This search is expensive in terms of complexity, but
             * as we are just testing, this does not matter very much.
             */
            int present = 0;
            for(t2=T; t2<t; ++t2) {
                if(pA->i[t2] == I && pA->j[t2] == J) {
                    present = 1;
                    break;
                }
            }
            if(present) {
                continue;
            }
            
            pA->i[t] = I;
            pA->j[t] = J;
            ++t;
            ++nnzFilled;
        }
    }
    
    if(symmetric) {
        /* Update nonzero count, as we now only have counted the elements in the lower-left triangle */
        NNZ = 0;
        long nonDummyDiagNonZeros;
        for(cmpnnt=0; cmpnnt<numComponents; ++cmpnnt) {
            nonDummyDiagNonZeros = component_m[cmpnnt] - component_dummies[cmpnnt];
            component_nnz[cmpnnt] = (component_nnz[cmpnnt]-nonDummyDiagNonZeros)*2 + nonDummyDiagNonZeros;
            NNZ += component_nnz[cmpnnt];
        }
    }
    
    /* Return value pointers */
    *pNumComponents = numComponents;
    *pComponent_m = component_m;
    *pComponent_n = component_n;
    *pComponent_weights = component_nnz;
    *p_i_to_I = i_to_I;
    *p_j_to_J = j_to_J;
    
    
    if(colWeights) {
        /* Reassign weights based on columns */
        long *component_weights = *pComponent_weights;
        pA->NrColWeights = N;
        pA->ColWeights = malloc( N * sizeof( long ) );
        if (pA->ColWeights == NULL) {
            printf("Error\n");
            exit(1);
        }
        for(J=0; J<N; ++J) {
            pA->ColWeights[J] = Random1(20, 50);
        }
        
        for(cmpnnt=0; cmpnnt<numComponents; ++cmpnnt) {
            component_weights[cmpnnt] = 0;
            for(j=0; j<component_n[cmpnnt]; ++j) {
                component_weights[cmpnnt] += pA->ColWeights[j_to_J[cmpnnt][j]];
            }
            
            /* Make sure each component has a unique weight, to be able to use it as identification */
            while(TRUE) {
                int found = FALSE;
                for(i=0; i<cmpnnt; ++i) {
                    if(component_weights[cmpnnt] == component_weights[i]) {
                        found = TRUE;
                        ++pA->ColWeights[j_to_J[cmpnnt][0]];
                        ++component_weights[cmpnnt];
                        break;
                    }
                }
                if(!found) {
                    break;
                }
            }
        }
        
    }
    
    free(component_dummies);
    
} /* end ConstructDisconnectedMatrix */


/**
 * Free memory for a matrix constructed with ConstructDisconnectedMatrix().
 * See ConstructDisconnectedMatrix() for details on the parameters.
 */
void DestructDisconnectedMatrix(struct sparsematrix *pA, int symmetric, int dummies, int colWeights,
                                long **pComponent_m, long **pComponent_n, long **pComponent_weights,
                                long ***p_i_to_I, long ***p_j_to_J) {
    free(*pComponent_weights);
    free(*pComponent_m);
    free(*p_i_to_I[0]);
    free(*p_i_to_I);
    if(!symmetric && !dummies) {
        free(*pComponent_n);
        free(*p_j_to_J[0]);
        free(*p_j_to_J);
    }
    
    MMSparseMatrixFreeMemory(pA);
    
} /* end DestructDisconnectedMatrix */
