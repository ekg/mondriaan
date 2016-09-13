#include "Matalloc.h"

long **matallocl(int m, int n) {

    /* This function allocates an m x n matrix of longs */

    int i;
    long *pl, **ppl;

    if (m == 0) {
        ppl = NULL;  
    } else { 
        ppl = (long **)malloc(m*sizeof(long *));
        
        if (ppl == NULL) {
            fprintf(stderr, "matallocl(): Not enough memory!\n");
            return NULL;
        }
        
        if (n == 0) {
            for (i = 0; i < m; i++)
                ppl[i] = NULL;
        } else {  
            pl = (long *)malloc(m*n*sizeof(long));
            
            if (pl == NULL) {
                fprintf(stderr, "matallocl: Not enough memory!\n");
                return NULL;
            }
            
            ppl[0] = pl;
            
            for (i = 1; i < m; i++)
                ppl[i] = ppl[i-1]+n;
        }
    }
    return ppl;

} /* end matallocl */


void matfreel(long **ppl) {

    /* This function frees a matrix of longs */

    if (ppl != NULL) {
        if (ppl[0] != NULL)
            free(ppl[0]);
        free(ppl);
    }

} /* end matfreel */
