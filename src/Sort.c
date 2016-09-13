#include "Sort.h"

long Random1(long lo, long hi) {

    /* This functions returns a randomly chosen index i in the range lo..hi
       if lo <= hi. Otherwise, it returns LONG_MAX. */

    double range;
    long i;

    if (lo > hi)
        return LONG_MAX;
    if (lo == hi) 
        return lo;

    range = (double) (hi - lo + 1);
    i = lo + (long) ((range*rand()) / (RAND_MAX+1.0));
    if (i < lo)
        i = lo;
    if (i > hi)
        i = hi;

    return i;

} /* end Random1 */

int SetRandomSeed(long seed) {

    /* This function sets the seed of the random number generator 
       to the given value if seed >= 0, and to a random value
       in the range 0-999999 determined by the UNIX wallclock
       if seed = -1. */

    unsigned int seed1;

#ifdef UNIX
    struct timeval timeval;
#endif

    if (seed < -1) {
        fprintf(stderr, "SetRandomSeed(): Warning, seed < -1!\n");
    }

    if (seed >= 0)
        /* set the seed to the given value */
        seed1 = (unsigned int) seed;
    else {
        /* set the seed to a random value using the UNIX wallclock */
#ifdef UNIX
        gettimeofday(&timeval, NULL);
        seed1 = (unsigned int) timeval.tv_usec; /* clock time in microseconds */
#endif
#ifndef UNIX
        /* fprintf(stderr, "SetRandomSeed(): random seed unavailable on non-UNIX systems!\n"); */
        seed1 = time(0);
#endif
    }

    srand(seed1);
#ifdef INFO2
    printf("Random number seed = %u\n", seed1);
#endif
    
    return TRUE;
} /* end SetRandomSeed */


void SwapLong(long *x, const long j0, const long j1) {
    /* This function swaps x[j0] and x[j1] in an array x of type long */
    const long tmp = x[j1];

    x[j1] = x[j0];
    x[j0] = tmp;
} /* end SwapLong */


void SwapDouble(double *x, long j0, long j1) { 
 
    /* This function swaps x[j0] and x[j1] in an array x of type double */ 
 
    double tmp; 
 
    tmp = x[j1]; 
    x[j1] = x[j0]; 
    x[j0] = tmp; 
 
} /* end SwapDouble*/ 


long *QSort(long *X, long LengthX) {

    /* This function sorts array X in decreasing order by
       randomised quicksort. It allocates and returns a bookkeeping 
       array I, which stores the original indices of the array items. */

    register long t; 
    long *I;
    
    if (!X) return 0;
   
    I  = (long *) malloc(LengthX * sizeof(long));
    
    if (I == NULL) {
        fprintf(stderr, "QSort(): Not enough memory for sort array!\n");
        
        return 0;
    }
  
    /*## Initialise index array at original values */
    for (t = 0; t < LengthX; t++)
        I[t] = t;
  
    quicksort(X, I, 0, LengthX - 1);
  
    return I;
}   /* end Qsort */


void quicksort(long *item, long *Index, long lo, long hi) {
 
    /* This function recursively sorts the items lo..hi in decreasing order
       and moves the corresponding index values in the same way.  */
 
    register long i, j;
    long i0, i1, i2, tmpidx, split;
 
    if (lo >= hi)
        return;
 
    /* Choose three splitters randomly */
    i0 = Random1(lo, hi);
    i1 = Random1(lo, hi);
    i2 = Random1(lo, hi);
 
    /* Sort i0, i1, i2 in increasing order of item value */
    if (item[i0] < item[i1]) {
        /* swap i0 and i1 */
        tmpidx = i0; i0 = i1; i1 = tmpidx;
    }

    if (item[i1] <  item[i2]) {
        /* swap i1 and i2 */
        tmpidx = i1; i1 = i2; i2 = tmpidx;
    }
    /* item[i2] is minimum value of three */

    if (item[i0] < item[i1]) {
        /* swap i0 and i1 */
        tmpidx = i0; i0 = i1; i1 = tmpidx;
    }
    /* item[i1] is median value of three */
 
    /* Use item[i1] as a splitter */
    split = item[i1];
    i = lo;
    j = hi;
    while (i <= j) {
        /* Loop Invariant (LI):
              for all r: lo <= r < i : item[r] >= split 
              for all r: j < r <= hi : item[r] <= split 
        */
        while (item[i] > split && i < hi)
            i++;
        while (item[j] < split && j > lo)
            j--;
        /* LI and lo <= i,j <= hi */
        if (i < j) {
            /* item[i] <= split <= item[j], so swap item i and j */
            SwapLong(item, i, j);
            SwapLong(Index, i, j);
            i++;
            j--;
        } else if (i == j) {
            if (item[i] >= split)
                i++;
            else
                j--;
        }
        /* LI */
    }
    /* i > j and LI */
 
    quicksort(item, Index, lo, j);
    quicksort(item, Index, i, hi);
    return;

} /* end quicksort */


int CSort(long *J, long *val, long maxval, long lo, long hi) {

    /* This function sorts the items J[lo..hi] by increasing value val,
       using a counting sort.
       Items with the same value retain the original order.
       maxval >= 0 is the maximum value that can occur;
       0 is the minimum value */

    long t, j, r, total, tmp, *start, *C;
    
    if (lo >= hi)
        return TRUE;

    C = (long *)malloc((hi-lo+1)*sizeof(long));
    start = (long *)malloc((maxval+1)*sizeof(long));
    
    if (C == NULL || start == NULL) {
        fprintf(stderr, "CSort(): Not enough memory!\n");
        return FALSE;
    }

    for (r=0; r<=maxval; r++)
        start[r] = 0;

    /* First pass. Count the number of items for each value. */
    for (t=lo; t<=hi; t++) {
        j = J[t];
        
        if (val[j] > maxval) {
            fprintf(stderr, "CSort(): val > maxval");
            return FALSE;
        }
        
        start[val[j]]++;
    }

    /* Make start cumulative */
    total = 0;
    for (r=0; r<=maxval; r++) {
        tmp = total;
        total += start[r];
        start[r] = tmp;
    }

    /* Second pass. Copy the items into C. */
    for (t=lo; t<=hi; t++) {
        j = J[t]; 
        C[start[val[j]]]= j;
        start[val[j]]++;
    }

    /* Third pass. Copy the items from C back into J. */
    for (t=lo; t<=hi; t++)
        J[t]= C[t-lo];

    free(start);
    free(C);
    
    return TRUE;
} /* end CSort */


void RandomPermute(long *X, long *Y, double *D, double *E, long lo, long hi) {

    /* This function randomly permutes the contents of the arrays
           X[lo..hi], Y[lo..hi] of type long and
           D[lo..hi], E[lo..hi] of type double, 
       all by the same random permutation.
       The function can also be used for < 4 arrays:
       e.g., if X=NULL on input, then X is not changed. */
       
    register long i;
  
    while (hi > lo) {
        i = Random1(lo, hi);
  
        if (X != NULL) {
            /*## swap: X[i] <-> X[hi] ##*/
            SwapLong(X, i, hi);
        }
        if (Y != NULL) {
            /*## swap: Y[i] <-> Y[hi] ##*/
            SwapLong(Y, i, hi);  
        }
        if (D != NULL) {
            /*## swap: D[i] <-> D[hi] ##*/
            SwapDouble(D, i, hi);  
        }

        if (E != NULL) { 
            /*## swap: E[i] <-> E[hi] ##*/
            SwapDouble(E, i, hi);   
        } 

        hi--;
    }      

} /* end RandomPermute */


