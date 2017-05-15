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
  
    quicksort3way(X, I, LengthX);
  
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
 
    /* Sort i0, i1, i2 in decreasing order of item value */
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




/**
 * Below, an implementation of 3-way quicksort is included.
 * This implementation combines the code from
 * https://sourceware.org/git/?p=glibc.git;a=blob;f=stdlib/qsort.c;h=12a5a7506a3337ecbebc6b1778425208ef8439c3;hb=HEAD
 * with knowledge taken from
 *     Engineering a sort function; Jon Bentley and M. Douglas McIlroy;
 *     Software - Practice and Experience; Vol. 23 (11), 1249-1265, 1993.
 *
 * In essence, the implementation below is an altered version of the glibc implementation, specialized for long values,
 * with additional logic to provide 3-way sorting functionality. Furthermore, pivots are chosen at random instead of
 * the median of 3 values at fixed positions, to provide better randomization of equal elements.
 * Random numbers are generated directly useing rand() instead of Random1(), to reduce the number of function calls.
 */


#define LT_GT >

/* Swap two items. */
#define SWAP(A, B)		{const long C = *(A); *(A) = *(B); *(B) = C; const long D = indices[A-list]; indices[A-list] = indices[B-list]; indices[B-list] = D;}

/* Stack node declarations used to store unfulfilled partition obligations. */
typedef struct
{
	long *lo;
	long *hi;
} stack_node;

/* The next 4 #defines implement a very fast in-line stack abstraction. */
/* The stack needs log (total_elements) entries.
 * Since total_elements has type size_t, we get as
 * upper bound for log (total_elements):
 * bits per byte (CHAR_BIT) * sizeof(size_t).
 */
#define STACK_SIZE	(CHAR_BIT * sizeof(size_t))
#define PUSH(low, high)	((void) ((top->lo = (low)), (top->hi = (high)), ++top))
#define	POP(low, high)	((void) (--top, (low = top->lo), (high = top->hi)))
#define	STACK_NOT_EMPTY	(stack < top)

/**
 * In each iteration, we execute a part of a recursion of the quicksort algorithm.
 * The list in each iteration can be viewed as follows:
 * +-----+-----+-----+-----+-----+
 * [  =  [  <  [  ?  ]  >  ]  =  ]
 * +-----+-----+-----+-----+-----+
 * a     b     c     d     e     f
 * 
 * After processing [c,d], we have:
 * +-----+-----++-----+-----+
 * [  =  [  <  ][  >  ]  =  ]
 * +-----+-----++-----+-----+
 * a     b     dc     e     f
 * 
 * After moving the equal elements, we have:
 * +-----+-----+-----+
 * [  <  ]  =  [  >  ]
 * +-----+-----+-----+
 * b=a   d     c   f=e
 * 
 * Note that this (and all comments below) rely on LT_GT == < . In practice, Mondriaan
 * chooses LT_GT to be >, so 'higher' and 'lower' should be swapped in all comments below.
 * 
 * Here:
 *  - [a,b) denotes all elements equal to the pivot, as found by c,
 *  - [b,c) denotes all elements lower than the pivot,
 *  - [c,d] denotes all elements to be processed,
 *  - (d,e] denotes all elements higher than the pivot,
 *  - (e,f] denotes all elements equal to the pivot, as found by d.
 * 
 * In each iteration, three random elements are selected from the then to be sorted set [c,d], of which the
 * median is determined. This element serves as pivot. Any elements equal to the pivot are sorted into
 * [a,b[ or ]e,f], while elements lower than or higher than the pivot go into respectively [b,c[ and ]d,e].
 * After [c,d] is empty (d < c), [a,b[ and ]e,f] are moved in between [b,d] and [c,e].
 * 
 */

void quicksort3way (long *list, long *indices, size_t total_elems)
{
	if (total_elems <= 1)
		return;

	long *a, *b, *c, *d, *e, *f;
	long n, el1, el2, el3, pivot;
	
	/* Initialise stack */
	stack_node stack[STACK_SIZE];
	stack_node *top = stack;
	PUSH (NULL, NULL);
	
	/* Initialise first iteration */
	c = list;
	d = &list[total_elems - 1];

	while (STACK_NOT_EMPTY)
	{
		/* Process [c,d] */
		a = b = c;
		f = e = d;
		n = d-c+1;
		
		if(n == 2) {
			/* Two elements are easily sorted */
			if(*d LT_GT *c) {
				SWAP(c, d);
			}
			++c;
			--d;
		}
		else if(n > 2) {
			
			/*** DETERMINE PIVOT ***/
			
			el1 = b[(long)((rand() / (RAND_MAX+1.0))*n)];
			el2 = b[(long)((rand() / (RAND_MAX+1.0))*n)];
			el3 = b[(long)((rand() / (RAND_MAX+1.0))*n)];
			if(el1 < el2) {
				if(el2 < el3)
					pivot = el2;
				else if(el1 < el3)
					pivot = el3;
				else
					pivot = el1;
			}
			else {
				if(el1 < el3)
					pivot = el1;
				else if(el2 < el3)
					pivot = el3;
				else
					pivot = el2;
			}
			
			/*** PARTITION ***/
			
			/* Check whether the start/end contains elements equal to the pivot */
			while(b <= e && *b == pivot) {
				++b;
				++c;
			}
			while(b <= e && *e == pivot) {
				--d;
				--e;
			}
			
			/* The `collapse the walls' section of quicksort. */
			while (c <= d)
			{
				/* Check for equals */
				while(c <= d && pivot == *c) {
					SWAP(b, c)
					++b;
					++c;
				}
				
				while(c <= d && pivot == *d) {
					SWAP(e, d)
					--e;
					--d;
				}
				
				/* Move c and d until they have to be swapped */
				while (c <= d && *c LT_GT pivot)
					++c;
				
				while (c <= d && pivot LT_GT *d)
					--d;
				
				/* With LT_GT == <, we now have *d <= pivot <= *c */
				
				if (c < d)
				{
					/* Now, *d <= pivot <= *c and c<d, hence swap. */
					SWAP (c, d);
					
					/* Move elements equal to pivot */
					if(pivot == *c) {
						SWAP(b, c)
						++b;
					}
					if(pivot == *d) {
						SWAP(e, d)
						--e;
					}
					
					/* Mark c and d as processed */
					++c;
					--d;
				}
				else if (c == d)
				{
					/* By the while-loops, we have *d <= pivot <= *c.
					 * Hence, if c == d, we have *d == *c == pivot.
					 * As the next step will be swapping back the elements equal to the pivot,
					 * we do not swap this element to an equal area.
					 * We move both c and d, effectively marking this element as equal to the pivot.
					 */
					++c;
					--d;
					break;
				}
				/* else, c > d and the loop will break */
			}
			/* Now, c-d is either 1 (most of times) or 2 (if (c==d) occurred above).
			 * Anyhow, d<c and the set of elements to be processed [c,d] is empty.
			 */
			
			/* If there exists a left/right partition (i.e., [b,c) or (d,e] is non-empty),
			 * move the left/right equal-partition (i.e., [a,b) or (e,f]) between c and d. */
			if(b <= d) {
				while(b > a) {
					--b;
					SWAP(d, b);
					--d;
				}
			} /* else, there exists no left partition, so moving equal elements is unnecessary. */
			
			if(e >= c) {
				while(e < f) {
					++e;
					SWAP(c, e);
					++c;
				}
			} /* else, there exists no right partition, so moving equal elements is unnecessary. */
		
		}
		else {
			fprintf(stderr, "quicksort3way(): Single element in set\n");
			exit(-1);
		}

		/*** RECURSION ***/
		
		/* Set up pointers for next iteration. First determine whether one or
		 * both partitions are completely sorted. If so, ignore one or both.
		 * Otherwise, push the larger partition's bounds on the stack and
		 * continue sorting the smaller one.
		 * The next set to be partitioned is assigned to [c,d].
		 */

		if (d <= b)
		{
			if (c >= e)
				/* 'Greater than' [c,e] and 'smaller than' [b,d] partition both contain less than two elements.
				 * These are sorted, so ignore both.
				 */
				POP (c, d);
			else
				/* 'Smaller than' partition [b,d] contains less than two elements.
				 * It is sorted, so ignore it, and process [c,e].
				 */
				d = e;
		}
		else if (c >= e)
			/* 'Greater than' partition [c,e] contains less than two elements.
			 * It is sorted, so ignore it, and process [b,d].
			 */
			c = b;
		else if ((d - b) > (e - c))
		{
			/* Push larger 'Smaller than' partition [b,d] and process [c,e]. */
			PUSH (b, d);
			d = e;
		}
		else
		{
			/* Push larger 'Greater than' partition [c,e] and process [b,d]. */
			PUSH (c, e);
			c = b;
		}
	}

}


#undef LT_GT
#undef SWAP
#undef STACK_SIZE
#undef PUSH
#undef POP
#undef STACK_NOT_EMPTY


