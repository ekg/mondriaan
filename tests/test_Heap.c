#include "Heap.h"
#include "Sort.h"
#include <string.h>

#define max(x,y) ( ( (x)>(y) ) ? (x) : (y) )


void test_Heapify();
void test_HeapPush();

int main(int argc, char **argv) {

    printf("Test Heap: ");
    
    SetRandomSeed(111);
    
    test_Heapify();
    test_HeapPush();

    printf("OK\n");
    exit(0);

} /* end main */


int compareLong (const void *a, const void *b) {
    long diff = *((long*)a) - *((long*)b);
    if(diff == 0)
        return 0;
    return (diff < 0) ? -1 : 1;
} /* end compareLong */

int compare2ndLong (const void *a, const void *b) {
    long diff = ((long*)a)[1] - ((long*)b)[1];
    if(diff == 0)
        return 0;
    return (diff < 0) ? -1 : 1;
} /* end compare2ndLong */



/**
 * Test Heapify()
 */
void test_Heapify() {
    long N = 100, i;
    
    long *items = (long*)malloc(2*N*sizeof(long));
    long *heapItems = (long*)malloc(2*N*sizeof(long));
    long tmp[2];
    
    if (items == NULL || heapItems  == NULL) {
        printf("Error\n");
        exit(1);
    }
    
    for(i=0; i<N; ++i) {
        items[2*i  ] = i;
        items[2*i+1] = Random1(0, 50);
    }
    
    memcpy(heapItems, items, 2*N*sizeof(long));
    
    struct heap Heap;
    HeapInit(&Heap, 2*sizeof(long), compare2ndLong);
    Heapify(&Heap, heapItems, N);
    
    qsort(items, N, 2*sizeof(long), compare2ndLong);
    
    for(i=0; i<N; ++i) {
        HeapPop(&Heap, tmp);
        if(tmp[1] != items[2*N-(2*i)-1]) {
            printf("Error\n");
            
        }
    }
    
    if(Heap.numItems > 0) {
        printf("Error\n");
        exit(1);
    }
    
    HeapDestroy(&Heap);
    free(items);
    
} /* end Heapify */

/**
 * Test HeapPush()
 */
void test_HeapPush() {
    
    struct heap Heap;
    HeapInit(&Heap, sizeof(long), compareLong);
    long tmp[1], tmp2[1], i, N = 100, max = -1;
    long *items = (long*)malloc(N*sizeof(long));
    
    if (items == NULL) {
        printf("Error\n");
        exit(1);
    }
    
    for(i=0; i<N; ++i) {
        *tmp = items[i] = Random1(0,N);
        max = max(max, *tmp);
        HeapPush(&Heap, tmp);
        if(!HeapPeek(&Heap, tmp2) || *tmp2 != max) {
            printf("Error\n");
            exit(1);
        }
    }
    
    qsort(items, N, sizeof(long), compareLong);
    
    for(i=0; i<N; ++i) {
        HeapPop(&Heap, tmp);
        if(*tmp != items[N-i-1]) {
            printf("Error\n");
            
        }
    }
    
    if(Heap.numItems > 0) {
        printf("Error\n");
        exit(1);
    }
    
    HeapDestroy(&Heap);
    free(items);
    
} /* end HeapPush */
