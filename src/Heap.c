
#include <string.h>
#include "Heap.h"

/* This is our own implementation of a heap, with some ideas inspired by
 *  -> https://gist.github.com/martinkunev/1365481
 *  -> https://github.com/robin-thomas/max-heap/blob/master/maxHeap.c
 */

#define LCHILD(x) (2*(x) + 1)
#define RCHILD(x) (2*(x) + 2)
#define PARENT(x) (((x) - 1)/2)


/**
 * Initialize a heap struct
 * A heap consists of items, and offers the ability to retrieve
 * the maximum/minimum item efficiently, depending on the comparison
 * function given.
 * 
 * Input:
 * itemSize         : The number of bytes one item consists of
 * compare          : Comparison function to compare items
 * 
 * Input/Output:
 * pHeap            : The heap struct
 */
void HeapInit(struct heap *pHeap, long itemSize, int (*compare)(const void *, const void*)) {
    pHeap->numItems = 0;
    pHeap->size = 0;
    pHeap->itemSize = itemSize;
    pHeap->items = NULL;
    pHeap->compare = compare;
}

/**
 * Free a heap struct
 * 
 * Input/Output:
 * pHeap            : The heap struct
 */
void HeapDestroy(struct heap *pHeap) {
    pHeap->numItems = 0;
    pHeap->size = 0;
    if(pHeap->items != NULL) {
        free(pHeap->items);
    }
    pHeap->items = NULL;
}

/**
 * Check the size of a heap. If the number of items in it
 * is equal to the maximum size, extend the maximum size.
 * 
 * Input/Output:
 * pHeap            : The heap struct
 */
void HeapCheckSize(struct heap *pHeap) {
    if(pHeap->items == NULL) {
        pHeap->size = 8;
        pHeap->items = malloc(pHeap->size * pHeap->itemSize);
        if(pHeap->items == NULL) {
            fprintf(stderr, "HeapPush(): Out of memory!\n");
            exit(1);
        }
    }
    else if(pHeap->numItems >= pHeap->size) {
        pHeap->size *= 2;
        char *newAlloc = realloc(pHeap->items, pHeap->size * pHeap->itemSize);
        if(newAlloc == NULL) {
            fprintf(stderr, "HeapPush(): Out of memory!\n");
            exit(1);
        }
        pHeap->items = newAlloc;
    }
}

/**
 * Copy the top element to item.
 * 
 * Input:
 * pHeap            : The heap struct
 * 
 * Output:
 * item             : Memory space to copy the top item to
 * 
 * Returns FALSE on error, TRUE otherwise
 */
int HeapPeek(struct heap *pHeap, void * const item) {
    if(pHeap->numItems <= 0 || item == NULL)
        return FALSE;
    
    memcpy(item, &pHeap->items[0], pHeap->itemSize);
    return TRUE;
}

/**
 * Copy the top element to item, and remove it from the heap.
 * 
 * Input/Output:
 * pHeap            : The heap struct
 * 
 * Output:
 * item             : Memory space to copy the top item to
 * 
 * Returns FALSE on error, TRUE otherwise
 */
int HeapPop(struct heap *pHeap, void * const item) {
    if(pHeap->numItems <= 0)
        return FALSE;
    
    /* Special case: the heap becomes empty */
    if(pHeap->numItems == 1) {
        --pHeap->numItems;
        
        if(item != NULL) {
            memcpy(item, &pHeap->items[0], pHeap->itemSize);
        }
        return TRUE;
    }
    
    /* The heap consists of at least two items; we can safely call HeapReplace */
    return HeapReplace(pHeap, &pHeap->items[(--pHeap->numItems) * pHeap->itemSize], item);
}

/**
 * Insert a new item into the heap
 * 
 * Input/Output:
 * pHeap            : The heap struct
 * 
 * Input:
 * new              : The item to insert into the heap
 */
void HeapPush(struct heap *pHeap, void * const item) {
    HeapCheckSize(pHeap);
    
    memcpy(&pHeap->items[pHeap->numItems * pHeap->itemSize], item, pHeap->itemSize);
    HeapSiftUp(pHeap, pHeap->numItems);
    
    ++pHeap->numItems;
}

/**
 * Copy the top element to item, and replace it by another item.
 * 
 * Input/Output:
 * pHeap            : The heap struct
 * 
 * Input:
 * new              : The item to replace the top item with
 * 
 * Output:
 * old              : Memory space to copy the old top item to
 * 
 * Returns FALSE on error, TRUE otherwise
 */
int HeapReplace(struct heap *pHeap, void * const new, void * const old) {
    if(pHeap->numItems <= 0)
        return FALSE;
    
    if(old != NULL) {
        memcpy(old, &pHeap->items[0], pHeap->itemSize);
    }
    
    if(&pHeap->items[0] != new) {
        memcpy(&pHeap->items[0], new, pHeap->itemSize);
    
        HeapSiftDown(pHeap, 0);
    }
    
    return TRUE;
}

/**
 * Sift up an item to its right position. This is used when an element at a leaf is modified.
 * 
 * Input:
 * item             : The item to sift up
 * 
 * Input/Output:
 * pHeap            : The heap struct
 */
void HeapSiftUp(struct heap *pHeap, long item) {
    
    char *base_ptr = pHeap->items;
    size_t itemSize = pHeap->itemSize;
    
    long index,     /* The item we are currently bubbling at */
         swap;      /* The item to swap the bubble with */
    
    char value[itemSize];
    
    memcpy(value, &base_ptr[itemSize * item], itemSize);
    
    for(index = item; TRUE; index = swap) {
        if(index == 0)
            break;
        
        swap = PARENT(index);
        
        if(pHeap->compare(value, &base_ptr[itemSize * swap]) <= 0) {
            break;
        }
        
        memcpy(&base_ptr[itemSize * index], &base_ptr[itemSize * swap], itemSize);
    }
    
    if(index != item) {
        memcpy(&base_ptr[itemSize * index], value, itemSize);
    }
    
}

/**
 * Sift down an item to its right position. This is used when the element at the root is modified.
 * 
 * Input:
 * item             : The item to sift down
 * 
 * Input/Output:
 * pHeap            : The heap struct
 */
void HeapSiftDown(struct heap *pHeap, long item) {
    
    char *base_ptr = pHeap->items;
    size_t itemSize = pHeap->itemSize;
    
    long index,     /* The item we are currently bubbling at */
         swap,      /* The item to swap the bubble with */
         childL,    /* Children of the bubble */
         childR;
    
    char value[itemSize];
    
    memcpy(value, &base_ptr[itemSize * item], itemSize);
    
    for(index = item; TRUE; index = swap) {
        
        childL = LCHILD(index);
        if(childL >= pHeap->numItems)
            break; /* This node has no children */
        childR = RCHILD(index);
        
        swap = childL;
        if(childR < pHeap->numItems && pHeap->compare(&base_ptr[itemSize * childR], &base_ptr[itemSize * childL]) >= 0) {
            swap = childR;
        }
        
        if(pHeap->compare(value, &base_ptr[itemSize * swap]) >= 0) {
            break;
        }
        
        memcpy(&base_ptr[itemSize * index], &base_ptr[itemSize * swap], itemSize);
        
    }
    
    if(index != item) {
        memcpy(&base_ptr[itemSize * index], value, itemSize);
    }
    
    
}

/**
 * Turn an array into a heap
 * 
 * Input:
 * numItems         : The number of items in the array
 * 
 * Input/Output:
 * items            : The array (will likely be permutated)
 * pHeap            : The heap struct
 */
void Heapify(struct heap *pHeap, void * const items, long numItems) {
    pHeap->items = (char *)items;
    pHeap->numItems = numItems;
    pHeap->size = numItems;
    
    long item; /* The item we want to move into place */
    
    for(item = PARENT(pHeap->numItems-1); item >= 0; --item) {
        HeapSiftDown(pHeap, item);
    }
    
}


#undef LCHILD
#undef RCHILD
#undef PARENT
