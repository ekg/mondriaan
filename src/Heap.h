#ifndef __Heap_h__
#define __Heap_h__

#include "Sort.h"

struct heap {
    long numItems;
    long size;
    
    char *items;
    size_t itemSize;
    int (*compare)(const void *, const void*);
};

/* Function declarations for Heap.c */
void HeapInit(struct heap *pHeap, long itemSize, int (*compare)(const void *, const void*));
void HeapDestroy(struct heap *pHeap);
void HeapCheckSize(struct heap *pHeap);
int HeapPeek(struct heap *pHeap, void * const item);
int HeapPop(struct heap *pHeap, void * const item);
void HeapPush(struct heap *pHeap, void * const item);
int HeapReplace(struct heap *pHeap, void * const new, void * const old);
void HeapSiftUp(struct heap *pHeap, long item);
void HeapSiftDown(struct heap *pHeap, long item);
void Heapify(struct heap *pHeap, void * const items, long numItems);

#endif /* __Heap_h__ */

