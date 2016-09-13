#include "GainBucket.h"

struct bucketentry *BucketInsert(struct gainbucket *pGB, long key, long VtxNr) {

    /* This function inserts vertex number VtxNr in the gainbucket GB,
       in the appropriate bucket with key value. The vertex must not be
       present already. A bucketentry representing VtxNr is allocated
       and a pointer to this bucketentry is returned. */

    struct bucket *pB, **ppB;
    struct bucketentry *pE;
    
    if (!pGB) {
        fprintf(stderr, "BucketInsert(): Null argument!\n");
        return NULL;
    }
    
    ppB = &(pGB->Root);
    pB = *ppB;
  
    /*## While the key is smaller than the value in the buckets, 
         go through the bucket list: ##*/
    while (*ppB != NULL) {
        pB = *ppB; /* pB points to current bucket */
  
        if (key < pB->value)
            ppB = &(pB->next); /* points to next bucket */
        else 
            break;
    }
  
    if (pB == NULL) {
        /*## Bucket list is empty. Create first bucket: ##*/
  
        *ppB = (struct bucket *) malloc(sizeof(struct bucket));
        
        if (*ppB == NULL) {
            fprintf(stderr, "BucketInsert(): Not enough memory for first bucket!\n");
            return NULL;
        }
        
        (*ppB)->value = key;
        (*ppB)->entry = NULL;
        (*ppB)->prev = (*ppB)->next = NULL;
  
        pGB->NrBuckets++;
    } else if (key < pB->value) { 
        /*## Create new bucket at the end of the current list: ##*/
        /* *ppB == NULL */
  
        *ppB = (struct bucket *) malloc(sizeof(struct bucket));
        
        if (*ppB == NULL) {
            fprintf(stderr, "BucketInsert(): Not enough memory for new bucket!\n");
            return NULL;
        }
  
        (*ppB)->value = key;
        (*ppB)->entry = NULL;
        (*ppB)->prev = pB;
        (*ppB)->next = NULL;
  
        pGB->NrBuckets++;  
    } else if (key > pB->value) {
        /*## Create new bucket between the previous and the current
             bucket in the list: ##*/
  
        *ppB = (struct bucket *) malloc(sizeof(struct bucket));
        
        if (*ppB == NULL) {
            fprintf(stderr, "BucketInsert(): Not enough memory for new bucket!\n");
            return NULL;
        }
      
        (*ppB)->value = key;
        (*ppB)->entry = NULL;
        (*ppB)->prev = pB->prev;
        (*ppB)->next = pB;
     
        if (pB->prev != NULL)
            (pB->prev)->next = *ppB;
        else
            pGB->Root = *ppB;
      
        pB->prev = *ppB;
        pGB->NrBuckets++;
    }
        /* otherwise key == pB->value (and hence pB == *ppB) */


    /*## Insert the vertex number into a new bucketentry at
         the start of the bucket: ##*/
    pE = (struct bucketentry *) malloc(sizeof(struct bucketentry));
    
    if (pE == NULL) {
        fprintf(stderr, "BucketInsert(): Not enough memory for bucket entry!\n");
        return NULL;
    }
  
    pE->vtxnr = VtxNr;
    pE->bucket = *ppB;  
    pE->prev = NULL;
    pE->next = (*ppB)->entry;
    
    if (pE->next != NULL)
        (pE->next)->prev = pE;
    
    (*ppB)->entry = pE;

    return pE;
} /* end BucketInsert */


struct bucketentry *BucketMove(struct gainbucket *pGB, struct bucketentry *pE,
                               long key) {
  
   /* This function moves the vertex number represented by bucketentry pE
      in the gainbucket GB to the appropriate bucket with key value. 
      The key value must differ from the current value.  
      A pointer to the new bucketentry is returned.
      The memory space of the original bucketentry is freed.
      If its bucket has become empty, the space of the bucket is also freed. */
  
    long VtxNr;
    struct bucket *pB;
    
    if (!pGB || !pE) {
        fprintf(stderr, "BucketMove(): Null arguments!\n");
        return NULL;
    }
  
    VtxNr = pE->vtxnr;
    pB = pE->bucket;
  
    if (pB->value == key) {
        fprintf(stderr, "BucketMove(): Destination bucket equals source!\n");
        return NULL;
    }
  
    /*## Remove this entry from the bucket: ##*/
    /* Adjust forward links */
    if (pE->prev != NULL)
        (pE->prev)->next = pE->next;  
    else
        pB->entry = pE->next;
  
    /* Adjust backward links */ 
    if (pE->next != NULL)
        (pE->next)->prev = pE->prev;
  
    if (pE != NULL) 
        free(pE);
  
    /*## If this bucket is now empty, remove it from the list: ##*/
    if (pB->entry == NULL) {      
        if (pB->prev != NULL)
             (pB->prev)->next = pB->next;
        else
             pGB->Root = pB->next;
  
        if (pB->next != NULL)
             (pB->next)->prev = pB->prev;
        
        if (pB != NULL)
          free(pB);
        pGB->NrBuckets--;
    }
    
    /*## Insert the vertex in the bucket of the new key: ##*/
    pE = BucketInsert(pGB, key, VtxNr);
  
    return pE;
} /* end BucketMove */
  
  
long BucketDeleteMax(struct gainbucket *pGB) {
  
   /* This function deletes the first vertex from the first bucket.
      This vertex has the maximum gain value. 
      The function returns the vertex number.
      The first bucket must exist and it should not be empty.
      The memory space of the original bucketentry is freed.
      If its bucket has become empty, the space of the bucket is also freed. */
 

    long VtxNr;
    struct bucket *pB;
    struct bucketentry *pE;
    
    if (!pGB) {
        fprintf(stderr, "BucketDeleteMax(): Null arguments!\n");
        return -1;
    }
    
    pB = pGB->Root ;
    pE = pB->entry;  
    
    if (!pB || !pE) {
        fprintf(stderr, "BucketDeleteMax(): Internal error!\n");
        return -1;
    }
  
    VtxNr = pE->vtxnr;
  
    /*## Remove this entry from the bucket: ##*/
    pB->entry = pE->next;
   
    if (pE->next != NULL)
        (pE->next)->prev = NULL;
  
    free(pE);
  
    /*## If this bucket is now empty, remove it from the list: ##*/
    if (pB->entry == NULL) {      
        pGB->Root = pB->next;
  
        if (pB->next != NULL)
            (pB->next)->prev = NULL;
        
        free(pB);
        pGB->NrBuckets--;
    }
  
    return(VtxNr);
} /* end BucketDeleteMax */


long GainBucketGetMaxVal(struct gainbucket *pGB) {

    /* This function gives the value of the vertex with the maximum
       gain value, if the gainbucket data structure is not empty.
       Otherwise, it returns LONG_MIN. */
    if (!pGB) return LONG_MIN;

    if (pGB->NrBuckets > 0)
        return((pGB->Root)->value);
    else
        return(LONG_MIN);

} /* end GainBucketGetMaxVal */


long GainBucketGetMaxValVertexNr(struct gainbucket *pGB) {

    /* This function gives the number of the vertex with the maximum
       gain value, if the gainbucket data structure is not empty.
       Otherwise, it returns LONG_MIN. */
    
    if (!pGB) return LONG_MIN;
    
    if (pGB->NrBuckets > 0)
        return(((pGB->Root)->entry)->vtxnr);
    else
        return(LONG_MIN);

} /* end GainBucketGetMaxValVertexNr */


int ClearGainBucket(struct gainbucket *pGB) {

   /* This function deletes all vertices and buckets
      and frees the corresponding memory space.
      As a result, pGB->Root = NULL and pGB->NrBuckets = 0. */

    struct bucket *pB;
    struct bucketentry *pE;
    
    if (!pGB) {
        fprintf(stderr, "ClearGainBucket(): Null argument!\n");
        return FALSE;
    }
    
    while ((pB = pGB->Root) != NULL) {      
        pGB->Root = pB->next;
  
        /*## Remove all entries from this bucket: ##*/  
        while ((pE = pB->entry) != NULL) {
            pB->entry = pE->next;
            free(pE);
        }
  
        /*## Remove this bucket from the list: ##*/
        free(pB);
        pGB->NrBuckets--;
    }

    return TRUE;
} /* end ClearGainBucket */
