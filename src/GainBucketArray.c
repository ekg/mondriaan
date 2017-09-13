#include "GainBucket.h"

int InitGainBucket(struct gainbucket *pGB, long MaxValue) {
    /* Initialize empty GainBucket structure.
       Must be called before using any other gainbucket method. */
    
    if (!pGB) {
        fprintf(stderr, "InitGainBucket(): Null parameter!\n");
        return FALSE;
    }
    
    pGB->MaxValue = MaxValue;
    pGB->MaxPresentValue = LONG_MIN;
    pGB->Root = (struct bucket*)calloc(2*MaxValue+1, sizeof(struct bucket));
    if (pGB->Root == NULL) {
        fprintf(stderr, "InitGainBucket(): Not enough memory!\n");
        return FALSE;
    }
    
    long k;
    for(k=-MaxValue; k<=MaxValue; ++k) {
        pGB->Root[pGB->MaxValue+k].value = k;
        pGB->Root[pGB->MaxValue+k].entry = NULL;
    }
    return TRUE;
} /* end InitGainBucket */


struct bucketentry *BucketInsert(struct gainbucket *pGB, long key, long VtxNr) {
    /* This function inserts vertex number VtxNr in the gainbucket GB,
       in the appropriate bucket with key value. The vertex must not be
       present already. A bucketentry representing VtxNr is allocated
       and a pointer to this bucketentry is returned. */

    struct bucket *pB;
    struct bucketentry *pE;
    
    if (!pGB || !pGB->Root) {
        fprintf(stderr, "BucketInsert(): Null parameter!\n");
        return NULL;
    }
    if(key > pGB->MaxValue || key < -pGB->MaxValue) {
        fprintf(stderr, "BucketInsert(): Invalid key!\n");
        return NULL;
    }
    
    pB = &(pGB->Root[pGB->MaxValue+key]);

    /*## Insert the vertex number into a new bucketentry at
         the start of the bucket: ##*/
    
    pE = (struct bucketentry *) malloc(sizeof(struct bucketentry));
    if (pE == NULL) {
        return NULL;
    }
  
    pE->vtxnr = VtxNr;
    pE->bucket = pB;  
    pE->prev = NULL;
    pE->next = pB->entry;
    
    if (pE->next != NULL) {
        (pE->next)->prev = pE;
    }
    else {
        pGB->NrBuckets++;
        if(key > pGB->MaxPresentValue)
            pGB->MaxPresentValue = key;
    }
    
    pB->entry = pE;
    
    return pE;
} /* end BucketInsert */


struct bucketentry *BucketMove(struct gainbucket *pGB, struct bucketentry *pE,
                               long key) {
    
   /* This function moves the vertex number represented by bucketentry pE
      in the gainbucket GB to the appropriate bucket with key value. 
      The key value must differ from the current value.  
      A pointer to the bucketentry is returned. */
   
    struct bucket *pB;
    
    if (!pGB || !pE || !pGB->Root) {
        fprintf(stderr, "BucketMove(): Null parameter!\n");
        return NULL;
    }
  
    pB = pE->bucket;
  
    if (pB->value == key) {
        return NULL;
    }
    if(key > pGB->MaxValue || key < -pGB->MaxValue) {
        fprintf(stderr, "BucketMove(): Invalid key!\n");
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
  
    /*## Check empty bucket: ##*/
    if (pB->entry == NULL) {
        
        pGB->NrBuckets--;
        /* If this bucket was the max bucket and the new value is lower, find the new max bucket */
        if(pB->value == pGB->MaxPresentValue && key < pGB->MaxPresentValue) {
            pGB->MaxPresentValue = LONG_MIN;
            for(pB--;pB>=pGB->Root;pB--) {
                if(pB->entry != NULL) {
                    pGB->MaxPresentValue = pB->value;
                    break;
                }
                if(pB->value < key)
                    break;
            }
        }
    }
    
    pB = &(pGB->Root[pGB->MaxValue+key]);

    /*## Reassign the bucketentry to the new bucket: ##*/
  
    pE->bucket = pB;  
    pE->prev = NULL;
    pE->next = pB->entry;
    
    if (pE->next != NULL) {
        (pE->next)->prev = pE;
    }
    else {
        pGB->NrBuckets++;
        if(key > pGB->MaxPresentValue)
            pGB->MaxPresentValue = key;
    }
    
    pB->entry = pE;
    
    return pE;
} /* end BucketMove */
  
  
long BucketDeleteMax(struct gainbucket *pGB) {
    
   /* This function deletes the first vertex from the first bucket.
      This vertex has the maximum gain value. 
      The function returns the vertex number.
      The first bucket must exist and it should not be empty.
      The memory space of the original bucketentry is freed. */

    long VtxNr;
    struct bucket *pB;
    struct bucketentry *pE;
    
    if (!pGB || !pGB->Root) {
        fprintf(stderr, "BucketDeleteMax(): Null parameter!\n");
        return -1;
    }
    
    pB = &(pGB->Root[pGB->MaxValue+pGB->MaxPresentValue]) ;
    pE = pB->entry;  
    
    if (!pB || !pE) {
        return -1;
    }
  
    VtxNr = pE->vtxnr;
  
    /*## Remove this entry from the bucket: ##*/
    pB->entry = pE->next;
   
    if (pE->next != NULL)
        (pE->next)->prev = NULL;
    
    free(pE);
  
    /*## If this bucket is now empty, find the new max bucket: ##*/
    if (pB->entry == NULL) {
        pGB->NrBuckets--;
        pGB->MaxPresentValue = LONG_MIN;
        for(pB--;pB>=pGB->Root;pB--) {
            if(pB->entry != NULL) {
                pGB->MaxPresentValue = pB->value;
                break;
            }
        }
    }
    
    return(VtxNr);
} /* end BucketDeleteMax */


long GainBucketGetMaxVal(struct gainbucket *pGB) {

    /* This function gives the value of the vertex with the maximum
       gain value, if the gainbucket data structure is not empty.
       Otherwise, it returns LONG_MIN. */
    
    if (!pGB || !pGB->Root) {
        return LONG_MIN;
    }

    if (pGB->NrBuckets > 0)
        return((pGB->Root[pGB->MaxValue+pGB->MaxPresentValue]).value);
    else
        return(LONG_MIN);

} /* end GainBucketGetMaxVal */


long GainBucketGetMaxValVertexNr(struct gainbucket *pGB) {
    /* This function gives the number of the vertex with the maximum
       gain value, if the gainbucket data structure is not empty.
       Otherwise, it returns LONG_MIN. */
    
    if (!pGB || !pGB->Root) {
        return LONG_MIN;
    }
    
    if (pGB->NrBuckets > 0)
        return(((pGB->Root[pGB->MaxValue+pGB->MaxPresentValue]).entry)->vtxnr);
    else
        return(LONG_MIN);

} /* end GainBucketGetMaxValVertexNr */


int ClearGainBucket(struct gainbucket *pGB) {
   /* This function deletes all bucketentries
      and frees the corresponding memory space.
      As a result, pGB->NrBuckets = 0. */

    struct bucket *pB;
    struct bucketentry *pE;
    
    if (!pGB) {
        fprintf(stderr, "ClearGainBucket(): Null parameter!\n");
        return FALSE;
    }
    
    if(pGB->Root == NULL) {
        return TRUE;
    }
    
    for(pB=&(pGB->Root[2*pGB->MaxValue]); pB>=pGB->Root; pB--) {
        if(pB->entry == NULL) {
            continue;
        }
        
        /*## Remove all entries from this bucket: ##*/  
        while ((pE = pB->entry) != NULL) {
            pB->entry = pE->next;
            free(pE);
        }
        
        pGB->NrBuckets--;
    }
    pGB->MaxPresentValue = LONG_MIN;
    
    return TRUE;
} /* end ClearGainBucket */

int DeleteGainBucket(struct gainbucket *pGB) {
   /* This function deletes the GainBucket structure
      and frees all corresponding memory space.
      As a result, pGB->Root = NULL and the structure
      cannot be used any more until InitGainBucket() is
      called. */
    if (!pGB) {
        fprintf(stderr, "DeleteGainBucket(): Null parameter!\n");
        return FALSE;
    }
    
    if(!ClearGainBucket(pGB)) {
        return FALSE;
    }
    
    if(pGB->Root == NULL)
        return TRUE;
    
    free(pGB->Root);
    pGB->Root = NULL;
    
    return TRUE;
} /* end DeleteGainBucket */
