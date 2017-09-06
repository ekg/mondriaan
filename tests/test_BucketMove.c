#include "GainBucket.h"

int main(int argc, char **argv) {

    void CheckBuckets2(struct gainbucket GB, long n, long bs);

    long n, bs, i, j, newkey;
    struct gainbucket GB ;
    struct bucket *pB ;
    struct bucketentry *pE ;

    printf("Test BucketMove: ");
    n= 100 ; /* number of buckets */
    bs = 5; /* bucket size */

    GB.NrBuckets =0 ;
    GB.Root = NULL ;
    
    if ( !InitGainBucket(&GB, n) ) {
        printf("Error\n");
        exit(1);
    }

    /* Insert vertices in increasing order of gain value */
    for (j=0; j<bs; j++)
        for (i=0; i<n; i++)
            BucketInsert(&GB, i, j*n + i) ;

    /* Move vertices to new buckets -1, -2, ... -n, while reversing the order 
       inside the buckets. These buckets are created at the end of the
       current bucket list */

#ifdef GAINBUCKET_ARRAY
    pB = &(GB.Root[GB.MaxValue+GB.MaxPresentValue]) ;
#else
    pB = GB.Root ;
#endif
    while ( pB != NULL && pB->value >= 0 ) {
        pE = pB->entry ; /* all buckets are nonempty */
        newkey = pB->value - n ;
        BucketMove(&GB, pE, newkey) ;
#ifdef GAINBUCKET_ARRAY
        pB = &(GB.Root[GB.MaxValue+GB.MaxPresentValue]) ;
#else
        pB = GB.Root ;
#endif
    }

    CheckBuckets2(GB, n, bs) ;
    
    if ( !DeleteGainBucket(&GB) ) {
        printf("Error\n");
        exit(1);
    }

    printf("OK\n") ;
    exit(0);

} /* end main */

  
void CheckBuckets2(struct gainbucket GB, long n, long bs) {
    long nrentries = 0, i, j ;
    struct bucket *pB ;
    struct bucketentry *pE ;
 
    if (GB.NrBuckets != n){
        printf("Error\n") ;
        exit(1); 
    }
 
#ifdef GAINBUCKET_ARRAY
    pB = &(GB.Root[GB.MaxValue+GB.MaxPresentValue]) ;
#else
    pB = GB.Root ;
#endif
 
    i = 0;
    while ( pB != NULL ) {
        j = 0;
        pE = pB->entry ; 
        while ( pE != NULL ) {   
            if (  (pE->bucket)->value != -1-i  /* = n-1-i - n */
                   || pE->vtxnr != n-1-i + j*n ) {
                 printf("Error\n") ;
                 exit(1);
            }

            nrentries++ ;
            j++ ;
            pE = pE->next ;
        }
#ifdef GAINBUCKET_ARRAY
        pB = (pB->value > -GB.MaxValue) ? pB-1 : NULL ;
#else
        pB = pB->next ;
#endif
        if (j!= bs){  
            printf("Error\n") ;
            exit(1);
        }
        i++ ;
    }
 
    if (nrentries != n*bs){
        printf("Error\n") ;
        exit(1);
    }
 
} /* end CheckBuckets2 */

