#include "GainBucket.h"

int main(int argc, char **argv) {

    void CheckBuckets(struct gainbucket GB, long n, long bs);

    long n, bs, i, j, offset;
    struct gainbucket GB ;

    printf("Test BucketInsert: ");
    n= 100 ; /* number of buckets */
    bs = 5; /* bucket size */
    offset = n/2 ;

    /* Insert vertices in increasing order of gain value */
    GB.NrBuckets =0 ;
    GB.Root = NULL ;
    for (j=0; j<bs; j++)
        for (i=0; i<n; i++)
            BucketInsert(&GB, i-offset, j*n + i) ;

    CheckBuckets(GB, n, bs) ;
    
    printf("OK\n") ;
    exit(0);

} /* end main */
  
void CheckBuckets(struct gainbucket GB, long n, long bs) {
    long nrentries = 0, i, j, offset ;
    struct bucket *pB ;
    struct bucketentry *pE ;
 
    /* printf("Number of Buckets: %ld \n", GB.NrBuckets) ; */
    if (GB.NrBuckets != n){
        printf("Error\n") ;
        exit(1);
    }

    pB = GB.Root ;

    offset = n/2 ; 
    i = 0;
    while ( pB != NULL ) {
        j = 0;
        pE = pB->entry ;
        while ( pE != NULL ) {
            /* printf("B[% ld] vntnr: %ld \n", (pE->bucket)->value,
                                                pE->vtxnr) ;  */
            if ( (pE->bucket)->value != n-1-i  - offset ||
                 pE->vtxnr != n-1-i + (bs-1-j) *n ) {
                 printf("Error\n") ;
                 exit(1);
            }
 
            nrentries++ ;
            j++ ;
            pE = pE->next ;
        }
        pB = pB->next ;
        if (j!= bs){  
            printf("Error\n") ;
            exit(1); 
        }
        i++ ;
    }
 
    /* printf("TotalNrEntries: %ld\n", nrentries) ; */
    if (nrentries != n*bs){
        printf("Error\n") ; 
        exit(1);
    } 

} /* end CheckBuckets */

