/* This file defines Gainbucket, a data structure
   which contains numbers of data items and their values. 
   The data item can be a vertex and the value its gain in a move.
   The numbers are integers >=0, and the values are integers
   without restriction. We call the items vertices and their values gains.

   The vertices are sorted in order of decreasing gain.
   Vertices with the same gain value are stored together in a bucket,
   implemented as a doubly linked list. The entries of a bucket,
   representing vertices, are called bucketentries.
   The list is terminated by NULL at both ends.

   The buckets themselves are also linked in a doubly linked list.
   This list is also terminated by NULL at both ends. */


#ifndef __GainBucket_h__
#define __GainBucket_h__

#include "Options.h" 

struct bucketentry {
    long vtxnr;            /* vertex number for this bucketentry */
    struct bucket *bucket; /* pointer to the bucket containing this entry */
    struct bucketentry *prev; /* pointer to previous bucketentry
                                  in the list of entries*/
    struct bucketentry *next; /* pointer to next bucketentry */
};
  
struct bucket {
    long value;           /* value for the entries in this bucket */
    struct bucket *prev;  /* pointer to previous bucket
                              in the list of buckets */
    struct bucket *next;  /* pointer to next bucket */
    struct bucketentry *entry; /* pointer to first bucketentry in this bucket */
};
  
struct gainbucket {
  long NrBuckets;       /* number of buckets in the gainbucket */
  struct bucket *Root;  /* pointer to the first bucket
                            in the list of buckets */
};


struct bucketentry *BucketInsert(struct gainbucket *pGB, long key, long VtxNr);
struct bucketentry *BucketMove(struct gainbucket *pGB, struct bucketentry *pE,
                               long key);
long BucketDeleteMax(struct gainbucket *pGB);

long GainBucketGetMaxVal(struct gainbucket *pGB);
long GainBucketGetMaxValVertexNr(struct gainbucket *pGB);

int ClearGainBucket(struct gainbucket *pGB);

#endif /* __GainBucket_h__ */
