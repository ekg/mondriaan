/* This file defines Gainbucket, a data structure
   which contains numbers of data items and their values. 
   This file is a placeholder that redirects to either
   the Linked List implementation or the Array implementation,
   as may be chosen in mondriaan.mk */

#if defined(GAINBUCKET_LIST)
#include "GainBucketList.h"
#elif defined(GAINBUCKET_ARRAY)
#include "GainBucketArray.h"
#else
#error "A gain bucket structure should be selected in mondriaan.mk"
#endif
