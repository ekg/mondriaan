#include "GainBucket.h"

#if defined(GAINBUCKET_LIST)
#include "GainBucketList.c"
#elif defined(GAINBUCKET_ARRAY)
#include "GainBucketArray.c"
#else
#error "A gain bucket structure should be selected in mondriaan.mk"
#endif
