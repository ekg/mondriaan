#ifndef __Sort_h__
#define __Sort_h__

#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include <time.h>

#ifdef UNIX
#include <sys/time.h>
#endif

/* Macros */
#define MIN(a, b)       ( (a) < (b) ? (a) : (b) )
#define MAX(a, b)       ( (a) > (b) ? (a) : (b) )

/* Logical constants */
#ifndef FALSE
#define FALSE 0
#endif
#ifndef TRUE
#define TRUE  1
#endif

/* Function declarations for Sort.c */
long Random1(long lo, long hi);
int SetRandomSeed(long seed);

void SwapLong(long *x, const long, const long);
void SwapDouble(double *x, long, long);

long *QSort(long *X, long LengthX);
void quicksort(long *item, long *Index, long lo, long hi);
int CSort(long *J, long *val, long maxval, long lo, long hi);

void RandomPermute(long *X, long *Y, double *D, double *E, long lo, long hi);

#endif /* __Sort_h__ */
