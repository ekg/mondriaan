#ifndef OPTIONS_H_DONE

#define OPTIONS_H_DONE

#include <sys/time.h>

#define MAXFNSIZE 64
#define CUR_REQ_OPTIONS 3
#define TRUE 1
#define FALSE 0

/* period of printing current depth and checking whether to stop */
#define PERIOD 8388608 /* 2^23 */

struct options {
    /* Structure holding the input options */

    /* Input filename */
    char fn[MAXFNSIZE];

    /* Load imbalance: epsset = TRUE means it is set as a fraction eps,
       epsset = FALSE means it is set as a maximum number of nonzeros k */
    double eps;
    unsigned int k;
    unsigned char epsset;
    struct timeval start, end;

    unsigned long long int nbranches;
    unsigned int maxvol;

    unsigned char resume;
    char resfn[MAXFNSIZE];

    double time;

    double maxruntime;
    
    char convert; /* 0 = no convert, 1 = convert after calculation, 2 = only convert */

};

char readoptions(struct options *o, int c, char **v);

#endif
