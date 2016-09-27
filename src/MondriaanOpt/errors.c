#include <stdio.h>
#include <stdlib.h>

void exitwitherror(unsigned int err){
    /* Exit the program immediately with error code 'err' */
    
    if(err==0){
        /* Options error, output program usage */
        fprintf(stderr,"MondriaanOpt - written by Daan Pelt and Rob Bisseling (2015)\n");
        fprintf(stderr,"\nInvalid options are given...\nUsage:\n\n");
        fprintf(stderr,"  ./program -m \"matrix file\" -e \"load imbalance\" -v \"volume\"\n");
        fprintf(stderr," or \n");
        fprintf(stderr,"  ./program -m \"matrix file\" -k \"number of nonzeros\" -v \"volume\"\n\n");
        fprintf(stderr,"Use -h option for more help\n");
        fprintf(stderr,"\n");
        exit(EXIT_SUCCESS);
    }else if(err==1){
        /* Help asked */
        fprintf(stderr,"MondriaanOpt - written by Daan Pelt and Rob Bisseling (2015)\n");
        printf("\nProgram to find optimal matrix bipartitioning.\n\nRequired options:\n");
        printf("  -m \"matrix file\"        : the input matrix file in Matrix Market (.mtx) format\n");
        printf("  -v \"volume\"             : the starting upper bound volume\n");
        printf("\nAlso, it is required to pass either of the following:\n");
        printf("  -e \"load imbalance\"     : the allowed load imbalance\n");
        printf("  -k \"number of nonzeros\" : the maximum allowed number of nonzeros per part\n");
        printf("\nFurther options:\n");
        printf("  -t \"seconds\"            : max running time in seconds\n");
        printf("  -r \"dumpfile\"           : resume with given dumpfile\n");
        printf("  -h                      : show this help\n");
        printf("\n");
        exit(EXIT_SUCCESS);
    }else if(err==2){
        /* Not enough memory */
        fprintf(stderr,"MondriaanOpt - written by Daan Pelt and Rob Bisseling (2015)\n");
        fprintf(stderr,"\nNot enough memory\n");
        fprintf(stderr,"\n");
        exit(EXIT_SUCCESS);
    }
} /* end exitwitherror */
