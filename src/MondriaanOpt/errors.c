#include <stdio.h>
#include <stdlib.h>

void exitwitherror(unsigned int err){
    /* Exit the program immediately with error code 'err' */
    
    if(err==0){
        /* Options error, output program usage */
        fprintf(stderr,"MondriaanOpt - written by Daan Pelt and Rob Bisseling (2015)\n");
        fprintf(stderr,"\nInvalid options are given...\nUsage:\n");
        fprintf(stderr,"  ./tools/MondriaanOpt matrix [P [eps]] [options]\n\n");
        
        fprintf(stderr,"Either [eps], -e or -k must be passed, and it is advised to pass -v.\n");
        fprintf(stderr,"Use -h option for more help\n");
        fprintf(stderr,"\n");
        exit(EXIT_SUCCESS);
    }else if(err==1){
        /* Help asked */
        fprintf(stderr,"MondriaanOpt - written by Daan Pelt and Rob Bisseling (2015)\n");
        printf("\nProgram to find optimal matrix bipartitioning.\n\nUsage:\n");
        printf("  ./tools/MondriaanOpt matrix [P [eps]] [options]\n\n");
        
        printf("One, two or three parameters may be passed, after which further options may be given.\n");
        printf("Either [eps], -e or -k must be passed, and it is advised to pass -v (see below).\n\n");
        
        printf("Parameters:\n");
        printf("  matrix                  : the input matrix file in Matrix Market (.mtx) format (may not start with a dash (-))\n");
        printf("  P                       : the number of processors. This must equal 2; this option is present for consistency with the other Mondriaan* commands\n");
        printf("  eps                     : the maximum allowed load imbalance\n\n");
        
        printf("Further options:\n");
        printf("  -v \"volume\"             : the starting upper bound volume\n");
        printf("  -e \"load imbalance\"     : the maximum allowed load imbalance\n");
        printf("  -k \"number of nonzeros\" : the maximum allowed number of nonzeros per part\n");
        printf("  -t \"seconds\"            : max running time in seconds\n");
        printf("  -h                      : show this help\n");
        printf("  -svg                    : Write visualisations of the partitioning to .svg files\n\n");
        
        printf("Apart from the matrix, at least one of [eps], -e or -k must be given, defining the maximum allowed load imbalance.\n");
        printf("The default value for the initial upper bound on the communication volume is m+n (m and n being the dimensions of the matrix),\n");
        printf("but it is strongly recommended to pass a better upper bound (-v) if available, to reduce computing time.\n\n");
        
        printf("Equivalent examples:\n");
        printf("  ./tools/MondriaanOpt tests/arc130.mtx 2 0.03 -v 17\n");
        printf("  ./tools/MondriaanOpt tests/arc130.mtx -e 0.03 -v 17\n");
        printf("  ./tools/MondriaanOpt tests/arc130.mtx -k 660 -v 17\n");
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
