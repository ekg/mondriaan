#include <stdlib.h>
#include <sys/time.h>

#include "options.h"
#include "matrix.h"
#include "solution.h"
#include "solver.h"

int main(int argc,char **argv){

    /* Solution, options, matrix structure */
    struct solution sol;
    struct options opt;
    struct mat a;

    /* Read options from command line */
    readoptions(&opt,argc,argv);

    /* Read matrix from the file specified on command line */
    readmatrixfromfile(&a,opt.fn);
    
    /* Initialise solution */
    initsolution(&a,&sol,&opt);

    /* Get starting time */
    gettimeofday(&(opt.start),NULL);

    /* Start branch and bound algorithm */
    solve(&a,&sol);

    /* Get ending time */
    gettimeofday(&(opt.end),NULL);

    /* Output final number of branches considered */
    printf("DONE!\nNUMBER OF BRANCHES CONSIDERED: %lu\n",(unsigned long)sol.nbranches);

    /* Get time used */
    opt.time += (double)(opt.end.tv_sec - opt.start.tv_sec) + (double)(opt.end.tv_usec - opt.start.tv_usec)/1000000.;

    /* Output time used */
    printf("Time taken: %lf s\n",opt.time);
    
    if(sol.maxvol < opt.maxvol) {
        /* A solution has been found */
        /* Output final volume upper bound to stdout */
        printf("Final volume: %d\n",sol.maxvol);
        
        /* Convert to other file formats */
        printConverted(&opt);
    }
    else {
        /* No solution found. */
        /* Subtract 1 from maxvol, as maxvol contains 'the volume we want to improve', while we want to print 'the upper bound'. */
        printf("No solution with a volume at most %d exists!\n", opt.maxvol-1);
    }
    
    exit (EXIT_SUCCESS);
}
