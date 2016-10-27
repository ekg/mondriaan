#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "options.h"
#include "errors.h"

char readoptions(struct options *o, int c, char **v){
    /* Reads options to 'o' structure from the command line,
       where c=argc and v=argv.
       The function returns FALSE if there is a problem, 
       and TRUE if everything is ok */

    /* Variables to check whether required options are given */
    char req[CUR_REQ_OPTIONS];

    int i;

    /* Initialise required options */
    for(i=0; i<CUR_REQ_OPTIONS; i++)
        req[i]=FALSE;

    /* First option to read */
    i=1;

    o->resume=FALSE;
    o->time=0.;
    o->maxruntime=0.;
    o->nbranches=0; /* number of branches traversed */
    o->SVG=SVGNo;
    o->maxvol = -1;
    
    while(i<c && v[i][0] != '-') {
        if(i == 1) {
            /* Read matrix filename */
            sprintf(o->fn,"%s",v[i]);
            req[0]=TRUE;
        }
        if(i == 2) {
            if(strcmp(v[i], "2") != 0) {
                exitwitherror(0);
            }
        }
        if(i == 3) {
            o->eps = atof(v[i]);
            o->epsset = TRUE;
            req[1]=TRUE;
        }
        
        i++;
    }
    

    /* Read while there are options left */
    while(i<c){
        if(strcmp(v[i],"-t")==0){
            /* Max running time in seconds */
            if(i==c-1) exitwitherror(0);
            o->maxruntime = atof(v[i+1]);
            i++;
        }else if(strcmp(v[i],"-e")==0){
            /* Load imbalance */
            if(i==c-1) exitwitherror(0);
            o->eps = atof(v[i+1]);
            o->epsset = TRUE;
            req[1]=TRUE;
            i++;
        }else if(strcmp(v[i],"-k")==0){
            /* Number of nonzeros */
            if(i==c-1) exitwitherror(0);
            o->k = atoi(v[i+1]);
            o->epsset = FALSE;
            req[1]=TRUE;
            i++;
        }else if(strcmp(v[i],"-v")==0){
            /* Starting upper bound on volume, e.g. MIN(m,n)+1 for an m by n matrix  */
            if(i==c-1) exitwitherror(0);
            /* We add +1, to change from 'upper bound' to 'the value we want to improve upon' */
            o->maxvol = atoi(v[i+1])+1;
            i++;
        }else if(strcmp(v[i],"-h")==0){
            /* Get help */
            exitwitherror(1);
        }else if(strcmp(v[i],"-r")==0){
            /* Resume with dumpfile of solution */
            if(i==c-1) exitwitherror(0);
            o->resume=TRUE;
            sprintf(o->resfn,"%s",v[i+1]);
            i++;
        }else if(strcmp(v[i],"-svg")==0){
            o->SVG=SVGYes;
        }
        i++;
    }
    
    /* Check if required options are given */
    for(i=0; i<CUR_REQ_OPTIONS; i++){
        if(!req[i])
            exitwitherror(0);
    }
    
    /* Everything ok */
    return TRUE;

} /* end readoptions */
