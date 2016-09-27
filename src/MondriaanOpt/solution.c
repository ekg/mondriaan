#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <libgen.h>

#include "options.h"
#include "solution.h"
#include "matrix.h"


void initsolution(const struct mat *a, struct solution *s, struct options *o){

    /* This function initialises the solution s for use with matrix a,
       and allocates the memory for the arrays of s */

    /* Allocate arrays */
    s->vals[ROW] = malloc(sizeof(unsigned char) * a->m);
    s->vals[COL] = malloc(sizeof(unsigned char) * a->n);

    s->n[TO0][ROW] = malloc(sizeof(unsigned int) * a->m);
    s->n[TO0][COL] = malloc(sizeof(unsigned int) * a->n);
    s->n[TO1][ROW] = malloc(sizeof(unsigned int) * a->m);
    s->n[TO1][COL] = malloc(sizeof(unsigned int) * a->n);

    s->startUnAssigned[ROW] = malloc(sizeof(unsigned int)*a->m);
    s->startUnAssigned[COL] = malloc(sizeof(unsigned int)*a->n);

    s->partialstate[ROW] = malloc(sizeof(unsigned char)*a->m);
    s->partialstate[COL] = malloc(sizeof(unsigned char)*a->n);

    s->partial[TO0] = malloc(sizeof(unsigned int)*a->cmax);
    s->partial[TO1] = malloc(sizeof(unsigned int)*a->cmax);

    s->justpartial = malloc(sizeof(unsigned int)*a->cmax);
    s->justpartialoff = malloc(sizeof(unsigned int)*a->cmax);

    s->prior = malloc(sizeof(unsigned int)*(a->m+a->n));

    s->nexts[ROW] = malloc(sizeof(struct priority)*(a->cmax+1));
    s->nexts[COL] = malloc(sizeof(struct priority)*(a->cmax+1));

    s->idxi[ROW] = malloc(sizeof(struct indexedint) * a->m);
    s->idxi[COL] = malloc(sizeof(struct indexedint) * a->n);

    /* Initialise graph structure */
    s->mGraph.match = malloc(sizeof(int)*(a->m+a->n));
    s->mGraph.nmatch = 0;
    s->mGraph.adj = malloc(sizeof(int)*(2*a->nnz));
    s->mGraph.start = malloc(sizeof(int)*(a->m+a->n));
    s->mGraph.length = malloc(sizeof(int)*(a->m+a->n));
    s->mGraph.visited = malloc(sizeof(int)*(a->m+a->n));
    s->mGraph.vertex1 = malloc(sizeof(int)*(a->m+a->n));
    s->mGraph.pred = malloc(sizeof(int)*(a->m+a->n));
    s->mGraph.vertex2 = malloc(sizeof(int)*(a->m+a->n));
    s->mGraph.dstart = malloc(sizeof(int)*(a->m+a->n + 1));
    s->mGraph.inGraph = malloc(sizeof(unsigned char)*(a->m+a->n));

    /* Set maximum number of nonzeros that can be assigned to a single part */
    if(o->epsset){
        if ((a->nnz)%2 == 0)
            o->k = (unsigned int)((1.+o->eps)*(a->nnz/2));
        else
            o->k = (unsigned int)((1.+o->eps)*((a->nnz + 1)/2));
    }
    s->max1 = o->k;

    /* Set volume upper bound */
    s->maxvol = o->maxvol;

    /* Set matrix name */
    s->matname = basename(o->fn);

    /* Set start time of solution process */
    time(&(s->starttime));

    /* Set the maximum runtime */
    s->maxruntime = o->maxruntime;

    return;

} /* end initsolution */


void insert_priority(struct solution *s, unsigned char dir, unsigned int c, unsigned int idx){

    /* This function inserts a new entry with index idx into the linked list
       for rows of length c (if dir = ROW) or columns of length c (if dir = COL).
       The entry is inserted at the header of the list. */

    struct priority *new = malloc(sizeof(struct priority));

    new->idx = idx;
    new->next = s->nexts[dir][c].next;
    s->nexts[dir][c].next = new;

} /* end insert_priority */


void delete_priority(struct solution *s, unsigned char dir, unsigned int c, unsigned int idx){

    /* This function searches for an entry with index idx in the linked list  
       for rows of length c (if dir = ROW) or columns of length c (if dir = COL).
       If such an entry exists in the list, it is deleted and its memory freed. */

    struct priority *prev = &(s->nexts[dir][c]);
    struct priority *cur = prev->next;

    while(cur!=NULL){
        /* end of list not yet reached */
        if(cur->idx==idx){
            /* entry found */
            prev->next = cur->next;
            free(cur);
            break;
        }
        prev = cur;
        cur = prev->next;
    }

} /* delete_priority */


void resetsolution(const struct mat *a, struct solution *s){

    /* This function resets the solution s for use with matrix a,
       setting all its variables to 0 or another default. */

    int i, j, c;

    /* Current volume */
    s->vol = 0;

    /* Number of nonzeros just partially assigned or unassigned */
    s->njustpartial = 0;
    s->njustpartialoff = 0;

    /* Initialise row information */
    for(i=0; i<a->m; i++){
        s->vals[ROW][i] = UNAS;
        for(j=0; j<2; j++)
            s->n[j][ROW][i] = 0; /* row i, part j */
        s->startUnAssigned[ROW][i] = 0;
        s->partialstate[ROW][i] = UNAS;
    }

    /* Initialise column information */
    for(i=0; i<a->n; i++){
        s->vals[COL][i] = UNAS;
        for(j=0; j<2; j++) 
            s->n[j][COL][i] = 0; /* column i assigned to part j */
        s->startUnAssigned[COL][i] = 0;
        s->partialstate[COL][i] = UNAS;
    }

    for(c=0; c<a->cmax; c++){
        s->partial[TO0][c] = 0;
        s->partial[TO1][c] = 0;
    }

    s->npartial[TO0] = 0;
    s->npartial[TO1] = 0;

    s->nprior = 0;

    s->tot[TO0] = 0;
    s->tot[TO1] = 0;
    s->depth = 0;
    s->nbranches = 0;
    s->prevdir = COL;
    s->FirstIdxAssignedTo0 = FALSE;

    /* Initialise linked list headers */
    for(c=0; c<=a->cmax; c++){
        s->nexts[ROW][c].next = NULL;
        s->nexts[COL][c].next = NULL;
    }

    /* Insert each row i in the linked list of its length */
    for(i=0; i<a->m; i++){
        s->idxi[ROW][i].val = i;
        c = a->starts[ROW][i+1] - a->starts[ROW][i];
        s->idxi[ROW][i].key = c ;
        insert_priority(s,ROW,c,i);
    }

    /* Insert columns in the linked list of their length */
    for(i=0; i<a->n; i++){
        s->idxi[COL][i].val = i+a->m; /* number the rows and columns with index 0 <= idx < m+n */
        c = a->starts[COL][i+1] - a->starts[COL][i];
        s->idxi[COL][i].key = c;
        insert_priority(s,COL,c,i+a->m);
    }

    /* Initialise graph arrays */
    for(i=0; i<a->m+a->n; i++){
        s->mGraph.match[i] = i;
        s->mGraph.length[i] = 0;
        s->mGraph.visited[i] = 0;
        s->mGraph.vertex1[i] = 0;
        s->mGraph.vertex2[i] = 0;
        s->mGraph.pred[i] = 0;
        s->mGraph.inGraph[i] = 0;
     }
    for(i=0; i<=a->m+a->n; i++)
        s->mGraph.dstart[i] = 0; /* note: one item extra */
    for(i=0; i<a->m; i++)
        s->mGraph.start[i] = a->starts[ROW][i];
    for(i=0; i<a->n; i++)
        s->mGraph.start[i+a->m] = a->starts[COL][i] + a->nnz;

    return;

} /* end resetsolution */


void printToSVG(const char *fn, struct solution *s, const struct mat *a){

    /* This function plots the bipartitioning of sparse matrix a
       computed as solution s. The output is to a file with filename fn.
       The format is Scalable Vector Graphics (SVG).

       The colour of a nonzero is:
          red if it is assigned to part 0,
          blue if it is assigned to part 1,
          yellow if it can be assigned to both, since both its row and its column are cut. */

    
    FILE *fp;
    int i, j, k;

    /* Open file */
    fp = fopen(fn,"w");
    if(fp==NULL)
        return ;

    /* Output standard SVG file intro */
    fprintf(fp,"<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\"?>\n");
    fprintf(fp,"<svg xmlns:svg=\"http://www.w3.org/2000/svg\" xmlns=\"http://www.w3.org/2000/svg\" version=\"1.0\"\n");

    /* Plot has width 10*n and height 10*m, where m = number of rows and n = number of columns */
    fprintf(fp,"width=\"%d\"\nheight=\"%d\"\nid=\"svg2\">",10*a->n,10*a->m);
    fprintf(fp,"<g>\n");

    /* Plot white background box */
    fprintf(fp,"<rect width=\"%d\" height=\"%d\" x=\"0\" y=\"0\" id=\"bkg\" style=\"fill:#ffffff;fill-opacity:1;\" />\n",10*a->n,10*a->m);

    for(i=0; i<a->m; i++){
        /* Plot the nonzeros of row i */
        for(k=a->starts[ROW][i]; k<a->starts[ROW][i+1]; k++){
            j = a->indices[ROW][k];
            /* Plot nonzero a(i,j) as a 10 by 10 box at position
                   x = 10*j (from the left margin),
                   y = 10*i (from the top margin) */
            if(s->vals[ROW][i]==TO0 || s->vals[COL][j]==TO0){
                /* Plot red nonzero */
                fprintf(fp,"<rect width=\"10\" height=\"10\" x=\"%d\" y=\"%d\" id=\"r%d-%d\" style=\"fill:#ff0000;fill-opacity:1;\" />\n",10*j,10*i,i,j);
            }else if(s->vals[ROW][i]==TO1 || s->vals[COL][j]==TO1){
                /* Plot blue nonzero */
                fprintf(fp,"<rect width=\"10\" height=\"10\" x=\"%d\" y=\"%d\" id=\"r%d-%d\" style=\"fill:#0000ff;fill-opacity:1;\" />\n",10*j,10*i,i,j);
            }else if(s->vals[ROW][i]==CUT && s->vals[COL][j]==CUT){
                /* Plot yellow nonzero */
                fprintf(fp,"<rect width=\"10\" height=\"10\" x=\"%d\" y=\"%d\" id=\"r%d-%d\" style=\"fill:#ffff00;fill-opacity:1;\" />\n",10*j,10*i,i,j);
            }else{
                /* Plot black nonzero. This should not happen. */
                fprintf(fp,"<rect width=\"10\" height=\"10\" x=\"%d\" y=\"%d\" id=\"r%d-%d\" style=\"fill:#000000;fill-opacity:1;\" />\n",10*j,10*i,i,j);
            }
        }
    }

    fprintf(fp,"</g>\n</svg>\n");

    /* Close file */
    fclose(fp);

} /* end printToSVG */


void printToMM(const char *fn, struct solution *s, const struct mat *a){

    /* This function outputs the bipartitioning of sparse matrix a
       computed as solution s. The output is to a file with filename fn. 
       The format is Matrix Market, representing a sparse integer matrix 
       with nonzeros a(i,j), where the indices are in the range 1 <= i <= m 
       and 1 <= j <= n.
     
       The integer value of a nonzero is:
          0 if it is assigned to part 0,
          1 if it is assigned to part 1,
          2 if it can be assigned to both, since both its row and its column are cut. */

    FILE *fp;
    int i, j, k;

    /* Open file */
    fp = fopen(fn,"w");
    if(fp==NULL)
        return ;

    /* Print the Matrix Market header */
    fprintf(fp,"%%%%MatrixMarket matrix coordinate integer general\n");
    
    /* Print the number of rows, columns, nonzeros */
    fprintf(fp,"%d %d %d\n", a->m, a->n, a->nnz);

    for(i=0; i<a->m; i++){
        /* Print the nonzeros of row i */
        for(k=a->starts[ROW][i]; k<a->starts[ROW][i+1]; k++){
            j = a->indices[ROW][k];
            if(s->vals[ROW][i]==TO0 || s->vals[COL][j]==TO0){
                /* Add 1 to the indices, because Matrix Market is 1-based */
                fprintf(fp,"%d %d %d\n", i+1, j+1, 0);
            }else if(s->vals[ROW][i]==TO1 || s->vals[COL][j]==TO1){
                fprintf(fp,"%d %d %d\n", i+1, j+1, 1);
            }else if(s->vals[ROW][i]==CUT && s->vals[COL][j]==CUT){
                fprintf(fp,"%d %d %d\n", i+1, j+1, 2);
            }
        }
    }

    /* Close file */
    fclose(fp);

    return ;

} /* end printToMM */

void printConverted(const char *fn) {
    
    struct opts Options; /* The Mondriaan options */
    struct sparsematrix A;
    FILE *File = NULL;
    
    SetDefaultOptions(&Options);
    
    /* Remove ".mtx" from matrix name */
    char fn0[MAXFNSIZE], fnI[MAXFNSIZE];
    int len;
    len = strlen(fn);
    memcpy(fn0, fn, len-3);
    fn0[len-4] = 0;

    /* Construct processor number file name */
    sprintf(fnI,"%s_P2.mtx",fn0);
    
    /* Combine values & distribution A */
    if(!SpMatReadIndexAndValueMatrixFiles(fn, fnI, &A)) {
        exit(-1);
    }
    
    /* Write the distributed matrix to file */
    char output[MAX_WORD_LENGTH];
    sprintf(output, "%s-P%d", fn, A.NrProcs);
    File = fopen(output, "w");
    if (!File) fprintf(stderr, "printConverted(): Unable to open '%s' for writing!\n", output);
    else {
        MMWriteSparseMatrix(&A, File, NULL, &Options);
        fclose(File);
    }
    
    /* Write out matrix the entries of which are processor indices. */
    if (!MMInsertProcessorIndices(&A)) {
        fprintf(stderr, "printConverted(): Unable to write processor indices!\n");
        exit(-1);
    }

    A.MMTypeCode[0] = 'M';
    sprintf(output, "%s-I%d", fn, A.NrProcs);

    File = fopen(output, "w");

    if (!File) fprintf(stderr, "printConverted(): Unable to open '%s' for writing!\n", output);
    else {
        MMWriteSparseMatrix(&A, File, NULL, &Options);
        fclose(File);
    }
    A.MMTypeCode[0] = 'D';
    
    /* Write the index sets of the Cartesian submatrices to file */
    sprintf(output, "%s-C%d", fn, A.NrProcs);
    File = fopen(output, "w");
    if (!File) fprintf(stderr, "printConverted(): Unable to open '%s' for writing!\n", output);
    else {
        MMWriteCartesianSubmatrices(&A, File);
        fclose(File);
    }
    
}
