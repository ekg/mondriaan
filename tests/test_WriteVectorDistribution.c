#include <stdlib.h>
#include <stdio.h>

#include "DistributeVecLib.h"

int main(int argc, char **argv) {

    FILE *fp; 
    char filename[MAX_WORD_LENGTH];
    /*char b1[MAX_WORD_LENGTH];
    char b2[MAX_WORD_LENGTH];
    char b3[MAX_WORD_LENGTH];
    char b4[MAX_WORD_LENGTH];
    char b5[MAX_WORD_LENGTH];
    char b6[MAX_WORD_LENGTH];*/
    long P, P1, n, n1, j, j1, s1;
    long int *owner;
    struct opts Options;
    SetDefaultOptions( &Options );

    printf("Test WriteVectorDistribution: ");
    P = 5;  /* number of  processors >= 1 */
    n = 100; /* length of the owner array */

    owner = (long int *)malloc(n*sizeof(long int));

    if (owner == NULL){
        printf("Error (allocation)\n");
        exit(1);
    }

    /* Initialise owners */
    for (j=0; j<n; j++) 
        owner[j]= (j+1)%P;

    /* Write the distribution to file */
    strcpy(filename,"outWriteVectorDistribution");
    
    fp = fopen(filename, "w");
    WriteVectorDistribution(owner, NULL, n, P, fp, &Options);
    fclose(fp);

    /* Read the values back and check whether they are OK */
    fp = fopen(filename, "r");
    
    if (!fp) {
        printf("Error (opening output)\n");
        exit(1);
    }
   
    /* Below only for Extended Matrix-Market format, DMM uses no header at all
    if (fscanf(fp,"%s %s %s %s %s %s",b1,b2,b3,b4,b5,b6) != 6) {
        printf("Error (reading banner)\n");
        exit(1);
    }

    if (strcmp(b1,"%%Extended-MatrixMarket") || strcmp(b2,"distributed-vector") 
        || strcmp(b4,"integer") || strcmp(b3,"array") || strcmp(b5,"general") 
        || strcmp(b6,"global")) {
        printf("Error (banner)\n");
        exit(1);
    }*/
 
    if (fscanf(fp,"%ld %ld",&n1,&P1) != 2) {
        printf("Error (reading header)\n");
        exit(1);
    }
    
    if (n1 != n || P1 != P) {
        printf("Error (header)\n");
        exit(1);
    }

    for (j=0; j<n; j++){ 
        if (fscanf(fp,"%ld %ld",&j1,&s1) != 2) {
            printf("Error (data field)\n");
            exit(1);
        }
        /* Values read from file are 1-based */
        if (j1 != j+1 || s1 != ((j+1)%P) + 1) {
            printf("Error (data)\n"); 
            exit(1); 
        } 
    }
    
    fclose(fp);

    printf("OK\n");
    exit(0);

} /* end main */
