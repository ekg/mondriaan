#include "DistributeVecLib.h"
#include "Sort.h"

int main(int argc, char **argv) {

    long P, n, i, j, s, nz, nztot,  *procstart, *Ns, *Nr, *Nv, *Nv1, *J;    
    int *procindex;
    long int *owner;

    printf("Test AssignRemainingNonemptyColumns: ");
    P = 17 ; /* P is the number of processors, must be odd, >= 3 */
    n = P ; /* P by P communication matrix C with processors
               in positions 0,1,...,P-1-j in column j for 0 <= j <= (P-1)/2 . */

    Ns = (long *)malloc(P*sizeof(long));
    Nr = (long *)malloc(P*sizeof(long));
    Nv = (long *)malloc(P*sizeof(long));
    Nv1 = (long *)malloc(P*sizeof(long));

    nz = ((P+1)*(3*P+1) ) / 8 ;
    procstart = (long *)malloc((n+1)*sizeof(long));
    procindex = (int *)malloc(nz*sizeof(int));
    owner = (long int *)malloc(n*sizeof(long int));

    if ( Ns == NULL || Nr  == NULL ||  Nv == NULL || Nv1 == NULL ||
         procstart == NULL || procindex == NULL || owner == NULL){
        printf("Error\n") ;
        exit(1);
    }

    /* Initialise Ns, Nr, Nv, Nv1 */
    for (i=0; i<P; i++){
        Ns[i] = 0;
        Nr[i] = 0;
        Nv[i] = 0;
        Nv1[i] = 0;
    }

    /* Initialise owner, procstart */
    nztot = 0;
    for (j=0; j <= (P-1)/2; j++){ 
        owner[j] = -1; 
        procstart[j] = nztot ;
        nztot += P-j ;
    } 
    for (j= (P+1)/2; j <n; j++){ 
        owner[j] = -1; 
        procstart[j] = nztot ;
    } 
    procstart[n] = nztot ;
    if ( nztot != nz) {
        printf("Error\n") ;
        exit(1);
    }


    /* Fill procindex. Only first (P+1)/2 columns. */
    for (j=0; j<= (P-1)/2 ; j++) 
            for (s=procstart[j]; s<procstart[j+1]; s++) 
                procindex[s] = s- procstart[j] ;

    AssignRemainingNonemptyColumns(n, P, owner, procstart, procindex, Ns, Nr, Nv) ;

    /* Check owners */
    /* First half should be owned */
    for (j=0; j <= (P-1)/2; j++){ 
        if (owner[j] < 0 || owner[j] >= P) {
            printf("Error\n") ;
            exit(1);
        }
    }
    /* Second half should be unowned */
    for (j= (P+1)/2; j<n; j++){ 
        if (owner[j] != -1) {
            printf("Error\n") ;
            exit(1);
        }
    } 

    /* Don't trust Nv. Count again. */
    for (j=0; j <= (P-1)/2; j++) 
        Nv1[owner[j]]++ ;
        
 
    /* Check maximum number of sends, receives, components.
       No processor should own more than 1 component */

    for (i=0; i<P; i++){  
        if (Ns[i] >= P || Nr[i] > (P+1)/2 || Nv[i] > 1 || Nv[i] != Nv[i] ){  
            printf("Error\n") ;
            exit(1);
        } 
    }
    
    /* Sort the Ns values in increasing order */
    J  = (long *) malloc(n * sizeof(long)) ;
    if (J == NULL){
        printf("Error\n") ;
        exit(1);
    }
    for (j=0; j < P  ; j++) 
        J[j] = j ;
 
    CSort(J, Ns, P-1, 0, P-1) ;
    
    /* The Ns values must be 0,0,...,0, (P-1)/2,...,P-1 */ 
    for (i=0; i<(P-1)/2;  i++){  
        if (Ns[J[i]] != 0 ){
            printf("Error\n") ;
            exit(1);
        } 
    }
    for (i=(P-1)/2; i < P; i++){  
        if (Ns[J[i]] != i ){
            printf("Error\n") ;
            exit(1);
        } 
    }
    
    printf("OK\n") ;
    exit(0);

} /* end main */
