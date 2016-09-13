#include "DistributeVecLib.h"

int main(int argc, char **argv) {

    long P, P2, k, n, i, j, *procstart, *Ns, *Nr, *Nv, *Sums;    
    int q, *procindex;
    long int *owner ;

    printf("Test RemoveColumnFromProc: ");
    P = 60 ; /* P is the number of processors, must be even, >=4 */
    P2 = P/2 ;
    k= 7 ;
    n= P*k ; /* P by n communication matrix C with processors
                in positions 0,1,...,P/2-1 in column j if j is even,
                and in positions P/2,...,P-1 in column j if j is odd. */

    Ns = (long *)malloc(P*sizeof(long));
    Nr = (long *)malloc(P*sizeof(long));
    Nv = (long *)malloc(P*sizeof(long));
    Sums = (long *)malloc(P*sizeof(long));
    procstart = (long *)malloc((n+1)*sizeof(long));
    procindex = (int *)malloc(n*P2*sizeof(int));
    owner = (long int *)malloc(n*sizeof(long int));

    if ( Ns == NULL || Nr  == NULL ||  Nv == NULL || Sums == NULL ||
         procstart == NULL || procindex == NULL || owner == NULL){
        printf("Error\n") ;
        exit(1);
    }

    /* Initialise procstart */
    for (j=0; j<n; j++){ 
        procstart[j]=P2*j ;
    } 
    procstart[n] = P2*n ;

    /* Fill procindex cyclically */
    for (i=0; i<P2*n; i++)
        procindex[i] = i%P ;

    /* Initialise owners, by the pattern
       0, P2, 1, P2+1, ... P2-1, P-1, 0, ... */
    for (j=0; j<n; j++)
        owner[j] = ((j/2)%P2) + (j%2)*P2 ;
 
    /* Initialise Ns, Nr, Nv */
    for (i=0; i<P; i++){  
        Ns[i] = n*(P2-1)/P;
        Nr[i] = Ns[i] ;
        Nv[i] = n/P ;
        Sums[i] = Ns[i] + Nr[i] ;
    }

    for (j=0; j<n; j++){  
        q= ((j/2)%P2) + (j%2)*P2 ;
        RemoveColumnFromProc(owner, procstart, procindex, Ns, Nr, Nv, Sums,
                             j, q);
    }  

    /* Check owners */
    for (j=0; j<n; j++){ 
        if (owner[j] != -1) {
            printf("Error\n") ;
            exit(1);
        }
    }
 
    /* Check number of sends, receives, components */
    for (i=0; i<P; i++){  
        if (Ns[i] != 0 || Nr[i] != 0 || Nv[i] != 0 || Sums[i] != n/2){  
            printf("Error\n") ;
            exit(1);
        } 
    }

    printf("OK\n") ;
    exit(0);

} /* end main */
