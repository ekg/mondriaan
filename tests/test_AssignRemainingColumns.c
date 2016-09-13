#include "DistributeVecLib.h"

int main(int argc, char **argv) {

    long P, k, n, i, j, proc, *Nv, *Nv1;    
    int q;
    long int *owner;

    printf("Test AssignRemainingColumns: ");
    P = 6 ; /* P is the number of processors, must be even */
    k= 10 ; /* k is the number of columns per processor,
               must be even */
    n= P*k ; /* n is the number of columns */

    Nv = (long *)malloc(P*sizeof(long));
    Nv1 = (long *)malloc(P*sizeof(long));
    owner = (long int *)malloc(n*sizeof(long int));

    if ( Nv == NULL || Nv1 == NULL || owner == NULL){
        printf("Error\n") ;
        exit(1);
    }

    /* Initialise Nv */
    for (i=0; i<P/2; i++)
        Nv[i] = k/2;
    for (i=P/2; i<P; i++)
        Nv[i] = 0;

    /* Initialise owner cyclically.
       Only first half of columns.  Only to procs 0..P/2-1.  */
    for (j=0; j<n; j++) 
        owner[j] = -1 ; 

    for (j=0; j<n/2; j++){ 
        proc = j%P ;
        if (proc < P/2)
            owner[j] = proc; 
    } 

    AssignRemainingColumns(n, P, owner, Nv);

    /* Compute Nv1 corresponding to owner */
    for (i=0; i<P; i++)
        Nv1[i] = 0;   

    for (j=0; j<n; j++){ 
        q=owner[j] ;
        if (q >= 0 && q < P) 
            Nv1[q]++ ;
    }

    /* Check that Nv[i] = Nv1[i] = k */
    for (i=0; i<P; i++) {
        if (Nv[i] != Nv1[i] || Nv[i] != k) {
            printf("Error\n") ;
            exit(1);
        }
    }
 
    printf("OK\n") ;
    exit(0);

} /* end main */
