#include "Sort.h"

int main(int argc, char **argv) {

    long random ;


    printf("Test SetRandomSeed: ");
    SetRandomSeed(-1) ;
    random = rand() ;
    if (random < 0 || random > (long)(RAND_MAX) ) {
        printf("Error\n") ;
        exit(1);
    } 
    printf("OK\n") ;

    exit(0);

} /* end main */
  
