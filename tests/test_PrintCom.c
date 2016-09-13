#include "DistributeVecLib.h"

int main(int argc, char **argv) {

    long l, ComVol, MaxOut, MaxIn, MaxCompnts, TotCompnts;
    int P ;

    printf("Test PrintCom:\n");
    P = 8 ; /* number of  processors >= 1 */
    l = 400 ;
    ComVol = 1000 ;
    MaxOut = 250 ;
    MaxIn = 125 ;
    MaxCompnts = 100 ;
    TotCompnts = 400;

    PrintCom(P, l, ROW, ComVol, MaxOut, MaxIn, MaxCompnts, TotCompnts) ;

    printf("OK\n") ;
    exit(0);

} /* end main */
