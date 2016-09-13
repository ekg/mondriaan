#include "Graph.h"
#include "HKLFM.h"
#include "Match.h"

struct opts Options ;

extern int ClearMoveLog(struct biparthypergraph *pHG);

int main(int argc, char **argv) {

    struct biparthypergraph HG ;

    printf("Test ClearMoveLog: ");

    HG.CurVtxLog = 1 ;
    HG.MinVtxLog = 1 ;
    ClearMoveLog(&HG) ;

    /* Check Log values */
    if (HG.CurVtxLog != 0 || HG.MinVtxLog != 0 ) {
        printf("Error\n") ;
        exit(1);
    }

    printf("OK\n") ;
    exit(0);

} /* end main */
