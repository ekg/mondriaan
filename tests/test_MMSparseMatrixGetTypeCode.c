#include "SparseMatrix.h"

int main(int argc, char **argv) {

    char line[MAX_LINE_LENGTH] ;
    struct sparsematrix A ;
    struct opts Options;
    SetDefaultOptions( &Options );

    printf("Test MMSparseMatrixGetTypeCode: ");
    
    A.MMTypeCode[0] = 'D' ; /* distributed matrix */
    A.MMTypeCode[1] = 'A' ; /* dense matrix (array) */
    A.MMTypeCode[2] = 'C' ; /* complex matrix */
    A.MMTypeCode[3] = 'K' ; /* skew-symmetric matrix */
    A.ViewType      = ViewTypeOriginal ; /* original view */

    MMSparseMatrixGetTypeCode(&A, line, NULL, &Options) ;

    /* Check banner line and 5 strings */
    if (strcmp(line,"%%MatrixMarket distributed-matrix array complex skew-symmetric\n") != 0 )
        printf("Literal error\n") ;
    if (strcmp(A.Banner, MM_Banner) != 0 ||
        strcmp(A.Object,"distributed-matrix") != 0 ||
        strcmp(A.Format,"array") != 0 ||
        strcmp(A.Field,"complex") != 0 ||
        strcmp(A.Symmetry,"skew-symmetric") != 0 ){
        printf("Typecode error\n") ;
        exit(1);
    }

    printf("OK\n") ;
    exit(0);

} /* end main */
