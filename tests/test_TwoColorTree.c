#include "SubsetSum.h"

void test_TwoColorTree();

int main(int argc, char **argv) {

    printf("Test TwoColorTree: ");
    
    SetRandomSeed(111);
    /*struct timeval tv;
    gettimeofday(&tv,NULL);
    unsigned long time_in_micros = 1000000 * tv.tv_sec + tv.tv_usec;
    SetRandomSeed(time_in_micros);*/
    
    test_TwoColorTree();

    printf("OK\n");
    exit(0);

} /* end main */


/**
 * Test TwoColorTree()
 */
void test_TwoColorTree() {
    long u, v, i, N = 100;
    
    struct Graph Tree;
    GraphInit(&Tree, N);
    
    for(v=1; v<N; ++v) {
        u = Random1(0, v-1);
        GraphAddEdge(&Tree, u, v);
    }
    
    TwoColorTree(&Tree);
    
    struct Node U;
    for(u=0; u<N; ++u) {
        U = Tree.nodes[u];
        
        if(U.val != 0 && U.val != 1) {
            printf("Error\n");
            exit(1);
        }
        
        for(i=0; i<U.deg; ++i) {
            if(Tree.nodes[U.adjList[i]].val == U.val) {
                printf("Error\n");
                exit(1);
            }
        }
    }
    
    GraphDestroy(&Tree);
    
} /* end test_TwoColorTree */
