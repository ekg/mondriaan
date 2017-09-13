#include "SubsetSum.h"

int test_SubsetSum(int KK_diff, double eps1, double eps2, long minN, long maxN);
void time_comparison();

int main(int argc, char **argv) {

    printf("Test SubsetSum: ");
    
    SetRandomSeed(111);
    /*struct timeval tv;
    gettimeofday(&tv,NULL);
    unsigned long time_in_micros = 1000000 * tv.tv_sec + tv.tv_usec;
    SetRandomSeed(time_in_micros);*/
    
    double epses[3] = {0.0, 0.03, 0.3};
    int i, e1, e2;
    
    for(i=0; i<2; ++i) {
        for(e1=0; e1<3; ++e1) {
            for(e2=0; e2<3; ++e2) {
                if(!test_SubsetSum(i, epses[e1], epses[e2], 2, 16)) {
                    printf("Error\n");
                    exit(1);
                }
            }
        }
    }
    
    /*time_comparison();*/
    
    printf("OK\n");
    exit(0);

} /* end main */


/**
 * Test SubsetSum() or KarmarkarKarp()
 * 
 * Input:
 * KK_diff           : Whether the Karmarkar-Karp differencing algorithm or naive enumeration should be used
 * eps1              : Maximum load imbalance of part 1
 * eps2              : Maximum load imbalance of part 2
 * minN              : Minimum N
 * maxN              : Maximum N
 */
int test_SubsetSum(int KK_diff, double eps1, double eps2, long minN, long maxN) {
    
    long i;
    
    long N = Random1(minN,maxN);
    long n1 = Random1(1, N-1);
    long n2 = N-n1;
    
    long w1 = 0, w2 = 0;
    long *weights = (long *)malloc(N*sizeof(long));
    
    /* Fill weights of subset 1 */
    for(i=0; i<n1; ++i) {
        weights[i] = Random1(1, 50);
        w1 += weights[i];
    }
    
    /* Fill weights of subset 2 */
    for(i=0; i<n2; ++i) {
        weights[n1+i] = Random1(1, 50);
        w2 += weights[n1+i];
    }
    
    RandomPermute(weights, NULL, NULL, NULL, 0, N-1);
    
    long *permutation = NULL;
    long *selected = NULL;
    
    long w_tot = w1+w2;
    long w_small = w1<w2 ? w1 : w2;
    long w_large = w1>w2 ? w1 : w2;
    long w_lo = w_small * (1+eps1);
    long w_hi = w_large * (1+eps2);
    long w_max = w_lo;
    long w_min = w_tot - w_hi;
    int found;
    if(KK_diff)
        found = KarmarkarKarp(weights, N, w_lo, w_hi, &permutation, &selected);
    else
        found = SubsetSumExp(weights, N, w_lo, w_hi, &permutation, &selected);
    
    /* By construction, a subset exists */
    if(!found) {
        if(KK_diff) {
            /* However, KK may not be able to find it */
            free(weights);
            return TRUE;
        }
        return FALSE;
    }
    
    /* Check weights of the found subsets */
    long sum[2] = {0, 0};
    for(i=0; i<N; ++i) {
        sum[selected[i]] += weights[i];
    }
    
    if(!(w_min <= sum[1] && sum[1] <= w_max) || !(w_tot - w_max <= sum[0] && sum[0] <= w_tot - w_min)) {
        return FALSE;
    }
    
    free(permutation);
    free(selected);
    free(weights);
    
    return TRUE;
    
} /* end test_SubsetSum */

/*
 * Benchmark utility to compare the SubsetSum and KarmarkarKarp methods.
 */
void time_comparison() {
    long i, N, Nmax = 500, numRuns = 100;
    
    clock_t starttime, endtime;
    double cputime;
    
    for(N=2; N<=Nmax; ++N) {
        
        printf("%ld ", N);
        
        starttime = clock();
        for(i=0; i<numRuns; ++i) {
            test_SubsetSum(0, 0.0, 0.0, N, N);
        }
        endtime = clock();
        cputime = ((double) (endtime - starttime)) / CLOCKS_PER_SEC;
        printf("%lf ", cputime/numRuns);
        
        
        starttime = clock();
        for(i=0; i<numRuns; ++i) {
            test_SubsetSum(1, 0.0, 0.0, N, N);
        }
        endtime = clock();
        cputime = ((double) (endtime - starttime)) / CLOCKS_PER_SEC;
        printf("%lf\n", cputime/numRuns);
        
    }
    
    
}
