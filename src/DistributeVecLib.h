
#ifndef __DistributeVecLib_h__
#define __DistributeVecLib_h__

#include "Options.h" 
#include "SparseMatrix.h" 

int AssignColumnToProc(long int *X, long *procstart, int *procindex,
                        long *Ns, long *Nr, long *Nv, long *Sums,
                        long j, int q);
int RemoveColumnFromProc(long int *X, long *procstart, int *procindex,
                          long *Ns, long *Nr, long *Nv, long *Sums, 
                          long j, int q);
int AssignRemainingColumns(long l, int P, long int *X, long *Nv);
int AssignRemainingNonemptyColumns(long l, int P, long int *X, 
                                    long *procstart, int *procindex,
                                    long *Ns, long *Nr, long *Nv);
  
int InitNprocs(const struct sparsematrix *pM, int dir, int *Nprocs);
int InitProcindex(const struct sparsematrix *pM, int dir, int *Nprocs,
                   long *procstart, int *procindex);

int GenerateHistogram(int *X, long l, int lo, int hi, long *Histogram );
int PrintHistogram(int lo, int hi, long *Histogram);

int CalcCom(const struct sparsematrix *pM, long int *X, int dir, long *ComVol, 
             long *MaxOut, long *MaxIn, long *MaxCompnts, long *TotCompnts);
int PrintCom(int P, long l, int dir, long ComVol, long MaxOut, long MaxIn,
              long MaxCompnts, long TotCompnts);
int CalcLocalLowerBound(const struct sparsematrix *pM, int dir, long *LB, int *Pact);

void PrintVecStatistics(int P, long *Ns, long *Nr, long *Nv);
int WriteVector(const long int *X, const char base, const char *name, long l, int P, FILE *fp, const struct opts *pOptions);
int WriteVectorDistribution(const long int *X, const char *name, long l, int P, FILE *fp, const struct opts *pOptions);
int WriteVectorCollection(long int **X, const char *name, const long i, const long *j, FILE *fp);
 
#endif /* __DistributeVecLib_h__ */
