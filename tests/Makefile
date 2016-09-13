# Mondriaan Tests Makefile

include ../mondriaan.mk

#scrap default rule:
% :: %.c ;

#targets
TESTTARGETS:=test_QSort test_CSort test_quicksort \
test_RandomPermute test_Random1 test_SetRandomSeed test_SwapLong \
test_SwapDouble test_BucketInsert test_BucketMove \
test_BucketDeleteMax test_GainBucketGetMaxVal test_ClearGainBucket \
test_InitNprocs test_InitProcindex test_InitIntersection \
test_FindProcLowestSum test_AssignColumnToProc test_RemoveColumnFromProc \
test_AssignRemainingColumns test_AssignRemainingNonemptyColumns \
test_CalcCom test_CalcLocalLowerBound test_matallocl test_matfreel \
test_FirstNzInRow test_DistributeVecOpt2 test_DistributeVecLocal \
test_PrintVecLocalStatistics test_PrintVecStatistics test_PrintHistogram \
test_PrintCom test_WriteVectorDistribution \
test_DistributeVecGreedyImprove test_InitSums test_GenerateHistogram \
test_DistributeVecOrig test_DistributeVecOrigEq test_DistributeVec \
test_BalanceParts test_DetermineSplit test_logb2 \
test_ComputeWeight test_SplitMatrixSimple test_SplitMatrixKLFM \
test_DistributeMatrixMondriaan \
test_CreateNewBiPartHyperGraph test_DeleteBiPartHyperGraph \
test_SparseMatrix2BiPartHyperGraph test_BiPartHyperGraph2SparseMatrix \
test_SparseMatrixSymmetric2Full test_SparseMatrixSymmetricLower2Random \
test_SparseMatrixSymmetricRandom2Lower test_AddDummiesToSparseMatrix \
test_RemoveDummiesFromSparseMatrix test_MMSparseMatrixInit \
test_MMDeleteSparseMatrix test_MMSparseMatrixAllocateMemory \
test_MMSparseMatrixFreeMemory test_MMSparseMatrixReadEntries \
test_MMSparseMatrixReadPstart test_MMSparseMatrixReadWeights \
test_MMSparseMatrixReadTail test_MMSparseMatrixReadHeader \
test_MMSparseMatrixPrintHeader test_MMSparseMatrixPrintTail \
test_MMSparseMatrixPrintPstart test_MMSparseMatrixPrintWeights \
test_MMSparseMatrixPrintEntries test_MMWriteSparseMatrix \
test_MMSparseMatrixGetTypeCode test_MMSparseMatrixSetTypeCode \
test_MMReadSparseMatrix test_SparseMatrixRemoveDuplicates \
test_SparseMatrixStructurallySymmetric test_CRS2CCS \
test_MMWriteCartesianSubmatrices test_GetParameters test_SetDefaultOptions \
test_SetOption test_PrintHelp test_MoveVtxInNetAdjncy \
test_FindMatchArbitrary test_FindMatchInprod \
test_CreateInitialBalancedPartition test_ComputeInitialGains \
test_ClearMoveLog test_MoveVertex test_UpdateGains test_HKLFM test_RunHKLFM \
test_LambdaLambdaMinusOneMetric \
test_FindOptimalPathMatching \
test_MatchATA \
test_CoarsenGraph test_RunMLGraphPart

.SECONDARY: ${TESTTARGETS:%=%.o}

test: ../tools/Mondriaan ${TESTTARGETS}
	./runtest

all: test


build: ../tools/Mondriaan ${TESTTARGETS}
	

../tools/Mondriaan:
	

%.target: ${OBJDEPS}
	${CC} ${CFLAGS} -o ${@:%.target=%} ${OBJDEPS} ${LIBS}

%.o: %.c %.h
	${CC} ${CFLAGS} -c -o $@ ${filter %.c,$^}

%.o: %.c
	${CC} ${CFLAGS} -c -o $@ $^

test_MoveVtxInNetAdjncy: test_MoveVtxInNetAdjncy.o ${GRAPH:%=%.o} ${MATCH:%=%.o}
	make -r OBJDEPS='$^' $@.target

test_FindMatchArbitrary: test_FindMatchArbitrary.o ${GRAPH:%=%.o} ${MATCH:%=%.o}
	make -r OBJDEPS='$^' $@.target

test_FindMatchInprod: test_FindMatchInprod.o ${GRAPH:%=%.o} ${MATCH:%=%.o}
	make -r OBJDEPS='$^' $@.target

test_CreateInitialBalancedPartition: test_CreateInitialBalancedPartition.o ${GRAPH:%=%.o} ${MATCH:%=%.o} ${HKLFM:%=%.o}
	make -r OBJDEPS='$^' $@.target

test_ComputeInitialGains: test_ComputeInitialGains.o ${GRAPH:%=%.o} ${MATCH:%=%.o} ${HKLFM:%=%.o}
	make -r OBJDEPS='$^' $@.target

test_HKLFM: test_HKLFM.o ${GRAPH:%=%.o} ${MATCH:%=%.o} ${HKLFM:%=%.o}
	make -r OBJDEPS='$^' $@.target

test_RunHKLFM: test_RunHKLFM.o ${GRAPH:%=%.o} ${MATCH:%=%.o} ${HKLFM:%=%.o}
	make -r OBJDEPS='$^' $@.target

test_CoarsenGraph: test_CoarsenGraph.o ${GRAPH:%=%.o} ${MATCH:%=%.o} ${HKLFM:%=%.o}
	make -r OBJDEPS='$^' $@.target

test_RunMLGraphPart: test_RunMLGraphPart.o ${GRAPH:%=%.o} ${MATCH:%=%.o} ${HKLFM:%=%.o}
	make -r OBJDEPS='$^' $@.target

test_UpdateGains: test_UpdateGains.o ${GRAPH:%=%.o} ${MATCH:%=%.o} ${HKLFM:%=%.o}
	make -r OBJDEPS='$^' $@.target

test_MoveVertex: test_MoveVertex.o ${GRAPH:%=%.o} ${MATCH:%=%.o} ${HKLFM:%=%.o}
	make -r OBJDEPS='$^' $@.target

test_ClearMoveLog: test_ClearMoveLog.o ${GRAPH:%=%.o} ${MATCH:%=%.o} ${HKLFM:%=%.o}
	make -r OBJDEPS='$^' $@.target

test_QSort: test_QSort.o ${SORT:%=%.o}
	make -r OBJDEPS='$^' $@.target

test_quicksort: test_quicksort.o ${SORT:%=%.o}
	make -r OBJDEPS='$^' $@.target

test_CSort: test_CSort.o ${SORT:%=%.o}
	make -r OBJDEPS='$^' $@.target

test_RandomPermute: test_RandomPermute.o ${SORT:%=%.o}
	make -r OBJDEPS='$^' $@.target

test_Random1: test_Random1.o ${SORT:%=%.o}
	make -r OBJDEPS='$^' $@.target

test_SetRandomSeed: test_SetRandomSeed.o ${SORT:%=%.o}
	make -r OBJDEPS='$^' $@.target

test_SwapLong: test_SwapLong.o ${SORT:%=%.o}
	make -r OBJDEPS='$^' $@.target

test_SwapDouble: test_SwapDouble.o ${SORT:%=%.o}
	make -r OBJDEPS='$^' $@.target

test_BucketInsert: test_BucketInsert.o ${GAINBUCKET:%=%.o} ${IO:%=%.o}
	make -r OBJDEPS='$^' $@.target

test_BucketMove: test_BucketMove.o ${GAINBUCKET:%=%.o} ${IO:%=%.o}
	make -r OBJDEPS='$^' $@.target

test_BucketDeleteMax: test_BucketDeleteMax.o ${GAINBUCKET:%=%.o} ${IO:%=%.o}
	make -r OBJDEPS='$^' $@.target

test_GainBucketGetMaxVal: test_GainBucketGetMaxVal.o ${GAINBUCKET:%=%.o} ${IO:%=%.o}
	make -r OBJDEPS='$^' $@.target

test_ClearGainBucket: test_ClearGainBucket.o ${GAINBUCKET:%=%.o} ${IO:%=%.o}
	make -r OBJDEPS='$^' $@.target

test_InitNprocs: test_InitNprocs.o ${DISTRIBUTEVECLIB:%=%.o}
	make -r OBJDEPS='$^' $@.target

test_InitProcindex: test_InitProcindex.o ${DISTRIBUTEVECLIB:%=%.o}
	make -r OBJDEPS='$^' $@.target

test_InitIntersection: test_InitIntersection.o ${DISTRIBUTEVECORIGEQ:%=%.o}
	make -r OBJDEPS='$^' $@.target

test_FindProcLowestSum: test_FindProcLowestSum.o ${DISTRIBUTEVECORIGEQ:%=%.o}
	make -r OBJDEPS='$^' $@.target

test_DistributeVecOrigEq: test_DistributeVecOrigEq.o ${DISTRIBUTEVECORIGEQ:%=%.o}
	make -r OBJDEPS='$^' $@.target

test_DistributeVec: test_DistributeVec.o ${DISTRIBUTEVEC:%=%.o}
	make -r OBJDEPS='$^' $@.target

test_AssignColumnToProc: test_AssignColumnToProc.o ${DISTRIBUTEVECLIB:%=%.o}
	make -r OBJDEPS='$^' $@.target

test_RemoveColumnFromProc: test_RemoveColumnFromProc.o ${DISTRIBUTEVECLIB:%=%.o}
	make -r OBJDEPS='$^' $@.target

test_AssignRemainingColumns: test_AssignRemainingColumns.o ${DISTRIBUTEVECLIB:%=%.o}
	make -r OBJDEPS='$^' $@.target

test_AssignRemainingNonemptyColumns: test_AssignRemainingNonemptyColumns.o ${DISTRIBUTEVECLIB:%=%.o}
	make -r OBJDEPS='$^' $@.target

test_CalcCom: test_CalcCom.o ${DISTRIBUTEVECLIB:%=%.o}
	make -r OBJDEPS='$^' $@.target

test_CalcLocalLowerBound: test_CalcLocalLowerBound.o ${DISTRIBUTEVECLIB:%=%.o}
	make -r OBJDEPS='$^' $@.target

test_matallocl: test_matallocl.o ${MATALLOC:%=%.o}
	make -r OBJDEPS='$^' $@.target

test_matfreel: test_matfreel.o ${MATALLOC:%=%.o}
	make -r OBJDEPS='$^' $@.target

test_FirstNzInRow: test_FirstNzInRow.o ${MATALLOC:%=%.o} ${DISTRIBUTEVECOPT2:%=%.o}
	make -r OBJDEPS='$^' $@.target

test_DistributeVecOpt2: test_DistributeVecOpt2.o ${MATALLOC:%=%.o} ${DISTRIBUTEVECOPT2:%=%.o}
	make -r OBJDEPS='$^' $@.target

test_DistributeVecLocal: test_DistributeVecLocal.o ${DISTRIBUTEVECLOCAL:%=%.o}
	make -r OBJDEPS='$^' $@.target

test_PrintVecLocalStatistics: test_PrintVecLocalStatistics.o ${DISTRIBUTEVECLOCAL:%=%.o}
	make -r OBJDEPS='$^' $@.target

test_PrintVecStatistics: test_PrintVecStatistics.o ${DISTRIBUTEVECLIB:%=%.o}
	make -r OBJDEPS='$^' $@.target

test_PrintHistogram: test_PrintHistogram.o ${DISTRIBUTEVECLIB:%=%.o}
	make -r OBJDEPS='$^' $@.target

test_PrintCom: test_PrintCom.o ${DISTRIBUTEVECLIB:%=%.o}
	make -r OBJDEPS='$^' $@.target

test_WriteVectorDistribution: test_WriteVectorDistribution.o ${DISTRIBUTEVECLIB:%=%.o}
	make -r OBJDEPS='$^' $@.target

test_DistributeVecGreedyImprove: test_DistributeVecGreedyImprove.o ${DISTRIBUTEVECGREEDY:%=%.o}
	make -r OBJDEPS='$^' $@.target

test_InitSums: test_InitSums.o ${DISTRIBUTEVECLIB:%=%.o} ${DISTRIBUTEVECORIG:%=%.o}
	make -r OBJDEPS='$^' $@.target

test_GenerateHistogram: test_GenerateHistogram.o ${DISTRIBUTEVECLIB:%=%.o}
	make -r OBJDEPS='$^' $@.target

test_DistributeVecOrig: test_DistributeVecOrig.o ${DISTRIBUTEVECORIG:%=%.o}
	make -r OBJDEPS='$^' $@.target

test_BalanceParts: test_BalanceParts.o ${DISTRIBUTEMAT:%=%.o}
	make -r OBJDEPS='$^' $@.target

test_DetermineSplit: test_DetermineSplit.o ${DISTRIBUTEMAT:%=%.o}
	make -r OBJDEPS='$^' $@.target

test_logb2: test_logb2.o ${DISTRIBUTEMAT:%=%.o}
	make -r OBJDEPS='$^' $@.target

test_ComputeWeight: test_ComputeWeight.o ${DISTRIBUTEMAT:%=%.o}
	make -r OBJDEPS='$^' $@.target

test_SplitMatrixSimple: test_SplitMatrixSimple.o ${DISTRIBUTEMAT:%=%.o}
	make -r OBJDEPS='$^' $@.target

test_SplitMatrixKLFM: test_SplitMatrixKLFM.o ${DISTRIBUTEMAT:%=%.o} ${DISTRIBUTEVECLIB:%=%.o}
	make -r OBJDEPS='$^' $@.target

test_DistributeMatrixMondriaan: test_DistributeMatrixMondriaan.o ${DISTRIBUTEMAT:%=%.o}
	make -r OBJDEPS='$^' $@.target

test_LambdaLambdaMinusOneMetric: test_LambdaLambdaMinusOneMetric.o ${DISTRIBUTEMAT:%=%.o}
	make -r OBJDEPS='$^' $@.target

test_FindOptimalPathMatching: test_FindOptimalPathMatching.o ${DISTRIBUTEMAT:%=%.o}
	make -r OBJDEPS='$^' $@.target

test_MatchATA: test_MatchATA.o ${SPARSEMATRIX:%=%.o} ${GRAPH:%=%.0} ${MATCH:%=%.o}
	make -r OBJDEPS='$^' $@.target

test_CreateNewBiPartHyperGraph: test_CreateNewBiPartHyperGraph.o ${GRAPH:%=%.o}
	make -r OBJDEPS='$^' $@.target

test_DeleteBiPartHyperGraph: test_DeleteBiPartHyperGraph.o ${GRAPH:%=%.o}
	make -r OBJDEPS='$^' $@.target

test_SparseMatrix2BiPartHyperGraph: test_SparseMatrix2BiPartHyperGraph.o ${GRAPH:%=%.o}
	make -r OBJDEPS='$^' $@.target

test_BiPartHyperGraph2SparseMatrix: test_BiPartHyperGraph2SparseMatrix.o ${GRAPH:%=%.o} ${MATALLOC:%=%.o}
	make -r OBJDEPS='$^' $@.target

test_SparseMatrixSymmetric2Full: test_SparseMatrixSymmetric2Full.o ${SPARSEMATRIX:%=%.o} ${MATALLOC:%=%.o}
	make -r OBJDEPS='$^' $@.target

test_SparseMatrixSymmetricLower2Random: test_SparseMatrixSymmetricLower2Random.o ${SPARSEMATRIX:%=%.o} ${MATALLOC:%=%.o}
	make -r OBJDEPS='$^' $@.target

test_SparseMatrixSymmetricRandom2Lower: test_SparseMatrixSymmetricRandom2Lower.o ${SPARSEMATRIX:%=%.o} ${MATALLOC:%=%.o}
	make -r OBJDEPS='$^' $@.target

test_AddDummiesToSparseMatrix: test_AddDummiesToSparseMatrix.o ${SPARSEMATRIX:%=%.o} ${MATALLOC:%=%.o}
	make -r OBJDEPS='$^' $@.target

test_RemoveDummiesFromSparseMatrix: test_RemoveDummiesFromSparseMatrix.o ${SPARSEMATRIX:%=%.o} ${MATALLOC:%=%.o}
	make -r OBJDEPS='$^' $@.target

test_MMSparseMatrixInit: test_MMSparseMatrixInit.o ${SPARSEMATRIX:%=%.o} ${MATALLOC:%=%.o}
	make -r OBJDEPS='$^' $@.target

test_MMDeleteSparseMatrix: test_MMDeleteSparseMatrix.o ${SPARSEMATRIX:%=%.o} ${MATALLOC:%=%.o}
	make -r OBJDEPS='$^' $@.target

test_MMSparseMatrixAllocateMemory: test_MMSparseMatrixAllocateMemory.o ${SPARSEMATRIX:%=%.o} ${MATALLOC:%=%.o}
	make -r OBJDEPS='$^' $@.target

test_MMSparseMatrixFreeMemory: test_MMSparseMatrixFreeMemory.o ${SPARSEMATRIX:%=%.o} ${MATALLOC:%=%.o}
	make -r OBJDEPS='$^' $@.target

test_MMSparseMatrixReadEntries: test_MMSparseMatrixReadEntries.o ${SPARSEMATRIX:%=%.o} ${MATALLOC:%=%.o}
	make -r OBJDEPS='$^' $@.target

test_MMSparseMatrixReadPstart: test_MMSparseMatrixReadPstart.o ${SPARSEMATRIX:%=%.o} ${MATALLOC:%=%.o}
	make -r OBJDEPS='$^' $@.target

test_MMSparseMatrixReadWeights: test_MMSparseMatrixReadWeights.o ${SPARSEMATRIX:%=%.o} ${MATALLOC:%=%.o}
	make -r OBJDEPS='$^' $@.target

test_MMSparseMatrixReadTail: test_MMSparseMatrixReadTail.o ${SPARSEMATRIX:%=%.o} ${MATALLOC:%=%.o}
	make -r OBJDEPS='$^' $@.target

test_MMSparseMatrixReadHeader: test_MMSparseMatrixReadHeader.o ${SPARSEMATRIX:%=%.o} ${MATALLOC:%=%.o}
	make -r OBJDEPS='$^' $@.target

test_MMReadSparseMatrix: test_MMReadSparseMatrix.o ${SPARSEMATRIX:%=%.o} ${MATALLOC:%=%.o}
	make -r OBJDEPS='$^' $@.target

test_MMSparseMatrixPrintHeader: test_MMSparseMatrixPrintHeader.o ${SPARSEMATRIX:%=%.o} ${MATALLOC:%=%.o}
	make -r OBJDEPS='$^' $@.target

test_MMSparseMatrixPrintTail: test_MMSparseMatrixPrintTail.o ${SPARSEMATRIX:%=%.o} ${MATALLOC:%=%.o}
	make -r OBJDEPS='$^' $@.target

test_MMSparseMatrixPrintPstart: test_MMSparseMatrixPrintPstart.o ${SPARSEMATRIX:%=%.o} ${MATALLOC:%=%.o}
	make -r OBJDEPS='$^' $@.target

test_MMSparseMatrixPrintWeights: test_MMSparseMatrixPrintWeights.o ${SPARSEMATRIX:%=%.o} ${MATALLOC:%=%.o}
	make -r OBJDEPS='$^' $@.target

test_MMSparseMatrixPrintEntries: test_MMSparseMatrixPrintEntries.o ${SPARSEMATRIX:%=%.o} ${MATALLOC:%=%.o}
	make -r OBJDEPS='$^' $@.target

test_MMWriteSparseMatrix: test_MMWriteSparseMatrix.o ${SPARSEMATRIX:%=%.o} ${MATALLOC:%=%.o}
	make -r OBJDEPS='$^' $@.target

test_MMSparseMatrixGetTypeCode: test_MMSparseMatrixGetTypeCode.o ${SPARSEMATRIX:%=%.o} ${MATALLOC:%=%.o}
	make -r OBJDEPS='$^' $@.target

test_MMSparseMatrixSetTypeCode: test_MMSparseMatrixSetTypeCode.o ${SPARSEMATRIX:%=%.o} ${MATALLOC:%=%.o}
	make -r OBJDEPS='$^' $@.target

test_SparseMatrixRemoveDuplicates: test_SparseMatrixRemoveDuplicates.o ${SPARSEMATRIX:%=%.o} ${MATALLOC:%=%.o}
	make -r OBJDEPS='$^' $@.target

test_SparseMatrixStructurallySymmetric: test_SparseMatrixStructurallySymmetric.o ${SPARSEMATRIX:%=%.o} ${MATALLOC:%=%.o}
	make -r OBJDEPS='$^' $@.target

test_CRS2CCS: test_CRS2CCS.o ${CARTESIAN:%=%.o}
	make -r OBJDEPS='$^' $@.target

test_MMWriteCartesianSubmatrices: test_MMWriteCartesianSubmatrices.o ${CARTESIAN:%=%.o}
	make -r OBJDEPS='$^' $@.target

test_GetParameters: test_GetParameters.o ${OPTIONS:%=%.o}
	make -r OBJDEPS='$^' $@.target

test_SetDefaultOptions: test_SetDefaultOptions.o ${OPTIONS:%=%.o}
	make -r OBJDEPS='$^' $@.target

test_SetOption: test_SetOption.o ${OPTIONS:%=%.o}
	make -r OBJDEPS='$^' $@.target

test_PrintHelp: test_PrintHelp.o ${OPTIONS:%=%.o}
	make -r OBJDEPS='$^' $@.target

%.cstub:
	rm -f ${@:%.cstub=%}.o

%.vcstub:
	rm -f ${@:%.vcstub=%}

veryclean: clean ${TESTTARGETS:%=%.vcstub} ;

clean: ${TESTTARGETS:%=%.cstub}
	rm -f `find . -name 'core'`

