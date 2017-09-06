#include <Mondriaan.h>

void PrintStatsHelp(int argc, char **argv);

long *ReadVectorFile(char *filename, int dir, struct sparsematrix *pA, int inspecting) {
	FILE *file = fopen(filename, "r");
	if (file) {
		long vecLen;
		int vecP;
		long *vec = ReadVector(0, &vecLen, &vecP, file);
		if(vecLen != ((dir==COL)?pA->m:pA->n) || vecP != pA->NrProcs) {
			fprintf(stderr, "main(): Invalid vector size!\n");
			exit(-1);
		}
		fclose(file);
		return vec;
	}
	else if(!inspecting) {
		fprintf(stderr, "main(): Unable to open '%s' for reading!\n", filename);
		exit(-1);
	}
	return NULL;
}


/* This function calculates the total communication volume for
 * a given distributed sparse matrix A, as produced by MondriaanOpt.
 * Typically these files end with -I2f, and contain a distribution
 * over two processors with a third virtual processor containing the
 * free nonzeros.
 * Note that these numbers are for only one direction (vector).

 * Input:
 *    A distributed matrix,
 *    dir = ROW (for distribution of v) or dir = COL (for u).
 * Output:
 *    ComVol is total communication volume for direction dir,
 */
int CalcComI2f(const struct sparsematrix *pM, int dir, long *ComVol) {
	
	int q, *procs;
	long l=0, j, t, *index;
	
	if (!pM || !ComVol) {
		fprintf(stderr, "CalcComI2f(): Null arguments!\n");
		return FALSE;
	}
	
	if (dir == COL) {
		l = pM->m;
		index = pM->i;
	} else if (dir == ROW) {
		l = pM->n;
		index = pM->j;
	} else {
		fprintf(stderr, "CalcComI2f(): Unknown direction!\n");
		return FALSE;
	}
	
	procs = (int *)calloc(l, sizeof(int));
	if (procs == NULL) {
		fprintf(stderr, "CalcComI2f(): Not enough memory!\n");
		return FALSE;
	}

	q=0;
	for (t=0; t<pM->NrNzElts; t++) {
		while(t == pM->Pstart[q+1])
			q++;  /* this terminates because t != pM->Pstart[P] = pM->NrNzElts */

		j = index[t];
		procs[j] |= (q+1);
	}
	
	*ComVol = 0;
	for (j=0; j<l; j++) {
		if(procs[j] == 3) {
			++(*ComVol);
		}
	}
	
	free(procs);
	
	return TRUE;
} /* end CalcComI2f */

/* This function calculates the imbalance for
 * a given distributed sparse matrix A, as produced by MondriaanOpt.
 * Typically these files end with -I2f, and contain a distribution
 * over two processors with a third virtual processor containing the
 * free nonzeros.

 * Input:
 *    A distributed matrix,
 * Output:
 *    Imbalance is the best achievable imbalance
 *    ImbalanceOpt is the best achievable imbalance as defined by MondriaanOpt
 */
int CalcImbalanceI2f(const struct sparsematrix *pM, double *Imbalance, double *ImbalanceOpt) {
	
	int q;
	
	if (!pM || !Imbalance || !ImbalanceOpt) {
		fprintf(stderr, "CalcComI2f(): Null arguments!\n");
		return FALSE;
	}
	
	long weights[pM->NrProcs];
	for(q=0; q<pM->NrProcs; ++q) {
		weights[q] = pM->Pstart[q+1] - pM->Pstart[q];
	}
	
	long smallWeight = (weights[0]<weights[1])?weights[0]:weights[1];
	long largeWeight = (weights[0]>weights[1])?weights[0]:weights[1];
	long freeWeight = (pM->NrProcs==3)?weights[2]:0;
	
	long diffWeight = largeWeight - smallWeight;
	long addWeight = (diffWeight<freeWeight)?diffWeight:freeWeight;
	smallWeight += addWeight;
	freeWeight -= addWeight;
	
	smallWeight += freeWeight/2;
	largeWeight += (freeWeight+1)/2;
	
	*Imbalance = (2*largeWeight)/(double)pM->NrNzElts - 1;
	*ImbalanceOpt = (2*largeWeight)/(double)(pM->NrNzElts + (pM->NrNzElts%2)) - 1;
	
	return TRUE;
} /* end CalcImbalanceI2f */


int main(int argc, char **argv) {
	
	if(argc < 2 || !strcmp(argv[1], "-h")) {
		PrintStatsHelp(argc, argv);
		exit((argc < 2)? 1 : 0);
	}
	
	/* Read options */
	int valuesAsProcs = FALSE, I2f = FALSE;
	int auto_uv = TRUE, auto_I2f = TRUE;
	int i;
	char *filename_u = NULL, *filename_v = NULL;
	for(i=2; i<argc; ++i) {
		if(!strcmp(argv[i], "-I")) {
			valuesAsProcs = TRUE;
		}
		if(!strcmp(argv[i], "-h")) {
			PrintStatsHelp(argc, argv);
			exit(0);
		}
		if(!strcmp(argv[i], "-I2f")) {
			valuesAsProcs = TRUE;
			I2f = TRUE;
		}
		if(!strcmp(argv[i], "-no-auto-uv")) {
			auto_uv = FALSE;
		}
		if(!strcmp(argv[i], "-no-auto-I2f")) {
			auto_I2f = FALSE;
		}
		if(!strcmp(argv[i], "-u") && i+1<argc) {
			filename_u = argv[i+1];
			++i;
		}
		if(!strcmp(argv[i], "-v") && i+1<argc) {
			filename_v = argv[i+1];
			++i;
		}
	}
	
	
	char *filename = argv[1];
	FILE *File;
	struct sparsematrix A;
	
	/* Read matrix file from disk or standard input. */
	if (!strcmp(filename, "-") || !strcmp(filename, "stdin")) {
		if (!MMReadSparseMatrix(stdin, &A)) {
			fprintf(stderr, "main(): Could not read matrix from standard input!\n");
			exit(-1);
		}
		
		filename = "stdin";
	}
	else {
		
		File = fopen(filename, "r");
		
		if (!File) {
			fprintf(stderr, "main(): Unable to open '%s' for reading!\n", filename);
			exit(-1);
		}
		
		if (!MMReadSparseMatrix(File, &A)) {
			fprintf(stderr, "main(): Could not read matrix!\n");
			exit(-1);
		}
		
		fclose(File);
	}
	
	/* If required, interpret matrix values as processor indices */
	if(valuesAsProcs) {
		if(!SpMatValuesToProcessorIndices(&A)) {
			fprintf(stderr, "main(): Could not interpret values as processor indices!\n");
			exit(-1);
		}
		
		if(!I2f && auto_I2f && !strcmp(strrchr(filename, '-'), "-I2f")) {
			fprintf(stderr, " -- Detected I2f format.\n");
			I2f = TRUE;
		}
	}
	
	if(A.MMTypeCode[0]!='D') {
		fprintf(stderr, "main(): Input matrix must be distributed!\n");
		exit(-1);
	}
	
	/* If matrix is symmetric, turn it to full form */
	if(A.m == A.n && (A.MMTypeCode[3] == 'S' || A.MMTypeCode[3] == 'K' || A.MMTypeCode[3] == 'H')) {
		if(!SparseMatrixSymmetric2Full(&A)) {
			fprintf(stderr, "main(): Could not transform symmetric matrix to full format!\n");
			exit(-1);
		}
	}
	
	/* If needed, try to automatically determine vector distribution files */
	long *u = NULL, *v = NULL;
	if(auto_uv && !I2f) {
		char basename[MAX_WORD_LENGTH], filename_vec[MAX_WORD_LENGTH];
		strcpy(basename, filename);
		char *suffix = strrchr(basename, '-');
		if(suffix != NULL) {
			*suffix = '\0';
		}
		
		if(filename_u == NULL) {
			sprintf(filename_vec, "%s-u%d", basename, A.NrProcs);
			u = ReadVectorFile(filename_vec, COL, &A, TRUE);
			if(u != NULL) {
				fprintf(stderr, " -- Found %s\n", filename_vec);
			}
		}
		
		if(filename_v == NULL) {
			sprintf(filename_vec, "%s-v%d", basename, A.NrProcs);
			v = ReadVectorFile(filename_vec, ROW, &A, TRUE);
			if(v != NULL) {
				fprintf(stderr, " -- Found %s\n", filename_vec);
			}
		}
	}
	
	/* Read manually passed u/v files */
	if(filename_u != NULL) {
		u = ReadVectorFile(filename_u, COL, &A, FALSE);
	}
	if(filename_v != NULL) {
		v = ReadVectorFile(filename_v, ROW, &A, FALSE);
	}
	
	/* Print general statistics */
	printf("Matrix:\n");
	printf("  m             : %ld\n", A.m);
	printf("  n             : %ld\n", A.n);
	printf("  nnz           : %ld\n", A.NrNzElts);
	printf("  P             : %d\n", I2f?2:A.NrProcs);
	printf("\n");
	
	if(I2f) {
		/* Print weights */
		printf("Weights:\n");
		printf("  proc 0        : %ld\n", A.Pstart[1]-A.Pstart[0]);
		printf("  proc 1        : %ld\n", A.Pstart[2]-A.Pstart[1]);
		printf("  free          : %ld\n", (A.NrProcs==3)?A.Pstart[3]-A.Pstart[2]:0);
		printf("\n");
		
		/* Print superstep v */
		long ComVolRow;
		CalcComI2f(&A, ROW, &ComVolRow);
		printf("Communication for vector v:\n");
		printf("  vol           : %ld\n", ComVolRow);
		printf("  avg = vol/P   : %g\n", ComVolRow/2.0);
		printf("\n");
		
		/* Print superstep u */
		long ComVolCol;
		CalcComI2f(&A, COL, &ComVolCol);
		printf("Communication for vector u:\n");
		printf("  vol           : %ld\n", ComVolCol);
		printf("  avg = vol/P   : %g\n", ComVolCol/2.0);
		printf("\n");
		
		/* Print total */
		double Imbalance, ImbalanceOpt;
		CalcImbalanceI2f(&A, &Imbalance, &ImbalanceOpt);
		printf("Total:\n");
		printf("  vol           : %ld\n", ComVolRow+ComVolCol);
		printf("  avg = vol/P   : %g\n", (ComVolRow+ComVolCol)/2.0);
		printf("  imbalance     : %g\n", Imbalance);
		printf("  imbalance(opt): %g\n", ImbalanceOpt);
		printf("\n");
		
	}
	else {
		/* Print weights */
		printf("Weights & imbalances:\n");
		long maxWeight = 0, weight;
		for(int q=0; q<A.NrProcs; ++q) {
			weight = A.Pstart[q+1] - A.Pstart[q];
			if(weight > maxWeight) {
				maxWeight = weight;
			}
			printf("  proc %-9d: %ld (%g)\n", q, weight, weight/(A.NrNzElts/(double)A.NrProcs) -1);
		}
		printf("\n");
		
		/* Print superstep v */
		long ComVolRow, MaxOutRow, MaxInRow, MaxCompntsRow, TotCompntsRow;
		CalcCom(&A, v, ROW, &ComVolRow, &MaxOutRow, &MaxInRow, &MaxCompntsRow, &TotCompntsRow);
		PrintCom(A.NrProcs, A.n, ROW, ComVolRow, MaxOutRow, MaxInRow, MaxCompntsRow, TotCompntsRow);
		printf("\n");
		
		/* Print superstep u */
		long ComVolCol, MaxOutCol, MaxInCol, MaxCompntsCol, TotCompntsCol;
		CalcCom(&A, u, COL, &ComVolCol, &MaxOutCol, &MaxInCol, &MaxCompntsCol, &TotCompntsCol);
		PrintCom(A.NrProcs, A.m, COL, ComVolCol, MaxOutCol, MaxInCol, MaxCompntsCol, TotCompntsCol);
		printf("\n");
		
		/* Print total */
		long MaxRow = (MaxInRow>MaxOutRow)?MaxInRow:MaxOutRow;
		long MaxCol = (MaxInCol>MaxOutCol)?MaxInCol:MaxOutCol;
		double Imbalance = maxWeight/(A.NrNzElts/(double)A.NrProcs) - 1;
		double ImbalanceOpt = maxWeight/ceil(A.NrNzElts/(double)A.NrProcs) - 1;
		printf("Total:\n");
		printf("  vol           : %ld\n", ComVolRow+ComVolCol);
		printf("  avg = vol/P   : %g\n", (ComVolRow+ComVolCol)/(double)A.NrProcs);
		printf("  max           : %ld\n", MaxRow+MaxCol);
		printf("  imbalance     : %g\n", Imbalance);
		printf("  imbalance(opt): %g\n", ImbalanceOpt);
		printf("\n");
	}
	
	/* Finish */
	if(u != NULL) {
		free(u);
	}
	if(v != NULL) {
		free(v);
	}

	MMDeleteSparseMatrix(&A);
	
	exit(0);
	
}

void PrintStatsHelp(int argc, char **argv) {

	printf("\nMondriaan version %s.\n", MONDRIAANVERSION);
	printf("Usage: %s [matrix] <options>\n\n", (argv && argv[0])?argv[0]:"./MondriaanStats");
	printf(" [matrix] - the partitioned matrix\n");
	printf(" <options>:\n");
	printf(" -I            - Interpret the (real) values of the nonzeros as processor indices.\n");
	printf("                 Typically used with matrix.mtx-Ix files.\n");
	printf(" -I2f          - Input file contains a bipartitioning, possibly with free nonzeros.\n");
	printf("                 Typically used with matrix.mtx-I2f files. -I2f implies -I\n");
	printf(" -no-auto-uv   - Do not determine the names of the u and v vector files automatically.\n");
	printf(" -no-auto-I2f  - Do not automatically determine -I2f suffix if -I is passed.\n");
	printf(" -u <file>     - Use <file> as distribution for u.\n");
	printf(" -v <file>     - Use <file> as distribution for v.\n");
	printf(" -h            - Show this help.\n");
	printf("\n");
	printf("When no options are given, the input matrix is assumed to be a distributed-matrix.\n");
	printf("\n");
	fflush(stdout);

} /* end PrintHelp */


