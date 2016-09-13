/*
Profile.c

Created by Bas Fagginger Auer.

This program profiles the performance of the Mondriaan library by measuring communication volume, calculation time and imbalance.

Usage:

./Profile 4 0.1 16 -SplitStrategy=hybrid a.mtx b.mtx c.mtx ...

  -- profiles a.mtx, b.mtx, c.mtx, ... by distributing them over 4 processors with imbalance 0.1 where we average over 16 attempts.
  To profile for multiple processor counts, list them as 2,4,6 (without spaces).
*/
#include <stdlib.h>
#include <stdio.h>
#include <memory.h>
#include <math.h>
#include <time.h>
#include <Mondriaan.h>

#define NUM_ATTRIBUTES 4

struct sMatrixAttribute {
	double Average, Variance, *Data;
};

struct sMatrixData {
	char File[256];
	
	long NumNz, Rows, Cols, P;
	struct sMatrixAttribute Attributes[NUM_ATTRIBUTES];
};

int NumProcs[256], NumNumProcs = 0;
struct sMatrixData *Matrices = NULL;
int NumMatrices = 0, NumSamples = 1;

int SetupAttributes(struct sMatrixData *Matrix) {
	/* Sets up the matrix attributes and allocates memory. */
	int i;
	
	if (Matrix == NULL || NumSamples <= 0) return FALSE;
	
	Matrix->NumNz = Matrix->Rows = Matrix->Cols = 0;
	
	for (i = 0; i < NUM_ATTRIBUTES; ++i) {
		Matrix->Attributes[i].Average = 0.0;
		Matrix->Attributes[i].Variance = 0.0;
		
		if ((Matrix->Attributes[i].Data = (double *)malloc(NumSamples*sizeof(double))) == NULL) return FALSE;
	}
	
	return TRUE;
}

int AverageAndFreeAttributes(struct sMatrixData *Matrix) {
	/* Calculates averages and variances of measured attributes and frees their data. */
	int i, j;
	
	if (Matrix == NULL) return FALSE;
	
	for (i = 0; i < NUM_ATTRIBUTES; ++i) {
		double Avg = 0.0, Var = 0.0;
		
		if (Matrix->Attributes[i].Data == NULL) return FALSE;
		
		for (j = 0; j < NumSamples; ++j) Avg += Matrix->Attributes[i].Data[j];
		
		Avg /= (double)NumSamples;
		
		for (j = 0; j < NumSamples; ++j) Var += (Matrix->Attributes[i].Data[j] - Avg)*(Matrix->Attributes[i].Data[j] - Avg);
		
		Var /= (double)NumSamples;
		
		free(Matrix->Attributes[i].Data);
		
		Matrix->Attributes[i].Average = Avg;
		Matrix->Attributes[i].Variance = Var;
		Matrix->Attributes[i].Data = NULL;
	}
	
	return TRUE;
}

int ParseCommandLineOptions(struct opts *pOptions, const int argc, char **argv) {
	int i, j;
	
	/* Check whether or not the given options are empty. */
	if (argc <= 0 || !argv || !pOptions) {
		fprintf(stderr, "ParseCommandLineOptions(): null arguments!\n");
		return FALSE;
	}
   
	/* Print Help if no option is given or a help option */
	if (argc == 1 || !strcmp(argv[1], "-h") || !strcmp(argv[1], "-help") || ! strcmp(argv[1], "--help")) {
		fprintf(stderr, "Usage: %s <nr procs.> <imbalance> <nr samples> a.mtx b.mtx c.mtx ...\n\n", argv[0]);
		return FALSE;
	}
  	
	/* Check the number of arguments */
	if (argc < 4) {
		fprintf(stderr, "ParseCommandLineOptions(): invalid number of arguments!\n");
		return FALSE;
	}
  
	/* Read the required arguments */
	pOptions->matrix = 0;
	
	if (TRUE) {
		char *P = argv[1];
		
		NumNumProcs = 0;
		
		while (strchr(P, ',') && NumNumProcs < 255) {
			if (sscanf(P, "%d,", &NumProcs[NumNumProcs++]) != 1) break;
			P = strchr(P, ',') + 1;
		}
		
		NumProcs[NumNumProcs++] = atoi(P);
		
		for (i = 0; i < NumNumProcs; ++i) {
			if (NumProcs[i] <= 0) fprintf(stderr, "ParseCommandLineOPtions(): invalid number of processors!\n");
		}
	}
	
	if (atof(argv[2]) < 0.0) fprintf(stderr, "ParseCommandLineOptions(): imbalance out of range!\n");
	else pOptions->eps = atof(argv[2]);
	
	if (atoi(argv[3]) <= 0) fprintf(stderr, "ParseCommandLineOptions(): invalid number of samples!\n");
	else NumSamples = atoi(argv[3]);
	
	/* Read the optional command line arguments, e.g.
		   -SplitStrategy=onedimrow
	   to overrule the corresponding defaults */
	i = 3;
	
	while (++i < argc) {
		if (argv[i][0] != '-') {
			NumMatrices++;
		}
		else {
			const char *Option = argv[i] + 1, *Value = 0;

			/* Determine the form of the option: -option=value or -option value. */
			if (strchr(Option, '=')) {
				Value = strchr(Option, '=') + 1;
				*(strchr(Option, '=')) = 0;
			}
			else if (++i < argc) {
				Value = argv[i];
			}
			else {
				fprintf(stderr, "ParseCommandLineOptions(): warning, missing value for '%s'!\n", Option);
			}

			/* Set the read option. */
			SetOption(pOptions, Option, Value);
		}
	}
	
	Matrices = (struct sMatrixData *)malloc(NumMatrices*NumNumProcs*sizeof(struct sMatrixData));
	
	memset(Matrices, 0, NumMatrices*NumNumProcs*sizeof(struct sMatrixData));
	
	if (Matrices == NULL) {
		fprintf(stderr, "ParseCommandLineOptions(): Not enough memory!\n");
		return FALSE;
	}
	
	NumMatrices = 0;
	i = 3;
	
	while (++i < argc) {
		if (argv[i][0] != '-') {
			for (j = 0; j < NumNumProcs; ++j) {
				struct sMatrixData *Data = &Matrices[NumMatrices*NumNumProcs + j];
				
				Data->P = NumProcs[j];
				strncpy(Data->File, argv[i], 256);
				Data->File[255] = '\0';
			}
			++NumMatrices;
		}
		else if (!strchr(argv[i], '=')) {
			++i;
		}
	}
	
	return TRUE;
}

void DoubleToLaTeX(char *Dest, const double Val, int NumDigits) {
	/* Small function that turns 1.423e-04 into 1.423 \cdot 10^{-4}. */
	char Tmp1[64], Tmp2[1024], *dp = Dest;
	int i;
	
	sprintf(Tmp1, "%%.%dle", NumDigits);
	sprintf(Tmp2, Tmp1, Val);
	
	for (i = 0; Tmp2[i] != '\0' && Tmp2[i] != 'e'; ++i) *dp++ = Tmp2[i];
	++i;
	*dp++ = ' ';
	*dp++ = '\\';
	*dp++ = 'c';
	*dp++ = 'd';
	*dp++ = 'o';
	*dp++ = 't';
	*dp++ = ' ';
	*dp++ = '1';
	*dp++ = '0';
	*dp++ = '^';
	*dp++ = '{';
	*dp = 'z';
	for ( ; Tmp2[i] != '\0'; ++i) if (Tmp2[i] != '0' && Tmp2[i] != '+') *dp++ = Tmp2[i];
	if (*dp == 'z') *dp++ = '0';
	*dp++ = '}';
	*dp++ = '\0';
}

int main(int argc, char **argv) {
	struct opts Options;
	int i;
	time_t CurrentTime;
	clock_t TotalClock = clock();
	
	/* Apply Mondriaan options. */
	SetDefaultOptions(&Options);
	
	if (!SetOptionsFromFile(&Options, "Mondriaan.defaults")) {
		fprintf(stderr, "main(): warning, cannot set options from 'Mondriaan.defaults', using default options!\n");
	}
	
	if (!ParseCommandLineOptions(&Options, argc, argv)) {
		fprintf(stderr, "main(): invalid command line parameters!\n");
		exit(EXIT_FAILURE);
	}
	
	if (!ApplyOptions(&Options)) {
		fprintf(stderr, "main(): could not apply given options!\n");
		exit(EXIT_FAILURE);
	}
	
	if (NumMatrices <= 0 || Matrices == NULL) {
		fprintf(stderr, "main(): Invalid number of supplied matrices or samples!\n");
		exit(EXIT_FAILURE);
	}
	
	/* Start profiling ... */
	fprintf(stderr, "Profiling Mondriaan for %d matrices, %d samples, %s processors, and %f imbalance.\n", NumMatrices, NumSamples, argv[1], Options.eps);
	
	for (i = 0; i < NumMatrices*NumNumProcs; ++i) {
		int j;
		
		fprintf(stderr, "[% 4d/%d] (% 4ld) %s ", i + 1, NumMatrices*NumNumProcs, Matrices[i].P, Matrices[i].File);
		fflush(stderr);
		
		if (!SetupAttributes(&Matrices[i])) {
			fprintf(stderr, "main(): Cannot setup attributes!\n");
			exit(EXIT_FAILURE);
		}
		
		Options.P = Matrices[i].P;
		
		/* Take the requested number of samples. */
		for (j = 0; j < NumSamples; ++j) {
			struct sparsematrix A;
			long int *UAssign, *VAssign, Symmetric;
			FILE *File;
			long l;
			int k;
			clock_t Clock;
			
			double Duration;
			long MaxNz, MinNz;
			long MaxComU, MaxComV, ComVolU, ComVolV;
			
			fprintf(stderr, ".");
			fflush(stderr);
			
			/* Read matrix from disk. */
			File = fopen(Matrices[i].File, "r");
			
			if (!File) {
				fprintf(stderr, "main(): Could not open '%s' for reading!\n", Matrices[i].File);
				exit(EXIT_FAILURE);
			}
				
			if (!MMReadSparseMatrix(File, &A)) {
				fprintf(stderr, "main(): Could not read matrix!\n");
				exit(EXIT_FAILURE);
			}
			
			fclose(File);
			
			/* Remove double zeroes. */
			if (!SparseMatrixRemoveDuplicates(&A)) exit(EXIT_FAILURE);
			
			/* Check symmetry. */
			if (A.m == A.n && (A.MMTypeCode[3] == 'S' || A.MMTypeCode[3] == 'K' || A.MMTypeCode[3] == 'H')) Symmetric = TRUE;
			else Symmetric = FALSE;
			
			if (Symmetric)
			{
				if (Options.SymmetricMatrix_UseSingleEntry == SingleEntNo) SparseMatrixSymmetric2Full(&A);
				else if (Options.SymmetricMatrix_SingleEntryType == ETypeRandom) SparseMatrixSymmetricLower2Random(&A);
			}
			
			/* Add dummies if requested. */
			if (A.m == A.n && Options.SquareMatrix_DistributeVectorsEqual == EqVecYes && Options.SquareMatrix_DistributeVectorsEqual_AddDummies == DumYes) AddDummiesToSparseMatrix(&A);
			
			/* Initialise processor array. */
			A.NrProcs = Options.P;
			A.Pstart = (long *)malloc((A.NrProcs + 1)*sizeof(long));
			
			if (A.Pstart == NULL) {
				fprintf(stderr, "main(): Cannot allocate processor array!\n");
				exit(EXIT_FAILURE);
			}
			
			A.Pstart[0] = 0;
			
			for (k = 1; k <= A.NrProcs; ++k) {
				A.Pstart[k] = A.NrNzElts;
			}
		
			/* Distribute the processors among the matrix entries. */
			SetRandomSeed(Options.Seed = 137*j + 12345);
			
			/* ==== Start Mondriaan */
			Clock = clock();
			
			if (!DistributeMatrixMondriaan(&A, Options.P, Options.eps, &Options, NULL)) {
				fprintf(stderr, "main(): Unable to distribute matrix!\n");
				exit(EXIT_FAILURE);
			}
			
			/* Remove dummies. */
			if (A.m == A.n && Options.SquareMatrix_DistributeVectorsEqual == EqVecYes && Options.SquareMatrix_DistributeVectorsEqual_AddDummies == DumYes) RemoveDummiesFromSparseMatrix(&A);
			
			/* Convert randomly represented matrix to lower triangular form. */
			if (Symmetric && Options.SymmetricMatrix_UseSingleEntry == SingleEntYes && Options.SymmetricMatrix_SingleEntryType == ETypeRandom) SparseMatrixSymmetricRandom2Lower(&A);
			
			/* Distribute vectors. */
			UAssign = (long int *)malloc(A.m*sizeof(long int));
			VAssign = (long int *)malloc(A.n*sizeof(long int));
			
			if (UAssign == NULL || VAssign == NULL) {
				fprintf(stderr, "main(): Cannot allocate vertex assign arrays!\n");
				exit(EXIT_FAILURE);
			}
			
			/* Convert symmetrically partitioned matrix to full form. */
			if (Symmetric && Options.SymmetricMatrix_UseSingleEntry == SingleEntYes) SparseMatrixSymmetric2Full(&A);
			
			if (A.m == A.n && Options.SquareMatrix_DistributeVectorsEqual == EqVecYes) {
				if (Symmetric && Options.SymmetricMatrix_UseSingleEntry == SingleEntYes) {
					MaxComV = DistributeVec(&A, VAssign, ROW, &Options);
					
					if (MaxComV < 0) {
						fprintf(stderr, "main(): Unable to distribute vector!\n");
						exit(EXIT_FAILURE);
					}
					
					for (k = 0; k < A.m; k++) {
						UAssign[k] = VAssign[k];
					}
						
					MaxComU = MaxComV;
				}
				else {
					MaxComU = DistributeVecOrigEq(&A, UAssign, VAssign, &Options);
					
					if (MaxComU < 0) {
						fprintf(stderr, "main(): Unable to distribute vector!\n");
						exit(EXIT_FAILURE);
					}
					
					MaxComV = 0;
				}
			}
			else {
				MaxComV = DistributeVec(&A, VAssign, ROW, &Options);
				MaxComU = DistributeVec(&A, UAssign, COL, &Options);
				
				if (MaxComV < 0 || MaxComU < 0) {
					fprintf(stderr, "main(): Unable to distribute vector!\n");
					exit(EXIT_FAILURE);
				}
			}
			
			/* ==== Stop Mondriaan */
			/* Calculate duration. */
			Duration = (double)(clock() - Clock)/(double)CLOCKS_PER_SEC;
			
			/* Determine minimum and maximum number of assigned nonzeroes. */
			MaxNz = MinNz = A.Pstart[1] - A.Pstart[0];
			
			for (k = 1; k < A.NrProcs; ++k) {
				l = A.Pstart[k + 1] - A.Pstart[k];
				
				if (l > MaxNz) MaxNz = l;
				if (l < MinNz) MinNz = l;
			}
			
			/* Calculate communication volume. */
			if (!CalcCom(&A, VAssign, ROW, &ComVolV, &l, &l, &l, &l) ||
			    !CalcCom(&A, UAssign, COL, &ComVolU, &l, &l, &l, &l)) {
				fprintf(stderr, "main(): Unable to calculate communication volume!\n");
				exit(EXIT_FAILURE);
			}
			
			/* Store attributes. */
			Matrices[i].NumNz = A.NrNzElts;
			Matrices[i].Rows = A.m;
			Matrices[i].Cols = A.n;
			
			Matrices[i].Attributes[0].Data[j] = Duration;
			Matrices[i].Attributes[1].Data[j] = (double)(Options.P*MaxNz - A.NrNzElts)/(double)A.NrNzElts;
			Matrices[i].Attributes[2].Data[j] = (double)(MaxComV + MaxComU);
			Matrices[i].Attributes[3].Data[j] = (double)(ComVolV + ComVolU);
			
			/* Free memory. */
			MMDeleteSparseMatrix(&A);
			free(UAssign);
			free(VAssign);
		}
		
		/* Average attributes. */
		if (!AverageAndFreeAttributes(&Matrices[i])) {
			fprintf(stderr, "main(): Cannot setup attributes!\n");
			exit(EXIT_FAILURE);
		}
		
		fprintf(stderr, "\n");
	}
	
	/* Write accumulated data to stdout. */
	fprintf(stderr, "Finished profiling, writing data ...\n");
	
	printf("%% Profiled Mondriaan for %d matrices, %d samples, and %f imbalance.\n", NumMatrices, NumSamples, Options.eps);
	printf("\\documentclass[a4paper, 10pt]{article}\n\n");
	printf("\\usepackage{lscape}\n");
	printf("\\usepackage{longtable}\n\n");
	printf("\\author{\\texttt{Profile.c}}\n");
	CurrentTime = time(NULL);
	printf("\\date{%s}\n", asctime(localtime(&CurrentTime)));
	printf("\\title{Profiling Mondriaan %s with %d %s}\n\n", MONDRIAANVERSION, NumMatrices, NumMatrices > 1 ? "matrices" : "matrix");
	printf("\\begin{document}\n\n");
	
	printf("\\maketitle\n\n");
	
	printf("\\section{Results}\n\n");
	printf("Used Mondriaan version %s to distribute %d matrices (listed in table \\ref{MondriaanMatrices}) over %d processors with maximum imbalance %f, taking the average of %d samples. The used options can be found in table \\ref{MondriaanSettings} en the numerical results in table \\ref{MondriaanResults}.\n", MONDRIAANVERSION, NumMatrices, Options.P, Options.eps, NumSamples);
	printf("This took %.1f minutes in total.\n\n", (double)(clock() - TotalClock)/(60.0*(double)CLOCKS_PER_SEC));
	
	/* Export options. */
	printf("\\begin{table}\n");
	printf("\\caption{Mondriaan configuration.}\n");
	printf("\\label{MondriaanSettings}\n");
	printf("\\begin{center}\n");
	
	if (!ExportOptionsToLaTeX(stdout, &Options)) {
		fprintf(stderr, "main(): Unable to create option table!\n");
		exit(EXIT_FAILURE);
	}
	
	printf("\\end{center}\n");
	printf("\\end{table}\n\n");
	
	/* Export list of test matrices. */
	printf("\\begin{table}\n");
	printf("\\caption{%d tested matrices.}\n", NumMatrices);
	printf("\\label{MondriaanMatrices}\n");
	printf("\\begin{center}\n");
	printf("\\begin{tabular}{l|lll}\nFile & $\\mathit{nz}$ & $m$ & $n$ \\\\\n\\hline\n");
	
	for (i = 0; i < NumMatrices*NumNumProcs; i += NumNumProcs) {
		printf("\\verb|%s| & %ld & %ld & %ld \\\\\n", Matrices[i].File, Matrices[i].NumNz, Matrices[i].Rows, Matrices[i].Cols);
	}
	
	printf("\\hline\n");
	
	printf("\\end{tabular}\n");
	printf("\\end{center}\n");
	printf("\\end{table}\n\n");
	
	/* Export test data. */
	printf("\\begin{landscape}\n");
	printf("\\begin{longtable}{lrrrrr}");
	printf("\\caption[Profile results]{Profile results for %d matrices.}\n", NumMatrices);
	printf("\\label{MondriaanResults} \\\\\n\n");
	
	printf("\\multicolumn{1}{c}{File} & \\multicolumn{1}{c}{$p$} & \\multicolumn{1}{c}{Time (s)} & \\multicolumn{1}{c}{$\\varepsilon$} & \\multicolumn{1}{c}{Max. com.} & \\multicolumn{1}{c}{Com. vol.} \\\\ \\hline\n");
	printf("\\endfirsthead\n\n");
	
	printf("\\multicolumn{6}{c}{\\tablename\\ \\thetable{} -- continued from previous page.} \\\\\n");
	printf("\\multicolumn{1}{c}{File} & \\multicolumn{1}{c}{$p$} & \\multicolumn{1}{c}{Time (s)} & \\multicolumn{1}{c}{$\\varepsilon$} & \\multicolumn{1}{c}{Max. com.} & \\multicolumn{1}{c}{Com. vol.} \\\\ \\hline\n");
	printf("\\endhead\n\n");
	
	printf("\\multicolumn{6}{c}{Continued on next page.} \\\\\n");
	printf("\\endfoot\n\n");
	
	printf("\\hline\n\\endlastfoot\n\n");
	
	for (i = 0; i < NumMatrices*NumNumProcs; i += NumNumProcs) {
		int j;
		
		for (j = 0; j < NumNumProcs; ++j) {
			char Tmp[256];
			const struct sMatrixData *Mat = &Matrices[i + j];
			
			if (j == 0) printf("\\verb|%s| & %ld", Mat->File, Mat->P);
			else printf(" & %ld", Mat->P);
			
			/*
			int k;
			for (k = 0; k < NUM_ATTRIBUTES; ++k) {
				char Tmp[256];
				
				DoubleToLaTeX(Tmp, Mat->Attributes[k].Average, 3);
				printf(" & $%s \\pm ", Tmp);
				DoubleToLaTeX(Tmp, sqrt(Mat->Attributes[k].Variance), 1);
				printf("%s$", Tmp);
			}
			*/
			
			DoubleToLaTeX(Tmp, Mat->Attributes[0].Average, 3);
			printf(" & $%s \\pm ", Tmp);
			DoubleToLaTeX(Tmp, sqrt(Mat->Attributes[0].Variance), 1);
			printf("%s$", Tmp);
			
			DoubleToLaTeX(Tmp, Mat->Attributes[1].Average, 3);
			printf(" & $%s$", Tmp);
			
			printf(" & $%ld \\pm %ld$", (long)Mat->Attributes[2].Average, (long)sqrt(Mat->Attributes[2].Variance));
			printf(" & $%ld \\pm %ld$", (long)Mat->Attributes[3].Average, (long)sqrt(Mat->Attributes[3].Variance));
			
			printf(" \\\\\n");
		}
		printf("\\hline\n");
	}
	
	printf("\n\\end{longtable}\n");
	printf("\\end{landscape}\n\n");
	
	printf("\\end{document}\n\n");
	
	/* Append raw data. */
	printf("Raw data:\n");
	
	for (i = 0; i < NumMatrices*NumNumProcs; ++i) {
		int j;
		
		printf("%s\t%ld\t%ld\t%ld\t%ld\t", Matrices[i].File, Matrices[i].NumNz, Matrices[i].Rows, Matrices[i].Cols, Matrices[i].P);
		
		for (j = 0; j < NUM_ATTRIBUTES - 1; ++j) printf("%e\t%e\t", Matrices[i].Attributes[j].Average, sqrt(Matrices[i].Attributes[j].Variance));
		printf("%e\t%e\n", Matrices[i].Attributes[j].Average, sqrt(Matrices[i].Attributes[j].Variance));
	}
	printf("\n");
	
	free(Matrices);

	fprintf(stderr, "Done!\n");
	
	exit(EXIT_SUCCESS);
}

