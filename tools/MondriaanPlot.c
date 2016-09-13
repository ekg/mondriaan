/*
MondriaanPlot.c

Created by Bas Fagginger Auer.

This program draws a series of pictures which illustrate the Mondriaan splitting process.

NOTA BENE: Symmetry settings are currently discarded!
*/
#include <stdlib.h>
#include <stdio.h>
#include <memory.h>
#include <math.h>
#include <Mondriaan.h>

struct rgb {
	float r;
	float g;
	float b;
};

/* Available colours for the image, colours are selected randomly at each split. */
/* The first colour is the default background colour, the format is {blue0, green0, red0, blue1, green1, red1, ...}. */
const unsigned char Colours[] = {0xff, 0xff, 0xff,
0xff, 0x00, 0x00,
0x00, 0xff, 0x00,
0x00, 0x00, 0xff,
0xff, 0xff, 0x00,
0xff, 0x00, 0xff,
0x00, 0xff, 0xff};
/* The total number of colours in the Colours[] array. */
const int NumColours = 7;
int CurrentColour = 0;
int Symmetric = FALSE;

/* Colour assigned to each processor. */
const unsigned char **PColours;

/* Width and height of the made snapshots. */
const int sWidth = 512, sHeight = 512;
/* RGB pixel data of the snapshot. */
unsigned char *sData;
/* For compatibility, we write shorts in the little-endian format. */
int IsLittleEndian = TRUE;

size_t fwriteushort(const unsigned short Value, FILE *Out) {
	unsigned short LEValue = Value;
	
	if (Out == NULL) return 0;
	
	if (!IsLittleEndian) {
		const unsigned char Lower = Value & 255, Upper = (Value >> 8) & 255;
		
		LEValue = (Lower << 8) | Upper;
	}
	
	return fwrite(&LEValue, sizeof(unsigned short), 1, Out);
}

/* Writes a 24-bit (RGBRGBRGB...) Truevision TGA file to a given file pointer. */
int WriteTGA(const char *Filename, const int Width, const int Height, const unsigned char *Data) {
	unsigned char IDLength = 0, ColourMapType = 0, ImageType = 2, ColourMapEntrySize = 24, BitsPerPixel = 24, ImageDescriptor = 0;
	unsigned short ColourMapOffset = 0, ColourMapLength = 0;
	unsigned short XOrig = 0, YOrig = 0, XSize = Width, YSize = Height;
	FILE *Out;
	
	if (!Filename || Width <= 0 || Height <= 0 || !Data) {
		fprintf(stderr, "WriteTGA(): Invalid parameters!\n");
		return FALSE;
	}
	
	if ((Out = fopen(Filename, "wb")) == NULL) {
		fprintf(stderr, "WriteTGA(): Unable to open '%s' for writing!\n", Filename);
		return FALSE;
	}
	
	fwrite(&IDLength, sizeof(unsigned char), 1, Out);
	fwrite(&ColourMapType, sizeof(unsigned char), 1, Out);
	fwrite(&ImageType, sizeof(unsigned char), 1, Out);
	
	fwriteushort(ColourMapOffset, Out);
	fwriteushort(ColourMapLength, Out);
	fwrite(&ColourMapEntrySize, sizeof(unsigned char), 1, Out);
	
	fwriteushort(XOrig, Out);
	fwriteushort(YOrig, Out);
	fwriteushort(XSize, Out);
	fwriteushort(YSize, Out);
	fwrite(&BitsPerPixel, sizeof(unsigned char), 1, Out);
	fwrite(&ImageDescriptor, sizeof(unsigned char), 1, Out);
	
	fwrite(Data, 3*sizeof(unsigned char), Width*Height, Out);
	
	fclose(Out);
	
	return TRUE;
}

/* Callback function called by DistributeMat() for each split. */
int Callback(int Iteration, int LastSplit, const struct sparsematrix *A) {
	char Filename[] = "img0000.tga";
	long i, j;
	const float DivX = (float)sWidth/(float)A->n, DivY = (float)sHeight/(float)A->m;
	
	printf("Writing image %d... ", Iteration);
	fflush(stdout);
	
	/* Add colour for the most recently added processor. */
	if (LastSplit + 1 < A->NrProcs) {
		memmove(&PColours[LastSplit + 1], &PColours[LastSplit], (A->NrProcs - 1 - LastSplit)*sizeof(char *));
	}
	
	if (++CurrentColour >= NumColours) CurrentColour = 1;
	
	/* Try to avoid using the same colour in adjacent processors. */
	if (PColours[LastSplit + 1] == &Colours[3*CurrentColour]) {
		if (++CurrentColour >= NumColours) CurrentColour = 1;
	}
	
	PColours[LastSplit] = &Colours[3*CurrentColour];
	
	/* Clear initial data buffer to background colour. */
	if (TRUE) {
		unsigned char *cp = sData;
		int x, y;
		
		for (y = 0; y < sHeight; y++) {
			for (x = 0; x < sWidth; x++) {
				*cp++ = Colours[0];
				*cp++ = Colours[1];
				*cp++ = Colours[2];
			}
		}
	}
	
	if (DivX > 1.0 || DivY > 1.0) {
		/* Blocks per processor are potentially bigger than one pixel. */
		for (i = 0; i < A->NrProcs; i++) {
			const unsigned char *Colour = PColours[i];
			
			for (j = A->Pstart[i]; j < A->Pstart[i + 1]; j++) {
				float Row = (A->row_perm_inv != NULL ? A->row_perm_inv[A->i[j]] : A->i[j]);
				float Col = (A->col_perm_inv != NULL ? A->col_perm_inv[A->j[j]] : A->j[j]);
				int x0 = (int)floor(DivX*Col), y0 = sHeight - (int)floor(DivY*(Row + 1));
				int x1 = (int)floor(DivX*(Col + 1.0)), y1 = sHeight - (int)floor(DivY*Row);
				int x, y;
				unsigned char *cp = &sData[3*(x0 + y0*sWidth)];
				
				if (x1 >= sWidth) x1 = sWidth - 1;
				if (y1 >= sHeight) y1 = sHeight - 1;
				
				for (y = y0; y <= y1; y++, cp += 3*(sWidth + x0 - x1 - 1)) {
					for (x = x0; x <= x1; x++) {
						*cp++ = Colour[0];
						*cp++ = Colour[1];
						*cp++ = Colour[2];
					}
				}
				
				if (Symmetric)
				{
					Col = (A->row_perm_inv != NULL ? A->row_perm_inv[A->i[j]] : A->i[j]);
					Row = (A->col_perm_inv != NULL ? A->col_perm_inv[A->j[j]] : A->j[j]);
					x0 = (int)floor(DivX*Col), y0 = sHeight - (int)floor(DivY*(Row + 1));
					x1 = (int)floor(DivX*(Col + 1.0)), y1 = sHeight - (int)floor(DivY*Row);
					cp = &sData[3*(x0 + y0*sWidth)];
					
					if (x1 >= sWidth) x1 = sWidth - 1;
					if (y1 >= sHeight) y1 = sHeight - 1;
					
					for (y = y0; y <= y1; y++, cp += 3*(sWidth + x0 - x1 - 1)) {
						for (x = x0; x <= x1; x++) {
							*cp++ = Colour[0];
							*cp++ = Colour[1];
							*cp++ = Colour[2];
						}
					}
				}
			}
		}
	}
	else {
		/* Blocks per processor are always smaller than or equal to one pixel. */
		for (i = 0; i < A->NrProcs; i++) {
			const unsigned char *Colour = PColours[i];
			
			for (j = A->Pstart[i]; j < A->Pstart[i + 1]; j++) {
				float Row = (A->row_perm_inv != NULL ? A->row_perm_inv[A->i[j]] : A->i[j]);
				float Col = (A->col_perm_inv != NULL ? A->col_perm_inv[A->j[j]] : A->j[j]);
				int x0 = (int)floor(DivX*Col), y0 = sHeight - 1 - (int)floor(DivY*Row);
				unsigned char *cp = &sData[3*(x0 + y0*sWidth)];
				
				*cp++ = Colour[0];
				*cp++ = Colour[1];
				*cp++ = Colour[2];
				
				if (Symmetric)
				{
					Col = (A->row_perm_inv != NULL ? A->row_perm_inv[A->i[j]] : A->i[j]);
					Row = (A->col_perm_inv != NULL ? A->col_perm_inv[A->j[j]] : A->j[j]);
					x0 = (int)floor(DivX*Col), y0 = sHeight - 1 - (int)floor(DivY*Row);
					cp = &sData[3*(x0 + y0*sWidth)];
					
					*cp++ = Colour[0];
					*cp++ = Colour[1];
					*cp++ = Colour[2];
				}
			}
		}
		
	}
	
	/* Write snapshot to disk as TGA image file. */
	sprintf(Filename, "img%04d.tga", Iteration);
	
	if (!WriteTGA(Filename, sWidth, sHeight, sData)) {
		fprintf(stderr, "Callback(): Unable to write TGA '%s'!\n", Filename);
                return FALSE;
	}
	
	/* Write two extra frames at the end. */
	if (Iteration == A->NrProcs) {
		sprintf(Filename, "img%04d.tga", Iteration + 1);
		WriteTGA(Filename, sWidth, sHeight, sData);
		sprintf(Filename, "img%04d.tga", Iteration + 2);
		WriteTGA(Filename, sWidth, sHeight, sData);
	}
	
	printf("done.\n");
        
        return TRUE;
}

int main(int argc, char **argv) {
	FILE *File;
	struct opts Options;
	struct sparsematrix A;
	int i;
	
	/* Check endian-ness. This may trigger a compiler warning. */
	if (TRUE)
	{
		unsigned char Test[2] = {1, 0};
		
		if (*((short *)Test) != 1) {
			fprintf(stderr, "main(): This is a big-endian system.\n");
			IsLittleEndian = FALSE;
		}
		else {
			fprintf(stderr, "main(): This is a little-endian system.\n");
			IsLittleEndian = TRUE;
		}
	}
	
	/* Apply Mondriaan options. */
	SetDefaultOptions(&Options);
	
	if (!SetOptionsFromFile(&Options, "Mondriaan.defaults")) {
		fprintf(stderr, "main(): warning, cannot set options from 'Mondriaan.defaults', using default options!\n");
	}
	
	if (!GetParameters(&Options, argc, argv)) {
		fprintf(stderr, "main(): invalid command line parameters!\n");
		exit(-1);
	}
	
	if (!ApplyOptions(&Options)) {
		fprintf(stderr, "main(): could not apply given options!\n");
		exit(-1);
	}
	
	/* Initialise image data buffer. */
	sData = (unsigned char *)malloc(3*sWidth*sHeight*sizeof(unsigned char));
	
	if (sData == NULL) {
		fprintf(stderr, "main(): not enough memory!\n");
		exit(-1);
	}
	else {
		/* Clear initial data buffer to background colour. */
		unsigned char *cp = sData;
		int x, y;
		
		for (y = 0; y < sHeight; y++) {
			for (x = 0; x < sWidth; x++) {
				*cp++ = Colours[0];
				*cp++ = Colours[1];
				*cp++ = Colours[2];
			}
		}
	}
	
	/* Read matrix from disk. */
	File = fopen(Options.matrix, "r");
	
	if (!File) {
		fprintf(stderr, "main(): Could not open '%s' for reading!\n", Options.matrix);
		exit(-1);
	}
		
	if (!MMReadSparseMatrix(File, &A)) {
		fprintf(stderr, "main(): Could not read matrix!\n");
		exit(-1);
	}
	
	fclose(File);
	
	/* Remove double zeroes. */
	if (!SparseMatrixRemoveDuplicates(&A)) exit(-1);
	
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
			
	/* Set up parameters of A. */
	A.NrProcs = Options.P;
	A.Pstart = (long *)malloc((A.NrProcs + 1)*sizeof(long));
	
	if (A.Pstart == NULL) {
		fprintf(stderr, "main(): Cannot allocate processor array!\n");
		exit(-1);
	}
	
	PColours = (const unsigned char **)malloc((A.NrProcs + 1)*sizeof(const unsigned char *));
	
	if (PColours == NULL) {
		fprintf(stderr, "main(): Not enough memory!\n");
		exit(-1);
	}
	
	/* Initialise processor array. */
	A.Pstart[0] = 0;
	
	for (i = 1; i <= A.NrProcs; i++) {
		A.Pstart[i] = A.NrNzElts;
		PColours[i] = &Colours[3];
	}
	
	/* Write first frame. */
	Callback(1, 0, &A);
	
	/* Distribute the processors among the matrix entries. */
	DistributeMatrixMondriaan (&A, A.NrProcs, Options.eps, &Options, Callback);
	
	free(sData);
	free(PColours);
	
	/* Remove dummies.
	if (A.m == A.n && Options.SquareMatrix_DistributeVectorsEqual == EqVecYes && Options.SquareMatrix_DistributeVectorsEqual_AddDummies == DumYes) RemoveDummiesFromSparseMatrix(&A);
	*/
	
	/* Free memory. */
	MMDeleteSparseMatrix(&A);
	
	exit(0);
}

