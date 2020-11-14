#ifndef TEST_H
#define TEST_H


#include "solver.h"
#include "kernel.h"

void matrix_vectorTestEHYB(matrixEHYB *inputMatrix, const double* vector_in, double* vector_out,
		int testpoint){
	FILE *tFile;
	int track = 0;	
	if(testpoint >= 0) {
		if ((tFile = fopen("testResultEHYB", "w")) == NULL){ 
			printf("file open error\n");
			exit(1);
		}
	}
	for(int part = 0; part < inputMatrix->nParts; ++part){
		int vecStart = inputMatrix->partBoundary[part];
		int vecEnd = vecStart + vectorCacheSize;
		int blockPartStart = part*blockPerPart;
		for(int block = 0; block < blockPerPart; ++block){
			int blockBias = inputMatrix->biasVecBlockELL[blockPartStart + block];
			int width = inputMatrix->widthVecBlockELL[blockPerPart + block];
			int row = vecStart + block*warpSize;
			for(int i = 0; i < 32; ++i){
				if(row + i < vecEnd && row + i < inputMatrix->dimension){
					for(int n = 0; n < width; ++n){
						int dataIdx = blockBias + i + n*warpSize;
						int col = inputMatrix->colBlockELL[dataIdx];
						double val = inputMatrix->valBlockELL[dataIdx];
						vector_out[row + i] += val*vector_in[col];
						if(row + i == testpoint && testpoint >= 0){ fprintf(tFile, 
								"row is %d V is %f vec is %f accum is %f\n", 
								testpoint, val, vector_in[col], vector_out[row + i]);
						}
					}
				}
			}
		}
	}

	int dataIdx, col;
	double val;
	for(int i = 0; i < inputMatrix->numOfRowER; ++i){
		int width = inputMatrix->widthVecER[i/warpSize];
		int bias = inputMatrix->biasVecER[i/warpSize];
		int row = inputMatrix->rowVecER[i];
		for(int n = 0; n < width; ++n){
			dataIdx = bias + i - ((i >> 5)<<5) + n*warpSize; 
			col = inputMatrix->colER[dataIdx];
			val = inputMatrix->valER[dataIdx]; 
			vector_out[row] += val*vector_in[col];
			if(row  == testpoint && testpoint >= 0){ fprintf(tFile, 
					"row is %d V is %f vec is %f accum is %f\n", 
					testpoint, val, vector_in[col], vector_out[row]);
			}
		}
	}
	if(testpoint >= 0)
		fclose(tFile);
}

void matrix_vectorTest(matrixCOO *matrix, const double* vector_in, double* vector_out, 
		int testpoint){
	FILE *tFile;
	if(testpoint >= 0) {
		if ((tFile = fopen("testResult", "w")) == NULL){ 
			printf("file open error\n");
			exit(1);
		}
	}
	for(unsigned int i = 0; i < matrix->dimension; ++i) vector_out[i] = 0;
	for(unsigned int i = 0; i < matrix->totalNum; ++i){
		vector_out[matrix->I[i]] += matrix->V[i]*vector_in[matrix->J[i]];	
		if(matrix->I[i] == testpoint && testpoint >= 0){ fprintf(tFile, 
				"row is %d V is %f vec is %f accum is %f\n", 
				testpoint, matrix->V[i], vector_in[matrix->J[i]], vector_out[matrix->I[i]]);
		}
		//if(vector_out[matrix->I[i]] > 0 && vector_out[matrix->I[i]] < 1)
		//	printf("strange vector out\n");

	}
	if(testpoint >= 0)
		fclose(tFile);
}



/*void calUnbalance(unsigned int part_boundary,
		unsigned int * numInRow, unsigned int * numInRowL, unsigned int * numInRowLP, 
	  	unsigned int nparts,
		int* numOfLines, int* numOfLinesL, numOfLinesLP,
		int* extraCOO, int* extraCOO_L, int* extraCOO_LP){
	*numOfLines = 0;	
	*extraCOO_LP = 0;
	*extraCOO_L = 0;
	*extraCOO = 0;
	for(unsigned int i = 0; i < nparts - 1; ++i){
		numOfLines += abs((int)(part_boundary[i+1] - part_boundary[i]) - shared_per_block);
		if((part_boundary[i+1] - part_boundary[i]) > shared_per_block){
			for(unsigned int j = shared_per_block+part_boundary[i]; j < part_boundary[i+1]; ++j){
				extraCOO += numInRow[j]; 
				extraCOO_L += numInRow_L[j]; 
				extraCOO_LP += numInRow_LP[j]; 
			}	
		}
	}

}*/
#endif
