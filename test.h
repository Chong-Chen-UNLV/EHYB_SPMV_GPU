#ifndef TEST_H
#define TEST_H


#include "solver.h"
#include "kernel.h"

void matrix_vectorTest(matrixCOO_S *matrix, const double* vector_in, double* vector_out, 
		int testpoint){
	FILE *tFile;
	if(testpoint > 0) {
		if ((tFile = fopen("testResult", "w")) == NULL){ 
			printf("file open error\n");
			exit(1);
		}
	}
	for(unsigned int i = 0; i < matrix->dimension; ++i) vector_out[i] = 0;
	for(unsigned int i = 0; i < matrix->totalNum; ++i){
		vector_out[matrix->I[i]] += matrix->V[i]*vector_in[matrix->J[i]];	
		
		if(matrix->I[i] == testpoint && testpoint > 0){ fprintf(tFile, 
				"row is %d V is %f vec is %f accum is %f\n", 
				testpoint, matrix->V[i], vector_in[matrix->J[i]], vector_out[matrix->I[i]]);
		}
	}
	if(testpoint > 0)
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
