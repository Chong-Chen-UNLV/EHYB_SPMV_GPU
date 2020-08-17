#ifndef TEST_H
#define TEST_H


#include "solver.h"
#include "kernel.h"

void mimicHYB(int* ELL_block_cols_vec,
		int* ELL_block_bias_vec,
		const int* part_boundary,
		int* colELL,
		double* matrixELL,
		int* I_COO,
		int* J_COO,
		double* V_COO,
		int size_COO,
		const double* vector_in,
		const int testpoint){

	FILE *tFile;
	double shared_data[vector_cache_size ];
	memset(shared_data, 0, sizeof(double)*vector_cache_size);
	if(testpoint > 0){
		if ((tFile = fopen("HYBResult", "w")) == NULL){ 
			printf("file open error\n");
			exit(1);
		}
		int part_idx = 0;
		int block_rowSize = 0;
		while(part_boundary[part_idx + 1] < testpoint) part_idx++;
		int block_idx = part_idx*block_per_part + (testpoint - part_boundary[part_idx])/ELL_threadSize; 
		int part_start = part_boundary[part_idx];
		int part_end = part_boundary[part_idx + 1];
		int boundary_start = part_boundary[part_idx] + (block_idx%block_per_part)*ELL_threadSize;
		int block_end = part_boundary[part_idx] + block_per_part*ELL_threadSize;	
		int row_local = testpoint-(block_idx%block_per_part)*ELL_threadSize - part_start;
		int block_data_bias = ELL_block_bias_vec[block_idx];
		int boundary_end;
		
		if((block_idx + 1)%block_per_part == 0){
			if(part_end > block_end){
				boundary_end = block_end;
			}
			else boundary_end = part_end;
		} else {
			boundary_end = boundary_start + ELL_threadSize;
		}
		block_rowSize = boundary_end - boundary_start;
		for(int i = part_start; i < part_start + vector_cache_size; ++ i){
			shared_data[i - part_start] = vector_in[i];	
		}
		int col_num = ELL_block_cols_vec[block_idx];
		double accum = 0;
		for(int i = 0; i < col_num; ++i){
			int data_idx = block_data_bias + block_rowSize*i + row_local;
			double data = matrixELL[data_idx];
			int col = colELL[data_idx];
			double vec = shared_data[col - part_start];
			accum += data*vec;
			fprintf(tFile, "row is %d V is %f vec is %f accum is %f\n", 
					testpoint, data, vec, accum);
		}
		fprintf(tFile, "-------start COO format------------\n" );	
					
		for(int i = 0; i < size_COO; ++i){
			if(I_COO[i] == testpoint){
				double vec = vector_in[J_COO[i]];
				accum += V_COO[i]*vec;
				fprintf(tFile, "row is %d V is %f vec is %f accum is %f\n", 
					testpoint, V_COO[i], vec, accum);
			}
		}
		fclose(tFile);
	}
}
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
