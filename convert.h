#ifndef CONVERT_H 
#define CONVERT_H 
#include "solver.h"
#include "kernel.h"

static inline int get_blocks_rodr(const int* boundary, 
		const int block_size){
	int block_num = 0;
	for(int i = 0; i < block_size; ++i){
		int partRows = boundary[i + 1] - boundary[i];
		block_num += ceil((float) partRows/ELL_threadSize);	
	}	
	return block_num;
}


void COO2EHYB(matrixCOO* inputMatrix, 
		matrixEHYB* outputMatrix,
		int* sizeBlockELL, 
		int* sizeER);

//void matrix_vectorHYB(struct abc inputMatrix);
#endif
