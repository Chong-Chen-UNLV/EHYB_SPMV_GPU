#ifndef CONVERT_H 
#define CONVERT_H 
#include "solver.h"
#include "kernel.h"

static inline unsigned int get_blocks_rodr(const unsigned int* boundary, 
		const unsigned int block_size){
	unsigned int block_num = 0;
	for(unsigned int i = 0; i < block_size; ++i){
		unsigned int partRows = boundary[i + 1] - boundary[i];
		block_num += ceil((float) partRows/ELL_threadSize);	
	}	
	return block_num;
}



void COO2ELL_block(unsigned int *size_COO, 
		unsigned int* ELL_block_cols_vec, 
		unsigned int* ELL_block_bias_vec,
		unsigned int **colELL, 
		double **matrixELL, 
		unsigned int **I_COO, 
		unsigned int **J_COO, 
		double **V_COO,
		const unsigned int *row_local, 
		const unsigned int *col_local, 
		const double* matrix_local, 
		const unsigned int *row_idx, 
		const unsigned int *numInRow, 
		const unsigned int max_col, 
		const unsigned int localMatrixSize, 
		const unsigned int loc_num_of_row, 
		const unsigned int part_size,
		const unsigned int block_num, 
		const unsigned int* boundary, 
		bool RODR,
		bool CACHE);

void COO2ELL(const unsigned int *row_local, const unsigned int *col_local, const double* matrix_local, 
		unsigned int **colELL,
	double **matrixELL, unsigned int **I_COO, unsigned int **J_COO, double **V_COO,const unsigned int *numInRow, 
	const unsigned int *row_idx, const unsigned int localMatrixSize, const unsigned int loc_num_of_row, 
	unsigned int *sizeOut, unsigned max_in, unsigned int *max_out);

//void matrix_vectorHYB(struct abc inputMatrix);
#endif
