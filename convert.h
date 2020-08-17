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



void COO2ELL_block(int *size_COO, 
		int* ELL_block_cols_vec, 
		int* ELL_block_bias_vec,
		int **colELL, 
		double **matrixELL, 
		int **I_COO, 
		int **J_COO, 
		double **V_COO,
		const int *row_local, 
		const int *col_local, 
		const double* matrix_local, 
		const int *row_idx, 
		const int *numInRow, 
		const int max_col, 
		const int localMatrixSize, 
		const int loc_num_of_row, 
		const int part_size,
		const int block_num, 
		const int* boundary, 
		bool RODR,
		bool CACHE);

void COO2ELL(const int *row_local, const int *col_local, const double* matrix_local, 
		int **colELL,
	double **matrixELL, int **I_COO, int **J_COO, double **V_COO,const int *numInRow, 
	const int *row_idx, const int localMatrixSize, const int loc_num_of_row, 
	int *sizeOut, max_in, int *max_out);

//void matrix_vectorHYB(struct abc inputMatrix);
#endif
