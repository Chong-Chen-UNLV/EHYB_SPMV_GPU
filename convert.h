#ifndef CONVERT_H 
#define CONVERT_H 

#include "kernel.h"

extern "C" {

void COO2ELL_block(unsigned int *size_COO, unsigned int* ELL_block_cols_vec, unsigned int* ELL_block_bias_vec,
		unsigned int **colELL, double **matrixELL, unsigned int **I_COO, unsigned int **J_COO, double **V_COO,
		const unsigned int *row_local, const unsigned int *col_local, const double* matrix_local, 
		const unsigned int *row_idx, const unsigned int *numInRow, 
		const unsigned int localMatrixSize, const unsigned int loc_num_of_row);

void COO2ELL(const unsigned int *row_local, const unsigned int *col_local, const double* matrix_local, unsigned int **colELL,
	double **matrixELL, unsigned int **I_COO, unsigned int **J_COO, double **V_COO,const unsigned int *numInRow, 
	const unsigned int *row_idx, const unsigned int localMatrixSize, const unsigned int loc_num_of_row, 
	unsigned int *sizeOut, unsigned max_in, unsigned int *max_out);
}

#endif
