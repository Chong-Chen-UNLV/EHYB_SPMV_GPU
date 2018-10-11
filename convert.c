#include "kernel.h"
#include "solver.h"

void COO2ELL(const unsigned int *row_local, const unsigned int *col_local, const double* matrix_local, unsigned int **colELL,
	double **matrixELL, unsigned int **I_COO, unsigned int **J_COO, double **V_COO,const unsigned int *nonzero_vec, 
	const unsigned int *row_idx, const unsigned int localMatrixSize, const unsigned int loc_num_of_row, 
	unsigned int *sizeOut, unsigned max_in, unsigned int *max_out){

	unsigned int maxRowNum;

	unsigned int size_COO=0, pointCOO=0;
	unsigned int row_bias=row_local[0];
	unsigned int num_bias=row_idx[row_bias];
	
	maxRowNum=ceil(((double) localMatrixSize)*1.2/((double) loc_num_of_row));
	if (maxRowNum > max_in){
		maxRowNum = max_in;	
	}
	else{
		for (unsigned int i=0;i<loc_num_of_row;i++){
			unsigned int row_idx=i+row_local[0];
			if (nonzero_vec[i+row_bias] > maxRowNum)
				size_COO+=nonzero_vec[i+row_bias]- maxRowNum;
		}	
	}
	if(size_COO > 0){	
		*I_COO=(unsigned int *)malloc(size_COO*sizeof(unsigned int));
		*J_COO=(unsigned int *)malloc(size_COO*sizeof(unsigned int));
		*V_COO=(double *)malloc(size_COO*sizeof(double));
	}
	*colELL=(unsigned int *)malloc(loc_num_of_row*maxRowNum*sizeof(unsigned int));
	*matrixELL=(double *)malloc(loc_num_of_row*maxRowNum*sizeof(double));
	
	unsigned int irregular=0;
	for (unsigned int i=0;i <loc_num_of_row; i++){
		unsigned int row_val = i+row_local[0];
		//goto COO format
		if (nonzero_vec[i+row_bias]>maxRowNum) {
			for (unsigned int j=0;j<nonzero_vec[i+row_bias];j++){

				//the ELL value should still be set as zero
				if (j<maxRowNum){
					(*colELL)[i+j*loc_num_of_row] = col_local[row_idx[row_val]+j-num_bias];
					(*matrixELL)[i+j*loc_num_of_row]=matrix_local[row_idx[row_val]+j-num_bias];
				}
				else{
					//assign COO value
					(*I_COO)[pointCOO]=row_local[row_idx[row_val]+j-num_bias];
					(*J_COO)[pointCOO]=col_local[row_idx[row_val]+j-num_bias];
					(*V_COO)[pointCOO]=matrix_local[row_idx[row_val]+j-num_bias];
					pointCOO++;
					if(pointCOO > size_COO)
						printf("error at pointCOO %d\n", pointCOO);
				}
			}
			irregular=irregular+1;
		}
		//goto ELL format
		else {
			for (unsigned int j=0;j<maxRowNum;j++){
				//write the ELL data
				if (j<nonzero_vec[i+row_bias]){
					(*colELL)[i+j*loc_num_of_row]=col_local[row_idx[row_val]+j-num_bias];
					(*matrixELL)[i+j*loc_num_of_row]=matrix_local[row_idx[row_val]+j-num_bias];
				}
				//write zero
				else{
					(*colELL)[i+j*loc_num_of_row]=0;
					(*matrixELL)[i+j*loc_num_of_row]=0;
				}
			}
		}
	}
	*max_out = maxRowNum;
	*sizeOut = size_COO;
	//printf("irregular row Number is %d\n",irregular);
	
}

static void ELL_block_wid_vec_gen(unsigned int* num_cols_vec, unsigned int* size_COO,
				const int* I, const int* J,
				const int dimension, const unsigned int *row_idx,
				unsigned int* ELL_block_wid_vec){
	
	unsigned int max_col = 0;
	unsigned int avg_local = 0;
	unsigned int block_nonz = 0;
	unsigned int block_idx = 0;
	unsigned int num_cols = 0;
	unsigned int local_end;
	*size_COO = 0;
	int block_nonZ;
	unsigned int ELL_blocks = ceil(((float) (dimension))/ELL_rodr_thread_size);
	for(unsigned int row = 0; row < dimension; row += ELL_blocks){
		block_idx = row/ELL_rodr_thread_size;
		if (row + ELL_rodr_thread_size <= dimension){
			local_end = row + ELL_rodr_thread_size;	
			block_nonZ = row_idx[row + ELL_rodr_thread_size] - row_idx[row];
			if(block_nonz > max_col) max_col = block_nonz;
			avg_local = ceil(((float) block_nonz))/ELL_rodr_thread_size;
		}else{
			local_end = dimension;
			block_nonz = row_idx[dimension] - row_idx[row];
			avg_local = ceil(((float) block_nonz)/(dimension - row));
		}
		if(max_col > avg_local*1.2){
			num_cols = avg_local*1.2;
		}
		else{
			num_cols = max_col;
		}

		for(int row_val = row; row_val < local_end; ++row_val){
			unsigned int nonzeros = row_idx[row_val + 1] - row_idx[row_val];
			if(nonzeros > num_cols){
				*size_COO += nonzeros - num_cols;	
			} 
		}
		ELL_block_wid_vec[block_idx] = num_cols;	
	}
}

static void COO2ELL_block_core(int* colELL, double* matrixELL,
		int* I_COO, int* J_COO, double* V_COO, 
		unsigned int* num_cols_vec, 
		const int* row_local, const int* col_local, const double* matrix_local){ 

	unsigned int irregular=0;
	unsigned int ELL_blocks = ceil( ((double)row_local)/ELL_threadSize);
	unsigned int *ELL_bias_vec = (unsigned int*) malloc(sizeof(unsigned int)*ELL_blocks);
	ELL_bias_vec[0] = 0;
	for (unsigned int blk = 1; blk < ELL_block_size; ++blk ){
		ELL_bias_vec[blk] = num_cols_vec[blk]*ELL_threadSize + ELL_bias_vec[blk - 1]
	}
	for (unsigned int row = 0; row < loc_num_of_row; row += ELL_blockSize){
		unsigned int block_idx = row/ELL_blockSize;
		unsigned int num_cols = num_cols_vec[block_idx];
		unsigned int ELL_bias = block_ELL_bias_vec[block_idx]; 
		
		for(i = row; i < row + ELL_blockSize; ++i){
			unsigned int row_val = i + row_bias;
			if (nonzero_vec[i+row_bias] > num_cols) {//goto COO format
				for (unsigned int j = 0; j < nonzero_vec[row_val + row_bias]; ++j){
					//the ELL value should still be set 
					if (j < num_cols){
						(*colELL)[ELL_bias + j*ELL_treadSize] = col_local[row_idx[row_val]+j-num_bias];
						(*matrixELL)[ELL_bias + j*ELL_treadSize]
							= matrix_local[row_idx[row_val]+j-num_bias];
					}
					else{
						//assign COO value
						(*I_COO)[pointCOO]=row_local[row_idx[row_val]+j-num_bias];
						(*J_COO)[pointCOO]=col_local[row_idx[row_val]+j-num_bias];
						(*V_COO)[pointCOO]=matrix_local[row_idx[row_val]+j-num_bias];
						pointCOO++;
						if(pointCOO > size_COO)
							printf("error at pointCOO %d\n", pointCOO);
					}
				}
				irregular=irregular+1;
			}else {//goto ELL format
				for (unsigned int j=0; j < num_cols; ++j){
					//write the ELL data
					if (j<nonzero_vec[i+row_bias]){
						(*colELL)[ELL_bias + j*ELL_treadSize] = col_local[row_idx[row_val]+j-num_bias];
						(*matrixELL)[ELL_bias + j*ELL_treadSize]
								= matrix_local[row_idx[row_val]+j-num_bias];
					}
					//write zero
					else{
						(*colELL)[ELL_bias + j*ELL_treadSize] = 0;
						(*matrixELL)[ELL_bias + j*ELL_treadSize] = 0;
					}
				}
			}
		}
	}
	free(ELL_bias_vec);
}

static void update_ELL_block_bias_vec(unsigned int block_num, unsigned int* ELL_block_wid_vec, unsigned int* ELL_block_bias_vec){
	for(unsigned int i = 0; i < block_num; ++i){

	}	
}

void COO2ELL_block(unsigned int* num_cols_vec, unsigned int *size_COO, unsigned int* block_ELL_bias_vec, 
		double **matrixELL, unsigned int **I_COO, unsigned int **J_COO, double **V_COO,
		const unsigned int *row_local, const unsigned int *col_local, const double* matrix_local, unsigned int **colELL,
		const unsigned int *nonzero_vec, 
		const unsigned int *row_idx, const unsigned int localMatrixSize, const unsigned int loc_num_of_row, 
		unsigned int *sizeOut, unsigned max_in, unsigned int *max_out){


	unsigned int point_COO = 0;
	unsigned int row_bias = row_local[0];
	unsigned int num_bias = row_num_accum[row_bias];
	unsigned int block_num = ceil(( float(loc_num_of_row))/ELL_block_size);
	*size_COO=0;	
	unsigned int* ELL_block_wid_vec = (unsigned int*) malloc(block_num*sizeof(unsigned int));
	ELL_block_wid_vec_gen(num_cols_vec, size_COO,
			row_local, col_local,
			loc_num_of_row, row_idx, ELL_block_wid_vec); 

	if(size_COO > 0){	
		*I_COO=(unsigned int *)malloc(size_COO*sizeof(unsigned int));
		*J_COO=(unsigned int *)malloc(size_COO*sizeof(unsigned int));
		*V_COO=(double *)malloc(size_COO*sizeof(double));
	}
	*colELL=(unsigned int *)malloc(loc_num_of_row*maxRowNum*sizeof(unsigned int));
	*matrixELL=(double *)malloc(loc_num_of_row*maxRowNum*sizeof(double));

	row_idx += ELL_block_size;
	COO2ELL_block_core(colELL, matrixELL,
			I_COO, J_COO, V_COO, 
			num_cols_vec, block_ELL_bias, 
			row_local, col_local, matrix_local);
	
	if (maxRowNum > max_in){
		maxRowNum = max_in;	
	}else{
		for (unsigned int i=0;i<loc_num_of_row;i++){
			unsigned int row_idx=i+row_local[0];
			if (nonzero_vec[i+row_bias] > maxRowNum)
				size_COO+=nonzero_vec[i+row_bias]- maxRowNum;
		}	
	}
	
	
	
	
	*max_out = maxRowNum;
	*sizeOut = size_COO;
	//printf("irregular row Number is %d\n",irregular);
	
}
