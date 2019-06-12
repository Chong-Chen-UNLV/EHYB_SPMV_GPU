#include "kernel.h"
#include "solver.h"
#include "convert.h"

void COO2ELL(const unsigned int *row_local, const unsigned int *col_local, 
		const double* matrix_local, unsigned int **colELL,
		double **matrixELL, unsigned int **I_COO, unsigned int **J_COO, 
		double **V_COO,const unsigned int *numInRow, 
		const unsigned int *row_idx, const unsigned int localMatrixSize, 
		const unsigned int loc_num_of_row, 
		unsigned int *sizeOut, unsigned max_in, unsigned int *max_out)
{
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
			if (numInRow[i+row_bias] > maxRowNum)
				size_COO+=numInRow[i+row_bias]- maxRowNum;
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
		if (numInRow[i+row_bias]>maxRowNum) {
			for (unsigned int j=0;j<numInRow[i+row_bias];j++){

				//the ELL value should still be set as zero
				if (j<maxRowNum){
					(*colELL)[i+j*loc_num_of_row] = col_local[row_idx[row_val]+j-num_bias];
					(*matrixELL)[i+j*loc_num_of_row]=matrix_local[row_idx[row_val]+j-num_bias];
				}
				else if (j >= maxRowNum){
					//assign COO value
					(*I_COO)[pointCOO]=row_local[row_idx[row_val]+j-num_bias];
					(*J_COO)[pointCOO]=col_local[row_idx[row_val]+j-num_bias];
					(*V_COO)[pointCOO]=matrix_local[row_idx[row_val]+j-num_bias];
					pointCOO++;
					if(pointCOO > size_COO)
						printf("error at pointCOO %d\n", pointCOO);
				}
				else;
			}
			irregular=irregular+1;
		}
		//goto ELL format
		else {
			for (unsigned int j=0;j<maxRowNum;j++){
				//write the ELL data
				if (j<numInRow[i+row_bias]){
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

static void ELL_block_cols_vec_gen_rodr(unsigned int* ELL_block_cols_vec, 
		unsigned int* size_COO, const unsigned int* part_boundary, 
		const int dimension, const int max_col, const int block_num,
	   	const unsigned int* row_idx)
{
	unsigned int avg_local = 0;
	//unsigned int block_nonz = 0;
	unsigned int block_idx = 0;
	unsigned int num_cols = 0;
	unsigned int local_end, local_start;
    unsigned int nonzero;
	unsigned int col_th, col_accum; 
	unsigned int *col_hist = (unsigned int*)calloc(max_col, sizeof(unsigned int));

	*size_COO = 0;
	for(unsigned int block_idx = 0; block_idx < block_num; ++block_idx){
		memset(col_hist, 0, max_col*sizeof(unsigned int));
		col_accum = 0;
		local_start = part_boundary[block_idx];
		local_end = part_boundary[block_idx + 1];
		col_th = ceil( float((local_end - local_start))*0.667);
		for(unsigned int row_val = local_start; row_val < local_end; ++row_val){
            nonzero = (row_idx[row_val + 1] - row_idx[row_val]); 
			if(nonzero > max_col)
				printf("error with large nonzero\n");
			col_hist[nonzero - 1] += 1;
        }
		num_cols = 0;
		for(unsigned int i = 0; i < max_col; ++i){
			col_accum += col_hist[i];
			if(col_accum > col_th){
				num_cols = i + 1;
				break;
			}
		}
		if(num_cols == 0){
		   	printf("error at ELL_block convert\n");
			exit(0);
		}		
   		for(unsigned int row_val = local_start; row_val < local_end; ++row_val){
			unsigned int nonzeros = row_idx[row_val + 1] - row_idx[row_val];
			if(nonzeros > num_cols){
				*size_COO += nonzeros - num_cols;	
			} 
		}
		ELL_block_cols_vec[block_idx] = num_cols;	
	}
	free(col_hist);
}

static void ELL_block_cols_vec_gen(unsigned int* ELL_block_cols_vec, unsigned int* size_COO,
				const int dimension, const int max_col, const unsigned int *row_idx)
{
	unsigned int avg_local = 0;
	//unsigned int block_nonz = 0;
	unsigned int block_idx = 0;
	unsigned int num_cols = 0;
	unsigned int local_end;
    unsigned int nonzero;
	unsigned int col_th, col_accum; 
	unsigned int *col_hist = (unsigned int*)calloc(max_col, sizeof(unsigned int));
	
	*size_COO = 0;

	for(unsigned int row = 0; row < dimension; row += ELL_threadSize){
		memset(col_hist, 0, max_col*sizeof(unsigned int));
		col_accum = 0;
		block_idx = row/ELL_threadSize;
		if (row + ELL_threadSize <= dimension){
			local_end = row + ELL_threadSize;	
			//block_nonz = row_idx[row + ELL_threadSize] - row_idx[row];
			//avg_local = ceil(((float) block_nonz))/ELL_threadSize;
		}else{
			local_end = dimension;
			//block_nonz = row_idx[dimension] - row_idx[row];
			//avg_local = ceil(((float) block_nonz)/(dimension - row));
		}
		col_th = ceil( float((local_end - row))*0.667);
        for(int row_val = row; row_val < local_end; ++row_val){
            nonzero = (row_idx[row_val + 1] - row_idx[row_val]); 
			if(nonzero > max_col)
				printf("error with large nonzero\n");
			col_hist[nonzero - 1] += 1;
        }
		num_cols = 0;
		for(unsigned int i = 0; i < max_col; ++i){
			col_accum += col_hist[i];
			if(col_accum > col_th){
				num_cols = i + 1;
				break;
			}
		}
		if(num_cols == 0){
		   	printf("error at ELL_block convert\n");
			exit(0);
		}		
        /*if(max_col > avg_local*1.2){
			num_cols = avg_local*1.2;
		}
		else{
			num_cols = max_col;
		}*/

		for(int row_val = row; row_val < local_end; ++row_val){
			unsigned int nonzeros = row_idx[row_val + 1] - row_idx[row_val];
			if(nonzeros > num_cols){
				*size_COO += nonzeros - num_cols;	
			} 
		}
		ELL_block_cols_vec[block_idx] = num_cols;	
	}
	free(col_hist);
}

static void COO2ELL_block_core(unsigned int* colELL, double* matrixELL,
		unsigned int* I_COO, unsigned int* J_COO, double* V_COO, const unsigned int size_COO,
		const unsigned int* row_idx, const unsigned int* numInRow, 
		const unsigned int loc_num_of_row, 
        const unsigned int* ELL_block_bias_vec, const unsigned int* ELL_block_cols_vec, 
		const unsigned int* row_local, const unsigned int* col_local, 
		const double* matrix_local,
		unsigned int block_num, const unsigned int* boundary, volatile bool RODR)
{ 
	unsigned int block_rowSize;
	unsigned int irregular=0;
	block_rowSize = ELL_threadSize;
	unsigned int row_bias = row_local[0]; 
	unsigned pointCOO = 0;
	for (unsigned int block_idx = 0; block_idx < block_num; ++block_idx){
		if(RODR){
			block_rowSize = boundary[block_idx + 1] - boundary[block_idx]; 
		}/*else{
			block_rowSize = ELL_threadSize;
		}*/
		unsigned int num_cols = ELL_block_cols_vec[block_idx];
		unsigned int block_bias = ELL_block_bias_vec[block_idx];
		unsigned int num_bias=row_idx[row_bias];
		unsigned int boundary_start, boundary_end;
		if(RODR){
			boundary_start = boundary[block_idx];
			boundary_end = boundary[block_idx + 1];	
		} else{
			boundary_start = block_idx*ELL_threadSize;
			boundary_end = (block_idx + 1)*ELL_threadSize;
		}
		
		for(unsigned int i = 0; i < block_rowSize; ++i){
			unsigned int row_val = i + boundary_start + row_bias;
			if(row_val >= loc_num_of_row) break;
			if (numInRow[row_val] > num_cols) {//goto COO format
				for (unsigned int j = 0; j < numInRow[row_val]; ++j){
					//the ELL value should still be set 
					if (j < num_cols){
						colELL[block_bias + i + j*block_rowSize] = 
							col_local[row_idx[row_val] + j - num_bias];
						matrixELL[block_bias + i + j*block_rowSize]
							= matrix_local[row_idx[row_val] + j - num_bias];
					}
					else{
						//assign COO value
						I_COO[pointCOO]=row_local[row_idx[row_val]+j-num_bias];
						J_COO[pointCOO]=col_local[row_idx[row_val]+j-num_bias];
						V_COO[pointCOO]=matrix_local[row_idx[row_val]+j-num_bias];
						pointCOO++;
						if(pointCOO > size_COO){
							printf("error at pointCOO %d\n", pointCOO);
							exit(0);
						}
					}
				}
				irregular=irregular+1;
			}else {//goto ELL format
				for (unsigned int j=0; j < num_cols; ++j){
					//write the ELL data
					if (j < numInRow[row_val]){
						colELL[block_bias + i + j*block_rowSize] = col_local[row_idx[row_val] + j - num_bias];
						matrixELL[block_bias + i + j*block_rowSize]
								= matrix_local[row_idx[row_val] + j - num_bias];
					}
					//write zero
					else{
						colELL[block_bias + i + j*block_rowSize] = 0;
						matrixELL[block_bias + i + j*block_rowSize] = 0;
					}
				}
			}
		}
	}
}

static void update_ELL_block_bias_vec(unsigned int block_num, unsigned int* ELL_block_cols_vec, unsigned int* ELL_block_bias_vec){
	ELL_block_bias_vec[0] = 0;
	for(unsigned int i = 1; i < block_num; ++i){
		ELL_block_bias_vec[i] = ELL_block_bias_vec[i-1] + ELL_block_cols_vec[i-1]*ELL_threadSize;	
	}	
}

/*parameters, first line: output of HYB related parameters, 2nd line: output of HYB matrices, 3rd line:input of local COO matrix 
4th line: CSR indeces (didn't include values), 5th line: input variables*/

void COO2ELL_block(unsigned int *size_COO, 
		unsigned int* ELL_block_cols_vec, unsigned int* ELL_block_bias_vec,
		unsigned int **colELL, double **matrixELL, unsigned int **I_COO, 
		unsigned int **J_COO, double **V_COO,
		const unsigned int *row_local, const unsigned int *col_local, const double* matrix_local, 
		const unsigned int *row_idx, const unsigned int *numInRow, const unsigned int max_col,
		const unsigned int localMatrixSize, const unsigned int loc_num_of_row, 
		const unsigned int block_num, const unsigned int* boundary, bool RODR){ 

	unsigned int point_COO = 0;
	unsigned int row_bias = row_local[0];
	unsigned int num_bias = row_idx[row_bias];
	*size_COO=0;	
	if(RODR){
		ELL_block_cols_vec_gen_rodr(ELL_block_cols_vec, 
		size_COO, boundary, 
		loc_num_of_row, max_col, block_num,
	   	row_idx);
	} else {
		ELL_block_cols_vec_gen(ELL_block_cols_vec, size_COO,
			loc_num_of_row, max_col, row_idx); 
	}

	if(*size_COO > 0){	
		*I_COO=(unsigned int *)malloc(*size_COO*sizeof(unsigned int));
		*J_COO=(unsigned int *)malloc(*size_COO*sizeof(unsigned int));
		*V_COO=(double *)malloc(*size_COO*sizeof(double));
	}
	unsigned int ELL_matrixSize = 0;
	for(unsigned int i = 0; i < block_num; ++i){
		ELL_block_bias_vec[i] = ELL_matrixSize;
		if(RODR){
			ELL_matrixSize += (boundary[i+1] - boundary[i])*ELL_block_cols_vec[i];
		} else {
			ELL_matrixSize += ELL_threadSize*ELL_block_cols_vec[i];
		}
	}
	ELL_block_bias_vec[block_num] = ELL_matrixSize;
	*colELL=(unsigned int *)malloc(ELL_matrixSize*sizeof(unsigned int));
	*matrixELL=(double *)malloc(ELL_matrixSize*sizeof(double));

	COO2ELL_block_core(*colELL, *matrixELL,
			*I_COO, *J_COO, *V_COO, *size_COO,
			row_idx, numInRow, loc_num_of_row, 
			ELL_block_bias_vec, ELL_block_cols_vec, 
			row_local, col_local, matrix_local, block_num, boundary, RODR);
	
}


