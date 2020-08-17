#include "kernel.h"
#include "solver.h"
#include "convert.h"
#include "Partition.h"

static void sortRordrListFull(unsigned int dimension,
		int* rodr_list, 
		int* numInRow){

	rowS* rodrSVec = (rowS*) malloc(dimension*sizeof(rowS));
	for(unsigned int i = 0; i < dimension; ++i){
		//the location of elements at SVec will be changed 
		//record the original idx since it will be used for
		//new rodr_list
		rodrSVec[rodr_list[i]].idx = i; 
		rodrSVec[rodr_list[i]].nonzeros = numInRow[i];	
	}	
	for(unsigned int i = 0; i < nparts; ++i){
		qsort(rodrSVec, dimension, sizeof(rowS), rowSCompare);	
	}	
	for(unsigned int i = 0; i < dimension; ++i){
		rodr_list[rodrSVec[i].idx] = i;
	}
	free(rodrSVec);

}
static void colsVecGenEHYB( int* widthVecBlockELL, 
		int* sizeCOO, 
		int* numOfRowER,
		int* rowVecER,
		int* widthVecER,
		const int* col,
		int* blockBoundary, 
		const int* partBoundary, 
		const int dimension, 
		const int maxCol, 
		const int partNum,
		const int* rowIdx,
		const bool CACHED)
{
	/*Reoder functions did an extra iteration 
 	* which provides numInRowER without ER reordering
 	* */
	int avg_local = 0;
	// int block_nonz = 0;
	int block_idx = 0;
	int num_cols = 0;
	int local_end, local_start;
	int nonzero;
	int col_th, col_accum; 
	int *col_hist = ( int*)calloc(max_col, sizeof( int));

	*sizeCOO = 0;
	block_boundary[0] = 0;
	int ELL_val = 0;
	int extraRows = 0;

	for( int part_idx = 0; part_idx < part_num; ++part_idx){
		int block_end = part_boundary[part_idx] + ELL_threadSize*block_per_part;
		int part_end = part_boundary[part_idx + 1];
		int part_start = part_boundary[part_idx];

		if(block_end < part_end)
			extraRows = part_end - block_end;
		else
			extraRows = 0;

		for( int i = 0; i < block_per_part; ++i){
			memset(col_hist, 0, max_col*sizeof( int));
			col_accum = 0;
			int local_start = part_boundary[part_idx]+i*ELL_threadSize;
			if(local_start > dimension) break;
			if (local_start + ELL_threadSize <= part_boundary[part_idx + 1]){
				local_end = local_start + ELL_threadSize;	
			}else{
				local_end = part_boundary[part_idx + 1];
			}
			block_boundary[block_idx + 1] = local_end;
			col_th = ceil( float(local_end - local_start)*0.75);
			for(int row_val = local_start; row_val < local_end; ++row_val){
				nonzero = (row_idx[row_val + 1] - row_idx[row_val]); 
				if(nonzero > max_col)
					printf("error with large nonzero\n");
				col_hist[nonzero - 1] += 1;
			}
			num_cols = 0;
			for( int i = 0; i < max_col; ++i){
				col_accum += col_hist[i];
				if(col_accum > col_th){
					num_cols = i + 1;
					break;
				}
			}
			if(num_cols == 0){
				printf("error at ELL_block convert rodr\n");
				exit(0);
			}		
			for(int row_val = local_start; row_val < local_end; ++row_val){
				ELL_val = 0;
				for(int j = row_idx[row_val]; j < row_idx[row_val + 1]; ++j){
					if(col[j] < part_start || col[j] >= part_start+ vector_cache_size){
						*size_COO += 1;
						ER_numInRow[ER_numOfRow] += 1;
						ER_line = true;
					} else 
						ELL_val += 1; 
				}
				if(ELL_val > num_cols){
					*size_COO += ELL_val - num_cols;
					ER_line = true;
				}
				if(ER_line){
					ER_rowVec[ER_numOfRow] = row_val;
					ER_numOfRow += 1;
					ER_line = false;
				}

			}
		}
		ELL_block_cols_vec[block_idx] = num_cols;	
		block_idx += 1;
	}
	if(extraRows > 0){
			for(int row_val = block_end; row_val < part_end; ++row_val){
				*size_COO += row_idx[row_val + 1] - row_idx[row_val];
			}	
		}
	}
	free(col_hist);
}

static void colsVecGenBlockELL( int* ELL_block_cols_vec, int* size_COO,
		const int dimension,const int max_col, const int *row_idx)
{
	 int avg_local = 0;
	// int block_nonz = 0;
	 int block_idx = 0;
	 int num_cols = 0;
	 int local_end;
	 int nonzero;
	 int col_th, col_accum; 
	 int *col_hist = ( int*)calloc(max_col, sizeof( int));

	*size_COO = 0;

	for( int row = 0; row < dimension; row += ELL_threadSize){
		memset(col_hist, 0, max_col*sizeof( int));
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
		for( int i = 0; i < max_col; ++i){
			col_accum += col_hist[i];
			if(col_accum > col_th){
				num_cols = i + 1;
				break;
			}
		}
		if(num_cols == 0){
			printf("error at ELL_block convert gen\n");
			exit(0);
		}		

		for(int row_val = row; row_val < local_end; ++row_val){
			 int nonzeros = row_idx[row_val + 1] - row_idx[row_val];
			if(nonzeros > num_cols){
				*size_COO += nonzeros - num_cols;	
			} 
		}
		ELL_block_cols_vec[block_idx] = num_cols;	
	}
	free(col_hist);
}


static void vecsGenER(int* I, 
	int* J,
	double* V,
	const int sizeCOO,
	const int numOfRowsER,
	const int* rowIdx,
	const int* numInRow,
	int* rodrList,
	int** rowVecER,
	int** biasVecER){
	
	numOfBlock = partSize*	
	
	*widthVecER = (int *) malloc(*sizeof(int));
	*biasVecER = (int *) malloc(*sizeof(int));
	*rowVecER = (int *) malloc(numOfRowER*sizeof(int));

	sortRordrListFull(numOfRows, rodr_list, numInRow);

	for(int i = 0; i< ER_numOfRows; ++i){
		origRowVal = ER_rowVec[i]; 
		//the Row location is not the value of row
		//it is the value of the ER_rowVec index, which is not equal to row
		//for example, row 100 is actually 10'th row of ER format submatrix
		rowLocation = rodr_list[i]; 
		ER_rowVec[rowLocation] = origRowVal;
		warpIdx = rowLocation/warpSize;		
		if(numInRow[i] > ER_colsVec[warpIdx]){
			ER_colsVec[warpBlockId] = numInRow[i];
		}
	}
	for(int i = 1; i < ER_numOfRows; ++i){
		ER_biasVec[i] = ER_biasVec[i-1] + warp_size*ER_colsVec[warpBlockId]; 
	}
}

static void COO2ER(int* I,
		int* J,
		double* V,
		const int sizeCOO,
		const int numInRowCOO,
		const int* numOfRowER,
		int* rowVecER,
		int* biasVecER,
		int* widthVecER,
		int* colER,
		double* valER){

	COO2ErVecsGen(I, J, V, sizeCOO, 
			numOfRowER, rowIdx, numInRow,
			rodrList, 
			&rowVecER, &biasVecER, &widthVecER);
	memset(ER_col, 0, sizeof(int)*size_ER);
	memset(ER_V, 0, sizeof(int)*size_ER);
	COOIdx = 0;
	for(int i = 0; i < numOfRow; ++i){
		rowLocation = rodr_list[i]; //again the Row location is not the value of row
		warpIdx = rowLocation/warpSize; 
		warpLane = rowLocation - warpIdx*warpSize;
		block_data_bias = ER_vec_bias[warpIdx]; 
		width = ER_vec_col_num[warpIdx];
		for(int j = 0; j < numInRow[i]; j++){
			dataIdx = block_data_bias + warpLane + warp_size*j; 
			ER_col[dataIdx] = J_COO[COOIdx]; 
			ER_V[dataIdx] = V_COO[COOIdx];
			COOIdx+=1;
		}	
	}
}

static void COO2EHYbCore( int* colELL, double* matrixELL,
		int* I_COO, int* J_COO, double* V_COO, 
		int* ER_rowVec, int* ER_numInRow, 
		const int ER_numOfRow,
		const int size_COO,
		const int* row_idx, const int* numInRow, 
		const int loc_num_of_row, 
		const int* ELL_block_bias_vec, const int* ELL_block_cols_vec, 
		const int* row_local, const int* col_local, 
		const double* matrix_local,
		int block_num, 
		const int* part_boundary, 
		volatile bool RODR, bool CACHED)
{ 
	int block_rowSize;
	int block_end = 0, part_start = 0, part_end = 0, extraRows = 0;
	int irregular=0;
	block_rowSize = ELL_threadSize;
	pointCOO = 0;
	int ELL_idx = 0;
	int boundary_start, boundary_end;
	int part_idx;
	int padding = 0;
	for ( int block_idx = 0; block_idx < block_num; ++block_idx){
		int num_cols = widthVecBlockELL[block_idx];
		int block_bias = biasVecBLockELL[block_idx];
		extraRows = 0;
		if(RODR){
			part_idx = block_idx/block_per_part;
			part_start = part_boundary[part_idx];
			part_end = part_boundary[part_idx + 1];
			boundary_start = part_boundary[part_idx] + (block_idx%block_per_part)*ELL_threadSize;
			block_end = part_boundary[part_idx] + block_per_part*ELL_threadSize;	

			if((block_idx + 1)%block_per_part == 0){
				if(part_end > block_end){
					extraRows = part_end - block_end;
					boundary_end = block_end;
				}
				else boundary_end = part_end;
			} else {
				boundary_end = boundary_start + ELL_threadSize;
			}
			block_rowSize = boundary_end - boundary_start;
		} else {
			boundary_start = block_idx*ELL_threadSize;
			boundary_end = (block_idx + 1)*ELL_threadSize;
			block_rowSize = ELL_threadSize;
		}

		block_rowSize = boundary_end - boundary_start;
		for( int i = 0; i < block_rowSize; ++i){
			int row_val = i + boundary_start;
			if(row_val >= loc_num_of_row) 
				break;
			ELL_idx = 0;
			for(int j = 0; j < numInRow[row_val]; ++j ){
				int col_val = col_local[row_idx[row_val] + j] ;
				if(ELL_idx < num_cols ){
					if(CACHED == true && (col_val < part_start || col_val >= part_start+vector_cache_size)){
						I_COO[pointCOO]=row_local[row_idx[row_val]+j];
						J_COO[pointCOO]=col_local[row_idx[row_val]+j];
						V_COO[pointCOO]=matrix_local[row_idx[row_val]+j];
						pointCOO++;
						if(pointCOO > size_COO){
							printf("error at pointCOO %d\n", pointCOO);
							exit(0);
						}
					} else {
						//goto ELL format
						colELL[block_bias + i + ELL_idx*block_rowSize] = col_val;
						matrixELL[block_bias + i + ELL_idx*block_rowSize]
							= matrix_local[row_idx[row_val] + j];
						ELL_idx += 1;
					}
				} else {
					//assign COO value
					I_COO[pointCOO]=row_local[row_idx[row_val]+j];
					J_COO[pointCOO]=col_local[row_idx[row_val]+j];
					V_COO[pointCOO]=matrix_local[row_idx[row_val]+j];
					pointCOO++;
					if(pointCOO > size_COO){
						printf("error at pointCOO %d at point 1\n", pointCOO);
						exit(0);
					}

				}

			}
			if(ELL_idx < num_cols ){
				for(int j = ELL_idx; j < num_cols; ++j){
					colELL[block_bias + i + j*block_rowSize] = part_start;
					matrixELL[block_bias + i + j*block_rowSize] = 0;
					if(part_start > 0 && (ELL_idx == 0 || colELL[block_bias] == 0))
						printf("debug poit at %d\n", row_val);
					padding += 1;
				}
			}
		}
		if(extraRows > 0){
			for( int row_val = block_end; row_val < part_end; ++row_val){
				for( int j = 0; j < numInRow[row_val]; ++j){
					I_COO[pointCOO]=row_local[row_idx[row_val]+j];
					J_COO[pointCOO]=col_local[row_idx[row_val]+j];
					V_COO[pointCOO]=matrix_local[row_idx[row_val]+j];
					pointCOO++;
					irregular+=1;
					if(pointCOO > size_COO){
						printf("error at pointCOO %d at point 2\n", pointCOO);

						exit(0);
					}
				}
			}	
		}
	}
	printf("irregular is %d, padding is %d pointCOO is %d, sizeCOO is %d\n", 
			irregular, padding, pointCOO, size_COO);
}

/*parameters, first line: output of HYB related parameters, 2nd line: output of HYB matrices, 3rd line:input of local COO matrix 
 4th line: ER index (didn't include values), 5th line: input variables*/

void COO2EHYB(matrixCOO* inputMatrix, 
		matrixEHYB* outputMatrix)
{ 

	int point_COO = 0;
	int ER_numOfRow;
	int* ER_numInRow = (int*)malloc(sizeof(int)*dimension);
	int* ER_rowVec = (int*)malloc(sizeof(int)*dimension);
	/*block_boundary will not be used outside this function*/
	int* block_boundary = ( int*)malloc((block_num + 1) * sizeof( int));
	int sizeCOO=0;	
	colsVecGenBlockELL(inputMatrix->widthVecBlockELL, 
			inputMatrix->biasVecBLockELL,
			&sizeCOO,
			);
	ELL_block_cols_vec_gen_rodr(ELL_block_cols_vec, 
		size_COO, 
		&ER_numOfRow,
		ER_rowVec,
		ER_widthVec,
		col_local,
		block_boundary,
		part_boundary, 
		loc_num_of_row, 
		max_col, 
		part_num,
		row_idx,
		CACHED);
	}

	if(*size_COO > 0){	
		*I_COO=( int *)malloc((*size_COO)*sizeof( int));
		*J_COO=( int *)malloc((*size_COO)*sizeof( int));
		*V_COO=(double *)malloc((*size_COO)*sizeof(double));
	}
	int ELL_matrixSize = 0;
	for( int i = 0; i < block_num; ++i){
		ELL_block_bias_vec[i] = ELL_matrixSize;
		if(RODR){
			ELL_matrixSize += (block_boundary[i+1] - block_boundary[i])*ELL_block_cols_vec[i];
		} else {
			ELL_matrixSize += ELL_threadSize*ELL_block_cols_vec[i];
		}
	}

	ELL_block_bias_vec[block_num] = ELL_matrixSize;
	*colELL=( int *)malloc(ELL_matrixSize*sizeof( int));
	*matrixELL=(double *)malloc(ELL_matrixSize*sizeof(double));

	COO2EHYB_block_core(*colELL, *matrixELL,
			*I_COO, *J_COO, *V_COO, *size_COO, COO_numInRow,
			row_idx, numInRow, loc_num_of_row, 
			ELL_block_bias_vec, ELL_block_cols_vec, 
			row_local, col_local, matrix_local, block_num, 
			part_boundary, RODR, CACHED);

	COO2ER(I_COO, J_COO, size_COO, COO_numInRow,
			ER_numOfRow, ER_rowVec, 
			ER_biasVec, ER_colsVec,
			ER_col, ER_val);

	free(*I_COO);
	free(*J_COO);
	free(*V_COO);
	free(COO_numInRow);
	free(block_boundary);
}

