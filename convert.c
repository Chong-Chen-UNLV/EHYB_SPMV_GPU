#include "kernel.h"
#include "solver.h"

void COO2ELL(const unsigned int *rowLocal, const unsigned int *colLocal, const double* matrixLocal, unsigned int **colELL,
	double **matrixELL, unsigned int **I_COO, unsigned int **J_COO, double **V_COO,const unsigned int *numInRow, 
	const unsigned int *rowNumAccum, const unsigned int localMatrixSize, const unsigned int localNumofRow, 
	unsigned int *sizeOut, unsigned max_in, unsigned int *max_out){

	unsigned int maxRowNum;

	unsigned int sizeCOO=0, pointCOO=0;
	unsigned int row_bias=rowLocal[0];
	unsigned int num_bias=rowNumAccum[row_bias];
	
	maxRowNum=ceil(((double) localMatrixSize)*1.2/((double) localNumofRow));
	if (maxRowNum > max_in){
		maxRowNum = max_in;	
	}
	else{
		for (unsigned int i=0;i<localNumofRow;i++){
			unsigned int row_idx=i+rowLocal[0];
			if (numInRow[i+row_bias] > maxRowNum)
				sizeCOO+=numInRow[i+row_bias]- maxRowNum;
		}	
	}
	if(sizeCOO > 0){	
		*I_COO=(unsigned int *)malloc(sizeCOO*sizeof(unsigned int));
		*J_COO=(unsigned int *)malloc(sizeCOO*sizeof(unsigned int));
		*V_COO=(double *)malloc(sizeCOO*sizeof(double));
	}
	*colELL=(unsigned int *)malloc(localNumofRow*maxRowNum*sizeof(unsigned int));
	*matrixELL=(double *)malloc(localNumofRow*maxRowNum*sizeof(double));
	
	unsigned int irregular=0;
	for (unsigned int i=0;i<localNumofRow;i++){
		unsigned int row_idx=i+rowLocal[0];
		//goto COO format
		if (numInRow[i+row_bias]>maxRowNum) {
			for (unsigned int j=0;j<numInRow[i+row_bias];j++){

				//the ELL value should still be set as zero
				if (j<maxRowNum){
					(*colELL)[i+j*localNumofRow]=colLocal[rowNumAccum[i]+j-num_bias];
					(*matrixELL)[i+j*localNumofRow]=matrixLocal[rowNumAccum[row_idx]+j-num_bias];
				}
				else{
					//assign COO value
					(*I_COO)[pointCOO]=rowLocal[rowNumAccum[row_idx]+j-num_bias];
					(*J_COO)[pointCOO]=colLocal[rowNumAccum[row_idx]+j-num_bias];
					(*V_COO)[pointCOO]=matrixLocal[rowNumAccum[row_idx]+j-num_bias];
					pointCOO++;
					if(pointCOO > sizeCOO)
						printf("error at pointCOO %d\n", pointCOO);
				}
			}
			irregular=irregular+1;
		}
		//goto ELL format
		else {
			for (unsigned int j=0;j<maxRowNum;j++){
				//write the ELL data
				if (j<numInRow[i+row_bias]){
					(*colELL)[i+j*localNumofRow]=colLocal[rowNumAccum[row_idx]+j-num_bias];
					(*matrixELL)[i+j*localNumofRow]=matrixLocal[rowNumAccum[row_idx]+j-num_bias];
				}
				//write zero
				else{
					(*colELL)[i+j*localNumofRow]=0;
					(*matrixELL)[i+j*localNumofRow]=0;
				}
			}
		}
	}
	*max_out = maxRowNum;
	*sizeOut = sizeCOO;
	//printf("irregular row Number is %d\n",irregular);
	
}

static void num_cols_vec_gen(unsigned int* num_cols_vec, unsigned int* sizeCOO
				const int* I, const int* J,
				const int* dimension, const unsigned int *rowNumAccum){
	
	unsigned int maxCol = 0;
	unsigned int avgLocal = 0;
	unsigned int blockNonZ = 0;
	unsigned int block_idx = 0;
	*sizeCOO = 0;

	for(int row = 0; row < dimension; row += ELL_blockSize){
		block_idx = row/ELL_blockSize;
		if (bolck_row + ELL_blockSize <= dimension){
			blockNonZ = rowNumAccum[row + ELL_blockSize] - rowNumAccum[row];
			if(blockNonZ > maxCol) maxCol = blockNonZ;
			avgPerLn = ceil( ((double) blockNonZ)/ELL_blockSize);
		}else{
			blockNonZ = rowNumAccum[dimension] - rowNumAccum[row];
			avgPerLn = ceil( ((double) blockNonZ)/(dimension - row));
		}
		if(maxCol > avgPerLin*1.2)
			num_cols = avgPerLin*1.2;
		else
			num_cols = maxCol;
		
		for(int row_idx = row; row_idx < row + ELL_blockSize; ++row_idx){
			if(numInRow[row_idx] > num_cols){
				*sizeCOO += numInRow[row_idx] - num_cols;	
			} 
		}
		num_cols_vec[block_idx] = num_cols;	
	}
}

static void COO2ELL_block_core(int* colELL, double* matrixELL,
		int* I_COO, int* J_COO, double* V_COO, 
		unsigned int* num_cols_vec, unsigned int* block_ELL_bias, 
		const int* row_local, const int* col_local, const double* matrixLocal){ 

	unsigned int irregular=0;
	unsigned int ELL_blocks = ceil( ((float)row_local)/ELL_threadSize);
	unsigned int *ELL_bias_vec = (unsigned int*) malloc(sizeof(unsigned int)*ELL_blocks);
	ELL_bias_vec[0] = 0;
	for (unsigned int blk = 1; blk < ELL_block_size; ++blk ){
		ELL_bias_vec[blk] = num_cols_vec[blk]*ELL_threadSize + ELL_bias_vec[blk - 1]
	}
	for (unsigned int row = 0; row < localNumofRow; row += ELL_blockSize){
		unsigned int block_idx = row/ELL_blockSize;
		unsigned int num_cols = num_cols_vec[block_idx];
		unsigned int ELL_bias = block_ELL_bias_vec[block_idx]; 
		
		for(i = row; i < row + ELL_blockSize; ++i){
			row_idx = i + row_bias;
			if (numInRow[i+row_bias] > num_cols) {//goto COO format
				for (unsigned int j = 0; j < numInRow[row_idx+row_bias]; ++j){
					//the ELL value should still be set 
					if (j < num_cols){
						(*colELL)[ELL_bias + j*ELL_treadSize] = colLocal[rowNumAccum[row_idx]+j-num_bias];
						(*matrixELL)[ELL_bias + j*ELL_treadSize]
							= matrixLocal[rowNumAccum[row_idx]+j-num_bias];
					}
					else{
						//assign COO value
						(*I_COO)[pointCOO]=rowLocal[rowNumAccum[row_idx]+j-num_bias];
						(*J_COO)[pointCOO]=colLocal[rowNumAccum[row_idx]+j-num_bias];
						(*V_COO)[pointCOO]=matrixLocal[rowNumAccum[row_idx]+j-num_bias];
						pointCOO++;
						if(pointCOO > sizeCOO)
							printf("error at pointCOO %d\n", pointCOO);
					}
				}
				irregular=irregular+1;
			}else {//goto ELL format
				for (unsigned int j=0; j < num_cols; ++j){
					//write the ELL data
					if (j<numInRow[i+row_bias]){
						(*colELL)[ELL_bias + j*ELL_treadSize] = colLocal[rowNumAccum[row_idx]+j-num_bias];
						(*matrixELL)[ELL_bias + j*ELL_treadSize]
								= matrixLocal[rowNumAccum[row_idx]+j-num_bias];
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
}


void COO2ELL_block(unsigned int* num_cols_vec, unsigned int *sizeCOO, unsigned int* block_ELL_bias_vec, 
		double **matrixELL, unsigned int **I_COO, unsigned int **J_COO, double **V_COO,
		const unsigned int *rowLocal, const unsigned int *colLocal, const double* matrixLocal, unsigned int **colELL,
		const unsigned int *numInRow, 
		const unsigned int *rowNumAccum, const unsigned int localMatrixSize, const unsigned int localNumofRow, 
		unsigned int *sizeOut, unsigned max_in, unsigned int *max_out){


	unsigned int sizeCOO=0, pointCOO=0;
	unsigned int row_bias=rowLocal[0];
	unsigned int num_bias=rowNumAccum[row_bias];
	

	num_cols_vec_gen(num_cols_vec, &sizeCOO,
			rowLocal, colLocal,
			localNumofRow, rowNumAccum); 

	if(sizeCOO > 0){	
		*I_COO=(unsigned int *)malloc(sizeCOO*sizeof(unsigned int));
		*J_COO=(unsigned int *)malloc(sizeCOO*sizeof(unsigned int));
		*V_COO=(double *)malloc(sizeCOO*sizeof(double));
	}
	*colELL=(unsigned int *)malloc(localNumofRow*maxRowNum*sizeof(unsigned int));
	*matrixELL=(double *)malloc(localNumofRow*maxRowNum*sizeof(double));

	row_idx += ELL_block_size;
	COO2ELL_block_core(rowLocal, colLocal, matrixLocal, row_idx, *colELL, 
			*matrixELL, *I_COO, *J_COO, *V_COO, row_bias, num_bias);	
	
	if (maxRowNum > max_in){
		maxRowNum = max_in;	
	}else{
		for (unsigned int i=0;i<localNumofRow;i++){
			unsigned int row_idx=i+rowLocal[0];
			if (numInRow[i+row_bias] > maxRowNum)
				sizeCOO+=numInRow[i+row_bias]- maxRowNum;
		}	
	}
	
	
	
	
	*max_out = maxRowNum;
	*sizeOut = sizeCOO;
	//printf("irregular row Number is %d\n",irregular);
	
}
