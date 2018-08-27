#include "kernel.h"
#include "solver.h"

void COO2ELL(const unsigned int *rowLocal, const unsigned int *colLocal, const float* matrixLocal, unsigned int **colELL,
	float **matrixELL, unsigned int **I_COO, unsigned int **J_COO, float **V_COO,const unsigned int *numInRow, 
	const unsigned int *rowNumAccum, const unsigned int localMatrixSize, const unsigned int localNumofRow, 
	unsigned int *sizeOut, unsigned max_in, unsigned int *max_out){

	unsigned int maxRowNum;

	unsigned int sizeCOO=0, pointCOO=0;
	unsigned int rowBias=rowLocal[0];
	unsigned int numBias=rowNumAccum[rowBias];
	
	maxRowNum=ceil(((float) localMatrixSize)*1.2/((float) localNumofRow));
	if (maxRowNum > max_in){
		maxRowNum = max_in;	
	}
	else{
		for (unsigned int i=0;i<localNumofRow;i++)
		{
			unsigned int rowIndex=i+rowLocal[0];
			if (numInRow[i+rowBias] > maxRowNum)
				sizeCOO+=numInRow[i+rowBias]- maxRowNum;
		}	
	}
	if(sizeCOO > 0){	
		*I_COO=(unsigned int *)malloc(sizeCOO*sizeof(unsigned int));
		*J_COO=(unsigned int *)malloc(sizeCOO*sizeof(unsigned int));
		*V_COO=(float *)malloc(sizeCOO*sizeof(float));
	}
	*colELL=(unsigned int *)malloc(localNumofRow*maxRowNum*sizeof(unsigned int));
	*matrixELL=(float *)malloc(localNumofRow*maxRowNum*sizeof(float));
	
	unsigned int irregular=0;
	for (unsigned int i=0;i<localNumofRow;i++)
	{
		unsigned int rowIndex=i+rowLocal[0];
		//goto COO format
		if (numInRow[i+rowBias]>maxRowNum) {
			for (unsigned int j=0;j<numInRow[i+rowBias];j++){

				//the ELL value should still be set as zero
				if (j<maxRowNum){
					(*colELL)[i+j*localNumofRow]=colLocal[rowNumAccum[i]+j-numBias];
					(*matrixELL)[i+j*localNumofRow]=matrixLocal[rowNumAccum[rowIndex]+j-numBias];
				}
				else{
					//assign COO value
					(*I_COO)[pointCOO]=rowLocal[rowNumAccum[rowIndex]+j-numBias];
					(*J_COO)[pointCOO]=colLocal[rowNumAccum[rowIndex]+j-numBias];
					(*V_COO)[pointCOO]=matrixLocal[rowNumAccum[rowIndex]+j-numBias];
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
				if (j<numInRow[i+rowBias]){
					(*colELL)[i+j*localNumofRow]=colLocal[rowNumAccum[rowIndex]+j-numBias];
					(*matrixELL)[i+j*localNumofRow]=matrixLocal[rowNumAccum[rowIndex]+j-numBias];
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

static const unsigned int num_cols find_num_cols_by_block(const int* I, const int* J, unsigned int* ELL_deepth_vec,
				const int* dimension, const unsigned int *numInRow)
	
	unsigned int maxRow = 0;
	unsigned int avgLocal = 0;
	unsigned int blockNonZ = 0;

	for(int block_row =0; block_row < dimension; block_row += ELL_blockSize){
		if (bolck_row + ELL_blockSize <= dimension){
			blockNonZ = rowNumAccum[block_row + ELL_blockSize] - rowNumAccum[block_row];
			avgPerLn = ceil( ((float) blockNonZ)/ELL_blockSize);
		}
		else{ 
			blockNonZ = rowNumAccum[dimension] - rowNumAccum[block_row];
			avgPerLn = ceil( ((float) blockNonZ)/(dimension - block_row));
		}
		num_cols = avgPerLin*1.2;
		for(int row_idx = block_row; row_idx < block_row + ELL_blockSize; ++row_idx){
			if
		}
		for(int row_idx = block_row; row_idx < block_row + ELL_blockSize; ++row_idx){
			sizeCOO += ;
	
		}	
	}

	num_cols = ceil(((float) localMatrixSize)*1.2/((float) localNumofRow));
	if(num_cols < maxRow)
		num_cols = maxRow;

	return num_cols;
}

static void COO2ELL_block_core(const int* row_local, const int* col_local, const float* matrixLocal, 
			const int row_idx, const int* colELL, const float* matrixELL,
			const int* I_COO, const int* J_COO, const float* V_COO,
			){

	unsigned int irregular=0;
	
	for (unsigned int row = 0; row < localNumofRow; row += ELL_blockSize){
			
		unsigned int row_idx = row + row_local[0];
		unsigned num_cols = num_cols_per_row_vec[block_idx];
		
		//goto COO format
		if (numInRow[i+rowBias] > num_cols) {
			for (unsigned int j=0;j<numInRow[i+rowBias];j++){

				//the ELL value should still be set as zero
				if (j < maxRowNum){
					(*colELL)[i+j*localNumofRow]=colLocal[rowNumAccum[i]+j-numBias];
					(*matrixELL)[i+j*localNumofRow]=matrixLocal[rowNumAccum[rowIndex]+j-numBias];
				}
				else{
					//assign COO value
					(*I_COO)[pointCOO]=rowLocal[rowNumAccum[rowIndex]+j-numBias];
					(*J_COO)[pointCOO]=colLocal[rowNumAccum[rowIndex]+j-numBias];
					(*V_COO)[pointCOO]=matrixLocal[rowNumAccum[rowIndex]+j-numBias];
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
				if (j<numInRow[i+rowBias]){
					(*colELL)[i+j*localNumofRow]=colLocal[rowNumAccum[rowIndex]+j-numBias];
					(*matrixELL)[i+j*localNumofRow]=matrixLocal[rowNumAccum[rowIndex]+j-numBias];
				}
				//write zero
				else{
					(*colELL)[i+j*localNumofRow]=0;
					(*matrixELL)[i+j*localNumofRow]=0;
				}
			}
		}
	}
}


void COO2ELL_block(const unsigned int *rowLocal, const unsigned int *colLocal, const float* matrixLocal, unsigned int **colELL,
		float **matrixELL, unsigned int **I_COO, unsigned int **J_COO, float **V_COO,const unsigned int *numInRow, 
		const unsigned int *rowNumAccum, const unsigned int localMatrixSize, const unsigned int localNumofRow, 
		unsigned int *sizeOut, unsigned max_in, unsigned int *max_out){

	unsigned int maxRowNum;

	unsigned int sizeCOO=0, pointCOO=0;
	unsigned int rowBias=rowLocal[0];
	unsigned int numBias=rowNumAccum[rowBias];
	
	for (unsigned int row = 0; row < localNumofRow; row += ELL_blockSize){
		
	}
	maxRowNum = find_num_cols_by_block(); 
	if(sizeCOO > 0){	
		*I_COO=(unsigned int *)malloc(sizeCOO*sizeof(unsigned int));
		*J_COO=(unsigned int *)malloc(sizeCOO*sizeof(unsigned int));
		*V_COO=(float *)malloc(sizeCOO*sizeof(float));
	}
	*colELL=(unsigned int *)malloc(localNumofRow*maxRowNum*sizeof(unsigned int));
	*matrixELL=(float *)malloc(localNumofRow*maxRowNum*sizeof(float));

	row_idx += ELL_block_size;
	COO2ELL_block_core(rowLocal, colLocal, matrixLocal, row_idx, *colELL, 
			*matrixELL, *I_COO, *J_COO, *V_COO, rowBias, numBias);	
	
	if (maxRowNum > max_in){
		maxRowNum = max_in;	
	}
	else{
		for (unsigned int i=0;i<localNumofRow;i++)
		{
			unsigned int rowIndex=i+rowLocal[0];
			if (numInRow[i+rowBias] > maxRowNum)
				sizeCOO+=numInRow[i+rowBias]- maxRowNum;
		}	
	}
	
	
	
	
	*max_out = maxRowNum;
	*sizeOut = sizeCOO;
	//printf("irregular row Number is %d\n",irregular);
	
}
