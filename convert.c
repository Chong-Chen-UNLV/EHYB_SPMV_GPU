#include "kernel.h"
#include "solver.h"
#include "convert.h"
#include "Partition.h"


static void sortRordrListFull(unsigned int dimension,
		int* reorderList, 
		int* numInRow)
{
	//This function generated the reorderList according to the
	//numInRow, since it is a wrapper of qsort function, and 
	//the logic is a bit different from similar function at
	//reordering.c, we make a seperate function in this file
	//instead of reuse the function in reordering.c
	rowS* rodrSVec = (rowS*) malloc(dimension*sizeof(rowS));
	for(unsigned int i = 0; i < dimension; ++i){
		//the location of elements at SVec will be changed 
		//record the original idx since it will be used for
		//new rodr_list
		rodrSVec[reorderList[i]].idx = i; 
		rodrSVec[reorderList[i]].nonzeros = numInRow[i];	
	}	
	for(unsigned int i = 0; i < nparts; ++i){
		qsort(rodrSVec, dimension, sizeof(rowS), rowSCompare);	
	}	
	for(unsigned int i = 0; i < dimension; ++i){
		reorderList[rodrSVec[i].idx] = i;
	}
	free(rodrSVec);
}

static void vecsGenBlockELL(matrixCOO* inputMatrix,
		matrixEHYB* outputMatrix,
		int* reorderListER,
		int* numInRowER)
{
	/*Reoder functions did an extra iteration 
	 * which provides numInRowER without ER reordering
	 * */
	// int block_nonz = 0;
	int blockIdx = 0;
	int numCols = 0;
	int nonZeros;
	int colTh, colAccum; 
	int *colHist = ( int*)calloc(max_col, sizeof( int));
	int *numInRow = inputMatrix->numInRow;
	int *numInRow2 = inputMatrix->numInRow2; 

	int extraRows = 0;
	int numOfRowER = 0;

	for(int partIdx = 0; partIdx < inputMatrix->nParts; ++partIdx){
		int checkER = 0;	
		int partStart = partBoundary[partIdx];
		int partEnd = partBoundary[partIdx + 1];
		if(partEnd > partStart + blockPerPart*warpSize)
			extraRows = partEnd - partStart - blockPerPart*warpSize;
		else
			extraRows = 0;

		for(int iter = 0; iter < blockPerPart; ++iter){
			blockIdx = iter + blockPerPart*partIdx;
			int blockStart = partStart + iter*warpSize;
			memset(colHist, 0, maxCol*sizeof(int));
			actualRow = 0;
			for(int row = blockStart; row < blockStart + warpSize; ++row){
				if(row < partEnd){
					if(numInRow2[row] > 0){
						colHist[numInRow2[row]- 1] += 1;
						actualRow += 1;	
					}
					if(numInRow2[row] == numInRow[row])
						numOfRowER += 1;
				} 
			}
			colTh = ceil(float(actualRow)*0.75);
			width = 0;
			for( int i = 0; i < max_col; ++i){
				colAccum += colHist[i];
				if(colAccum > colTh){
					numCols = i + 1;
					if(checkER == 0) checkER = 1;
					break;
				}
			}
			if(numCols == 0){
				printf("error at ELL_block convert rodr\n");
				exit(0);
			}
			widthVecBlockELL[blockIdx] = numCols;
		}
		if(extraRows > 0){
			if(checkER == 0) checkER = 1;
			numOfRowER += 1;
			for(int row = partStart + blockPerPart*warpSize; row < partEnd; ++row){
				numInRowER[row] += numInRow[row];
			}	
		}
	}
	if(checkER == 0){
		printf("perfect matrix, not in the scope of this study\n");		
		exit(0);
	}
	outputMatrix->numOfRowER = numOfRowER;
	outputMatrix->rowVecER = (int*)malloc(numOfRowER*sizeof(int));
	outputMatrix->biasVecER = (int*)malloc(ceil((float numOfRowER)/warpSize)*sizeof(int));
	free(col_hist);
}

static void vecsGenER(matrixEHYB* inputMatrix, int* numInRowER, int* reorderListER)
{
	int warpIdx;

	sortRordrListFull(dimension, reorderListER, numInRowER);

	for(int i = 0; i < dimension; ++i){
		//the Row location is not the value of row
		//it is the value of the ER_rowVec index, which is not equal to row
		//for example, row 100 is actually 10'th row of ER format submatrix
		if(numInRowER[i] > 0){
			rowLocation = reorderList[i]; 
			rowVecER[rowLocation] = i;
			warpIdx = rowLocation/warpSize;		
			if(numInRowER[i] > widthVecER[warpIdx]){
				widthVecER[warpIdx] = numInRowER[i];
			}
		}
	}

	for(int i = 1; i < inputMatrix->numOfRowER; ++i){
		biasVecER[i] = biasVecER[i-1] + warpSize*widthVecER[warpBlockId]; 
	}
}

static void COO2EHYBCore(matrixCOO* inputMatrix, 
		matrixEHYB* outputMatrix,
		int* reorderListER)
{ 
	int* partBoundary = inputMatrix->partBoundary;
	int* widthVecBlockELL = inputMatrix->widthVecBlockELL;
	int* biasVecBLockELL = inputMatrix->biasVecBLockELL;
	int blockEnd = 0, partStart = 0, partEnd = 0, extraRows = 0;
	int irregular=0;
	int writedInRowELL = 0;
	int writedInRowER = 0;
	int partIdx;
	int padding = 0;
	int nBlocks = inputMatrix->nParts*blockPerPart;

	for(int blockIdx = 0; blockIdx < nBlocks; ++blockIdx){
		int widthBlockELL = widthVecBlockELL[blockIdx];
		int biasBlockELL = biasVecBLockELL[blockIdx]; 
		partIdx = blockIdx/blockPerPart;
		partStart = partBoundary[partIdx];
		partEnd = partBoundary[partIdx + 1];
		int blockStart = partStart + warpSize*(blockIdx%blockPerPart);
		if(blockIdx%blockPerPart == 0)
			extraRows = partEnd - (partStart + warpSize*blockPerPart);

		for(int i = 0; i < warpSize; ++i){
			writedInRowER = 0;
			writedInRowELL= 0;
			rowVal = blockStart + i;
			if(reorderListER[rowVal] > 0){
				rowLocER = reorderListER[rowVal];
				blockIdxER = rowLocER/warpSize;
				rowLocInBlockER = rowLocER - blockIdxER*warpSize;
				biasER = biasVecER[blockIdxER];
				widthER = widthVecER[blockIdxER];
			}

			if(rowVal < partEnd){
				for(int j = 0; j < numInRow[rowVal]; ++j){
					int tmpIdx = j + rowIdx[rowVal];	
					if(J[tmpIdx] < partEnd && J[tmpIdx] > partStart && writedInRowELL < widthBlockELL){
						colBlockELL[biasBlockELL+i+writedInRowELL*warpSize] = J[tmpIdx];
						valBlockELL[biasBlockELL+i+writedInRowELL*warpSize] = V[tmpIdx];
						writedInRowELL += 1;	
					} else {
						colER[biasER+rowLocInBlockER+writedInRowER*warpSize] = J[tmpIdx];
						valER[biasER+rowLocInBlockER+writedInRowER*warpSize] = V[tmpIdx];
						writedInRowER += 1;
					}
				}
			} 
		}	

		if(extraRows > 0){
			for(int rowVal = blockEnd; rowVal < partEnd; ++rowVal){
				if(rowVecER[reorderListER[rowVal]] != rowVal){
					printf("error at rowVecER\n");
					exit(0);
				}
				rowLocER = reorderListER[rowVal];
				blockIdxER = rowLocER/warpSize;
				rowLocInBlockER = rowLocER - blockIdxER*warpSize;
				biasER = biasVecER[blockIdxER];
				widthER = widthVecER[blockIdxER];
				for(int j = 0; j < numInRow[rowVal]; ++j){
					colER[biasER+rowLocInBlockER+j*warpSize]=J[rowIdx[rowVal]+j];
					valER[biasER+rowLocInBlockER+j*warpSize]=V[rowIdx[row_val]+j];
					irregular+=1;
					if(pointER > sizeER){
						printf("error at pointER %d\n", pointER);
						exit(0);
					}
				}
			}
		}
		extraRow = 0;
	}

}

/*parameters, first line: output of HYB related parameters, 2nd line: output of HYB matrices, 3rd line:input of local COO matrix 
  4th line: ER index (didn't include values), 5th line: input variables*/

void COO2EHYB(matrixCOO* inputMatrix, 
		matrixEHYB* outputMatrix,
		int* sizeBlockELL, 
		int* sizeER)
{ 

	int* reorderListER = (int*)calloc(sizeof(int)*dimension);
	int* numInRowER = (int*)calloc(sizeof(int)*dimension);
	outputMatrix->widthVecBlockELL = (int*)malloc(outputMatrix->nParts*blockPerPart*sizeof(int));
	outputMatrix->biasVecBlockELL = (int*)malloc(outputMatrix->nParts*blockPerPart*sizeof(int));
	/*block_boundary will not be used outside this function*/
	vecsGenBlockELL(inputMatrix, outputMatrix, numInRowER);

	*sizeBlockELL = 0;
	for(int i = 0; i < blockPerPart*outputMatrix->nParts; ++i){
		*sizeBlockELL += warpSize*outputMatrix->widthVecBlockELL[i]; 
	}
	outputMatrix->valBlockELL = (double*) malloc(sizeof(double)*(*sizeBlockELL));
	outputMatrix->colBlockELL = (int*) malloc(sizeof(int)*(*sizeBlockELL));

	vecsGenER(inputMatrix, outputMatrix, numInRowER, reorderListER);
	*sizeER = 0;	
	for(int i = 0; i < ceil((float outputMatrix->numOfRowER)/warpSize); ++i){
		*sizeER += warpSize*outputMatrix->widthVecER[i];
	}

	outputMatrix->valER = (double*) malloc(sizeof(double)*(*sizeER));
	outputMatrix->colER = (int*) malloc(sizeof(int)*(*sizeER));

	COO2EHYBCore(inputMatrix, 
			outputMatrix,
			reorderListER);

}

