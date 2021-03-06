#include "kernel.h"
#include "spmv.h"
#include "convert.h"
#include "Partition.h"

#define warpSize 32

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
		rodrSVec[i].idx = i; 
		rodrSVec[i].nonzeros = numInRow[i];	
	}	
	qsort(rodrSVec, dimension, sizeof(rowS), rowSCompare);	
	for(unsigned int i = 0; i < dimension; ++i){
		reorderList[rodrSVec[i].idx] = i;
	}
	free(rodrSVec);
}

static void longVecCopy(
		const int* longVecBoundary,
		const int* longVecRow,
		matrixCOO* inputMatrix,
		matrixEHYB* outputMatrix){
		
	int nLongVec = outputMatrix->nLongVec;
	if(nLongVec > 0){
		int vecEleSize = longVecBoundary[nLongVec]; 
		outputMatrix->longVecBoundary = (int*)malloc((nLongVec+1)*sizeof(int));
		outputMatrix->longVecRow = (int*)malloc(nLongVec*sizeof(int));
		outputMatrix->longVecCol = (int*)malloc(vecEleSize*sizeof(int));
		outputMatrix->longVecVal = (double*)malloc(vecEleSize*sizeof(double));
		memcpy(outputMatrix->longVecBoundary, longVecBoundary, sizeof(int)*(nLongVec+1) );
		memcpy(outputMatrix->longVecRow, longVecRow, sizeof(int)*nLongVec );
		for(int i = 0; i < nLongVec; ++i){
			int row = longVecRow[i];				
			for(int j = 0; j < inputMatrix->numInRow[row]; ++j){
				int idxL = longVecBoundary[i] + j;
				int idxM = inputMatrix->rowIdx[row] + j;
				outputMatrix->longVecCol[idxL] = inputMatrix->I[idxM];
				outputMatrix->longVecVal[idxL] = inputMatrix->V[idxM];
			}
		}
	}
	
}

static void vecsGenBlockELL(matrixCOO* inputMatrix,
		matrixEHYB* outputMatrix,
		int* tmpBoundary,
		int* reorderListER,
		int* numInRowER)
{
	/*Reoder functions did an extra iteration 
	 * which provides numInRowER without ER reordering
	 * */
	int blockIdx = 0;
	int16_t numCols = 0;
	int maxCol = inputMatrix->maxCol;
	int nonZeros;
	int *numInRow = inputMatrix->numInRow;
	int *numInRow2 = inputMatrix->numInRow2; 

	int extraRows = 0;
	int numOfRowER = 0;
	int toER = 0;
	int blockPerPart = (inputMatrix->vectorCacheSize/warpSize);
	int longVecId = 1;
	const int buffSize1 = inputMatrix->dimension/20;//assume that number of long row is less than 5%, it is a reasonable assumption
	int* longVecBoundary = (int*)malloc(buffSize1*sizeof(int));
	longVecBoundary[0] = 0;
	int* longVecMap = (int*)malloc(buffSize1*sizeof(int));//map from real row to the longVec part id
	int* longVecRow = (int*)malloc(buffSize1*sizeof(int));//reverse of the "map", longVec part id back to real row
	for(int partIdx = 0; partIdx < inputMatrix->nParts; ++partIdx){
		int partStart = inputMatrix->partBoundary[partIdx];
		int realPartStart = partStart;
		int partEnd = inputMatrix->partBoundary[partIdx + 1];
		//the long vec exist at the begining of partition
		while(numInRow2[realPartStart] > threadLongVec){
			//inputMatrix->pureVecId[pureVecNum] = row;
			//inputMatrix->pureVecId[pureVecNum] += 1;
			longVecBoundary[longVecId] = numInRow[realPartStart] + longVecBoundary[longVecId - 1]; 
			longVecMap[realPartStart] = longVecId;	
			longVecRow[longVecId] = realPartStart;	
			longVecId += 1;
			realPartStart += 1;
			outputMatrix->nLongVec+=1;
		}
		tmpBoundary[partIdx] = realPartStart;
		if(partEnd > partStart + blockPerPart*warpSize)
			extraRows = partEnd - partStart - blockPerPart*warpSize;
		//if(partIdx == 2)
		//	printf("start fined error on 0 numCols\n");	
		for(int iter = 0; iter < blockPerPart; ++iter){
			blockIdx = iter + blockPerPart*partIdx;
			int blockStart = partStart + iter*warpSize;
			numCols = 0;
			for(int row = blockStart; (row < blockStart + warpSize && row < partEnd); ++row){
				if(row >=realPartStart){
					if(numInRow2[row] > numCols)
						numCols = numInRow2[row];
					if(numInRow2[row] != numInRow[row]){
						numOfRowER += 1;
						numInRowER[row] = numInRow[row] - numInRow2[row];
						toER += numInRowER[row];
					}
				} 
			}
			if(numCols == 0 && blockStart < partEnd && blockStart + warpSize > realPartStart){
				printf("error at ELL_block convert rodr\n");
				exit(0);
			}
			outputMatrix->widthVecBlockELL[blockIdx] = numCols;
		}
		if(extraRows > 0){
			for(int row = partStart + blockPerPart*warpSize; row < partEnd; ++row){
				numOfRowER += 1;
				numInRowER[row] += numInRow[row];
				toER += numInRowER[row];
			}	
		}
	}
	if(numOfRowER == 0){
		printf("perfect matrix, not in the scope of this study\n");		
		exit(0);
	}
	printf("toER is %d, kernel calculation is %d\n", toER, inputMatrix->totalNum - toER);
	outputMatrix->numOfRowER = numOfRowER;
	outputMatrix->rowVecER = (int*)malloc(numOfRowER*sizeof(int));
	outputMatrix->biasVecER = (int*)malloc(ceil(((double) numOfRowER)/warpSize)*sizeof(int));
	//----start longVec transfer-----------
	longVecCopy(longVecBoundary, longVecRow, inputMatrix, outputMatrix);
}

static void vecsGenER(matrixEHYB* inputMatrix, int* numInRowER, int* reorderListER)
{
	int warpIdx, rowLocation;

	sortRordrListFull(inputMatrix->dimension, reorderListER, numInRowER);

	for(int i = 0; i < inputMatrix->dimension; ++i){
		//the Row location is not the value of row
		//it is the value of the ER_rowVec index, which is not equal to row
		//for example, row 100 is actually 10'th row of ER format submatrix
		if(numInRowER[i] > 0){
			rowLocation = reorderListER[i]; 
			inputMatrix->rowVecER[rowLocation] = i;
			warpIdx = rowLocation/warpSize;		
			if(numInRowER[i] > inputMatrix->widthVecER[warpIdx]){
				inputMatrix->widthVecER[warpIdx] = numInRowER[i];
			}
		}
	}

}

static void COO2EHYBCore(matrixCOO* inputMatrix, 
		matrixEHYB* outputMatrix,
		int* tmpBoundary,
		int* reorderListER,
		int* numInRowER)
{ 
	int* partBoundary = inputMatrix->partBoundary;
	
	int* rowIdx = inputMatrix->rowIdx;
	int* numInRow = inputMatrix->numInRow;
	int* I = inputMatrix->I;
	int* J = inputMatrix->J;
	double* V = inputMatrix->V;

	int16_t* widthVecBlockELL = outputMatrix->widthVecBlockELL;
	int* biasVecBlockELL = outputMatrix->biasVecBlockELL;
	int16_t* colBlockELL = outputMatrix->colBlockELL; 
	double* valBlockELL = outputMatrix->valBlockELL; 

	int16_t* widthVecER = outputMatrix->widthVecER;
	int* rowVecER = outputMatrix->rowVecER; 
	int* biasVecER = outputMatrix->biasVecER;
	int* colER = outputMatrix->colER; 
	double* valER = outputMatrix->valER; 

	int partIdx, partStart, realPartStart, partEnd, fetchEnd, extraRows;
	int rowLocER, blockIdxER, rowLocInBlockER, biasER; 
	int16_t widthER;     	
	int biasBlockELL, widthBlockELL;     	
	int irregular=0;
	int writedInRowELL = 0;
	int writedInRowER = 0;
	int wasteElement = 0;
	int blockPerPart = (inputMatrix->vectorCacheSize/warpSize);

	int nBlocks = inputMatrix->nParts*blockPerPart;
	
	for(int blockIdx = 0; blockIdx < nBlocks; ++blockIdx){
		widthBlockELL = widthVecBlockELL[blockIdx];
		biasBlockELL = biasVecBlockELL[blockIdx]; 
		partIdx = blockIdx/blockPerPart;
		partStart = partBoundary[partIdx];
		realPartStart = tmpBoundary[partIdx];	
		partEnd = partBoundary[partIdx + 1];
		fetchEnd =  partStart + inputMatrix->vectorCacheSize;
		int blockStart = partStart + warpSize*(blockIdx%blockPerPart);
		if(blockIdx%blockPerPart == 0)
			extraRows = partEnd - (partStart + warpSize*blockPerPart);

		for(int i = 0; i < warpSize; ++i){
			writedInRowER = 0;
			writedInRowELL= 0;
			int rowVal = blockStart + i;
			if(rowVal >= realPartStart && rowVal < partEnd){
				if(reorderListER[rowVal] < outputMatrix->numOfRowER){
					rowLocER = reorderListER[rowVal];
					if(rowVecER[rowLocER] != rowVal){
						printf("error at rowVecER\n");
						exit(0);
					}
					blockIdxER = rowLocER/warpSize;
					rowLocInBlockER = rowLocER - blockIdxER*warpSize;
					biasER = biasVecER[blockIdxER];
					widthER = widthVecER[blockIdxER];
				} else {
					if(numInRowER[rowVal] > 0){
						printf("error: got to technical not ER line\n");	
						exit(0);
					}
					widthER = -1;
				}
				for(int j = 0; j < numInRow[rowVal]; ++j){
					int tmpIdx = j + rowIdx[rowVal];	
					if(I[tmpIdx] != rowVal){
						printf("row val check failed\n");
						exit(0);
					}
					if(J[tmpIdx] < fetchEnd && J[tmpIdx] >= partStart){
						colBlockELL[biasBlockELL+i+writedInRowELL*warpSize] = (int16_t)(J[tmpIdx] - partStart);
						valBlockELL[biasBlockELL+i+writedInRowELL*warpSize] = V[tmpIdx];
						writedInRowELL += 1;	
						if(writedInRowELL > widthBlockELL){
							printf("write more elements to blockELL than its width %d\n", widthBlockELL);
							exit(1);
						}
					} else {
						if(widthER < 0){
							printf("go to noexist irrigular line\n");
							exit(0);
						}	
						if(writedInRowER >= widthER){
							printf("error at rowVecER inner with rowVal %d\n", rowVal);
							exit(0);
						}
						colER[biasER+rowLocInBlockER+writedInRowER*warpSize] = J[tmpIdx];
						valER[biasER+rowLocInBlockER+writedInRowER*warpSize] = V[tmpIdx];
						writedInRowER += 1;
					}
				}
				while(writedInRowELL < widthBlockELL){
					colBlockELL[biasBlockELL+i+writedInRowELL*warpSize] = 0;
					valBlockELL[biasBlockELL+i+writedInRowELL*warpSize] = 0;
					wasteElement += 1;
					writedInRowELL+=1;
				}
				double val1, val2;
			} else {
				for(int j = 0; j < widthBlockELL; ++j){
					colBlockELL[biasBlockELL+i+j*warpSize] = 0;
					valBlockELL[biasBlockELL+i+j*warpSize] = 0;
					wasteElement += 1;
				}
			}	
			
		}	
		if(extraRows > 0){
			for(int rowVal = (partStart + warpSize*blockPerPart); rowVal < partEnd; ++rowVal){
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
					valER[biasER+rowLocInBlockER+j*warpSize]=V[rowIdx[rowVal]+j];
					irregular+=1;
					if(j > widthER){
						printf("error at ER width%d\n");
						exit(0);
					}
				}
			}
		}
		extraRows = 0;
	}

	printf("wasteElement is %d\n", wasteElement);
}

/*parameters, first line: output of HYB related parameters, 2nd line: output of HYB matrices, 3rd line:input of local COO matrix 
  4th line: ER index (didn't include values), 5th line: input variables*/

void COO2EHYB(matrixCOO* inputMatrix, 
		matrixEHYB* outputMatrix,
		int* sizeBlockELL, 
		int* sizeER) 
		
{ 

	int* reorderListER = (int*)calloc(inputMatrix->dimension, sizeof(int));
	int* numInRowER = (int*)calloc(inputMatrix->dimension, sizeof(int));
	int* tmpBoundary = (int*)calloc(inputMatrix->nParts, sizeof(int));
	outputMatrix->dimension = inputMatrix->dimension;
	outputMatrix->nParts = inputMatrix->nParts;
	outputMatrix->vectorCacheSize = inputMatrix->vectorCacheSize;
	outputMatrix->kernelPerPart = inputMatrix->kernelPerPart;
	outputMatrix->partBoundary = inputMatrix->partBoundary;
	outputMatrix->widthVecBlockELL = (int16_t*)calloc(outputMatrix->nParts*(inputMatrix->vectorCacheSize/warpSize), sizeof(int16_t));
	outputMatrix->biasVecBlockELL = (int*)calloc(outputMatrix->nParts*(inputMatrix->vectorCacheSize/warpSize), sizeof(int));
	/*block_boundary will not be used outside this function*/
	vecsGenBlockELL(inputMatrix, outputMatrix, tmpBoundary, reorderListER, numInRowER);

	*sizeBlockELL = 0;
	for(int i = 0; i < (inputMatrix->vectorCacheSize/warpSize)*outputMatrix->nParts; ++i){
		outputMatrix->biasVecBlockELL[i] = *sizeBlockELL; 
		*sizeBlockELL += warpSize*outputMatrix->widthVecBlockELL[i]; 
	}
	outputMatrix->valBlockELL = (double*)calloc((*sizeBlockELL), sizeof(double));
	outputMatrix->colBlockELL = (int16_t*)calloc((*sizeBlockELL), sizeof(int16_t));
	int blockNumER = ceil((double (outputMatrix->numOfRowER))/warpSize);
	outputMatrix->biasVecER = (int*)calloc(blockNumER, sizeof(int));
	outputMatrix->widthVecER = (int16_t *)calloc(blockNumER, sizeof(int16_t));

	vecsGenER(outputMatrix, numInRowER, reorderListER);
	for(int i = 1; i < blockNumER; ++i){
		outputMatrix->biasVecER[i] = outputMatrix->biasVecER[i-1] + warpSize*(outputMatrix->widthVecER[i-1]); 
	}
	*sizeER = 0;	
	for(int i = 0; i < blockNumER; ++i){
		*sizeER += warpSize*outputMatrix->widthVecER[i];
	}

	outputMatrix->valER = (double*) calloc((*sizeER), sizeof(double));
	outputMatrix->colER = (int*) calloc((*sizeER), sizeof(int));

	COO2EHYBCore(inputMatrix, 
			outputMatrix,
			tmpBoundary,
			reorderListER,
			numInRowER);
	
	//checkBlockELL(outputMatrix);
	outputMatrix->reorderListER = reorderListER;
	free(numInRowER);
	free(tmpBoundary);
}

