#include "kernel.h"
#include "solver.h"
#include "convert.h"
#include "Partition.h"

#define warpSize 32

static bool checkRowErr(double* V, double* valBlockELL, double* valER, 
		int numInRowCOO, int widthBlockELL, int widthER, 
		int biasCOO, int biasBlockELL, int biasER,
		double *val1, double *val2){
	*val1 = 0; 
	*val2 = 0;
	for(int i = 0; i < numInRowCOO; ++i){
		*val1 += V[biasCOO + i];
	}	
	for(int i = 0; i < widthBlockELL; ++i){
		*val2 += valBlockELL[biasBlockELL + i*warpSize];
	}
	if(widthER > 0){
		for(int i = 0; i < widthER; ++i){
			*val2 += valER[biasER + i*warpSize];
		}
	}
	if (*val1 - *val2 < .001 && *val1 - *val2 > -.001)
		return false;
	else 
		return true;
}
static void checkBlockELL(matrixEHYB* inputMatrix){
	
	for(int part = 0; part < inputMatrix->nParts; ++part){
		int vecStart = inputMatrix->partBoundary[part];
		int vecEnd = vecStart + vectorCacheSize;
		int blockPartStart = part*blockPerPart;
		for(int block = 0; block < blockPerPart; ++block){
			int blockBias = inputMatrix->biasVecBlockELL[blockPartStart + block];
			int width = inputMatrix->widthVecBlockELL[blockPerPart + block];
			int row = vecStart + block*warpSize;
			for(int i = 0; i < 32; ++i){
				if(row + i < vecEnd){
					for(int n = 0; n < width; ++n){
						int dataIdx = blockBias + i + n*warpSize;
						int col = inputMatrix->colBlockELL[dataIdx];
						double val = inputMatrix->valBlockELL[dataIdx];
						if(val == 1)
							printf("error found!\n");
					}
				}
			}
		}
	}
}
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
	int maxCol = inputMatrix->maxCol;
	int nonZeros;
	int colTh = 0, colAccum = 0; 
	int *colHist = ( int*)calloc(maxCol, sizeof( int));
	int *numInRow = inputMatrix->numInRow;
	int *numInRow2 = inputMatrix->numInRow2; 

	int extraRows = 0;
	int numOfRowER = 0;
	int actualRow;
	int toER = 0;

	for(int partIdx = 0; partIdx < inputMatrix->nParts; ++partIdx){
		int partStart = inputMatrix->partBoundary[partIdx];
		int partEnd = inputMatrix->partBoundary[partIdx + 1];
		if(partEnd > partStart + blockPerPart*warpSize)
			extraRows = partEnd - partStart - blockPerPart*warpSize;
		else
			extraRows = 0;

		for(int iter = 0; iter < blockPerPart; ++iter){
			blockIdx = iter + blockPerPart*partIdx;
			int blockStart = partStart + iter*warpSize;
			memset(colHist, 0, maxCol*sizeof(int));
			actualRow = 0;
			for(int row = blockStart; (row < blockStart + warpSize && row < partEnd); ++row){
				if(numInRow2[row] > 0){
					colHist[numInRow2[row]- 1] += 1;
					actualRow += 1;	
				}
				if(numInRow2[row] != numInRow[row]){
					numOfRowER += 1;
					numInRowER[row] = numInRow[row] - numInRow2[row];
					toER += numInRowER[row];
				}
			}
			colTh = ceil(float(actualRow)*0.75);
			colAccum = 0;
			for( int i = 0; i < maxCol; ++i){
				colAccum += colHist[i];
				if(colAccum > colTh){
					numCols = i + 1;
					break;
				}
			}
			if(numCols == 0){
				printf("error at ELL_block convert rodr\n");
				exit(0);
			}
			for(int row = blockStart; (row < blockStart + warpSize && row < partEnd); ++row){
				if(numInRow2[row] > numCols){
					if(numInRowER[row] == 0)
						numOfRowER += 1;
					numInRowER[row] += numInRow2[row] - numCols;
					toER += numInRow2[row] - numCols;
				}
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
	printf("toER is %d\n", toER);
	outputMatrix->numOfRowER = numOfRowER;
	outputMatrix->rowVecER = (int*)malloc(numOfRowER*sizeof(int));
	outputMatrix->biasVecER = (int*)malloc(ceil(((float) numOfRowER)/warpSize)*sizeof(int));
	free(colHist);

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
		int* reorderListER,
		int* numInRowER)
{ 
	int* partBoundary = inputMatrix->partBoundary;
	
	int* rowIdx = inputMatrix->rowIdx;
	int* numInRow = inputMatrix->numInRow;
	int* I = inputMatrix->I;
	int* J = inputMatrix->J;
	double* V = inputMatrix->V;

	int* widthVecBlockELL = outputMatrix->widthVecBlockELL;
	int* biasVecBlockELL = outputMatrix->biasVecBlockELL;
	int* colBlockELL = outputMatrix->colBlockELL; 
	double* valBlockELL = outputMatrix->valBlockELL; 

	int* widthVecER = outputMatrix->widthVecER;
	int* rowVecER = outputMatrix->rowVecER; 
	int* biasVecER = outputMatrix->biasVecER;
	int* colER = outputMatrix->colER; 
	double* valER = outputMatrix->valER; 

	int partIdx, partStart, partEnd, fetchEnd, extraRows;
	int rowLocER, blockIdxER, rowLocInBlockER, biasER, widthER;     	
	int biasBlockELL, widthBlockELL;     	
	int irregular=0;
	int writedInRowELL = 0;
	int writedInRowER = 0;
	int padding = 0;
	int nBlocks = inputMatrix->nParts*blockPerPart;
	

	for(int blockIdx = 0; blockIdx < nBlocks; ++blockIdx){
		widthBlockELL = widthVecBlockELL[blockIdx];
		biasBlockELL = biasVecBlockELL[blockIdx]; 
		partIdx = blockIdx/blockPerPart;
		partStart = partBoundary[partIdx];
		partEnd = partBoundary[partIdx + 1];
		fetchEnd =  partStart + vectorCacheSize;
		int blockStart = partStart + warpSize*(blockIdx%blockPerPart);
		if(blockIdx%blockPerPart == 0)
			extraRows = partEnd - (partStart + warpSize*blockPerPart);

		for(int i = 0; i < warpSize; ++i){
			writedInRowER = 0;
			writedInRowELL= 0;
			int rowVal = blockStart + i;
			if(rowVal < partEnd){
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
					}
					widthER = -1;
				}
				for(int j = 0; j < numInRow[rowVal]; ++j){
					int tmpIdx = j + rowIdx[rowVal];	
					if(I[tmpIdx] != rowVal){
						printf("row val check failed\n");
						exit(0);
					}
					if(J[tmpIdx] < fetchEnd && J[tmpIdx] >= partStart && writedInRowELL < widthBlockELL){
						colBlockELL[biasBlockELL+i+writedInRowELL*warpSize] = J[tmpIdx];
						valBlockELL[biasBlockELL+i+writedInRowELL*warpSize] = V[tmpIdx];
						writedInRowELL += 1;	
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
					colBlockELL[biasBlockELL+i+writedInRowELL*warpSize] = partStart;
					valBlockELL[biasBlockELL+i+writedInRowELL*warpSize] = 0;
					writedInRowELL+=1;
				}
				double val1, val2;
			} else {
				for(int j = 0; j < widthBlockELL; ++j){
					colBlockELL[biasBlockELL+i+j*warpSize] = partStart;
					valBlockELL[biasBlockELL+i+j*warpSize] = 0;
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
						printf("error at ER%d\n");
						exit(0);
					}
				}
			}
		}
		extraRows = 0;
	}

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
	outputMatrix->dimension = inputMatrix->dimension;
	outputMatrix->nParts = inputMatrix->nParts;
	outputMatrix->partBoundary = inputMatrix->partBoundary;
	outputMatrix->widthVecBlockELL = (int*)calloc(outputMatrix->nParts*blockPerPart, sizeof(int));
	outputMatrix->biasVecBlockELL = (int*)calloc(outputMatrix->nParts*blockPerPart, sizeof(int));
	/*block_boundary will not be used outside this function*/
	vecsGenBlockELL(inputMatrix, outputMatrix, reorderListER, numInRowER);

	*sizeBlockELL = 0;
	for(int i = 0; i < blockPerPart*outputMatrix->nParts; ++i){
		outputMatrix->biasVecBlockELL[i] = *sizeBlockELL; 
		*sizeBlockELL += warpSize*outputMatrix->widthVecBlockELL[i]; 
	}
	outputMatrix->valBlockELL = (double*)calloc((*sizeBlockELL), sizeof(double));
	outputMatrix->colBlockELL = (int*)calloc((*sizeBlockELL), sizeof(int));
	int blockNumER = ceil((float (outputMatrix->numOfRowER))/warpSize);
	outputMatrix->biasVecER = (int*)calloc(blockNumER, sizeof(int));
	outputMatrix->widthVecER = (int*)calloc(blockNumER, sizeof(int));

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
			reorderListER,
			numInRowER);
	
	//checkBlockELL(outputMatrix);
	outputMatrix->reorderListER = reorderListER;
	free(numInRowER);
}

