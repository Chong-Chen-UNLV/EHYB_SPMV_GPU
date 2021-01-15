#include "kernel.h"
#include "spmv.h"
#include "convert.h"
#include "reordering.h"

static void cudaMallocTransDataEHYB(matrixEHYB* localMatrix, matrixEHYB* localMatrix_d, 
		const int sizeBlockELL, const int sizeER){


	localMatrix_d->dimension = localMatrix->dimension;
	localMatrix_d->numOfRowER = localMatrix->numOfRowER;
	localMatrix_d->nParts = localMatrix->nParts;
	int blockNumER = ceil(((float) localMatrix->numOfRowER)/warpSize);

    cudaMalloc((void **) &(localMatrix_d->biasVecBlockELL), localMatrix->nParts*blockPerPart*sizeof(int));
    cudaMalloc((void **) &(localMatrix_d->widthVecBlockELL), localMatrix->nParts*blockPerPart*sizeof(int));
    cudaMalloc((void **) &(localMatrix_d->partBoundary), (localMatrix->nParts+1)*sizeof(int));
    cudaMalloc((void **) &(localMatrix_d->valBlockELL), sizeBlockELL*sizeof(float));
    cudaMalloc((void **) &(localMatrix_d->colBlockELL), sizeBlockELL*sizeof(int16_t));

    cudaMalloc((void **) &(localMatrix_d->rowVecER), localMatrix_d->numOfRowER*sizeof(int));
    cudaMalloc((void **) &(localMatrix_d->biasVecER), blockNumER*sizeof(int));
    cudaMalloc((void **) &(localMatrix_d->widthVecER), blockNumER*sizeof(int));
    cudaMalloc((void **) &(localMatrix_d->colER), sizeER*sizeof(int));
    cudaMalloc((void **) &(localMatrix_d->valER), sizeER*sizeof(float));

    cudaMemcpy(localMatrix_d->biasVecBlockELL, localMatrix->biasVecBlockELL, localMatrix->nParts*blockPerPart*sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(localMatrix_d->widthVecBlockELL, localMatrix->widthVecBlockELL, localMatrix->nParts*blockPerPart*sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(localMatrix_d->partBoundary, localMatrix->partBoundary, (localMatrix->nParts+1)*sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(localMatrix_d->valBlockELL, localMatrix->valBlockELL, sizeBlockELL*sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(localMatrix_d->colBlockELL, localMatrix->colBlockELL, sizeBlockELL*sizeof(int16_t), cudaMemcpyHostToDevice);

    cudaMemcpy(localMatrix_d->rowVecER, localMatrix->rowVecER, localMatrix_d->numOfRowER*sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(localMatrix_d->biasVecER, localMatrix->biasVecER, blockNumER*sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(localMatrix_d->widthVecER, localMatrix->widthVecER, blockNumER*sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(localMatrix_d->colER, localMatrix->colER, sizeER*sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(localMatrix_d->valER, localMatrix->valER, sizeER*sizeof(float), cudaMemcpyHostToDevice);
}
extern "C"
void spmvGPuEHYB(matrixCOO* localMatrix, 
		const float *vectorIn, float *vectorOut,  
		const int MAXIter, int *realIter)
{
	//This function treat y as input and x as output, (solve the equation Ax=y) y is the vector we already known, x is the vector we are looking for

	int sizeBlockELL, sizeER;
	int dimension = localMatrix->dimension;
	int totalNum = localMatrix->totalNum;
	matrixEHYB localMatrixEHYB, localMatrixEHYB_d;

	COO2EHYB(localMatrix, 
			&localMatrixEHYB,
			&sizeBlockELL,
			&sizeER);

	cudaMallocTransDataEHYB(&localMatrixEHYB,
			&localMatrixEHYB_d, 
			sizeBlockELL,
			sizeER);
	printf("sizeER is %d\n", sizeER);

	float *vectorIn_d, *vectorOut_d;
	int16_t *biasIdxBlock_d;
	size_t size1 = dimension*sizeof(float);
	//cudaMalloc((void **) &biasIdxBlock_d, sizeof(uint16_t)*localMatrix->nParts);
	cudaMalloc((void **) &vectorOut_d,size1);
	cudaMalloc((void **) &vectorIn_d,size1);
	//float *x=(float *) malloc(size1);
	int iter=0;
	struct timeval start1, end1;

	//float *x=(float *) malloc(size1);
	//initialize
	
	gettimeofday(&start1, NULL);
    cudaMemcpy(vectorIn_d, vectorIn, dimension*sizeof(float), cudaMemcpyHostToDevice);
	while (iter<MAXIter){
		matrixVectorEHYB(&localMatrixEHYB_d, vectorIn_d, vectorOut_d, -1);
		iter++;
	}
	cudaMemcpy(vectorOut, vectorOut_d, dimension*sizeof(float), cudaMemcpyDeviceToHost);
	gettimeofday(&end1, NULL);	
	float timeByMs=((end1.tv_sec * 1000000 + end1.tv_usec)-(start1.tv_sec * 1000000 + start1.tv_usec))/1000;
	printf("iter is %d, time is %f ms, GPU Gflops is %f\n ",iter, timeByMs, 
			(1e-9*(totalNum*2)*1000*iter)/timeByMs );
	cudaFree(localMatrixEHYB_d.valER);
	cudaFree(localMatrixEHYB_d.colER);
	cudaFree(localMatrixEHYB_d.biasVecER);
	cudaFree(localMatrixEHYB_d.widthVecER);
	cudaFree(localMatrixEHYB_d.rowVecER);
	cudaFree(localMatrixEHYB_d.biasVecBlockELL);
	cudaFree(localMatrixEHYB_d.widthVecBlockELL);
	cudaFree(localMatrixEHYB_d.colBlockELL);
	cudaFree(localMatrixEHYB_d.valBlockELL);
	cudaFree(localMatrixEHYB_d.partBoundary);
}

void solverGPuUnprecondCUSPARSE(matrixCOO* localMatrix, 
		const float *vector_in, float *vector_out,  
		const int MAXIter)
{
	//exampine the performance using cusparse library functions with
	//CSR format
	//float dotp0,dotr0,dotr1,doth;
	int dimension, totalNum; 
    int *rowIdx, *J; 
    float* V;
    dimension = localMatrix->dimension; 
    totalNum = localMatrix->totalNum; 

	rowIdx = localMatrix->rowIdx; 
    J = localMatrix->J;
    V = localMatrix->V;
	
	int* col_d;
	int* rowIdx_d;
	float *V_d;

	float *vector_in_d, *vector_out_d;
	size_t size1=dimension*sizeof(float);
	
	cudaMalloc((void **) &rowIdx_d, (dimension+1)*sizeof(float));
	cudaMalloc((void **) &vector_out_d,size1);
	cudaMalloc((void **) &vector_in_d,size1);
	cudaMalloc((void **) &col_d,totalNum*sizeof(int));
	cudaMalloc((void **) &V_d,totalNum*sizeof(float));
	//float *x=(float *) malloc(size1);
	int iter=0;
	float const1 = 1.0;
	//initialize
   	if(cudaSuccess != cudaMemcpy(rowIdx_d, rowIdx, (dimension+1)*sizeof(int), cudaMemcpyHostToDevice)) printf("error1\n");
    if(cudaSuccess !=cudaMemcpy(col_d, J, totalNum*sizeof(int), cudaMemcpyHostToDevice)) printf("error2\n");
    if(cudaSuccess !=cudaMemcpy(V_d, V, totalNum*sizeof(float), cudaMemcpyHostToDevice)) printf("error3\n");
	
	struct timeval start1, end1;
	
	//if BSR doing the format change 
	//cusparseStatus_tcusparseDcsr2gebsr_bufferSize(handle, dir, m, n, descrA, csrValA, csrRowPtrA, 
	//		csrColIndA, rowBlockDim, colBlockDim, pBufferSize);
	cusparseHandle_t handleSparse;
	cusparseCreate(&handleSparse);
	cusparseOperation_t transA = CUSPARSE_OPERATION_NON_TRANSPOSE;
	cusparseMatDescr_t descr = 0;
	int status = cusparseCreateMatDescr(&descr);
	if (status != CUSPARSE_STATUS_SUCCESS ) {
		exit(0);	
	}
	cusparseSetMatType (descr, CUSPARSE_MATRIX_TYPE_GENERAL);
	cusparseSetMatIndexBase (descr, CUSPARSE_INDEX_BASE_ZERO);
	gettimeofday(&start1, NULL);
	float one = 1.0;
	float zero = 0.0;
    if(cudaSuccess != cudaMemcpy(vector_in_d, vector_in, dimension*sizeof(float), cudaMemcpyHostToDevice)) printf("error4\n");
	while (iter<MAXIter){
		//int errorIdx = 0;
		//float compareError;
		
		cusparseStatus_t smpvStatus = 
		cusparseScsrmv(handleSparse,
				transA,
				dimension,
				dimension,
				totalNum,
				&one,
				descr,
				V_d,
				rowIdx_d,
				col_d,
				vector_in_d,
				&zero,
				vector_out_d);
		iter++;
	}
	cudaMemcpy(vector_out, vector_out_d, dimension*sizeof(float), cudaMemcpyDeviceToHost);
	gettimeofday(&end1, NULL);	
	float timeByMs=((end1.tv_sec * 1000000 + end1.tv_usec)-(start1.tv_sec * 1000000 + start1.tv_usec))/1000;
	printf("iter is %d, time is %f ms, GPU csrmv Gflops is %f\n ",iter, timeByMs, (1e-9*(totalNum*2)*1000*iter)/timeByMs);
			

}

