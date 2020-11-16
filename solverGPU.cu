#include "kernel.h"
#include "solver.h"
#include "convert.h"
#include "reordering.h"
#include "test.h"

static void cudaMallocTransDataEHYB(matrixEHYB* localMatrix, matrixEHYB* localMatrix_d, 
		const int sizeBlockELL, const int sizeER){


	localMatrix_d->dimension = localMatrix->dimension;
	localMatrix_d->numOfRowER = localMatrix->numOfRowER;
	localMatrix_d->nParts = localMatrix->nParts;
	int blockNumER = ceil(((float) localMatrix->numOfRowER)/warpSize);

    cudaMalloc((void **) &(localMatrix_d->biasVecBlockELL), localMatrix->nParts*blockPerPart*sizeof(int));
    cudaMalloc((void **) &(localMatrix_d->widthVecBlockELL), localMatrix->nParts*blockPerPart*sizeof(int));
    cudaMalloc((void **) &(localMatrix_d->partBoundary), (localMatrix->nParts+1)*sizeof(int));
    cudaMalloc((void **) &(localMatrix_d->valBlockELL), sizeBlockELL*sizeof(double));
    cudaMalloc((void **) &(localMatrix_d->colBlockELL), sizeBlockELL*sizeof(int));

    cudaMalloc((void **) &(localMatrix_d->rowVecER), localMatrix_d->numOfRowER*sizeof(int));
    cudaMalloc((void **) &(localMatrix_d->biasVecER), blockNumER*sizeof(int));
    cudaMalloc((void **) &(localMatrix_d->widthVecER), blockNumER*sizeof(int));
    cudaMalloc((void **) &(localMatrix_d->colER), sizeER*sizeof(int));
    cudaMalloc((void **) &(localMatrix_d->valER), sizeER*sizeof(double));

    cudaMemcpy(localMatrix_d->biasVecBlockELL, localMatrix->biasVecBlockELL, localMatrix->nParts*blockPerPart*sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(localMatrix_d->widthVecBlockELL, localMatrix->widthVecBlockELL, localMatrix->nParts*blockPerPart*sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(localMatrix_d->partBoundary, localMatrix->partBoundary, (localMatrix->nParts+1)*sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(localMatrix_d->valBlockELL, localMatrix->valBlockELL, sizeBlockELL*sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(localMatrix_d->colBlockELL, localMatrix->colBlockELL, sizeBlockELL*sizeof(int), cudaMemcpyHostToDevice);

    cudaMemcpy(localMatrix_d->rowVecER, localMatrix->rowVecER, localMatrix_d->numOfRowER*sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(localMatrix_d->biasVecER, localMatrix->biasVecER, blockNumER*sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(localMatrix_d->widthVecER, localMatrix->widthVecER, blockNumER*sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(localMatrix_d->colER, localMatrix->colER, sizeER*sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(localMatrix_d->valER, localMatrix->valER, sizeER*sizeof(double), cudaMemcpyHostToDevice);
}
extern "C"
void solverGPuUnprecondEHYB(matrixCOO* localMatrix, 
		const double *vectorIn, double *vectorOut,  
		const int MAXIter, int *realIter)
{
	//This function treat y as input and x as output, (solve the equation Ax=y) y is the vector we already known, x is the vector we are looking for
	double dotp0,dotr0,dotr1,doth;

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
	cublasHandle_t handle;
	cublasCreate(&handle);

	double *bp_d, *pk_d, *rk_d, *vectorOut_d;
	size_t size1 = dimension*sizeof(double);
	cudaMalloc((void **) &bp_d,size1);
	cudaMalloc((void **) &pk_d,size1);
	cudaMalloc((void **) &rk_d,size1);
	cudaMalloc((void **) &vectorOut_d,size1);
	//double *x=(double *) malloc(size1);
	double threshold=0.0000001;
	int iter=0;
	double const1 = 1.0;
	double error, alphak, _alphak, gamak;
	error=1000;
	//initialize
	doth=0;
    cudaMemcpy(pk_d, vectorIn, dimension*sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(rk_d, vectorIn, dimension*sizeof(double), cudaMemcpyHostToDevice);
	for (int i=0;i<dimension;i++) {
		doth=doth+vectorIn[i]*vectorIn[i];
	}
	struct timeval start1, end1;
	double *bp=(double *) malloc(size1);
	double *bp_g =(double *) malloc(size1);
	double *pk=(double *) malloc(size1);
	double *rk=(double *) malloc(size1);

	//double *x=(double *) malloc(size1);
	error=1000;
	//initialize
	for (int i=0;i<dimension;i++)
	{
		pk[i]=vectorIn[i];
		rk[i]=vectorIn[i];
		vectorOut[i]=0;
		bp[i]=0;
		bp_g[i]=0;
	}
	gettimeofday(&start1, NULL);
	while (error>threshold&&iter<MAXIter){
		dotp0=0;
		dotr0=0;
		dotr1=0;
		int errorIdx = 0;
		double compareError;
		cudaMemset(bp_d, 0, size1);
		matrixVectorEHYB(&localMatrixEHYB_d, pk_d, bp_d, -1);
		cublasDdot(handle,dimension,bp_d,1,pk_d,1,&dotp0);
		cublasDdot(handle,dimension,rk_d,1,rk_d,1,&dotr0);
			
		alphak=dotr0/dotp0;
		_alphak = -alphak;
		
		cublasDaxpy(handle,dimension,&alphak,pk_d,1,vectorOut_d,1);
		cublasDaxpy(handle,dimension,&_alphak,bp_d,1,rk_d,1);
		cublasDdot(handle,dimension,rk_d,1,rk_d,1,&dotr1);
		
		gamak=dotr1/dotr0;

		cublasDscal(handle,dimension,&gamak, pk_d,1);
		cublasDaxpy(handle,dimension,&const1, rk_d, 1, pk_d, 1);
		
		//printf("at iter %d, alphak is %f, gamak is %f\n",iter, alphak,gamak);
		error=sqrt(dotr1)/sqrt(doth);
		//error_track[iter]=error;
		//printf("error at %d is %f\n",iter, error);
		iter++;
	}
	cudaMemcpy(vectorOut, vectorOut_d, dimension*sizeof(double), cudaMemcpyDeviceToHost);
	gettimeofday(&end1, NULL);	
	double timeByMs=((end1.tv_sec * 1000000 + end1.tv_usec)-(start1.tv_sec * 1000000 + start1.tv_usec))/1000;
	printf("iter is %d, time is %f ms, GPU Gflops is %f, under estimate flops is %f\n ",iter, timeByMs, 
			(1e-9*(totalNum*2+13*dimension)*1000*iter)/timeByMs, (1e-9*(totalNum*2)*1000*iter)/timeByMs);
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
		const double *vector_in, double *vector_out,  
		const int MAXIter, int *realIter,  const cb_s cb,
		const int partSize, const int* partBoundary)
{
	//exampine the performance using cusparse library functions with
	//CSR format
	//double dotp0,dotr0,dotr1,doth;
	double dotp0,dotr0,dotr1,doth;
	int dimension, totalNum; 
    int *rowIdx, *J; 
    double* V;
    dimension = localMatrix->dimension; 
    totalNum = localMatrix->totalNum; 

	rowIdx = localMatrix->rowIdx; 
    J = localMatrix->J;
    V = localMatrix->V;
	
	int* col_d;
	int* rowIdx_d;
	double *V_d;
	cublasHandle_t handleBlas;
	cublasCreate(&handleBlas);
	cusparseHandle_t handleSparse;
	cusparseCreate(&handleSparse);

	double *bp_d, *pk_d, *rk_d, *vector_out_d;
	size_t size1=dimension*sizeof(double);
	cudaMalloc((void **) &bp_d,size1);
	cudaMalloc((void **) &pk_d,size1);
	cudaMalloc((void **) &rk_d,size1);
	cudaMalloc((void **) &rowIdx_d,size1);
	cudaMalloc((void **) &vector_out_d,size1);
	cudaMalloc((void **) &col_d,totalNum*sizeof(int));
	cudaMalloc((void **) &V_d,totalNum*sizeof(double));
	//double *x=(double *) malloc(size1);
	double threshold=0.0000001;
	int iter=0;
	double const1 = 1.0;
	double error, alphak, _alphak, gamak;
	error=1000;
	//initialize
	doth=0;
    cudaMemcpy(pk_d, vector_in, dimension*sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(rk_d, vector_in, dimension*sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(rowIdx_d, rowIdx, dimension*sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(col_d, J, totalNum*sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(V_d, V, totalNum*sizeof(double), cudaMemcpyHostToDevice);
	for (int i=0;i<dimension;i++) {
		doth=doth+vector_in[i]*vector_in[i];
	}
	struct timeval start1, end1;
	double *bp=(double *) malloc(size1);
	double *bp_g =(double *) malloc(size1);
	double *pk=(double *) malloc(size1);
	double *rk=(double *) malloc(size1);

	//double *bp_dt = (double *) malloc(size1);
	//double *pk_dt = (double *) malloc(size1);
	//double *rk_dt = (double *) malloc(size1);
	//double *x=(double *) malloc(size1);
	error=1000;
	//initialize
	for (int i=0;i<dimension;i++)
	{
		pk[i]=vector_in[i];
		rk[i]=vector_in[i];
		vector_out[i]=0;
		bp[i]=0;
	}
	//if BSR doing the format change 
	//cusparseStatus_tcusparseDcsr2gebsr_bufferSize(handle, dir, m, n, descrA, csrValA, csrRowPtrA, 
	//		csrColIndA, rowBlockDim, colBlockDim, pBufferSize);
	cusparseOperation_t transA = CUSPARSE_OPERATION_NON_TRANSPOSE;
	cusparseMatDescr_t descr = 0;
	int status = cusparseCreateMatDescr(&descr);
	if (status != CUSPARSE_STATUS_SUCCESS ) {
		exit(0);	
	}
	cusparseSetMatType (descr, CUSPARSE_MATRIX_TYPE_GENERAL);
	cusparseSetMatIndexBase (descr, CUSPARSE_INDEX_BASE_ZERO);
	gettimeofday(&start1, NULL);
	double one = 1.0;
	double zero = 0.0;
	while (error>threshold&&iter<MAXIter){
		dotp0=0;
		dotr0=0;
		dotr1=0;
		//int errorIdx = 0;
		//double compareError;
		
		cudaMemset(bp_d, 0, size1);
		cusparseStatus_t smpvStatus = 
		cusparseDcsrmv(handleSparse,
				transA,
				dimension,
				dimension,
				totalNum,
				&one,
				descr,
				V_d,
				rowIdx_d,
				col_d,
				pk_d,
				&zero,
				bp_d);

		cublasDdot(handleBlas,dimension,bp_d,1,pk_d,1,&dotp0);
		cublasDdot(handleBlas,dimension,rk_d,1,rk_d,1,&dotr0);
			
		alphak=dotr0/dotp0;
		_alphak = -alphak;
		
		cublasDaxpy(handleBlas,dimension,&alphak,pk_d,1,vector_out_d,1);
		cublasDaxpy(handleBlas,dimension,&_alphak,bp_d,1,rk_d,1);
		cublasDdot(handleBlas,dimension,rk_d,1,rk_d,1,&dotr1);
		
		gamak=dotr1/dotr0;

		cublasDscal(handleBlas,dimension,&gamak,pk_d,1);
		cublasDaxpy(handleBlas,dimension,&const1, rk_d, 1, pk_d, 1);
		
		//printf("at iter %d, alphak is %f, gamak is %f\n",iter, alphak,gamak);
		error=sqrt(dotr1)/sqrt(doth);
		//error_track[iter]=error;
		//printf("error at %d is %f\n",iter, error);
		iter++;
	}
	cudaMemcpy(vector_out, vector_out_d, dimension*sizeof(double), cudaMemcpyDeviceToHost);
	gettimeofday(&end1, NULL);	
	double timeByMs=((end1.tv_sec * 1000000 + end1.tv_usec)-(start1.tv_sec * 1000000 + start1.tv_usec))/1000;
	printf("iter is %d, time is %f ms, GPU Gflops is %f, under estimate flops is %f\n ",iter, timeByMs, 
			(1e-9*(totalNum*2+13*dimension)*1000*iter)/timeByMs, (1e-9*(totalNum*2)*1000*iter)/timeByMs);

}

