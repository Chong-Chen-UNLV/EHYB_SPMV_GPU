#include "kernel.h"
#include "solver.h"
#include "convert.h"
#include "test.h"

static inline int calPadding(int dimension,
		const int* numInRow, const int* part_boundary, 	
		int* ELL_block_cols_vec, const int part_size){

	int padded = 0;
	for(int i = 0; i < part_size; ++i){
		for(int j = part_boundary[i]; j < part_boundary[i+1]; ++j){
			if(j < part_boundary[i] + block_per_part*ELL_threadSize){
				int blockIdx = i*block_per_part + (j-part_boundary[i])/ELL_threadSize;
				if(numInRow[j] < ELL_block_cols_vec[blockIdx])	
					padded += (ELL_block_cols_vec[blockIdx] - numInRow[j]);
			}
		}
	}
	return padded;
}

static void cudaMallocTransDataBlockELL(int** col_d, double** V_d, 
				int** ELL_block_cols_vec_d, int** ELL_block_bias_vec_d,
                int* colELL, double* matrixELL,
				int* ELL_block_cols_vec, int* ELL_block_bias_vec,
				int dimension, int blocks){
                
	int ELL_size = ELL_block_bias_vec[blocks];
	//cudaSetDevice(cuda_device);
    cudaMalloc((void **) col_d, ELL_size*sizeof(int));
    cudaMalloc((void **) V_d, ELL_size*sizeof(double));
    cudaMalloc((void **) ELL_block_cols_vec_d, blocks*sizeof(int));
    cudaMalloc((void **) ELL_block_bias_vec_d, (blocks + 1)*sizeof(int));

    cudaMemcpy(*col_d,colELL,ELL_size*sizeof(int),cudaMemcpyHostToDevice);
    cudaMemcpy(*V_d,matrixELL,ELL_size*sizeof(double),cudaMemcpyHostToDevice);
    cudaMemcpy(*ELL_block_cols_vec_d, ELL_block_cols_vec, blocks*sizeof(int),cudaMemcpyHostToDevice);
    cudaMemcpy(*ELL_block_bias_vec_d, ELL_block_bias_vec, (blocks + 1)*sizeof(int),cudaMemcpyHostToDevice);
}

void solverGPuUnprecondEHYB(matrixCOO_S* localMatrix, 
		const double *vector_in, double *vector_out,  
		const int MAXIter, int *realIter,  const cb_s cb,
		const int part_size, const int* part_boundary)
{
	//This function treat y as input and x as output, (solve the equation Ax=y) y is the vector we already known, x is the vector we are looking for
	double dotp0,dotr0,dotr1,doth;
	int dimension, totalNum, maxRowNum, totalNumCOO; 
    int *row_idx, *numInRow, *I, *J; 
    double *V;
    dimension = localMatrix->dimension; 
    totalNum = localMatrix->totalNum; 
    maxRowNum = localMatrix->maxRowNum; 

	int ELL_blocks;
    int *widthVecBlockELL, *widthVecBlockELL_d;    
    int *biasVecBLockELL, *biasVecBLockELL_d;	
	
    row_idx = localMatrix->rowIdx; 
    numInRow = localMatrix->numInRow; 
    I = localMatrix->I; 
    J = localMatrix->J;
    V = localMatrix->V;
	int *colBlockELL *colER;
	int *rowVecER, *biasVecER, *widthVecER; 
	int* ER_numOfRow;
	double *ER_val;
	double* matrixELL, *V_COO;
	int *col_d;
	double *V_d;
	volatile bool RODR, BLOCK, CACHE;
	int *I_COO_d, *J_COO_d;
	double *V_COO_d;
	BLOCK = cb.BLOCK;
	CACHE = cb.CACHE;
	RODR = cb.RODR;

	matrixEHYbS localMatrixEHYB, localMatrixEHYB_d;
	int *part_boundary_d;
	
        //RODR will change the number of blocks since it provides cache option
	ELL_blocks = part_size*block_per_part; 

    ELL_block_cols_vec = (int*)malloc(ELL_blocks*sizeof(int));  
    ELL_block_bias_vec = (int*)malloc((ELL_blocks +1)*sizeof(int));  
	
	COO2EHYB(&localMatrix, 
		&localMatrixEHYB,
		biasVecBLockELL,
		RODR,
		CACHE);

	cudaMallocTransDataEHYB(localMatrixEHYB,
			&localMatrixEHYB_d);

	int padding = calPadding(dimension, numInRow, part_boundary, 	
		ELL_block_cols_vec, part_size);
	printf("totalNumCOO is %d padding is %d, waste rate on ELL is %f\n", 
			totalNumCOO, padding, ((float)padding)/totalNum);

	cublasHandle_t handle;
	cublasCreate(&handle);

	
	localMatrixEHYB_d.dimension = dimension;
	localMatrixEHYB_d.col_d = col_d;
	localMatrixEHYB_d.V_d = V_d;
	localMatrixEHYB_d.V_COO_d = V_COO_d;
	localMatrixEHYB_d.ELL_block_bias_vec_d = ELL_block_bias_vec_d;
	localMatrixEHYB_d.ELL_block_cols_vec_d = ELL_block_cols_vec_d;

	double *bp_d, *pk_d, *rk_d, *vector_out_d;
	size_t size0=dimension*sizeof(int);
	size_t size1=dimension*sizeof(double);
	cudaMalloc((void **) &bp_d,size1);
	cudaMalloc((void **) &pk_d,size1);
	cudaMalloc((void **) &rk_d,size1);
	cudaMalloc((void **) &vector_out_d,size1);
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
	for (int i=0;i<dimension;i++) {
		doth=doth+vector_in[i]*vector_in[i];
	}
	struct timeval start1, end1;
	double *bp=(double *) malloc(size1);
	double *bp_g =(double *) malloc(size1);
	double *pk=(double *) malloc(size1);
	double *rk=(double *) malloc(size1);

	double *bp_dt = (double *) malloc(size1);
	double *pk_dt = (double *) malloc(size1);
	double *rk_dt = (double *) malloc(size1);
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
	printf("ELL_block_cols_vec[0] is %d\n", ELL_block_cols_vec[0]);	
	gettimeofday(&start1, NULL);
	while (error>threshold&&iter<MAXIter){
		dotp0=0;
		dotr0=0;
		dotr1=0;
		int errorIdx = 0;
		double compareError;
		cudaMemset(bp_d, 0, size1);

		//matrix_vectorTest(localMatrix, vector_in, bp, 0);
		
		matrix_vectorEHYB(&localMatrixEHYB_d, pk_d, bp_d, cb, 0);
		//cudaMemcpy(bp_g, bp_d, dimension*sizeof(double), cudaMemcpyDeviceToHost);
		//for(int i = 0; i < dimension; ++i){
		//	compareError = 	(bp_g[i] - bp[i])/bp[i];
		//	if(errorIdx == 0 && (compareError > 0.0000001 || compareError < -0.0000001)){ 	
		//		printf("bp[%d] of GPU result is %f, test value is %f, difference is %f\n", i, bp_g[i], bp[i], compareError);
		//		if(errorIdx == 0) errorIdx = i;	
		//	}
		//}

		//matrix_vectorTest(localMatrix, vector_in, bp, errorIdx);
		//mimicHYB(ELL_block_cols_vec,
		//		ELL_block_bias_vec,
		//		part_boundary,
		//		colELL,
		//		matrixELL,
		//		I_COO,
		//		J_COO,
		//		V_COO,
		//		totalNumCOO,
		//		vector_in,
		//		errorIdx);
		//
		//cudaMemset(bp_d, 0, size1);
		//matrix_vectorHYB(&localMatrixHYB_d, pk_d, bp_d, cb, errorIdx,
		//	   part_size, part_boundary_d, errorIdx);
		//
		//for (int k=0;k<5;k++) printf("pk %d at iter %d is %f\n",k, iter, pk[k]);
		//printf("CPU start\n");
		
		cublasDdot(handle,dimension,bp_d,1,pk_d,1,&dotp0);
		cublasDdot(handle,dimension,rk_d,1,rk_d,1,&dotr0);
			
		alphak=dotr0/dotp0;
		_alphak = -alphak;
		
		cublasDaxpy(handle,dimension,&alphak,pk_d,1,vector_out_d,1);
		cublasDaxpy(handle,dimension,&_alphak,bp_d,1,rk_d,1);
		cublasDdot(handle,dimension,rk_d,1,rk_d,1,&dotr1);
		
		gamak=dotr1/dotr0;

		cublasDscal(handle,dimension,&gamak,pk_d,1);
		cublasDaxpy(handle,dimension,&const1, rk_d, 1, pk_d, 1);
		
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

void solverGPU_unprecondCUSPARSE(matrixCOO_S* localMatrix, 
		const double *vector_in, double *vector_out,  
		const int MAXIter, int *realIter,  const cb_s cb,
		const int part_size, const int* part_boundary)
{
	//exampine the performance using cusparse library functions with
	//CSR format
	//double dotp0,dotr0,dotr1,doth;
	double dotp0,dotr0,dotr1,doth;
	int dimension, totalNum, maxRowNum, totalNumCOO; 
    int *row_idx, *numInRow, *I, *J; 
    double *V;
    dimension = localMatrix->dimension; 
    totalNum = localMatrix->totalNum; 
    maxRowNum = localMatrix->maxRowNum; 

	row_idx = localMatrix->row_idx; 
    numInRow = localMatrix->numInRow; 

    I = localMatrix->I; 
    J = localMatrix->J;
    V = localMatrix->V;
	
	int *col_d;
	double *V_d;
	volatile bool RODR, BLOCK, CACHE;
	int *I_COO_d, *J_COO_d;
	double *V_COO_d;
	cublasHandle_t handle;
	cublasCreate(&handle);

	double *bp_d, *pk_d, *rk_d, *vector_out_d;
	size_t size0=dimension*sizeof(int);
	size_t size1=dimension*sizeof(double);
	cudaMalloc((void **) &bp_d,size1);
	cudaMalloc((void **) &pk_d,size1);
	cudaMalloc((void **) &rk_d,size1);
	cudaMalloc((void **) &vector_out_d,size1);
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
	for (int i=0;i<dimension;i++) {
		doth=doth+vector_in[i]*vector_in[i];
	}
	struct timeval start1, end1;
	double *bp=(double *) malloc(size1);
	double *bp_g =(double *) malloc(size1);
	double *pk=(double *) malloc(size1);
	double *rk=(double *) malloc(size1);

	double *bp_dt = (double *) malloc(size1);
	double *pk_dt = (double *) malloc(size1);
	double *rk_dt = (double *) malloc(size1);
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
	cusparseOperation_t transA = cusparseOperation_t.CUSPARSE_OPERATION_NON_TRANSPOSE;

	gettimeofday(&start1, NULL);

	while (error>threshold&&iter<MAXIter){
		dotp0=0;
		dotr0=0;
		dotr1=0;
		int errorIdx = 0;
		double compareError;
		
		cudaMemset(bp_d, 0, size1);
		cusparseStatus_t smpvStatus = 
		cusparseDcsrmv(handle,
				transA,
				dimension,
				dimension,
				totalNum,
				1.0,
				CUSPARSE_MATRIX_TYPE_GENERAL,
				V_d,
				rowIdx_d,
				col_d,
				pk_d,
				0.0,
				bp_d);

		cublasDdot(handle,dimension,bp_d,1,pk_d,1,&dotp0);
		cublasDdot(handle,dimension,rk_d,1,rk_d,1,&dotr0);
			
		alphak=dotr0/dotp0;
		_alphak = -alphak;
		
		cublasDaxpy(handle,dimension,&alphak,pk_d,1,vector_out_d,1);
		cublasDaxpy(handle,dimension,&_alphak,bp_d,1,rk_d,1);
		cublasDdot(handle,dimension,rk_d,1,rk_d,1,&dotr1);
		
		gamak=dotr1/dotr0;

		cublasDscal(handle,dimension,&gamak,pk_d,1);
		cublasDaxpy(handle,dimension,&const1, rk_d, 1, pk_d, 1);
		
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

