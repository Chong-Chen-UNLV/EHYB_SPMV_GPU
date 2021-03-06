#include "kernel.h"

#define FULL_MASK 0xffffffff
#define BASE 262144 //1024*1024

#define block_size 512	
#define thread_size 512
#define block_size2 16
#define thread_size2 512
#define warpSize 32

/*kernel function for initialize*/
__global__ void kernelInitialize(const int num, double *x)
{
	int idx=blockDim.x * blockIdx.x+ threadIdx.x;
	
	for (int n=idx;n<num;n+=BASE) x[n]=0;
}

__global__ void kernelInitializeAll(const int num, double *pk, double *bp, double *x, double *zk, const double *vector_in)
{
	int idx=blockDim.x * blockIdx.x+ threadIdx.x;
	double temp;
	for (int n=idx;n<num;n+=BASE) 
	{
		temp=zk[n];
		pk[n]=temp;
		bp[n]=0;
		x[n]=0;
	}
}

__global__ void kernelInitializeR(const int num,double *rk, const double *vector_in)
{
	int idx=blockDim.x * blockIdx.x+ threadIdx.x;
	double temp;
	for (int n=idx;n<num;n+=BASE) 
	{
		temp=vector_in[n];
		rk[n]=temp;
	}
}
__global__ void longRowKernel(const int* vecBoundary, 
		const int* row, 
		const int* col,
		const double* val,
		const double* x,
		double* y)
{
	int partIdx = blockIdx.x;
	int xIdx = threadIdx.x;
	int vecStart = vecBoundary[partIdx];
	int vecEnd = vecBoundary[partIdx + 1];
	int rowIdx = row[partIdx];
	uint32_t warpLane = threadIdx.x- ((threadIdx.x>>5)<<5);
	double result = 0;	
	for(int colIdx = vecStart + xIdx; colIdx < vecEnd; colIdx+= threadLongVec){
		result += val[colIdx]*x[col[colIdx]];	
	}
	__syncwarp();
	for (int offset = 16; offset > 0; offset /= 2)
		result += __shfl_down_sync(FULL_MASK, result, offset);
	__syncwarp();
	if(warpLane == 0){
		atomicAdd(&y[rowIdx], result); 
	}
}

__global__ void vecReorderER(const int numOfRowER,
			const int* rowVecER,
			const double* yER,
			double* y){
	uint32_t idx = blockDim.x*blockIdx.x+threadIdx.x;
	if(idx < numOfRowER){
		y[rowVecER[idx]]+=yER[idx];
	}
}
//first version implementation,
//concern about performance loss from inbalance between blocks 
__global__ void kernelER(const int numOfRowER,
			const int* rowVecER,
			const int* biasVecER,  
			const int16_t* widthVecER, 
			const int* colER, 
			const double *valER, const double * x, double * y)
{
	int width;  
	int bias;
	int dataIdx;
	int row, col;	
	double val;
	uint32_t idx = blockDim.x*blockIdx.x+threadIdx.x;
	uint32_t warpIdx = idx>>5;
	uint32_t warpLane = threadIdx.x- ((threadIdx.x>>5)<<5);
	if(idx < numOfRowER){
		row = rowVecER[idx];
		width = widthVecER[warpIdx];//cache will work when every threads read same global address
		bias = biasVecER[warpIdx];
		double dot = 0;
		for(int n=0; n < width; ++n){
			dataIdx = bias + warpLane + warpSize*n ;
			col=colER[dataIdx];
			val=valER[dataIdx];
			dot += val* x[col];
		}
		y[row]+=dot;
	}
}

__global__ void kernelCachedBlockedELL(
		const int16_t vectorCacheSize,
		const int16_t* widthVecBlockELL,
		const int* biasVecBlockELL,  
		const int16_t *colBlockELL, 
		const double *valBlockELL, 
		const int numOfRowER,
		const int* biasVecER, 
		const int16_t* widthVecER, 
		const int* colER, const double* valER,
		int* warpIdxER,
		const double * x,
		double * y,
		double * yER,
		const int* partBoundary)
{
	int partIdx = blockIdx.x; 
	int xIdx = threadIdx.x;
	__shared__ int biasIdxBlock; 
	extern __shared__ double cachedVec[];
	int vecStart = partBoundary[blockIdx.x];
	int vecEnd = partBoundary[blockIdx.x + 1];
	int warpLane = xIdx - ((xIdx>>5)<<5); //xIdx%32 = xIdx - (xIdx/32)*32)
	int warpIdx = (xIdx/32);
	int row = 0;
	int biasIdxWarp;
	int blockStartIdx = (vectorCacheSize/32)*partIdx;	
	for (int i = xIdx; i < vectorCacheSize; i += threadELL){
		cachedVec[i] = x[i + vecStart];
		if(i + vecStart + vectorCacheSize < vecEnd)
			y[i + vecStart +  vectorCacheSize] = 0;
	}
	if(xIdx == 0) biasIdxBlock = warpPerBlock; 

	biasIdxWarp = warpIdx;
	__syncthreads();
	double val, dot;
	int dataIdx; 
	int col;
	int bias, width;
	while(biasIdxWarp < vectorCacheSize/32){//the thread is step with stride threadELL
		dot = 0;
		row = warpLane + biasIdxWarp*warpSize + vecStart;
		if(row < vecEnd){
			bias = biasVecBlockELL[biasIdxWarp + blockStartIdx]; 
			width = widthVecBlockELL[biasIdxWarp + blockStartIdx];
			for(int n=0; n< width; ++n){
				dataIdx = bias + warpSize*n + warpLane;//however the data storage is stride with block_rowSize
				val= valBlockELL[dataIdx];
				col = colBlockELL[dataIdx];
				dot += val*cachedVec[col];
			}
			y[row] = dot;
		}
		if(warpLane == 0)
			biasIdxWarp = atomicAdd(&biasIdxBlock, 1); 
		biasIdxWarp = __shfl_sync(FULL_MASK, biasIdxWarp, 0);
	 	__syncwarp();	
	}
	biasIdxWarp = 0;
	__syncwarp();
	if(warpLane == 0)
		biasIdxWarp = atomicAdd(warpIdxER, 1); 
	__syncwarp();	
	biasIdxWarp = __shfl_sync(FULL_MASK, biasIdxWarp, 0);

	while(biasIdxWarp < (numOfRowER/32 +1)){
		row = biasIdxWarp*32 + warpLane;
		if(row < numOfRowER){
			width = widthVecER[biasIdxWarp];//cache will work when every threads read same global address
			bias = biasVecER[biasIdxWarp];
			double dot = 0;
			for(int n=0; n < width; ++n){
				dataIdx = bias + warpLane + warpSize*n ;
				col=colER[dataIdx];
				val=valER[dataIdx];
				dot += val* x[col];
			}
			yER[row]=dot;
		}
		if(warpLane == 0)
			biasIdxWarp = atomicAdd(warpIdxER, 1); 
		__syncwarp();	
		biasIdxWarp = __shfl_sync(FULL_MASK, biasIdxWarp, 0);
	}
}

__global__ void kernelCachedBlockedELL_small(
		const int16_t vectorCacheSize,
		const int16_t kernelPerPart,
		int* biasIdxBlock,
		const int16_t* widthVecBlockELL,
		const int* biasVecBlockELL,  
		const int16_t *colBlockELL, 
		const double *valBlockELL, 
		const int numOfRowER,
		const int* biasVecER, 
		const int16_t* widthVecER, 
		const int* colER, const double* valER,
		int* warpIdxER,
		const double * x,
		double * y,
		double * yER,
		const int* partBoundary)
{
	int partIdx = (blockIdx.x)/kernelPerPart; 
	int xIdx = threadIdx.x;
	extern __shared__ double cachedVec[];  
	//__shared__ int biasIdxBlock; 
	int vecStart = partBoundary[partIdx];
	int vecEnd = partBoundary[partIdx + 1];
	int warpLane = xIdx - ((xIdx>>5)<<5); //xIdx%32 = xIdx - (xIdx/32)*32)
	int warpIdx = (xIdx/warpSize);
	int row = 0;
	int biasIdxWarp;
	int blockStartIdx = (vectorCacheSize/warpSize)*partIdx;	
	for (int i = xIdx; i < vectorCacheSize; i += threadELL){
		cachedVec[i] = x[i + vecStart];
		if((i + vecStart + vectorCacheSize < vecEnd) && (blockIdx.x)%kernelPerPart == 0)
			y[i + vecStart +  vectorCacheSize] = 0;
	}
	
	if(xIdx == 0 && (blockIdx.x)%kernelPerPart == 0) biasIdxBlock[partIdx] = warpPerBlock*kernelPerPart; 
	biasIdxWarp = warpIdx*kernelPerPart+((blockIdx.x)%kernelPerPart);
	__syncthreads();
	double val, dot;
	int dataIdx; 
	int col;
	int bias, width;
	while(biasIdxWarp < vectorCacheSize/warpSize){//the thread is step with stride threadELL
		dot = 0;
		row = warpLane + biasIdxWarp*warpSize + vecStart;
		if(row < vecEnd){
			bias = biasVecBlockELL[biasIdxWarp + blockStartIdx]; 
			width = widthVecBlockELL[biasIdxWarp + blockStartIdx];
			for(int n=0; n< width; ++n){
				dataIdx = bias + warpSize*n + warpLane;//however the data storage is stride with block_rowSize
				val= valBlockELL[dataIdx];
				col = colBlockELL[dataIdx];
				dot += val*cachedVec[col];
			}
			y[row] = dot;
		}
		if(warpLane == 0)
			biasIdxWarp = atomicAdd(&(biasIdxBlock[partIdx]), 1); 
		__syncwarp();	
		biasIdxWarp = __shfl_sync(FULL_MASK, biasIdxWarp, 0);
	}
	biasIdxWarp = 0;
	__syncwarp();
	if(warpLane == 0)
		biasIdxWarp = atomicAdd(warpIdxER, 1); 
	__syncwarp();	
	biasIdxWarp = __shfl_sync(FULL_MASK, biasIdxWarp, 0);

	while(biasIdxWarp < (numOfRowER/32 +1)){
		row = biasIdxWarp*32 + warpLane;
		if(row < numOfRowER){
			width = widthVecER[biasIdxWarp];//cache will work when every threads read same global address
			bias = biasVecER[biasIdxWarp];
			double dot = 0;
			for(int n=0; n < width; ++n){
				dataIdx = bias + warpLane + warpSize*n ;
				col=colER[dataIdx];
				val=valER[dataIdx];
				dot += val* x[col];
			}
			yER[row]=dot;
		}
		if(warpLane == 0)
			biasIdxWarp = atomicAdd(warpIdxER, 1); 
		__syncwarp();	
		biasIdxWarp = __shfl_sync(FULL_MASK, biasIdxWarp, 0);
	}
}


//y=x+gamak*y
__global__ void kernelMyxpy(const int dimension, double gamak, const double *x, double *y)
{
	int idx=blockDim.x*blockIdx.x+threadIdx.x;
	int n=idx;
	while(n<dimension){
		y[n]=x[n]+gamak*y[n];
		n=n+BASE;
	}
}

extern "C"
void initialize_all(const int dimension, double *pk_d, double *bp_d, double *x, double *zk, const double *vector_in_d)
{
	kernelInitializeAll<<<block_size,thread_size>>>(dimension, pk_d, bp_d, x, zk, vector_in_d);
}

void initialize_bp(int num, double *x)
{
	kernelInitialize<<<block_size,thread_size>>>(num,x);
}

void initialize_r(int num, double *rk, double *vector_in)
{
	kernelInitializeR<<<block_size,thread_size>>>(num,rk,vector_in);
}
void myxpy(const int dimension, double gamak, const double *x, double *y)
{
	kernelMyxpy<<<block_size,thread_size>>>(dimension,gamak,x,y);
}

void initialDeviceArray(int num, double *x)
{
	kernelInitialize<<<512,512>>>(num,x);
}


void matrixVectorBlockELL(const int nParts, 
		const int16_t vectorCacheSize,
		const int16_t* widthVecBlockELL_d, 
		const int* biasVecBlockELL_d,    
		const int16_t* colBlockELL_d,
		const double* valBlockELL_d, 
		const int* partBoundary_d,
		const int numOfRowER,
		const int* rowVecER_d, 
		const int* biasVecER_d, 
		const int16_t* widthVecER_d, 
		const int* colER_d, 
		const double* valER_d,
		int* warpIdxER_d,
		const double *x_d, 
		double *y_d,
		double *yER_d,
		const int nLongVec,
		const int* longVecBoundary,
		const int* longVecRow,
		const int* longVecCol,
		double* longVecVal)
{
 	//int maxbytes = 65536; // 64 KB
 	//int maxbytes = 73728; // 72 KB
 	//int maxbytes = 81920; // 80 KB
 	int maxbytes = 96256; // 94 KB
	cudaFuncSetAttribute(kernelCachedBlockedELL, cudaFuncAttributeMaxDynamicSharedMemorySize, maxbytes);
	kernelCachedBlockedELL<<<nParts, threadELL, maxbytes>>>(
			vectorCacheSize,
			widthVecBlockELL_d,
			biasVecBlockELL_d,  
			colBlockELL_d, valBlockELL_d, 
			numOfRowER,
			biasVecER_d, 
			widthVecER_d, 
			colER_d, valER_d,
			warpIdxER_d,
			x_d,
			y_d,
			yER_d,
			partBoundary_d);
	int threadSizeER, blockSizeER;
	if(numOfRowER/smSize> 1023){
		threadSizeER = 1024;
		blockSizeER = ceil(((double) numOfRowER)/1024);
	} else {
		blockSizeER = smSize;
		threadSizeER = 	ceil(((double) numOfRowER)/smSize);
	}
	vecReorderER<<<blockSizeER, threadSizeER>>>(numOfRowER, rowVecER_d, yER_d, y_d);

	if(nLongVec > 0){
		longRowKernel<<<nLongVec, threadLongVec>>>(longVecBoundary,longVecRow,
				longVecCol, longVecVal, x_d, y_d);
	}
}

void matrixVectorBlockELL_small(const int nParts,	
		int16_t vectorCacheSize,
		const int16_t kernelPerPart,
		int* biasIdxBlock_d,
		const int16_t* widthVecBlockELL_d, 
		const int* biasVecBlockELL_d,    
		const int16_t* colBlockELL_d,
		const double* valBlockELL_d, 
		const int* partBoundary_d,
		const int numOfRowER,
		const int* rowVecER_d, 
		const int* biasVecER_d, 
		const int16_t* widthVecER_d, 
		const int* colER_d, 
		const double* valER_d,
		int* warpIdxER_d,
		const double *x_d, 
		double *y_d,
		double *yER_d,
		const int nLongVec,
		const int* longVecBoundary,
		const int* longVecRow,
		const int* longVecCol,
		double* longVecVal)
{
 	//int maxbytes = 65536; // 64 KB
 	//int maxbytes = 73728; // 72 KB
 	//int maxbytes = 81920; // 80 KB
 	int maxbytes = 96256; // 94 KB
    cudaFuncSetAttribute(kernelCachedBlockedELL_small, cudaFuncAttributeMaxDynamicSharedMemorySize, maxbytes);
	uint16_t kernelNum = kernelPerPart*nParts;
	kernelCachedBlockedELL_small<<<kernelNum, threadELL, maxbytes>>>(
			vectorCacheSize,
			kernelPerPart,
			biasIdxBlock_d,
			widthVecBlockELL_d,
			biasVecBlockELL_d,  
			colBlockELL_d, valBlockELL_d, 
			numOfRowER,
			biasVecER_d, 
			widthVecER_d, 
			colER_d, valER_d,
			warpIdxER_d,
			x_d,
			y_d,
			yER_d,
			partBoundary_d);
	int threadSizeER, blockSizeER;
	if(numOfRowER/smSize> 1023){
		threadSizeER = 1024;
		blockSizeER = ceil(((double) numOfRowER)/1024);
	} else {
		blockSizeER = smSize;
		threadSizeER = 	ceil(((double) numOfRowER)/smSize);
	}
	vecReorderER<<<blockSizeER, threadSizeER>>>(numOfRowER, rowVecER_d, yER_d, y_d);
	if(nLongVec > 0){
		longRowKernel<<<nLongVec, threadLongVec>>>(longVecBoundary,longVecRow,
				longVecCol, longVecVal, x_d, y_d);
	}
}

void matrixVectorER(const int numOfRowER, 
		const int* rowVecER_d, const int* biasVecER_d, 
		const int16_t* widthVecER_d, 
		const int* colER_d, const double* valER_d,
		const double* vectorIn_d, double* vectorOut_d)
{

	int blockSizeLocal;
	blockSizeLocal=ceil(((double) numOfRowER)/threadELL);//for data with 2 million elements, we have interval size 200
	kernelER<<<blockSizeLocal, threadELL>>>(numOfRowER, 
			rowVecER_d, 
			biasVecER_d, 
			widthVecER_d,
			colER_d, 
			valER_d, 
			vectorIn_d, 
			vectorOut_d);

}

//void matrixVectorEHYB_NC(matrixEHYB* inputMatrix_d, 
//		//int16_t* biasIdxBlock_d, 
//		double* vectorIn_d,
//		double* vectorOut_d)
//{
//
//	kernelCachedBlockedELL_NC<<<inputMatrix_d->nParts, threadELL>>>(
//			//biasIdxBlock_d,
//			inputMatrix_d->widthVecBlockELL,
//			inputMatrix_d->biasVecBlockELL,  
//			inputMatrix_d->colBlockELL, 
//			inputMatrix_d->valBlockELL, 
//			vectorIn_d,
//			vectorOut_d,
//			inputMatrix_d->partBoundary);
//
//	
//	matrixVectorER(inputMatrix_d->numOfRowER, inputMatrix_d->rowVecER, 
//			inputMatrix_d->biasVecER,
//			inputMatrix_d->widthVecER,
//			inputMatrix_d->colER, 
//			inputMatrix_d->valER,
//			vectorIn_d, vectorOut_d);
//
//}

void matrixVectorEHYB(matrixEHYB* inputMatrix_d, 
		//int16_t* biasIdxBlock_d, 
		double* vectorIn_d,
		double* vectorOut_d)
{
	matrixVectorBlockELL(inputMatrix_d->nParts, 
			inputMatrix_d->vectorCacheSize,
			inputMatrix_d->widthVecBlockELL,
			inputMatrix_d->biasVecBlockELL,  
			inputMatrix_d->colBlockELL, 
			inputMatrix_d->valBlockELL, 
			inputMatrix_d->partBoundary,
			inputMatrix_d->numOfRowER, 
			inputMatrix_d->rowVecER, 
			inputMatrix_d->biasVecER,
			inputMatrix_d->widthVecER,
			inputMatrix_d->colER, 
			inputMatrix_d->valER,
			inputMatrix_d->warpIdxER_d,
			vectorIn_d,
			vectorOut_d,
			inputMatrix_d->outER,
			inputMatrix_d->nLongVec,
			inputMatrix_d->longVecBoundary,
			inputMatrix_d->longVecRow,
			inputMatrix_d->longVecCol,
			inputMatrix_d->longVecVal);

}

void matrixVectorEHYB_small(matrixEHYB* inputMatrix_d, 
		//const int16_t kernelPerPart,
		int* biasIdxBlock_d, 
		double* vectorIn_d,
		double* vectorOut_d)
{

	matrixVectorBlockELL_small(inputMatrix_d->nParts, 
			inputMatrix_d->vectorCacheSize,
			inputMatrix_d->kernelPerPart,
			biasIdxBlock_d,
			inputMatrix_d->widthVecBlockELL,
			inputMatrix_d->biasVecBlockELL,  
			inputMatrix_d->colBlockELL, 
			inputMatrix_d->valBlockELL, 
			inputMatrix_d->partBoundary,
			inputMatrix_d->numOfRowER, 
			inputMatrix_d->rowVecER, 
			inputMatrix_d->biasVecER,
			inputMatrix_d->widthVecER,
			inputMatrix_d->colER, 
			inputMatrix_d->valER,
			inputMatrix_d->warpIdxER_d,
			vectorIn_d,
			vectorOut_d,
			inputMatrix_d->outER,
			inputMatrix_d->nLongVec,
			inputMatrix_d->longVecBoundary,
			inputMatrix_d->longVecRow,
			inputMatrix_d->longVecCol,
			inputMatrix_d->longVecVal);
	
}
