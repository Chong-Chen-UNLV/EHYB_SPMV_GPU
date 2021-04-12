#include "kernel.h"

#define FULL_MASK 0xffffffff
#define BASE 262144 //1024*1024

#define block_size 512	
#define thread_size 512
#define block_size2 16
#define thread_size2 512
#define warpSize 32

/*kernel function for initialize*/
__global__ void kernelInitialize(const int num, float *x)
{
	int idx=blockDim.x * blockIdx.x+ threadIdx.x;
	
	for (int n=idx;n<num;n+=BASE) x[n]=0;
}

__global__ void kernelInitializeAll(const int num, float *pk, float *bp, float *x, float *zk, const float *vector_in)
{
	int idx=blockDim.x * blockIdx.x+ threadIdx.x;
	float temp;
	for (int n=idx;n<num;n+=BASE) 
	{
		temp=zk[n];
		pk[n]=temp;
		bp[n]=0;
		x[n]=0;
	}
}

__global__ void kernelInitializeR(const int num,float *rk, const float *vector_in)
{
	int idx=blockDim.x * blockIdx.x+ threadIdx.x;
	float temp;
	for (int n=idx;n<num;n+=BASE) 
	{
		temp=vector_in[n];
		rk[n]=temp;
	}
}
__global__ void vecReorderER(const int numOfRowER,
			const int* rowVecER,
			const float* yER,
			float* y){
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
			const float *valER, const float * x, float * y)
{
	int width;  
	int bias;
	int dataIdx;
	int row, col;	
	float val;
	uint32_t idx = blockDim.x*blockIdx.x+threadIdx.x;
	uint32_t warpIdx = idx>>5;
	uint32_t warpLane = threadIdx.x- ((threadIdx.x>>5)<<5);
	if(idx < numOfRowER){
		row = rowVecER[idx];
		width = widthVecER[warpIdx];//cache will work when every threads read same global address
		bias = biasVecER[warpIdx];
		float dot = 0;
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
		const int16_t* widthVecBlockELL,
		const int* biasVecBlockELL,  
		const int16_t *colBlockELL, 
		const float *valBlockELL, 
		const int numOfRowER,
		const int* biasVecER, 
		const int16_t* widthVecER, 
		const int* colER, const float* valER,
		int* warpIdxER,
		const float * x,
		float * y,
		float * yER,
		const int* partBoundary)
{
	int partIdx = blockIdx.x; 
	int xIdx = threadIdx.x;
	extern __shared__ float cachedVec[];  
	__shared__ int biasIdxBlock; 
	extern __shared__ float cachedVec[];
	//__shared__ volatile int sharedBias[blockPerPart];  
	//__shared__ volatile int sharedWidth[blockPerPart];  
	int vecStart = partBoundary[blockIdx.x];
	int vecEnd = partBoundary[blockIdx.x + 1];
	int warpLane = xIdx - ((xIdx>>5)<<5); //xIdx%32 = xIdx - (xIdx/32)*32)
	int warpIdx = (xIdx/32);
	int row = 0;
	int biasIdxWarp;
	int blockStartIdx = blockPerPart*partIdx;	
	for (int i = xIdx; i < vectorCacheSize; i += threadELL){
		cachedVec[i] = x[i + vecStart];
		if(i + vecStart + vectorCacheSize < vecEnd)
			y[i + vecStart +  vectorCacheSize] = 0;
	}
	//if(xIdx < blockPerPart){
	//	sharedBias[xIdx] = biasVecBlockELL[blockStartIdx + xIdx];	
	//	sharedWidth[xIdx] = widthVecBlockELL[blockStartIdx+ xIdx];	
	//}
	if(xIdx == 0) biasIdxBlock = warpPerBlock; 
	biasIdxWarp = warpIdx;
	__syncthreads();
	float val, dot;
	int dataIdx; 
	int col;
	int bias, width;
	while(biasIdxWarp < blockPerPart){//the thread is step with stride threadELL
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
			float dot = 0;
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
		const int16_t kernelPerPart,
		int* biasIdxBlock,
		const int16_t* widthVecBlockELL,
		const int* biasVecBlockELL,  
		const int16_t *colBlockELL, 
		const float *valBlockELL, 
		const int numOfRowER,
		const int* biasVecER, 
		const int16_t* widthVecER, 
		const int* colER, const float* valER,
		int* warpIdxER,
		const float * x,
		float * y,
		float * yER,
		const int* partBoundary)
{
	int partIdx = (blockIdx.x)/kernelPerPart; 
	int xIdx = threadIdx.x;
	extern __shared__ float cachedVec[];  
	//__shared__ int biasIdxBlock; 
	//__shared__ volatile int sharedBias[blockPerPart];  
	//__shared__ volatile int sharedWidth[blockPerPart];  
	int vecStart = partBoundary[partIdx];
	int vecEnd = partBoundary[partIdx + 1];
	int warpLane = xIdx - ((xIdx>>5)<<5); //xIdx%32 = xIdx - (xIdx/32)*32)
	int warpIdx = (xIdx/32);
	int row = 0;
	int biasIdxWarp;
	int blockStartIdx = blockPerPart*partIdx;	
	for (int i = xIdx; i < vectorCacheSize; i += threadELL){
		cachedVec[i] = x[i + vecStart];
		if((i + vecStart + vectorCacheSize < vecEnd) && (blockIdx.x)%kernelPerPart == 0)
			y[i + vecStart +  vectorCacheSize] = 0;
	}
	
	if(xIdx == 0 && (blockIdx.x)%kernelPerPart == 0) biasIdxBlock[partIdx] = warpPerBlock*kernelPerPart; 
	biasIdxWarp = warpIdx*kernelPerPart+((blockIdx.x)%kernelPerPart);
	__syncthreads();
	float val, dot;
	int dataIdx; 
	int col;
	int bias, width;
	while(biasIdxWarp < blockPerPart){//the thread is step with stride threadELL
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
			float dot = 0;
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

__global__ void kernelCachedBlockedELL_NC(
		//int16_t* biasIdxBlock,
		const int16_t* widthVecBlockELL,
		const int* biasVecBlockELL,  
		const int16_t *colBlockELL, 
		const float *valBlockELL, 
		const float * x,
		float * y,
		const int* partBoundary)
{
	int partIdx = blockIdx.x; 
	int xIdx = threadIdx.x;
	//__shared__ volatile float cachedVec[vectorCacheSize];  
	__shared__ int biasIdxBlock; 
	//__shared__ volatile int sharedBias[blockPerPart];  
	//__shared__ volatile int sharedWidth[blockPerPart];  
	int vecStart = partBoundary[blockIdx.x];
	int vecEnd = partBoundary[blockIdx.x + 1];
	int warpLane = xIdx - ((xIdx>>5)<<5); //xIdx%32 = xIdx - (xIdx/32)*32)
	int warpIdx = (xIdx/32);
	int row = 0;
	int biasIdxWarp;
	int blockStartIdx = blockPerPart*partIdx;	
	//for (int i = xIdx; i < vectorCacheSize; i += threadELL){
	//	cachedVec[i] = x[i + vecStart];
	//}
	//if(xIdx < blockPerPart){
	//	sharedBias[xIdx] = biasVecBlockELL[blockStartIdx + xIdx];	
	//	sharedWidth[xIdx] = widthVecBlockELL[blockStartIdx+ xIdx];	
	//}
	if(xIdx == 0) biasIdxBlock = warpPerBlock; 
	biasIdxWarp = warpIdx;
	__syncthreads();
	float val, dot;
	int dataIdx; 
	int col;
	int bias, width;
	for(int i = 0; i < loopInKernel; ++i){//the thread is step with stride threadELL
		dot = 0;
		row = warpLane + biasIdxWarp*warpSize + vecStart;
		if(row < vecEnd){
			bias = biasVecBlockELL[biasIdxWarp + blockStartIdx]; 
			width = widthVecBlockELL[biasIdxWarp + blockStartIdx];
			for(int n=0; n< width; ++n){
				dataIdx = bias + warpSize*n + warpLane;//however the data storage is stride with block_rowSize
				val= valBlockELL[dataIdx];
				col = colBlockELL[dataIdx];
				dot += val*x[vecStart + col];
			}
			y[row] = dot;
		}
		if(warpLane == 0)
			biasIdxWarp = atomicAdd(&biasIdxBlock, 1); 
		biasIdxWarp = __shfl_sync(FULL_MASK, biasIdxWarp, 0);
	 	__syncwarp();	
	}
}


//y=x+gamak*y
__global__ void kernelMyxpy(const int dimension, float gamak, const float *x, float *y)
{
	int idx=blockDim.x*blockIdx.x+threadIdx.x;
	int n=idx;
	while(n<dimension){
		y[n]=x[n]+gamak*y[n];
		n=n+BASE;
	}
}

extern "C"
void initialize_all(const int dimension, float *pk_d, float *bp_d, float *x, float *zk, const float *vector_in_d)
{
	kernelInitializeAll<<<block_size,thread_size>>>(dimension, pk_d, bp_d, x, zk, vector_in_d);
}

void initialize_bp(int num, float *x)
{
	kernelInitialize<<<block_size,thread_size>>>(num,x);
}

void initialize_r(int num, float *rk, float *vector_in)
{
	kernelInitializeR<<<block_size,thread_size>>>(num,rk,vector_in);
}
void myxpy(const int dimension, float gamak, const float *x, float *y)
{
	kernelMyxpy<<<block_size,thread_size>>>(dimension,gamak,x,y);
}

void initialDeviceArray(int num, float *x)
{
	kernelInitialize<<<512,512>>>(num,x);
}



void matrixVectorBlockELL(const int nParts, 
		const int16_t* widthVecBlockELL_d, 
		const int* biasVecBlockELL_d,    
		const int16_t* colBlockELL_d,
		const float* valBlockELL_d, 
		const int* partBoundary_d,
		const int numOfRowER,
		const int* rowVecER_d, 
		const int* biasVecER_d, 
		const int16_t* widthVecER_d, 
		const int* colER_d, 
		const float* valER_d,
		int* warpIdxER_d,
		const float *x_d, 
		float *y_d,
		float *yER_d)
{
 	//int maxbytes = 65536; // 64 KB
 	//int maxbytes = 73728; // 72 KB
 	//int maxbytes = 81920; // 80 KB
 	int maxbytes = 96256; // 94 KB
	cudaFuncSetAttribute(kernelCachedBlockedELL, cudaFuncAttributeMaxDynamicSharedMemorySize, maxbytes);
	kernelCachedBlockedELL<<<nParts, threadELL, sharedPerBlock>>>(
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
	if(numOfRowER/80 > 1023){
		threadSizeER = 1024;
		blockSizeER = ceil(((float) numOfRowER)/1024);
	} else {
		blockSizeER = 80;
		threadSizeER = 	ceil(((float) numOfRowER)/80);
	}
	vecReorderER<<<blockSizeER, threadSizeER>>>(numOfRowER, rowVecER_d, yER_d, y_d);
}

void matrixVectorBlockELL_small(const int nParts,	
		const int16_t kernelPerPart,
		int* biasIdxBlock_d,
		const int16_t* widthVecBlockELL_d, 
		const int* biasVecBlockELL_d,    
		const int16_t* colBlockELL_d,
		const float* valBlockELL_d, 
		const int* partBoundary_d,
		const int numOfRowER,
		const int* rowVecER_d, 
		const int* biasVecER_d, 
		const int16_t* widthVecER_d, 
		const int* colER_d, 
		const float* valER_d,
		int* warpIdxER_d,
		const float *x_d, 
		float *y_d,
		float *yER_d)
{
 	//int maxbytes = 65536; // 64 KB
 	//int maxbytes = 73728; // 72 KB
 	//int maxbytes = 81920; // 80 KB
 	//int maxbytes = 96256; // 94 KB
    //cudaFuncSetAttribute(kernelCachedBlockedELL_small, cudaFuncAttributeMaxDynamicSharedMemorySize, maxbytes);
	uint16_t kernelNum = kernelPerPart*nParts;
	kernelCachedBlockedELL_small<<<kernelNum, threadELL, sharedPerBlock>>>(
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
	if(numOfRowER/80 > 1023){
		threadSizeER = 1024;
		blockSizeER = ceil(((float) numOfRowER)/1024);
	} else {
		blockSizeER = 80;
		threadSizeER = 	ceil(((float) numOfRowER)/80);
	}
	vecReorderER<<<blockSizeER, threadSizeER>>>(numOfRowER, rowVecER_d, yER_d, y_d);
}

void matrixVectorER(const int numOfRowER, 
		const int* rowVecER_d, const int* biasVecER_d, 
		const int16_t* widthVecER_d, 
		const int* colER_d, const float* valER_d,
		const float* vectorIn_d, float* vectorOut_d)
{

	int blockSizeLocal;
	blockSizeLocal=ceil(((float) numOfRowER)/threadELL);//for data with 2 million elements, we have interval size 200
	kernelER<<<blockSizeLocal, threadELL>>>(numOfRowER, 
			rowVecER_d, 
			biasVecER_d, 
			widthVecER_d,
			colER_d, 
			valER_d, 
			vectorIn_d, 
			vectorOut_d);

}

void matrixVectorEHYB_NC(matrixEHYB* inputMatrix_d, 
		//int16_t* biasIdxBlock_d, 
		float* vectorIn_d,
		float* vectorOut_d)
{

	kernelCachedBlockedELL_NC<<<inputMatrix_d->nParts, threadELL>>>(
			//biasIdxBlock_d,
			inputMatrix_d->widthVecBlockELL,
			inputMatrix_d->biasVecBlockELL,  
			inputMatrix_d->colBlockELL, 
			inputMatrix_d->valBlockELL, 
			vectorIn_d,
			vectorOut_d,
			inputMatrix_d->partBoundary);

	
	matrixVectorER(inputMatrix_d->numOfRowER, inputMatrix_d->rowVecER, 
			inputMatrix_d->biasVecER,
			inputMatrix_d->widthVecER,
			inputMatrix_d->colER, 
			inputMatrix_d->valER,
			vectorIn_d, vectorOut_d);

}

void matrixVectorEHYB(matrixEHYB* inputMatrix_d, 
		//int16_t* biasIdxBlock_d, 
		float* vectorIn_d,
		float* vectorOut_d)
{
	matrixVectorBlockELL(inputMatrix_d->nParts, 
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
			inputMatrix_d->outER);
			
	//matrixVectorER(inputMatrix_d->numOfRowER, inputMatrix_d->rowVecER, 
	//		inputMatrix_d->biasVecER,
	//		inputMatrix_d->widthVecER,
	//		inputMatrix_d->colER, 
	//		inputMatrix_d->valER,
	//		vectorIn_d, vectorOut_d);

}

void matrixVectorEHYB_small(matrixEHYB* inputMatrix_d, 
		//const int16_t kernelPerPart,
		int* biasIdxBlock_d, 
		float* vectorIn_d,
		float* vectorOut_d)
{

	matrixVectorBlockELL_small(inputMatrix_d->nParts, 
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
			inputMatrix_d->outER);
	
	//matrixVectorER(inputMatrix_d->numOfRowER, inputMatrix_d->rowVecER, 
	//		inputMatrix_d->biasVecER,
	//		inputMatrix_d->widthVecER,
	//		inputMatrix_d->colER, 
	//		inputMatrix_d->valER,
	//		vectorIn_d, vectorOut_d);

}
