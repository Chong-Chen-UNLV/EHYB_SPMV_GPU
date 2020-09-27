#include "kernel.h"

#define FULL_MASK 0xffffffff
#define BASE 262144 //1024*1024

#define block_size 512	
#define thread_size 512
#define block_size2 16
#define thread_size2 512
#define WARP_SIZE 32

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

//first version implementation,
//concern about performance loss from inbalance between blocks 
__global__ void kernelER(const int numOfRowER,
			const int* widthVecER_d, 
			const int* biasVecER_d,  
			const int* rowVecER_d;
			const int* colER_d, 
			const double *valER_d, const double * x, double * y)
{
	int width;  
	int bias;
	int dataIdx;
	int row, col;	
	double val;
	uint32_t idx = blockDim.x*blockIdx.x+threadIdx.x;
	uint32_t warpIdx = idx>>5;
	row = rowVecER_d[idx];
	if(row < numOfRowER){
		width = widthVecER_d[warpIdx];//cache will work when every threads read same global address
		bias = biasVecER_d[warpIdx];
		double dot = 0;
		for(int n=0; n < width; ++n){
			dataIdx = bias + warpSize*n + threadIdx.x;
			col=colER_d[dataIdx];
			val=valER_d[dataIdx];
			dot += val* x[col];
		}
		y[row]+=dot;
	}
}

__global__ void kernelCachedBlockedELL_test(const int* widthVecBlockELL_d,
		const int* vecBias,  
		const int *indices, const double *data, const double * x,
		double * y,
		const int* partBoundary,
		const int testPoint)
{
	int partIdx = blockIdx.x; 
	int xIdx = threadIdx.x;
	__shared__ volatile double cachedVec[vectorCacheSize];  
	__shared__ volatile int sharedBias[warpPerBlock];  
	__shared__ volatile int sharedWidth[warpPerBlock];  
	int vecStart = partBoundary[blockIdx.x];
	int vecEnd = partBoundary[blockIdx.x + 1];
	int warpLane = xIdx - (xIdx<<5)>>5; //xIdx%32 = xIdx - (xIdx/32)*32)
	int row = 0;

	for (int i = xIdx; i < vectorCacheSize; i += threadELL){
		cachedVec[i] = x[i + vecStart];
	}
	if(xIdx < blockPerPart){
		sharedBias[xIdx] = vecBias[blockPerPart*partIdx + xIdx];	
		sharedWidth[xIdx] = vecWidth[blockPerPart*partIdx + xIdx];	
	}
	double val, dot;
	int blockId, dataIdx; 
	int col;
	int bias, width;

	#pragma unroll
	for(int i = 0; i < memPerThread/elementSize; ++i){//the thread is step with stride ELL_threadSize
		dot = 0;
		//each iteration go through (1024/warpSize)=32 blocks in blockELL format, which is i >> 5
		//the warpIdx is xIdx<<4
		bias = sharedBias[i>>5 + xIdx<<5]; 
		width = sharedWidth[i>>5 + xIdx<<5];
		row = i*threadELL + vecStart + xIdx;
		if(row < vec_end){
			for(int n=0; n< width; ++n){
				dataIdx = bias + warpSize*n + warpLane;//however the data storage is stride with block_rowSize
				val=data[dataIdx];
				col = indices[dataIdx]- vecStart;
				//if(col < 0 || col >= vector_cache_size) 
				//	dot += val*cached_vec[col];
				//else
				if(row == testPoint)
					dot += val*cached_vec[col] - 1 + 0.999;
				else
					dot += val*cachedVec[col];
			}
			//if(row == testPoint)
			//	y[row] = dot+0.01;
			//else 
			y[row] = dot;
		}
	}		
}

__global__ void kernelCachedBlockedELL(const int* widthVecBlockELL_d,
		const int* vecBias,  
		const int *indices, const double *data, const double * x,
		double * y,
		const int* partBoundary)
{
	int partIdx = blockIdx.x; 
	int xIdx = threadIdx.x;
	__shared__ volatile double cachedVec[vectorCacheSize];  
	__shared__ volatile int sharedBias[warpPerBlock];  
	__shared__ volatile int sharedWidth[warpPerBlock];  
	int vecStart = partBoundary[blockIdx.x];
	int vecEnd = partBoundary[blockIdx.x + 1];
	int warpLane = xIdx - (xIdx<<5)>>5; //xIdx%32 = xIdx - (xIdx/32)*32)
	int row = 0;

	for (int i = xIdx; i < vectorCacheSize; i += threadELL){
		cachedVec[i] = x[i + vecStart];
	}
	if(xIdx < blockPerPart){
		sharedBias[xIdx] = vecBias[blockPerPart*partIdx + xIdx];	
		sharedWidth[xIdx] = vecWidth[blockPerPart*partIdx + xIdx];	
	}
	double val, dot;
	int blockId, dataIdx; 
	int col;
	int bias, width;

	#pragma unroll
	for(int i = 0; i < memPerThread/elementSize; ++i){//the thread is step with stride ELL_threadSize
		dot = 0;
		//each iteration go through (1024/warpSize)=32 blocks in blockELL format, which is i >> 5
		//the warpIdx is xIdx<<4
		bias = sharedBias[i>>5 + xIdx<<5]; 
		width = sharedWidth[i>>5 + xIdx<<5];
		row = i*threadELL + vecStart + xIdx;
		if(row < vec_end){
			for(int n=0; n< width; ++n){
				dataIdx = bias + warpSize*n + warpLane;//however the data storage is stride with block_rowSize
				val=data[dataIdx];
				col = indices[dataIdx]- vecStart;
				//if(col < 0 || col >= vector_cache_size) 
				//	dot += val*cached_vec[col];
				//else
				//if(row == testPoint)
				//	dot += val*cached_vec[col] - 1 + 0.999;
				//else
				dot += val*cachedVec[col];
			}
			//if(row == testPoint)
			//	y[row] = dot+0.01;
			//else 
			y[row] = dot;
		}
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



void matrixVectorBlockELL(const int nParts, const int testPoint, 
		const int* widthVecBlockELL_d, 
		const int* biasVecBLockELL_d,    
		const int* colBlockELL_d,
		const double* valBlockELL_d, 
		const int* partBoundary_d,
		const double *x_d, double *y_d)
{


	
	//printf("ELL_blocks is %d\n", ELL_blocks);
	if(!CACHE){	
		printf("only examine cached implementation now\n");
		exit(0);	
		//gpuErrchk( cudaPeekAtLastError() );
	} else {
		if(testPoint >= 0){
			kernelCachedBlockedELL_test<<<nParts, threadELL>>>(widthVecBlockELL_d,
					biasVecBLockELL_d,  
					colBlockELL_d, valBlockELL_d, 
					x_d,
					y_d,
					partBoundary_d,
					testPoint);
		} else {
			kernelCachedBlockedELL<<<nParts, threadELL>>>(widthVecBlockELL_d,
					biasVecBLockELL_d,  
					colBlockELL_d, valBlockELL_d, 
					x_d,
					y_d,
					partBoundary_d);
		}

	}
}

void matrixVectorER(const int numOfRowER, 
		const int* rowVecER_d, const int* biasVecER_d, 
		const int* widthVecER_d, 
		const int* colER_d, const double* valER_d,
		const double* vectorIn_d, const double* vectorOut_d)
{

	int blockSizeLocal;
	blockSizeLocal=ceil(((double) numOfRowER)/threadSizeELL);//for data with 2 million elements, we have interval size 200
	kernelER<<<blockSizeLocal, threadSizeER>>>(numOfRowER, 
			rowVecER_d, 
			biasVecER_d, 
			widthVecER_d,
			colER_d, 
			valER_d, 
			vectorIn_d, 
			vectorOut_d);

}

void matrixVectorEHYB(matrixEHYB* inputMatrix_d, double* vectorIn_d,
		double* vectorOut_d, const int testPoint)
{

	matrixVectorBlockELL(inputMatrix_d->nParts, 
			testPoint,
			inputMatrix_d->widthVecBlockELL,
			inputMatrix_d->biasVecBLockELL,  
			inputMatrix_d->colBlockELL, 
			inputMatrix_d->valBlockELL, 
			inputMatrix_d->partBoundary,
			vectorIn_d,
			vectorOut_d);
	
	matrixVectorER(numOfRowER, inputMatrix_d->rowVecER, 
			inputMatrix_d->widthVecER,
			inputMatrix_d->biasVecER,
			inputMatrix_d->colER, 
			inputMatrix_d->valER,
			vectorIn_d, vectorOut_d);

}
