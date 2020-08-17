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

//for ELL format matrix, output y=data*x
__global__ void ELL_kernel(const int num_rows, const int cal_rows, const int num_cols_per_row,
			const int *indices, const double *data, const double * x, double * y) 
{
	int row= blockDim.x*blockIdx.x+threadIdx.x;
	if (row<cal_rows){
		double dot =0;
		
		for (int n=0; n< num_cols_per_row; n++){
			int col=indices[num_rows * n + row];
			double val=data[num_rows*n+row];

			if (val != 0)
				dot += val* x[col];
		}
			y[row]=dot;
	}
}


//first version implementation,
//concern about performance loss from inbalance between blocks 
__global__ void ER_kernel(const int num_rows,
			const int* num_cols_per_row_vec, 
			const int* block_data_bias_vec,  
			const int* row_vec;
			const int* indices, 
			const double *data, const double * x, double * y)
{
	
	int block_idx = blockIdx.x; 
	int thread_idx = threadIdx.x; 
	int num_cols_per_row;  
	int block_data_bias;
	int data_idx;
	int row, col;	
	double val;
	uint32_t idx = blockDim.x*blockIdx.x+threadIdx.x;
	uint32_t warpIdx = idx>>4;
	row = row_vec[idx];
	if(row < num_rows){
		num_cols_per_row = num_cols_per_row_vec[warpIdx];//cache will work when every threads read same global address
		block_data_bias = block_data_bias_vec[warpIdx];
		double dot = 0;
		for(int n=0; n< num_cols_per_row; n++){
			data_idx = block_data_bias + warp_size*n + thread_idx;
			col=indices[data_idx];
			val=data[data_idx];
			dot += val* x[col];
		}
		y[row]+=dot;
	}
}

__global__ void ELL_kernel_block(const int num_rows,
			const int* num_cols_per_row_vec, 
			const int* block_data_bias_vec,  
			const int *indices, 
			const double *data, const double * x, double * y,
			bool tex = false) 
{
	
	int block_idx = blockIdx.x; 
	int thread_idx = threadIdx.x; 
	int num_cols_per_row;  
	int block_data_bias;
	int data_idx;
	int col;	
	double val;
	int row= blockDim.x*blockIdx.x+threadIdx.x;
	if(row < num_rows){
		num_cols_per_row = num_cols_per_row_vec[block_idx];//cache will work when every threads read same global address
		block_data_bias = block_data_bias_vec[block_idx];
		double dot =0;
		for (int n=0; n< num_cols_per_row; n++){
			data_idx = block_data_bias + ELL_threadSize*n + thread_idx;
			col=indices[data_idx];
			val=data[data_idx];

			if (val != 0){
				dot += val* x[col];
			}
		}
		y[row]+=dot;
	}
}

__global__ void ELL_kernel_rodr(const int* num_cols_per_row_vec,
		const int* block_data_bias_vec,  
		const int *indices, const double *data, const double * x,
		double * y,
		const int* part_boundary,
		const bool tex)
{
	int part_idx = blockIdx.x; 
	int x_idx = threadIdx.x;
	int vec_start = part_boundary[blockIdx.x];
	int vec_end = part_boundary[blockIdx.x + 1];
	int row = 0;
	double val, dot;
	int block_idx, data_idx, col;
	int block_rowSize, block_data_bias, num_cols_per_row;
	int block_base = part_idx * block_per_part;
		
	int endBlockRow;
	if(vec_end >= vec_start+block_per_part*ELL_threadSize)
		endBlockRow = ELL_threadSize; 
	else 
		endBlockRow = vec_end - (vec_start+(block_per_part - 1)*ELL_threadSize);

	for(int i = 0; i < block_per_part; ++i){//the thread is step with stride ELL_threadSize
		dot =0;
		block_rowSize = ELL_threadSize;
		block_idx = block_base + i;
		if(i == (block_per_part -1)){
			block_rowSize = endBlockRow;
		}
		block_data_bias = block_data_bias_vec[block_idx];
		num_cols_per_row = num_cols_per_row_vec[block_idx];//cache will work if it is shared by all threads
		row = i*ELL_threadSize + vec_start + x_idx;
		if(row < vec_end){
			for(int n=0; n< num_cols_per_row; ++n){
				data_idx = block_data_bias + block_rowSize*n + x_idx;//however the data storage is stride with block_rowSize
				col=indices[data_idx];
				val=data[data_idx];
				dot += val* x[col];

			}
			y[row] = dot;
		}
		block_idx += 1;
	}		
}

__global__ void ELL_kernel_rodr_test(const int* num_cols_per_row_vec,
		const int* block_data_bias_vec,  
		const int *indices, const double *data, const double * x,
		double * y,
		const int* part_boundary,
		const int testPoint)
{
	int part_idx = blockIdx.x; 
	int x_idx = threadIdx.x;
	int vec_start = part_boundary[blockIdx.x];
	int vec_end = part_boundary[blockIdx.x + 1];
	int row = 0;
	double val, dot;
	int block_idx, data_idx, col;
	int block_rowSize, block_data_bias, num_cols_per_row;
	int block_base = part_idx * block_per_part;
		
	int endBlockRow;
	if(vec_end >= vec_start+block_per_part*ELL_threadSize)
		endBlockRow = ELL_threadSize; 
	else 
		endBlockRow = vec_end - (vec_start+(block_per_part - 1)*ELL_threadSize);


	for(int i = 0; i < block_per_part; ++i){//the thread is step with stride ELL_threadSize
		dot =0;
		block_rowSize = ELL_threadSize;
		block_idx = block_base + i;
		if(i == (block_per_part -1)){
			block_rowSize = endBlockRow;
		}
		block_data_bias = block_data_bias_vec[block_idx];
		num_cols_per_row = num_cols_per_row_vec[block_idx];
		row = i*ELL_threadSize + vec_start + x_idx;
		if(row < vec_end){
			for(int n=0; n< num_cols_per_row; ++n){
				data_idx = block_data_bias + block_rowSize*n + x_idx;//however the data storage is stride with block_rowSize
				col=indices[data_idx];
				val=data[data_idx];
				if(val != 0){
					//if(row == testPoint)
					//	dot = dot + val* x[col];
					//else
						dot += val* x[col];
				}
			}
			//if(row == testPoint)
			//	y[row] = dot + 1 - 0.999;
			//else
				y[row] = dot;
		}
		block_idx += 1;
	}		
}

__global__ void ELL_cached_kernel_rodr(const int* vecWidth,
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
	int row = 0;

	for (int i = xIdx; i < vectorCacheSize; i += EllThreadSize){
		cachedVec[i] = x[i + vecStart];
	}
	if(xIdx < blockPerPart){
		sharedBias[xIdx] = vecBias[blockPerPart*partIdx+xIdx];	
		sharedWidth[xIdx] = vecWidth[blockPerPart*partIdx+xIdx];	
	}
	double val, dot;
	int block_idx, data_idx; 
	int col;
	int block_rowSize, block_data_bias, num_cols_per_row;
	int block_base = part_idx * block_per_part;

	for(int i = 0; i < block_per_part; ++i){//the thread is step with stride ELL_threadSize
		dot =0;
		block_rowSize = ELL_threadSize;
		block_idx = block_base + i;
		if(i == (block_per_part -1)){
			block_rowSize = endBlockRow;
		}
		block_data_bias = block_data_bias_vec[block_idx];
		num_cols_per_row = num_cols_per_row_vec[block_idx];
		row = i*ELL_threadSize + vec_start + x_idx;
		if(row < vec_end){
			for(int n=0; n< num_cols_per_row; ++n){
				data_idx = block_data_bias + block_rowSize*n + x_idx;//however the data storage is stride with block_rowSize
				val=data[data_idx];
				col = indices[data_idx]- vec_start;
				//if(col < 0 || col >= vector_cache_size) 
				//	dot += val*cached_vec[col];
				//else
				//if(row == testPoint)
				//	dot += val*cached_vec[col] - 1 + 0.999;
				//else
					dot += val*cached_vec[col];
			}
			//if(row == testPoint)
			//	y[row] = dot+0.01;
			//else 
				y[row] = dot;
		}
		block_idx += 1;
	}		
}


__global__ void COO_shared(const int num_nozeros, const int interval_size,
				const int *I, const int *J, const double *V,
				const double *x, double *y)
{
	__shared__ volatile int rows[48*thread_size/WARP_SIZE];  //why using 48? because we need 16 additional junk elements
	__shared__ volatile double vals[thread_size];

	int thread_id = blockDim.x*blockIdx.x + threadIdx.x;
	int thread_lane= threadIdx.x & (WARP_SIZE-1); //great idea! think about it
	int warp_id = thread_id / WARP_SIZE;

	int interval_begin=warp_id*interval_size;
	int interval_end =min(interval_begin+interval_size, num_nozeros);
	/*how about the interval is not the multiple of warp_size?*/
	//int iteration_end=((interval_end)/WARP_SIZE)*WARP_SIZE;

	int idx=16*(threadIdx.x/32+1) + threadIdx.x;//every warp has 16 "junk" rows elements

	rows[idx-16]=-1;
	
	int n;
	n=interval_begin+thread_lane;
	while (n< interval_end)
	{
		int row =I[n];
		//double val=V[n]*fetch_x(J[n], x);
		double val=V[n]*x[J[n]];

		
		rows[idx] =row;
		vals[threadIdx.x] =val;

        if(row == rows[idx -  1]) { vals[threadIdx.x] = val = val + vals[threadIdx.x -  1]; }
        if(row == rows[idx -  2]) { vals[threadIdx.x] = val = val + vals[threadIdx.x -  2]; }
        if(row == rows[idx -  4]) { vals[threadIdx.x] = val = val + vals[threadIdx.x -  4]; }
        if(row == rows[idx -  8]) { vals[threadIdx.x] = val = val + vals[threadIdx.x -  8]; }
        if(row == rows[idx - 16]) { vals[threadIdx.x] = val = val + vals[threadIdx.x - 16]; }

		if(thread_lane == 31 || n == interval_end -1){
			atomicAdd(&y[row],val);  
		}else{
			if(row != rows[idx + 1]){
					//y[row] += val;
				atomicAdd(&y[row],val);  
				
			}
		}
        
		n+=WARP_SIZE;
	}
	

}
__global__ void COO_shared2(const int num_nozeros,
		const int *I, const int *J, const double *V,
		const double *x, double *y, int testPoint)
{
	__shared__ volatile int rows[48*thread_size/WARP_SIZE];  //why using 48? because we need 16 additional junk elements
	__shared__ volatile double vals[thread_size];

	int thread_lane= threadIdx.x & (WARP_SIZE-1); //great idea! think about it
	int warp_id = threadIdx.x / WARP_SIZE;
	/*how about the interval is not the multiple of warp_size?*/
	//int iteration_end=((interval_end)/WARP_SIZE)*WARP_SIZE;

	int idx=16*(threadIdx.x/32+1) + threadIdx.x;//every warp has 16 "junk" rows elements
	rows[idx-16]=-1;
	rows[idx]=0;

	int n;
	int row;
	double val;

	for(uint16_t it = 0;  it < step_p_blk; ++it)
	{
		n = blockDim.x*blockIdx.x*step_p_blk + warp_id*WARP_SIZE*step_p_blk 
			+ it*WARP_SIZE + thread_lane;
		if(n < num_nozeros){
			row =I[n];
			//double val=V[n]*fetch_x(J[n], x);
			val=V[n]*x[J[n]];

			rows[idx] =row;
			vals[threadIdx.x] =val;

			if(row == rows[idx -  1]) { vals[threadIdx.x] = val = val + vals[threadIdx.x -  1]; __syncwarp();}
			if(row == rows[idx -  2]) { vals[threadIdx.x] = val = val + vals[threadIdx.x -  2]; __syncwarp();}
			if(row == rows[idx -  4]) { vals[threadIdx.x] = val = val + vals[threadIdx.x -  4]; __syncwarp();}
			if(row == rows[idx -  8]) { vals[threadIdx.x] = val = val + vals[threadIdx.x -  8]; __syncwarp();}
			if(row == rows[idx - 16]) { vals[threadIdx.x] = val = val + vals[threadIdx.x - 16]; __syncwarp();}

			if(thread_lane == 31 || n == num_nozeros -1){
				if(row == testPoint){
					atomicAdd(&y[row],val);  
				}else
					atomicAdd(&y[row],val);  
			}else{
				if(row != rows[idx + 1]){
					if(row == testPoint ){
						atomicAdd(&y[row],val);  
					} else
						atomicAdd(&y[row],val);  
				}
			}
		}
	}
}



__global__ void COO_atomic(const int num_nozeros, const int interval_size, 
				const int *I, const int *J, const double *V, 
				const double *x, double *y, bool tex, int testPoint) 
{
	
	int thread_id = blockDim.x*blockIdx.x + threadIdx.x;
	int thread_lane= threadIdx.x & (WARP_SIZE-1); //great idea! think about it
	int warp_id = thread_id / WARP_SIZE;
	
	int interval_begin=warp_id*interval_size;
	int interval_end =min(interval_begin+interval_size, num_nozeros);
	/*how about the interval is not the multiple of warp_size?*/
	//int iteration_end=((interval_end)/WARP_SIZE)*WARP_SIZE;
	int row;
	double val;
	int row_tmp;
	double val_tmp;
	
	int n;
	n=interval_begin+thread_lane;
	while (n< interval_end)
	{
		row = I[n];
		val = V[n]*x[J[n]];
		
		//double val=V[n]*x[J[n]];
		val_tmp = __shfl_up_sync(FULL_MASK, val, 1);
		row_tmp = __shfl_up_sync(FULL_MASK, row, 1);
		if(thread_lane > 0 && row == row_tmp) { 
			val += val_tmp; 
		} 
		val_tmp = __shfl_up_sync(FULL_MASK, val, 2);
		row_tmp = __shfl_up_sync(FULL_MASK, row, 2);
		if(thread_lane > 1 && row == row_tmp) { 
			val += val_tmp; 
		}
		val_tmp = __shfl_up_sync(FULL_MASK, val, 4);
		row_tmp = __shfl_up_sync(FULL_MASK, row, 4);
		if(thread_lane > 3 && row == row_tmp) { 
			val += val_tmp; 
		}
		val_tmp = __shfl_up_sync(FULL_MASK, val, 8);
		row_tmp = __shfl_up_sync(FULL_MASK, row, 8);
		if(thread_lane > 7 && row == row_tmp) { 
			val += val_tmp; 
		}
		val_tmp = __shfl_up_sync(FULL_MASK, val, 16);
		row_tmp = __shfl_up_sync(FULL_MASK, row, 16);
		if(thread_lane > 15 && row == row_tmp) { 
			val += val_tmp; 
		}
		row_tmp = __shfl_down_sync(FULL_MASK, row, 1);
		if(thread_lane == 31 || n == interval_end -1){
			//if(row == testPoint){
			//	y[row] += val;
			//}else
				atomicAdd(&y[row],val);  
		}else{
			if(row != row_tmp){
			//	if(row == testPoint ){
			//		y[row] += val;
			//	} else
					y[row] += val;
					//atomicAdd(&y[row],val);  
				
			}
		}	
		n+=WARP_SIZE;
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


void matrix_vectorELL(const int num_rows, const int cal_rows, 
			const int num_cols_per_row,  const int *J,
 			const double *V, const double *x, double *y,
			const bool RODR, const int rodr_blocks, const int* part_boundary_d)
{
	int ELL_blocks = ceil((double) num_rows/ELL_threadSize);
	//printf("ELL_blocks is %d\n", ELL_blocks);
	//bind_x(x);
	ELL_kernel<<<ELL_blocks, ELL_threadSize>>>(num_rows, cal_rows, num_cols_per_row, J, V, x,y);
	//unbind_x(x);
	
}

void matrix_vectorELL_block(const int num_rows, const int testPoint, 
			const int* num_cols_per_row_vec, 
			const int* block_data_bias_vec,    
			const int *J,
 			const double *V, const double *x, double *y,
			const bool CACHE, const int rodr_blocks, const int* part_boundary_d,
			const bool tex=false)
{
	int ELL_blocks = ceil((double) num_rows/ELL_threadSize);
	//printf("ELL_blocks is %d\n", ELL_blocks);
	if(rodr_blocks > 0){
		if(CACHE){	
			ELL_cached_kernel_rodr<<<rodr_blocks, ELL_threadSize>>>(num_cols_per_row_vec, 
					block_data_bias_vec,
					J, V, x, y, part_boundary_d, testPoint, tex);
			//gpuErrchk( cudaPeekAtLastError() );
		} else {
			if(testPoint > 0){
				ELL_kernel_rodr_test<<<rodr_blocks, ELL_threadSize>>>(num_cols_per_row_vec, 
						block_data_bias_vec,
						J, V, x, y, part_boundary_d, testPoint);
			} else {
				ELL_kernel_rodr<<<rodr_blocks, ELL_threadSize>>>(num_cols_per_row_vec, 
						block_data_bias_vec,
						J, V, x, y, part_boundary_d, tex);
			}
			
		}

	}else{
		ELL_kernel_block<<<ELL_blocks, ELL_threadSize>>>(num_rows, num_cols_per_row_vec, 
			block_data_bias_vec, J, V, x,y);
	}
}

//using warp interval for computing on different warps
void matrix_vectorCOO_warp(const int num_nozeros_compensation, int *I, int *J, double *V, double *x_d, double *y_d, int testPoint, bool tex=false)
{
	int interval_size2;
	interval_size2=ceil(((double) num_nozeros_compensation)/(512*512/WARP_SIZE));//for data with 2 million elements, we have interval size 200	
	//COO_atomic<<<512, 512>>>(num_nozeros_compensation, interval_size2, I, J, V, x_d, y_d, tex, testPoint);
	COO_shared<<<512, 512>>>(num_nozeros_compensation, interval_size2, I, J, V, x_d, y_d);
}

void matrix_vectorCOO(const int num_nozeros_compensation, int *I, int *J, double *V, double *x_d, double *y_d, int testPoint, bool tex=false)
{
	int blockSizeLocal;

	blockSizeLocal=ceil(((double) num_nozeros_compensation)/(step_p_blk*threadSizeCOO));//for data with 2 million elements, we have interval size 200
	//COO_atomic<<<512, 512>>>(num_nozeros_compensation, interval_size2, I, J, V, x_d, y_d, tex, testPoint);
	COO_shared2<<<blockSizeLocal, threadSizeCOO>>>(num_nozeros_compensation, I, J, V, x_d, y_d, testPoint);

}

void matrix_vectorER(const int ER_numOfRow, 
		const int* ER_rowVec_d, const int* ER_biasVec_d, 
		const int* ER_widthVec_d, 
		const int* ER_col_d, const double*  ER_val_d,
		const double* vector_in_d, const double* vector_out_d)
{

	int blockSizeLocal;
	blockSizeLocal=ceil(((double) ER_numOfRow)/threadSizeER);//for data with 2 million elements, we have interval size 200
	ER_kernel<<<blockSizeLocal, threadSizeER>>>(ER_numOfRow, ER_rowVec_d, ER_biasVec_d, ER_widthVec_d,
			ER_col_d, ER_val_d, vector_in_d, vector_out_d);

}

void matrix_vectorEHYB(matrixEHYbS* inputMatrix_d, double* vector_in_d,
		double* vector_out_d, cb_s cb, const int testPoint)
{
	int dimension = inputMatrix_d->dimension;
	int partSize = inputMatrix_d->partSize;
	int* widthVecBlockELL_d = inputMatrix_d->widthVecBlockELL;
	int* biasVecBLockELL_d = inputMatrix_d->;
	int* colBlockELL_d = inputMatrix_d->colBlockELL;
	double* valBlockELL_d = inputMatrix->valBlockELL;
	int* rowVecER_d = inputMatrix_d->rowVecER;
	int* widthVecER_d = inputMatrix_d->widthVecER;
	int* biasVecER_d = inputMatrix_d->biasVecER;
	int* part_boundary_d = inputMatrix_d->part_boundary;
	const int numOfRowER = inputMatrix->numOfRowER; 
	int* colER_d= inputMatrix->colER; 
	double* valER_d = inputMatrix->valER;
	
	matrixVectorBlockELL(dimension, testPoint, widthVecBlockELL_d, 
			biasVecBLockELL_d,
			colBlockELL_d, valBlockELL_d, vectorIn_d, vectorOut_d,
			cb.CACHE, partSize, partBoundary_d);

	if(numER > 0){
		matrixVectorER(numOfRowsER, rowVecER_d, 
				widthVecER_d,
				biasVecER_d,
				colER_d, valER_d,
				vectorIn_d, vectorOut_d);
	} 

}
