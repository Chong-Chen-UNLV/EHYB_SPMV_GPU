#include "kernel.h"

#define FULL_MASK 0xffffffff
#define BASE 262144 //1024*1024

#define block_size 512	
#define thread_size 512
#define block_size2 16
#define thread_size2 512
#define WARP_SIZE 32


texture<double, 1, cudaReadModeElementType> texInput;

texInput.addressMode[0] = cudaAddressModeBorder;
texInput.addressMode[1] = cudaAddressModeBorder;
texInput.filterMode = cudaFilterModePoint;
texInput.normalized = false;


/*the device function for level 2 reduce*/
__device__ void segreduce_block(const int * idx, double * val)
{
    double left = 0;
    if( threadIdx.x >=   1 && idx[threadIdx.x] == idx[threadIdx.x -   1] ) { left = val[threadIdx.x -   1]; } __syncthreads(); val[threadIdx.x] += left; left = 0; __syncthreads();  
    if( threadIdx.x >=   2 && idx[threadIdx.x] == idx[threadIdx.x -   2] ) { left = val[threadIdx.x -   2]; } __syncthreads(); val[threadIdx.x] += left; left = 0; __syncthreads();
    if( threadIdx.x >=   4 && idx[threadIdx.x] == idx[threadIdx.x -   4] ) { left = val[threadIdx.x -   4]; } __syncthreads(); val[threadIdx.x] += left; left = 0; __syncthreads();
    if( threadIdx.x >=   8 && idx[threadIdx.x] == idx[threadIdx.x -   8] ) { left = val[threadIdx.x -   8]; } __syncthreads(); val[threadIdx.x] += left; left = 0; __syncthreads();
    if( threadIdx.x >=  16 && idx[threadIdx.x] == idx[threadIdx.x -  16] ) { left = val[threadIdx.x -  16]; } __syncthreads(); val[threadIdx.x] += left; left = 0; __syncthreads();
    if( threadIdx.x >=  32 && idx[threadIdx.x] == idx[threadIdx.x -  32] ) { left = val[threadIdx.x -  32]; } __syncthreads(); val[threadIdx.x] += left; left = 0; __syncthreads();  
    if( threadIdx.x >=  64 && idx[threadIdx.x] == idx[threadIdx.x -  64] ) { left = val[threadIdx.x -  64]; } __syncthreads(); val[threadIdx.x] += left; left = 0; __syncthreads();
    if( threadIdx.x >= 128 && idx[threadIdx.x] == idx[threadIdx.x - 128] ) { left = val[threadIdx.x - 128]; } __syncthreads(); val[threadIdx.x] += left; left = 0; __syncthreads();
    if( threadIdx.x >= 256 && idx[threadIdx.x] == idx[threadIdx.x - 256] ) { left = val[threadIdx.x - 256]; } __syncthreads(); val[threadIdx.x] += left; left = 0; __syncthreads();
	if( threadIdx.x >= 512 && idx[threadIdx.x] == idx[threadIdx.x - 512] ) { left = val[threadIdx.x - 512]; } __syncthreads(); val[threadIdx.x] += left; left = 0; __syncthreads();
}

__device__ 
double get_val(const uint32_t idx, const uint32_t scope1, const uint32_t scope2,  const double *vec, volatile double* cached_vec){
	if(idx > scope1 && idx < scope2)
		return cached_vec[idx - scope1];
	else
		return vec[idx];
}

/*kernel function for initialize*/
__global__ void kernelInitialize(const uint32_t num, double *x)
{
	uint32_t idx=blockDim.x * blockIdx.x+ threadIdx.x;
	
	for (uint32_t n=idx;n<num;n+=BASE) x[n]=0;
}

__global__ void kernelInitializeAll(const uint32_t num, double *pk, double *bp, double *x, double *zk, const double *vector_in)
{
	uint32_t idx=blockDim.x * blockIdx.x+ threadIdx.x;
	double temp;
	for (uint32_t n=idx;n<num;n+=BASE) 
	{
		temp=zk[n];
		pk[n]=temp;
		bp[n]=0;
		x[n]=0;
	}
}

__global__ void kernelInitializeR(const uint32_t num,double *rk, const double *vector_in)
{
	uint32_t idx=blockDim.x * blockIdx.x+ threadIdx.x;
	double temp;
	for (uint32_t n=idx;n<num;n+=BASE) 
	{
		temp=vector_in[n];
		rk[n]=temp;
	}
}

//for ELL format matrix, output y=data*x
__global__ void ELL_kernel(const uint32_t num_rows, const uint32_t cal_rows, const uint32_t num_cols_per_row,
			const uint32_t *indices, const double *data, const double * x, double * y) 
{
	uint32_t row= blockDim.x*blockIdx.x+threadIdx.x;
	if (row<cal_rows){
		double dot =0;
		
		for (uint32_t n=0; n< num_cols_per_row; n++){
			uint32_t col=indices[num_rows * n + row];
			double val=data[num_rows*n+row];

			if (val != 0)
				dot += val* x[col];
		}
			y[row]=dot;
	}
}
/*
__global__ void ELL_kernel_float(const uint32_t num_rows, const uint32_t cal_rows, const uint32_t num_cols_per_row,
			const uint32_t *indices, const float* data, const double * x, double * y) 
			
{
	uint32_t row= blockDim.x*blockIdx.x+threadIdx.x;
	if (row<cal_rows){
		double dot =0;
		for (uint32_t n=0; n< num_cols_per_row; n++){
			uint32_t col=indices[num_rows * n + row];
			double val=data[num_rows*n+row];

			if (val != 0)
				dot += val* x[col];
		}
		y[row]+=dot;
	}
}*/
__global__ void ELL_kernel_block(const uint32_t num_rows,
			const uint32_t* num_cols_per_row_vec, 
			const uint32_t* block_data_bias_vec,  
			const uint32_t *indices, 
			const double *data, const double * x, double * y) 
{
	
	uint32_t block_idx = blockIdx.x; 
	uint32_t thread_idx = threadIdx.x; 
	uint32_t num_cols_per_row;  
	uint32_t block_data_bias;
	uint32_t data_idx;
	uint32_t col;	
	double val;
	uint32_t row= blockDim.x*blockIdx.x+threadIdx.x;
	if(row < num_rows){
		num_cols_per_row = num_cols_per_row_vec[block_idx];//cache will work when every threads read same global address
		block_data_bias = block_data_bias_vec[block_idx];
		double dot =0;
		for (uint32_t n=0; n< num_cols_per_row; n++){
			data_idx = block_data_bias + ELL_threadSize*n + thread_idx;
			col=indices[data_idx];
			val=data[data_idx];

			if (val != 0)
				dot += val* x[col];
		}
		y[row]+=dot;
	}
}

__global__ void ELL_kernel_rodr(const uint32_t* num_cols_per_row_vec,
		const uint32_t* block_data_bias_vec,  
		const uint32_t *indices, const double *data, const double * x,
		double * y,
		const uint32_t* part_boundary,
		const bool tex)
{
	uint32_t part_idx = blockIdx.x; 
	uint32_t x_idx = threadIdx.x;
	uint32_t vec_start = part_boundary[blockIdx.x];
	uint32_t vec_end = part_boundary[blockIdx.x + 1];
	uint32_t row = 0;
	double val, dot;
	uint32_t block_idx, data_idx, col;
	uint32_t block_rowSize, block_data_bias, num_cols_per_row;
	uint32_t block_base = part_idx * block_per_part;
		
	uint32_t endBlockRow;
	if(vec_end >= vec_start+block_per_part*ELL_threadSize)
		endBlockRow = ELL_threadSize; 
	else 
		endBlockRow = vec_end - (vec_start+(block_per_part - 1)*ELL_threadSize);

	for(uint32_t i = 0; i < block_per_part; ++i){//the thread is step with stride ELL_threadSize
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
			for(uint32_t n=0; n< num_cols_per_row; ++n){
				data_idx = block_data_bias + block_rowSize*n + x_idx;//however the data storage is stride with block_rowSize
				col=indices[data_idx];
				val=data[data_idx];
				if(tex == false)
					dot += val* x[col];
				else
					dot += tex1Dfetch(texInput, col);

			}
			y[row] = dot;
		}
		block_idx += 1;
	}		
}

__global__ void ELL_kernel_rodr_test(const uint32_t* num_cols_per_row_vec,
		const uint32_t* block_data_bias_vec,  
		const uint32_t *indices, const double *data, const double * x,
		double * y,
		const uint32_t* part_boundary,
		const uint32_t testPoint)
{
	uint32_t part_idx = blockIdx.x; 
	uint32_t x_idx = threadIdx.x;
	uint32_t vec_start = part_boundary[blockIdx.x];
	uint32_t vec_end = part_boundary[blockIdx.x + 1];
	uint32_t row = 0;
	double val, dot;
	uint32_t block_idx, data_idx, col;
	uint32_t block_rowSize, block_data_bias, num_cols_per_row;
	uint32_t block_base = part_idx * block_per_part;
		
	uint32_t endBlockRow;
	if(vec_end >= vec_start+block_per_part*ELL_threadSize)
		endBlockRow = ELL_threadSize; 
	else 
		endBlockRow = vec_end - (vec_start+(block_per_part - 1)*ELL_threadSize);


	for(uint32_t i = 0; i < block_per_part; ++i){//the thread is step with stride ELL_threadSize
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
			for(uint32_t n=0; n< num_cols_per_row; ++n){
				data_idx = block_data_bias + block_rowSize*n + x_idx;//however the data storage is stride with block_rowSize
				col=indices[data_idx];
				val=data[data_idx];
				if(val != 0){
					if(row == testPoint)
						dot = dot + val* x[col];
					else
						dot += val* x[col];
				}
			}
			if(row == testPoint)
				y[row] = dot + 1 - 0.999;
			else
				y[row] = dot;
		}
		block_idx += 1;
	}		
}

__global__ void ELL_cached_kernel_rodr(const uint32_t* num_cols_per_row_vec,
		const uint32_t* block_data_bias_vec,  
		const uint32_t *indices, const double *data, const double * x,
		double * y,
		const uint32_t* part_boundary,
		const bool tex)
{
	uint32_t part_idx = blockIdx.x; 
	uint32_t x_idx = threadIdx.x;
	__shared__ volatile double cached_vec[vector_cache_size];  
	uint32_t vec_start = part_boundary[blockIdx.x];
	uint32_t vec_end = part_boundary[blockIdx.x + 1];
	uint32_t row = 0;

	for (uint32_t i = x_idx; i < vector_cache_size; i += ELL_threadSize){
		if(i < vec_end) cached_vec[i] = x[i + vec_start];
		else cached_vec[i] = 0;
	}
	double val, fetched, dot;
	uint32_t block_idx, data_idx, col;
	uint32_t block_rowSize, block_data_bias, num_cols_per_row;
	uint32_t block_base = part_idx * block_per_part;
		
	uint32_t endBlockRow;
	if(vec_end >= vec_start+block_per_part*ELL_threadSize)
		endBlockRow = ELL_threadSize; 
	else 
		endBlockRow = vec_end - (vec_start+(block_per_part - 1)*ELL_threadSize);

	for(uint32_t i = 0; i < block_per_part; ++i){//the thread is step with stride ELL_threadSize
		dot =0;
		block_rowSize = ELL_threadSize;
		block_idx = block_base + i;
		if(i == block_per_part -1){
			block_rowSize = endBlockRow;
		}
		block_data_bias = block_data_bias_vec[block_idx];
		num_cols_per_row = num_cols_per_row_vec[block_idx];
		row = i*ELL_threadSize + vec_start + x_idx;
		if(row < vec_end){
			for(uint32_t n=0; n< num_cols_per_row; ++n){
				data_idx = block_data_bias + block_rowSize*n + x_idx;//however the data storage is stride with block_rowSize
				col=indices[data_idx];
				val=data[data_idx];
				if(val != 0){
					if(col > vec_start && col < vec_start + vector_cache_size)
						fetched = cached_vec[col - vec_start];
					else{
						if(tex == false)
							fetched = x[col];
						else
							fetched = tex1Dfetch(x, col);
					}
					dot += val*fetched;
				}
			}
			y[row] = dot;
		}
		block_idx += 1;
	}		
}

/*
__global__ void ELL_cached_kernel_float(const uint32_t num_rows,  
				const uint32_t num_cols_per_row, 
				const uint32_t *indices, const float* data, const double * x,
				double * y,
				const uint32_t* part_boundary)
{
	uint32_t x_idx = blockDim.x*blockIdx.x+threadIdx.x;
	__shared__ volatile double cached_vec[vector_cache_size];  
	uint32_t vec_start = part_boundary[blockIdx.x];
	uint32_t vec_end = part_boundary[blockIdx.x + 1];
	double val, dot;
	uint32_t col;

	for (uint32_t i = x_idx; i < vector_cache_size; i += ELL_threadSize){
		cached_vec[i] = x[i + vec_start];
	}
	for(uint32_t row = x_idx; row < vec_end; row += ELL_threadSize){
		dot =0;
		for (uint32_t n=0; n< num_cols_per_row; n++){
			col=indices[num_rows*n + row];
			val=data[num_rows*n + row];
			if (val != 0)
				dot += val*get_val(col, vec_start, vec_start + vector_cache_size, x, cached_vec);
		}
		y[row] = dot;
	}		
}*/

__global__ void COO_atomic(const uint32_t num_nozeros, const uint32_t interval_size, 
				const uint32_t *I, const uint32_t *J, const double *V, 
				const double *x, double *y) 
{
	__shared__ volatile int rows[48*thread_size/WARP_SIZE];  //why using 48? because we need 16 additional junk elements
	__shared__ volatile double vals[thread_size];
	
	uint32_t thread_id = blockDim.x*blockIdx.x + threadIdx.x;
	uint32_t thread_lane= threadIdx.x & (WARP_SIZE-1); //great idea! think about it
	uint32_t warp_id = thread_id / WARP_SIZE;
	
	uint32_t interval_begin=warp_id*interval_size;
	uint32_t interval_end =min(interval_begin+interval_size, num_nozeros);
	/*how about the interval is not the multiple of warp_size?*/
	//uint32_t iteration_end=((interval_end)/WARP_SIZE)*WARP_SIZE;
	int row;
	double val;
	
	uint32_t n;
	n=interval_begin+thread_lane;
	while (n< interval_end)
	{
		uint32_t row =I[n];
		double val=V[n]*tex1Dfetch(texInput, I[n]);
		//double val=V[n]*x[J[n]];

        if(row == __shfl_up_sync(FULL_MASK, row, 1)) { val + __shfl_up_sync(FULL_MASK, val, 1); } 
        if(row == __shfl_up_sync(FULL_MASK, row, 2)) { val + __shfl_up_sync(FULL_MASK, val, 2); } 
        if(row == __shfl_up_sync(FULL_MASK, row, 4)) { val + __shfl_up_sync(FULL_MASK, val, 4); } 
        if(row == __shfl_up_sync(FULL_MASK, row, 8)) { val + __shfl_up_sync(FULL_MASK, val, 8); } 
        if(row == __shfl_up_sync(FULL_MASK, row, 16)) { val + __shfl_up_sync(FULL_MASK, val, 16); } 

        if(row != shfl_down_sync(FULL_MASK, row, 1 ) && n<interval_end-1)
           atomicAdd(&y[row],vals);  
		n+=WARP_SIZE;
	}
	
}

//for COO format matrix, output y=data*x
//the basic idea is come from
__global__ void COO_level1(const uint32_t num_nozeros, const uint32_t interval_size, 
				const uint32_t *I, const uint32_t *J, const double *V, 
				const double *x, double *y, int *temp_rows, 
				double *temp_vals)
{
	__shared__ volatile int rows[48*thread_size/WARP_SIZE];  //why using 48? because we need 16 additional junk elements
	__shared__ volatile double vals[thread_size];
	
	uint32_t thread_id = blockDim.x*blockIdx.x + threadIdx.x;
	uint32_t thread_lane= threadIdx.x & (WARP_SIZE-1); //great idea! think about it
	uint32_t warp_id = thread_id / WARP_SIZE;
	
	uint32_t interval_begin=warp_id*interval_size;
	uint32_t interval_end =min(interval_begin+interval_size, num_nozeros);
	/*how about the interval is not the multiple of warp_size?*/
	//uint32_t iteration_end=((interval_end)/WARP_SIZE)*WARP_SIZE;
	
	uint32_t idx=16*(threadIdx.x/32+1) + threadIdx.x;//every warp has 16 "junk" rows elements
	
	rows[idx-16]=-1;
	
	if(interval_begin >= interval_end)
	{
		temp_rows[warp_id] = -1;
		return;
	}	
	if (thread_lane ==31)
	{
		// initialize the cary in values
		rows[idx]=I[interval_begin];
		vals[threadIdx.x]=0;
	}
	uint32_t n;
	n=interval_begin+thread_lane;
	while (n< interval_end)
	{
		uint32_t row =I[n];
		//double val=V[n]*fetch_x(J[n], x);
		double val=V[n]*x[J[n]];
		
		if (thread_lane==0)
		{
			if (row==rows[idx+31])
				val+=vals[threadIdx.x+31]; //don't confused by the "plus" 31, because the former end is the new start
			else 
				y[rows[idx+31]] += vals[threadIdx.x+31];//try to fix the bug from orignial library functions
		}
		rows[idx] =row;
		vals[threadIdx.x] =val;
		
        if(row == rows[idx -  1]) { vals[threadIdx.x] = val = val + vals[threadIdx.x -  1]; } 
        if(row == rows[idx -  2]) { vals[threadIdx.x] = val = val + vals[threadIdx.x -  2]; }
        if(row == rows[idx -  4]) { vals[threadIdx.x] = val = val + vals[threadIdx.x -  4]; }
        if(row == rows[idx -  8]) { vals[threadIdx.x] = val = val + vals[threadIdx.x -  8]; }
        if(row == rows[idx - 16]) { vals[threadIdx.x] = val = val + vals[threadIdx.x - 16]; }

        if(thread_lane < 31 && row < rows[idx + 1] && n<interval_end-1)
            y[row] += vals[threadIdx.x];  
		n+=WARP_SIZE;
	}
	
	/*now we consider the reminder of interval_size/warp_size*/

    /*program at one warp is automatically sychronized*/
	if(n==(interval_end+WARP_SIZE-1))
    {
        // write the carry out values
        temp_rows[warp_id] = rows[idx];
        temp_vals[warp_id] = vals[threadIdx.x];
    }	
	
}



//y=x+gamak*y
__global__ void kernelMyxpy(const uint32_t dimension, double gamak, const double *x, double *y)
{
	uint32_t idx=blockDim.x*blockIdx.x+threadIdx.x;
	uint32_t n=idx;
	while(n<dimension){
		y[n]=x[n]+gamak*y[n];
		n=n+BASE;
	}
}
extern "C"
void initialize_all(const uint32_t dimension, double *pk_d, double *bp_d, double *x, double *zk, const double *vector_in_d)
{
	kernelInitializeAll<<<block_size,thread_size>>>(dimension, pk_d, bp_d, x, zk, vector_in_d);
}

void initialize_bp(uint32_t num, double *x)
{
	kernelInitialize<<<block_size,thread_size>>>(num,x);
}

void initialize_r(uint32_t num, double *rk, double *vector_in)
{
	kernelInitializeR<<<block_size,thread_size>>>(num,rk,vector_in);
}
void myxpy(const uint32_t dimension, double gamak, const double *x, double *y)
{
	kernelMyxpy<<<block_size,thread_size>>>(dimension,gamak,x,y);
}

void initialDeviceArray(uint32_t num, double *x)
{
	kernelInitialize<<<512,512>>>(num,x);
}


void matrix_vectorELL(const uint32_t num_rows, const uint32_t cal_rows, 
			const uint32_t num_cols_per_row,  const uint32_t *J,
 			const double *V, const double *x, double *y,
			const bool RODR, const uint32_t rodr_blocks, const uint32_t* part_boundary_d)
{
	uint32_t ELL_blocks = ceil((double) num_rows/ELL_threadSize);
	//printf("ELL_blocks is %d\n", ELL_blocks);
	//bind_x(x);
	ELL_kernel<<<ELL_blocks, ELL_threadSize>>>(num_rows, cal_rows, num_cols_per_row, J, V, x,y);
	//unbind_x(x);
	
}

void matrix_vectorELL_block(const uint32_t num_rows, const uint32_t testPoint, 
			const uint32_t* num_cols_per_row_vec, 
			const uint32_t* block_data_bias_vec,    
			const uint32_t *J,
 			const double *V, const double *x, double *y,
			const bool CACHE, const uint32_t rodr_blocks, const uint32_t* part_boundary_d,
			const bool tex=false)
{
	uint32_t ELL_blocks = ceil((double) num_rows/ELL_threadSize);
	//printf("ELL_blocks is %d\n", ELL_blocks);
	//bind_x(x);
	
	if(rodr_blocks > 0){
		if(CACHE){	
			ELL_cached_kernel_rodr<<<rodr_blocks, ELL_threadSize>>>(num_cols_per_row_vec, 
					block_data_bias_vec,
					J, V, x, y, part_boundary_d);
			gpuErrchk( cudaPeekAtLastError() );
		} else {
			if(testPoint > 0){
				ELL_kernel_rodr_test<<<rodr_blocks, ELL_threadSize>>>(num_cols_per_row_vec, 
						block_data_bias_vec,
						J, V, x, y, part_boundary_d, testPoint);
			} else {
				ELL_kernel_rodr<<<rodr_blocks, ELL_threadSize>>>(num_cols_per_row_vec, 
						block_data_bias_vec,
						J, V, x, y, part_boundary_d);
			}
			
		}

	}else{
		ELL_kernel_block<<<ELL_blocks, ELL_threadSize>>>(num_rows, num_cols_per_row_vec, 
			block_data_bias_vec, J, V, x,y);
	}
	if(tex == true)

		
	//unbind_x(x);
	
}

void matrix_vectorCOO(const uint32_t num_nozeros_compensation, uint32_t *I, uint32_t *J, double *V, double *x_d, double *y_d)
{
	uint32_t interval_size2;
	interval_size2=ceil(((double) num_nozeros_compensation)/(block_size*thread_size/WARP_SIZE));//for data with 2 million elements, we have interval size 200	
	COO_atomic(num_nozeros_compensation, interval_size2, I, J, V, x_d, y_d);
}

void matrix_vectorHYB(matrixHYB_S_d* inputMatrix, double* vector_in_d,
		double* vector_out_d, cb_s cb, const uint32_t testPoint,
		const uint32_t part_size, const uint32_t* part_boundary_d, 
		bool tex=false)
{
	uint32_t dimension = inputMatrix->dimension;
	uint32_t ELL_width = inputMatrix->ELL_width;
	uint32_t totalNumCOO = inputMatrix->totalNumCOO;
	uint32_t* col_d = inputMatrix->col_d;
	uint32_t* I_COO_d = inputMatrix->I_COO_d;
	uint32_t* J_COO_d = inputMatrix->J_COO_d;
	double* V_d = inputMatrix->V_d;
	double* V_COO_d = inputMatrix->V_COO_d;
	uint32_t* ELL_block_bias_vec_d = inputMatrix->ELL_block_bias_vec_d;
	uint32_t* ELL_block_cols_vec_d = inputMatrix->ELL_block_cols_vec_d;
	size_t offset = 0;
	if(tex==true)
		cudaBindTexture(&offset, texInput, vector_in_d, sizeof(double)*dimension);	
	if(!cb.BLOCK){
		matrix_vectorELL(dimension, dimension, ELL_width, col_d,V_d,
				vector_in_d, vector_out_d, false, 0, NULL);
	} else {
		if(cb.RODR){
			matrix_vectorELL_block(dimension, testPoint, ELL_block_cols_vec_d, 
					ELL_block_bias_vec_d,
					col_d,V_d, vector_in_d, vector_out_d,
					cb.CACHE, part_size, part_boundary_d);
		}
		else{
			matrix_vectorELL_block(dimension, testPoint, ELL_block_cols_vec_d, 
					ELL_block_bias_vec_d,
					col_d, V_d, 
					vector_in_d, vector_out_d,
					false, 0, NULL);
		}
	}

	if (totalNumCOO > 0) matrix_vectorCOO(totalNumCOO, I_COO_d, J_COO_d, V_COO_d, 
			vector_in_d, vector_out_d);

	if(tex==true)
		cudaUnbindTexture(texInput);
}
