#include "kernel.h"

#define BASE 262144 //1024*1024

#define block_size 512	
#define thread_size 512
#define block_size2 16
#define thread_size2 512
#define WARP_SIZE 32


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
double get_val(const unsigned int idx, const unsigned int scope1, const unsigned int scope2,  const double *vec, volatile double* cached_vec){
	if(idx > scope1 && idx < scope2)
		return cached_vec[idx - scope1];
	else
		return vec[idx];
}

/*kernel function for initialize*/
__global__ void kernelInitialize(const unsigned int num, double *x)
{
	unsigned int idx=blockDim.x * blockIdx.x+ threadIdx.x;
	
	for (unsigned int n=idx;n<num;n+=BASE) x[n]=0;
}

__global__ void kernelInitializeAll(const unsigned int num, double *pk, double *bp, double *x, double *zk, const double *vector_in)
{
	unsigned int idx=blockDim.x * blockIdx.x+ threadIdx.x;
	double temp;
	for (unsigned int n=idx;n<num;n+=BASE) 
	{
		temp=zk[n];
		pk[n]=temp;
		bp[n]=0;
		x[n]=0;
	}
}

__global__ void kernelInitializeR(const unsigned int num,double *rk, const double *vector_in)
{
	unsigned int idx=blockDim.x * blockIdx.x+ threadIdx.x;
	double temp;
	for (unsigned int n=idx;n<num;n+=BASE) 
	{
		temp=vector_in[n];
		rk[n]=temp;
	}
}

//for ELL format matrix, output y=data*x
__global__ void ELL_kernel(const unsigned int num_rows, const unsigned int cal_rows, const unsigned int num_cols_per_row,
			const unsigned int *indices, const double *data, const double * x, double * y, 
			const unsigned int bias0, const unsigned int bias1)
{
	unsigned int row= blockDim.x*blockIdx.x+threadIdx.x;
	if (row<cal_rows){
		double dot =0;
		for (unsigned int n=0; n< num_cols_per_row; n++){
			unsigned int col=indices[num_rows * n + row];
			double val=data[num_rows*n+row];

			if (val != 0)
				dot += val* x[col-bias0];
		}
		y[row+bias1]+=dot;
	}
}
/*
__global__ void ELL_kernel_float(const unsigned int num_rows, const unsigned int cal_rows, const unsigned int num_cols_per_row,
			const unsigned int *indices, const float* data, const double * x, double * y, 
			const unsigned int bias0, const unsigned int bias1)
{
	unsigned int row= blockDim.x*blockIdx.x+threadIdx.x;
	if (row<cal_rows){
		double dot =0;
		for (unsigned int n=0; n< num_cols_per_row; n++){
			unsigned int col=indices[num_rows * n + row];
			double val=data[num_rows*n+row];

			if (val != 0)
				dot += val* x[col-bias0];
		}
		y[row+bias1]+=dot;
	}
}*/
__global__ void ELL_kernel_block(const unsigned int num_rows, const unsigned int cal_rows, 
			const unsigned int* num_cols_per_row_vec, const unsigned int* block_data_bias_vec,  
			const unsigned int *indices, const double *data, const double * x, double * y, 
			const unsigned int bias0, const unsigned int bias1){
	
	unsigned int block_idx = blockIdx.x; 
	unsigned int thread_idx = threadIdx.x; 
	unsigned int num_cols_per_row;  
	unsigned int block_data_bias;
	unsigned int data;
	unsigned int col;	
	double val;
	unsigned int row= blockDim.x*blockIdx.x+threadIdx.x;
	if(row < cal_rows){
		num_cols_per_row = num_cols_per_row_vec[block_idx];//cache will work when every threads read same global address
		block_data_bias = block_data_bias_vec[block_idx];
		double dot =0;
		for (unsigned int n=0; n< num_cols_per_row; n++){
			data_idx = block_data_bias + ELL_threadSize*n + thread_idx;
			col=indices[data_idx];
			val=data[data_idx];

			if (val != 0)
				dot += val* x[col-bias0];
		}
		y[row+bias1]+=dot;
	}
}

/* bias0 and bias1 is reserved for future distributed version*/
__global__ void ELL_cached_kernel(const unsigned int num_rows,  
				const unsigned int num_cols_per_row, 
				const unsigned int *indices, const double *data, const double * x,
				double * y, const unsigned int bias0, 
				const unsigned int bias1, const unsigned int* part_boundary)
{
	unsigned int x_idx = blockDim.x*blockIdx.x+threadIdx.x;
	__shared__ volatile double cached_vec[vector_cache_size];  
	unsigned int vec_start = part_boundary[blockIdx.x] + bias0;
	unsigned int vec_end = part_boundary[blockIdx.x + 1] + bias0;
	double val, dot;
	unsigned int col;

	for (unsigned int i = x_idx; i < vector_cache_size; i += ELL_threadSize){
		if(i < vec_end) cached_vec[i] = x[i + vec_start];
		else cached_vec[i] = 0;
	}
	for(unsigned int row = x_idx; row < vec_end; row += ELL_threadSize){
		dot =0;
		for (unsigned int n=0; n< num_cols_per_row; n++){
			col=indices[num_rows*n + row];
			val=data[num_rows*n + row];
			if (val != 0)
				dot += val*get_val(col, vec_start, vec_start + vector_cache_size, x, cached_vec);
		}
		y[row+bias1] = dot;
	}		
}

__global__ void ELL_cached_kernel_block(const unsigned int* num_cols_per_row_vec, 
		const unsigned int* block_data_bias_vec,  
		const unsigned int *indices, const double *data, const double * x,
		double * y, const unsigned int bias0, 
		const unsigned int bias1, const unsigned int* part_boundary){

	unsigned int block_idx = blockIdx.x; 
	unsigned int x_idx = blockDim.x*blockIdx.x+threadIdx.x;
	__shared__ volatile double cached_vec[vector_cache_size];  
	unsigned int vec_start = part_boundary[blockIdx.x] + bias0;
	unsigned int vec_end = part_boundary[blockIdx.x + 1] + bias0;
	double val, dot;
	unsigned int col;
	num_cols_per_row = num_cols_per_row_vec[block_idx];//cache will work when every threads read same global address
	//vec_start + vector_cache_size will be slightly different from vec_end
	for (unsigned int i = x_idx; i < vector_cache_size; i += ELL_threadSize){
		if(i < vec_end) cached_vec[i] = x[i + vec_start];
		else cached_vec[i] = 0;
	}

	block_data_bias = block_data_bias_vec[block_idx];
	for(unsigned int row = x_idx; row < vec_end; row += ELL_threadSize){
		dot =0;
		for (unsigned int n=0; n< num_cols_per_row; n++){
			data_idx = block_data_bias + ELL_threadSize*n + row;
			col=indices[data_idx];
			val=data[data_idx];
			if (val != 0)
				dot += val*get_val(col, vec_start, vec_start + vector_cache_size, x, cached_vec);
		}
		y[row+bias1] = dot;
	}		
}
/*
__global__ void ELL_cached_kernel_float(const unsigned int num_rows,  
				const unsigned int num_cols_per_row, 
				const unsigned int *indices, const float* data, const double * x,
				double * y, const unsigned int bias0, 
				const unsigned int bias1, const unsigned int* part_boundary)
{
	unsigned int x_idx = blockDim.x*blockIdx.x+threadIdx.x;
	__shared__ volatile double cached_vec[vector_cache_size];  
	unsigned int vec_start = part_boundary[blockIdx.x] + bias0;
	unsigned int vec_end = part_boundary[blockIdx.x + 1] + bias0;
	double val, dot;
	unsigned int col;

	for (unsigned int i = x_idx; i < vector_cache_size; i += ELL_threadSize){
		cached_vec[i] = x[i + vec_start];
	}
	for(unsigned int row = x_idx; row < vec_end; row += ELL_threadSize){
		dot =0;
		for (unsigned int n=0; n< num_cols_per_row; n++){
			col=indices[num_rows*n + row];
			val=data[num_rows*n + row];
			if (val != 0)
				dot += val*get_val(col, vec_start, vec_start + vector_cache_size, x, cached_vec);
		}
		y[row+bias1] = dot;
	}		
}*/

//for COO format matrix, output y=data*x
//the basic idea is come from
__global__ void COO_level1(const unsigned int num_nozeros, const unsigned int interval_size, 
				const unsigned int *I, const unsigned int *J, const double *V, 
				const double *x, double *y, int *temp_rows, 
				double *temp_vals, const unsigned int xp,const unsigned int yp)
{
	__shared__ volatile int rows[48*thread_size/WARP_SIZE];  //why using 48? because we need 16 additional junk elements
	__shared__ volatile double vals[thread_size];
	
	unsigned int thread_id = blockDim.x*blockIdx.x + threadIdx.x;
	unsigned int thread_lane= threadIdx.x & (WARP_SIZE-1); //great idea! think about it
	unsigned int warp_id = thread_id / WARP_SIZE;
	
	unsigned int interval_begin=warp_id*interval_size;
	unsigned int interval_end =min(interval_begin+interval_size, num_nozeros);
	/*how about the interval is not the multiple of warp_size?*/
	//unsigned int iteration_end=((interval_end)/WARP_SIZE)*WARP_SIZE;
	
	unsigned int idx=16*(threadIdx.x/32+1) + threadIdx.x;//every warp has 16 "junk" rows elements
	
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
	unsigned int n;
	n=interval_begin+thread_lane;
	while (n< interval_end)
	{
		unsigned int row =I[n];
		//double val=V[n]*fetch_x(J[n], x);
		double val=V[n]*x[J[n]-xp];
		
		if (thread_lane==0)
		{
			if (row==rows[idx+31])
				val+=vals[threadIdx.x+31]; //don't confused by the "plus" 31, because the former end is the new start
			else 
				y[rows[idx+31]-yp] += vals[threadIdx.x+31];//try to fix the bug from orignial library functions
		}
		rows[idx] =row;
		vals[threadIdx.x] =val;
		
        if(row == rows[idx -  1]) { vals[threadIdx.x] = val = val + vals[threadIdx.x -  1]; } 
        if(row == rows[idx -  2]) { vals[threadIdx.x] = val = val + vals[threadIdx.x -  2]; }
        if(row == rows[idx -  4]) { vals[threadIdx.x] = val = val + vals[threadIdx.x -  4]; }
        if(row == rows[idx -  8]) { vals[threadIdx.x] = val = val + vals[threadIdx.x -  8]; }
        if(row == rows[idx - 16]) { vals[threadIdx.x] = val = val + vals[threadIdx.x - 16]; }

        if(thread_lane < 31 && row < rows[idx + 1] && n<interval_end-1)
            y[row-yp] += vals[threadIdx.x];  
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

/* The second level of the segmented reduction operation
Why we need second level of reduction? because the program running at different block can not be sychronized 
Notice the number of input elements is fixed, and the number is relatively much small than the dimension of matrixs
consider block_size=512, thread_size=512 (i.e. 512 block, each block has 512 threads) and wrapSize=32 the dimension of
temp_rows will be 512*512/32=8192, this is a fix number
So, we should set this device function's block_size=512/32=16, thread_size=512*/

__global__ void COO_level2(const int * temp_rows,
                       	const double * temp_vals,
			int * temp_rows2,
			double * temp_vals2,
                       	double * y, const unsigned int p)
/*The bias is */									
{
    __shared__ int rows[thread_size2 + 1];    
    __shared__ double vals[thread_size2 + 1];
	unsigned int idx_t=threadIdx.x;
	unsigned int idx_g=blockDim.x*blockIdx.x+threadIdx.x;
	
    if (threadIdx.x == 0)
    {
        rows[thread_size2] =  -1;
        vals[thread_size2] =   0;
	temp_rows2[blockIdx.x]=-1;
    }
    	
	rows[idx_t]=temp_rows[idx_g];
	vals[idx_t]=temp_vals[idx_g];
	__syncthreads();
	
	segreduce_block(rows, vals);
	
	if (rows[threadIdx.x] != rows[threadIdx.x + 1])
	{
		if (threadIdx.x!=(thread_size2-1))
		{
			if (rows[threadIdx.x]>=0 && rows[threadIdx.x]>=p) y[rows[threadIdx.x]-p] += vals[threadIdx.x];
		}
		else
		{
			temp_rows2[blockIdx.x]=rows[threadIdx.x];
			temp_vals2[blockIdx.x]=vals[threadIdx.x];
		}
	}		
	
}

//no sychronize between blocks, so we need to restart another kernel function
__global__ void COO_level3(const unsigned int num,
                            const int * temp_rows,
                            double * temp_vals,
                            double * y,const unsigned int p)
{
	/*only 16 elements, single thread is enough*/
	unsigned int i=0;
	for (i=0;i<num-1;i++)
	{
		if (temp_rows[i]!=temp_rows[i+1])
			if (temp_rows[i]>=0 && temp_rows[i]>=p) y[temp_rows[i]-p] +=temp_vals[i];
		else if (temp_rows[i]==temp_rows[i+1])
			temp_vals[i+1]=temp_vals[i+1]+temp_vals[i]; //don't forget update the values! also at most situation (sparse matrix) it is unnecessary	
	}
	/*the last elements of input data will not disturb by any other elements, so update the output directly*/
	if (temp_rows[i]>=0 && temp_rows[i]>=p) 
		y[temp_rows[i]-p] +=temp_vals[i];
}

/*The single thread version of reduction*/
__global__ void COO_level2_serial(const int * temp_rows,
                              const double * temp_vals,
                                    double * y,const unsigned int p)
{
	unsigned int i=0;
	for (i=0;i<(block_size*thread_size/WARP_SIZE);i++)
	{
		if (temp_rows[i]>=0 && temp_rows[i]>=p) 
			y[temp_rows[i]-p]+=temp_vals[i];
	}
}

__global__ void COO_level2_serial2(const int * temp_rows,
                              const double * temp_vals,
                                    double * y, const unsigned int p)
{
	unsigned int i=0;
	for (i=0;i<(block_size2*thread_size2/WARP_SIZE);i++)
	{
		
	if (temp_rows[i]>=0 && temp_rows[i]>=p) 
 		y[temp_rows[i]-p]+=temp_vals[i];
	}
}

__global__ void COO_level2_serial3(const unsigned int num, const int * temp_rows,
                              const double * temp_vals,
                                    double * y,const unsigned int p)
{
	unsigned int i=0;
	for (i=0;i<num;i++)
	{
		
	if (temp_rows[i]>=0 && temp_rows[i]>=p) 
 		y[temp_rows[i]-p]+=temp_vals[i];
	}
}

__global__ void COO_level1_serial(const unsigned int num, unsigned int *I, unsigned int *J, double *V, double *x, double *y, const unsigned int xp, const unsigned int yp)
{
	unsigned int i=0;
	for (i=0;i<num;i++)
	{
		if (I[i]>=yp) y[I[i]-yp]+=V[i]*x[J[i]-xp];
	}
}

//y=x+gamak*y
__global__ void kernelMyxpy(const unsigned int dimension, double gamak, const double *x, double *y)
{
	unsigned int idx=blockDim.x*blockIdx.x+threadIdx.x;
	unsigned int n=idx;
	while(n<dimension){
		y[n]=x[n]+gamak*y[n];
		n=n+BASE;
	}
}
extern "C"
void initialize_all(const unsigned int dimension, double *pk_d, double *bp_d, double *x, double *zk, const double *vector_in_d)
{
	kernelInitializeAll<<<block_size,thread_size>>>(dimension, pk_d, bp_d, x, zk, vector_in_d);
}

void initialize_bp(unsigned int num, double *x)
{
	kernelInitialize<<<block_size,thread_size>>>(num,x);
}

void initialize_r(unsigned int num, double *rk, double *vector_in)
{
	kernelInitializeR<<<block_size,thread_size>>>(num,rk,vector_in);
}
void myxpy(const unsigned int dimension, double gamak, const double *x, double *y)
{
	kernelMyxpy<<<block_size,thread_size>>>(dimension,gamak,x,y);
}

void initialDeviceArray(unsigned int num, double *x)
{
	kernelInitialize<<<512,512>>>(num,x);
}


void matrix_vectorELL(const unsigned int num_rows, const unsigned int cal_rows, 
			const unsigned int num_cols_per_row,  const unsigned int *J,
 			const double *V, const double *x, double *y, const unsigned int bias0, const unsigned int bias1, 
			const bool RODR, const unsigned int rodr_blocks, const unsigned int* part_boundary_d)
{
	/*bias0 is for x and bias1 is for y, in precond solver, x, y may have different start point, 
		bias0 is "absolut bias", bias 1 is relative bias*/
	unsigned int ELL_blocks = ceil((double) num_rows/ELL_threadSize);
	//printf("ELL_blocks is %d\n", ELL_blocks);
	//bind_x(x);
	if(RODR){
		
		ELL_cached_kernel<<<rodr_blocks, ELL_threadSize>>>(num_rows, num_cols_per_row, J,
 			V, x,y, bias0, bias1, part_boundary_d);	
	} else { 
		ELL_kernel<<<ELL_blocks, ELL_threadSize>>>(num_rows, cal_rows, num_cols_per_row, J, V, x,y, bias0, bias1);
	}
	//unbind_x(x);
	
}

/*void matrix_vectorELL_float(const unsigned int num_rows, const unsigned int cal_rows, 
			const unsigned int num_cols_per_row,  const unsigned int *J,
 			const float* V, const double *x, double *y, const unsigned int bias0, const unsigned int bias1, 
			const bool RODR, const unsigned int rodr_blocks, const unsigned int* part_boundary_d)
{
	/*bias0 is for x and bias1 is for y, in precond solver, x, y may have different start point,
		bias0 is "absolut bias", bias 1 is relative bias*
	unsigned int ELL_blocks = ceil((double) num_rows/ELL_threadSize);
	//printf("ELL_blocks is %d\n", ELL_blocks);
	//bind_x(x);
	if(RODR){
		
		ELL_cached_kernel_float<<<rodr_blocks, ELL_threadSize>>>(num_rows, num_cols_per_row, J,
 			V, x,y, bias0, bias1, part_boundary_d);	
	} else { 
		ELL_kernel_float<<<ELL_blocks, ELL_threadSize>>>(num_rows, cal_rows, num_cols_per_row, J, V, x,y, bias0, bias1);
	}
	//unbind_x(x);
	
}*/

void matrix_vectorELL_block(const unsigned int num_rows, const unsigned int cal_rows, 
			const unsigned int* num_cols_per_row_vec, 
			const unsigned int* block_data_bias_vec,    
			const unsigned int *J,
 			const double *V, const double *x, double *y, const unsigned int bias0, const unsigned int bias1, 
			const bool RODR, const unsigned int rodr_blocks, const unsigned int* part_boundary_d)
{
	/*bias0 is for x and bias1 is for y, in precond solver, x, y may have different start point, 
		bias0 is "absolut bias", bias 1 is relative bias*/
	unsigned int ELL_blocks = ceil((double) num_rows/ELL_threadSize);
	//printf("ELL_blocks is %d\n", ELL_blocks);
	//bind_x(x);
	if(RODR){
		
		ELL_cached_kernel_block<<<rodr_blocks, ELL_threadSize>>>(num_rows, num_cols_per_row_vec, 
			block_data_bias_vec,
			J, V, x, y, bias0, bias1, part_boundary_d);
	}
	else
		ELL_kernel_block<<<ELL_blocks, ELL_threadSize>>>(num_rows, cal_rows, num_cols_per_row_vec, 
			block_data_bias_vec, J, V, x,y, bias0, bias1);
	//unbind_x(x);
	
}

void matrix_vectorCOO(const unsigned int num_nozeros_compensation, unsigned int *I, unsigned int *J, double *V, double *x_d, double *y_d, unsigned int bias0, unsigned int bias1)
{
	/*bias0 is for input vector, bias1 is for output vector, different from the ELL format both bias0 and bias1 is absolut bias*/
	unsigned int interval_size2;
	interval_size2=ceil(((double) num_nozeros_compensation)/(block_size*thread_size/WARP_SIZE));//for data with 2 million elements, we have interval size 200	
	//printf("num_nozeros_compensation is %d, intervalSize is %d\n",num_nozeros_compensation, interval_size2 );
	if (interval_size2>2*32)
	{
		//512*512
		size_t sizeKernel0=(block_size*thread_size/WARP_SIZE)*sizeof(unsigned int);
		size_t sizeKernel1=(block_size*thread_size/WARP_SIZE)*sizeof(double);		
		int *temp_rows1;
		double *temp_vals1;
		int *temp_rows2;
		double *temp_vals2;
		cudaMalloc((void**)&temp_rows1, sizeKernel0);
		cudaMalloc((void**)&temp_vals1, sizeKernel1);
		cudaMalloc((void**)&temp_rows2, block_size2*sizeof(unsigned int));
		cudaMalloc((void**)&temp_vals2, block_size2*sizeof(double));
		COO_level1<<<block_size,thread_size>>>(num_nozeros_compensation,interval_size2, 
					I, J, V, x_d, y_d, temp_rows1, temp_vals1, bias0, bias1);
		COO_level2<<<block_size2,thread_size2>>>(temp_rows1,temp_vals1,temp_rows2,temp_vals2,y_d, bias1);
		COO_level3<<<1,1>>>(block_size2,temp_rows2,temp_vals2,y_d, bias1);
		//COO_level2_serial<<<1,1>>>(temp_rows1,temp_vals1,y_d, bias1);
	}
	else if (interval_size2>1)
	//if (interval_size2>32)
	{
		//16*512
		//printf("situation 2 happened!\n");
		size_t sizeKernel2=(block_size2*thread_size2/WARP_SIZE)*sizeof(unsigned int);
		size_t sizeKernel3=(block_size2*thread_size2/WARP_SIZE)*sizeof(double);
		int *temp_rows3;
		double *temp_vals3;
		cudaMalloc((void**)&temp_rows3, sizeKernel2);
		cudaMalloc((void**)&temp_vals3, sizeKernel3);
		COO_level1<<<block_size2,thread_size2>>>(num_nozeros_compensation, interval_size2, 
				I, J, V, x_d, y_d, temp_rows3, temp_vals3, bias0, bias1);
		//512 calculation excuted serially
		COO_level2_serial2<<<1,1>>>(temp_rows3,temp_vals3,y_d, bias1);		
	}
	/*else if (interval_size2>4)
	{
		//16*32
		unsigned int iterval_size3=ceil((double) num_nozeros_compensation/(512*2/32));
		unsigned int *temp_rows4;
		double *temp_vals4;
		cudaMalloc((void**)&temp_rows4, 32*sizeof(unsigned int));
		cudaMalloc((void**)&temp_vals4, 32*sizeof(double));		
		COO_level1<<<2,512>>>(num_nozeros_compensation, iterval_size3, I,J,V,x_d,y_d,temp_rows4,temp_vals4,bias0,bias1);
		//16 calculation excuted serially
		COO_level2_serial3<<<1,1>>>(32, temp_rows4,temp_vals4,y_d, bias1);		
	}*/
	else
	{
		//less than 512, all calculation excuted serially
		//printf("situation 3 happen\n");
		COO_level1_serial<<<1,1>>>(num_nozeros_compensation, I,J,V,x_d,y_d, bias0,bias1);
	}
	
}

/*void matrix_vectorHYP(const unsigned int num_rows, const unsigned int max, const unsigned int *J, const double *V, unsigned int NumCOO, unsigned int *I_COO, unsigned int *J_COO, double *V_COO, )
{
	matrix_vectorELL();
	matrix_vectorCOO();
}*/
