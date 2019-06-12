#ifndef KERNEL_H
#define KERNEL_H

#include <stdlib.h>
#include <stdbool.h>
#include <stdio.h>
#include <math.h>
#include <sys/time.h>
#include <time.h>
#include <unistd.h>
//#include <mpi.h>

#include <cublas_v2.h>
#include <cuda_runtime.h>
#include <cuda_runtime_api.h>
#include <cuda.h>

#define ELL_threadSize 1024 

#define shared_per_block (28*1024)
#define element_size 8 //if single precision, 4, if double precision, 8 
#define vector_cache_size shared_per_block/element_size

//#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
//inline void gpuAssert(cudaError_t code, const char *file, unsigned int line, bool abort=true)
//{
//   abort=true;
//   if (code != cudaSuccess) 
//   {
//      fprunsigned intf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
//      if (abort) exit(code);
//   }
//}
extern "C"
void initialize_all(const unsigned int dimension, double *pk_d, double *bp_d, double *x, double *zk, const double *vector_in_d);
void initialize_bp(unsigned int num, double *x);
void initialize_r(unsigned int num, double *rk, double *vector_in);
void myxpy(const unsigned int dimension, double gamak, const double *x, double *y);
void matrix_vectorCOO(const unsigned int num_nozeros_compensation, unsigned int *I, 
		unsigned int *J, double *V, double *x, double *y, unsigned int bias0, unsigned int bias1);
void matrix_vectorELL(const unsigned int num_rows, const unsigned int cal_rows, 
		const unsigned int num_cols_per_row,  const unsigned int *J,
		const double *V, const double *x, double *y, const unsigned int bias0, const unsigned int bias1, 
		const bool RODR, const unsigned int rodr_blocks, const unsigned int* part_boundary_d);

void matrix_vectorELL_block(const unsigned int num_rows, const unsigned int cal_rows, 
			const unsigned int* num_cols_per_row_vec, 
			const unsigned int* block_data_bias_vec,    
			const unsigned int *J,
 			const double *V, const double *x, double *y, const unsigned int bias0, const unsigned int bias1, 
			const bool RODR, const unsigned int rodr_blocks, const unsigned int* part_boundary_d);

#endif
