#ifndef KERNEL_H
#define KERNEL_H

#include <stdlib.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <sys/time.h>
#include <time.h>
#include <unistd.h>
//#include <mpi.h>

#include <cublas_v2.h>
#include <cuda_runtime.h>
#include <cuda_runtime_api.h>
#include <cuda.h>
#include "solver.h"

const int ELL_threadSize = 512; 

const int shared_per_block = 24*1024;
const int element_size = 8; //if single precision, 4, if double precision, 8 
const int vector_cache_size = shared_per_block/element_size;
const int block_per_part = shared_per_block/(ELL_threadSize*element_size); 

#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
{
   if (code != cudaSuccess) 
   {
      fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
      if (abort) exit(code);
   }
}


extern "C"
void initialize_all(const uint32_t dimension, double *pk_d, double *bp_d, double *x, double *zk, const double *vector_in_d);
void initialize_bp(uint32_t num, double *x);
void initialize_r(uint32_t num, double *rk, double *vector_in);
void myxpy(const uint32_t dimension, double gamak, const double *x, double *y);
void matrix_vectorCOO(const uint32_t num_nozeros_compensation, uint32_t *I, 
		uint32_t *J, double *V, double *x, double *y);

void matrix_vectorELL(const uint32_t num_rows, const uint32_t cal_rows, 
		const uint32_t num_cols_per_row,  const uint32_t *J,
		const double *V, const double *x, double *y, 
		const bool RODR, const uint32_t rodr_blocks, const uint32_t* part_boundary_d);

void matrix_vectorELL_block(const uint32_t num_rows, const uint32_t cal_rows, 
			const uint32_t* num_cols_per_row_vec, 
			const uint32_t* block_data_bias_vec,    
			const uint32_t *J,
 			const double *V, const double *x, double *y,
			const bool RODR, const uint32_t rodr_blocks, const uint32_t* part_boundary_d);
void matrix_vectorHYB(matrixHYB_S_d* inputMatrix, double* vector_in_d,
		double* vector_out_d, cb_s cb, const uint32_t testPoint,
		const uint32_t part_size, const uint32_t* part_boundary);

#endif
