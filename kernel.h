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

const int ELLthreadSize = 512; 

const int sharedPerBlock = 32*1024;
const int elementSize = 8; //if single precision, 4, if double precision, 8 
const int vectorCacheSize = sharedPerBlock/elementSize;
const int blockPerPart = sharedPerBlock/(EllThreadSize*elementSize); 
const int stepPerblk = 8;
const int threadSizeCOO= 256;

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
void initialize_all(const int dimension, double *pk_d, double *bp_d, double *x, double *zk, const double *vector_in_d);
void initialize_bp(int num, double *x);
void initialize_r(int num, double *rk, double *vector_in);
void myxpy(const int dimension, double gamak, const double *x, double *y);
void matrix_vectorCOO(const int num_nozeros_compensation, int *I, 
		int *J, double *V, double *x, double *y, int testPoint, const bool tex);

void matrix_vectorELL(const int num_rows, const int cal_rows, 
		const int num_cols_per_row,  const int *J,
		const double *V, const double *x, double *y, 
		const bool RODR, const int rodr_blocks, const int* part_boundary_d);

void matrix_vectorELL_block(const int num_rows, const int cal_rows, 
			const int* num_cols_per_row_vec, 
			const int* block_data_bias_vec,    
			const int *J,
 			const double *V, const double *x, double *y,
			const bool RODR, const int rodr_blocks, const int* part_boundary_d, const bool tex);
void matrix_vectorHYB(matrixHYB_S_d* inputMatrix, double* vector_in_d,
		double* vector_out_d, cb_s cb, const int testPoint,
		const int part_size, const int* part_boundary, const bool tex);

void matrix_vectorEHYB(matrixHYB_S_d* inputMatrix, double* vector_in_d,
		double* vector_out_d, cb_s cb, const int testPoint,
		const int part_size, const int* part_boundary, const bool tex);

#endif
