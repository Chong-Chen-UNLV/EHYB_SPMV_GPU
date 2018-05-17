#ifndef KERNEL_H
#define KERNEL_H

#include <stdlib.h>
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

#define shared_per_block 24*1024
#define element_size 4 //if single precision, 4, if double precision, 8 
#define vector_cache_size shared_per_block/element_size

#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
{
   if (code != cudaSuccess) 
   {
      fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
      if (abort) exit(code);
   }
}

extern "C"{

void initialize_all(const int dimension, float *pk_d, float *bp_d, float *x, float *zk, const float *vector_in_d);
void initialize_bp(int num, float *x);
void initialize_r(int num, float *rk, float *vector_in);
void myxpy(const int dimension, float gamak, const float *x, float *y);
void matrix_vectorCOO(const int num_nozeros_compensation, int *I, int *J, float *V, float *x, float *y, int bias0, int bias1);
void matrix_vectorELL(const int num_rows, const int cal_rows, const int num_cols_per_row,  const int *J,
		const float *V, const float *x, float *y, const int bias0, const int bias1, 
		const bool RODR, const int rodr_blocks, const int* part_boundary_d);

void COO2ELL(const int *rowLocal, const int *colLocal, const float* matrixLocal, int **colELL,
	float **matrixELL, int **I_COO, int **J_COO, float **V_COO,const int *numInRow, 
	const int *rowNumAccum, const int localMatrixSize, const int localNumofRow, 
	int *sizeOut, int *max);

}
#endif
