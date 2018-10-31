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

#define shared_per_block (24*1024)
#define element_size 4 //if single precision, 4, if double precision, 8 
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
void initialize_all(const unsigned int dimension, float *pk_d, float *bp_d, float *x, float *zk, const float *vector_in_d);
void initialize_bp(unsigned int num, float *x);
void initialize_r(unsigned int num, float *rk, float *vector_in);
void myxpy(const unsigned int dimension, float gamak, const float *x, float *y);
void matrix_vectorCOO(const unsigned int num_nozeros_compensation, unsigned int *I, 
		unsigned int *J, float *V, float *x, float *y, unsigned int bias0, unsigned int bias1);
void matrix_vectorELL(const unsigned int num_rows, const unsigned int cal_rows, 
		const unsigned int num_cols_per_row,  const unsigned int *J,
		const float *V, const float *x, float *y, const unsigned int bias0, const unsigned int bias1, 
		const bool RODR, const unsigned int rodr_blocks, const unsigned int* part_boundary_d);

void COO2ELL(const unsigned int *rowLocal, const unsigned int *colLocal, const float* matrixLocal, unsigned int **colELL,
	float **matrixELL, unsigned int **I_COO, unsigned int **J_COO, float **V_COO,const unsigned int *numInRow, 
	const unsigned int *rowNumAccum, const unsigned int localMatrixSize, const unsigned int localNumofRow, 
	unsigned int *sizeOut, unsigned int *Max);
#endif
