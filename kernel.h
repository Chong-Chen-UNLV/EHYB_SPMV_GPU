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
#include <cusparse_v2.h>
#include <cuda_runtime.h>
#include <cuda_runtime_api.h>
#include <cuda.h>
#include "spmv.h"
#define warpSize  32

const int memPerThread = 72;
const int threadELL = 1024;
const int warpPerBlock = threadELL/warpSize;
const int sharedPerBlock = memPerThread*threadELL;//1024 is the maximum threads per block
const int elementSize = 8; //if single precision, 4, if double precision, 8 
const int loopInKernel =  memPerThread/elementSize;
const int vectorCacheSize = sharedPerBlock/elementSize;
const int blockPerPart = sharedPerBlock/(warpSize*elementSize); 

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

void matrixVectorEHYB_NC(matrixEHYB* inputMatrix, 
		//int16_t* biasIdxBlock_d,
		double* vector_in_d,
		double* vector_out_d, const int testPoint);

void matrixVectorEHYB(matrixEHYB* inputMatrix, 
		//int16_t* biasIdxBlock_d,
		double* vector_in_d,
		double* vector_out_d, const int testPoint);

#endif
