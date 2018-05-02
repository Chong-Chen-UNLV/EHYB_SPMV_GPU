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

extern "C"{

void initialize_all(const int dimension, float *pk_d, float *bp_d, float *x, float *zk, const float *vector_in_d);
void initialize_bp(int num, float *x);
void initialize_r(int num, float *rk, float *vector_in);
void myxpy(const int dimension, float gamak, const float *x, float *y);
void matrix_vectorCOO(const int num_nozeros_compensation, int *I, int *J, float *V, float *x, float *y, int bias0, int bias1);
void matrix_vectorELL(const int num_rows, const int cal_rows, const int num_cols_per_row,  const int *J,
		const float *V, const float *x, float *y, const int bias0, const int bias1, 
		const bool RODR, const float* part_boundary_d);

void COO2ELL(const int *rowLocal, const int *colLocal, const float* matrixLocal, int **colELL,
	float **matrixELL, int **I_COO, int **J_COO, float **V_COO,const int *numInRow, 
	const int *rowNumAccum, const int localMatrixSize, const int localNumofRow, 
	int *sizeOut, int *max);

void solverCPU_precond(const int dimension, const int totalNum, const int *I, const int *J,
		 	const float *V, const int totalNumPrecond, const int *I_precond,
			const int *J_precond, const float *V_precond, const float *vector_in,
			float *vector_out, float *error_track, const int MAXIter);
void solverCPU(const int dimension, const int totalNum, const int *I, const int *J, 
		const float *V, const int totalNumPrecond, const int *I_precond, 
		const int *J_precond, const float *V_precond, const float *vector_in, 
		float *vector_out, float *error_track);
void solverGPU_HYB(const int dimension, const int totalNum, 
		const int *I_accum, const int* numInRow, 
		const int* I, const int *J, const float *V, 
		const int totalNumPrecond, 
		const int *I_precond_accum, const int *numInRowL, 
		const int *I_precond, const int *J_precond, const float *V_precond, 
		const int totalNumPrecondP,
		const int *I_precondP_accum, const int *numInRowLP, 
		const int *I_precondP, const int *J_precondP, const float *V_precondP, 
		const float *y, float *x,  
		const int MAXIter, int *realIter, const bool RODR, 
		const int partition_size, const int* part_boundary);
}
#endif
