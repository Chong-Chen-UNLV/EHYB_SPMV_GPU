#ifndef SOLVER_H
#define SOLVER_H

#include <omp.h>
void solverPrecondCPU(const int procNum, const int dimension, const int totalNum, const int *row_idx, const int *J, 
		const float *V, const int totalNumPrecond, const int *row_idxL, 
		const int *J_precond, const float *V_precond, const int totalNumPrecondP,
		const int *row_idxLP, const int *J_precondP, const float *V_precondP, 
		const float *vector_in, float *vector_out, const int MAXIter, int *realIter);

void solverGPU_HYB(const int dimension, 
		const int totalNum, const int* numInRow, 
		const int *row_idx,  const int* I, const int *J, const float *V, 
		const int totalNumPrecond, const int *numInRowL,
		const int *row_idxL,  
		const int *I_precond, const int *J_precond, const float *V_precond, 
		const int totalNumPrecondP, const int *numInRowLP,
		const int *row_idxLP,  
		const int *I_precondP, const int *J_precondP, const float *V_precondP, 
		const float *vector_in, float *vector_out,  
		const int MAXIter, int *realIter, const bool RODR, 
		const int partition_size, const int* part_boundary);

void solver(const int dimension, const int totalNum, const int *I, const int *J, const double *V, const double *vector_in, 
			double *vector_out, double *error_track, int MAXIter);

#endif
