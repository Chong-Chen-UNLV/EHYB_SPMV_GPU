#ifndef SOLVER_H
#define SOLVER_H

#include <omp.h>
void solverPrecondCPU(const unsigned int procNum, const unsigned int dimension, 
		const unsigned int totalNum, const unsigned int *row_idx, const unsigned int *J, 
		const float *V, const unsigned int totalNumPrecond, const unsigned int *row_idxL, 
		const unsigned int *J_precond, const float *V_precond, const unsigned int totalNumPrecondP,
		const unsigned int *row_idxLP, const unsigned int *J_precondP, const float *V_precondP, 
		const float *vector_in, float *vector_out, const unsigned int MAXIter, unsigned int *realIter);

void solverGPU_HYB(const unsigned int dimension, 
		const unsigned int totalNum, const unsigned int* numInRow, const unsigned int maxRowNum, 
		const unsigned int *row_idx,  const unsigned int* I, const unsigned int *J, const float *V, 
		const unsigned int totalNumPrecond, const unsigned int *numInRowL, const unsigned int maxRowNumPrecond,
		const unsigned int *row_idxL,  
		const unsigned int *I_precond, const unsigned int *J_precond, const float *V_precond, 
		const unsigned int totalNumPrecondP, const unsigned int *numInRowLP, const unsigned int maxRowNumPrecondP,
		const unsigned int *row_idxLP,  
		const unsigned int *I_precondP, const unsigned int *J_precondP, const float *V_precondP, 
		const float *vector_in, float *vector_out,  
		const unsigned int MAXIter, unsigned int *realIter, const bool RODR, const bool BLOCK,
		const unsigned int partition_size, const unsigned int* part_boundary);

#endif
