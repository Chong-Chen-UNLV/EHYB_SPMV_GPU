#ifndef SOLVER_H
#define SOLVER_H

void solverPrecondCPU(const int dimension, const int totalNum, const int *I_accum, const int *J, 
		const float *V, const int totalNumPrecond, const int *I_precond_accum, 
		const int *J_precond, const float *V_precond, const int totalNumPrecondP,
		const int *I_precondP_accum, const int *J_precondP, const float *V_precondP, 
		const float *y, float *x, const int MAXIter, int *realIter);

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

void solver(const int dimension, const int totalNum, const int *I, const int *J, const double *V, const double *vector_in, 
			double *vector_out, double *error_track, int MAXIter);

#endif
