#ifndef SOLVER_H
#define SOLVER_H

#include <omp.h>
#include "convert.h"

typedef struct _cb{
	bool GPU;
	bool RODR;
	bool CACHE;
	bool BLOCK;
}cb_s;

typedef struct _matrixCOO_S{
    unsigned int totalNum;
    unsigned int dimension;                 
    unsigned int maxRowNum;
    unsigned int* row_idx;                 
    unsigned int* numInRow;
    unsigned int* I;
    unsigned int* J;
    double* V;
}matrixCOO_S;

inline void init_cb(cb_s* in_s)
{
    in_s->GPU = false;
    in_s->RODR = false; 
    in_s->BLOCK = false;
    in_s->CACHE = false;
}

inline void init_matrixCOO_S(matrixCOO_S* matrix, unsigned int dimension,
        unsigned int totalNum, unsigned int maxRowNum, unsigned int* row_idx, 
        unsigned int* numInRow, unsigned int* I,  unsigned int* J, double* V ){
    matrix->totalNum = totalNum;
    matrix->dimension = dimension;
    matrix->maxRowNum = maxRowNum;
    matrix->row_idx = row_idx;
    matrix->numInRow = numInRow;
    matrix->I = I;
    matrix->J = J;
    matrix->V = V;
}

void solverPrecondCPU(const unsigned int procNum, const unsigned int dimension, 
		const unsigned int totalNum, const unsigned int *row_idx, const unsigned int *J, 
		const double *V, const unsigned int totalNumPrecond, const unsigned int *row_idxL, 
		const unsigned int *J_precond, const double *V_precond, const unsigned int totalNumPrecondP,
		const unsigned int *row_idxLP, const unsigned int *J_precondP, const double *V_precondP, 
		const double *vector_in, double *vector_out, const unsigned int MAXIter, unsigned int *realIter);

void solverGPU_HYB(matrixCOO_S* localMatrix, matrixCOO_S* localMatrixPrecond, 
                matrixCOO_S* localMatrixPrecondP,
        		const double *vector_in, double *vector_out,  
                const unsigned int MAXIter, unsigned int *realIter,  const cb_s cb,
                const unsigned int partition_size, const unsigned int* part_boundary);

#endif
