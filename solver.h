#ifndef SOLVER_H
#define SOLVER_H

#include <omp.h>
#include <stdint.h>

struct abc{ int a; int b; int c;};

typedef struct _cb{
	bool PRECOND;
	bool GPU;
	bool RODR;
	bool CACHE;
	bool BLOCK;
	bool FACT;
	bool SORT;
}cb_s;

typedef struct _matrixCOO{
    int totalNum;
    int dimension;                 
    int maxRowNum;
	int nparts;
    int* rowIdx;                 
    int* numInRow;
    int* numInRow2;//used for ER part for EHYB
    int* I;
    int* J;
    double* V;
	int* diag;
}matrixCOO;

typedef struct _matrixEHYB{
	int dimension;
	int	 blocks;
	int	 numOfRowER;
	int* widthVecBlockELL;
	int* biasVecBLockELL;
	int* colBlockELL;
	double* valBlockELL;
	int* partBoundary;
	int* widthVecER;
	int* rowVecER;
	int* biasVecER; 
	int* colER;
	double* valER;
}matrixEHYB;

inline void init_cb(cb_s* in_s)
{
	in_s->PRECOND = false;
    in_s->GPU = false;
    in_s->RODR = true; 
    in_s->BLOCK = true;
    in_s->CACHE = true;
	in_s->FACT = true;
}

inline void init_matrixCOO_S(matrixCOO_S* matrix, uint32_t dimension,
        uint32_t totalNum, uint32_t maxRowNum, uint32_t* row_idx, 
        uint32_t* numInRow, uint32_t* I,  uint32_t* J, double* V ){
    matrix->totalNum = totalNum;
    matrix->dimension = dimension;
    matrix->maxRowNum = maxRowNum;
    matrix->row_idx = row_idx;
    matrix->numInRow = numInRow;
    matrix->I = I;
    matrix->J = J;
    matrix->V = V;
}

void solverPrecondCPU(const uint32_t procNum, const uint32_t dimension, 
		const uint32_t totalNum, const uint32_t *row_idx, const uint32_t *J, 
		const double *V, const uint32_t totalNumPrecond, const uint32_t *row_idxL, 
		const uint32_t *J_precond, const double *V_precond, const uint32_t totalNumPrecondP,
		const uint32_t *row_idxLP, const uint32_t *J_precondP, const double *V_precondP, 
		const double *vector_in, double *vector_out, const uint32_t MAXIter, uint32_t *realIter);

void solverGPU_UprecondEHYB(matrixCOO_S* localMatrix, 
		const double *vector_in, double *vector_out,  
		const uint32_t MAXIter, uint32_t *realIter,  const cb_s cb,
		const uint32_t part_size, const uint32_t* part_boundary);

void solverGPU_EHYB(matrixCOO_S* localMatrix, matrixCOO_S* localMatrixPrecond, 
                matrixCOO_S* localMatrixPrecondP,
        		const double *vector_in, double *vector_out,  
                const uint32_t MAXIter, uint32_t *realIter,  const cb_s cb,
                const uint32_t partition_size, const uint32_t* part_boundary);

#endif
