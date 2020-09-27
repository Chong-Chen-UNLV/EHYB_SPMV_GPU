#ifndef SOLVER_H
#define SOLVER_H

#include <omp.h>
#include <stdint.h>

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
    int* numInRow2;//"real" numInRow for vetex goto blockELL
    int* I;
    int* J;
    double* V;
    double* diag;
	int* reorderList;
}matrixCOO;

typedef struct _matrixEHYB{
	int dimension;
	int	nParts;
	int	numOfRowER;
	int* reorderList;
	int* reorderListER;
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

void solverGPuUprecondEHYB(matrixEHYB* localMatrix, 
		const double* vectorIn, double* vectorOut,  
		const int MAXIter, int* realIter);

void solverGPuUnprecondCUSPARSE(matrixCOO* localMatrix, 
		const double *vector_in, double *vector_out,  
		const int MAXIter, int *realIter,  const cb_s cb,
		const int partSize, const int* partBoundary);

#endif
