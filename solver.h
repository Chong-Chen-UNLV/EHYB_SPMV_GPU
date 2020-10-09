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
    int maxCol;
	int nParts;
    int* rowIdx;                 
    int* numInRow;
    int* numInRow2;//"real" numInRow for vetex goto blockELL
    int* I;
    int* J;
    double* V;
    double* diag;
	int* partBoundary;
	int* reorderList;
}matrixCOO;

typedef struct _matrixEHYB{
	int dimension;
	int	nParts;
	int	numOfRowER;
	int* reorderList;
	int* reorderListER;
	int* widthVecBlockELL;
	int* biasVecBlockELL;
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


void solverGPuUnprecondEHYB(matrixCOO* localMatrix, 
		const double* vectorIn, double* vectorOut,  
		const int MAXIter, int* realIter);

void solverGPuUnprecondCUSPARSE(matrixCOO* localMatrix, 
		const double *vector_in, double *vector_out,  
		const int MAXIter, int *realIter);

#endif
