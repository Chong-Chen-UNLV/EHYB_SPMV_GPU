#ifndef SPMV_H
#define SPMV_H

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
	uint16_t vectorCacheSize;
	int16_t kernelPerPart;
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
	int16_t vectorCacheSize;
	int	kernelPerPart;
	int	numOfRowER;
	int* warpIdxER_d;
	//int	threadSizeER;
	//int	blockSizeER;
	int* reorderList;
	int* reorderListER;
	int16_t* widthVecBlockELL;
	int* biasVecBlockELL;
	int16_t* colBlockELL;
	double* valBlockELL;
	int* partBoundary;
	int16_t* widthVecER;
	int* rowVecER;
	int* biasVecER; 
	int* colER;
	double* valER;
	double* outER;
	//the long vectors 
	int nLongVec = 0;
	int* longVecBoundary;
	int* longVecRow;
	int* longVecCol;
	double* longVecVal;
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

extern "C"
void spmvGPuEHYB(matrixCOO* localMatrix, 
		const double* vectorIn, double* vectorOut,  
		const int MAXIter, int* realIter);

void solverGPuUnprecondCUSPARSE(matrixCOO* localMatrix, 
		const double *vector_in, double *vector_out,  
		const int MAXIter);

void spmvGeneric(matrixCOO* localMatrix, 
		const double *vector_in, double *vector_out,  
		const int MAXIter);
#endif

