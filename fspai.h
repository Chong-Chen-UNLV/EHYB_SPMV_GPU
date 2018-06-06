#ifndef FSPAI_H
#define FSPAI_H

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <sys/time.h>
#include <time.h>
#include <unistd.h>
#include <omp.h>

#define bzero(b,len) (memset((b), '\0', (len)), (void) 0)
#define padding 1

typedef struct _S
{
	unsigned int *I;
	unsigned int *J;
	float *V;
	unsigned int *I_precond;
	unsigned int *J_precond;
	float *V_precond;
	unsigned int maxRowNum;
	unsigned int *numInRow;	
	unsigned int *row_idx;
	unsigned int *numInRowPrecond;
	unsigned int *row_idxPrecond; 
	float *diag;
	unsigned int colStart;
	unsigned int colEnd;
	unsigned int id;
}S;

typedef struct _Sort_S{
	unsigned int idx;
	float val; 
}Sort_S;

inline int qs_compare(Sort_S *A, Sort_S *B);

void solver(const unsigned int dimension, const unsigned int totalNum, 
		const unsigned int *I, const unsigned int *J, 
		const float *V, float *tempCSR,
		const float *vector_in, float *vector_out,
		float *bp, float *pk, float *rk, int MAXIter);
			
void insertSort(unsigned int *J, float *V, unsigned int num, unsigned int *outJ, float *outV);

void solverPrecondCOO(const unsigned int dimension, const unsigned int totalNum,
			const unsigned int *I, const unsigned int *J, 
			const float *V, const unsigned int totalNumPrecond, 
			const unsigned int *I_precond,
			const unsigned int *J_precond, const float *V_precond, 
			const unsigned int totalNumPrecondP, const unsigned int *I_precondP, 
			const unsigned int *J_precondP, const float *V_precondP,
			const float *y, float *x, const unsigned int MAXIter, unsigned int *realIter);

void solverPrecondPhi(const unsigned int dimension, const unsigned int totalNum, 
			const unsigned int *I, const unsigned int *J, 
			const float *V, const unsigned int totalNumPrecond, 
			const unsigned int *I_precond,
			const unsigned int *J_precond, const float *V_precond, 
			const unsigned int totalNumPrecondP, const unsigned int *I_precondP, 
			const unsigned int *J_precondP, const float *V_precondP,
			const float *y, float *x, const int MAXIter, 
			int *realIter, int Rank);
				
void formatChange(unsigned int dimension, unsigned int *numInRow, unsigned int *totalNum_1, unsigned int *I, unsigned int *J, float *V, unsigned int **I_1, unsigned int **J_1, float **V_1);				

void solverCPU(const unsigned int dimension, const unsigned int totalNum, 
		const unsigned int *I, const unsigned int *J, 
		const float *V, const float *vector_in, 
		float *vector_out, float *error_track, int MAXIter, int *realIter);	
	
#endif
