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
	int *I;
	int *J;
	float *V;
	int *I_precond;
	int *J_precond;
	float *V_precond;
	int maxRowNum;
	int *numInRow;	
	int *row_idx;
	int *numInRowPrecond;
	int *row_idxPrecond; 
	float *diag;
	int colStart;
	int colEnd;
	int id;
}S;

void solver(const int dimension, const int totalNum, const int *I, const int *J, const float *V, float *tempCSR,
			const float *vector_in,float *vector_out, float *bp, float *pk, float *rk, int MAXIter);
			
void insertSort(int *J, float *V, int num, int *outJ, float *outV);

void solverPrecondCOO(const int dimension, const int totalNum, const int *I, const int *J, const float *V, const int totalNumPrecond, const int *I_precond,
				const int *J_precond, const float *V_precond, const int totalNumPrecondP, const int *I_precondP, const int *J_precondP, const float *V_precondP,
				const float *y, float *x, const int MAXIter, int *realIter);

void solverPrecondPhi(const int dimension, const int totalNum, const int *I, const int *J, const float *V, const int totalNumPrecond, const int *I_precond,
				const int *J_precond, const float *V_precond, const int totalNumPrecondP, const int *I_precondP, const int *J_precondP, const float *V_precondP,
				const float *y, float *x, const int MAXIter, int *realIter, int rank);
				
void formatChange(int dimension, int *numInRow, int *totalNum_1, int *I, int *J, float *V, int **I_1, int **J_1, float **V_1);				

void solverCPU(const int dimension, const int totalNum, const int *I, const int *J, const float *V, const float *vector_in, 
			float *vector_out, float *error_track, int MAXIter, int *realIter);	
	
#endif
