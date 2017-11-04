#ifndef FSPAI_H
#define FSPAI_H

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <sys/time.h>
#include <time.h>
#include <unistd.h>
#include <omp.h>


#define padding 1

typedef struct _S
{
	int *I;
	int *J;
	double *V;
	int *I_precond;
	int *J_precond;
	double *V_precond;
	int maxRowNum;
	int *numInRow;	
	int *rowNumAccum;
	int *numInRowPrecond;
	int *rowNumAccumPrecond; 
	double *diag;
	int colStart;
	int colEnd;
	int id;
}S;

void solver(const int dimension, const int totalNum, const int *I, const int *J, const double *V, double *tempCSR,
			const double *vector_in,double *vector_out, double *bp, double *pk, double *rk, int MAXIter);
			
void insertSort(int *J, double *V, int num, int *outJ, double *outV);

void solverPrecondCOO(const int dimension, const int totalNum, const int *I, const int *J, const double *V, const int totalNumPrecond, const int *I_precond,
				const int *J_precond, const double *V_precond, const int totalNumPrecondP, const int *I_precondP, const int *J_precondP, const double *V_precondP,
				const double *y, double *x, const int MAXIter, int *realIter);

void solverPrecondPhi(const int dimension, const int totalNum, const int *I, const int *J, const double *V, const int totalNumPrecond, const int *I_precond,
				const int *J_precond, const double *V_precond, const int totalNumPrecondP, const int *I_precondP, const int *J_precondP, const double *V_precondP,
				const double *y, double *x, const int MAXIter, int *realIter, int rank);
				
void formatChange(int dimension, int *numInRow, int *totalNum_1, int *I, int *J, double *V, int **I_1, int **J_1, double **V_1);				

void solverCPU(const int dimension, const int totalNum, const int *I, const int *J, const float *V, const float *vector_in, 
			float *vector_out, float *error_track, int MAXIter, int *realIter);	
	
#endif
