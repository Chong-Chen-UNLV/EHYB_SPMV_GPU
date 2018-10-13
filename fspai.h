#ifndef FSPAI_H
#define FSPAI_H

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <sys/time.h>
#include <time.h>
#include <unistd.h>
#include <omp.h>

//#define bzero(b,len) (memset((b), '\0', (len)), (void) 0)
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
	double val; 
}Sort_S;

inline int qs_compare(Sort_S *A, Sort_S *B);

void solver(const unsigned int dimension, const unsigned int totalNum, 
		const unsigned int *I, const unsigned int *J, 
		const double *V, double *tempCSR,
		const double *vector_in, double *vector_out,
		double *bp, double *pk, double *rk, int MAXIter);
			
void insertSort(unsigned int *J, double *V, unsigned int num, unsigned int *outJ, double *outV);

/*void solverPrecondCOO(const unsigned int dimension, const unsigned int totalNum,
			const unsigned int *I, const unsigned int *J, 
			const double *V, const unsigned int totalNumPrecond, 
			const unsigned int *I_precond,
			const unsigned int *J_precond, const double *V_precond, 
			const unsigned int totalNumPrecondP, const unsigned int *I_precondP, 
			const unsigned int *J_precondP, const double *V_precondP,
			const double *y, double *x, const unsigned int MAXIter, unsigned int *realIter);

void solverPrecondPhi(const unsigned int dimension, const unsigned int totalNum, 
			const unsigned int *I, const unsigned int *J, 
			const double *V, const unsigned int totalNumPrecond, 
			const unsigned int *I_precond,
			const unsigned int *J_precond, const double *V_precond, 
			const unsigned int totalNumPrecondP, const unsigned int *I_precondP, 
			const unsigned int *J_precondP, const double *V_precondP,
			const double *y, double *x, const int MAXIter, 
			int *realIter, int Rank);*/
				
void formatChange(unsigned int dimension, unsigned int *numInRow, unsigned int *totalNum_1, unsigned int *I, unsigned int *J, double *V, unsigned int **I_1, unsigned int **J_1, double **V_1);				

void solverCPU(const unsigned int dimension, const unsigned int totalNum, 
		const unsigned int *I, const unsigned int *J, 
		const double *V, const double *vector_in, 
		double *vector_out, double *error_track, int MAXIter, int *realIter);	
	
#endif
