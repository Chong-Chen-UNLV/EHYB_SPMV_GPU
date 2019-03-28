#include "solver.h"

void matrix_vectorTest(matrixCOO_S *matrix, const double* vector_in, double* vector_out){
	for(unsigned int i = 0; i < matrix->dimension; ++i) vector_out[i] = 0;
	for(unsigned int i = 0; i < matrix->totalNum; ++i){
		vector_out[matrix->I[i]] += matrix->V[i]*vector_in[matrix->J[i]];	
	}
}
