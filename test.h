#include "solver.h"

void matrix_vectorTest(matrixCOO_S *matrix, const double* vector_in, double* vector_out, 
		unsigned int testpoint){
	FILE *tFile;
	if(testpoint > 0) {
		if ((tFile = fopen("testResult", "w")) == NULL){ 
			printf("file open error\n");
			exit(1);
		}
	}
	for(unsigned int i = 0; i < matrix->dimension; ++i) vector_out[i] = 0;
	for(unsigned int i = 0; i < matrix->totalNum; ++i){
		vector_out[matrix->I[i]] += matrix->V[i]*vector_in[matrix->J[i]];	
		
		if(matrix->I[i] == testpoint && testpoint > 0){ fprintf(tFile, 
				"row is %d V is %f vec is %f accum is %f\n", 
				testpoint, matrix->V[i], vector_in[matrix->J[i]], vector_out[matrix->I[i]]);
		}
	}
	if(testpoint > 0)
		fclose(tFile);
}
