#ifndef REORDERING_H
#define REORDERING_H
#include <mtmetis.h>
#include "solver.h"


void matrixReorder(const unsigned int* dimension_in, const unsigned int totalNum, 
		const cb_s cb,
		const unsigned int* I, 
		const unsigned int* J, const double* V, 
		unsigned int* numInRow, unsigned int* row_idx, 
		unsigned int* I_rodr, unsigned int* J_rodr, 
		double* V_rodr, unsigned int* rodr_list, unsigned int* part_boundary,
		const unsigned int nparts);


void vectorReorder(const unsigned int dimension, const double* v_in, double* v_rodr, const unsigned int* rodr_list);
void vectorRecover(const unsigned int dimension, const double* v_rodr, double* v, const unsigned int* rodr_list);

#endif
