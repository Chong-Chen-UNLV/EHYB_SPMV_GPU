#ifndef REORDERING_H
#define REORDERING_H
#include <mtmetis.h>

void matrix_reorder(const unsigned int* dimension_in, const unsigned int totalNum, const unsigned int* I,
		const unsigned int* J, const double* V, 
		unsigned int* numInRow, unsigned int* row_idx, 
		unsigned int* I_rodr, unsigned int* J_rodr, 
		double* V_rodr, unsigned int* rodr_list, unsigned int* part_boundary,
		const unsigned int nparts);

void update_numInRowL(const unsigned int totalNum, 
			const unsigned int dimension, 
			unsigned int* I_rodr, 
			unsigned int* J_rodr, 
			double* V_rodr, 
			unsigned int* numInRowL,
			unsigned int* maxL,
			unsigned int* maxLP,
			unsigned int* row_idxL, 
			unsigned int* row_idxLP, 
			double* diag);

void vector_reorder(const unsigned int dimension, const double* v_in, double* v_rodr, const unsigned int* rodr_list);
void vector_recover(const unsigned int dimension, const double* v_rodr, double* v, const unsigned int* rodr_list);

#endif
