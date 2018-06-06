#ifndef REORDERING_H
#define REORDERING_H
#include <mtmetis.h>

void matrix_reorder(const unsigned int* dimension_in, const unsigned int totalNum, const unsigned int* I,
		const unsigned int* J, const float* V, 
		unsigned int* numInRow, unsigned int* row_idx, 
		unsigned int* I_rodr, unsigned int* J_rodr, 
		float* V_rodr, unsigned int* rodr_list, unsigned int* part_boundary,
		const unsigned int nparts);

void update_numInRowL(const unsigned int totalNum, 
			const unsigned int dimension, 
			unsigned int* I_rodr, 
			unsigned int* J_rodr, 
			float* V_rodr, 
			unsigned int* numInRowL,
			unsigned int* row_idxL, 
			float* diag);

void vector_reorder(const unsigned int dimension, const float* v_in, float* v_rodr, const unsigned int* rodr_list);
void vector_recover(const unsigned int dimension, const float* v_rodr, float* v, const unsigned int* rodr_list);

#endif
