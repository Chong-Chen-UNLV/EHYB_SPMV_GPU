#ifndef REORDERING_H
#define REORDERING_H
#include <mtmetis.h>
void matrix_reorder(const int* dimension_in, const int totalNum, const int* I, const int* J, const float* V, 
		int* numInRow, int* I_rodr, int* J_rodr, float* V_rodr, int* rodr_list, int* part_boundary,
		const int rodr_option, const unsigned int nparts);
void update_numInRow(const int totalNum, const int dimension, int* I_rodr, 
		int* J_odr, float* V_odr, int* numInRow, int* numInRowL,
		int* row_idx, int* row_idxL, float* diag)
void vector_reorder(const int dimension, const int* rodr_list, const float* v_in, float* v_rodr);
void vector_recover(const int dimension, const int* rodr_list, const float* v_rodr, float* v);

#endif
