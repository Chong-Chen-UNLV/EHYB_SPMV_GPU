#ifndef REORDERING_H
#define REORDERING_H
#include "solver.h"


void matrixReorder(matrixCOO* localMatrixCOO);

void vectorReorder(const int dimension, const float* v_in, float* v_rodr, const int* rodr_list);
void vectorRecover(const int dimension, const float* v_rodr, float* v, const int* rodr_list);

#endif
