#ifndef REORDERING_H
#define REORDERING_H
#include "solver.h"


void matrixReorder(matrixCOO* localMatrixCOO);

void vectorReorder(const int dimension, const double* v_in, double* v_rodr, const int* rodr_list);
void vectorRecover(const int dimension, const double* v_rodr, double* v, const int* rodr_list);

#endif
