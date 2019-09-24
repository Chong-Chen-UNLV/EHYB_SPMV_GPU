#ifndef PARTITION_H
#define PARTITION_H

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <sys/time.h>
#include <time.h>
#include <unistd.h>
#include <mtmetis.h>

typedef _rowS{
	unsigned int idx;
	unsigned int nonzeros;
}rowS;

inline int row_compare(const void *A, const void *B){
	Sort_S *A_ = (Sort_S *) A;
	Sort_S *B_ = (Sort_S *) B;
	if(A_->idx > B_->idx) return 1;
	if(A_->idx == B_->idx) return 0;
	if(A_->idx < B_->idx) return -1;
}

#endif
