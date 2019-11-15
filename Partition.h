#ifndef PARTITION_H
#define PARTITION_H

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <sys/time.h>
#include <time.h>
#include <unistd.h>
#include <mtmetis.h>

typedef struct _rowS{
	unsigned int idx;
	unsigned int nonzeros;
}rowS;

inline int rowSCompare(const void *A, const void *B){
	rowS *A_ = (rowS*) A;
	rowS *B_ = (rowS*) B;
	if(A_->nonzeros > B_->nonzeros) return -1;
	if(A_->nonzeros == B_->nonzeros) return 0;
	if(A_->nonzeros < B_->nonzeros) return 1;
}

#endif
