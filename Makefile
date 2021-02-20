################################################################################
#
## Build script for projct
#
#################################################################################
CC=g++

CU_INC=-I/cm/shared/apps/cuda10.2/toolkit/10.2.89/include
CUDA_LIB =-L/cm/shared/apps/cuda10.2/toolkit/10.2.89/lib64 -lcublas -lcusparse -lcudart 

CFLAGS= -g -Wall -gdwarf-2 -O3 -funroll-loops -fopenmp -I./mt-metis-0.6.0/include -I./include -L./ -L./lib -lmtmetis #-lopenblas -std=gnu99 
	CUFLAGS = -O3 -arch=sm_70 
	CFILES = convert.c mmio.c reordering.c solver_test.c      
	CUFILES = $(wildcard *.cu)
	OBJECTS = $(CFILES:.c=.o)
CU_OBJECTS = $(CUFILES:.cu=.o)

all : $(OBJECTS) $(CU_OBJECTS)
	$(CC) -m64 $^ $(CFLAGS) $(CU_INC) $(CUDA_LIB)  -o spmv.out

$(OBJECTS) : $(CFILES) *.h
	$(CC) -m64 $(CFILES) $(CFLAGS) $(CU_INC) $(CUDA_LIB) -c

$(CU_OBJECTS) : $(CUFILES) *.h
	nvcc -c $(CUFILES) $(CUFLAGS)

clean :
	rm *.o


