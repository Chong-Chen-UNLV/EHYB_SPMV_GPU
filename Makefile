################################################################################
#
## Build script for project
#
#################################################################################
CC=g++

CU_INC=-I/usr/local/cuda/include
CUDA_LIB =-L/usr/local/cuda/lib64 -lcublas -lcudart -lcuda

CFLAGS= -g -gdwarf-2 -O3 -funroll-loops -fopenmp -I./mt-metis-0.6.0/include -I./include -L./ -L./lib -lmtmetis #-lopenblas -std=gnu99 
	CUFLAGS = -O3 --use_fast_math -gencode=arch=compute_61,code=compute_61 
	CFILES          =       $(wildcard *.c)
	CUFILES         =       $(wildcard *.cu)
	OBJECTS         =       $(CFILES:.c=.o)
CU_OBJECTS      =       $(CUFILES:.cu=.o)

all : $(OBJECTS) $(CU_OBJECTS)
	$(CC) -m64 $^ $(CFLAGS) $(CU_INC) $(CUDA_LIB)  -o Solver_test.out

$(OBJECTS) : $(CFILES)
	$(CC) -m64 $^ $(CFLAGS) $(CU_INC) $(CUDA_LIB) -c

$(CU_OBJECTS) : $(CUFILES)
	nvcc -c -g -G $^ $(CUFLAGS)

clean :
	rm *.o


