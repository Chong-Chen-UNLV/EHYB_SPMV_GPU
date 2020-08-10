################################################################################
#
## Build script for project
#
#################################################################################
CC=g++

CU_INC=-I/opt/packages/cuda/10.1/include
CUDA_LIB =-L/opt/packages/cuda/10.1/lib64 -lcublas -lcudart 

CFLAGS= -g -gdwarf-2 -O0 -funroll-loops -fopenmp -I./mt-metis-0.6.0/include -I./include -L./ -L./lib -lmtmetis #-lopenblas -std=gnu99 
	CUFLAGS = -O0 --use_fast_math -gencode=arch=compute_60,code=compute_60 
	CFILES          =       $(wildcard *.c)
	CUFILES         =       $(wildcard *.cu)
	OBJECTS         =       $(CFILES:.c=.o)
CU_OBJECTS      =       $(CUFILES:.cu=.o)

all : $(OBJECTS) $(CU_OBJECTS)
	$(CC) -m64 $^ $(CFLAGS) $(CU_INC) $(CUDA_LIB)  -o Solver_test.out

$(OBJECTS) : $(CFILES) *.h
	$(CC) -m64 $(CFILES) $(CFLAGS) $(CU_INC) $(CUDA_LIB) -c

$(CU_OBJECTS) : $(CUFILES) *.h
	nvcc -c -g -G $(CUFILES) $(CUFLAGS)

clean :
	rm *.o


