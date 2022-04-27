################################################################################
#
## Build script for cuda projct
#
#################################################################################
CC=g++
CUDA_HOME=/home/chong/anaconda3/envs/cuda
CU_INC=-I $(CUDA_HOME)/include
CUDA_LIB =-L $(CUDA_HOME)/lib64 -lcublas -lcusparse -lcudart 

CFLAGS= -g -w -gdwarf-2 -O3 -funroll-loops -fopenmp -L./ -L./lib -lmtmetis #-lopenblas -std=gnu99 
CUFLAGS = -O3 -arch=sm_86 
CFILES = convert.c mmio.c reordering.c solver_test.c    
CUFILES = $(wildcard *.cu)
OBJECTS = $(CFILES:.c=.o)
CU_OBJECTS = $(CUFILES:.cu=.o)

all : $(OBJECTS) $(CU_OBJECTS)
	$(CC) -m64 $^ $(CFLAGS) $(CU_INC) $(CUDA_LIB)  -o spmvAlg1.out

$(OBJECTS) : $(CFILES) *.h 
	$(CC) -m64 $(CFILES) $(CFLAGS) $(CU_INC) $(CUDA_LIB) -c

$(CU_OBJECTS) : $(CUFILES) *.h 
	nvcc -g -c $(CUFILES) $(CUFLAGS)

clean :
	rm *.o


