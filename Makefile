################################################################################
#
# Build script for project
#
################################################################################
CC=icc
CFLAGS= -funroll-loops -openmp -std=c99 -O3 -par-affinity=scatter -mkl #-vec-report=3 
CFILES		=	$(wildcard *.c)
CUFILES		=	$(wildcard *.cu)
OBJECTS		=	$(CFILES:.c=.o)
CU_OBJECTS	= 	$(CUFILES:.cu=.o)

Phi1 : $(OBJECTS) $(CU_OBJECTS)
	$(CC) -m64 $^ $(CFLAGS) -o Solver_test.out 

$(OBJECTS) : $(CFILES)
	$(CC) -m64 $^ $(CFLAGS) -c

kernel.o : $(CUFILES)
	nvcc -c -g -G $^

clean :

	rm *.o


