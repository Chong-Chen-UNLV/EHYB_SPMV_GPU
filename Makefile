################################################################################
#
# Build script for project
#
################################################################################
CC=icc
CFLAGS= -funroll-loops -openmp -std=c99 -O3 -par-affinity=scatter -mkl #-vec-report=3 
CFILES		=	$(wildcard *.c)
OBJECTS		=	$(CFILES:.c=.o)

Phi1 : $(OBJECTS)
	$(CC) -m64 $^ $(CFLAGS) -o Solver_test.out 

%.o : %.c
	$(CC) -m64 $^ $(CFLAGS) -c


clean :

	rm *.o


