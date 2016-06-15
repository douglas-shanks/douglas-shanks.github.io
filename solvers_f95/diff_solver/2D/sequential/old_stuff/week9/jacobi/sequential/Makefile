
F95 = gfortran
OPT = -O

all: iterative

header.o: header.f90
	$(F95) $(OPT) -c header.f90

laplace.o: laplace.f90 header.o
	$(F95) $(OPT) -c laplace.f90

matmult.o: matmult.f90 header.o
	$(F95) $(OPT) -c matmult.f90

jacobi.o: jacobi.f90 header.o
	$(F95) $(OPT) -c jacobi.f90

error.o: error.f90 header.o
	$(F95) $(OPT) -c error.f90

iterative.o: iterative.f90 header.o
	$(F95) $(OPT) -c iterative.f90

iterative: matmult.o laplace.o jacobi.o error.o iterative.o
	$(F95) $(OPT) -o iterative matmult.o laplace.o jacobi.o error.o iterative.o

clean: 
	rm -rf *~ *.o *.mod core  