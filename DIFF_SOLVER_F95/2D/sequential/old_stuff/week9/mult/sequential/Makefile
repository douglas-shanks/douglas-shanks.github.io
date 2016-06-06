
F95 = gfortran
OPT = -O

all: test_mult

header.o: header.f90
	$(F95) $(OPT) -c header.f90

laplace.o: laplace.f90 header.o
	$(F95) $(OPT) -c laplace.f90

matmult.o: matmult.f90 header.o
	$(F95) $(OPT) -c matmult.f90

test_mult.o: test_mult.f90 header.o
	$(F95) $(OPT) -c test_mult.f90

test_mult: matmult.o laplace.o test_mult.o
	$(F95) $(OPT) -o test_mult matmult.o laplace.o test_mult.o

clean:
	rm -rf *~ *.o *.mod core
