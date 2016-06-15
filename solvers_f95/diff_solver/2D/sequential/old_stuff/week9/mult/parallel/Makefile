
F95       = mpif90
FLAGS     =
LIB_LIST  =

OPT = -O

all: test_mult

greetings.o: greetings.f90
	$(F95) $(OPT) $(FLAGS) -c greetings.f90

header.o: header.f90
	$(F95) $(OPT) -c header.f90

laplace.o: laplace.f90 header.o
	$(F95) $(OPT) $(FLAGS) -c laplace.f90

matmult.o: matmult.f90 header.o
	$(F95) $(OPT) $(FLAGS) -c matmult.f90

test_mult.o: test_mult.f90 header.o
	$(F95) $(OPT) $(FLAGS) -c test_mult.f90

test_mult: greetings.o matmult.o laplace.o test_mult.o
	 $(F95) $(OPT) $(FLAGS) -o test_mult greetings.o matmult.o laplace.o test_mult.o $(LIB_LIST)

clean: 
	rm -rf *~ *.o *.mod core  
