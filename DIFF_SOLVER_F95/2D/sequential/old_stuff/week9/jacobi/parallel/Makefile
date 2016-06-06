
F95       = mpif90
FLAGS     =
LIB_LIST  =

OPT = -O

all: iterative

greetings.o: greetings.f90
	$(F95) $(OPT) $(FLAGS) -c greetings.f90

header.o: header.f90
	$(F95) $(OPT) -c header.f90

laplace.o: laplace.f90 header.o
	$(F95) $(OPT) $(FLAGS) -c laplace.f90

matmult.o: matmult.f90 header.o
	$(F95) $(OPT) $(FLAGS) -c matmult.f90

iterative.o: iterative.f90 header.o
	$(F95) $(OPT) $(FLAGS) -c iterative.f90

iterative: greetings.o matmult.o laplace.o iterative.o
	 $(F95) $(OPT) $(FLAGS) -o iterative greetings.o matmult.o laplace.o iterative.o $(LIB_LIST)

clean: 
	rm -rf *~ *.o *.mod core  
