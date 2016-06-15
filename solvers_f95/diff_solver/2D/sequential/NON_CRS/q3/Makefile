LIB_LIST = -llapack -lblas

F95 = gfortran
#OPT = -Wall -C
OPT = -O

all: testsolvers

laplace.o: laplace.f90
	$(F95) $(OPT) -c laplace.f90

rhs.o: rhs.f90
	$(F95) $(OPT) -c rhs.f90

testsolvers.o: testsolvers.f90
	$(F95) $(OPT) -c testsolvers.f90   

testsolvers: testsolvers.o laplace.o rhs.o 
	$(F95) -o testsolvers $(OPT) testsolvers.o laplace.o rhs.o $(LIB_LIST)

# Note that the lapack library is only used in the linking step

clean: 
	rm -rf *.o  core*   # creates a command to remove old .o and core files




















