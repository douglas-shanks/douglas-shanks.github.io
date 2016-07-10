
F95       = mpif90 
FLAGS     =
LIB_LIST  =

all: test_steady

greetings.o: greetings.f90
	$(F95) $(OPT) $(FLAGS) -c greetings.f90

header.o: header.f90
	$(F95) $(OPT) -c header.f90

sparsegather.o: sparsegather.f90 header.o
	$(F95) $(OPT) $(FLAGS) -c sparsegather.f90

sysmatrix.o: sysmatrix.f90 header.o
	$(F95) $(OPT) -c sysmatrix.f90

rhs.o: rhs.f90 header.o
	$(F95) $(OPT) -c rhs.f90
		
alpha.o: alpha.f90 header.o
	$(F95) $(OPT) -c alpha.f90	
	
error.o: error.f90 header.o
	$(F95) $(OPT) $(FLAGS) -c error.f90

matmult.o: matmult.f90 header.o
	$(F95) $(OPT) $(FLAGS) -c matmult.f90

jacobi.o: jacobi.f90 header.o
	$(F95) $(OPT) $(FLAGS) -c jacobi.f90

test_steady.o: test_steady.f90 header.o
	$(F95) $(OPT) $(FLAGS) -c test_steady.f90

test_steady: greetings.o sparsegather.o header.o error.o matmult.o sysmatrix.o rhs.o jacobi.o test_steady.o alpha.o
	 $(F95) $(OPT) -o test_steady greetings.o header.o sparsegather.o error.o matmult.o sysmatrix.o rhs.o jacobi.o test_steady.o alpha.o$(LIB_LIST)

clean: 
	rm -rf *~ *.o *.mod core  
