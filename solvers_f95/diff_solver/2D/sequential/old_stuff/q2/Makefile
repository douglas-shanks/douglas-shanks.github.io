LIB_LIST    = -lblas

F95 = gfortran
OPT = -Wall -fimplicit-none

all: testblas3 testlupp

testblas3.o: testblas3.f90
	$(F95) -c $(OPT)  testblas3.f90    

testlupp.o: testlupp.f90
	$(F95) -c $(OPT)  testlupp.f90    

blublas2.o: blublas2.f90
	$(F95) -c $(OPT) blublas2.f90     # simple lu (without pivoting)
                                         # of a subblock of a matrix
                                         # using blas level 2

luppblas2.o: luppblas2.f90
	$(F95) -c $(OPT) luppblas2.f90    # lu with partial pivoting
                                         # using blas level 2
lublas3.o: lublas3.f90
	$(F95) -c $(OPT) lublas3.f90      # simple lu (without pivoting)
                                         # using blas level 3

bsblas2.o: bsblas2.f90 
	$(F95) -c  $(OPT) bsblas2.f90     # bs using blas level 2

rsymm.o: rsymm.f90 
	$(F95) -c $(OPT) rsymm.f90

testblas3: testblas3.o blublas2.o lublas3.o bsblas2.o rsymm.o
	$(F95) -o testblas3 testblas3.o blublas2.o lublas3.o bsblas2.o rsymm.o $(LIB_LIST)

testlupp: testlupp.o luppblas2.o bsblas2.o
	$(F95) -o testlupp testlupp.o luppblas2.o bsblas2.o $(LIB_LIST)

# Note that the blas are only used in the linking step

clean: 
	rm -rf *.o  core*   # creates a command to remove old .o and core files




















