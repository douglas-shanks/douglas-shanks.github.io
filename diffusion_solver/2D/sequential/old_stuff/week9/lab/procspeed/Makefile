LIB_LIST  = $(BLASDIR)/lib/libf77blas.a $(BLASDIR)/lib/libatlas.a 
# Above is the location on aquila of the BLAS library compatible with the 
# ifort compiler. This is called by blublas2 and lublas3

F95 = gfortran
OPT = -u -O

all: testblas3

testblas3.o: testblas3.f90
	$(F95) -c $(OPT)  testblas3.f90    

blublas2.o: blublas2.f90
	$(F95) -c $(OPT) blublas2.f90

lublas3.o: lublas3.f90
	$(F95) -c $(OPT) lublas3.f90

bsblas2.o: bsblas2.f90 
	$(F95) -c  $(OPT) bsblas2.f90

rsymm.o: rsymm.f90 
	$(F95) -c $(OPT) rsymm.f90

testblas3: testblas3.o blublas2.o lublas3.o bsblas2.o rsymm.o
	$(F95) -o testblas3 testblas3.o blublas2.o lublas3.o bsblas2.o rsymm.o $(LIB_LIST)

clean: 
	rm -rf *.o  core*



















