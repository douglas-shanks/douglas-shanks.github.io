LIB_LIST = -llapack -lblas

# Above is the location of the blas library which is called by lublas2
# The filename with full path for this library is /usr/lib/libblas.a
# Since it is in /usr/lib we do not need to specify a path.
# The compiler looks there automatically.
F95 = gfortran
OPT = -O

all:test2D

exactsoln2D.o: exactsoln2D.f90
	$(F95) -c exactsoln2D.f90	
	
test2D.o: test2D.f90
	$(F95) -c test2D.f90
	
FTCSmatrix.o: FTCSmatrix.f90
	$(F95) -c FTCSmatrix.f90
	
BTCSmatrix.o: BTCSmatrix.f90
	$(F95) -c BTCSmatrix.f90

CNmatrix.o: CNmatrix.f90
	$(F95) -c CNmatrix.f90	

CNrhs.o: CNrhs.f90
	$(F95) -c CNrhs.f90
	
norm.o: norm.f90
	$(F95) -c norm.f90	

test2D: test2D.o exactsoln2D.o FTCSmatrix.o BTCSmatrix.o CNmatrix.o CNrhs.o norm.o
	 $(F95) -o test2D test2D.o exactsoln2D.o FTCSmatrix.o BTCSmatrix.o CNmatrix.o CNrhs.o norm.o  $(LIB_LIST)

# Note that the blas are only used in the linking step
clean: 
	rm -rf *~ *.o *.mod core  
