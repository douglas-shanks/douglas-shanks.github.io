LIB_LIST = -llapack -lblas
F95 = gfortran
OPT = -Wall -fimplicit-none

all: testlu

rhs.o: rhs.f90
	$(F95) $(OPT) -c rhs.f90

laplace.o: laplace.f90
	$(F95) $(OPT) -c laplace.f90

lu1.o: lu1.f90
	$(F95) $(OPT) -c lu1.f90

bs1.o: bs1.f90
	$(F95) $(OPT) -c bs1.f90

testlu.o: testlu.f90 
	$(F95)  $(OPT) -c testlu.f90   

testlu: testlu.o rhs.o laplace.o lu1.o bs1.o
	$(F95) $(OPT) -o testlu rhs.o laplace.o testlu.o  lu1.o  bs1.o $(LIB_LIST)

clean: 
	rm -rf *.o  core* 



