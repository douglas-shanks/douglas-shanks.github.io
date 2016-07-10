LIB_LIST  =  $(BLASDIR)/lib/libf77blas.a $(BLASDIR)/lib/libatlas.a

F95   = mpif90
FLAGS =

OPT = -O

all: lsquares test_comm

lsquares.o: lsquares.f90
	gfortran $(FLAGS) $(OPT) -c lsquares.f90

lu1.o: lu1.f90
	gfortran $(FLAGS) $(OPT) -c lu1.f90

bs1.o: bs1.f90
	gfortran $(FLAGS) $(OPT) -c bs1.f90

nrmeq.o: nrmeq.f90
	gfortran $(FLAGS) $(OPT) -c nrmeq.f90

lsquares: lsquares.o lu1.o bs1.o nrmeq.o
	gfortran -o lsquares lsquares.o lu1.o bs1.o nrmeq.o $(LIB_LIST)

test_comm.o: test_comm.f90
	$(F95) $(FLAGS) $(OPT) -c test_comm.f90

test_comm: test_comm.o
	$(F95) $(OPT) -o test_comm test_comm.o

clean: 
	rm -rf *.o  core




















