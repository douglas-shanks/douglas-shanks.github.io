program test_comm

 include "mpif.h"

 integer :: imin,imax
 parameter (imin = 1, imax = 7)
 integer :: myid, numprocs, ierr, rc
 real (kind=8), dimension(:), allocatable :: buffer    
 real (kind=8) :: t1,t2
 real (kind=8), dimension(imax) :: time  
 integer :: i,stat(MPI_STATUS_SIZE),dest
 integer, dimension(imax) :: n

 call MPI_INIT(ierr)
 call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)
 call MPI_COMM_SIZE(MPI_COMM_WORLD, numprocs, ierr)

 allocate(buffer(2**imax))

 if (myid == 0) then

    call random_number(buffer)

    do i=imin,imax
       
       n(i) = 2**i 

       t1 = MPI_WTIME()

       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       !!!  For each value of n = 2**i send the first n     !!!
       !!!  entries of buffer (i.e. an array of length n)   !!!
       !!!  50 times to process 1 and back and measure the  !!!
       !!!  elapsed time. Use this to calculate the time    !!!
       !!!  for one send.                                   !!! 
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


       !!!   ADD NECCESSARY SENDS AND RECEIVES HERE  !!!


       t2 = MPI_WTIME() 

       ! time(i) = ??
           
    end do

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!  Now write the data (i.e. the value of n and the    !!!
    !!!  time for one send) to a file lsdata.txt in the     !!!
    !!!  format needed for lsquares.f90                     !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    open(unit=2,file='lsdata.txt')
    write(2,*) imax-imin+1
    write(2,*) n(imin:imax)
    write(2,*) time(imin:imax)
    close(2)

 else if (myid == 1) then

    do i=imin,imax
       
       n = 2**i 

       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       
       !!!  Processor 1 has to do the opposite, i.e.        !!!
       !!!  receive the message 50 times and send it back   !!!
       !!!  to Processor 0.                                 !!!
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


       !!!   ADD NECCESSARY RECEIVES AND SENDS HERE  !!!


    end do

 else
    print*,' Please start on 2 processors only! '
 end if    

 deallocate(buffer)

 call MPI_FINALIZE(rc)

end program test_comm
