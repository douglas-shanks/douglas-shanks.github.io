!
!  Description: Testing the Matrix-Vector multiplication for compressed 
!               row storage format (parallel version)
!
! -----------------------------------------------------------------------

program test_mult

  use header
  include "mpif.h"

  type(Matrix)  :: A
  type(Vector)  :: u, b
  real(kind=8)  :: t1, t2
  integer       :: m, n, flag, i
  integer       :: myid, nprocs, nrows, ibeg, iend
  integer       :: ierr

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!      Beginning of program - Initialisation of MPI context
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, myid, ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, nprocs, ierr)

  call greetings(myid,nprocs)

  if ( myid == 0 ) then
     print*,'Number of discretisation points per coordinate direction:'
     open(unit=2,file="input.dat")
     read(2,*) m
     print*, 'Value of m =',  m

     if (mod((m-1)*(m-1),nprocs) /= 0) then
        print*,'(m-1)*(m-1) has to be a multiple of nprocs =',nprocs
        m = 0
     end if
  end if

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!      Broadcast m to the other processes
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

  call MPI_Bcast(m,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  if (m==0) stop
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!      Calculate the start and end indices of the rows to be held locally
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

  n     = (m-1)*(m-1)

  nrows = n / nprocs
  ibeg  = myid * nrows + 1
  iend  = (myid+1) * nrows

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!     Allocate memory for A, u, and b and set dimensions
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

  allocate(A%aa(5*nrows))
  allocate(A%jj(5*nrows))
  allocate(A%ii(n+1))

  A%n    = n
  A%ibeg = ibeg
  A%iend = iend

  allocate(u%xx(n))
  allocate(b%xx(n))

  b%n    = n
  b%ibeg = ibeg
  b%iend = iend

  u%n    = n
  u%ibeg = ibeg
  u%iend = iend

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!     Construct the matrix A
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

  call Laplace(A,m,ibeg,iend) 

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!     Set u to the vector of all ones and then multiply it by A
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

  u%xx(ibeg:iend) = 1.0_8

  call MPI_Barrier(MPI_COMM_WORLD,ierr)
  t1 = MPI_WTime()

  call Mat_Mult(A,u,b)

  call MPI_Barrier(MPI_COMM_WORLD,ierr)
  t2 = MPI_WTime()

  if (myid == 0) write(*,*) 'Elapsed time on ',nprocs,' processes:',t2-t1

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!     Print out the solution to verify code
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

  if (myid == 0) then 

     print*,'Do you want to print the solution? (0 ... no, 1 ... yes)'
     read(2,*) flag
     close(2)

  end if

  call MPI_Bcast(flag,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

  if (flag == 1) then
        
     n_loc = b%iend - b%ibeg + 1
     call MPI_Gather(b%xx(u%ibeg),n_loc,MPI_DOUBLE_PRECISION,b%xx,  &
          &          n_loc,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

     if (myid == 0) then 
        do i=0,m-2
           print*, real(b%xx(i*(m-1)+1:(i+1)*(m-1)))
        end do
     end if

  end if

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!     Deallocate memory
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

  deallocate(A%aa)
  deallocate(A%jj)
  deallocate(A%ii)
  deallocate(u%xx)
  deallocate(b%xx)

  call MPI_Finalize(ierr)

end program test_mult




