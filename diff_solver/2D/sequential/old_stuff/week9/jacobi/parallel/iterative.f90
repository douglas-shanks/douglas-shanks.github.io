!
!  Description: Solves Poisson's problem by Jacobi's method 
!               (parallel version)
!
! -----------------------------------------------------------------------

program iterative

  use header
  include "mpif.h"


  real (kind=8) :: eps
  integer :: kmax
  parameter (eps = 1.0d-8, kmax = 9999)

  type(Matrix)  :: A
  type(Vector)  :: u, u_ex, b
  integer       :: m
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
     close(2)
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
!     Allocate memory for A, u, u_ex and b and set dimensions
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

  allocate(A%aa(5*nrows))
  allocate(A%jj(5*nrows))
  allocate(A%ii(n+1))

  A%n    = n
  A%ibeg = ibeg
  A%iend = iend

  allocate(u%xx(n))
  allocate(b%xx(n))
  allocate(u_ex%xx(n))

  b%n    = n
  b%ibeg = ibeg
  b%iend = iend

  u%n    = n
  u%ibeg = ibeg
  u%iend = iend

  u_ex%n    = n
  u_ex%ibeg = ibeg
  u_ex%iend = iend

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!     Construct the matrix A
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

  call Laplace(A,m,ibeg,iend) 

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!     Set exact solution u_ex to a random vector and calculate rhs b = A*u
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

  call random_number(u_ex%xx(ibeg:iend))

  call Mat_Mult(A,u_ex,b)

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!     THIS IS WHERE THE CALL TO THE PARALLEL JACOBI METHOD SHOULD GO !!!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

! |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv





! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
! |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 


! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!     Deallocate memory
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

  deallocate(A%aa)
  deallocate(A%jj)
  deallocate(A%ii)
  deallocate(u_ex%xx)
  deallocate(u%xx)
  deallocate(b%xx)

  call MPI_Finalize(ierr)

end program iterative




