!
!  Description: Solves Poisson's problem by Jacobi's method 
!               (sequential version)
!
! -----------------------------------------------------------------------

program iterative

  use header

  real (kind=8) :: eps
  integer :: kmax
  parameter (eps = 1.0d-8, kmax = 9999)

  type(Matrix)  :: A
  type(Vector)  :: u, u_ex, b
  real (kind=8) :: norm
  integer       :: m, n, its

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Beginning of program
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

  print*,'type m'
  read*, m
  print*, 'value of m =',  m
 
  n = (m-1)*(m-1)
 
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!     Allocate memory for A, u, u_ex, and b and set dimensions
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

  allocate(A%aa(5*n))
  allocate(A%jj(5*n))
  allocate(A%ii(n+1))

  allocate(b%xx(n))
  allocate(u%xx(n))
  allocate(u_ex%xx(n))

  A%n    = n
  b%n    = n
  u%n    = n
  u_ex%n = n

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!     Construct the matrix A
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

  call Laplace(A,m) 

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!     Set exact solution to a random vector and then calculate rhs
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

  call random_number(u_ex%xx) 

  call Mat_Mult(A,u_ex,b)

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!     Apply Jacobi's method to solve the system
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  
  call Jacobi(A,u,b,eps,kmax,its)

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!     Check the error
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

  call Error(u,u_ex,norm)
  
  if (its > kmax) then
     write(*,100) kmax,norm
  else
     write(*,110) its,norm
  endif

100 format('Maximum number of iterations (',i5,          &
    &      ') reached, norm of the error is ',e10.4)
110 format('After ',i5,' iterations the norm of the error is ',e10.4)

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!     Deallocate memory
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

  deallocate(A%aa)
  deallocate(A%jj)
  deallocate(A%ii)
  deallocate(u_ex%xx)
  deallocate(u%xx)
  deallocate(b%xx)

end program iterative




