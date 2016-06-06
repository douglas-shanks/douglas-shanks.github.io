!
!  Description: Testing the Matrix-Vector multiplication for compressed 
!               row storage format
!
! -----------------------------------------------------------------------

program test_mult

  use header

  type(Matrix)  :: A
  type(Vector)  :: u, b
  integer       :: m, n, flag, i

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     Beginning of program
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

  print*,'Number of discretisation points per coordinate direction:'
  read*, m
  print*, 'value of m =',  m

  n = (m-1)*(m-1)
 
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!     Allocate memory for A, u, and b and set dimensions
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

  allocate(A%aa(5*n))
  allocate(A%jj(5*n))
  allocate(A%ii(n+1))

  allocate(u%xx(n))
  allocate(b%xx(n))

  A%n = n
  u%n = n
  b%n = n

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!     Construct the matrix A
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

  call Laplace(A,m) 

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!     Set u to the vector of all ones and then multiply it by A
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

  u%xx = 1.0_8

  call Mat_Mult(A,u,b)

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!     Print out the solution to verify code
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

  print*,'Do you want to print the solution? (0 ... no, 1 ... yes)'
  read*, flag

  if (flag == 1) then
     do i=0,m-2
        print*, real(b%xx(i*(m-1)+1:(i+1)*(m-1)))
     end do
  end if

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!     Deallocate memory
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

  deallocate(A%aa)
  deallocate(A%jj)
  deallocate(A%ii)
  deallocate(u%xx)
  deallocate(b%xx)

end program test_mult




