program testsolvers

! approximates the solution to Poisson's equation 
!
!    - Laplacian (u) = f
!
! on the unit square with zero dirichlet conditions, by finite 
! difference methods.
! This leads to a linear system:
!
!            A x = b
! 
! The matrix A is produces by the function laplace 
! The right hand side vector b is produced by the function rhs
!  
! Uses LAPACK to solve this linear system,
! use the standard drivers sgesv and spbsv
 
implicit none

integer :: m,n
real (kind=8),  dimension(:,:), allocatable :: A
real (kind=8),  dimension(:), allocatable :: b,x
integer :: i,j
real :: t1,t2

! extra storage needed by dgesv, dpbsv:
real (kind=8), dimension(:,:), allocatable::AB,AFB
integer, dimension(:), allocatable :: ipiv
integer :: info

print*,'type m'
read*, m
print*, 'value of m =',  m

n = (m-1)**2

! allocate storage

 allocate(A(n,n))
 allocate(x(n))
 allocate(b(n))
 allocate(AB(m,n))
 allocate(ipiv(n))

 call laplace(A,m)     ! construct the matrix A and rhs b
open(unit=2, file='A.txt', ACTION="write", STATUS="replace")
write(2,*) n
write(2,*) 
write(2,*) A

 call rhs(b,m)

 x = b              ! save the rhs b into x

! now solve, and produce output information
  
 call cpu_time(t1)
 call dgesv(n,1,A,n,ipiv,x,n,info)             
 call cpu_time(t2)

 print*,'  '
 print*,'results from general solver dgesv:'
 print*,  '    time for solve'   
 print'(f12.6)', t2-t1
 print*,'  '  
 print*, 'error in fd solution  ', abs(x((n+1)/2) - &
         (1.0_8/16.0_8)*exp(0.5_8))
 print*,'  '

! redo this with spd  solver  

 call laplace(A,m)     ! construct the matrix A and rhs b

 call rhs(b,m)

 x = b              ! save the rhs b into x

 call cpu_time(t1)
 call dposv('Lower',n,1,A,n,x,n,info)             
 call cpu_time(t2)
 
 print*,'  '
 print*,'results from dposv:'
 print*,  '    time for solve'   
 print'(f12.6)', t2-t1
 print*,'  '  
 print*, 'error in fd solution  ', abs(x((n+1)/2) - &
         (1.0_8/16.0_8)*exp(0.5_8))
 print*,'  '


! redo this with the banded solver  

 call laplace(A,m)     ! construct the matrix A and rhs b

 call rhs(b,m)

 x = b              ! save the rhs b into x

! assemble the array AB needed in the banded solver:
! see the argument list for dpbsv

 do j = 1,n
    do i = j, min(n,j+m-1)
       AB(1+i-j,j) = A(i,j)
    end do
 end do

 call cpu_time(t1)
 call dpbsv('Lower',n,m-1,1,AB,m,x,n,info)             
 call cpu_time(t2)

 print*,'  '
 print*,'results from dpbsv:'
 print*,  '    time for solve'   
 print'(f12.6)', t2-t1
 print*,'  '  
 print*, 'error in fd solution  ', abs(x((n+1)/2) - &
         (1.0_8/16.0_8)*exp(0.5_8))
 print*,'  '

! deallocate storage

 deallocate(A)
 deallocate(x)
 deallocate(b)
 deallocate(AB)
 deallocate(ipiv)


end program testsolvers












