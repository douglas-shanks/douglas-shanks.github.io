! main program to test lu

program testlu

real (kind=8), allocatable, dimension(:,:) :: A
real (kind=8), allocatable, dimension(:) :: x,b
real (kind=8) :: t1, t2
integer, allocatable, dimension(:) :: P
integer :: n, m

integer, dimension(:), allocatable :: ipiv
integer :: info

print*,'type m'
read*, m
print*, 'value of m =',  m

n = (m-1)**2

allocate(A(n,n))
allocate(P(n))
allocate(x(n))
allocate(b(n))
allocate(ipiv(n))

! use lu1.f90 and bs1.f90

 call laplace(A,m)     ! construct the matrix A and rhs b
 call rhs(b,m)

 call cpu_time(t1)
 call dgesv(n,1,A,n,ipiv,b,n,info)             
 call cpu_time(t2)

 print*,'  '
 print*,'results from general solver dgesv:'
 print*,  '    time for solve'   
 print'(f12.6)', t2-t1
 print*,'  ' 


 call cpu_time(t1)
call lu1(A,n)             ! get   LU factors of  A                       
call bs1(A,b,n,x)         ! solve LUx = b
call cpu_time(t2)

 print*,'  '
 print*,'results from bs1:'
 print*,  '    time for solve'   
 print'(f12.6)', t2-t1
 print*,'  ' 

!print*,' '
!print*,'computed solution '
!print*,x

end program testlu

