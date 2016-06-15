! main program to test luppblas2

program testlupp

real (kind=8), allocatable, dimension(:,:) :: A
real (kind=8), allocatable, dimension(:) :: x_e,x,b

integer, allocatable, dimension(:) :: p
integer :: n

print*,'type n'
read*, n
print*, 'value of n =',  n

allocate(A(n,n))
allocate(p(n))
allocate(x(n))
allocate(x_e(n))
allocate(b(n))

! use luppblas2.f90 and bsblas2.f90

  call random_number(A) 

  call random_number(x_e) 

! assemble RHS: so that Ax = b has the exact solution x_e
  b  = matmul(A,x_e)

  call luppblas2(A,n,p)             ! get   LU factors of  A

  b = b(p(:))
  
  call bsblas2(A,b,n,x)           ! solve LUx = b

  print*, 'error in solve'   ! check that the solver has worked
  print*, sqrt(dot_product(x-x_e,x-x_e))

deallocate(A)
deallocate(p)
deallocate(x)
deallocate(x_e)
deallocate(b)

end program testlupp

