program testblas3

! compares times for blublas2 and lublas3 on random well-conditioned
! test matrices
 
integer :: n
real (kind=8),  dimension(:,:), allocatable :: A,A1
real (kind=8),  dimension(:), allocatable :: b,x,x_e

real :: t1, t2
integer :: isolve
integer :: bsize,ierror



print*,'type n'
read*, n
print*, 'value of n =',  n

! allocate storage

 allocate(A(n,n))
 allocate(A1(n,n))
 allocate(x(n))
 allocate(x_e(n))   
 allocate(b(n))

! create a random well-conditioned symmetric positive definite matrix

 call rsymm(A,n)

 A1 = A                     ! save as A1 (only needed for the error check)
 
! do i = 1,n
! print*,'matrix A'
! print*,A(i,:)
! end do 

!  random solution vector
  call random_number(x_e) 

! assemble RHS: so that Ax = b has the exact solution x_e
  b  = matmul(A,x_e)

! now solve

 isolve = 1 

 do 
 
    print*,'which solver?  2  for blublas2, 3 for lublas3 '
    read*,isolve

    call cpu_time(t1)    ! start the clock
      
    if(isolve == 2) then  
       print*,'using solver blublas2'
       call blublas2(A,n,n,1)
    else
       print*,'using solver lublas3'
       print*,'type block size'
       read*,bsize
       call lublas3(A,n,bsize)
    end if
      
    call cpu_time(t2)    ! stop the clock

    call bsblas2(A,b,n,x)     ! backsolve

    print*, 'error in solve'   ! check that the solver has worked
    print*, sqrt(dot_product(x-x_e,x-x_e))

    print*, 'Time taken for solver', t2 - t1, 'seconds   ' 
                                 ! print the time 
    
    print*,'Type 1 for another solve, 0 to stop'
    read*, isolve

    if (isolve == 0) exit
    
    A = A1             ! reload A

 end do 

 deallocate(A)
 deallocate(A1)
 deallocate(x)
 deallocate(x_e)   
 deallocate(b)


end program testblas3











