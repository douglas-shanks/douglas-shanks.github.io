program lsquares 

! solves the least squares data fitting problem 

implicit none
integer :: n,d
real (kind=8),  dimension(:,:), allocatable :: A,C
real (kind=8),  dimension(:), allocatable :: y,t,x,b
integer :: i,j

! read in the dimensions

open(unit=2,file='lsdata.txt')

read(2,*) n
print*, 'number of data samples =',  n

allocate(t(n))
allocate(y(n))

read(2,*) t(:)          ! read data 
read(2,*) y(:)
close(2)

print*,' '
print*,'data points:'
do i = 1,n
print*, t(i), y(i)    ! print data to screen (optional)
end do 
print*,' '

print*, 'Type the degree of polynomial you want to fit'
read*, d
print*,' '
print*, 'degree of polynomial to fit the data = ', d

allocate(C(n,d+1))
allocate(x(d+1))
allocate(A(d+1,d+1))
allocate(b(d+1))

! assemble the matrix C:

C(:,1) = 1.0_8

do i = 1,n
   do j = 2,d+1
      C(i,j) = t(i)**(j-1)
   end do
end do 

! solve the least squares problem

call nrmeq(C,y,n,d+1,A,b)  ! Calculates A = C^TC and b = C^Ty

call lu1(A,d+1)     ! computes the LU decomposition of A = C^TC

call bs1(A,b,d+1,x) ! finds the solution to Ax = b by back subsitution

print*, ' ' 
print*,'The coefficients of the lest squares polynomial are:'
print*, ' ' 
print*,x(1:d+1)
print*, ' ' 

! output the  results in a format which matlab can read 
 
open(3,file='lsoutput.txt')
write(3,*) d
write(3,*) ' '
write(3,*) (x(i),i=1,d+1)
close(3)

! good practice: deallocate storage

deallocate(t)
deallocate(C)
deallocate(x)
deallocate(y)
deallocate(A)
deallocate(b)


end program lsquares













