subroutine bs1(A,b,n,x)

! Back substitution to solve A*x = b
! assuming that A contains its LU factors

implicit none
 
integer, intent(in) :: n
real (kind=8), intent(in), dimension(n,n) :: A
real (kind=8), intent(in), dimension(n) :: b
real (kind=8), intent(out), dimension(n) :: x

real (kind=8), dimension(n) :: y

integer i,j

! first solve L*y = b

y(1)  = b(1)

do i = 2,n
   y(i) = b(i)
   do j= 1,i-1
      y(i) = y(i) - A(i,j)*y(j)
   end do
end do 

! Now solve U*x = y 

x(n) = y(n)/A(n,n);

do i = n-1,1,-1

  x(i) = y(i)
  do j = i+1,n
     x(i) = x(i) - A(i,j)*x(j)
  end do
  x(i) = x(i)/A(i,i)

end do 


end subroutine bs1








