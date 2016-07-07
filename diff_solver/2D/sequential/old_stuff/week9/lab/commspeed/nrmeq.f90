subroutine nrmeq(C,y,n,m,A,b)
implicit none 

integer, intent(in) :: n,m
real (kind=8), intent(in),dimension(n,m):: C
real (kind=8), intent(in),dimension(n):: y
real (kind=8), intent(out),dimension(m,m):: A
real (kind=8), intent(out),dimension(m):: b

integer i,j,k

do i = 1,m

  do j = 1,m
    A(i,j) = 0.0_8
    do k = 1,n
      A(i,j) = A(i,j) + C(k,i)*C(k,j)
    end do
  end do

  b(i) = 0.0_8
  do k = 1,n
    b(i) = b(i) + C(k,i)*y(k)
  end do 

end do

end subroutine nrmeq

