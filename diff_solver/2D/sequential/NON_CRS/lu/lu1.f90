subroutine lu1(A,n)
! Gaussian Elimination without pivoting
! overwrites A with the L,U factors
implicit none 

integer, intent(in) :: n
real (kind=8), intent(inout), dimension(n,n):: A

integer i,j,k

do i = 1,n-1

  do j = i+1,n
    A(j,i) = A(j,i)/A(i,i)
  end do

  do j = i+1,n
    do k = i+1,n
      A(j,k) = A(j,k) - A(j,i)*A(i,k)
    end do
  end do

end do

end subroutine lu1

