subroutine norm(error,U1,U2,m)

! This actually computes the error not the norm

integer, intent(in) :: m
real(kind=8), intent(in), dimension((m-1)**2) :: U1
real(kind=8), intent(in), dimension((m-1)**2) :: U2
real(kind=8), intent(out) :: error

error = 0
do i=1,(m-1)
	error = error + SQRT((U1(i)- U2(i))**2) / SQRT(real((m-1)**2,8))
enddo

end subroutine norm
