subroutine rhs(b,m)

! finds the right-hand side vector in the approximation of Poisson's 
! equation on the unit square

integer, intent(in) :: m
real(kind=8), intent(out), dimension((m-1)**2) :: b
integer :: i,j
real (kind=8) :: x1,x2,h

h = 1.0_8/real(m,8)

do  i = 1,m-1

  do j = 1,m-1

    x1 = i*h
    x2 = j*h

    b(i+(j-1)*(m-1)) = (3.0_8*x1+x1**2)*exp(x1)*x2*(1-x2) &
                    + 2*(x1*(1-x1))*exp(x1)  
  end do
end do
 

end subroutine rhs
