subroutine rsymm(A,n)

! assembles a random well-conditioned n times n symmetric matrix A
! It is just a small symmetric random perturbation of the identity
! and so should be well-conditioned

integer, intent(in) :: n
real (kind=8), dimension(n,n), intent(out) :: A
integer :: i

call random_number(A)

A = 0.01*(A*transpose(A))    ! small symmetric matrix

do i = 1,n
  A(i,i) = A(i,i) + 1.0        
end do

end subroutine rsymm

