subroutine bsblas2(A,b,n,x)

! Back substitution to solve A*x = b
! assuming that A contains its LU factors 
! using level 2 blas

implicit none  
integer, intent(in) :: n
real (kind=8), intent(in), dimension(n,n) :: A
real (kind=8), intent(in), dimension(n) :: b
real (kind=8), intent(out), dimension(n) :: x

! avoid overwriting b:

x = b

! first solve L*y = b  Writes the answer into x

call dtrsv('L','N','U',n,A,n,x,1)

! Now solve U*x = y   again writes the answer into x

call dtrsv('U','N','N',n,A,n,x,1)

 

end subroutine bsblas2








