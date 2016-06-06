subroutine blublas2(A,n,bsize,k)

! Gaussian Elimination without pivoting
! using blas level 2 for a square bsize X bsize submatrix
! with top left entry A(k,k) of a square n x n matrix A.

integer, intent(in) :: n,bsize,k
real (kind=8), intent(inout), dimension(n,n) :: A
integer i

 do i = 1,bsize-1

   call dscal(bsize-i,1.0/A(k+i-1,k+i-1),A(k+i,k+i-1),1)
   call dger(bsize-i,bsize-i,-1.0d0,A(k+i,k+i-1),1, &
             A(k+i-1,k+i),n,A(k+i,k+i),n)
 end do 

end subroutine blublas2


