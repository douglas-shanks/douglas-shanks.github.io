subroutine lublas3(A,n,bsize)

! Gaussian Elimination without pivoting
! using blocking and blas level 3 for a square n X n matrix
! Assumes the blocksize divides n 
  
integer, intent(in) :: n, bsize
real (kind=8), intent(inout), dimension(n,n) :: A

integer :: nb
integer ::  i, j

nb = n/bsize

print*,'number of blocks', nb

do j = 1,nb-1

  i = (j-1)*bsize + 1

  call blublas2(A,n,bsize,i)    ! LU decomposition of bsize x bsize 
                                ! block with leading entry A(i,i)  
                                ! Note that this block is overwritten 
                                ! with the LU factors 

  call dtrsm('Left','Lower','No Transpose','Unit Diagonal', &
             bsize,n-j*bsize,1.0_8, &
             A(i,i),n,A(i,i+bsize),n)    ! Step 2 of the algorithm: 
                                         ! a unit lower triangular solve 
                                         ! with many right-hand sides

  call dtrsm('Right','Upper','No Transpose','Non Unit Diagonal', &
             n-j*bsize,bsize,1.0_8, &
             A(i,i),n,A(i+bsize,i),n)    ! Step 3 of the algorithm: 
                                         ! an upper  triangular solve 
                                         ! with many right-hand sides

  call dgemm('No transpose','No transpose',  &
             n-j*bsize,n-j*bsize, bsize, & 
             -1.0_8,A(i+bsize,i),n,A(i,i+bsize),n, &
              1.0_8,A(i+bsize,i+bsize),n)  ! Step 4 of the algorithm:
                                         ! form the remainder matrix 
                                         ! which is still to be factorised. 

 end do 

 ! LU decomposition of final nb x nb  block
 
 call blublas2(A,n,n-(nb-1)*bsize,(nb-1)*bsize+1)
 
end subroutine lublas3










