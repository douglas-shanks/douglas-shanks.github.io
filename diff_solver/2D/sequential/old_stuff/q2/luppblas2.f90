subroutine luppblas2(A,n,p)

! Gaussian Elimination without pivoting
! using blas level 2 for a square n X n submatrix
  
integer, intent(in) :: n
real (kind=8), intent(inout), dimension(n,n) :: A
integer, intent(out), dimension(n) :: p

integer :: ipivot, ptemp, i, idamax

external idamax

! initialise the permutation vector p

 do i = 1,n
   p(i) = i
 end do

 do i = 1,n-1

! first find the row index of the new pivot and put it in ipivot

  ipivot = i - 1 + idamax(n-i+1,A(i,i),1)

! Swap row i and row ipivot in A and in p (only if i /= ipivot)

  if(ipivot /= i) then 

     call dswap(n,A(i,1),n,A(ipivot,1),n)

     ptemp = p(ipivot)
     p(ipivot) = p(i)
     p(i) = ptemp

  end if

! now carry on as before

  call dscal(n-i,1.0/A(i,i),A(i+1,i),1)

  call dger(n-i,n-i,-1.0d0,A(i+1,i),1,A(i,i+1),n,A(i+1,i+1),n)

end do 

end subroutine luppblas2



