subroutine laplace(A,m) 

! fortran95  function  to assemble the finite  matrix      
! arising from the finite difference  discretisation of the Laplacian
! partial differential operator. This is not the same as the Laplacian 
! of a graph 

integer, intent(in) :: m
real(kind=8),  intent(out),  dimension((m-1)**2,(m-1)**2) :: A    
real(kind=8), dimension(m-1,m-1) :: B,ID
integer :: i

! initialise with zeros

 A = 0.0_8
 B = 0.0_8
 ID = 0.0_8
           
 do i = 1,m-1
   B(i,i) = 4.0_8
   ID(i,i) = -1.0_8
   if (i>1) then  
     B(i,i-1) = -1.0_8
   end if
   if (i<m-1) then 
     B(i,i+1) = -1.0_8
   end if      
 end do  


 do i=1,m-1
   A((m-1)*(i-1)+1:(m-1)*i,(m-1)*(i-1)+1:(m-1)*i) = B
   if (i>1) then 
     A((m-1)*(i-1)+1:(m-1)*i,(m-1)*(i-2)+1:(m-1)*(i-1))= ID
   end if

   if (i<m-1) then 
     A((m-1)*(i-1)+1:(m-1)*i,(m-1)*i+1:(m-1)*(i+1)) = ID
   end if
 end do 

 A = (m**2)*A

end subroutine laplace












