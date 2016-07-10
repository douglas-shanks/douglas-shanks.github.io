subroutine CNmatrix(A,alpha,tmax,m,Nt)

! Subroutine to form matrix arriving from 
! CN for the heat equation
! Non CRS version

 integer, intent(in) :: m,Nt
 real(kind=8), intent(in) :: alpha, tmax
 real(kind=8),  intent(out),  dimension((m-1)**2,(m-1)**2) :: A
 real(kind=8), dimension(m-1,m-1) :: B,ID
 real(kind=8):: hx, ht, R1

! Mesh spacing for time and space

 hx = 1.0_8/real(m,8)
 ht = tmax/real(Nt-1,8)
 R1 = alpha*ht/(hx**2)

! initialise with zeros

 A = 0.0_8
 B = 0.0_8
 ID = 0.0_8
           
 do i = 1,m-1
   B(i,i) = 1.0_8 + 2.0_8*R1
   ID(i,i) = -R1/(2.0_8)
   if (i>1) then  
     B(i,i-1) = -R1/(2.0_8)
   end if
   if (i<m-1) then 
     B(i,i+1) = -R1/(2.0_8)
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
  
end subroutine CNmatrix

