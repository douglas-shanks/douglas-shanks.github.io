!==================================================================
!
!   Function evaluating the material parameter at the point (x,y) in 
! 	the PDE
! 		-\grad a(x,y) \grad U + sigma U = g,  in D
!									  U = 0,  on D
!   (At the moment alpha(x,y) = 1.0 everywhere, i.e. the Laplacian)
!
!==================================================================

function alpha(x,y) result (val)

use header

real(kind=8), intent(in) :: x,y
real(kind=8) :: val

val = 1.0d0

end function alpha
