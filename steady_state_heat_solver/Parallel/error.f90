! -----------------------------------------------------------------------
subroutine Error(u,v,norm)
! -----------------------------------------------------------------------
!
!  Description: Subroutine to calculate || u - v ||_2 
!
! -----------------------------------------------------------------------

  use header
  include "mpif.h"

  type(Vector), intent(inout) :: u, v
  real(kind=8), intent(out) :: norm

  real(kind=8) :: mynrm
  integer :: i

  mynrm = 0.0_8
  
  do i=u%ibeg,u%iend
     mynrm = mynrm + (u%xx(i) - v%xx(i)) * (u%xx(i) - v%xx(i))
  end do

  call MPI_Allreduce(mynrm,norm,1,MPI_DOUBLE_PRECISION,MPI_SUM, &
                     MPI_COMM_WORLD,ierr)
  norm = sqrt(norm)

end subroutine Error
