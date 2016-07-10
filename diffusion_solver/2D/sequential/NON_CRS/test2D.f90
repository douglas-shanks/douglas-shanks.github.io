program test2D

!---------------------------------------------------------------------
!
!	A test program to check that 2D FTCS, BTCS and CN converges
!
!---------------------------------------------------------------------

implicit none

 integer:: N,m, Nt
 integer i,j
 real(kind=8), dimension(:,:), allocatable ::   Aftcs,Acn,Acnrhs,Abtcs
 real(kind=8), dimension(:), allocatable ::  x,xcn,y,Uexact,Uftcs1,Uftcs,Ucn,Ubtcs
 real(kind=8):: alpha,tmax, hx, ht, xx, yy, R1, alpha1, beta1, errorftcs, errorcn,errorbtcs

! Define Pi for use

 real (kind=8) :: PI = 4.0_8*ATAN(1.0_8)
 real :: t1,t2

! extra storage needed by dgesv, dpbsv:

 integer, dimension(:), allocatable :: ipiv1, ipiv2
 integer :: info1, info2

! Constants needed for dgmev

 alpha1=1.d0; beta1=0.d0

 WRITE(*,*)
 WRITE(*,*) '------------------------------------------------------------------'
 WRITE(*,*)
 WRITE(*,*) ' dU/dt = \alpha D^2 U in \Omega, U = 0 on d\Omega '
 WRITE(*,*)
 WRITE(*,*) ' This is a test program where we solve the 2D heat equation '
 WRITE(*,*) ' We compare the exact solution to that with the FTCS method '
 WRITE(*,*) ' BTCS method and Crank Nicholson method.'
 WRITE(*,*)
 WRITE(*,*) '------------------------------------------------------------------'
 WRITE(*,*)
 WRITE(*,*) ' Inputs:'
 WRITE(*,*)
 WRITE(*,*) ' What is alpha? '
 READ(*,*) 	alpha
 WRITE(*,*) ' What is tmax? '
 READ(*,*) 	tmax
 WRITE(*,*)
 WRITE(*,*) ' How many grid points spatially? '
 READ(*,*) 	m
 WRITE(*,*) ' How many time steps? '
 READ(*,*) 	Nt
 WRITE(*,*) 

! Total grid points, not including Dirichlet boundary

 N = (m-1)**2

 if (m<0 .or. Nt<0) then
 	print *, "N and Nt must be positive "
 	stop
 endif

! Allocate storage

 allocate(y(N))
 allocate(x(N))
 allocate(xcn(N))
 allocate(Uftcs1(N))
 allocate(Uexact(N))
 allocate(Uftcs(N))
 allocate(Ubtcs(N))
 allocate(Ucn(N))
 allocate(Aftcs(N,N))
 allocate(Acn(N,N))
 allocate(Abtcs(N,N))
 allocate(Acnrhs(N,N))
 allocate(ipiv1(N))
 allocate(ipiv2(N)) 

! Mesh spacing for time and space

 hx = 1.0_8/real(m-1,8)
 ht = tmax/real(Nt-1,8)
 R1 = alpha*ht/(hx**2)

 print*, 'Value of CFL coefficient'
 print*, 				R1
 print*, ' '
 if (R1 >= 0.25_8) then
	print *, "Cannot guarantee stability for FTCS "
 endif

!-----------------------------------------------------------
!					EXACT SOLUTION
!-----------------------------------------------------------

 call exactsoln2D(Uexact,alpha,tmax,m)

 open(unit=2, file='uexact.txt', ACTION="write", STATUS="replace")
 write(2,*)m
 write(2,*) 
 write(2, *)( Uexact)

!-----------------------------------------------------------
!						FTCS
!-----------------------------------------------------------

! Initial solution, u_0

 call exactsoln2D(Uftcs1,alpha,0.0_8,m)

 !open(unit=2, file='uftcs1.txt', ACTION="write", STATUS="replace")
 !write(2,*) m
 !write(2,*) 
 !write(2, *)( Uftcs1)

 Uftcs = Uftcs1

 call FTCSmatrix(Aftcs,alpha,tmax,m,Nt)
 
 !open(unit=2, file='Aftcs.txt', ACTION="write", STATUS="replace")
 !write(2,*) N
 !write(2,*) 
 !write(2,*) Aftcs

! Then it should just be a case of U^{j} = A*U{j-1}

 call cpu_time(t1)
 do j=2,Nt

 	x = Uftcs
	call dgemv('n',N,N,alpha1,Aftcs,N,x,1,beta1,y,1)
	Uftcs = y
	
 end do
 call cpu_time(t2)

 write(*,*)
 write(*,*) '------------------------------------------------------------------'
 print*,    ' FTCS Solution '
 write(*,*) '------------------------------------------------------------------'
 write(*,*) 

 open(unit=2, file='uftcs.txt', ACTION="write", STATUS="replace")
 write(2,*) m
 write(2,*) 
 write(2, *)( Uftcs)

! Compute the error

 call norm(errorftcs,Uexact,Uftcs,m)

 write(*,*) ' Error '
 write(*,*)
 print*, errorftcs	 
 write(*,*)
 print*, ' time for solve'   
 print'(f12.6)', t2-t1

!-----------------------------------------------------------
!						BTCS
!-----------------------------------------------------------

! Initial solution, u_0

 call exactsoln2D(Ubtcs,alpha,0.0_8,m)

 call cpu_time(t1)
 do j=2,Nt

	! Make BTCS matrix and RHS (This has to be done at each
	! iterate as apparetnly dgesv changes the matrix!

	call BTCSmatrix(Abtcs,alpha,tmax,m,Nt)

	! now solve, and produce output information 

	call dgesv(N,1,Abtcs,N,ipiv2,Ubtcs,N,info2)

 end do
 call cpu_time(t2)

 write(*,*)
 write(*,*) '------------------------------------------------------------------'
 print*,    ' BTCS Solution '
 write(*,*) '------------------------------------------------------------------'
 write(*,*) 

 open(unit=2, file='ubtcs.txt', ACTION="write", STATUS="replace")
 write(2,*) m
 write(2,*) 
 write(2, *)( Ubtcs)

! Compute the error

 call norm(errorbtcs,Uexact,Ubtcs,m)

 write(*,*) ' Error '
 write(*,*)
 print*, errorbtcs	 
 write(*,*)
 print*, ' time for solve'   
 print'(f12.6)', t2-t1

!-----------------------------------------------------------
!				Crank Nicholson
!-----------------------------------------------------------

! Initial solution, u_0
 
 call exactsoln2D(Ucn,alpha,0.0_8,m)

 call cpu_time(t1)
 do j=2,Nt

	! Make CN matrix and RHS (This has to be done at each
	! iterate as apparetnly dgesv changes the matrix!

	call CNmatrix(Acn,alpha,tmax,m,Nt)
	! make the rhs

	call CNrhs(Acnrhs,alpha,tmax,m,Nt)
	call dgemv('n',N,N,alpha1,Acnrhs,N,Ucn,1,beta1,xcn,1)	

	! now solve, and produce output information 

	call dgesv(N,1,Acn,N,ipiv1,xcn,N,info1)
	Ucn = xcn

 end do
 call cpu_time(t2)

 write(*,*)
 write(*,*) '------------------------------------------------------------------'
 print*,    ' Crank Nicholson Solution '
 write(*,*) '------------------------------------------------------------------'
 write(*,*)

 open(unit=2, file='ucn.txt', ACTION="write", STATUS="replace")
 write(2,*) m
 write(2,*) 
 write(2, *)( Ucn)

! Compute the error

 call norm(errorcn,Uexact,Ucn,m)

 write(*,*) ' Error '
 write(*,*)
 print*, errorcn	 
 write(*,*)
 print*, ' time for solve'   
 print'(f12.6)', t2-t1

!---------------------------------------------------------------
! deallocate storage
!---------------------------------------------------------------

 deallocate(x)
 deallocate(xcn)
 deallocate(y)
 deallocate(Uexact)
 deallocate(Uftcs1)
 deallocate(Uftcs,Ubtcs)
 deallocate(Ucn)
 deallocate(Aftcs,Abtcs)
 deallocate(Acn)
 deallocate(Acnrhs)
 deallocate(ipiv1)
 deallocate(ipiv2) 


 WRITE(*,*) '------------------------------------------------------------------' 
 WRITE(*,*) ' Program finished '
 WRITE(*,*) '------------------------------------------------------------------'

end program test2D
