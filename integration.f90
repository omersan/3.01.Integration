!-----------------------------------------------------------------------------!
! MAE 5093 DOCUMENTATION | Engineering Numerical Analysis
!-----------------------------------------------------------------------------!
! >>> Numerical Integration 
!-----------------------------------------------------------------------------!
! References: 
! * Fundamentals of Engineering Numerical Analysis by P. Moin (2012) -- Ch.1
! * Numerical Recipes: The Art of Scientific Computing, 2nd Edition (1992) 
!-----------------------------------------------------------------------------!
! Written by Omer San
!            CFDLab, Oklahoma State University, cfdlab.osu@gmail.com
!            www.cfdlab.org
! 
! Last updated: Sep. 3, 2015
!-----------------------------------------------------------------------------!


program integration

implicit none
integer::nx
real*8 ::dx,x0,xN,pi,eI,gd0,gdN,s1,s2,s3
integer::i
real*8,allocatable :: x(:), g(:)

nx = 32

pi = 4.0d0*datan(1.0d0)

x0 = 1.0d0
xN = pi

dx = (xN-x0)/dfloat(nx)

!grid
allocate(x(0:nx))
do i=0,nx
x(i) = x0 + dfloat(i)*dx
end do

!function
allocate(g(0:nx))
do i=0,nx
g(i) = dsin(x(i))/(2.0d0*x(i)**3)
end do

!exact integration for the function above:
eI = 0.1985572988d0

!derivative data at the boundaries (only needed for end correction in trap method)
gd0 = (x(0)*dcos(x(0)) -3.0d0*dsin(x(0)))/(2.0d0*x(0)**4)
gdN = (x(nx)*dcos(x(nx)) -3.0d0*dsin(x(nx)))/(2.0d0*x(nx)**4)

!compute integrals with different methods

call trap1D(nx,dx,g,s1)

call simp1D(nx,dx,g,s2)

call trapend1D(nx,dx,g,gd0,gdN,s3)


!write solutions
write(*,100)"exact: ", eI, 100.0d0*dabs(eI-eI)/eI
write(*,100)"trapezoidal: ", s1, 100.0d0*dabs(s1-eI)/eI
write(*,100)"Simpson: ", s2, 100.0d0*dabs(s2-eI)/eI
write(*,100)"trapezoidal+end: ", s3, 100.0d0*dabs(s3-eI)/eI



100 format(A20,2F20.6)

end


!---------------------------------------------------------------------------!
!Numerical integral routines
!---------------------------------------------------------------------------!

!----------------------------------------------------------!
!trapezoidal rule for numerical integration of g(i)
!for equally distributed mesh with interval dx 
!second-order accurate
!----------------------------------------------------------!
subroutine trap1D(nx,dx,g,s)
implicit none
integer::nx,i
real*8 ::dx,g(0:nx),s,ds,th

	th = 0.5d0*dx
    
	s  = 0.0d0
	do i=0,nx-1
	ds = th*(g(i)+g(i+1))
	s  = s + ds
	end do


return
end

!----------------------------------------------------------!
!trapezoidal rule for numerical integration of g(i)
!plus end-correction
!for equally distributed mesh with interval dx 
!fourth-order accurate
!----------------------------------------------------------!
subroutine trapend1D(nx,dx,g,gd0,gdN,s)
implicit none
integer::nx,i
real*8 ::dx,g(0:nx),s,ds,th,gd0,gdN

	th = 0.5d0*dx
    
	s  = 0.0d0
	do i=0,nx-1
	ds = th*(g(i)+g(i+1))
	s  = s + ds
	end do

    !end correction
    s = s - dx*dx/12.0d0*(gdN-gd0)


return
end

!----------------------------------------------------------!
!Simpson's 1/3 rule for numerical integration of g(i)
!for equally distributed mesh with interval dx 
!nx should be even number
!fourth-order accurate
!----------------------------------------------------------!
subroutine simp1D(nx,dx,g,s)
implicit none
integer::nx,i,nh
real*8 ::dx,g(0:nx),s,ds,th

	nh = int(nx/2)
	th = 1.0d0/3.0d0*dx
    
	s  = 0.0d0
	do i=0,nh-1
	ds = th*(g(2*i)+4.0d0*g(2*i+1)+g(2*i+2))
	s  = s + ds
	end do

return
end

!----------------------------------------------------------!
!high-order (fourth) integration
!doesn't assume even points
!5th order scheme for numerical integration of g(i)
!4th order on the boundaries
!for equally distributed mesh with interval dx 
!----------------------------------------------------------!
subroutine ho1D(nx,dx,g,s)
implicit none
integer::nx,i
real*8 ::dx,g(0:nx),s,ds,th


	th = dx/24.0d0

	s = 0.0d0
	do i=1,nx-2
	ds = th*(-g(i-1)+13.0d0*g(i)+13.0d0*g(i+1)-g(i+2))
	s = s + ds
	end do
    !boundaries
	i=0
	s = s + 2.0d0*th*(5.0d0*g(i)+8.0d0*g(i+1)-g(i+2)) 
	i=nx-1
	s = s + 2.0d0*th*(5.0d0*g(i+1)+8.0d0*g(i)-g(i-1)) 

	
return
end
