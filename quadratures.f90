!-----------------------------------------------------------------------------!
! MAE 5093 DOCUMENTATION | Engineering Numerical Analysis
!-----------------------------------------------------------------------------!
! >>> Numerical Integration with Quadratures
!-----------------------------------------------------------------------------!
! References: 
! * Fundamentals of Engineering Numerical Analysis by P. Moin (2012) -- Ch.1
! * Numerical Recipes: The Art of Scientific Computing, 2nd Edition (1992) 
!-----------------------------------------------------------------------------!
! Written by Omer San
!            CFDLab, Oklahoma State University, cfdlab.osu@gmail.com
!            www.cfdlab.org
! 
! Last updated: Sep. 8, 2015
!-----------------------------------------------------------------------------!


program integration
implicit none
real*8 ::f,a,b,eI,S1,k(5),w(5)
integer::j

!exact solution:
eI = 0.5d0*dlog(8.0d0)**2

a =  1.0d0 !lower bound
b =  8.0d0 !upper bound

!Gauss-Legendre (n=5)
k(1) = 0.0d0
k(2) = 1.0d0/3.0d0*dsqrt(5.0d0-2.0d0*dsqrt(10.0d0/7.0d0))
k(3) =-1.0d0/3.0d0*dsqrt(5.0d0-2.0d0*dsqrt(10.0d0/7.0d0))
k(4) = 1.0d0/3.0d0*dsqrt(5.0d0+2.0d0*dsqrt(10.0d0/7.0d0))
k(5) =-1.0d0/3.0d0*dsqrt(5.0d0+2.0d0*dsqrt(10.0d0/7.0d0))

w(1) = 128.0d0/225.0d0
w(2) = (322.0d0 + 13.0d0*dsqrt(70.0d0))/900.0d0
w(3) = (322.0d0 + 13.0d0*dsqrt(70.0d0))/900.0d0
w(4) = (322.0d0 - 13.0d0*dsqrt(70.0d0))/900.0d0
w(5) = (322.0d0 - 13.0d0*dsqrt(70.0d0))/900.0d0



S1 = 0.0d0

do j=1,5
S1 = S1 + (b-a)/2.0d0*f((b+a)/2.0d0 + (b-a)/2.0d0*k(j))*w(j)
end do

write(*,19)"exact: ", eI, 100.0d0*dabs((eI-eI)/eI)
write(*,19)"Gauss-Legendre: ", S1, 100.0d0*dabs((S1-eI)/eI)


19 format(A20,2F20.10)



end

!-----------------------------------------------------------------------------!
!Given function to integrate
!-----------------------------------------------------------------------------!
real*8 function f(x)
implicit none
real*8 :: x
f = dlog(x)/x
end function f
