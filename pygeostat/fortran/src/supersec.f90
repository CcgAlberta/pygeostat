
module supersecondary

    implicit none

    public :: supersec

    real*8, allocatable, public :: supvals(:)

    contains

    subroutine supersec(nxyz, nsec, secvals, rho, wts, tmin, tmax)
!--------------------------------------------------------------------------------------------------
!
!   This subroutine calculates a super-secondary value for each cell combining multiple secondary
!   variables.
!   
!   Code retrieved and modified from Ryan M. Barnett's supersec (March 2013) program.
!   
!   Parameters:
!     nxyz(int):       Number of cells
!     nsec(int):       Number of secondary variables to combine
!     secvals(double): 2-D array containing the secondary values to combine
!     rho(double):     Correlation coefficient of the final super-secondary variable and the
!                      primary variable
!     wts(double):     1-D array containing the weights for each of the secondary variables
!     tmin(double):    Minimum allowable super-secondary value
!     tmax(double):    Maximum allowable super-secondary value  
!   
!   Author: Warren Black -  May 26, 2016
!
!--------------------------------------------------------------------------------------------------

        implicit none

        integer, intent(in)   :: nxyz, nsec
        real*8,  intent(in)   :: tmin,tmax,rho
        real*8,  intent(in)   :: wts(nsec)
        real*8,  intent(in)   :: secvals(nxyz, nsec)

        real*8 :: ming,maxg,minp,maxp

        integer :: error   ! Error flag
        integer :: i

        if( allocated(supvals) ) deallocate(supvals)
        allocate( supvals(nxyz), stat=error )

        !Establish the maximum/minimum allowed Gaussian value
        minp=1/real(nxyz)
        call gauinv(minp,ming,error)
        maxp=1-minp
        call gauinv(maxp,maxg,error)

        ! Calculate super-secondary values
        do i=1,nxyz
            if( all(secvals(i,1:nsec) >= tmin ) .and. all( secvals(i,1:nsec) <= tmax ) )then
                supvals(i) = sum( wts * secvals(i,1:nsec) ) / rho
                if( supvals(i) > maxg ) supvals(i) = maxg
                if( supvals(i) < ming ) supvals(i) = ming
            endif
        enddo

    end subroutine supersec 

    subroutine gauinv(p,xp,ierr)
!--------------------------------------------------------------------------------------------------

! Computes the inverse of the standard normal cumulative distribution function with a numerical
! approximation from : Statistical Computing, by W.J. Kennedy, Jr. and James E. Gentle, 1980,
! p. 95.
!
! INPUT/OUTPUT:
!
!   p    = double precision cumulative probability value: dble(psingle)
!   xp   = G^-1 (p) in single precision
!   ierr = 1 - then error situation (p out of range), 0 - OK
!
!--------------------------------------------------------------------------------------------------
    real*8 p0,p1,p2,p3,p4,q0,q1,q2,q3,q4,y,pp,lim,p,xp
    save   p0,p1,p2,p3,p4,q0,q1,q2,q3,q4,lim
    integer ierr
    ! Coefficients of approximation:
    data lim/1.0e-10/
    data p0/-0.322232431088/,p1/-1.0/,p2/-0.342242088547/
    data p3/-0.0204231210245/,p4/-0.0000453642210148/
    data q0/0.0993484626060/,q1/0.588581570495/,q2/0.531103462366/
    data q3/0.103537752850/,q4/0.0038560700634/
    ! Check for an error situation:
    ierr = 1
    if(p.lt.lim) then
          xp = -1.0e10
          return
    end if
    if(p.gt.(1.0-lim)) then
          xp =  1.0e10
          return
    end if
    ierr = 0      
    ! Get k for an error situation:
    pp   = p
    if(p.gt.0.5) pp = 1 - pp
    xp   = 0.0
    if(p.eq.0.5) return
    !
    ! Approximate the function:
    !
    y  = dsqrt(dlog(1.0/(pp*pp)))
    xp = real( y + ((((y*p4+p3)*y+p2)*y+p1)*y+p0) / ((((y*q4+q3)*y+q2)*y+q1)*y+q0) )
    if(real(p).eq.real(pp)) xp = -xp
    ! Return with G^-1(p):
    return

    end subroutine gauinv

end module supersecondary