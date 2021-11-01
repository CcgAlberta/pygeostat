
module backtr
    
    implicit none

    public :: backtransformdata

    contains


    subroutine backtransformdata(nd,nt,vrgs,vr,vrg,zmin,zmax,ltail,ltpar,utail,utpar,ltrim, &
                                 utrim,values)
!--------------------------------------------------------------------------------------------------
!
! This subroutine backtransforms a standard normal array by a specified back transform table. Can 
! also specify an option for the tails of the distribution.
!
! Parameters:
!   nd               number of data considered 
!   nt               number of values in the back transform tbale
!   vrgs             normal score value to be back transformed
!   vr(nt)           original data values that were transformed
!   vrg(nt)          the corresponding transformed values
!   zmin,zmax        limits possibly used for linear or power model
!   ltail            option to handle values less than vrg(1):
!   ltpar            parameter required for option ltail
!   utail            option to handle values greater than vrg(nt):
!   utpar            parameter required for option utail
!   ltrim            lower trimming limit
!   utrim            upper trimming limit
!   backtr           back-transformed values
!
! Warren Black : Minor adjustments to backtr.f90 by Clayton V. Deutch (GSLIB) for wrapping in
!                python.
!----------------------------------------------------------------------------------------------

        implicit none

        integer, intent(in)  :: nd, nt
        real*8,  intent(in)  :: vrgs(:), vr(:), vrg(:)
        real*8,  intent(in)  :: zmin, zmax, ltpar, utpar, ltrim, utrim
        integer, intent(in)  :: ltail, utail
        real*8,  intent(out) :: values(nd)
        integer              :: i, test

        do i = 1, nd
            if (vrgs(i) < ltrim .or. vrgs(i) > utrim)then
                values(i) = -999
            else
                values(i) = backtrval(vrgs(i), nt, vr, vrg, zmin, zmax, ltail, ltpar, utail, utpar)
            endif
        end do

    end subroutine backtransformdata 

    real*8 function backtrval (vrgs,nt,vr,vrg,zmin,zmax,ltail,ltpar,utail,utpar) result(value)
!----------------------------------------------------------------------------------------------
!
! This subroutine backtransforms a standard normal deviate from a specified back transform
! table and option for the tails of the distribution.  Call once with "first" set to true then
! set to false unless one of the options for the tail changes.
!
! Parameters:
!   nt               number of values in the back transform tbale
!   vrgs             normal score value to be back transformed
!   vr(nt)           original data values that were transformed
!   vrg(nt)          the corresponding transformed values
!   zmin,zmax        limits possibly used for linear or power model
!   ltail            option to handle values less than vrg(1):
!   ltpar            parameter required for option ltail
!   utail            option to handle values greater than vrg(nt):
!   utpar            parameter required for option utail
!   value            value being returned
!
! Warren Black : Minor adjustments to backtr.f90 by Clayton V. Deutch (GSLIB) for wrapping in
!                python.
!----------------------------------------------------------------------------------------------
        real*8,    parameter   :: EPSLON=1.0e-20
        integer, intent(in)  :: nt
        real*8,  intent(in)  :: vr(nt), vrg(nt)
        real*8,  intent(in)  :: vrgs, zmin, zmax, ltpar, utpar
        integer, intent(in)  :: ltail, utail
        integer              :: j
        real*8               :: cdfhi, cdflo, cdfbt, cpow, lambda
        !
        ! Value in the lower tail?    1=linear, 2=power, (3 and 4 are invalid):
        !
        if(vrgs.le.vrg(1)) then
            value = vr(1)
            cdflo  = gcum(vrg(1))
            cdfbt  = gcum(vrgs)
            if(ltail.eq.1) then
                 value = powint(DBLE(0.0),cdflo,zmin,vr(1),cdfbt,DBLE(1.0))
            else if(ltail.eq.2) then
                cpow   = 1.0 / ltpar
                value = powint(DBLE(0.0),cdflo,zmin,vr(1),cdfbt,cpow)
            endif
        !
        ! Value in the upper tail?     1=linear, 2=power, 4=hyperbolic:
        !
        else if(vrgs.ge.vrg(nt)) then
            value = vr(nt)
            cdfhi  = gcum(vrg(nt))
            cdfbt  = gcum(vrgs)
            if(utail.eq.1) then
                value = powint(cdfhi,DBLE(1.0),vr(nt),zmax,cdfbt,DBLE(1.0))
            else if(utail.eq.2) then
                cpow   = 1.0 / utpar
                value = powint(cdfhi,DBLE(1.0),vr(nt),zmax,cdfbt,cpow)
            else if(utail.eq.4) then
                lambda = (vr(nt)**utpar)*(1.0-gcum(vrg(nt)))
                value = (lambda/(1.0-gcum(vrgs)))**(1.0/utpar)
            endif
        else
        !
        ! Value within the transformation table:
        !
            call locate(vrg,nt,1,nt,vrgs,j)
            j = max(min((nt-1),j),1)
            value = powint(vrg(j),vrg(j+1),vr(j),vr(j+1),vrgs,DBLE(1.0))
        endif
        return

    end function backtrval

    real*8 function gcum(x) result(value)
!-----------------------------------------------------------------------
!
! Evaluate the standard normal cdf given a normal deviate x.  gcum is
! the area under a unit normal curve to the left of x.  The results are
! accurate only to about 5 decimal places.
!
!-----------------------------------------------------------------------
        real*8, intent(in)  :: x
        real*8              :: z, t, e2
        z = x
        if(z.lt.0.) z = -z
        t    = 1./(1.+ 0.2316419*z)
        value = t*(0.31938153   + t*(-0.356563782 + t*(1.781477937 + &
               t*(-1.821255978 + t*1.330274429))))
        e2   = 0.
        !
        !  6 standard deviations out gets treated as infinity:
        !
        if(z.le.6.) e2 = exp(-z*z/2.)*0.3989422803
        value = 1.0- e2 * value
        if(x.ge.0.) return
        value = 1.0 - value
        return
    end function gcum

    real*8 function powint(xlow,xhigh,ylow,yhigh,xval,pow) result(value)
!-----------------------------------------------------------------------
!
! Power interpolate the value of y between (xlow,ylow) and (xhigh,yhigh)
!                 for a value of x and a power pow.
!
!-----------------------------------------------------------------------
        real*8,   parameter    :: EPSLON = 1.0e-20
        real*8, intent (in)  :: xlow, xhigh, ylow, yhigh, xval, pow

        if((xhigh-xlow).lt.EPSLON) then
            value = (yhigh+ylow)/2.0
        else
            value = ylow + (yhigh-ylow)*(((xval-xlow)/(xhigh-xlow))**pow)
        end if

        return

        print *, value

    end function powint

    subroutine locate(xx,n,is,ie,x,j)
!-----------------------------------------------------------------------
!
! Given an array "xx" of length "n", and given a value "x", this routine
! returns a value "j" such that "x" is between xx(j) and xx(j+1).  xx
! must be monotonic, either increasing or decreasing.  j=is-1 or j=ie is
! returned to indicate that x is out of range.
!
! Bisection Concept From "Numerical Recipes", Press et. al. 1986  pp 90.
!-----------------------------------------------------------------------
        implicit none
        integer, intent(in) :: n,ie,is
        real*8, dimension(n), intent(in) :: xx
        real*8, intent(in) :: x
        integer, intent(out) :: j
        integer jl,ju,jm,is_modif
        !
        ! Initialize lower and upper methods:
        !
        if (is.le.0) then
          is_modif = 1
        else
          is_modif = is
        end if
        jl = is_modif-1
        ju = ie
        if(xx(n).le.x) then
              j = ie
              return
        end if
        !
        ! If we are not done then compute a midpoint:
        !
     10 if(ju-jl.gt.1) then
            jm = (ju+jl)/2
            !
            ! Replace the lower or upper limit with the midpoint:
            !
            if((xx(ie).gt.xx(is_modif)).eqv.(x.gt.xx(jm))) then
                  jl = jm
            else
                  ju = jm
            endif
            go to 10
        endif
        !
        ! Return with the array index:
        !
        j = jl
        return
    end subroutine locate

end module backtr