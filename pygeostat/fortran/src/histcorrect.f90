! histcorrect is designed to be a very fast histogram correction
! subroutine. It is a post-processing method to correct the results
! of simulations which may have variances that are too small, or
! means that are incorrect or not correctly distributed. Code is based
! off of the trans_trend program used under the GSLIB License. 
!
! Compile as a python DLL with the command:
!      f2py -c -m --fcompiler=gnu95 histcorrect random.f90 normaldist.f90 histcorrect.f90
!
! (c) 2014, Jared L. Deutsch
!
! Note: C does not support assumed shape arrays, so length arguments
! (old style fortran) are required, even if they are implicitly
! defined by the wrapper. This is required for a useful DLL that
! can be used in languages other than Fortran. 
!
! Another note: a bunch of other GSLIB subroutines are included here.
! This is temporary and marks the start of a conversion of these routines
! to a format useful from a DLL authoring perspective. 

  subroutine histcorrect(refvals,refwts,nref,simvals,nsim, &
                         usemean,globmean,usekv,kvs,omega, &
                         niter,nsub,rseed)
!-----------------------------------------------------------------------
!  
!             histcorrect - Realization Histogram Correction
!             **********************************************
!
!  Parameters:
!    - A reference distribution
!       refvals(nref) - values from the reference distribution
!       refwts(nref) - weights for the values, could be all 1s for
!                      equal weighting
!   PLANNED BUT NOT CURRENTLY IMPLEMENTED FEATURE:
!       If nref = 1, then refvals(1) will be treated as a Gaussian mean
!       and refwts(1) as a Gaussian variance. The correction will then
!       transform values to have a normal distribution with
!             corrected ~ N(refvals(1),refwts(1))
!    - A simulated distribution to correct
!       simvals(nsim) - simulated values to be corrected. All trimming
!                       must be done BEFORE calling histcorrect
!    - Flag for correcting to a specific mean value and mean value
!       usemean - logical for using global mean correction or not
!       globmean - global mean to correct to
!    - Flag for usage of kriging variance and kriging variances
!       usekv - logical for using kriging variance correction or not
!       kvs(nsim) - estimation (kriging variances)
!       omega - kriging variance control parameter (0.33 <= omega <= 2.0)
!       niter - number of iterations (only if using kriging variances)
!    - Subsample size for fast CDF determination
!       nsub - size of random sample to use
!              Note if <= 0 or >= nsim, no randomness
!    - Random number seed for mersenne twister random number generator
!       rseed - odd, integer random number seed used for random despiking
!               and subsampling. Random number generation could be replaced
!               with the current state of the random number generator
!               in a wrapper program. 
!  Required/Returns (modifies in place):
!    - simvals(nsim) - corrected simulation values
!
!-----------------------------------------------------------------------
  use random ! Mersenne-Twister random number generator
  use normaldist, only : snorm_inv,norm_pdf,snorm_pdf ! Normal distribution routines
  use sortem
  implicit none
  real, parameter :: EPSLON=1.0e-6  
! Subroutine parameters
  integer, intent(in) :: nref,nsim,niter,nsub,rseed
  logical, intent(in) :: usemean, usekv
  real, intent(in) :: omega,globmean
  real, dimension(nref), intent(in) :: refvals,refwts
  real, dimension(nsim), intent(in) :: kvs
  real, dimension(nsim), intent(inout) :: simvals
! Internal variables
  integer :: i,iter,err,cdfidx
  real :: targtotalwt,oldcp,cp,maxkv,zmin,zmax,refmean,step,znew,quantile,getz
  logical :: subsampling
  real, allocatable :: currvals(:),currcdf(:)
! Copies of values to modify here without changing the source
  real, dimension(nref) :: targvals,targcdf
  real, dimension(nsim) :: kvfactors
  integer :: nloops,ncurrcdf
!
! Initialize the random number generator
!
  call init_genrand(rseed)
  ! Subsampling?
  if (nsub.gt.2 .and. nsub.lt.nsim) then
    subsampling = .true.
    ncurrcdf = nsub
  else
    subsampling = .false.
    ncurrcdf = nsim
  end if
!
! Setup the target reference distribution
!
  if (nref .gt. 1) then
    ! Copy over the target reference values with despiking + weights
    do i=1,nref
      targvals(i) = refvals(i) + real(grnd())*EPSLON
    end do
    targcdf = refwts
    ! Sort the target reference values for the CDF
    call modsortem(targvals,nref,1,targcdf)
    ! Get the CDF values
    oldcp  = 0.0
    cp     = 0.0
    targtotalwt = sum(targcdf)
    refmean = 0.0
    do i=1,nref
      refmean = refmean + targvals(i)*targcdf(i)/targtotalwt
      cp         = cp + targcdf(i) / targtotalwt
      targcdf(i) = (cp + oldcp)  * 0.5
      oldcp      = cp
    end do
    zmin = minval(targvals)
    zmax = maxval(targvals)
!
! Modify the target distribution by the global mean - staying within zmin and zmax
!
    if (usemean) then
      refmean = globmean/refmean
      do i=1,nref
        targvals(i) = targvals(i)*refmean
        if (targvals(i) .lt. zmin) targvals(i) = zmin
        if (targvals(i) .gt. zmax) targvals(i) = zmax
      end do
    end if
  else
    ! nref == 1 so we are using the Gaussian distribution
    ! This is not implemented yet, return -999s
    simvals = -999.0
    return
  end if
!
! Establish the kriging variance factors
!
  if (usekv) then
    ! Max kriging variance?
    maxkv = maxval(kvs)
    ! Kriging variance factor
    do i=1,nsim
      if (kvs(i).ge.0.0) then
        kvfactors(i) = (kvs(i)/maxkv)**omega
      else
        kvfactors(i) = 1.0
      end if
    end do
    ! Keep iterations since we are using kriging variance
    nloops = niter
  else
    ! Kriging factors are all 1
    kvfactors = 1.0
    ! No iterations if not using kriging variance!
    nloops = 1
  end if
!
! Allocate space for building the current CDF
!
  allocate(currvals(ncurrcdf),stat=err)
  if (err /= 0) print *, "array: Allocation request denied"
  allocate(currcdf(ncurrcdf),stat=err)
  if (err /= 0) print *, "array: Allocation request denied"
  ! Set up the CDF values now - these won't change and just depend on the length
  step = 1.0/real(ncurrcdf)
  currcdf(1) = step*0.5
  do i=2,ncurrcdf
    currcdf(i) = currcdf(1) + step*(i-1)
  end do
!
! Main loop to correct the histogram (well, loop only if using kriging variances)
!
  do iter = 1,nloops
!
! Establish the current CDF for the histogram, possibly by subsampling
!
    if (subsampling) then
      ! Draw a random subsample - this is completely random and there is the
      ! chance that the same sample will be used multiple times...but it is FAST
      ! Make sure to use the random despiking
      do i=1,ncurrcdf
        currvals(i) = simvals(1+int(real(grnd())*(nsim-1)))+real(grnd())*EPSLON
      end do
    else
      ! Copy over the values
      do i=1,ncurrcdf
        currvals(i) = simvals(i)+real(grnd())*EPSLON
      end do
    end if
    ! Sort the CDF values
    call modsortem(currvals,ncurrcdf,0)
!
! Perform a quantile-quantile transform to correct the CDF
!
    do i=1,nsim
      !
      ! Get the quantile
      !
      call locate(currvals,ncurrcdf,1,ncurrcdf,simvals(i),cdfidx)
      if (subsampling) then
        ! Linearly interpolate for the current quantile
        if (cdfidx .gt. 0 .and. cdfidx .lt. ncurrcdf) then
          quantile = currcdf(cdfidx)+(simvals(i)-currvals(cdfidx))* &
                        (step)/(currvals(cdfidx+1)-currvals(cdfidx))
        elseif (cdfidx .le. 0) then
          if (currvals(1)-zmin .eq. 0.0) then
            quantile = 0
          else
            quantile = (simvals(i)-zmin)* &
                        (currcdf(1))/(currvals(1)-zmin)
          end if
        elseif (cdfidx .ge. ncurrcdf) then
          if (zmax-currvals(ncurrcdf) .eq. 0.0) then
            quantile = 1
          else
            quantile = currcdf(ncurrcdf)+(simvals(i)-currvals(ncurrcdf))* &
                        (currcdf(1))/(zmax-currvals(ncurrcdf))
          end if
        end if
      else
        quantile = currcdf(cdfidx)
      end if
      !
      ! Get the new Z value from the target histogram
      !
      if (quantile.eq.0) then
        znew = zmin
      elseif (quantile.eq.1) then
        znew = zmax
      else
        znew = getz(quantile,nref,targvals,targcdf,zmin,zmax,1,0.0,1,0.0)
      end if
      !
      ! Adjust considering the kriging variance factor
      !
      simvals(i) = simvals(i)+kvfactors(i)*(znew-simvals(i))
      !
      ! Check we are between the bounds and assign
      !
      if (simvals(i) .lt. zmin) simvals(i) = zmin
      if (simvals(i) .gt. zmax) simvals(i) = zmax
    end do
! Go back and do another iteration, if iterating
  end do
! Done - Fortran deallocates automatically
  return
  end subroutine histcorrect

! GSLIB subroutines are used under the license distributed with GSLIB:
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!                                                                      %
! Copyright (C) 1996, The Board of Trustees of the Leland Stanford     %
! Junior University.  All rights reserved.                             %
!                                                                      %
! The programs in GSLIB are distributed in the hope that they will be  %
! useful, but WITHOUT ANY WARRANTY.  No author or distributor accepts  %
! responsibility to anyone for the consequences of using them or for   %
! whether they serve any particular purpose or work at all, unless he  %
! says so in writing.  Everyone is granted permission to copy, modify  %
! and redistribute the programs in GSLIB, but only under the condition %
! that this notice and the above copyright notice remain intact.       %
!                                                                      %
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  real function powint(xlow,xhigh,ylow,yhigh,xval,pow)
!-----------------------------------------------------------------------
!
! Power interpolate the value of y between (xlow,ylow) and (xhigh,yhigh)
!                 for a value of x and a power pow.
!
!-----------------------------------------------------------------------
  implicit none
  real, parameter :: EPSLON=1.0e-20
  real, intent(in) :: xlow,xhigh,ylow,yhigh,xval,pow

  if((xhigh-xlow).lt.EPSLON) then
        powint = (yhigh+ylow)/2.0
  else
        powint = ylow + (yhigh-ylow)*(((xval-xlow)/(xhigh-xlow))**pow)
  end if

  return
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
  real, dimension(n), intent(in) :: xx
  real, intent(in) :: x
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
  end

      real function getz(pval,nt,vr,cdf,zmin,zmax,ltail,ltpar,utail,utpar)
!-----------------------------------------------------------------------
!
!           Back Transform Univariate Data from Normal Scores
!           *************************************************
!
! This subroutine backtransforms from a
! specified back transform table and option for the tails of the
! distribution.  
!
! INPUT VARIABLES:
!
!   pval             probability value to use
!   nt               number of values in the back transform tbale
!   vr(nt)           original data values that were transformed
!   cdf(nt)          the corresponding transformed values
!   zmin,zmax        limits possibly used for linear or power model
!   ltail            option to handle values less than cdf(1)
!   ltpar            parameter required for option ltail
!   utail            option to handle values greater than cdf(nt)
!   utpar            parameter required for option utail
!-----------------------------------------------------------------------
  implicit none
  real, parameter :: EPSLON=1.0e-20
  integer, intent(in) :: nt
  real, dimension(nt), intent(in) :: vr,cdf
  real, intent(in) :: ltpar,utpar,pval,zmin,zmax
  integer, intent(in) :: ltail,utail
  real lambda,cdfhi,cdflo,powint,cpow
  integer j
!
! Value in the lower tail?    1=linear, 2=power, (3 and 4 are invalid):
!
  if(pval.le.cdf(1)) then
              getz = vr(1)
        if(ltail.eq.1) then
              getz = powint(0.0,cdf(1),zmin,vr(1),pval,1.0)
        else if(ltail.eq.2) then
              cpow = 1.0 / ltpar
              getz = powint(0.0,cdf(1),zmin,vr(1),pval,cpow)
        endif
!
! Value in the upper tail?     1=linear, 2=power, 4=hyperbolic:
!
  else if(pval.ge.cdf(nt)) then
              cdfhi  = cdf(nt)
              getz   = vr(nt)
        if(utail.eq.1) then
              getz   = powint(cdfhi,1.0,vr(nt),zmax,pval,1.0)
        else if(utail.eq.2) then
              cpow   = 1.0 / utpar
              getz   = powint(cdfhi,1.0,vr(nt),zmax,pval,cpow)
        else if(utail.eq.4) then
              lambda = (vr(nt)**utpar)*(1.0-cdf(nt))
              getz   = (lambda/(1.0-pval))**(1.0/utpar)
        endif
  else
!
! Value within the transformation table:
!
        call locate(cdf,nt,1,nt,pval,j)
        j    = max(min((nt-1),j),1)
        getz = powint(cdf(j),cdf(j+1),vr(j),vr(j+1),pval,1.0)
  endif
  if(getz.lt.zmin) getz = zmin
  if(getz.gt.zmax) getz = zmax
  return
  end function getz

  subroutine gauinv(p,xp,ierr)
!-----------------------------------------------------------------------
!
! Computes the inverse of the standard normal cumulative distribution
! function with a numerical approximation from : Statistical Computing,
! by W.J. Kennedy, Jr. and James E. Gentle, 1980, p. 95.
!
!
! INPUT/OUTPUT:
!
!   p    = double precision cumulative probability value: dble(psingle)
!   xp   = G^-1 (p) in single precision
!   ierr = 1 - then error situation (p out of range), 0 - OK
!
!-----------------------------------------------------------------------
  real(kind=8), intent(in) :: p
  real, intent(out) :: xp
  integer, intent(out) :: ierr
  real(kind=8) p0,p1,p2,p3,p4,q0,q1,q2,q3,q4,y,pp,lim
  save   p0,p1,p2,p3,p4,q0,q1,q2,q3,q4,lim
!
! Coefficients of approximation:
!
  data lim/1.0e-20/
  data p0/-0.322232431088/,p1/-1.0/,p2/-0.342242088547/, &
       p3/-0.0204231210245/,p4/-0.0000453642210148/
  data q0/0.0993484626060/,q1/0.588581570495/,q2/0.531103462366/, &
       q3/0.103537752850/,q4/0.0038560700634/
!
! Check for an error situation:
!
  ierr = 1
  if(p.lt.lim) then
        xp = -1.0e21
        return
  end if
  if(p.gt.(1.0d0-lim)) then
        xp =  1.0e21
        return
  end if
  ierr = 0      
!
! Get k for an error situation:
!
  pp   = p
  if(p.gt.0.5d0) pp = 1 - pp
  xp   = 0.0
  if(p.eq.0.5d0) return
!
! Approximate the function:
!
  y  = dsqrt(dlog(1.0d0/(pp*pp)))
  xp = real( y + ((((y*p4+p3)*y+p2)*y+p1)*y+p0) / &
                ((((y*q4+q3)*y+q2)*y+q1)*y+q0) )
  if(real(p).eq.real(pp)) xp = -xp
!
! Return with G^-1(p):
!
  return
  end subroutine gauinv

  real function gcum(x)
!-----------------------------------------------------------------------
!
! Evaluate the standard normal cdf given a normal deviate x.  gcum is
! the area under a unit normal curve to the left of x.  The results are
! accurate only to about 5 decimal places.
!
!-----------------------------------------------------------------------
  real, intent(in) :: x
  real z
  z = x
  if(z.lt.0.) z = -z
  t    = 1./(1.+ 0.2316419*z)
  gcum = t*(0.31938153   + t*(-0.356563782 + t*(1.781477937 + &
         t*(-1.821255978 + t*1.330274429))))
  e2   = 0.
!
!  6 standard deviations out gets treated as infinity:
!
  if(z.le.6.) e2 = exp(-z*z/2.)*0.3989422803
  gcum = 1.0- e2 * gcum
  if(x.ge.0.) return
  gcum = 1.0 - gcum
  return
  end function gcum

  double precision function acorni(idum)
!-----------------------------------------------------------------------
!
! Fortran implementation of ACORN random number generator of order less
! than or equal to 12 (higher orders can be obtained by increasing the
! parameter value MAXORD).
!
! NOTES: 1. The variable idum is a dummy variable. The common block
!           IACO is used to transfer data into the function.
!
!        2. Before the first call to ACORN the common block IACO must
!           be initialised by the user, as follows. The values of
!           variables in the common block must not subsequently be
!           changed by the user.
!
!             KORDEI - order of generator required ( must be =< MAXORD)
!
!             MAXINT - modulus for generator, must be chosen small
!                      enough that 2*MAXINT does not overflow
!
!             ixv(1) - seed for random number generator
!                      require 0 < ixv(1) < MAXINT
!
!             (ixv(I+1),I=1,KORDEI)
!                    - KORDEI initial values for generator
!                      require 0 =< ixv(I+1) < MAXINT
!
!        3. After initialisation, each call to ACORN generates a single
!           random number between 0 and 1.
!
!        4. An example of suitable values for parameters is
!
!             KORDEI   = 10
!             MAXINT   = 2**30
!             ixv(1)   = an odd integer in the (approximate) range 
!                        (0.001 * MAXINT) to (0.999 * MAXINT)
!             ixv(I+1) = 0, I=1,KORDEI
!
! Author: R.S.Wikramaratna,                           Date: October 1990
!-----------------------------------------------------------------------
  implicit double precision (a-h,o-z)
  parameter (KORDEI=12,MAXOP1=KORDEI+1,MAXINT=2**30)
  common/iaco/ ixv(MAXOP1)
  do i=1,KORDEI
        ixv(i+1)=(ixv(i+1)+ixv(i))
        if(ixv(i+1).ge.MAXINT) ixv(i+1)=ixv(i+1)-MAXINT
  end do
  acorni=dble(ixv(KORDEI+1))/MAXINT
  return
  end
