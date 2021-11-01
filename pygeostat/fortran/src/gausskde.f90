! gausskde.f90 - contains routines for Gaussian kernel density estimation
!
! (c) 2013-2014, Ryan Barnett

  subroutine kernel_gauss(ovals,ndata,nvar,dloc,cov,cov_inv,bandw,pdf)
      implicit none
      integer, intent (in) :: ndata  !Number of homotopic samples
      integer, intent (in) :: nvar  !Number of homotopic variables
      real*8, intent (in) ::  ovals(ndata,nvar) !Homotopic samples
      real*8, intent (in) ::  dloc(nvar) ! Multi-dimensional location to have density estimated
      real*8, intent (in) :: cov(nvar,nvar) !Data covariance matrix
      real*8 , intent(in) :: cov_inv(nvar,nvar) !Inverse of data covariance matrix
      real*8, intent (in) :: bandw !Kernel bandwidth
      real*8, intent (out) :: pdf !Kernel density
      integer i,j,k
      real*8 dist(nvar,1),normm(1,1),norm,Cinv(nvar,nvar)
      real*8, parameter :: PI = 3.141592653589793238462643383279502884197
      Cinv=cov_inv/(2*bandw) !Scale bandwidth by the covariance inverse
      pdf=0 
      ! Loop over all data to calculate the summation of all kernels at the specified location
      do i = 1,ndata
          !Calculate the distance between the location and the 
          !datapoint 
          do j=1,nvar
              dist(j,1)=dloc(j)-ovals(i,j)
          enddo
          ! Convert to mahalanobis distance (cov_inv also considers bandwidth)
          normm=matmul(transpose(dist),matmul(Cinv,dist))
          norm=normm(1,1)
          pdf=pdf+exp(-norm)
      enddo
      ! Apply the licit 'front' term
      do i=1,nvar
          norm = sum( cov(i,:)**2 )
      enddo
      norm=sqrt(norm)
      pdf=pdf*(1/(2*PI)**real(nvar/2)*norm)
  end subroutine kernel_gauss
  
  subroutine buildcdf(nd,pdf,pdf_x,cdf,cdf_x)
      implicit none
      integer, intent (in) :: nd  !Number of pdf discretizations
      real*8, intent (in) ::  pdf(nd) !PDF values
      real*8, intent (in) ::  pdf_x(nd) !X values associated with PDF
      real*8, intent (out) ::  cdf(nd+2) !CDF values
      real*8, intent (out) ::  cdf_x(nd+2) !X values associated with cdf (not the same length as pdf_x due to tail values
      integer i,j,k
      real*8 cdf1(nd+2),slope
      !Build the cdf and extrapolate the tails
      cdf1(1) = 0.D0; cdf1(nd+2) = 1.D0
      do j = 2,nd+1
          cdf1(j) = cdf1(j-1) + pdf(j-1) 
          cdf_x(j) = pdf_x(j-1)
      enddo
      cdf = cdf1
      do j = 2,nd+1
          cdf(j) = ( cdf1(j) + cdf1(j-1) )/2.D0
      enddo
      ! Linear extrapolation for the lower tail value
      slope = ( cdf(3) - cdf(2) ) / ( cdf_x(3) - cdf_x(2) )
      if( slope < 0 )then
          cdf_x(1) = cdf_x(2) -  ( cdf(2) - cdf(1) ) / slope
      else
          cdf_x(1) = cdf_x(2) - ( cdf_x(3) - cdf_x(2) )
      endif
      ! Linear extrapolation for the upper tail value
      slope = ( cdf(nd+1) - cdf(nd) ) / ( cdf_x(nd+1) - cdf_x(nd) )
      if( slope < 0 )then
          cdf_x(nd+2) = cdf_x(nd+1) + ( cdf(nd+2) - cdf(nd+1) ) / slope
      else
          cdf_x(nd+2) = cdf_x(nd+1) + ( cdf_x(nd+1) - cdf_x(nd) )
      endif
  end subroutine buildcdf
  
  subroutine covariance(datavals,nd,nv,cov)
      implicit none
      integer, intent(in) :: nd,nv !Number of data and variables
      real*8, intent(in) :: datavals(nd,nv) !Data
      real*8, intent(out) :: cov(nv,nv) !Output covariance
      integer i,j,k
      cov=0
      do i = 1,nv-1
          cov(i,i) = 1.0
          do j=i+1,nv
              cov(i,j) = sum( datavals(:,i) * datavals(:,j) )
              cov(i,j) = cov(i,j) / real(nd)
              cov(j,i) = cov(i,j) 
          enddo    
      enddo
      cov(nv,nv) = 1.0
  end subroutine covariance

  real(kind=8) function getz(pval,nt,vr,cdf,zmin,zmax,ltail,ltpar,utail,utpar)
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
  real(kind=8), parameter :: EPSLON=1.0d-20
  integer, intent(in) :: nt
  real(kind=8), dimension(nt), intent(in) :: vr,cdf
  real(kind=8), intent(in) :: ltpar,utpar,pval,zmin,zmax
  integer, intent(in) :: ltail,utail
  real(kind=8) lambda,cdfhi,cdflo,powint,cpow
  integer j
!
! Value in the lower tail?    1=linear, 2=power, (3 and 4 are invalid):
!
  if(pval.le.cdf(1)) then
              getz = vr(1)
        if(ltail.eq.1) then
              getz = powint(0.0,cdf(1),zmin,vr(1),pval,1.0d0)
        else if(ltail.eq.2) then
              cpow = 1.0d0 / ltpar
              getz = powint(0.0,cdf(1),zmin,vr(1),pval,cpow)
        endif
!
! Value in the upper tail?     1=linear, 2=power, 4=hyperbolic:
!
  else if(pval.ge.cdf(nt)) then
              cdfhi  = cdf(nt)
              getz   = vr(nt)
        if(utail.eq.1) then
              getz   = powint(cdfhi,1.0d0,vr(nt),zmax,pval,1.0d0)
        else if(utail.eq.2) then
              cpow   = 1.0d0 / utpar
              getz   = powint(cdfhi,1.0d0,vr(nt),zmax,pval,cpow)
        else if(utail.eq.4) then
              lambda = (vr(nt)**utpar)*(1.0d0-cdf(nt))
              getz   = (lambda/(1.0d0-pval))**(1.0/utpar)
        endif
  else
!
! Value within the transformation table:
!
        call locate(cdf,nt,1,nt,pval,j)
        j    = max(min((nt-1),j),1)
        getz = powint(cdf(j),cdf(j+1),vr(j),vr(j+1),pval,1.0d0)
  endif
  if(getz.lt.zmin) getz = zmin
  if(getz.gt.zmax) getz = zmax
  return
  end function getz

  real(kind=8) function powint(xlow,xhigh,ylow,yhigh,xval,pow)
!-----------------------------------------------------------------------
!
! Power interpolate the value of y between (xlow,ylow) and (xhigh,yhigh)
!                 for a value of x and a power pow.
!
!-----------------------------------------------------------------------
  implicit none
  real(kind=8), parameter :: EPSLON=1.0d-20
  real(kind=8), intent(in) :: xlow,xhigh,ylow,yhigh,xval,pow

  if((xhigh-xlow).lt.EPSLON) then
        powint = (yhigh+ylow)/2.0d0
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
  real(kind=8), dimension(n), intent(in) :: xx
  real(kind=8), intent(in) :: x
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

