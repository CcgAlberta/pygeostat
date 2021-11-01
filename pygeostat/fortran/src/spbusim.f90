! spbusim implements Semiparametric Bayesian Updating for Nonlinear Inference
!
! Compile as a python DLL with the command:
!    call f2py -c -m --fcompiler=gnu95 spbusim random.f90 covasubs.for sortem.f90 normaldist.f90 gausskde.f90 linearinterp.f90 spbusim.f90 ../../resource/lapack_solve.a
!
!-----------------------------------------------------------------------------
! Copyright (c) 2014-2015, Jared L. Deutsch
! All rights reserved.
!
! Redistribution and use in source and binary forms, with or without modification,
! are permitted provided that the following conditions are met:
!
! 1. Redistributions of source code must retain the above copyright notice,
!    this list of conditions and the following disclaimer.
!
! 2. Redistributions in binary form must reproduce the above copyright notice,
!    this list of conditions and the following disclaimer in the documentation
!    and/or other materials provided with the distribution.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
! "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
! TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
! PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS
! BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
! CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
! SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
! INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
! CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
! ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
! THE POSSIBILITY OF SUCH DAMAGE.
!-----------------------------------------------------------------------------


  subroutine ngspbusim(xzbiv,nxdisc,xmin,xsiz,nzdisc,zmin,zsiz, &
    xsec,compxyz,ncomps,nreal,rseed,tmin, &
    nst,c0,it,cc,azm,dip,tilt,ahmax,ahmin,avert,zsim)
!-----------------------------------------------------------------------
!
!  ngspbusim - Semiparametric Bayesian Updating for Nonlinear Inference
!  ********************************************************************
!
! This subroutine applies semiparametric Bayesian updating for nonlinear
! inference when the data are not nicely composited (ie: nested). It does
! not require a Gaussian transform of the bivariate distribution.
!
! 1. The prior distribution is inferred using multigaussian kriging of
!    previously simulated Z_v values
! 2. The likelihood distribution is inferred by calculating the conditional
!    distribution of Z_v|X_v=x with the colocated sample
! 3. The global distribution is calculated from the marginal of the input
!    bivariate distribution.
!
!  Parameters:
!    - Small 2D bivariate distribution
!       xzbiv(nxdisc*nzdisc) - gridded bivariate distribution data where
!                              Z and X are in original units.
!       nxdisc,xmin,xsiz - gridded x parameters for xzbiv
!       nzdisc,zmin,zsiz - gridded z parameters for xzbiv
!    - Small composite data (X_v) and locations to be simulated
!       xsec(ncomps) - secondary data - X_v
!       compxyz(3,ncomps) - x,y,z points corresponding to xsec
!       ncomps - number of small composites
!    - Simulation parameters
!       nreal - number of realizations of each comp
!       rseed - random number seed
!       tmin - lower trimming limit for secondary variables
!    - Variogram parameters (standardized)
!       nst,c0,it,cc,azm,dip,tilt,ahmax,ahmin,avert
!  Returns:
!    - Simulated Z values at composite locations
!       zsim(nreal*ncomps) - simulated Z values
!
! (C) Jared Deutsch 2015
!-----------------------------------------------------------------------
  use random ! Mersenne-Twister random number generator
  use normaldist, only : snorm_inv, snorm_pdf, norm_pdf, snorm_cdf ! Normal distribution routines
  use linearinterp, only: gridbilerp_d ! Bilinear interpolation on a grid
  implicit none

  real(kind=8), parameter :: EPSLON=1.0e-6
  real(kind=8), parameter :: BIGDBLE=1d21,SMALLDBLE=1d-6

  integer, intent(in) :: nxdisc,nzdisc,nst,nreal,ncomps,rseed
  real(kind=8), intent(in) :: xzbiv(nxdisc*nzdisc)
  real(kind=8), intent(in) :: xsec(ncomps)
  real(kind=8), intent(in) :: compxyz(3,ncomps)
  real(kind=8), intent(in) :: xmin,xsiz,zmin,zsiz,tmin
  integer, intent(in), dimension(nst) :: it
  real(kind=8), intent(in) :: c0
  real(kind=8), intent(in), dimension(nst) :: cc,azm,dip,tilt, &
                                              ahmax,ahmin,avert
  real(kind=8), dimension(nreal*ncomps), intent(out) :: zsim
  integer :: i,j,icomp,ireal,idisc,ixdisc,izdisc,idx
  real(kind=8), dimension(nzdisc) :: zdisc,globpdf,priorpdf,updatedpdf, &
                                     likelipdf
  real(kind=8), dimension(nzdisc+2) :: updatedcdf,updatedcdf_z,globcdf,globcdf_z, &
                                       priorcdf,priorcdf_z
  real(kind=8) :: zdiscmin,zdiscsiz,sumpdf
  integer, dimension(ncomps) :: ncomps_shuffled
  integer :: tmp_int,rand_int,test,nd
  real(kind=8), dimension(ncomps,ncomps) :: spatial_cov,lhs
  real(kind=8), dimension(ncomps) :: wts,rhs,simvector,simvectorns
  real(kind=8) :: prior_mean,prior_var,zval,powint,prior_stdev,globalpvalue, &
                  zval_lower,zval_upper,cdfval_lower,cdfval_upper
!
! GSLIB interface variables and setup
!
  integer :: ist
  integer :: MAXROT
  real(kind=8), allocatable :: rotmat(:,:,:)
  real(kind=8) :: cmax,maxcov,cova,x,y,z,c0gslib
  real(kind=8), dimension(nst) :: ccgslib,aagslib,anis1,anis2
! GSLIB interface allocation
  MAXROT = nst
  allocate(rotmat(MAXROT,3,3))
  c0gslib = c0
! GSLIB anisotropy ratios
  do ist=1,nst
    anis1(ist) = ahmin(ist)/max(ahmax(ist),SMALLDBLE)
    anis2(ist) = avert(ist)/max(ahmax(ist),SMALLDBLE)
    ccgslib(ist) = cc(ist)
    aagslib(ist) = ahmax(ist)
  end do
! GSLIB rotation matrix
  rotmat = 0
  do ist=1,nst
    call setrot(azm(ist),dip(ist),tilt(ist),anis1(ist),anis2(ist), &
      ist,MAXROT,rotmat)
  end do
! Maximum covariance for calculating variogram values from cova3
  call cova3(0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,1,nst,nst,c0gslib, &
                it,ccgslib(1:nst),aagslib,1,MAXROT,rotmat,cmax,maxcov)
! Check that the variogram is standardized
  if (maxcov .gt. 1.001 .or. maxcov .lt. 0.999) then
    write(*,*) 'ERROR: the input variogram must be standardized!',maxcov
    stop
  end if
! Initialize the random number generator
  call init_genrand(rseed)

!
! Z_v discretization
!
  zdisc(1) = zmin
  do idisc=2,nzdisc
    zdisc(idisc) = zdisc(idisc-1)+zsiz
  end do
!
! Calculate the marginal distribution f(Z_v)
!
  ! Initialize the global pdf to 0
  globpdf = 0d0
  ! For each discretized z value, integrate over the bivariate
  idx = 0
  sumpdf = 0
  do izdisc=1,nzdisc
    do ixdisc=1,nxdisc
      idx = idx + 1
      globpdf(izdisc) = globpdf(izdisc) + xzbiv(idx)
    end do
    sumpdf = sumpdf + globpdf(izdisc)
  end do
  do izdisc=1,nzdisc
    globpdf(izdisc) = globpdf(izdisc)/sumpdf
  end do
  ! Integrate for the global CDF
  call buildcdf(nzdisc, globpdf, zdisc, globcdf, globcdf_z)
!
! Pre-calculate the (ncomps x ncomps) spatial covariance matrix
!
  do i=1,ncomps
    do j=1,ncomps
      call cova3(compxyz(1,i),compxyz(2,i),compxyz(3,i), &
                 compxyz(1,j),compxyz(2,j),compxyz(3,j), &
                 1,nst,nst,c0gslib,it,ccgslib(1:nst),aagslib,1,MAXROT,rotmat,cmax,cova)
      spatial_cov(i,j) = cova
    end do
  end do
!
! For each realization
!
  do ireal=1,nreal
    ! Shuffle the simulation order using a Fisher-Yates shuffle
    do icomp=1,ncomps
      ncomps_shuffled(icomp) = icomp
    end do
    do icomp=ncomps,1,-1
      tmp_int = ncomps_shuffled(icomp)
      rand_int = floor(grnd()*ncomps)+1
      ncomps_shuffled(icomp) = ncomps_shuffled(rand_int)
      ncomps_shuffled(rand_int) = tmp_int
    end do
!
! For each small composite
!
    do tmp_int=1,ncomps
      icomp = ncomps_shuffled(tmp_int)
!
! Infer the prior distribution using previously simulated values
!
      if (tmp_int .eq. 1) then
        ! No prior distribution since we have no spatial distribution - set to global
        prior_mean = 0.0d0
        prior_var = 1.0d0
      else
        ! Build up local spatial covariance matrices
        nd = tmp_int-1
        do i=1,nd
          do j=1,nd
            lhs(i,j) = spatial_cov(ncomps_shuffled(i),ncomps_shuffled(j))
          end do
        end do
        do i=1,nd
          rhs(i) = spatial_cov(ncomps_shuffled(tmp_int),ncomps_shuffled(i))
        end do
        ! Solve the normal equations
        call solve(lhs(1:nd,1:nd), wts(1:nd), rhs(1:nd), nd, 1, test)
        ! Calculate the prior mean and variance
        prior_mean = 0.0d0
        prior_var = 1.0d0
        do i=1,nd
          ! Use the Gaussian transformed value for MG kriging
          prior_mean = prior_mean + wts(i)*simvectorns(ncomps_shuffled(i))
          prior_var = prior_var - wts(i)*rhs(i)
        end do
      end if
      ! Build a CDF with evenly discretized CDF values
      prior_stdev = dsqrt(prior_var)
      ! Endpoints
      priorcdf(1) = 0d0
      priorcdf_z(1) = globcdf_z(1)
      priorcdf(nzdisc+2) = 1d0
      priorcdf_z(nzdisc+2) = globcdf_z(nzdisc+2)
      ! Middle
      do izdisc=1,nzdisc
        priorcdf(izdisc+1) = dble(izdisc)/dble(nzdisc+1) !ie: priorcdf = 0.01,0.02,...,0.99
        globalpvalue = snorm_cdf(snorm_inv(priorcdf(izdisc+1))*prior_stdev+prior_mean)
        call locate(globcdf,nzdisc+2,1,nzdisc+2,globalpvalue,j)
        j = min(max(1,j),nzdisc+1)
        priorcdf_z(izdisc+1) = powint(globcdf(j),globcdf(j+1),globcdf_z(j),globcdf_z(j+1),globalpvalue,1.0d0)
      end do
      ! For each discrete Z interval - get the prior pdf
      do izdisc=1,nzdisc
        ! Z interval to calculate PDF on
        zval_lower = zdisc(izdisc) - 0.5d0*(zsiz)
        zval_upper = zdisc(izdisc) + 0.5d0*(zsiz)
        ! CDF values corresponding to Z interval
        call locate(priorcdf_z,nzdisc+2,1,nzdisc+2,zval_lower,j)
        j = min(max(1,j),nzdisc+1)
        cdfval_lower = powint(priorcdf_z(j),priorcdf_z(j+1),priorcdf(j),priorcdf(j+1),zval_lower,1.0d0)
        call locate(priorcdf_z,nzdisc+2,1,nzdisc+2,zval_upper,j)
        j = min(max(1,j),nzdisc+1)
        cdfval_upper = powint(priorcdf_z(j),priorcdf_z(j+1),priorcdf(j),priorcdf(j+1),zval_upper,1.0d0)
        ! Difference is the estimate of the prior PDF
        priorpdf(izdisc) = max(cdfval_upper-cdfval_lower,0d0)
      end do
      sumpdf = sum(priorpdf)
      do idisc=1,nzdisc
        priorpdf(idisc) = priorpdf(idisc)/sumpdf
      end do
!
! Infer the likelihood distribution using the bivariate distribution and secondary data
!
      if (xsec(icomp) .gt. tmin) then
        do idisc=1,nzdisc
          ! Bilinear interpolation of the bivariate distribution
          call gridbilerp_d(nxdisc,xsiz,nzdisc,zsiz,xzbiv,xsec(icomp)-xmin,zdisc(idisc)-zmin,likelipdf(idisc))
        end do
        sumpdf = sum(likelipdf)
        do idisc=1,nzdisc
          likelipdf(idisc) = likelipdf(idisc)/sumpdf
        end do
      else
        ! If secondary is missing, set to global
        do idisc=1,nzdisc
          likelipdf(idisc) = globpdf(idisc)
        end do
      end if
!
! Combine the distributions using Bayesian Updating
!
      updatedpdf = likelipdf*priorpdf/globpdf
      updatedpdf = updatedpdf/sum(updatedpdf)
      ! Integrate for the updated CDF
      call buildcdf(nzdisc, updatedpdf, zdisc, updatedcdf, updatedcdf_z)
      ! Draw a sample from the CDF
      zval = grnd()
      call locate(updatedcdf,nzdisc+2,1,nzdisc+2,zval,j)
      j = min(max(1,j),nzdisc+1)
      simvector(icomp) = powint(updatedcdf(j),updatedcdf(j+1),updatedcdf_z(j),updatedcdf_z(j+1),zval,1.0d0)
      ! Also normal score transform this value for use in the prior distribution
      zval = simvector(icomp)
      call locate(globcdf_z,nzdisc+2,1,nzdisc+2,zval,j)
      j = min(max(1,j),nzdisc+1)
      simvectorns(icomp) = powint(globcdf_z(j),globcdf_z(j+1),globcdf(j),globcdf(j+1),zval,1.0d0)
      simvectorns(icomp) = snorm_inv(simvectorns(icomp))
      ! Copy over the simulated value
      zsim((ireal-1)*ncomps+icomp) = simvector(icomp)
    end do
  end do

  deallocate(rotmat)
  return
  end subroutine ngspbusim





  subroutine ncspbusim(xzbiv,nxdisc,xmin,xsiz,nzdisc,zmin,zsiz, &
    xsec,compxyz,ncomps,nreal,rseed,tmin, &
    nst,c0,it,cc,azm,dip,tilt,ahmax,ahmin,avert,zsim)
!-----------------------------------------------------------------------
!
!  ncspbusim - Semiparametric Bayesian Updating for Nonlinear Inference
!  ********************************************************************
!
! This subroutine applies semiparametric Bayesian updating for nonlinear
! inference when the data are not nicely composited (ie: nested).
!
!  Parameters:
!    - Small 2D bivariate distribution
!       xzbiv(nxdisc*nzdisc) - gridded bivariate distribution data where
!                              at least Z is normal scored!
!       nxdisc,xmin,xsiz - gridded x parameters for xzbiv
!       nzdisc,zmin,zsiz - gridded z parameters for xzbiv
!    - Small composites data
!       xsec(ncomps) - secondary data - X_v
!       compxyz(3,ncomps) - x,y,z points corresponding to xsec
!       ncomps - number of small composites
!    - Simulation parameters
!       nreal - number of realizations of each comp
!       rseed - random number seed
!       tmin - lower trimming limit for secondary variables
!    - Variogram parameters (standardized)
!       nst,c0,it,cc,azm,dip,tilt,ahmax,ahmin,avert
!  Returns:
!    - Simulated Z values at composite locations
!       zsim(nreal*ncomps)
!
! (C) Jared Deutsch 2015
!-----------------------------------------------------------------------
  use random ! Mersenne-Twister random number generator
  use normaldist, only : norm_pdf, snorm_pdf ! Normal distribution routines
  use linearinterp, only: gridbilerp_d ! Bilinear interpolation on a grid
  implicit none

  real(kind=8), parameter :: EPSLON=1.0e-6
  real(kind=8), parameter :: BIGDBLE=1d21,SMALLDBLE=1d-6

  integer, intent(in) :: nxdisc,nzdisc,nst,nreal,ncomps,rseed
  real(kind=8), intent(in) :: xzbiv(nxdisc*nzdisc)
  real(kind=8), intent(in) :: xsec(ncomps)
  real(kind=8), intent(in) :: compxyz(3,ncomps)
  real(kind=8), intent(in) :: xmin,xsiz,zmin,zsiz,tmin
  integer, intent(in), dimension(nst) :: it
  real(kind=8), intent(in) :: c0
  real(kind=8), intent(in), dimension(nst) :: cc,azm,dip,tilt, &
                                              ahmax,ahmin,avert
  real(kind=8), dimension(nreal*ncomps), intent(out) :: zsim
  integer :: i,j,icomp,ireal,idisc
  real(kind=8), dimension(nzdisc) :: zdisc,globpdf,priorpdf,updatedpdf, &
                                     likelipdf
  real(kind=8), dimension(nzdisc+2) :: updatedcdf,updatedcdf_z
  real(kind=8) :: zdiscmin,zdiscsiz,sumpdf
  integer, dimension(ncomps) :: ncomps_shuffled
  integer :: tmp_int,rand_int,test,nd
  real(kind=8), dimension(ncomps,ncomps) :: spatial_cov,lhs
  real(kind=8), dimension(ncomps) :: wts,rhs,simvector
  real(kind=8) :: prior_mean,prior_var,zval,powint
!
! GSLIB interface variables and setup
!
  integer :: ist
  integer :: MAXROT
  real(kind=8), allocatable :: rotmat(:,:,:)
  real(kind=8) :: cmax,maxcov,cova,x,y,z,c0gslib
  real(kind=8), dimension(nst) :: ccgslib,aagslib,anis1,anis2
! GSLIB interface allocation
  MAXROT = nst
  allocate(rotmat(MAXROT,3,3))
  c0gslib = c0
! GSLIB anisotropy ratios
  do ist=1,nst
    anis1(ist) = ahmin(ist)/max(ahmax(ist),SMALLDBLE)
    anis2(ist) = avert(ist)/max(ahmax(ist),SMALLDBLE)
    ccgslib(ist) = cc(ist)
    aagslib(ist) = ahmax(ist)
  end do
! GSLIB rotation matrix
  rotmat = 0
  do ist=1,nst
    call setrot(azm(ist),dip(ist),tilt(ist),anis1(ist),anis2(ist), &
      ist,MAXROT,rotmat)
  end do
! Maximum covariance for calculating variogram values from cova3
  call cova3(0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,1,nst,nst,c0gslib, &
                it,ccgslib(1:nst),aagslib,1,MAXROT,rotmat,cmax,maxcov)
! Check that the variogram is standardized
  if (maxcov .gt. 1.001 .or. maxcov .lt. 0.999) then
    write(*,*) 'ERROR: the input variogram must be standardized!',maxcov
    stop
  end if
! Initialize the random number generator
  call init_genrand(rseed)

!
! Z_v discretization
!
  zdisc(1) = zmin
  do idisc=2,nzdisc
    zdisc(idisc) = zdisc(idisc-1)+zsiz
  end do
!
! Calculate global PDF for Z_v - it is a standard normal deviate
!
  sumpdf = 0
  do idisc=1,nzdisc
    globpdf(idisc) = snorm_pdf(zdisc(idisc))
    sumpdf = sumpdf + globpdf(idisc)
  end do
  do idisc=1,nzdisc
    globpdf(idisc) = globpdf(idisc)/sumpdf
  end do
!
! Pre-calculate the (ncomps x ncomps) spatial covariance matrix
!
  do i=1,ncomps
    do j=1,ncomps
      call cova3(compxyz(1,i),compxyz(2,i),compxyz(3,i), &
                 compxyz(1,j),compxyz(2,j),compxyz(3,j), &
                 1,nst,nst,c0gslib,it,ccgslib(1:nst),aagslib,1,MAXROT,rotmat,cmax,cova)
      spatial_cov(i,j) = cova
    end do
  end do
!
! For each realization
!
  do ireal=1,nreal
    ! Shuffle the simulation order using a Fisher-Yates shuffle
    do icomp=1,ncomps
      ncomps_shuffled(icomp) = icomp
    end do
    do icomp=ncomps,1,-1
      tmp_int = ncomps_shuffled(icomp)
      rand_int = floor(grnd()*ncomps)+1
      ncomps_shuffled(icomp) = ncomps_shuffled(rand_int)
      ncomps_shuffled(rand_int) = tmp_int
    end do
!
! For each small composite
!
    do tmp_int=1,ncomps
      icomp = ncomps_shuffled(tmp_int)
!
! Infer the prior distribution using previously simulated values
!
      if (tmp_int .eq. 1) then
        ! No prior distribution since we have no spatial distribution - set to global
        prior_mean = 0.0d0
        prior_var = 1.0d0
      else
        ! Build up local spatial covariance matrices
        nd = tmp_int-1
        do i=1,nd
          do j=1,nd
            lhs(i,j) = spatial_cov(ncomps_shuffled(i),ncomps_shuffled(j))
          end do
        end do
        do i=1,nd
          rhs(i) = spatial_cov(ncomps_shuffled(tmp_int),ncomps_shuffled(i))
        end do
        ! Solve the normal equations
        call solve(lhs(1:nd,1:nd), wts(1:nd), rhs(1:nd), nd, 1, test)
        ! Calculate the prior mean and variance
        prior_mean = 0.0d0
        prior_var = 1.0d0
        do i=1,nd
          prior_mean = prior_mean + wts(i)*simvector(ncomps_shuffled(i))
          prior_var = prior_var - wts(i)*rhs(i)
        end do
      end if
      ! Discretize the prior PDF
      do idisc=1,nzdisc
        priorpdf(idisc) = norm_pdf(zdisc(idisc), prior_mean, prior_var)
      end do
      sumpdf = sum(priorpdf)
      do idisc=1,nzdisc
        priorpdf(idisc) = priorpdf(idisc)/sumpdf
      end do
!
! Infer the likelihood distribution using the bivariate distribution and secondary data
!
      if (xsec(icomp) .gt. tmin) then
        do idisc=1,nzdisc
          ! Bilinear interpolation of the bivariate distribution
          call gridbilerp_d(nxdisc,xsiz,nzdisc,zsiz,xzbiv,xsec(icomp)-xmin,zdisc(idisc)-zmin,likelipdf(idisc))
        end do
        sumpdf = sum(likelipdf)
        do idisc=1,nzdisc
          likelipdf(idisc) = likelipdf(idisc)/sumpdf
        end do
      else
        ! If secondary is missing, set to global
        do idisc=1,nzdisc
          likelipdf(idisc) = globpdf(idisc)
        end do
      end if
!
! Combine the distributions using Bayesian Updating
!
      updatedpdf = likelipdf*priorpdf/globpdf
      updatedpdf = updatedpdf/sum(updatedpdf)
      ! Integrate for the updated CDF
      call buildcdf(nzdisc, updatedpdf, zdisc, updatedcdf, updatedcdf_z)
      ! Draw a sample from the CDF
      zval = grnd()
      call locate(updatedcdf,nzdisc+2,1,nzdisc+2,zval,j)
      j = min(max(1,j),nzdisc+1)
      simvector(icomp) = powint(updatedcdf(j),updatedcdf(j+1),updatedcdf_z(j),updatedcdf_z(j+1),zval,1.0d0)
      ! Copy over the simulated value
      zsim((ireal-1)*ncomps+icomp) = simvector(icomp)
    end do
  end do

  deallocate(rotmat)
  return
  end subroutine ncspbusim

  subroutine spbusim(xzbiv,wtsbiv,nbivdata,xsec,zcomp,compxyz, &
    nsmall,ncomps,ndisc,nreal,rseed,bandw, &
    nst,c0,it,cc,azm,dip,tilt,ahmax,ahmin,avert, &
    zsim,zcompout)
!-----------------------------------------------------------------------
!
!  spbusim - Semiparametric Bayesian Updating for Nonlinear Inference
!  ******************************************************************
!
! This
!
!  Parameters:
!    - Small 2D bivariate distribution data
!       xzbiv(nbivdata,2) - bivariate distribution data
!       wtsbiv(nbivdata) - weights for bivariate distribution data
!       nbivdata - number of points in the bivariate distribution
!    - Composites
!       xsec(nsmall,ncomps) - secondary data - X_v
!       zcomp(ncomps) - known composite samples - Z_V
!       compxyz(3,nsmall,ncomps) - x,y,z points corresponding to xsec
!    - Simulation parameters
!       ndisc - number of points to discretize Z_v from min(Z_v) to max(Z_v)
!       nreal - number of realizations of each comp
!       rseed - random number seed
!       bandw - Kernel Density Estimation bandwidth
!    - Variogram parameters (standard)
!       nst,c0,it,cc,azm,dip,tilt,ahmax,ahmin,avert
!
! (C) Jared Deutsch 2014
!-----------------------------------------------------------------------
  use sortem ! Quicksort module
  use random ! Mersenne-Twister random number generator
  use normaldist, only : snorm_inv, snorm_pdf, norm_pdf, snorm_cdf ! Normal distribution routines
  implicit none
  real(kind=8), parameter :: EPSLON=1.0e-6
  real(kind=8), parameter :: BIGDBLE=1d21,SMALLDBLE=1d-6
  real(kind=8), intent(in) :: xzbiv(nbivdata,2),wtsbiv(nbivdata)
  integer, intent(in) :: nbivdata,nst,ndisc,nreal,nsmall,ncomps,rseed
  real(kind=8), intent(in) :: xsec(nsmall,ncomps),zcomp(ncomps)
  real(kind=8), intent(in) :: compxyz(3,nsmall,ncomps)
  integer, intent(in), dimension(nst) :: it
  real(kind=8), intent(in) :: c0,bandw
  real(kind=8), intent(in), dimension(nst) :: cc,azm,dip,tilt, &
                                              ahmax,ahmin,avert
  real(kind=8), intent(out) :: zsim(nsmall,nreal*ncomps), &
                               zcompout(nreal*ncomps)

  integer :: i,j,icomp,ismall,ireal,idisc
  real(kind=8), allocatable :: nsxzbiv(:,:),globalcdf_z(:),globalcdf(:), &
                               nsxsec(:,:),updatedcdf(:),updatedcdf_z(:)
  real(kind=8), dimension(ndisc,nsmall) :: likelipdf
  real(kind=8), dimension(ndisc) :: zdisc,globpdf,priorpdf,updatedpdf
  real(kind=8) :: zdiscmin,zdiscsiz,sumpdf
  real(kind=8), dimension(2,2) :: cov,cov_inv
  real(kind=8), dimension(2) :: dloc
  integer, dimension(nsmall) :: nsmall_shuffled
  integer :: tmp_int,rand_int,test,nd
  real(kind=8), dimension(nsmall,nsmall) :: spatial_cov,lhs
  real(kind=8), dimension(nsmall) :: wts,rhs,simvector
  real(kind=8) :: prior_mean,prior_var,zval,powint

! GSLIB interface variables
  integer :: ist
  integer :: MAXROT
  real(kind=8), allocatable :: rotmat(:,:,:)
  real(kind=8) :: cmax,maxcov,cova,x,y,z,c0gslib
  real(kind=8), dimension(nst) :: ccgslib,aagslib,anis1,anis2
! GSLIB interface allocation
  MAXROT = nst
  allocate(rotmat(MAXROT,3,3))
  c0gslib = c0
! GSLIB anisotropy ratios
  do ist=1,nst
    anis1(ist) = ahmin(ist)/max(ahmax(ist),SMALLDBLE)
    anis2(ist) = avert(ist)/max(ahmax(ist),SMALLDBLE)
    ccgslib(ist) = cc(ist)
    aagslib(ist) = ahmax(ist)
  end do
! GSLIB rotation matrix
  rotmat = 0
  do ist=1,nst
    call setrot(azm(ist),dip(ist),tilt(ist),anis1(ist),anis2(ist), &
      ist,MAXROT,rotmat)
  end do
! Maximum covariance for calculating variogram values from cova3
  call cova3(0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,1,nst,nst,c0gslib, &
                it,ccgslib(1:nst),aagslib,1,MAXROT,rotmat,cmax,maxcov)

! Initialize the random number generator
  call init_genrand(rseed)


! Additional memory allocation
  !write(*,*) 'allocating additional memory'
  allocate(nsxzbiv(nbivdata,2),globalcdf_z(nbivdata+2),globalcdf(nbivdata+2), &
           nsxsec(nsmall,ncomps),updatedcdf(ndisc+2),updatedcdf_z(ndisc+2))

! Univariate normal score transform of Zv and Xv points for inferring
! normal score bivariate
  ! Note that Z is transformed after X so the variables can be re-used
  ! as there is no reason to keep the tables for X after transforming the secondaries
  call nscore(xzbiv(:,1),wtsbiv,nbivdata,.true.,EPSLON, &
    globalcdf_z,globalcdf,nsxzbiv(:,1))
  ! Normal score transform xsec secondary values
  do icomp=1,ncomps
    do ismall=1,nsmall
      zval = xsec(ismall,icomp)
      call locate(globalcdf_z,nbivdata+2,1,nbivdata+2,zval,j)
      j = min(max(1,j),nbivdata+1)
      nsxsec(ismall,icomp) = powint(globalcdf_z(j),globalcdf_z(j+1),globalcdf(j),globalcdf(j+1),zval,1.0d0)
      nsxsec(ismall,icomp) = snorm_inv(nsxsec(ismall,icomp))
    end do
  end do

  ! Normal score transform the Z values and keep the transformation table
  call nscore(xzbiv(:,2),wtsbiv,nbivdata,.true.,EPSLON, &
    globalcdf_z,globalcdf,nsxzbiv(:,2))
! Discretize global PDF for Z - now a standard normal deviate
  zdiscmin = minval(nsxzbiv(:,2))
  zdiscsiz = (maxval(nsxzbiv(:,2))-zdiscmin)/(ndisc-1)
  zdisc(1) = zdiscmin
  do idisc=2,ndisc
    zdisc(idisc) = zdisc(idisc-1)+zdiscsiz
  end do
  sumpdf = 0
  do idisc=1,ndisc
    globpdf(idisc) = snorm_pdf(zdisc(idisc))
    sumpdf = sumpdf + globpdf(idisc)
  end do
  do idisc=1,ndisc
    globpdf(idisc) = globpdf(idisc)/sumpdf
  end do
! Covariance matrix and inverse covariance of bivariate data
  call covariance(nsxzbiv,nbivdata,2,cov)
  cov_inv = cov
  !write(*,*) cov
  call invert(cov_inv,2)
  !write(*,*) cov_inv
! Loop over icomp=1,...,ncomps composites
  do icomp=1,ncomps
    !write(*,*) icomp
! Infer the likelihood distributions using Kernel Density Estimation
    !write(*,*) 'starting KDE'
    do ismall=1,nsmall
      ! Multivariate location - this is the secondary value X
      dloc(1) = nsxsec(ismall,icomp)
      sumpdf = 0.0d0
      do idisc=1,ndisc
        ! Multivariate location - this is the primary value Z
        dloc(2) = zdisc(idisc)
        call kernel_gauss(nsxzbiv,nbivdata,2,dloc,cov,cov_inv,bandw,likelipdf(idisc,ismall))
        sumpdf = sumpdf + likelipdf(idisc,ismall)
      end do
      if (sumpdf .lt. SMALLDBLE) then
        write(*,*) 'ERROR: No KDE possible for normal score secondary =',xsec(ismall,icomp)
        write(*,*) '       Assuming global distribution for likelihood!'
        write(*,*) 'If multiple errors of this type are encountered, the bandwidth may'
        write(*,*) 'be too small.'
        likelipdf(:,ismall) = globpdf
        sumpdf = sum(likelipdf(:,ismall))
      end if
      do idisc=1,ndisc
        likelipdf(idisc,ismall) = likelipdf(idisc,ismall)/sumpdf
      end do
    end do
    !write(*,*) 'pre-calculating spatial covariance matrix'
! Pre-calculate the (nsmall x nsmall) spatial covariance matrix
    do i=1,nsmall
      do j=1,nsmall
        call cova3(compxyz(1,i,icomp),compxyz(2,i,icomp),compxyz(3,i,icomp), &
                   compxyz(1,j,icomp),compxyz(2,j,icomp),compxyz(3,j,icomp), &
                   1,nst,nst,c0gslib,it,ccgslib(1:nst),aagslib,1,MAXROT,rotmat,cmax,cova)
        spatial_cov(i,j) = cova
      end do
    end do
    !write(*,*) 'looping over realizations'
! Loop over realizations ireal=1,...,nreal
    do ireal=1,nreal
! Loop over ismall=1,...,nsmall (SHUFFLED)
      ! Shuffle the array first using a Fisher-Yates shuffle
      do ismall=1,nsmall
        nsmall_shuffled(ismall) = ismall
      end do
      do ismall=nsmall,1,-1
        tmp_int = nsmall_shuffled(ismall)
        rand_int = floor(grnd()*nsmall)+1
        nsmall_shuffled(ismall) = nsmall_shuffled(rand_int)
        nsmall_shuffled(rand_int) = tmp_int
      end do
      do tmp_int=1,nsmall
        ismall = nsmall_shuffled(tmp_int)
        !write(*,*) 'tmp_int =',tmp_int
        !write(*,*) 'nsmall_shuffled(tmp_int) =',nsmall_shuffled(tmp_int)
! Infer the prior distribution using previously simulated values
        if (tmp_int .eq. 1) then
          ! No prior distribution since we have no spatial distribution - set to global
          prior_mean = 0.0d0
          prior_var = 1.0d0
        else
          ! Build up local spatial covariance matrices using pre-calculated values
          nd = tmp_int-1
          do i=1,nd
            do j=1,nd
              lhs(i,j) = spatial_cov(nsmall_shuffled(i),nsmall_shuffled(j))
            end do
          end do
          do i=1,nd
            rhs(i) = spatial_cov(nsmall_shuffled(tmp_int),nsmall_shuffled(i))
          end do
          ! Solve the normal equations
          call solve(lhs(1:nd,1:nd), wts(1:nd), rhs(1:nd), nd, 1, test)
          ! Calculate the prior mean and variance
          prior_mean = 0.0d0
          prior_var = 1.0d0
          do i=1,nd
            prior_mean = prior_mean + wts(i)*simvector(nsmall_shuffled(i))
            prior_var = prior_var - wts(i)*rhs(i)
          end do
        end if
        ! Discretize the prior PDF
        do idisc=1,ndisc
          priorpdf(idisc) = norm_pdf(zdisc(idisc), prior_mean, prior_var)
        end do
        sumpdf = sum(priorpdf)
        do idisc=1,ndisc
          priorpdf(idisc) = priorpdf(idisc)/sumpdf
        end do
! Combine the distributions using Bayesian Updating
        updatedpdf = likelipdf(:,ismall)*priorpdf/globpdf
        updatedpdf = updatedpdf/sum(updatedpdf)
! Integrate for the updated CDF
        call buildcdf(ndisc, updatedpdf, zdisc, updatedcdf, updatedcdf_z)
! Draw a sample from the CDF
        zval = grnd()
        call locate(updatedcdf,ndisc+2,1,ndisc+2,zval,j)
        j = min(max(1,j),ndisc+1)
        simvector(ismall) = powint(updatedcdf(j),updatedcdf(j+1),updatedcdf_z(j),updatedcdf_z(j+1),zval,1.0d0)
      end do
! Copy over the simulated vector and true values
      zsim(:,(ireal-1)*ncomps+icomp) = simvector
      zcompout((ireal-1)*ncomps+icomp) = zcomp(icomp)
    end do
  end do
! Back-transform all simulated values and return with the truth
  do icomp=1,ncomps
    do ireal=1,nreal
      do ismall=1,nsmall
        zval = snorm_cdf(zsim(ismall,(ireal-1)*ncomps+icomp))
        call locate(globalcdf,nbivdata+2,1,nbivdata+2,zval,j)
        j = min(max(1,j),nbivdata+1)
        zsim(ismall,(ireal-1)*ncomps+icomp) = powint(globalcdf(j),globalcdf(j+1),globalcdf_z(j),globalcdf_z(j+1),zval,1.0d0)
      end do
    end do
  end do
  return
  end subroutine spbusim

  subroutine nscore(zvar,weights,nvalues,rdespike,epslon,cdf_z,cdf,yvar)
!-----------------------------------------------------------------------
!
!                     nscore - Normal score transform
!                     *******************************
!
! This subroutine is based on the GSLIB program nscore. Weights, and
! random despiking can optionally be used. A CDF with values and
! quantile transformed Gaussian values are returned. No trimming
! is done here - trim before calling.
!
! (C) Jared Deutsch 2014
!-----------------------------------------------------------------------
  use sortem ! Quicksort module
  use random ! Mersenne-Twister random number generator
  use normaldist, only : snorm_inv ! Normal distribution routines
  implicit none
  integer, parameter :: rseed=31415
  integer, intent(in) :: nvalues
  real(kind=8), intent(in) :: zvar(nvalues)
  logical, intent(in) :: rdespike
  real(kind=8), intent(in) :: weights(nvalues)
  real(kind=8), intent(in) :: epslon
  real(kind=8), intent(out) :: cdf_z(nvalues+2),cdf(nvalues+2),yvar(nvalues)

  real(kind=8) :: vr(nvalues),wts(nvalues)
  real(kind=8) :: eps,totalwt,zval,powint
  integer :: i,j

  if (epslon .le. 0d0) then
    eps = 1.0d-6
  else
    eps = epslon
  end if
  if (rdespike) then
    ! Initialize the random number generator
    call init_genrand(rseed)
    ! Random despiking - add a small random component
    do i=1,nvalues
      vr(i) = zvar(i) + grnd()*eps
    end do
  else
    ! No random despiking - copy over the array
    vr = zvar
  end if
  ! Ensure weights sum to exactly 1.0
  totalwt = sum(weights)
  wts = weights/totalwt
  ! Sort the variable and weights
  call dblemodsortem(vr,nvalues,1,wts)
  ! Build the CDF
  call buildcdf(nvalues,wts,vr,cdf,cdf_z)
  ! Look up the Gaussian values corresponding to the CDF values
  ! for zvar
  do i=1,nvalues
    zval = zvar(i)
    call locate(cdf_z,nvalues+2,1,nvalues+2,zval,j)
    j = min(max(1,j),nvalues+1)
    yvar(i) = powint(cdf_z(j),cdf_z(j+1),cdf(j),cdf(j+1),zval,1.0d0)
    yvar(i) = snorm_inv(yvar(i))
  end do
  return
  end subroutine nscore
