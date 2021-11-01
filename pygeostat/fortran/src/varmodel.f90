! varmodel is a powerful, easy to use variogram modeling and fitting program.
!
! Compile as a python DLL with the command:
!  call f2py -c -m --fcompiler=gnu95 varmodel covasubs.for random.f90 varmodel.f90
!
!-----------------------------------------------------------------------------
! Copyright (c) 2014, Jared L. Deutsch, Matthew V. Deutsch and Clayton V. Deutsch
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

  subroutine varmodelfit(nst,c0,it,cc,azm,dip,tilt,ahmax,ahmin,avert,sill, &
    maxiter,rseed,minpairs,npairswt,invdistwt,fixhmaxvertanis,hmaxvertanis, &
    fixhminhmaxanis,hminhmaxanis,nvargs,varlagdist,varvalue,varnumpairs, &
    varazm,vardip,opt_c0,opt_it,opt_cc,opt_azm,opt_dip,opt_tilt, &
    opt_ahmax,opt_ahmin,opt_avert)
!-----------------------------------------------------------------------
!
!             varmodelfit - Semi-automatic variogram fitting
!             **********************************************
!
! This program fits 3D directional variograms given a set of constraints
! and partially optimized parameters.
!
!  Parameters:
!    - Partially optimized variogram parameters
!       nst - number of structures (cannot be optimized here, must be set)
!       c0(2) - min and max allowable values for the nugget effect
!       it(nst,2) - min and max structure types allowed
!       cc(nst,2) - min and max allowable values for structure contributions
!       azm(nst,2) - min and max allowable values for ang1/azimuth
!       dip(nst,2) - min and max allowable values for ang2/azimuth
!       tilt(nst,2) - min and max allowable values for ang3/azimuth
!       ahmax(nst,2) - min and max allowable values for major range
!       ahmin(nst,2) - min and max allowable values for minor range
!       avert(nst,2) - min and max allowable values for vertical range
!       sill(2) - min and max allowable sill value
!    - Optimization parameters
!       maxiter - maximum number of iterations
!       rseed - random number seed
!       minpairs - minimum number of pairs to consider
!       npairswt - logical indicating whether or not to optimize with # pairs
!       invdistwt - logical indicating whether or not to use inverse distance weighting
!       fixhmaxvertanis - logical indicating whether or not to fix hmax/vert anisotropy
!       hmaxvertanis - hmax/vert anisotropy
!       fixhminhmaxanis - logical indicating whether or not to fix hmin/hmax anisotropy
!       hminhmaxanis - hmin/hmax anisotropy
!    - Experimental variogram values
!       nvargs - number of experimental variogram points
!       varlagdist(nvargs) - lag distances
!       varvalue(nvargs) - variogram value
!       varnumpairs(nvargs) - number of pairs found for this variogram value
!       varazm(nvargs) - variogram azimuth
!       vardip(nvargs) - variogram dips
!    - Returns:
!       opt_c0 - nugget effect
!       opt_it(nst) - structure types
!       opt_cc(nst) - values for structure contributions
!       opt_azm(nst) - values for ang1/azimuth
!       opt_dip(nst) - values for ang2/azimuth
!       opt_tilt(nst) - values for ang3/azimuth
!       opt_ahmax(nst) - values for major range
!       opt_ahmin(nst) - values for minor range
!       opt_avert(nst) - values for vertical range
!
! (C) Jared Deutsch 2014
!-----------------------------------------------------------------------
  use random ! Mersenne-Twister random number generator
  implicit none
  real, parameter :: EPSLON=1.0e-6
  real(kind=8), parameter :: BIGDBLE=1d21,SMALLDBLE=1d-6,PI=4*atan(1.0d0)
  real(kind=8), parameter :: DEG2RAD=PI/180d0
! Parameters
  integer, intent(in) :: nst,maxiter,nvargs,minpairs,rseed
  real(kind=8), dimension(2), intent(in) :: c0,sill
  integer, dimension(nst,2), intent(in) :: it
  real(kind=8), dimension(nst,2), intent(in) :: cc,azm,dip,tilt,ahmax,ahmin,avert
  logical, intent(in) :: npairswt,invdistwt,fixhmaxvertanis,fixhminhmaxanis
  real(kind=8),intent(in) :: hmaxvertanis,hminhmaxanis
  integer, dimension(nvargs), intent(in) :: varnumpairs
  real(kind=8), dimension(nvargs), intent(in) :: varlagdist,varvalue,varazm,vardip
! Return variables
  real(kind=8), intent(out) :: opt_c0
  integer, dimension(nst), intent(out) :: opt_it
  real(kind=8), dimension(nst), intent(out) :: opt_cc,opt_azm,opt_dip,opt_tilt, &
    opt_ahmax,opt_ahmin,opt_avert
! Internal variables
  real(kind=8) :: old_c0,old_sill,old_obj,opt_obj
  integer, dimension(nst) :: old_it
  real(kind=8), dimension(nst) :: old_cc,old_azm,old_dip,old_tilt, &
    old_ahmax,old_ahmin,old_avert
  integer :: i,j,iter,ist,ivarg,nvargsused,ipert
  real(kind=8), dimension(nvargs) :: varwts,varfitvalue,suminvdist
  real(kind=8) :: npairstotal,rpert,randomshift
  real(kind=8) :: maxlagdist(3)
  logical, dimension(nvargs) :: usevar
  logical :: fixc0,fixsill,success
  logical, dimension(nst) :: fixcc
! Initialize the random number generator
  call init_genrand(rseed)
! Pre-calculate the sum inverse distance and total number of pairs
  npairstotal = 0
  suminvdist = SMALLDBLE
  do ivarg=1,nvargs
    if (varlagdist(ivarg) .gt. 0d0) then
      ! Only calculate sum inverse distance weight for the same direction
      do j=1,nvargs
        if ((varazm(ivarg) .eq. varazm(j)) .and. &
            (vardip(ivarg) .eq. vardip(j)) .and. &
            (varlagdist(j) .gt. 0d0)) then
          suminvdist(ivarg) = suminvdist(ivarg) + 1d0/varlagdist(j)
        end if
      end do
    end if
    if (varnumpairs(ivarg) .gt. 0) then
      npairstotal = npairstotal + varnumpairs(ivarg)
    end if
  end do
! Figure out what variograms we are using to speed up calculations
  nvargsused = 0
  maxlagdist = SMALLDBLE
  do ivarg=1,nvargs
    if (varnumpairs(ivarg) .ge. max(minpairs,1)) then
      usevar(ivarg) = .true.
      nvargsused = nvargsused + 1
    else
      usevar(ivarg) = .false.
    end if
  end do
! Have we got anything to work with?
  if (nvargsused .lt. 2) then
    write(*,*) 'Too few variogram values to even consider optimizing'
    write(*,*) 'Check that the minimum number of pairs specified is reasonable.'
    stop
  end if
! Weights for variogram points
  do ivarg=1,nvargs
    ! Initialize for equal weighting
    varwts(ivarg) = 1d0
    ! Inverse distance weighting of points
    if (invdistwt) then
      if (varlagdist(ivarg) .gt. 0d0) then
        varwts(ivarg) = varwts(ivarg)*1d0/varlagdist(ivarg)/suminvdist(ivarg)
      end if
    end if
    ! Is this a valid variogram point or are there too few pairs?
    if (varnumpairs(ivarg) .lt. max(minpairs,1)) then
      varwts(ivarg) = 0
    else
      ! Number of pairs weighting
      if (npairswt) then
        varwts(ivarg) = varwts(ivarg)*varnumpairs(ivarg)/npairstotal
      end if
    end if
  end do
!-----------------------------------------------------------------------
! Initial parameter estimates if not fixed
!-----------------------------------------------------------------------
! Sill
  if (sill(1) .ne. sill(2)) then
    ! Estimate a reasonable sill value to start with
    old_sill = 0
    do ivarg=1,nvargs
      if (usevar(ivarg)) then
        old_sill = old_sill + varvalue(ivarg)
      end if
    end do
    old_sill = old_sill/nvargsused
    ! If this is less than sill(1) and greater than sill(2), pick the closest
    if (old_sill .lt. sill(1)) old_sill = sill(1)
    if (old_sill .gt. sill(2)) old_sill = sill(2)
    fixc0 = .false.
  else
    old_sill = sill(1)
    fixsill = .true.
  end if
! Nugget effect
  if (c0(1) .ne. c0(2)) then
    opt_c0 = 0
    if (opt_c0 .lt. c0(1)) opt_c0 = c0(1)
    if (opt_c0 .gt. c0(2)) opt_c0 = c0(2)
    fixc0 = .false.
  else
    opt_c0 = c0(1)
    fixc0 = .true.
  end if
! Nested variance contributions
  ! Initial pass to determine fixed structures
  do ist=1,nst
    if (cc(ist,1) .ne. cc(ist,2)) then
      if (.not. fixsill) then
        ! No fixed sill, so set to fraction of variability
        opt_cc(ist) = old_sill/nst
      else
        opt_cc(ist) = 0
      end if
      fixcc(ist) = .false.
    else
      opt_cc(ist) = cc(ist,1)
      fixcc(ist) = .true.
    end if
  end do
  ! Fix the variance contributions (sills)
  call fixsills(nst,c0,sill,cc,fixc0,fixcc,opt_c0,opt_cc,success)
  if (.not. success) then
    write(*,*) 'WARNING: Unsuccessful sill optimization on first pass'
    write(*,*) '  Are your constraints consistent and possible?'
  end if
! Structure types
  do ist=1,nst
    if (it(ist,1) .ne. it(ist,2)) then
      ! Default to spherical, or lowest structure number
      if (1 .ge. it(ist,1)) then
        opt_it(ist) = 1
      else
        opt_it(ist) = it(ist,1)
      end if
    else
      opt_it(ist) = it(ist,1)
    end if
  end do
! Angles
  do ist=1,nst
    if (azm(ist,1) .ne. azm(ist,2)) then
      opt_azm = 0
      if (opt_azm(ist) .lt. azm(ist,1)) opt_azm(ist) = azm(ist,1)
      if (opt_azm(ist) .gt. azm(ist,2)) opt_azm(ist) = azm(ist,2)
    else
      opt_azm(ist) = azm(ist,1)
    end if
    if (dip(ist,1) .ne. dip(ist,2)) then
      opt_dip = 0
      if (opt_dip(ist) .lt. dip(ist,1)) opt_dip(ist) = dip(ist,1)
      if (opt_dip(ist) .gt. dip(ist,2)) opt_dip(ist) = dip(ist,2)
    else
      opt_dip(ist) = dip(ist,1)
    end if
    if (tilt(ist,1) .ne. tilt(ist,2)) then
      opt_tilt = 0
      if (opt_tilt(ist) .lt. tilt(ist,1)) opt_tilt(ist) = tilt(ist,1)
      if (opt_tilt(ist) .gt. tilt(ist,2)) opt_tilt(ist) = tilt(ist,2)
    else
      opt_tilt(ist) = tilt(ist,1)
    end if
  end do
! Ranges
  do ist=1,nst
    if (ahmax(ist,1) .ne. ahmax(ist,2)) then
      ! Find the combination of closest angle/furthest distnace
      do ivarg=1,nvargs
        if ((varazm(ivarg).eq.opt_azm(ist)) .and. &
            (vardip(ivarg).eq.opt_dip(ist))) then
          if (varlagdist(ivarg).gt.maxlagdist(1)) maxlagdist(1) = varlagdist(ivarg)
        end if
      end do
      if (maxlagdist(1).le.SMALLDBLE) then
        do ivarg=1,nvargs
          if (varlagdist(ivarg).gt.maxlagdist(1)) maxlagdist(1) = varlagdist(ivarg)
        end do
      end if
      opt_ahmax(ist) = real(ist)/real(nst)*maxlagdist(1)
      if (opt_ahmax(ist) .lt. ahmax(ist,1)) opt_ahmax(ist) = ahmax(ist,1)
      if (opt_ahmax(ist) .gt. ahmax(ist,2)) opt_ahmax(ist) = ahmax(ist,2)
    else
      opt_ahmax(ist) = ahmax(ist,1)
    end if
    if (ahmin(ist,1) .ne. ahmin(ist,2)) then
      ! Find the combination of closest angle/furthest distnace
      do ivarg=1,nvargs
        if ((varazm(ivarg)+90d0.eq.opt_azm(ist) .or. varazm(ivarg)-90d0.eq.opt_azm(ist)) .and. &
            (vardip(ivarg).eq.opt_dip(ist))) then
          if (varlagdist(ivarg).gt.maxlagdist(2)) maxlagdist(2) = varlagdist(ivarg)
        end if
      end do
      if (maxlagdist(2).le.SMALLDBLE) then
        do ivarg=1,nvargs
          if (varlagdist(ivarg).gt.maxlagdist(2)) maxlagdist(2) = varlagdist(ivarg)
        end do
      end if
      opt_ahmin(ist) = real(ist)/real(nst)*maxlagdist(2)
      if (opt_ahmin(ist) .lt. ahmin(ist,1)) opt_ahmin(ist) = ahmin(ist,1)
      if (opt_ahmin(ist) .gt. ahmin(ist,2)) opt_ahmin(ist) = ahmin(ist,2)
    else
      opt_ahmin(ist) = ahmin(ist,1)
    end if
    if (avert(ist,1) .ne. avert(ist,2)) then
      ! Find the combination of closest angle/furthest distnace
      do ivarg=1,nvargs
        if ((varazm(ivarg).eq.opt_azm(ist)) .and. &
            (vardip(ivarg)+90d0.eq.opt_dip(ist) .or. vardip(ivarg)-90d0.eq.opt_dip(ist))) then
          if (varlagdist(ivarg).gt.maxlagdist(3)) maxlagdist(3) = varlagdist(ivarg)
        end if
      end do
      if (maxlagdist(3).le.SMALLDBLE) then
        do ivarg=1,nvargs
          if (varlagdist(ivarg).gt.maxlagdist(3)) maxlagdist(3) = varlagdist(ivarg)
        end do
      end if
      opt_avert(ist) = real(ist)/real(nst)*maxlagdist(3)
      if (opt_avert(ist) .lt. avert(ist,1)) opt_avert(ist) = avert(ist,1)
      if (opt_avert(ist) .gt. avert(ist,2)) opt_avert(ist) = avert(ist,2)
    else
      opt_avert(ist) = avert(ist,1)
    end if
  end do
!-----------------------------------------------------------------------
! Optimization
!-----------------------------------------------------------------------
! Calculate initial objective function
  call varmodelpts(nst,opt_c0,opt_it,opt_cc,opt_azm,opt_dip, &
          opt_tilt,opt_ahmax,opt_ahmin,opt_avert, &
          nvargs,varlagdist,varazm,vardip,varfitvalue)
  opt_obj = 0
  do ivarg=1,nvargs
    opt_obj = opt_obj + varwts(ivarg)*(varfitvalue(ivarg)-varvalue(ivarg))**2
  end do
  write(*,*) 'Starting objective value = ',opt_obj
! Copy over the initial values
  old_obj = opt_obj
  old_c0 = opt_c0
  old_it = opt_it
  old_cc = opt_cc
  old_azm = opt_azm
  old_dip = opt_dip
  old_tilt = opt_tilt
  old_ahmax = opt_ahmax
  old_ahmin = opt_ahmin
  old_avert = opt_avert
  old_sill = opt_c0 + sum(opt_cc)
! Structure type, range, sill, angle optimization
  do iter=1,maxiter
! Draw a random value to pick what to change
    rpert = grnd()
! Nugget effect and variance contribution change
    if (rpert .lt. 0.35d0) then
      ipert = floor(grnd()*(nst+1)+1)
      ! Random shift is +/- 7.5% of old sill
      randomshift = 0.15d0*(grnd()-0.5d0)*old_sill
      if ((ipert .eq. nst+1 .and. .not. fixc0)) then
        ! Change the nugget effect
        opt_c0 = opt_c0+randomshift
      elseif (ipert .ne. nst+1) then
        if (.not. fixcc(ipert)) then
          opt_cc(ipert) = opt_cc(ipert)+randomshift
        end if
      end if
      ! Fix the variance contributions (sills)
      call fixsills(nst,c0,sill,cc,fixc0,fixcc,opt_c0,opt_cc,success)
    elseif ((rpert .ge. 0.35d0) .and. (rpert .lt. 0.75d0)) then
! Range change
      ! Pick a structure number
      ist = floor(grnd()*nst+1)
      ! Pick a range value - ahmax,ahmin or avert
      ipert = floor(grnd()*3+1)
      ! Random shift is +/- 7.5% of maximum range in this direction
      randomshift = 0.15d0*(grnd()-0.5d0)
      if (ipert.eq.3) then
        opt_avert(ist) = opt_avert(ist)+randomshift*maxval(opt_avert)
        if (ist.lt.nst) then
          if (opt_avert(ist).gt.opt_avert(ist+1)) opt_avert(ist) = opt_avert(ist+1)
        end if
        if (opt_avert(ist) .lt. avert(ist,1)) opt_avert(ist) = avert(ist,1)
        if (opt_avert(ist) .gt. avert(ist,2)) opt_avert(ist) = avert(ist,2)
        if (fixhmaxvertanis) then
          opt_ahmax(ist) = opt_avert(ist)*hmaxvertanis
        end if
      elseif (ipert.eq.2) then
        opt_ahmin(ist) = opt_ahmin(ist)+randomshift*maxval(opt_ahmin)
        if (ist.lt.nst) then
          if (opt_ahmin(ist).gt.opt_ahmin(ist+1)) opt_ahmin(ist) = opt_ahmin(ist+1)
        end if
        if (opt_ahmin(ist) .lt. ahmin(ist,1)) opt_ahmin(ist) = ahmin(ist,1)
        if (opt_ahmin(ist) .gt. ahmin(ist,2)) opt_ahmin(ist) = ahmin(ist,2)
        if (fixhminhmaxanis) then
          opt_ahmax(ist) = opt_ahmin(ist)/hminhmaxanis
        end if
      else
        opt_ahmax(ist) = opt_ahmax(ist)+randomshift*maxval(opt_ahmax)
        if (ist.lt.nst) then
          if (opt_ahmax(ist).gt.opt_ahmax(ist+1)) opt_ahmax(ist) = opt_ahmax(ist+1)
        end if
        if (opt_ahmax(ist) .lt. ahmax(ist,1)) opt_ahmax(ist) = ahmax(ist,1)
        if (opt_ahmax(ist) .gt. ahmax(ist,2)) opt_ahmax(ist) = ahmax(ist,2)
        if (fixhmaxvertanis) then
          opt_avert(ist) = opt_ahmax(ist)/hmaxvertanis
        end if
        if (fixhminhmaxanis) then
          opt_ahmin(ist) = opt_ahmax(ist)*hminhmaxanis
        end if
      end if
    elseif ((rpert .ge. 0.75d0) .and. (rpert .lt. 0.95d0)) then
! Angle change
      ! Pick a structure number
      ist = floor(grnd()*nst+1)
      ! Pick an angle - azm,dip,tilt
      ipert = floor(grnd()*3+1)
      ! Random shift is +/- 7.5 degrees for this direction
      randomshift = 15d0*(grnd()-0.5d0)
      if (ipert.eq.3) then
        opt_azm(ist) = opt_azm(ist)+randomshift
        if (opt_azm(ist) .lt. azm(ist,1)) opt_azm(ist) = azm(ist,1)
        if (opt_azm(ist) .gt. azm(ist,2)) opt_azm(ist) = azm(ist,2)
      elseif (ipert.eq.2) then
        opt_dip(ist) = opt_dip(ist)+randomshift
        if (opt_dip(ist) .lt. dip(ist,1)) opt_dip(ist) = dip(ist,1)
        if (opt_dip(ist) .gt. dip(ist,2)) opt_dip(ist) = dip(ist,2)
      else
        opt_tilt(ist) = opt_tilt(ist)+randomshift
        if (opt_tilt(ist) .lt. tilt(ist,1)) opt_tilt(ist) = tilt(ist,1)
        if (opt_tilt(ist) .gt. tilt(ist,2)) opt_tilt(ist) = tilt(ist,2)
        end if
    else
! Structure type change
      ist = floor(grnd()*nst+1)
      ipert = floor(grnd()*(it(ist,2)-it(ist,1))+it(ist,1))
      opt_it(ist) = ipert
    end if
! Calculate a new objective value
    call varmodelpts(nst,opt_c0,opt_it,opt_cc,opt_azm,opt_dip, &
          opt_tilt,opt_ahmax,opt_ahmin,opt_avert, &
          nvargs,varlagdist,varazm,vardip,varfitvalue)
    opt_obj = 0
    do ivarg=1,nvargs
      opt_obj = opt_obj + varwts(ivarg)*(varfitvalue(ivarg)-varvalue(ivarg))**2
    end do
    ! Accept?
    if (opt_obj .lt. old_obj) then
      ! Accept the iteration
      old_obj = opt_obj
      old_c0 = opt_c0
      old_it = opt_it
      old_cc = opt_cc
      old_azm = opt_azm
      old_dip = opt_dip
      old_tilt = opt_tilt
      old_ahmax = opt_ahmax
      old_ahmin = opt_ahmin
      old_avert = opt_avert
      old_sill = opt_c0 + sum(opt_cc)
    else
      ! Reject the iteration
      opt_obj = old_obj
      opt_c0 = old_c0
      opt_it = old_it
      opt_cc = old_cc
      opt_azm = old_azm
      opt_dip = old_dip
      opt_tilt = old_tilt
      opt_ahmax = old_ahmax
      opt_ahmin = old_ahmin
      opt_avert = old_avert
    end if
  end do
  write(*,*) 'Final objective value = ',opt_obj
! Final check to fix sills before returning
  call fixsills(nst,c0,sill,cc,fixc0,fixcc,opt_c0,opt_cc,success)
  if (.not. success) then
    write(*,*) 'ERROR: Unsuccessful sill optimization on final pass'
    write(*,*) '  Are your constraints consistent and possible?'
  end if
  return
  end subroutine varmodelfit

  subroutine varmodelpts(nst,c0,it,cc,azm,dip,tilt,ahmax,ahmin,avert, &
      nvargs,varlagdist,varazm,vardip,varmodelvals)
!-----------------------------------------------------------------------
! This subroutine calculates variogram model points at input lags with
! provided lag distances, azimuths and dips using the GSLIB cova3
! subroutine.
!
! (C) Jared Deutsch 2014
! NOTE : RM - 2019 - pygeostat cova3 modified for full double precision
!-----------------------------------------------------------------------
  implicit none
  real(kind=8), parameter :: SMALLDBLE=1d-6,PI=4*atan(1.0d0)
  real(kind=8), parameter :: DEG2RAD=PI/180d0
! Parameters
  integer, intent(in) :: nst,nvargs
  real(kind=8), intent(in) :: c0
  integer, dimension(nst), intent(in) :: it
  real(kind=8), dimension(nst), intent(in) :: cc,azm,dip,tilt,ahmax,ahmin,avert
  real(kind=8), dimension(nvargs), intent(in) :: varlagdist,varazm,vardip
! Return
  real(kind=8), dimension(nvargs), intent(out) :: varmodelvals
! Internal variables
  integer :: ist,ivarg,i,j
! GSLIB interface variables
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
! Calculate variogram values
  do ivarg=1,nvargs
    x = sin(DEG2RAD*varazm(ivarg))*cos(DEG2RAD*vardip(ivarg))*varlagdist(ivarg)
    y = cos(DEG2RAD*varazm(ivarg))*cos(DEG2RAD*vardip(ivarg))*varlagdist(ivarg)
    z = sin(DEG2RAD*vardip(ivarg))*varlagdist(ivarg)
    call cova3(0.0d0,0.0d0,0.0d0,x,y,z,1,nst,nst,c0gslib, &
                it,ccgslib(1:nst),aagslib,1,MAXROT,rotmat,cmax,cova)
    varmodelvals(ivarg) = maxcov - cova
  end do
  return
  end subroutine varmodelpts

  subroutine fixsills(nst,c0,sill,cc,fixc0,fixcc,opt_c0,opt_cc,success)
!-----------------------------------------------------------------------
! This subroutine enforces the constraints imposed by the user on sills
! if possible. As usual, garbage in = garbage out. No excessive error
! checking is performed, but maybe it should be.
!
! (C) Jared Deutsch 2014
!-----------------------------------------------------------------------
  implicit none
  integer, parameter :: maxiter = 100
  real(kind=8), parameter :: abstol = 1.0d-6
  integer, intent(in) :: nst
  real(kind=8), dimension(2), intent(in) :: c0,sill
  real(kind=8), dimension(nst,2), intent(in) :: cc
  logical, intent(in) :: fixc0
  logical, dimension(nst), intent(in) :: fixcc
  logical, intent(out) :: success
  real(kind=8), intent(inout) :: opt_c0
  real(kind=8), dimension(nst), intent(inout) :: opt_cc
  real(kind=8) :: unfixedvar,inputsill,inputunfixed,rtol,closestsill
  integer :: ist,niter
  logical sillsfixed
! Counter initialization
  sillsfixed = .false.
  rtol = abstol*sill(1)
  niter = 0
  success = .true.
! Main iterative loop to correct sills
  do while(.not. sillsfixed)
! Make sure that constraints are enforced and calculate input sill
    inputsill = 0d0
    if (opt_c0.lt.c0(1)) opt_c0 = c0(1)
    if (opt_c0.gt.c0(2)) opt_c0 = c0(2)
    inputsill = inputsill + opt_c0
    do ist=1,nst
      if (opt_cc(ist).lt.cc(ist,1)) opt_cc(ist) = cc(ist,1)
      if (opt_cc(ist).gt.cc(ist,2)) opt_cc(ist) = cc(ist,2)
      inputsill = inputsill + opt_cc(ist)
    end do
! Are we done?
    if ((sill(1) .eq. sill(2)) .and. &
        (abs(inputsill-sill(1)).lt.rtol)) then
      sillsfixed = .true.
    elseif ((inputsill .ge. sill(1)) .and. &
            (inputsill .le. sill(2))) then
      sillsfixed = .true.
    else
      ! Not done - try and fix - what sill constraint are we closest to?
      if (abs(sill(2)-inputsill).lt.abs(sill(1)-inputsill)) then
        closestsill = sill(2)
      else
        closestsill = sill(1)
      end if
      ! How much variability is not fixed?
      unfixedvar = closestsill
      inputunfixed = 0d0
      if (fixc0) then
        unfixedvar = unfixedvar - c0(1)
      else
        inputunfixed = inputunfixed + opt_c0
      end if
      do ist=1,nst
        if (fixcc(ist)) then
          unfixedvar = unfixedvar - cc(ist,1)
        else
          inputunfixed = inputunfixed + opt_cc(ist)
        end if
      end do
      ! Scale all the variable sills
      if (inputunfixed .gt. 0d0) then
        if (.not. fixc0) opt_c0 = opt_c0*unfixedvar/inputunfixed
        do ist = 1,nst
          if (.not. fixcc(ist)) then
            opt_cc(ist) = opt_cc(ist)*unfixedvar/inputunfixed
          end if
        end do
      end if
    end if
! Maximum iteration condition
    niter = niter + 1
    if (niter .gt. maxiter) then
      sillsfixed = .true.
      success = .false.
    end if
! End main iterative loop
  end do
! Seriously, double check these constraints again - these dominate
  if (opt_c0.lt.c0(1)) opt_c0 = c0(1)
  if (opt_c0.gt.c0(2)) opt_c0 = c0(2)
  do ist=1,nst
    if (opt_cc(ist).lt.cc(ist,1)) opt_cc(ist) = cc(ist,1)
    if (opt_cc(ist).gt.cc(ist,2)) opt_cc(ist) = cc(ist,2)
  end do
  return
  end subroutine fixsills
