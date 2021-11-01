! varcalc is a fast and flexible 3D experimental variogram calculation subroutine
! designed to replace the GSLIB gamv program. 
!
! Compile as a python DLL with the command:
!    call f2py -c -m --fcompiler=gnu95 varcalc sortem.f90 random.f90 varsubs.f90 varcalc.f90
!
!-----------------------------------------------------------------------------
! Copyright (c) 2014-2015, Jared L. Deutsch, Matthew V. Deutsch, and Clayton V. Deutsch
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

  subroutine varcalc(xyz,vals,ndata,nvars,tmin,tmax,strict, &
                invariotypes,variostd,variosills,varioheads,variotails,nvarios, &
                azm,dip,tilt,azmtol,diptol,nlags,lagdist,lagtol,bandhorz,bandvert, &
                variocutcat,varidx,varlagdist,varnumpairs,varvalue,varazm,vardip)
!-----------------------------------------------------------------------
!  
!             varcalc - Experimental variogram calculation
!             ********************************************
!
! This program calculates experimental variograms for 3D, irregularly
! spaced data using a set of directions, tolerances and variogram types. 
! The variogram calculation convention described in GSLIB (Geostatistical
! Software Library and User's Guide, Deutsch & Journel, 1998) is used. This
! program is designed to be used as a compiled library or wrapped with
! a parameter file typical of GSLIB programs. All parameters are described
! below:
!
!  Parameters:
!    - 3D data
!       xyz(ndata,3) - x,y,z coordinates
!       vals(ndata,nvars) - variable values
!       ndata - number of data
!       nvars - number of variables
!       tmin - lower trimming limit for variables
!       tmax - upper trimming limit for variables
!       strict - logical indicating whether to run some checks for common
!                errors. Setting to true will slightly slow execution but
!                could catch common user errors! Recommended to set to true.
!    - Variogram type parameters
!       variotypes(nvarios) - variogram types (see list below)
!       variostd(nvarios) - logical indicating standardize or not
!       variosills(nvarios) - variogram sill to standardize to if variostd
!       varioheads(nvarios) - index of head variogram variable
!       variotails(nvarios) - index of tail variogram variable
!       nvarios - number of variogram types to compute
!       variocutcat(nvarios,2) - cutoff/categories for variograms of type 9-12
!                                otherwise these values will be ignored
!    - Experimental variogram parameters
!       azm - azimuth angle measured positive clockwise from the +y axis (N-S)
!       dip - dip angle measured negative down from the xy plane along the azimuth
!       tilt - tilt (or plunge) angle measured positive counterclockwise about Y'
!       azmtol - azimuth angle tolerance
!       diptol - dip angle tolerance
!       nlags - number of experimental lags to compute
!       lagdist - experimental lag distance (h)
!       lagtol - experimental lag distance tolerance
!       bandhorz - horizontal bandwidth (perpendicular to azm-Z plane)
!       bandvert - vertical bandwidth (perpendicular to horizontal plane)       
!    - Returns:
!       varidx(nvarios*(nlags+1)) - variogram index (ie: i=1,...,nvarios)
!       varlagdist(nvarios*(nlags+1)) - average experimental lag distance
!       varnumpairs(nvarios*(nlags+1)) - number of pairs found
!       varvalue(nvarios*(nlags+1)) - variogram value
!       varazm(nvarios*(nlags+1)) = calculation azimuth direction
!       vardip(nvarios*(nlags+1)) = calculation dip direction
!
! Variogram types:
!   1 - traditional semivariogram
!   2 - cross semivariogram (full definition, see GSLIB)
!   3 - covariance (non-ergodic covariance)
!   4 - correlogram (standardized non-ergodic covariance)
!   5 - General relative variogram
!   6 - Pairwise relative variogram
!   7 - Semivariogram of logarithms
!   8 - Semimadogram
!   9 - Indicator semivariogram - continuous - with cutoff
!  10 - Indicator semivariogram - categorical - with category
!  11 - Indicator cross semivariogram - continuous - with 2 cutoffs
!  12 - Indicator cross semivariogram - categorical - with 2 categories
!
! (C) Jared Deutsch 2014
!-----------------------------------------------------------------------
  use sortem ! Quicksort module
  use varsubs ! Variogram subroutine module
  implicit none
  real, parameter :: EPSLON=1.0e-6
  real(kind=8), parameter :: BIGDBLE=1d21,SMALLDBLE=1d-6,PI=4*atan(1.0d0)
  real(kind=8), parameter :: DEG2RAD=PI/180d0
  integer, parameter :: lout=3
! Subroutine parameters
  integer, intent(in) :: ndata,nvars,nvarios,nlags
! Slicing (ndata,3) is fastest as :,i since the resulting slices are contiguous
! therefore xyz is constructed for rotations as:
  real(kind=8), dimension(ndata,3), intent(in)  :: xyz
  real(kind=8), dimension(ndata,nvars), intent(in)  :: vals
  real(kind=8), intent(in) :: tmin,tmax,azm,dip,tilt,diptol,lagdist,lagtol, &
                      azmtol,bandhorz,bandvert
  integer, dimension(nvarios), intent(in) :: invariotypes,varioheads,variotails
  real(kind=8), dimension(nvarios,2), intent(in) :: variocutcat
  real(kind=8), dimension(nvarios), intent(in) :: variosills
  logical, dimension(nvarios), intent(in) :: variostd
  logical, intent(in) :: strict
! Returns
  integer, dimension(nvarios*(nlags+1)), intent(out) :: varidx,varnumpairs
  real(kind=8), dimension(nvarios*(nlags+1)), intent(out) :: varlagdist,varvalue
  real(kind=8), dimension(nvarios*(nlags+1)), intent(out) :: varazm,vardip
! Internal variables
  logical :: omni
  real(kind=8), dimension(ndata,3) :: xyzrot
  real(kind=8), dimension(ndata) :: dblerotindex
  integer, dimension(ndata) :: introtindex
  real(kind=8), dimension(3,3) :: rotmat,forwardrotmat,reverserotmat
  real(kind=8) :: angle,azmtolrad,diptolrad,diptoldist,azmtoldist,farthest, &
                  ymin,ymax,zmin,zmax,xdist,h
  integer :: i,j,k,m,iv,outputidx,startidx
  real(kind=8), dimension(nvarios*(nlags+1)) :: headvariance,tailvariance, &
                                               headmean,tailmean
  logical, dimension(nvarios) :: flipped
  integer, dimension(nvarios) :: variotypes
! Initialize return variables
  varnumpairs = 0
  varlagdist = 0
  varvalue = 0
  flipped = .false.
! Initialize other accumulators
  headvariance = 0
  tailvariance = 0
  headmean = 0
  tailmean = 0
! Initialize flipped indicator - Used for variance - covariance, 1 - corr
  do iv=1,nvarios
    variotypes(iv) = invariotypes(iv)
    if (variotypes(iv).lt. 0) then
      variotypes(iv) = abs(variotypes(iv))
      flipped(iv) = .true.
    end if
  end do
! Initialize index tracking for values
  do i=1,ndata
    dblerotindex(i) = i
  end do
! If 'strict' is true, perform some checks for common errors
  if (strict) then
    ! Implement more of these - positivity constraints, weird cross etc
    do iv=1,nvarios
      if (variotypes(iv).ne.2 .and. variotypes(iv).ne.3 .and. variotypes(iv).lt.11) then
        ! Check that the head and tails match
        if (varioheads(iv).ne.variotails(iv)) then
          write(*,*) 'WARNING: you specified variogram type',variotypes(iv)
          write(*,*) '         but the head and tail variables differ!'
          write(*,*) '         This could lead to strange errors, such as'
          write(*,*) '         including values that should be trimmed and'
          write(*,*) '         meaningless values. '
        end if
      end if
    end do
  end if
! Correct the angle tolerances to lie between 0 and 90 degrees
  omni = .false.
  if (azmtol .ge. 90d0) then
    azmtolrad = 90d0*DEG2RAD
    omni = .true.
  elseif (azmtol .lt. 0d0) then
    azmtolrad =  0d0
  else
    azmtolrad = azmtol*DEG2RAD
  end if
  if (diptol .ge. 90d0) then
    diptolrad = 90d0*DEG2RAD
    omni = .true.
  elseif (diptol .lt. 0d0) then
    diptolrad =  0d0
  else
    diptolrad = diptol*DEG2RAD
  end if
! Precalculate the point where the angle tolerances no longer come in to play
  if (azmtol .gt. 90d0-SMALLDBLE) then
    azmtoldist = 0d0
  elseif (azmtol .lt. SMALLDBLE) then
    azmtoldist = BIGDBLE
  else
    azmtoldist = bandhorz/tan(azmtolrad)
  end if
  if (diptol .gt. 90d0-SMALLDBLE) then
    diptoldist = 0d0
  elseif (diptol .lt. SMALLDBLE) then
    diptoldist = BIGDBLE
  else
    diptoldist = bandvert/tan(diptolrad)
  end if
! Set up the rotation matrix using the GSLIB angles
  call setrotmat(azm,dip,tilt,forwardrotmat,reverserotmat)
! Rotate to align the new x' axis with the variogram direction
  xyzrot = matmul(xyz,forwardrotmat)
! Sort along the direction vector
  call dblemodsortem(xyzrot(:,1),ndata,3,xyzrot(:,2),xyzrot(:,3),dblerotindex)
! Track the sorted indices for computing values
  do i=1,ndata
    introtindex(i) = int(dblerotindex(i))
  end do
! Main variogram loop over all points indexed by i
  do i=1,ndata-1   ! This point, i, is the "tail" variable in GSLIB notation
    ! Farthest point to consider from the current location
    farthest = xyzrot(i,1) + nlags*lagdist + lagtol
    ! Pre-calculated bandwidth min/max y's and z's
    ymax = xyzrot(i,2) + bandhorz
    ymin = xyzrot(i,2) - bandhorz
    zmax = xyzrot(i,3) + bandvert
    zmin = xyzrot(i,3) - bandvert
    ! Loop over all possible pairs forward of the current point
    if (omni) then
      startidx = 1
    else
      startidx = i+1
    end if
    do j=startidx,ndata   ! This point, j, is the "head" variable in GSLIB notation
      ! If we are too far away, skip checking other pairs
      if (xyzrot(j,1) .gt. farthest) exit
      ! Check horizontal bandwidths
      if (xyzrot(j,2) .gt. ymax) cycle
      if (xyzrot(j,2) .lt. ymin) cycle
      ! Check vertical bandwidths
      if (xyzrot(j,3) .gt. zmax) cycle
      if (xyzrot(j,3) .lt. zmin) cycle
      ! Check if the points are at the same point on the x' axis
      xdist = abs(xyzrot(j,1)-xyzrot(i,1))
      if (i.eq.j) then
        cycle
      else
        ! Check azimuth tolerance if it matters at this distance
        if (xdist .lt. azmtoldist) then
          ! Check absolute azimuth tolerance
          angle = abs(atan((xyzrot(j,2)-xyzrot(i,2))/(max(xdist,SMALLDBLE))))
          if (angle .gt. azmtolrad) cycle
        end if
        ! Check dip tolerance if it matters at this distance
        if (xdist .lt. diptoldist) then
          ! Check absolute dip tolerance
          angle = abs(atan((xyzrot(j,3)-xyzrot(i,3))/(max(xdist,SMALLDBLE))))
          if (angle .gt. diptolrad) cycle
        end if
      end if
      ! At this point all tests have been passed so place the pair in the
      ! appropriate bin - due to possible overlapping bins each bin
      ! should be checked to see if pairs fall in it
      ! Find the distance between the values
      h = distance(xyzrot(i,:),xyzrot(j,:))
      do k=0,nlags
        ! Are we inside this bin
        if ((h .ge. (k*lagdist-lagtol)) .and. &
            (h .le. (k*lagdist+lagtol))) then
          ! We are inside the bin, so calculate each variogram value
          do iv=1,nvarios
            ! Are both values licit?
            if (variotypes(iv) .ne. 3) then
                if (vals(introtindex(i),varioheads(iv)) .le. tmin .or. &
                    vals(introtindex(i),varioheads(iv)) .gt. tmax) cycle
                if (vals(introtindex(j),variotails(iv)) .le. tmin .or. &
                    vals(introtindex(j),variotails(iv)) .gt. tmax) cycle
            else
                if (vals(introtindex(j),varioheads(iv)) .le. tmin .or. &
                  vals(introtindex(j),varioheads(iv)) .gt. tmax) cycle
                if (vals(introtindex(i),variotails(iv)) .le. tmin .or. &
                  vals(introtindex(i),variotails(iv)) .gt. tmax) cycle
            end if
            ! Check other pairs for cross variogram as well 2, 11, 12
            if (variotypes(iv) .eq. 2 .or. variotypes(iv) .ge. 11) then
              if (vals(introtindex(j),varioheads(iv)) .le. tmin .or. &
                  vals(introtindex(j),varioheads(iv)) .gt. tmax) cycle
              if (vals(introtindex(i),variotails(iv)) .le. tmin .or. &
                  vals(introtindex(i),variotails(iv)) .gt. tmax) cycle
            end if
            ! Not trimming, calculate and add the values on
            outputidx = (iv-1)*(nlags+1)+k+1
            varnumpairs(outputidx) = varnumpairs(outputidx) + 1
            varlagdist(outputidx) = varlagdist(outputidx) + h
            if (variotypes(iv) .ne. 2 .and. variotypes(iv) .lt. 9) then
              varvalue(outputidx) = varvalue(outputidx) + &
                  expgam(vals(introtindex(i),variotails(iv)), &
                         vals(introtindex(j),varioheads(iv)), &
                         variotypes(iv))
            elseif (variotypes(iv) .eq. 2) then ! Cross
              varvalue(outputidx) = varvalue(outputidx) + &
                  expgam(vals(introtindex(i),variotails(iv)), &
                         vals(introtindex(j),variotails(iv)), &
                         variotypes(iv), &
                         vals(introtindex(i),varioheads(iv)), &
                         vals(introtindex(j),varioheads(iv)))
            elseif (variotypes(iv) .ge. 9) then ! Indicator + indicator cross
              varvalue(outputidx) = varvalue(outputidx) + &
                  expgam(vals(introtindex(i),variotails(iv)), &
                         vals(introtindex(j),variotails(iv)), &
                         variotypes(iv), &
                         vals(introtindex(i),varioheads(iv)), &
                         vals(introtindex(j),varioheads(iv)), &
                         variocutcat(iv,1),variocutcat(iv,2))
            end if
            ! Should we accumulate the tail and head means?
            if ((variotypes(iv) .eq. 3) .or. &
                (variotypes(iv) .eq. 4) .or. &
                (variotypes(iv) .eq. 5)) then
              if (vals(introtindex(i),variotails(iv)) .gt. tmin .and. &
                    vals(introtindex(i),variotails(iv)) .le. tmax) then
                tailmean(outputidx) = tailmean(outputidx) + &
                            vals(introtindex(i),variotails(iv))
              end if
              if (vals(introtindex(j),varioheads(iv)) .gt. tmin .and. &
                    vals(introtindex(j),varioheads(iv)) .le. tmax) then
                headmean(outputidx) = headmean(outputidx) + &
                            vals(introtindex(j),varioheads(iv))
              end if
            end if
            ! Should we accumulate the tail and head variances?
            if (variotypes(iv) .eq. 4) then
              tailvariance(outputidx) = tailvariance(outputidx) + &
                            vals(introtindex(i),variotails(iv))**2
              headvariance(outputidx) = headvariance(outputidx) + &
                            vals(introtindex(j),varioheads(iv))**2
            end if
          end do
        end if
      end do
    end do
  end do
! Compute the average h distance and appropriate final variogram measure
  do k=0,nlags
    do iv=1,nvarios
      outputidx = (iv-1)*(nlags+1)+k+1
      varidx(outputidx) = iv
      varazm(outputidx) = azm
      vardip(outputidx) = dip
      ! Compute the variogram measures if any pairs were found, otherwise skip
      if (varnumpairs(outputidx) .gt. 0) then
        ! Average distance
        varlagdist(outputidx) = varlagdist(outputidx)/varnumpairs(outputidx)
        ! Compute the appropriate variogram measure
        varvalue(outputidx) = varvalue(outputidx)/varnumpairs(outputidx)
        ! Standardize the tail and head mean and variance if calculated
        if ((variotypes(iv) .eq. 3) .or. &
            (variotypes(iv) .eq. 4) .or. &
            (variotypes(iv) .eq. 5)) then
          tailmean(outputidx) = tailmean(outputidx)/varnumpairs(outputidx)
          headmean(outputidx) = headmean(outputidx)/varnumpairs(outputidx)
        end if
        if (variotypes(iv) .eq. 4) then
          tailvariance(outputidx) = tailvariance(outputidx)/varnumpairs(outputidx) &
                                      - tailmean(outputidx)**2
          headvariance(outputidx) = headvariance(outputidx)/varnumpairs(outputidx) &
                                      - headmean(outputidx)**2
        end if
        ! Compute the final measure
        select case (variotypes(iv))
          case(1) ! Traditional semivariogram
            varvalue(outputidx) = 0.5d0*varvalue(outputidx)
          case(2) ! Traditional cross semivariogram
            varvalue(outputidx) = 0.5d0*varvalue(outputidx)
          case(3) ! Covariance (non-ergodic covariance)
            varvalue(outputidx) = varvalue(outputidx) - &
                                tailmean(outputidx)*headmean(outputidx)
          case(4) ! Correlogram
            if ((tailvariance(outputidx) .gt. 0d0) .and. & 
                        (headvariance(outputidx) .gt. 0d0)) then
              varvalue(outputidx) = (varvalue(outputidx) - &
                                tailmean(outputidx)*headmean(outputidx))/ &
                      (sqrt(tailvariance(outputidx))*sqrt(headvariance(outputidx)))
            else
              varvalue(outputidx) = -999d0
            end if
          case(5) ! General relative semivariogram
            varvalue(outputidx) = varvalue(outputidx)/ &
                                ((tailmean(outputidx)+headmean(outputidx))/2)**2
          case(6) ! Pairwise relative semivariogram
            varvalue(outputidx) = 0.5d0*varvalue(outputidx)
          case(7) ! Semivariogram of logarithms
            varvalue(outputidx) = 0.5d0*varvalue(outputidx)
          case(8) ! Semimadogram
            varvalue(outputidx) = 0.5d0*varvalue(outputidx)
          case(9) ! Indicator semivariogram - continuous
            varvalue(outputidx) = 0.5d0*varvalue(outputidx)
          case(10) ! Indicator semivariogram - categorical
            varvalue(outputidx) = 0.5d0*varvalue(outputidx)
          case(11) ! Traditional cross indicator semivariogram - continuous
            write(*,*) 'WARNING: Cross indicator variogram results have not been verified!'
            varvalue(outputidx) = 0.5d0*varvalue(outputidx)
          case(12) ! Traditional cross indicator semivariogram - categorical
            write(*,*) 'WARNING: Cross indicator variogram results have not been verified!'
            varvalue(outputidx) = 0.5d0*varvalue(outputidx)
          case default
            write(*,*) 'ERROR: Invalid variogram type',variotypes(iv)
            stop
        end select
        ! 'Flip' the variogram if negative value was input
        if (flipped(iv)) then
          select case (variotypes(iv))
            case(3) ! Variance - covariance (non-ergodic covariance)
              if (varvalue(outputidx) .ne. -999d0) then
                varvalue(outputidx) = variosills(iv) - varvalue(outputidx)
              end if
            case(4) ! 1 - Correlogram
              if (varvalue(outputidx) .ne. -999d0) then
                varvalue(outputidx) = 1.0d0 - varvalue(outputidx)
              end if
            case default
              write(*,*) 'ERROR: Invalid flipped variogram type',-1*variotypes(iv)
              stop
          end select
        end if
        ! Standardize the variogram if standardizing
        if (variostd(iv)) then
          varvalue(outputidx) = varvalue(outputidx)/variosills(iv)
        end if
      else
        varlagdist(outputidx) = k*lagdist
        varvalue(outputidx) = -999d0
      end if       
    end do
  end do
  return
  end subroutine varcalc

  subroutine sillcalc(vals,ndata,nvars,tmin,tmax,strict,nvarios, &
                variotypes,varioheads,variotails,variocutcat,variosills)
!-----------------------------------------------------------------------
!  
!         sillcalc - Experimental variogram sill calculation
!         **************************************************
!
! This subroutine infers variogram sill values for different variogram
! types. 
!
!  Parameters:
!    - Data
!       vals(ndata,nvars) - variable values
!       ndata - number of data
!       nvars - number of variables
!       tmin - lower trimming limit for variables
!       tmax - upper trimming limit for variables
!       strict - logical indicating whether to run some checks for common
!                errors. Setting to true will slightly slow execution but
!                could catch common user errors! Recommended to set to true.
!    - Variogram type parameters
!       variotypes(nvarios) - variogram types (see list below)
!       varioheads(nvarios) - index of head variogram variable
!       variotails(nvarios) - index of tail variogram variable
!       variocutcat(nvarios,2) - cutoff/categories for variograms of type 9-12
!                                otherwise these values will be ignored
!       nvarios - number of variogram types to compute
!    - Returns
!       variosills(nvarios) - variogram sill to standardize to if variostd
!
! Variogram types:
!   1 - traditional semivariogram
!   2 - cross semivariogram (full definition, see GSLIB)
!   3 - covariance (non-ergodic covariance)
!   4 - correlogram (standardized non-ergodic covariance)
!   5 - General relative variogram
!   6 - Pairwise relative variogram
!   7 - Semivariogram of logarithms
!   8 - Semimadogram
!   9 - Indicator semivariogram - continuous - with cutoff
!  10 - Indicator semivariogram - categorical - with category
!  11 - Indicator cross semivariogram - continuous - with 2 cutoffs
!  12 - Indicator cross semivariogram - categorical - with 2 categories
!
! (C) Jared Deutsch 2014
!  Modified by Ryan Barnett - added a sill calculation for 
!     variotype = 10 (categorical) - Jan. 2016
!  Modified by Jared Deutsch - added a sill calculation for 
!     variotype = 9 (continuous categorical) - Aug. 2016
!  Modified by Jared Deutsch - added a sill calculation for 
!     variotype = 6 (pairwise relative) - Nov. 2016
!-----------------------------------------------------------------------
  use varsubs ! Variogram subroutine module
  use random ! Random sampling module, only used for pairwise relative sill
  implicit none
  real, parameter :: EPSLON=1.0e-6
  integer, parameter :: RSEED=31415, NPRSAMPLE=50000
  real(kind=8), parameter :: BIGDBLE=1d21,SMALLDBLE=1d-6,PI=4*atan(1.0d0)
  real(kind=8), parameter :: DEG2RAD=PI/180d0
  integer, parameter :: lout=3
! Subroutine parameters
  integer, intent(in) :: ndata,nvars,nvarios
  real(kind=8), dimension(ndata,nvars), intent(in)  :: vals
  real(kind=8), intent(in) :: tmin,tmax
  integer, dimension(nvarios), intent(in) :: variotypes,varioheads,variotails
  real(kind=8), dimension(nvarios,2), intent(in) :: variocutcat
  logical, intent(in) :: strict
  real(kind=8), dimension(nvarios), intent(out) :: variosills
! Internal variables
  integer :: iv,i,j,ndraw
  real(kind=8), dimension(nvarios,2) :: means,sumsqs
  integer, dimension(nvarios,2) :: nvalues
  real(kind=8), dimension(ndata) :: prob
  real(kind=8) :: prop, val1, val2
! Initialization
  variosills = 0
  means = 0
  sumsqs = 0
  nvalues = 0
! If 'strict' is true, perform some checks for common errors
  if (strict) then
    ! Reject types with no sill calculation implemented here
    do iv=1,nvarios
      if ( (variotypes(iv) .gt. 6 .and. variotypes(iv) .lt. 9) &
            .or. (variotypes(iv) .gt. 10) .or. (variotypes(iv) .eq. 5)) then
        write(*,*) 'WARNING: Variogram types with no implemented'
        write(*,*) '         sill calculation input!'
      end if
    end do
    ! Implement more checks!
  end if

  do iv=1,nvarios
! Sill of traditional semivariogram is the variance
    if (variotypes(iv).eq.1) then
      do i=1,ndata
        if ((vals(i,variotails(iv)).gt.tmin) .and. (vals(i,variotails(iv)).lt.tmax)) then
          means(iv,1) = means(iv,1) + vals(i,variotails(iv))
          sumsqs(iv,1) = sumsqs(iv,1) + vals(i,variotails(iv))*vals(i,variotails(iv))
          nvalues(iv,1) = nvalues(iv,1) + 1
        end if
      end do
      if (nvalues(iv,1).gt.0) then
        means(iv,1) = means(iv,1)/nvalues(iv,1)
        sumsqs(iv,1) = sumsqs(iv,1)/nvalues(iv,1)
        variosills(iv) = sumsqs(iv,1) - means(iv,1)*means(iv,1)
      else
        variosills(iv) = 0
      end if
    end if
! Sill of traditional cross semivariogram is the covariance
    if (variotypes(iv).eq.2) then
      do i=1,ndata
        if ((vals(i,variotails(iv)).gt.tmin) .and. (vals(i,variotails(iv)).lt.tmax) .and. &
            (vals(i,varioheads(iv)).gt.tmin) .and. (vals(i,varioheads(iv)).lt.tmax)) then
          means(iv,1) = means(iv,1) + vals(i,variotails(iv))
          means(iv,2) = means(iv,2) + vals(i,varioheads(iv))
          sumsqs(iv,1) = sumsqs(iv,1) + vals(i,variotails(iv))*vals(i,varioheads(iv))
          nvalues(iv,1) = nvalues(iv,1) + 1
        end if
      end do
      if (nvalues(iv,1).gt.0) then
        means(iv,1) = means(iv,1)/nvalues(iv,1)
        means(iv,2) = means(iv,2)/nvalues(iv,1)
        sumsqs(iv,1) = sumsqs(iv,1)/nvalues(iv,1)
        variosills(iv) = sumsqs(iv,1) - means(iv,1)*means(iv,2)
      else
        variosills(iv) = 0
      end if
    end if
! Sill of the covariance
    if (abs(variotypes(iv)).eq.3) then
      do i=1,ndata
        if ((vals(i,variotails(iv)).gt.tmin) .and. (vals(i,variotails(iv)).lt.tmax)) then
          means(iv,1) = means(iv,1) + vals(i,variotails(iv))
          sumsqs(iv,1) = sumsqs(iv,1) + vals(i,variotails(iv))*vals(i,variotails(iv))
          nvalues(iv,1) = nvalues(iv,1) + 1
        end if
        if ((vals(i,varioheads(iv)).gt.tmin) .and. (vals(i,varioheads(iv)).lt.tmax)) then
          means(iv,2) = means(iv,2) + vals(i,varioheads(iv))
          sumsqs(iv,2) = sumsqs(iv,2) + vals(i,varioheads(iv))*vals(i,varioheads(iv))
          nvalues(iv,2) = nvalues(iv,2) + 1
        end if
      end do
    if (nvalues(iv,1).gt.0 .and. nvalues(iv,2).gt.0) then
        means(iv,1) = means(iv,1)/nvalues(iv,1)
        means(iv,2) = means(iv,2)/nvalues(iv,2)
        sumsqs(iv,1) = sumsqs(iv,1)/nvalues(iv,1)
        sumsqs(iv,2) = sumsqs(iv,2)/nvalues(iv,2)
        variosills(iv) = sqrt(sumsqs(iv,1) - means(iv,1)*means(iv,1))*sqrt(sumsqs(iv,2) - means(iv,2)*means(iv,2))
      else
        variosills(iv) = 0
      end if
    end if
! Sill of the correlogram
    if (abs(variotypes(iv)).eq.4) then
      variosills(iv) = 1
    end if
! Sill of the pairwise relative variogram
    if (abs(variotypes(iv)).eq.6) then
      !Initialize the random number genersator
      call init_genrand(RSEED)
      ! Initialize tabulating the sill
      variosills(iv) = 0
      ! Note THIS VERSION DOES NOT INCORPORATE WEIGHTS!!!
      ! Monte Carlo Sample for PR Sill Determination
      ndraw = 0
      do i=1,NPRSAMPLE
        ! Draw a random sample
        j = max(min(nint(ndata*grnd())+1, ndata), 1)
        if ((vals(j,varioheads(iv)).gt.tmin) .and. (vals(j,varioheads(iv)).lt.tmax)) then
          val1 = vals(j,varioheads(iv))
        else
          cycle
        end if
        j = max(min(nint(ndata*grnd())+1, ndata), 1)
        if ((vals(j,varioheads(iv)).gt.tmin) .and. (vals(j,varioheads(iv)).lt.tmax)) then
          val2 = vals(j,varioheads(iv))
        else
          cycle
        end if
        if ((val1+val2).ge.EPSLON) then
          ndraw = ndraw + 1
          variosills(iv) = variosills(iv) + (val1-val2)**2 / ((val1+val2)*0.5)**2
        end if
      end do
      ! Report the sill
      if (ndraw.gt.0) then
        variosills(iv) = 0.5d0*variosills(iv)/ndraw
      else
        ! Error - no values drawn, report 0
        variosills(iv) = 0
      endif
    end if
! Sill of the continuous indicator variogram
    if (abs(variotypes(iv)).eq.9) then
      ! Determine the proportion of the indicator
      prop = 0.d0
      do i=1,ndata
        if ((vals(i,variotails(iv)).gt.tmin) .and. (vals(i,variotails(iv)).lt.tmax)) then
            if( vals(i,variotails(iv)) .lt. variocutcat(variotails(iv),2) ) prop = prop + 1.d0
          nvalues(iv,1) = nvalues(iv,1) + 1
        end if
      end do
      if (nvalues(iv,1).gt.0) then
        prop = prop / dble(nvalues(iv,1))
        variosills(iv) = prop*(1.d0-prop)
      else
        variosills(iv) = 0
      end if
    end if
! Sill of the categorical variogram
    if (abs(variotypes(iv)).eq.10) then
      ! Determine the proportion of the category
      prop = 0.d0
      do i=1,ndata
        if ((vals(i,variotails(iv)).gt.tmin) .and. (vals(i,variotails(iv)).lt.tmax)) then
            if( vals(i,variotails(iv)) == variocutcat(variotails(iv),2) ) prop = prop + 1.d0
          nvalues(iv,1) = nvalues(iv,1) + 1
        end if
      end do
      if (nvalues(iv,1).gt.0) then
        prop = prop / dble(nvalues(iv,1))
        variosills(iv) = prop*(1.d0-prop)
      else
        variosills(iv) = 0
      end if
    end if
  end do
  return
  end subroutine sillcalc

  subroutine varvisual(xyzorigin,azm,dip,tilt,azmtol,diptol,nlags,lagdist,lagtol, &
                bandhorz,bandvert,xyzvis,vtkfl)
!-----------------------------------------------------------------------
!  
!        varvisual - Experimental variogram tolerance plotting
!        *****************************************************
!
! This useful subroutine will generate a VTK file using the tolerance
! parameters specified for variogram calculation. This file can be viewed
! in the free ParaView software, or other free 3D visualizers. This is
! primarily aimed at directionaly, 3D variograms and may break if used
! for 3D omni or other strange cases. 
!
!  Parameters:
!    - 3D coordinate for drawing the origin of the variogram
!       xyzorigin(3) - x,y,z coordinates for variogram origin for drawing
!    - Experimental variogram parameters
!       azm - azimuth angle measured positive clockwise from the +y axis (N-S)
!       dip - dip angle measured negative down from the xy plane along the azimuth
!       tilt - tilt or plunge angle
!       azmtol - azimuth angle tolerance
!       diptol - dip angle tolerance
!       nlags - number of experimental lags to compute
!       lagdist - experimental lag distance (h)
!       lagtol - experimental lag distance tolerance
!       bandhorz - horizontal bandwidth (perpendicular to azm-Z plane)
!       bandvert - vertical bandwidth (perpendicular to horizontal plane)       
!    - Returns:
!       xyzvis(13,3)
!    - Writes:
!       vtfkfl - triangulated VTK file with variogram plotted
!
! (C) Jared Deutsch 2014
!-----------------------------------------------------------------------
  use varsubs ! Variogram subroutines including rotations
  implicit none
  integer, parameter :: lout=3
  real, parameter :: EPSLON=1.0e-6
  real(kind=8), parameter :: BIGDBLE=1d21,SMALLDBLE=1d-6,PI=4*atan(1.0d0)
  real(kind=8), parameter :: DEG2RAD=PI/180d0
! Subroutine parameters
  integer, intent(in) :: nlags
  real(kind=8), dimension(3), intent(in)  :: xyzorigin
  real(kind=8), intent(in) :: azm,dip,tilt,diptol,lagdist,lagtol, &
                      azmtol,bandhorz,bandvert
  character(len=512), intent(in) :: vtkfl
! Returns
  real(kind=8), dimension(13,3), intent(out) :: xyzvis
! Internal variables
  logical :: omni
  real(kind=8), dimension(4,3) :: xyzvissave
  real(kind=8), dimension(3,3) :: rotmat,forwardrotmat,reverserotmat
  real(kind=8) :: angle,azmtolrad,diptolrad,diptoldist,azmtoldist,farthest, &
                  ymin,ymax,zmin,zmax,xdist,h
  integer :: i,j,k,iv,outputidx
! Visualization variables
  real(kind=8), dimension(3) :: xyzmid
  real(kind=8) :: maxpltdist
  integer :: ios,colour,nsave,nstart
! Polygon indexing
  character(len=12), parameter, dimension(12) :: polypoints = & 
    (/'3 0 1 2     ', &
      '3 0 2 4     ', &
      '3 0 1 3     ', &
      '3 0 3 4     ', &
      '4 1 2 6 5   ', &
      '4 2 4 8 6   ', &
      '4 1 3 7 5   ', &
      '4 3 4 8 7   ', &
      '4 5 6 10 9  ', &
      '4 6 8 12 10 ', &
      '4 5 7 11 9  ', &
      '4 7 8 12 11 '/)

! Correct the angle tolerances to lie between 0 and 90 degrees
  omni = .false.
  if (azmtol .ge. 90d0) then
    azmtolrad = 90d0*DEG2RAD
    omni = .true.
  elseif (azmtol .lt. 0d0) then
    azmtolrad =  0d0
  else
    azmtolrad = azmtol*DEG2RAD
  end if
  if (diptol .ge. 90d0) then
    diptolrad = 90d0*DEG2RAD
    omni = .true.
  elseif (diptol .lt. 0d0) then
    diptolrad =  0d0
  else
    diptolrad = diptol*DEG2RAD
  end if
! Precalculate the point where the angle tolerances no longer come in to play
  if (azmtol .gt. 90d0-SMALLDBLE) then
    azmtoldist = 0d0
  elseif (azmtol .lt. SMALLDBLE) then
    azmtoldist = BIGDBLE
  else
    azmtoldist = bandhorz/tan(azmtolrad)
  end if
  if (diptol .gt. 90d0-SMALLDBLE) then
    diptoldist = 0d0
  elseif (diptol .lt. SMALLDBLE) then
    diptoldist = BIGDBLE
  else
    diptoldist = bandvert/tan(diptolrad)
  end if
! Set up the rotation matrix using the GSLIB angles
  call setrotmat(azm,dip,tilt,forwardrotmat,reverserotmat)
! Maximum plotting distance
  maxpltdist = nlags*lagdist+lagtol
! Initialize
  xyzvis(:,:) = 0
! Sets of points to plot (before translation):
!   xyzvis(1,:) - origin at 0,0,0
!   xyzvis(2:5,:) - 4 points corresponding to azmtoldist,?,? if azmtoldist < maxpltdist
!   xyzvis(6:9,:) - 4 points corresponding to diptoldist,?,? if diptoldist < maxpltdist
!   xyzvis(10:13,:) - 4 points corresponding to maxpltdist,?,?
  nsave = 1
  nstart = nsave+1
! Azimuth tolerance
  if (azmtoldist .lt. maxpltdist) then
      ! Plot the azimuth points
      nsave = nsave + 4
      ! x'
      xyzvis(nstart:nsave,1) = azmtoldist
      ! y'
      ymax = bandhorz
      ymin = -1d0*bandhorz
      ! z'
      if (azmtoldist .lt. diptoldist) then
        zmax = azmtoldist*tan(diptolrad)
        zmin = -1d0*zmax
      else
        zmax = bandvert
        zmin = -1d0*bandvert
      end if
      xyzvis(nstart,2) = ymax
      xyzvis(nstart+1,2) = ymin
      xyzvis(nstart+2,2) = ymax
      xyzvis(nstart+3,2) = ymin
      xyzvis(nstart,3) = zmax
      xyzvis(nstart+1,3) = zmax
      xyzvis(nstart+2,3) = zmin
      xyzvis(nstart+3,3) = zmin
  end if
  nstart = nsave+1
! Dip tolerance
  if (diptoldist .lt. maxpltdist) then
      ! Plot the dip points
      nsave = nsave + 4
      ! x'
      xyzvis(nstart:nsave,1) = diptoldist
      ! y'
      if (diptoldist .lt. azmtoldist) then
        ymax = diptoldist*tan(azmtolrad)
        ymin = -1d0*ymax
      else
        ymax = bandhorz
        ymin = -1d0*bandhorz
      end if
      ! z'
      zmax = bandvert
      zmin = -1d0*bandvert
      xyzvis(nstart,2) = ymax
      xyzvis(nstart+1,2) = ymin
      xyzvis(nstart+2,2) = ymax
      xyzvis(nstart+3,2) = ymin
      xyzvis(nstart,3) = zmax
      xyzvis(nstart+1,3) = zmax
      xyzvis(nstart+2,3) = zmin
      xyzvis(nstart+3,3) = zmin
  end if
! Swap dip and azimuth tolerances if dip < azm
  if ((diptoldist .lt. maxpltdist) .and. &
      (azmtoldist .lt. maxpltdist) .and. &
      (diptoldist .lt. azmtoldist)) then
      xyzvissave = xyzvis(2:5,:)
      xyzvis(2:5,:) = xyzvis(6:9,:)
      xyzvis(6:9,:) = xyzvissave
  end if
  nstart = nsave+1
! Max distance
  nsave = nsave + 4
  ! x'
  xyzvis(nstart:nsave,1) = maxpltdist
  ! y'
  if (maxpltdist .lt. azmtoldist) then
    ymax = maxpltdist*tan(azmtolrad)
    ymin = -1d0*ymax
  else
    ymax = bandhorz
    ymin = -1d0*bandhorz
  end if
  ! z'
  if (maxpltdist .lt. diptoldist) then
    zmax = maxpltdist*tan(diptolrad)
    zmin = -1d0*zmax
  else
    zmax = bandvert
    zmin = -1d0*bandvert
  end if
  xyzvis(nstart,2) = ymax
  xyzvis(nstart+1,2) = ymin
  xyzvis(nstart+2,2) = ymax
  xyzvis(nstart+3,2) = ymin
  xyzvis(nstart,3) = zmax
  xyzvis(nstart+1,3) = zmax
  xyzvis(nstart+2,3) = zmin
  xyzvis(nstart+3,3) = zmin
! Rotate back to original coordinates and translate
  xyzvis = matmul(xyzvis,reverserotmat)
  do i=1,nsave
    do j=1,3
      xyzvis(i,j) = xyzvis(i,j) + xyzorigin(j)
    end do
  end do
! Save out the vtk file
  open(lout,file=vtkfl,status='UNKNOWN')
  write(lout,'(a)') '# vtk DataFile Version 3.0'
  write(lout,'(a)') 'Variogram Spec Visualization'
  write(lout,'(a)') 'ASCII'
  write(lout,'(a)') 'DATASET POLYDATA'
  write(lout,'(a,i3,a)') 'POINTS',nsave,' float'
  do i=1,nsave
    write(lout,'(3f15.5)') (xyzvis(i,j),j=1,3)
  end do
  write(lout,*)
  write(lout,'(a,i3)') 'VERTICES 1',nsave+1
  write(lout,'(14i3)') nsave,(j,j=0,nsave-1)
  write(lout,*)
  if (nsave .eq. 5) then
    j = 16
  else
    j = 16+(nsave-5)*5
  end if
  write(lout,'(a,2i3)') 'POLYGONS',nsave-1,j
  do i=1,nsave-1
    write(lout,'(a)') polypoints(i)
  end do
  close(lout)
! Return with the variogram points
  return
  end subroutine varvisual

  subroutine varpairs(xyz,xyzorigin,ndata,azm,dip,tilt,azmtol,diptol, &
                nlags,lagdist,lagtol,bandhorz,bandvert,paired)
!-----------------------------------------------------------------------
!  
!             varpairs - Experimental variogram pairs
!             ***************************************
!
! This program is made to show what points would be paired given a
! set of variogram direction parameters and the origin 0,0,0
!
!  Parameters:
!    - 3D data
!       xyzorigin(3) - x,y,z coordinates of the origin
!       xyz(ndata,3) - x,y,z coordinates to pair with
!       ndata - number of data
!    - Experimental variogram parameters
!       azm - azimuth angle measured positive clockwise from the +y axis (N-S)
!       dip - dip angle measured negative down from the xy plane along the azimuth
!       tilt - tilt (or plunge) angle measured positive counterclockwise about Y'
!       azmtol - azimuth angle tolerance
!       diptol - dip angle tolerance
!       nlags - number of experimental lags to compute
!       lagdist - experimental lag distance (h)
!       lagtol - experimental lag distance tolerance
!       bandhorz - horizontal bandwidth (perpendicular to azm-Z plane)
!       bandvert - vertical bandwidth (perpendicular to horizontal plane)       
!    - Returns:
!       paired(ndata) - indicator of whether or not xyz(ndata,:) is paired
!
! (C) Jared Deutsch 2014
!-----------------------------------------------------------------------
  use sortem ! Quicksort module
  use varsubs ! Variogram subroutine module
  implicit none
  real, parameter :: EPSLON=1.0e-6
  real(kind=8), parameter :: BIGDBLE=1d21,SMALLDBLE=1d-6,PI=4*atan(1.0d0)
  real(kind=8), parameter :: DEG2RAD=PI/180d0
  integer, parameter :: lout=3
! Subroutine parameters
  integer, intent(in) :: ndata,nlags
  real(kind=8), dimension(ndata,3), intent(in)  :: xyz
  real(kind=8), intent(in) :: azm,dip,tilt,diptol,lagdist,lagtol, &
                      azmtol,bandhorz,bandvert
  real(kind=8), dimension(3), intent(in) :: xyzorigin
! Returns
  integer, dimension(ndata), intent(out) :: paired
! Internal variables
  logical :: omni
  real(kind=8), dimension(ndata,3) :: xyzrot
  real(kind=8), dimension(ndata) :: dblerotindex
  integer, dimension(ndata) :: introtindex
  real(kind=8), dimension(3,3) :: rotmat,forwardrotmat,reverserotmat
  real(kind=8) :: angle,azmtolrad,diptolrad,diptoldist,azmtoldist,farthest, &
                  ymin,ymax,zmin,zmax,xdist,h
  integer :: i,j,k,m,iv,outputidx,startidx
! Initialize return variables
  paired = 0
! Initialize index tracking for values
  do i=1,ndata
    dblerotindex(i) = i
  end do
! Correct the angle tolerances to lie between 0 and 90 degrees
  omni = .false.
  if (azmtol .ge. 90d0) then
    azmtolrad = 90d0*DEG2RAD
    omni = .true.
  elseif (azmtol .lt. 0d0) then
    azmtolrad =  0d0
  else
    azmtolrad = azmtol*DEG2RAD
  end if
  if (diptol .ge. 90d0) then
    diptolrad = 90d0*DEG2RAD
    omni = .true.
  elseif (diptol .lt. 0d0) then
    diptolrad =  0d0
  else
    diptolrad = diptol*DEG2RAD
  end if
! Precalculate the point where the angle tolerances no longer come in to play
  if (azmtol .gt. 90d0-SMALLDBLE) then
    azmtoldist = 0d0
  elseif (azmtol .lt. SMALLDBLE) then
    azmtoldist = BIGDBLE
  else
    azmtoldist = bandhorz/tan(azmtolrad)
  end if
  if (diptol .gt. 90d0-SMALLDBLE) then
    diptoldist = 0d0
  elseif (diptol .lt. SMALLDBLE) then
    diptoldist = BIGDBLE
  else
    diptoldist = bandvert/tan(diptolrad)
  end if
! Set up the rotation matrix using the GSLIB angles
  call setrotmat(azm,dip,tilt,forwardrotmat,reverserotmat)
! Rotate to align the new x' axis with the variogram direction
  xyzrot = matmul(xyz,forwardrotmat)
! Sort along the direction vector
  call dblemodsortem(xyzrot(:,1),ndata,3,xyzrot(:,2),xyzrot(:,3),dblerotindex)
! Track the sorted indices for computing values
  do i=1,ndata
    introtindex(i) = int(dblerotindex(i))
  end do
! Main variogram loop over all points indexed by i
  do i=1,1   ! This point, i, is the "tail" variable in GSLIB notation
    xyzrot(i,1) = xyzorigin(1)
    xyzrot(i,2) = xyzorigin(2)
    xyzrot(i,3) = xyzorigin(3)
    ! Farthest point to consider from the current location
    farthest = xyzrot(i,1) + nlags*lagdist + lagtol
    ! Pre-calculated bandwidth min/max y's and z's
    ymax = xyzrot(i,2) + bandhorz
    ymin = xyzrot(i,2) - bandhorz
    zmax = xyzrot(i,3) + bandvert
    zmin = xyzrot(i,3) - bandvert
    ! Loop over all possible pairs forward of the current point
    if (omni) then
      startidx = 1
    else
      startidx = i+1
    end if
    do j=startidx,ndata   ! This point, j, is the "head" variable in GSLIB notation
      ! For this visualization module only - check that we are pointed the right way!
      if (xyzrot(j,1) .le. xyzrot(i,1) .and. .not. omni) cycle
      ! If we are too far away, skip checking other pairs
      if (xyzrot(j,1) .gt. farthest) exit
      ! Check horizontal bandwidths
      if (xyzrot(j,2) .gt. ymax) cycle
      if (xyzrot(j,2) .lt. ymin) cycle
      ! Check vertical bandwidths
      if (xyzrot(j,3) .gt. zmax) cycle
      if (xyzrot(j,3) .lt. zmin) cycle
      ! Check if the points are at the same point on the x' axis
      xdist = abs(xyzrot(j,1)-xyzrot(i,1))
      if (i.eq.j) then
        cycle
      else
        ! Check azimuth tolerance if it matters at this distance
        if (xdist .lt. azmtoldist) then
          ! Check absolute azimuth tolerance
          angle = abs(atan((xyzrot(j,2)-xyzrot(i,2))/(max(xdist,SMALLDBLE))))
          if (angle .gt. azmtolrad) cycle
        end if
        ! Check dip tolerance if it matters at this distance
        if (xdist .lt. diptoldist) then
          ! Check absolute dip tolerance
          angle = abs(atan((xyzrot(j,3)-xyzrot(i,3))/(max(xdist,SMALLDBLE))))
          if (angle .gt. diptolrad) cycle
        end if
      end if
      ! At this point all tests have been passed so place the pair in the
      ! appropriate bin - due to possible overlapping bins each bin
      ! should be checked to see if pairs fall in it
      ! Find the distance between the values
      h = distance(xyzrot(i,:),xyzrot(j,:))
      do k=0,nlags
        ! Are we inside this bin
        if ((h .ge. (k*lagdist-lagtol)) .and. &
            (h .le. (k*lagdist+lagtol))) then
          ! We are inside the bin, so calculate each variogram value
          paired(introtindex(j)) = k+1
        end if
      end do
    end do
  end do
  return
  end subroutine varpairs
