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

module varsubs
  public :: distance,expgam,setrotmat
contains

  real(kind=8) function distance(xyz_a,xyz_b)
!-----------------------------------------------------------------------
!  
!                  Distance given 2 x,y,z vectors
!
! (c) Jared Deutsch, 2014
!-----------------------------------------------------------------------
    real(kind=8), dimension(3), intent(in) :: xyz_a,xyz_b
    distance = sqrt((xyz_a(1)-xyz_b(1))**2 + &
                    (xyz_a(2)-xyz_b(2))**2 + &
                    (xyz_a(3)-xyz_b(3))**2)
  end function distance

  real(kind=8) function expgam(variotail,variohead,variotype, &
                               variotail2,variohead2,variocutcat1,variocutcat2)
!-----------------------------------------------------------------------
!  
! Intermediate experimental variogram value given head, tail and type
! values
!
! Notes:
!   - The cross variogram definition used is the full definition
!     where head and tail values are required at u and u+h lags!
!     This is the definition used in GSLIB. 
!
! (c) Jared Deutsch, 2014
!-----------------------------------------------------------------------
    implicit none
    real(kind=8), intent(in) :: variotail,variohead
    ! Optional cross variogram and indicator variogram parameters
    real(kind=8), intent(in), optional :: variotail2,variohead2,variocutcat2,variocutcat1
    integer, intent(in) :: variotype
    integer :: tailind,headind,tail2ind,head2ind
    
    select case (variotype)
      case(1) ! Traditional semivariogram
        expgam = (variotail-variohead)*(variotail-variohead)
      case(2) ! Traditional cross semivariogram
        expgam = (variotail-variohead)*(variotail2-variohead2)
      case(3) ! Covariance (non-ergodic covariance)
        expgam = variotail*variohead
      case(4) ! Correlogram
        expgam = variotail*variohead
      case(5) ! General relative semivariogram
        expgam = (variotail-variohead)*(variotail-variohead)
      case(6) ! Pairwise relative semivariogram
        expgam = (variotail-variohead)*(variotail-variohead)/((variotail+variohead)/2)**2
      case(7) ! Semivariogram of logarithms (ln)
        expgam = (log(variotail)-log(variohead))**2
      case(8) ! Semimadogram
        expgam = abs(variotail-variohead)
      case(9) ! Indicator semivariogram - continuous
        if (variotail .lt. variocutcat1) then
          tailind = 0
        else
          tailind = 1
        end if
        if (variohead .lt. variocutcat1) then
          headind = 0
        else
          headind = 1
        end if
        expgam = (tailind-headind)*(tailind-headind)
      case(10) ! Indicator semivariogram - categorical
        if (int(variotail) .eq. int(variocutcat1)) then
          tailind = 1
        else
          tailind = 0
        end if
        if (int(variohead) .eq. int(variocutcat1)) then
          headind = 1
        else
          headind = 0
        end if
        expgam = (tailind-headind)*(tailind-headind)
      case(11) ! Traditional cross indicator semivariogram - continuous
        if (variotail .lt. variocutcat1) then
          tailind = 0
        else
          tailind = 1
        end if
        if (variohead .lt. variocutcat1) then
          headind = 0
        else
          headind = 1
        end if
        if (variotail2 .lt. variocutcat2) then
          tail2ind = 0
        else
          tail2ind = 1
        end if
        if (variohead2 .lt. variocutcat2) then
          head2ind = 0
        else
          head2ind = 1
        end if
        expgam = (tailind-headind)*(tail2ind-head2ind)
      case(12) ! Traditional cross indicator semivariogram - categorical
        if (int(variotail) .eq. int(variocutcat1)) then
          tailind = 1
        else
          tailind = 0
        end if
        if (int(variohead) .eq. int(variocutcat1)) then
          headind = 1
        else
          headind = 0
        end if
        if (int(variotail2) .eq. int(variocutcat2)) then
          tail2ind = 1
        else
          tail2ind = 0
        end if
        if (int(variohead2) .eq. int(variocutcat2)) then
          head2ind = 1
        else
          head2ind = 0
        end if
        expgam = (tailind-headind)*(tail2ind-head2ind)
      case default
        write(*,*) 'ERROR: Invalid variogram type',variotype
        stop
    end select

  end function expgam

  subroutine setrotmat(azm,dip,tilt,forwardrotmat,reverserotmat)
!-----------------------------------------------------------------------
!  
!        Sets up the forward rotation matrix and reverse matrix
!
! (c) Jared Deutsch, 2014
!-----------------------------------------------------------------------
  implicit none
  real(kind=8), parameter :: PI = 4*atan(1.0d0)
  real(kind=8), intent(in) :: azm,dip,tilt
  real(kind=8), dimension(3,3) :: azmrotmat,diprotmat,tiltrotmat
  real(kind=8), dimension(3,3), intent(out) :: forwardrotmat,reverserotmat
  real(kind=8) :: angle

! Get the angle to rotate about the Z axis in radians
  angle = (90d0-azm)*PI/180d0
! Rotation matrix for azimuth correction
  azmrotmat = 0d0
  azmrotmat(3,3) = 1d0
  azmrotmat(1,1) = cos(angle)
  azmrotmat(1,2) = -1d0*sin(angle)
  azmrotmat(2,1) = sin(angle)
  azmrotmat(2,2) = cos(angle)

! Get the angle to rotate about the new X' axis
  angle = -1d0*(dip)*PI/180d0
! Rotation matrix for dip correction
  diprotmat = 0d0
  diprotmat(2,2) = 1d0
  diprotmat(1,1) = cos(angle)
  diprotmat(1,3) = sin(angle)
  diprotmat(3,1) = -1d0*sin(angle)
  diprotmat(3,3) = cos(angle)

! Get the angle to rotate about the new Y' axis
  angle = -1d0*(tilt)*PI/180d0
! Rotation matrix for tilt correction
  tiltrotmat = 0d0
  tiltrotmat(1,1) = 1d0
  tiltrotmat(2,2) = cos(angle)
  tiltrotmat(2,3) = sin(angle)
  tiltrotmat(3,2) = -1d0*sin(angle)
  tiltrotmat(3,3) = cos(angle)

! Complete forward rotation matrix is the product of these 2 matrices
  forwardrotmat = matmul(matmul(azmrotmat,diprotmat),tiltrotmat)
! Reverse rotation matrix is the transpose of the forward matrix
! as these matrices are orthogonal
  reverserotmat = transpose(forwardrotmat)
  end subroutine setrotmat

end module varsubs
