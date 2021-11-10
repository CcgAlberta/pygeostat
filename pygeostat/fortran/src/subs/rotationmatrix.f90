! These subroutines are copyright John Manchuk (2013) and used
! with permission. 

!Module with types and functions for using rotation matrices
!in gslib type programs
!COPYRIGHT JOHN MANCHUK JANUARY 2013
module rotationmatrix
    implicit none
    
    public :: set_rmat, &
              dsqrd, &
              fix_angles
    private
    
    real*8, parameter :: PI=3.141592653589793238462643383279502884197D0
    real*8, parameter :: DEG2RAD = PI / 180.D0

contains

    !Populate a rotation matrix based on angles (uncorrected!), and
    !ranges that describe a variogram or search ellipse.
    function set_rmat ( angles, ranges ) result ( rm )
        real*8, intent (in) :: angles( 3 )
        real*8, intent (in) :: ranges( 3 )
        real*8 :: sina, sinb, sint, cosa, cosb, cost
        real*8, dimension ( 3 ) :: a
        real*8 :: r1, r2
        real*8 :: rm(3,3)
        
        ! Fix the angles and get the required sines and cosines:
        a = angles
        call fix_angles ( a )
        sina = dsin(a(1)) ; sinb = dsin(a(2)) ; sint = dsin(a(3))
        cosa = dcos(a(1)) ; cosb = dcos(a(2)) ; cost = dcos(a(3))

        ! Construct the rotation matrix:
        r1 = ranges(1) / ranges(2)
        r2 = ranges(1) / ranges(3)
        rm(1,:) = (/ cosb * cosa , cosb * sina , -sinb /)
        rm(2,:) = (/ r1 * (-cost * sina + sint * sinb * cosa) , &
                     r1 * ( cost * cosa + sint * sinb * sina) , &
                     r1 * ( sint * cosb) /)
        rm(3,:) = (/ r2 * ( sint * sina + cost * sinb * cosa) , &
                     r2 * (-sint * cosa + cost * sinb * sina) , &
                     r2 * ( cost * cosb) /)
        Return
    end function set_rmat
    
    !Convert three angles of anisotropy to those that are usable.  Note that
    !the vector of angles will be updated and returned
    subroutine fix_angles ( angles )
        real*8, intent (inout) :: angles( 3 )
        !1st - angle between the major axis of anisotropy and the
        !      E-W axis. Note: Counter clockwise is positive.
        !2nd - angle between major axis and the horizontal plane.
        !      (The dip of the ellipsoid measured positive down)
        !3rd - Angle of rotation of minor axis about the major axis
        !      of the ellipsoid.
        if(angles(1) >= 0.0D+00 .and. angles(1) < 270.0D+00) then
            angles(1) = ( 90.0D+00 - angles(1)) * DEG2RAD
        else
            angles(1) = (450.0D+00 - angles(1)) * DEG2RAD
        endif
        angles(2) = -1.0D+00 * angles(2) * DEG2RAD
        angles(3) =  1.0D+00 * angles(3) * DEG2RAD
        Return
    end subroutine fix_angles
    
    !Get the anisotropic distance between two points
    real*8 function dsqrd ( a, b, rm )
        real*8, intent (in) :: a(3), b(3)
        real*8, optional :: rm(3,3)
        real*8 :: ba(3)
        ba = b - a
        if(present(rm))then
            dsqrd = (rm(1,1)*ba(1) + rm(1,2)*ba(2) + rm(1,3)*ba(3))**2 + &
                    (rm(2,1)*ba(1) + rm(2,2)*ba(2) + rm(2,3)*ba(3))**2 + &
                    (rm(3,1)*ba(1) + rm(3,2)*ba(2) + rm(3,3)*ba(3))**2
        else
            dsqrd = ba(1)**2 + ba(2)**2 + ba(3)**2
        endif
    end function dsqrd

end module rotationmatrix