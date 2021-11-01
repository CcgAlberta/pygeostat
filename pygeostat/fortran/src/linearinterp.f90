! Linear and bilinear interpolation routines, and numerical gradient on a grid
!
! Compile as a python DLL with the command:
!     call f2py -c -m --fcompiler=gnu95 --opt='-O2' linearinterp linearinterp.f90
!
! Note that in {PYTHONINSTALLDIR}/Lib/site-packages/numpy/distutils/fcompiler/gnu.py
!   line 337 (approximately)
!      raise NotImplementedError("Only MS compiler supported with gfortran on win64")
!   must be changed to
!      pass #raise NotImplementedError("Only MS compiler supported with gfortran on win64")
!
!-----------------------------------------------------------------------------
! Copyright (c) 2015, Jared L. Deutsch
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

module linearinterp
    implicit none
    
    ! Linear interpolation functions
    PUBLIC ::  lerp_d, &      ! double precision 1D linear interpolation
               lerp_d_trim, & ! double precision 1D trimmed linear interpolation
               lerp_s, &      ! single precision 1D linear interpolation
               bilerp_d, &    ! double precision 2D bilinear interpolation
               bilerp_d_trim, & ! double precision trimmed bilinear interpolation
               bilerp_s,&     ! single precision 2D bilinear interpolation
               gridbilerp_d,& ! double precision 2D bilinear interpolation on a grid
               gridbilerp_s,& ! single precision 2D bilinear interpolation on a grid
               checkgradient_d, & ! double precision 2D gradient checker for ability to compute on a grid
               gradient_d,&   ! double precision 2D gradient using a grid
               gridgradient_d, & ! double precision 2D gradient on a grid
               discretizesegment, & ! discretizes a line segment
               extractpoints_d, & ! 1D grid extraction from a line segment
               unitperpvector, & ! Computes the unit perpendicular vector to a line given by 2 points
               discpolyline, & ! Discretize a polyline
               unfoldcoords ! Unfolds a set of 2D coordinates onto a 1D line with perpendicular components

contains

    ! 1D linear interpolation (double precision)
    real*8 function lerp_d(x1,y1,x2,y2,x) result(ylerp)
        real*8, intent(in) :: x1,y1,x2,y2,x ! (x1,y1), (x2,y2), x
        ylerp = y1 + (y2-y1)/(x2-x1)*(x-x1)
    end function lerp_d

    ! 1D linear interpolation (double precision, trimmed)
    real*8 function lerp_d_trim(x1,y1,x2,y2,x,tmin) result(ylerp)
        real*8, intent(in) :: x1,y1,x2,y2,x,tmin ! (x1,y1), (x2,y2), x, lower trimming limit
        if (y1 .gt. tmin .and. y2 .gt. tmin) then
            ylerp = lerp_d(x1,y1,x2,y2,x)
        elseif (y1 .gt. tmin .and. y2 .le. tmin) then
            ylerp = y1
        elseif (y1 .le. tmin .and. y2 .gt. tmin) then
            ylerp = y2
        else
            ylerp = tmin
        end if
    end function lerp_d_trim

    ! 1D linear interpolation (single precision)
    real function lerp_s(x1,y1,x2,y2,x) result(ylerp)
        real, intent(in) :: x1,y1,x2,y2,x ! (x1,y1), (x2,y2), x
        ylerp = y1 + (y2-y1)/(x2-x1)*(x-x1)
    end function lerp_s

    ! 2D bilinear interpolation (double precision)
    ! Note: q12 is the interpolating value at x1,y2 and so on...
    real*8 function bilerp_d(x1,y1,x2,y2,q11,q12,q21,q22,x,y) result(qlerp)
        real*8, intent(in) :: x1,y1,x2,y2,q11,q12,q21,q22,x,y
        real*8 :: r1,r2
        r1 = ((x2-x)/(x2-x1))*q11 + ((x-x1)/(x2-x1))*q21
        r2 = ((x2-x)/(x2-x1))*q12 + ((x-x1)/(x2-x1))*q22
        qlerp = ((y2-y)/(y2-y1))*r1 + ((y-y1)/(y2-y1))*r2
    end function bilerp_d

    ! 2D bilinear interpolation (double precision, trimmed)
    ! Note: q12 is the interpolating value at x1,y2 and so on...
    real*8 function bilerp_d_trim(x1,y1,x2,y2,q11,q12,q21,q22,x,y,tmin) result(qlerp)
        real*8, intent(in) :: x1,y1,x2,y2,q11,q12,q21,q22,x,y,tmin ! tmin = lower trimming limit
        real*8 :: r1,r2
        if (q11 .gt. tmin .and. q21 .gt. tmin) then
            r1 = ((x2-x)/(x2-x1))*q11 + ((x-x1)/(x2-x1))*q21
        elseif (q11 .gt. tmin .and. q21 .le. tmin) then
            r1 = q11
        elseif (q11 .le. tmin .and. q21 .gt. tmin) then
            r1 = q21
        else
            r1 = tmin
        end if
        if (q12 .gt. tmin .and. q22 .gt. tmin) then
            r2 = ((x2-x)/(x2-x1))*q12 + ((x-x1)/(x2-x1))*q22
        elseif (q12 .gt. tmin .and. q22 .le. tmin) then
            r2 = q12
        elseif (q12 .le. tmin .and. q22 .gt. tmin) then
            r2 = q22
        else
            r2 = tmin
        end if
        if (r1 .gt. tmin .and. r2 .gt. tmin) then
            qlerp = ((y2-y)/(y2-y1))*r1 + ((y-y1)/(y2-y1))*r2
        elseif (r1 .gt. tmin .and. r1 .le. tmin) then
            qlerp = r1
        elseif (q12 .le. tmin .and. r2 .gt. tmin) then
            qlerp = r2
        else
            qlerp = tmin
        end if
    end function bilerp_d_trim

    ! 2D bilinear interpolation (single precision)
    ! Note: q12 is the interpolating value at x1,y2 and so on...
    real function bilerp_s(x1,y1,x2,y2,q11,q12,q21,q22,x,y) result(qlerp)
        real, intent(in) :: x1,y1,x2,y2,q11,q12,q21,q22,x,y
        real :: r1,r2
        r1 = ((x2-x)/(x2-x1))*q11 + ((x-x1)/(x2-x1))*q21
        r2 = ((x2-x)/(x2-x1))*q12 + ((x-x1)/(x2-x1))*q22
        qlerp = ((y2-y)/(y2-y1))*r1 + ((y-y1)/(y2-y1))*r2
    end function bilerp_s

    ! 2D index in to grid
    integer function grididx(ix,iy,nx)
        integer, intent(in) :: ix,iy,nx
        grididx = ix + (iy-1)*nx
    end function grididx

    ! 2D bilinear interpolation on a regular 2D grid (double precision)
    ! Note: Requires that x and y already be standardized (ie: x-xmin,y-ymin)
    subroutine gridbilerp_d(nx,xsiz,ny,ysiz,qvalues,x,y,qlerp)
        integer, intent(in) :: nx,ny
        real*8, intent(in) :: xsiz,ysiz,x,y,qvalues(nx*ny)
        real*8, intent(out) :: qlerp
        integer :: ix,iy
        real*8 :: x1,y1,x2,y2
        ! Get the lower left hand ix,iy coordinate
        ix = floor((x-0.5d0*xsiz)/xsiz) + 1
        iy = floor((y-0.5d0*ysiz)/ysiz) + 1
        ! Clip ix,iy to fall within 1,nx and 1,ny
        ix = max(min(ix,nx),1)
        iy = max(min(iy,ny),1)
        ! Is this an edge case (falls on the perimeter) or outside the grid? If so just return that grid value
        ! Check the x edges and corners
        if (x <= (0.5d0*xsiz)) then ! left hand edge
            if (y <= (0.5d0*ysiz)) then ! bottom left hand corner
                !write(*,*) 'bottom left hand corner'
                qlerp = qvalues(grididx(ix,iy,nx))
                return
            elseif (y >= (ny*ysiz-0.5d0*ysiz)) then ! top left hand corner
                !write(*,*) 'top left hand corner'
                qlerp = qvalues(grididx(ix,iy,nx))
                return
            else
                !write(*,*) 'left edge somewhere in the middle'
                y1 = iy*ysiz - 0.5d0*ysiz
                y2 = y1 + ysiz
                qlerp = lerp_d(y1,qvalues(grididx(ix,iy,nx)), &
                               y2,qvalues(grididx(ix,iy+1,nx)), y)
                return
            end if
        elseif (x >= (nx*xsiz-(0.5d0*xsiz))) then ! right hand edge
            if (y <= (0.5d0*ysiz)) then ! bottom right hand corner
                !write(*,*) 'bottom right hand corner'
                qlerp = qvalues(grididx(ix,iy,nx))
                return
            elseif (y >= (ny*ysiz-0.5d0*ysiz)) then ! top right hand corner
                !write(*,*) 'top right hand corner'
                qlerp = qvalues(grididx(ix,iy,nx))
                return
            else
                !write(*,*) 'right edge somewhere in the middle'
                y1 = iy*ysiz - 0.5d0*ysiz
                y2 = y1 + ysiz
                qlerp = lerp_d(y1,qvalues(grididx(ix,iy,nx)), &
                               y2,qvalues(grididx(ix,iy+1,nx)), y)
                return
            end if
        end if
        ! Check the top and bottom edges, corners have already been checked
        if (y <= (0.5d0*ysiz)) then ! bottom right hand corner
            !write(*,*) 'bottom edge somewhere in the middle'
            x1 = ix*xsiz - 0.5d0*xsiz
            x2 = x1 + xsiz
            qlerp = lerp_d(x1,qvalues(grididx(ix,iy,nx)), &
                           x2,qvalues(grididx(ix+1,iy,nx)), x)
            return
        elseif (y >= (ny*ysiz-0.5d0*ysiz)) then 
            !write(*,*) 'top edge somewhere in the middle'
            x1 = ix*xsiz - 0.5d0*xsiz
            x2 = x1 + xsiz
            qlerp = lerp_d(x1,qvalues(grididx(ix,iy,nx)), &
                           x2,qvalues(grididx(ix+1,iy,nx)), x)
            return
        end if
        ! This is not an edge case, perform bilinear interpolation
        !write(*,*) 'somewhere in the middle'
        x1 = ix*xsiz - 0.5d0*xsiz
        x2 = x1 + xsiz
        y1 = iy*ysiz - 0.5d0*ysiz
        y2 = y1 + ysiz
        qlerp = bilerp_d(x1,y1,x2,y2, &
                         qvalues(grididx(ix,iy,nx)), &
                         qvalues(grididx(ix,iy+1,nx)), &
                         qvalues(grididx(ix+1,iy,nx)), &
                         qvalues(grididx(ix+1,iy+1,nx)), x, y)
        return
    end subroutine gridbilerp_d

    ! 2D trimmed bilinear interpolation on a regular 2D grid (double precision)
    ! Note: Requires that x and y already be standardized (ie: x-xmin,y-ymin)
    subroutine gridbilerp_d_trim(nx,xsiz,ny,ysiz,qvalues,x,y,qlerp,tmin)
        integer, intent(in) :: nx,ny
        real*8, intent(in) :: xsiz,ysiz,x,y,qvalues(nx*ny),tmin
        real*8, intent(out) :: qlerp
        integer :: ix,iy
        real*8 :: x1,y1,x2,y2
        ! Get the lower left hand ix,iy coordinate
        ix = floor((x-0.5d0*xsiz)/xsiz) + 1
        iy = floor((y-0.5d0*ysiz)/ysiz) + 1
        ! Clip ix,iy to fall within 1,nx and 1,ny
        ix = max(min(ix,nx),1)
        iy = max(min(iy,ny),1)
        ! Is this an edge case (falls on the perimeter) or outside the grid? If so just return that grid value
        ! Check the x edges and corners
        if (x <= (0.5d0*xsiz)) then ! left hand edge
            if (y <= (0.5d0*ysiz)) then ! bottom left hand corner
                !write(*,*) 'bottom left hand corner'
                qlerp = qvalues(grididx(ix,iy,nx))
                return
            elseif (y >= (ny*ysiz-0.5d0*ysiz)) then ! top left hand corner
                !write(*,*) 'top left hand corner'
                qlerp = qvalues(grididx(ix,iy,nx))
                return
            else
                !write(*,*) 'left edge somewhere in the middle'
                y1 = iy*ysiz - 0.5d0*ysiz
                y2 = y1 + ysiz
                qlerp = lerp_d_trim(y1,qvalues(grididx(ix,iy,nx)), &
                               y2,qvalues(grididx(ix,iy+1,nx)), y, tmin)
                return
            end if
        elseif (x >= (nx*xsiz-(0.5d0*xsiz))) then ! right hand edge
            if (y <= (0.5d0*ysiz)) then ! bottom right hand corner
                !write(*,*) 'bottom right hand corner'
                qlerp = qvalues(grididx(ix,iy,nx))
                return
            elseif (y >= (ny*ysiz-0.5d0*ysiz)) then ! top right hand corner
                !write(*,*) 'top right hand corner'
                qlerp = qvalues(grididx(ix,iy,nx))
                return
            else
                !write(*,*) 'right edge somewhere in the middle'
                y1 = iy*ysiz - 0.5d0*ysiz
                y2 = y1 + ysiz
                qlerp = lerp_d_trim(y1,qvalues(grididx(ix,iy,nx)), &
                               y2,qvalues(grididx(ix,iy+1,nx)), y, tmin)
                return
            end if
        end if
        ! Check the top and bottom edges, corners have already been checked
        if (y <= (0.5d0*ysiz)) then ! bottom right hand corner
            !write(*,*) 'bottom edge somewhere in the middle'
            x1 = ix*xsiz - 0.5d0*xsiz
            x2 = x1 + xsiz
            qlerp = lerp_d_trim(x1,qvalues(grididx(ix,iy,nx)), &
                           x2,qvalues(grididx(ix+1,iy,nx)), x, tmin)
            return
        elseif (y >= (ny*ysiz-0.5d0*ysiz)) then 
            !write(*,*) 'top edge somewhere in the middle'
            x1 = ix*xsiz - 0.5d0*xsiz
            x2 = x1 + xsiz
            qlerp = lerp_d_trim(x1,qvalues(grididx(ix,iy,nx)), &
                           x2,qvalues(grididx(ix+1,iy,nx)), x, tmin)
            return
        end if
        ! This is not an edge case, perform bilinear interpolation
        !write(*,*) 'somewhere in the middle'
        x1 = ix*xsiz - 0.5d0*xsiz
        x2 = x1 + xsiz
        y1 = iy*ysiz - 0.5d0*ysiz
        y2 = y1 + ysiz
        qlerp = bilerp_d_trim(x1,y1,x2,y2, &
                         qvalues(grididx(ix,iy,nx)), &
                         qvalues(grididx(ix,iy+1,nx)), &
                         qvalues(grididx(ix+1,iy,nx)), &
                         qvalues(grididx(ix+1,iy+1,nx)), x, y, tmin)
        return
    end subroutine gridbilerp_d_trim

    ! 2D bilinear interpolation on a regular 2D grid (single precision)
    ! Note: Requires that x and y already be standardized (ie: x-xmin,y-ymin)
    subroutine gridbilerp_s(nx,xsiz,ny,ysiz,qvalues,x,y,qlerp)
        integer, intent(in) :: nx,ny
        real, intent(in) :: xsiz,ysiz,x,y,qvalues(nx*ny)
        real, intent(out) :: qlerp
        integer :: ix,iy
        real :: x1,y1,x2,y2
        ! Get the lower left hand ix,iy coordinate
        ix = floor((x-0.5*xsiz)/xsiz) + 1
        iy = floor((y-0.5*ysiz)/ysiz) + 1
        ! Clip ix,iy to fall within 1,nx and 1,ny
        ix = max(min(ix,nx),1)
        iy = max(min(iy,ny),1)
        ! Is this an edge case (falls on the perimeter) or outside the grid? If so just return that grid value
        ! Check the x edges and corners
        if (x <= (0.5*xsiz)) then ! left hand edge
            if (y <= (0.5*ysiz)) then ! bottom left hand corner
                !write(*,*) 'bottom left hand corner'
                qlerp = qvalues(grididx(ix,iy,nx))
                return
            elseif (y >= (ny*ysiz-0.5*ysiz)) then ! top left hand corner
                !write(*,*) 'top left hand corner'
                qlerp = qvalues(grididx(ix,iy,nx))
                return
            else
                !write(*,*) 'left edge somewhere in the middle'
                y1 = iy*ysiz - 0.5*ysiz
                y2 = y1 + ysiz
                qlerp = lerp_s(y1,qvalues(grididx(ix,iy,nx)), &
                               y2,qvalues(grididx(ix,iy+1,nx)), y)
                return
            end if
        elseif (x >= (nx*xsiz-(0.5*xsiz))) then ! right hand edge
            if (y <= (0.5*ysiz)) then ! bottom right hand corner
                !write(*,*) 'bottom right hand corner'
                qlerp = qvalues(grididx(ix,iy,nx))
                return
            elseif (y >= (ny*ysiz-0.5*ysiz)) then ! top right hand corner
                !write(*,*) 'top right hand corner'
                qlerp = qvalues(grididx(ix,iy,nx))
                return
            else
                !write(*,*) 'right edge somewhere in the middle'
                y1 = iy*ysiz - 0.5*ysiz
                y2 = y1 + ysiz
                qlerp = lerp_s(y1,qvalues(grididx(ix,iy,nx)), &
                               y2,qvalues(grididx(ix,iy+1,nx)), y)
                return
            end if
        end if
        ! Check the top and bottom edges, corners have already been checked
        if (y <= (0.5*ysiz)) then ! bottom right hand corner
            !write(*,*) 'bottom edge somewhere in the middle'
            x1 = ix*xsiz - 0.5*xsiz
            x2 = x1 + xsiz
            qlerp = lerp_s(x1,qvalues(grididx(ix,iy,nx)), &
                           x2,qvalues(grididx(ix+1,iy,nx)), x)
            return
        elseif (y >= (ny*ysiz-0.5*ysiz)) then 
            !write(*,*) 'top edge somewhere in the middle'
            x1 = ix*xsiz - 0.5*xsiz
            x2 = x1 + xsiz
            qlerp = lerp_s(x1,qvalues(grididx(ix,iy,nx)), &
                           x2,qvalues(grididx(ix+1,iy,nx)), x)
            return
        end if
        ! This is not an edge case, perform bilinear interpolation
        !write(*,*) 'somewhere in the middle'
        x1 = ix*xsiz - 0.5*xsiz
        x2 = x1 + xsiz
        y1 = iy*ysiz - 0.5*ysiz
        y2 = y1 + ysiz
        qlerp = bilerp_s(x1,y1,x2,y2, &
                         qvalues(grididx(ix,iy,nx)), &
                         qvalues(grididx(ix,iy+1,nx)), &
                         qvalues(grididx(ix+1,iy,nx)), &
                         qvalues(grididx(ix+1,iy+1,nx)), x, y)
        return
    end subroutine gridbilerp_s

    ! Checks if a 2D gradient calculation would be valid at this location
    ! Returns 0 if invalid
    integer function checkgradient_d(nx,xsiz,ny,ysiz,qvalues,ix,iy,tmin) result(validgrad)
        integer, intent(in) :: nx,ny,ix,iy
        real*8, intent(in) :: xsiz,ysiz,qvalues(nx*ny),tmin
        integer testix,testiy
        validgrad = 1
        ! Valid center point?
        if (qvalues(grididx(ix,iy,nx)) .le. tmin) then
            validgrad = 0
            return
        end if
        ! Inside the grid?
        if (ix.lt.1 .or. iy.lt.1 .or. ix.gt.nx .or. iy.gt.ny) then
            validgrad = 0
            return
        end if
        ! Check around x
        testiy = iy
        testix = ix - 1
        if (testix.ge.1) then
            if (qvalues(grididx(testix,testiy,nx)) .le. tmin) then
                validgrad = 0
                return
            end if
        end if
        testix = ix + 1
        if (testix.le.nx) then
            if (qvalues(grididx(testix,testiy,nx)) .le. tmin) then
                validgrad = 0
                return
            end if
        end if
        ! Check around y
        testix = ix
        testiy = iy - 1
        if (testiy.ge.1) then
            if (qvalues(grididx(testix,testiy,nx)) .le. tmin) then
                validgrad = 0
                return
            end if
        end if
        testiy = iy + 1
        if (testiy.le.ny) then
            if (qvalues(grididx(testix,testiy,nx)) .le. tmin) then
                validgrad = 0
                return
            end if
        end if
        return
    end function checkgradient_d

    ! 2D gradient calculation at a point on a grid node (double precision)
    ! Note: Requires the 1-ordered grid node locations ix,iy for the derivatives
    ! Note: Uses the 3-point central difference method where possible
    subroutine gradient_d(nx,xsiz,ny,ysiz,qvalues,ix,iy,xderiv,yderiv)
        integer, intent(in) :: nx,ny,ix,iy
        real*8, intent(in) :: xsiz,ysiz,qvalues(nx*ny)
        real*8, intent(out) :: xderiv,yderiv
        real*8 :: x1,y1,x2,y2
        integer ixfixed,iyfixed
        ! Clip ix,iy to fall within 1,nx and 1,ny
        ixfixed = max(min(ix,nx),1)
        iyfixed = max(min(iy,ny),1)
        ! Compute the x derivative - prefer the 2 sided finite difference
        if (ixfixed .eq. 1) then ! 2-point forward difference
            xderiv = (qvalues(grididx(ixfixed+1,iyfixed,nx))-qvalues(grididx(ixfixed,iyfixed,nx)))/(xsiz)
        elseif (ixfixed .eq. nx) then ! 2-point backwards difference
            xderiv = (qvalues(grididx(ixfixed,iyfixed,nx))-qvalues(grididx(ixfixed-1,iyfixed,nx)))/(xsiz)
        else ! 3-point central difference
            xderiv = (qvalues(grididx(ixfixed+1,iyfixed,nx))-qvalues(grididx(ixfixed-1,iyfixed,nx)))/(2.0d0*xsiz)
        end if
        ! Compute the y derivative - prefer the 2 sided finite difference
        if (iyfixed .eq. 1) then ! 2-point forward difference
            yderiv = (qvalues(grididx(ixfixed,iyfixed+1,nx))-qvalues(grididx(ixfixed,iyfixed,nx)))/(ysiz)
        elseif (iyfixed .eq. ny) then ! 2-point backwards difference
            yderiv = (qvalues(grididx(ixfixed,iyfixed,nx))-qvalues(grididx(ixfixed,iyfixed-1,nx)))/(ysiz)
        else ! 3-point central difference
            yderiv = (qvalues(grididx(ixfixed,iyfixed+1,nx))-qvalues(grididx(ixfixed,iyfixed-1,nx)))/(2.0d0*ysiz)
        end if
        return
    end subroutine gradient_d

    ! 2D gradient calculation on a regular 2D grid (double precision)
    ! Note: Requires that x and y already be standardized (ie: x-xmin,y-ymin)
    subroutine gridgradient_d(nx,xsiz,ny,ysiz,qvalues,x,y,xderiv,yderiv)
        integer, intent(in) :: nx,ny
        real*8, intent(in) :: xsiz,ysiz,x,y,qvalues(nx*ny)
        real*8, intent(out) :: xderiv,yderiv
        integer :: ix,iy
        real*8 :: x1,y1,x2,y2
        ! Get the lower left hand ix,iy coordinate
        ix = floor((x-0.5d0*xsiz)/xsiz) + 1
        iy = floor((y-0.5d0*ysiz)/ysiz) + 1
        ! Clip ix,iy to fall within 1,nx and 1,ny
        ix = max(min(ix,nx),1)
        iy = max(min(iy,ny),1)
        ! Initialize derivatives
        xderiv = 0.0d0
        yderiv = 0.0d0
        ! Find out where we fall on the grid
        if (ix .lt. nx) then ! Not on x-boundary
            if (iy .lt. ny) then ! Not on y-boundary
                xderiv = 0.5d0*((qvalues(grididx(ix+1,iy,nx)) - qvalues(grididx(ix,iy,nx)))/xsiz + &
                                (qvalues(grididx(ix+1,iy+1,nx)) - qvalues(grididx(ix,iy+1,nx)))/xsiz)
                yderiv = 0.5d0*((qvalues(grididx(ix,iy+1,nx)) - qvalues(grididx(ix,iy,nx)))/ysiz + &
                                (qvalues(grididx(ix+1,iy+1,nx)) - qvalues(grididx(ix+1,iy,nx)))/ysiz)
            else ! On y-boundary
                xderiv = 0.5d0*((qvalues(grididx(ix+1,iy,nx)) - qvalues(grididx(ix,iy,nx)))/xsiz + &
                                (qvalues(grididx(ix+1,iy-1,nx)) - qvalues(grididx(ix,iy-1,nx)))/xsiz)
                yderiv = 0.5d0*((qvalues(grididx(ix,iy,nx)) - qvalues(grididx(ix,iy-1,nx)))/ysiz + &
                                (qvalues(grididx(ix+1,iy,nx)) - qvalues(grididx(ix+1,iy-1,nx)))/ysiz)
            end if
        else ! On x-boundary
            if (iy .lt. ny) then ! Not on y-boundary
                xderiv = 0.5d0*((qvalues(grididx(ix,iy,nx)) - qvalues(grididx(ix-1,iy,nx)))/xsiz + &
                                (qvalues(grididx(ix,iy+1,nx)) - qvalues(grididx(ix-1,iy+1,nx)))/xsiz)
                yderiv = 0.5d0*((qvalues(grididx(ix,iy+1,nx)) - qvalues(grididx(ix,iy,nx)))/ysiz + &
                                (qvalues(grididx(ix-1,iy+1,nx)) - qvalues(grididx(ix-1,iy,nx)))/ysiz)
            else ! On y-boundary
                xderiv = 0.5d0*((qvalues(grididx(ix,iy,nx)) - qvalues(grididx(ix-1,iy,nx)))/xsiz + &
                                (qvalues(grididx(ix,iy-1,nx)) - qvalues(grididx(ix-1,iy-1,nx)))/xsiz)
                yderiv = 0.5d0*((qvalues(grididx(ix,iy,nx)) - qvalues(grididx(ix,iy-1,nx)))/ysiz + &
                                (qvalues(grididx(ix-1,iy,nx)) - qvalues(grididx(ix-1,iy-1,nx)))/ysiz)
            end if
        end if
        return
    end subroutine gridgradient_d

    ! Number of cells in a line segment
    integer function numbercellsinline_d(x1,y1,x2,y2,cellsize) result(ncells)
        real*8, intent(in) :: x1,y1,x2,y2,cellsize
        ncells = nint(dsqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1))/cellsize) ! Use the nearest integer
    end function numbercellsinline_d

    ! 2D line segment discretization for a grid (double precision)
    subroutine discretizesegment(x1,y1,x2,y2,cellsize,ncells,cellcentersx,cellcentersy)
        real*8, parameter :: EPSLON = 1.0d-8
        real*8, intent(in) :: x1,y1,x2,y2,cellsize
        integer, intent(in) :: ncells
        real*8, intent(out) :: cellcentersx(ncells),cellcentersy(ncells)
        real*8 :: xinc,yinc,slope
        integer :: i
        ! Compute the increments for increasing x and y
        if (abs(x2-x1) > EPSLON*abs(x1)) then
            slope = (y2-y1)/(x2-x1)
            xinc = cellsize/dsqrt(slope*slope + 1.0d0)
            yinc = xinc*slope
        else
            xinc = 0d0
            yinc = cellsize*((y2-y1)/dsqrt((y2-y1)*(y2-y1)))
        end if
        ! write(*,*) dsqrt(xinc**2+yinc**2) <- should exactly equal cellsize
        ! First location is at x1+xinc/2,y1+yinc/2
        if (x1 .gt. x2) then
            cellcentersx(1) = x2 + xinc/2.0d0
            cellcentersy(1) = y2 + yinc/2.0d0
        else
            cellcentersx(1) = x1 + xinc/2.0d0
            cellcentersy(1) = y1 + yinc/2.0d0
        end if
        ! Fill in the rest of the locations
        do i=2,ncells
            cellcentersx(i) = cellcentersx(i-1) + xinc
            cellcentersy(i) = cellcentersy(i-1) + yinc
        end do
        return
    end subroutine discretizesegment

    ! Unit perpendicular vector
    subroutine unitperpvector(x1,y1,x2,y2,unitperpx,unitperpy)
        real*8, parameter :: EPSLON = 1.0d-8
        real*8, intent(in) :: x1,y1,x2,y2
        real*8, intent(out) :: unitperpx,unitperpy
        real*8 :: xinc,yinc,slope
        integer :: i
        ! Compute the increments for increasing x and y
        if (abs(x2-x1) > EPSLON*abs(x1)) then
            slope = (y2-y1)/(x2-x1)
            xinc = 1.0d0/dsqrt(slope*slope + 1.0d0)
            yinc = xinc*slope
        else
            xinc = 0d0
            yinc = (y2-y1)/dsqrt((y2-y1)*(y2-y1))
        end if
        ! Switch them around and multiply x by negative to get the perpendicular
        if (x1 .gt. x2) then
            unitperpx = yinc
            unitperpy = -1d0*xinc
        else
            unitperpx = -1d0*yinc
            unitperpy = xinc
        end if
        return
    end subroutine unitperpvector

    ! 1D grid extraction from a 2D line segment on a grid (double precision)
    ! Uses bilinear interpolation and computes a local gradient
    subroutine extractpoints_d(nx,xmin,xsiz,ny,ymin,ysiz,qvalues,ncells,pointsx,pointsy,values,gradientsx,gradientsy)
        integer, intent(in) :: nx,ny,ncells
        real*8, intent(in) :: xmin,xsiz,ymin,ysiz,qvalues(nx*ny),pointsx(ncells),pointsy(ncells)
        real*8, intent(out) :: values(ncells),gradientsx(ncells),gradientsy(ncells)
        integer :: icell
        real*8 :: x,y
        ! Get all the points and gradients
        do icell=1,ncells
            ! Normalize x and y for the subroutines
            x = pointsx(icell) - xmin
            y = pointsy(icell) - ymin
            ! Use bilinear interpolation to get a value
            call gridbilerp_d(nx,xsiz,ny,ysiz,qvalues,x,y,values(icell))
            ! Compute the gradient
            call gridgradient_d(nx,xsiz,ny,ysiz,qvalues,x,y,gradientsx(icell),gradientsy(icell))
        end do
        return
    end subroutine extractpoints_d

    ! Calculate the cumulative lengths of segments of an ordered polyline
    subroutine polycumulativelengths(nlinepts,linex,liney,linecumlengths)
        integer, intent(in) :: nlinepts
        real*8, intent(in) :: linex(nlinepts),liney(nlinepts)
        real*8, intent(out) :: linecumlengths(nlinepts-1)
        integer iline
        do iline=1,nlinepts-1
          linecumlengths(iline) = dsqrt((linex(iline+1)-linex(iline))*(linex(iline+1)-linex(iline)) + &
                                  (liney(iline+1)-liney(iline))*(liney(iline+1)-liney(iline)))
        end do
        do iline=2,nlinepts-1
          linecumlengths(iline) = linecumlengths(iline) + linecumlengths(iline-1)
        end do
    end subroutine polycumulativelengths

    ! Calculate the segment lengths of an ordered polyline
    subroutine polysegmentlengths(nlinepts,linex,liney,segmentlengths)
        integer, intent(in) :: nlinepts
        real*8, intent(in) :: linex(nlinepts),liney(nlinepts)
        real*8, intent(out) :: segmentlengths(nlinepts-1)
        integer iline
        do iline=1,nlinepts-1
          segmentlengths(iline) = dsqrt((linex(iline+1)-linex(iline))*(linex(iline+1)-linex(iline)) + &
                                  (liney(iline+1)-liney(iline))*(liney(iline+1)-liney(iline)))
        end do
    end subroutine polysegmentlengths

    ! Calculate the length of an ordered polyline
    real*8 function polylinelength(nlinepts,linex,liney) result(length)
        integer, intent(in) :: nlinepts
        real*8, intent(in) :: linex(nlinepts),liney(nlinepts)
        integer iline
        length = 0d0
        do iline=1,nlinepts-1
          length = length + dsqrt((linex(iline+1)-linex(iline))*(linex(iline+1)-linex(iline)) + &
                                  (liney(iline+1)-liney(iline))*(liney(iline+1)-liney(iline)))
        end do
    end function polylinelength

!-----------------------------------------------------------------------
! Discretize a polyline
!    - Input:
!       nlinepts - number of points composing the polyline (# of lines is 1 less than # of points)
!       linex(nlinepts) - x coordinates of line points
!       liney(nlinepts) - y coordinates of line points
!       ncells - number of grid cells (along x)
!       cellsize - grid cell size (along x)
!    - Returns:
!       foldedx(ncells) - folded x coordinates of cells
!       foldedy(ncells) - folded y coordinates of cells
!       segmentid(ncells) - segment id corresponding to polyline
!       unitperpx(ndata) - unit perpendicular vector x component
!       unitperpy(ndata) - unit perpendicular vector y component
!-----------------------------------------------------------------------
    subroutine discpolyline(nlinepts,linex,liney,ncells,cellsize,foldedx,foldedy,segmentid,unitperpx,unitperpy)
        real*8, parameter :: EPSLON = 1.0d-8
        integer, intent(in) :: nlinepts,ncells
        real*8, intent(in) :: linex(nlinepts),liney(nlinepts),cellsize
        real*8, intent(out) :: foldedx(ncells),foldedy(ncells),unitperpx(ncells),unitperpy(ncells)
        integer, intent(out) :: segmentid(ncells)
        ! Temporary variables
        real*8 :: linecumlengths(nlinepts-1),gridxprime,xprimemin,xprimenormed, &
                  segmentlengths(nlinepts-1),xinc,yinc,slope, &
                  lineperpx(nlinepts-1),lineperpy(nlinepts-1)
        integer iline,icell,nlines,currentline
        ! Get line cumulative lengths first
        call polycumulativelengths(nlinepts,linex,liney,linecumlengths)
        ! Segment lengths
        call polysegmentlengths(nlinepts,linex,liney,segmentlengths)
        ! Unit perpendiculars
        nlines = nlinepts-1
        do iline=1,nlines
            call unitperpvector(linex(iline),liney(iline),linex(iline+1),liney(iline+1), &
                                lineperpx(iline),lineperpy(iline))
        end do
        ! Offsetting
        xprimemin = cellsize*0.5d0
        currentline = 1
        ! Now step through the grid cells
        do icell=1,ncells
            ! Current location on xprime around polyline
            gridxprime = dble(icell)*cellsize - xprimemin
            ! Have we switched lines?
            if (gridxprime .gt. linecumlengths(currentline)) then
                do iline=currentline,nlines
                    if (gridxprime .lt. linecumlengths(iline)) then
                        exit
                    end if
                end do
                currentline = iline
            end if
            ! Where do we fall on the line?
            if (currentline .gt. 1) then
                xprimenormed = (gridxprime - linecumlengths(currentline-1))/segmentlengths(currentline)
            else
                xprimenormed = gridxprime/segmentlengths(currentline)
            end if
            ! Vector maths
            foldedx(icell) = xprimenormed*(linex(currentline+1) - linex(currentline)) + linex(currentline)
            foldedy(icell) = xprimenormed*(liney(currentline+1) - liney(currentline)) + liney(currentline)
            segmentid(icell) = currentline
            unitperpx(icell) = lineperpx(currentline)
            unitperpy(icell) = lineperpy(currentline)
        end do      
        return
    end subroutine discpolyline

!-----------------------------------------------------------------------
! Unfolds 2D point on to a line using the shortest perpendicular distance
!    - Input:
!       nlinepts - number of points composing the polyline (# of lines is 1 less than # of points)
!       linex(nlinepts) - x coordinates of line points
!       liney(nlinepts) - y coordinates of line points
!       ndata - number of data points to unfold
!       datax(ndata) - x coordinates of data
!       datay(ndata) - y coordinates of data
!    - Returns:
!       unfoldedx(ndata) - unfolded x' coordinates of data
!       unfoldedy(ndata) - unfolded y' coordinates of data
!       unfoldedvertx(nlinepts) - unfolded line vertex x' points
!       nearestx(ndata) - nearest x point on line for data
!       nearesty(ndata) - nearest y point on line for data
!-----------------------------------------------------------------------
    subroutine unfoldcoords(nlinepts,linex,liney,ndata,datax,datay,unfoldedx,unfoldedy,unfoldedvertx,nearestx,nearesty)
        integer, intent(in) :: nlinepts,ndata
        real*8, intent(in) :: linex(nlinepts),liney(nlinepts),datax(ndata),datay(ndata)
        real*8, intent(out) :: unfoldedx(ndata),unfoldedy(ndata),unfoldedvertx(nlinepts),nearestx(ndata),nearesty(ndata)
        ! Temporary variables
        integer iline,nline,idata,minidx
        real*8 :: unitperpx(nlinepts-1),unitperpy(nlinepts-1),perpdist(nlinepts-1), &
                  aline(nlinepts-1),bline(nlinepts-1),cline(nlinepts-1),sqrta2b2line(nlinepts-1), &
                  a2b2line(nlinepts-1),xptonline(nlinepts-1),yptonline(nlinepts-1), &
                  dotprod,sqlength,mindist,vertexdist(nlinepts),xprime,linelength(nlinepts-1)
        logical :: inbetween(nlinepts-1)
        ! Number of lines is 1 less than number of points
        nline = nlinepts-1
        ! Save the unfolded x coordinates
        unfoldedvertx(1) = 0d0
        ! Compute the unit perpendicular vectors and standard form
        do iline = 1,nline
            ! Unit perpendicular vector
            call unitperpvector(linex(iline),liney(iline),linex(iline+1),liney(iline+1),unitperpx(iline),unitperpy(iline))
            ! Standard form Ax + By + C = 0
            !     A    x +     B    y +       C       = 0
            ! (y1 – y2)x + (x2 – x1)y + (x1y2 – x2y1) = 0
            aline(iline) = liney(iline) - liney(iline+1)
            bline(iline) = linex(iline+1) - linex(iline)
            cline(iline) = linex(iline)*liney(iline+1) - linex(iline+1)*liney(iline)
            ! Precompute sqrt(a^2+b^2)
            sqrta2b2line(iline) = dsqrt(aline(iline)*aline(iline) + bline(iline)*bline(iline))
            ! Precompute a^2+b^2
            a2b2line(iline) = aline(iline)*aline(iline) + bline(iline)*bline(iline)
            ! Precompute the line length
            linelength(iline) = dsqrt((linex(iline+1)-linex(iline))*(linex(iline+1)-linex(iline)) + &
                               (liney(iline+1)-liney(iline))*(liney(iline+1)-liney(iline)))
            ! Save the unfolded x coordinates
            unfoldedvertx(iline+1) = unfoldedvertx(iline) + linelength(iline)
        end do
        ! Find the nearest line and compute the distance
        do idata = 1,ndata
            ! Compute all distances and nearest points
            do iline = 1,nline
                ! Compute the perpendicular distance
                perpdist(iline) = abs(aline(iline)*datax(idata) + &
                                      bline(iline)*datay(idata) + &
                                      cline(iline))/sqrta2b2line(iline)
                ! Find the closest point on the line
                xptonline(iline) = (bline(iline)*(bline(iline)*datax(idata) - &
                                    aline(iline)*datay(idata)) - &
                                    aline(iline)*cline(iline))/a2b2line(iline)
                yptonline(iline) = (aline(iline)*(-1d0*bline(iline)*datax(idata) + &
                                    aline(iline)*datay(idata)) - &
                                    bline(iline)*cline(iline))/a2b2line(iline)
                ! Check if the closest point is in between the line bounds
                inbetween(iline) = .true.
                dotprod = (xptonline(iline)-linex(iline))*(linex(iline+1)-linex(iline)) + &
                          (yptonline(iline)-liney(iline))*(liney(iline+1)-liney(iline))
                if (dotprod .lt. 0d0) then
                    inbetween(iline) = .false.
                else
                    sqlength = (linex(iline+1)-linex(iline))*(linex(iline+1)-linex(iline)) + &
                               (liney(iline+1)-liney(iline))*(liney(iline+1)-liney(iline))
                    if (dotprod .gt. sqlength) then
                        inbetween(iline) = .false.
                    end if
                end if
            end do
            ! Compute the distances to all vertices as well in case of convex line/edges
            do iline = 1,nlinepts
                vertexdist(iline) = dsqrt((datax(idata)-linex(iline))*(datax(idata)-linex(iline)) + &
                                          (datay(idata)-liney(iline))*(datay(idata)-liney(iline)))
            end  do
            mindist = 1.0d21
            ! Get the minimum distance that falls within the range of each line
            xprime = 0d0
            do iline = 1,nline
                if (inbetween(iline)) then
                    if (perpdist(iline) .lt. mindist) then
                        nearestx(idata) = xptonline(iline)
                        nearesty(idata) = yptonline(iline)
                        mindist = perpdist(iline)
                        ! Vertex x + distance along
                        unfoldedx(idata) = xprime + dsqrt((nearestx(idata)-linex(iline))*(nearestx(idata)-linex(iline)) + &
                                               (nearesty(idata)-liney(iline))*(nearesty(idata)-liney(iline)))
                        ! Test if dot product is positive for (+) or (-) y
                        if (((datax(idata)-xptonline(iline))*unitperpx(iline) + &
                            (datay(idata)-yptonline(iline))*unitperpy(iline)) .gt. 0d0) then
                            unfoldedy(idata) = perpdist(iline)
                        else
                            unfoldedy(idata) = -1d0*perpdist(iline)
                        end if
                    end if
                end if
                ! Increment the distance along the line
                xprime = xprime + linelength(iline)
            end do
            ! Now get the minimum distance to all vertices
            xprime = 0d0
            do iline = 1,nlinepts
                if (vertexdist(iline) .lt. mindist) then
                    nearestx(idata) = linex(iline)
                    nearesty(idata) = liney(iline)
                    mindist = vertexdist(iline)
                    unfoldedx(idata) = xprime ! Vertex x
                    ! Test if dot product is positive for (+) or (-) y
                    if (((datax(idata)-linex(iline))*unitperpx(min(iline,nline)) + &
                        (datay(idata)-liney(iline))*unitperpy(min(iline,nline))) .gt. 0d0) then
                        unfoldedy(idata) = vertexdist(iline)
                    else
                        unfoldedy(idata) = -1d0*vertexdist(iline)
                    end if
                end if
                ! Increment the length if not the first vertex
                if (iline .lt. nlinepts) then
                    xprime = xprime + linelength(iline)
                end if
            end  do
        end do
        return
    end subroutine unfoldcoords

end module linearinterp
