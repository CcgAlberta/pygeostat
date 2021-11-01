module sgvertices
!_____________________________________________________________________________
! Module for calculating the vertices a structured grid, given the 
! structured grid centroids. A simple interpolation/extrapolation 
! algorithm is used.
!   
! calc_sgvertices provides the public subroutine call, acting as an interface
! to the calculation of a 2-D structured surface(nz=1) or 3-D structured
! grid. The calc_sgvertices parameters are provided below.
!
! Input:
!   nxc: number of nodes in the x dimension
!   nyc: number of nodes in the y dimension
!   nzc: number of nodes in the z dimension (1 activates 2-D)
!   centroids(nxc*nyc*nzc, 3): coordinates of the centroids, which
!       follow a GSLIB-style convention in its ordering
!
! Output:
!   vertices((nx+1)*(ny+1)*(nz+1), 3)): coordinates of the calculated vertices
!   test: returns an error code (0 if no error)
!
! Version 1.0.0, Ryan M. Barnett, April 2018
!_____________________________________________________________________________
implicit none

! Public routines
public :: calc_sgvertices

private

integer*8 :: i, j 
! Dimensions of the structured grid cells
integer*8 :: nx, ny, nz
integer*8 :: ix, iy, iz
! Dimensions of the structured grid vertices
integer*8 :: nxv, nyv, nzv
integer*8 :: ixv, iyv, izv
! Working variables shared by the 2d and 3d subroutines
integer*8 :: ex, ey, ez, start_endp(2), start_endp1(2), start_end(2), start_end1(2), &
             extrap(2), extrap1(2)
 
contains

subroutine calc_sgvertices(nxc, nyc, nzc, i3d, centroids, vertices, test)
    ! Public interface for the 2-D and 3-D calculations
    ! Input:
    !   nxc: number of nodes in the x dimension
    !   nyc: number of nodes in the y dimension
    !   nzc: number of nodes in the z dimension (1 activates 2-D)
    !   centroids(nxc*nyc*nzc, 3): coordinates of the centroids, which
    !       follow a GSLIB-style convention in its ordering
    !   i3d: should be 0 if 2-D (nz = 1) or 1 if 3-D (nz>1). Must be provided 
    !       since the python interface requires immediate set dimensioning of the 
    !       output array.
    ! Output:
    !   vertices(:, 3)): coordinates of the calculated vertices, which are 
    !       (nx+1)*(ny+1) if 2-D, or (nx+1)*(ny+1)*(nz+1) if 3-D.
    !                    
    !   test: returns an error code (0 if no error)
    integer*8, intent(in) :: nxc, nyc, nzc, i3d
    real*8, dimension(nxc*nyc*nzc, 3), intent(in) :: centroids
    real*8, dimension((nxc+1)*(nyc+1)*(nzc+i3d), 3), intent(out) :: vertices
    integer, intent(out) :: test
    test = 0
    ! Checks
    if(nzc == 1)then
        if(i3d /= 0)then
            test = 1
            write(*,*) 'i3d must be 0 if nzc is 1'
            return
        end if
    else
        if(i3d /= 1)then
            test = 1
            write(*,*) 'i3d must be 1 if nzc greater than 1'
            return
        end if
    end if
    ! Iniitalize the global working variables
    nx = nxc
    nxv = nx + 1
    ny = nyc
    nyv = ny + 1
    if(nzc > 1)then
        nz = nzc
        nzv = nz + 1
        call sgvertices_3d(centroids, vertices, test)
    else
        call sgvertices_2d(centroids, vertices, test)
    end if
end subroutine calc_sgvertices
    
subroutine sgvertices_3d(centroids, vertices, test) 
    real*8, dimension(:, :), intent(in) :: centroids 
    real*8, dimension(:, :), intent(out) :: vertices
    integer, intent(out) :: test
    real*8, allocatable, dimension(:, :, :, :) :: cs  
    real*8, allocatable, dimension(:, :, :, :), target :: vs
    ! Pointer to invidual vertices
    real*8, dimension(:), pointer :: vt
    !
    ! Allocate the working arrays
    !
    test = 0
    allocate(vs(nxv, nyv, nzv, 3), cs(nx, ny, nz, 3), stat=test)
    if(test /= 0)then
        test = 97
        return
    end if
    i = 0
    do iz = 1,nz
        do iy = 1,ny
            do ix = 1,nx
                i = i + 1
                cs(ix, iy, iz, :) = centroids(i, :)
            end do
        end do
    end do

    !
    ! STEP ONE, vertices on the interior of the grid, which are shared
    ! by 8 grid cells - simply average the centroids
    !
    do izv = 2, nz
        do iyv = 2, ny
            do ixv = 2, nx
                vt => vs(ixv, iyv, izv, :)
                ! First, the four centroids in the plane below
                vt = cs(ixv-1, iyv-1, izv-1, :)
                vt = vt + cs(ixv, iyv-1, izv-1, :)
                vt = vt + cs(ixv-1, iyv, izv-1, :)
                vt = vt + cs(ixv, iyv, izv-1, :)
                ! Second, the four centroids in the plane above
                vt = vt + cs(ixv-1, iyv-1, izv, :)
                vt = vt + cs(ixv, iyv-1, izv, :)
                vt = vt + cs(ixv-1, iyv, izv, :)
                vt = vt + cs(ixv, iyv, izv, :)
                vt = vt/8.d0
            end do
        end do
    end do
    !
    ! STEP TWO, vertices along the planes (but not the edges/corners of the grid)
    ! Average the two varying values and extrapolate the constant value to its extent
    ! First, the x = 0 and nxv planes
    !
    start_endp = [int8(1), nxv]
    start_end = [int8(1), nx]
    extrap = [int8(2), nx-int8(1)]
    do i = 1, 2
        ixv = start_endp(i)
        ix = start_end(i)
        ex = extrap(i)
        do izv = 2, nz
            do iyv = 2, ny
                vt => vs(ixv, iyv, izv, :)
                vt = cs(ix, iyv-1, izv-1, :)
                vt = vt + (cs(ix, iyv-1, izv-1, :)-cs(ex, iyv-1, izv-1, :))*.5d0
                vt = vt + cs(ix, iyv, izv-1, :)
                vt = vt + (cs(ix, iyv, izv-1, :)-cs(ex, iyv, izv-1, :))*.5d0
                vt = vt + cs(ix, iyv-1, izv, :)
                vt = vt + (cs(ix, iyv-1, izv, :)-cs(ex, iyv-1, izv, :))*.5d0
                vt = vt + cs(ix, iyv, izv, :)
                vt = vt + (cs(ix, iyv, izv, :)-cs(ex, iyv, izv, :))*.5d0
                vt = vt/4.d0
            end do
        end do
    enddo
    ! Second, the y = 0 and nyv planes
    start_endp = [int8(1), nyv]
    start_end = [int8(1), ny]
    extrap = [int8(2), ny-int8(1)]
    do i = 1, 2
        iyv = start_endp(i)
        iy = start_end(i)
        ey = extrap(i)
        do izv = 2, nz
            do ixv = 2, nx
                vt => vs(ixv, iyv, izv, :)
                vt = cs(ixv-1, iy, izv-1, :)
                vt = vt + (cs(ixv-1, iy, izv-1, :)-cs(ixv-1, ey, izv-1, :))*.5d0
                vt = vt + cs(ixv, iy, izv-1, :)
                vt = vt + (cs(ixv, iy, izv-1, :)-cs(ixv, ey, izv-1, :))*.5d0
                vt = vt + cs(ixv-1, iy, izv, :)
                vt = vt + (cs(ixv-1, iy, izv, :)-cs(ixv-1, ey, izv, :))*.5d0
                vt = vt + cs(ixv, iy, izv, :)
                vt = vt + (cs(ixv, iy, izv, :)-cs(ixv, ey, izv, :))*.5d0
                vt = vt/4.d0
            end do
        end do
    enddo
    ! Third, the z = 0 and nzv planes
    start_endp = [int8(1), nzv]
    start_end = [int8(1), nz]
    extrap = [int8(2), nz-int8(1)]
    do i = 1, 2
        izv = start_endp(i)
        iz = start_end(i)
        ez = extrap(i)
        do iyv = 2, ny
            do ixv = 2, nx
                vt => vs(ixv, iyv, izv, :)
                vt = cs(ixv-1, iyv-1, iz, :)
                vt = vt + (cs(ixv-1, iyv-1, iz, :)-cs(ixv-1, iyv-1, ez, :))*.5d0
                vt = vt + cs(ixv, iyv-1, iz, :)
                vt = vt + (cs(ixv, iyv-1, iz, :)-cs(ixv, iyv-1, ez, :))*.5d0
                vt = vt + cs(ixv-1, iyv, iz, :)
                vt = vt + (cs(ixv-1, iyv, iz, :)-cs(ixv-1, iyv, ez, :))*.5d0
                vt = vt + cs(ixv, iyv, iz, :)
                vt = vt + (cs(ixv, iyv, iz, :)-cs(ixv, iyv, ez, :))*.5d0
                vt = vt/4.d0
            end do
        end do
    enddo
    !
    ! STEP THREE, vertices along the edges (but not the corners of the grid)
    ! Average the varying value and extrapolate the constant values to their extents
    ! First, the edges where only z varies (x/y = 0 and x/y= nxv/nyv)
    !
    start_endp = [int8(1), nxv]
    start_end = [int8(1), nx]
    extrap = [int8(2), nx-int8(1)]
    start_endp1 = [int8(1), nyv]
    start_end1 = [int8(1), ny]
    extrap1 = [int8(2), ny-int8(1)]
    do i = 1, 2
        ixv = start_endp(i)
        ix = start_end(i)
        ex = extrap(i)
        do j = 1, 2
            iyv = start_endp1(j)
            iy = start_end1(j)
            ey = extrap1(j)
            do izv = 2, nz
                vt => vs(ixv, iyv, izv, :)
                vt = cs(ix, iy, izv, :)
                vt = vt + (cs(ix, iy, izv, 1:2)-cs(ex, ey, izv, :))*.5d0
                vt = vt + cs(ix, iy, izv-1, :)
                vt = vt + (cs(ix, iy, izv-1, 1:2)-cs(ex, ey, izv-1, :))*.5d0
                vt = vt/2.d0
            end do
        end do
    enddo    
    ! Second, the edges where only y varies (x/z = 0 and x/z= nxv/nzv)
    start_endp = [int8(1), nxv]
    start_end = [int8(1), nx]
    extrap = [int8(2), nx-int8(1)]
    start_endp1 = [int8(1), nzv]
    start_end1 = [int8(1), nz]
    extrap1 = [int8(2), nz-int8(1)]
    do i = 1, 2
        ixv = start_endp(i)
        ix = start_end(i)
        ex = extrap(i)
        do j = 1,2
            izv = start_endp1(j)
            iz = start_end1(j)
            ez = extrap1(j)
            do iyv = 2,ny
                vt => vs(ixv, iyv, izv, :)
                vt = cs(ix, iyv-1, iz, :)
                vt = vt + (cs(ix, iyv-1, iz, :)-cs(ex, iyv-1, ez, :))*.5d0
                vt = vt + cs(ix, iyv, iz, :)
                vt = vt + (cs(ix, iyv, iz, :)-cs(ex, iyv, ez, :))*.5d0
                vt = vt/2.d0
            end do
        end do
    enddo
    ! Third, the edges where only x varies (y/z = 0 and y/z= nyv/nzv)
    start_endp = [int8(1), nyv]
    start_end = [int8(1), ny]
    extrap = [int8(2), ny-int8(1)]
    start_endp1 = [int8(1), nzv]
    start_end1 = [int8(1), nz]
    extrap1 = [int8(2), nz-int8(1)]
    do i = 1, 2
        iyv = start_endp(i)
        iy = start_end(i)
        ey = extrap(i)
        do j = 1, 2
            izv = start_endp1(j)
            iz = start_end1(j)
            ez = extrap1(j)
            do ixv = 2, nx
                vt => vs(ixv, iyv, izv, :)
                vt = cs(ixv-1, iy, iz, :)
                vt = vt + (cs(ixv-1, iy, iz, :)-cs(ixv-1, ey, ez, :))*.5d0
                vt = vt + cs(ixv, iy, iz, :)
                vt = vt + (cs(ixv, iy, iz, :)-cs(ixv, ey, ez, :))*.5d0
                vt = vt/2.d0
            end do
        end do
    enddo
    !
    ! STEP Four, vertices along the corners
    ! Simply extrapolate from the centroid
    !
    ! First, the z = 0 corners
    vs(1, 1, 1, :) = cs(1, 1, 1, :) + (cs(1, 1, 1, :)-cs(2, 2, 2, :))*.5d0
    vs(nxv, 1, 1, :) = cs(nx, 1, 1, :) + (cs(nx, 1, 1, :)-cs(nx-1, 2, 2, :))*.5d0
    vs(1, nyv, 1, :) = cs(1, ny, 1, :) + (cs(1, ny, 1, :)-cs(2, ny-1, 2, :))*.5d0
    vs(nxv, nyv, 1, :) = cs(nx, ny, 1, :) + &
                                 (cs(nx, ny, 1, :)-cs(nx-1, ny-1, 2, :))*.5d0
    ! Second, the z = nzv corners
    vs(1, 1, nzv, :) = cs(1, 1, nz, :) + (cs(1, 1, nz, :)-cs(2, 2, nz-1, :))*.5d0
    vs(nxv, 1, nzv, :) = cs(nx, 1, nz, :) + (cs(nx, 1, nz, :)-cs(nx-1, 2, nz-1, :))*.5d0
    vs(1, nyv, nzv, :) = cs(1, ny, nz, :) + (cs(1, ny, nz, :)-cs(2, ny-1, nz-1, :))*.5d0
    vs(nxv, nyv, nzv, :) = cs(nx, ny, nz, :) + &
                                 (cs(nx, ny, nz, :)-cs(nx-1, ny-1, nz-1, :))*.5d0
    ! Convert to the output dimensions
    i = 0
    do izv = 1,nzv
        do iyv = 1,nyv
            do ixv = 1,nxv
                i = i + 1
                vertices(i, :) = vs(ixv, iyv, izv, :)
            end do
        end do
    end do
end subroutine sgvertices_3d


subroutine sgvertices_2d(centroids, vertices, test) 
    real*8, dimension(:, :), intent(in) :: centroids 
    real*8, dimension(:, :), intent(out) :: vertices
    integer, intent(out) :: test    
    real*8, allocatable, dimension(:, :, :) :: cs 
    real*8, allocatable, dimension(:, :, :), target ::vs 
    ! Pointer to invidual vertices
    real*8, dimension(:), pointer :: vt
    !
    ! Allocate the working arrays
    !
    allocate(vs(nxv, nyv, 3),  cs(nx, ny, 3), stat=test)
    if(test /= 0)then
        test = 97
        return
    end if
    i = 0
    do iy = 1,ny
        do ix = 1,nx
            i = i + 1
            cs(ix, iy, :) = centroids(i, :)
        end do
    end do
    !
    ! STEP ONE, vertices on the interior of the grid, which are shared
    ! by 4 grid cells - simply average the centroids
    !
    do iyv = 2, ny
        do ixv = 2, nx
            vt => vs(ixv, iyv, :)
            vt = cs(ixv-1, iyv-1, :)
            vt = vt + cs(ixv, iyv-1, :)
            vt = vt + cs(ixv-1, iyv, :)
            vt = vt + cs(ixv, iyv, :)
            vt = vt/4.d0
        end do
    end do
    !
    ! STEP TWO, vertices along the edges (but not the corners of the grid)
    ! Average the varying value and extrapolate the constant values to their extents 
    ! First, the edges where only y varies (x = 0 and x= nxv)
    start_endp = [int8(1), nxv]
    start_end = [int8(1), nx]
    extrap = [int8(2), nx-int8(1)]
    do i = 1, 2
        ixv = start_endp(i)
        ix = start_end(i)
        ex = extrap(i)
        do iyv = 2,ny
            vt => vs(ixv, iyv, :)
            vt = cs(ix, iyv-1, :)
            vt = vt + (cs(ix, iyv-1, :)-cs(ex, iyv-1, :))*.5d0
            vt = vt + cs(ix, iyv, :)
            vt = vt + (cs(ix, iyv, :)-cs(ex, iyv, :))*.5d0
            vt = vt/2.d0
        end do
    enddo
    ! Second, the edges where only x varies (y = 0 and y= nyv)
    start_endp = [int8(1), nyv]
    start_end = [int8(1), ny]
    extrap = [int8(2), ny-int8(1)]
    do i = 1, 2
        iyv = start_endp(i)
        iy = start_end(i)
        ey = extrap(i)
        do ixv = 2, nx
            vt => vs(ixv, iyv, :)
            vt = cs(ixv-1, iy, :)
            vt = vt + (cs(ixv-1, iy, :)-cs(ixv-1, ey, :))*.5d0
            vt = vt + cs(ixv, iy, :)
            vt = vt + (cs(ixv, iy, :)-cs(ixv, ey, :))*.5d0
            vt = vt/2.d0
        end do
    enddo
    !
    ! STEP Three, vertices along the corners
    ! Simply extrapolate from the centroid
    !
    vs(1, 1, :) = cs(1, 1, :) + (cs(1, 1, :)-cs(2, 2, :))*.5d0
    vs(nxv, 1, :) = cs(nx, 1, :) + (cs(nx, 1, :)-cs(nx-1, 2, :))*.5d0
    vs(1, nyv, :) = cs(1, ny, :) + (cs(1, ny, :)-cs(2, ny-1, :))*.5d0
    vs(nxv, nyv, :) = cs(nx, ny, :) + &
                                 (cs(nx, ny, :)-cs(nx-1, ny-1, :))*.5d0
    ! Convert to output
    i = 0
    do iyv = 1,nyv
        do ixv = 1,nxv
            i = i + 1
            vertices(i, :) = vs(ixv, iyv, :)
        end do
    end do
end subroutine sgvertices_2d

end module sgvertices