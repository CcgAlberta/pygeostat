! public
      subroutine getcollocated(datafl,nx,xmn,xsiz,ny,ymn,ysiz,nz, &
                               zmn,zsiz,ndat,xyz,nsecv,seccols, &
                               secdat)
!---------------------------------------------------------------------
!
!               Get Collocated Exhaustive Secondary Data
!               ****************************************
!
!  Description
!  -----------
!  The subroutine grabs the collocated exhaustive secondary data based
!  on the xyz coordinates passed and returns it as an array.
!
!  Compile using the command:
!  f2py -c -m --compiler=cygwin getcollocated getcollocated.f90 sortem.f90
!
!  Parameters
!  ----------
!  datafl: Path to gridded data
!  nx, zsiz, ect.: Grid definition parameters
!  ndat: Number of locations to get data from
!  xyz: Array of the xyz coordinates
!  nsecv: Number of secondary variables to return
!  seccols: Columns which contain the secondary variables to retrieve
!  secdat: Array in which the retrieved secondary variables are
!          returned
!
!  Contributors
!  ------------
!  Author  - Warren E. Black                   DATE: February 22, 2016
!
!  Revised - Name                              DATE:
!            Description of revisions.
!
!  (c) 2016, Warren E. Black
!---------------------------------------------------------------------

      use sortem

      implicit none

      character, intent(in)     :: datafl*512
      integer, intent(in)       :: nx, ny, nz, ndat, nsecv
      real(kind=8), intent(in)  :: xmn, ymn, zmn, xsiz, ysiz, zsiz
      real(kind=8), intent(in)  :: seccols(:), xyz(:,:)
      real(kind=8), intent(out) :: secdat(ndat,nsecv)
      integer                   :: ix, iy, iz, loc, i, j, k, l, m, ncell, nvar, test
      logical                   :: inflag
      real(kind=8)              :: gindex(ndat, 2)
      real(kind=8),allocatable  :: var(:)

      ! Set the default value of the secondary data to -999
      secdat=-999

      ! Read over the header
      open(20,file=datafl,action='read',status='old')
      read(20,*)
      read(20,*) nvar
      j = 1
      do i=1, nvar
            read(20, *)
      end do

      ! Allocate a few arrays
      allocate(var(nvar),stat = test)
      if (test.ne.0) stop 'Error:  Allocation failed'

      ! Get an array of the indexes needed, -999 indicates that the
      ! xyz is outside of the grid
      do i=1, ndat
         call getindx(nx,xmn,xsiz,xyz(i, 1),ix,inflag)
         if (inflag .eqv. .False.) then
            if (nx .eq. 1) then
               ix = 1
            else
               ix = -999
            endif
         end if
         call getindx(ny,ymn,ysiz,xyz(i, 2),iy,inflag)
         if (inflag .eqv. .False.) then
            if (ny .eq. 1) then
               iy = 1
            else
               iy = -999
            endif
         end if
         call getindx(nz,zmn,zsiz,xyz(i, 3),iz,inflag)
         if (inflag .eqv. .False.) then
            if (nz .eq. 1) then
               iz = 1
            else
               iz = -999
            endif
         end if
         ! Get the 1-d index location
         loc = (iz-1)*nx*ny+(iy-1)*nx+ix
         gindex(i, 1) = loc
         gindex(i, 2) = i
      end do

      ! Sort the array
      call dblemodsortem(gindex(:,1),ndat,1,gindex(:,2))

      ! Grab the data
      j = 1
      ncell = nx*ny*nz
      do i=1, ncell
         ! If the xyz coord is outside of the grid, skip it
         if (gindex(j, 1) .lt. 0) then
            j = j + 1
            cycle
         endif
         if (i .eq. gindex(j, 1)) then
            read(20,*,err=99) (var(l),l=1,nvar)
            l=1
            ! Only record the desired columns
            do k=1,nsecv
               if (seccols(l) .eq. k) then
                  m = gindex(j, 2)
                  secdat(m,l) = var(k)
               endif
               l = l +1
            enddo
            j = j +1
         else
            read(20, *)
         endif
      enddo

      ! Wrap it up
      close(20)
      return
 99   stop ' ERROR in data file'
      end subroutine getcollocated

  subroutine getindx(n,min,siz,loc,index,inflag)
!---------------------------------------------------------------------
!
!     Gets the coordinate index location of a point within a grid
!     ***********************************************************
!
!
! n       number of "nodes" or "cells" in this coordinate direction
! min     origin at the center of the first cell
! siz     size of the cells
! loc     location of the point being considered
! index   output index within [1,n]
! inflag  true if the location is actually in the grid (false
!         otherwise e.g., if the location is outside then index will
!         be set to nearest boundary
!
!---------------------------------------------------------------------
      integer, parameter    :: dp=kind(0.d0)
      integer               :: n,index
      real(dp)              :: min,siz,loc
      logical               :: inflag
!
! Compute the index of "loc":
!
      index = int( (loc-min)/siz + 1.5 )
!
! Check to see if in or out:
!
      if(index.lt.1) then
            index  = 1
            inflag = .false.
      else if(index.gt.n) then
            index  = n
            inflag = .false.
      else
            inflag = .true.
      end if
!
! Return to calling program:
!
      return
      end subroutine getindx
