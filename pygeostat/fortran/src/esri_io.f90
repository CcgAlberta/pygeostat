! public
module esri_io

      implicit none
      real*8, allocatable, dimension(:,:), public :: data_double(:, :)

contains
    subroutine readgrid(ascfl, xmn, ymn, zmn, xsiz, ysiz, zsiz, nx, ny, nz)
    !--------------------------------------------------------------------!
    !               Convert ESRI ASCII GRID to GSLIB Grid                !
    !               *************************************                !
    !                                                                    !
    !        File: arcgistogslib.f90                                     !
    !      Author: Warren E. Black                                       !
    !     Created: 02/25/15                                              !
    ! Description: Converts ASCII grid output from ArcGIS to GSLIB grid  !
    !              format                                                !
    !                                                                    !
    !--------------------------------------------------------------------!

        character, intent(in)   :: ascfl*512
        real*8, intent(out)     :: xmn, ymn, zmn, xsiz, ysiz, zsiz
        integer*8, intent(out)  :: nx, ny, nz
        integer                 :: nan, nsim, test, i, j, l, err, lin
        character(len=70)       :: dump

        !---------------------------------------------------------------------
        !                            Global Variables
        !---------------------------------------------------------------------

        nz=1
        zmn=0
        zsiz=1
        nsim=1

        !---------------------------------------------------------------------
        !                              Read Files
        !---------------------------------------------------------------------

        !Read ArcGIS output
        lin = 1
        open(lin, file=ascfl)
              !Read header
              read(lin,*) dump,nx
              read(lin,*) dump,ny
              read(lin,*) dump,xmn
              read(lin,*) dump,ymn
              read(lin,*) dump,xsiz
              ysiz=xsiz
              xmn=xmn+(xsiz/2)
              ymn=ymn+(ysiz/2)
              read(lin,*) dump,nan
              !Allocate the memory
              if( allocated(data_double) ) deallocate(data_double)
              allocate( data_double(nx, ny), stat=err )
              if (err .ne. 0 ) stop ' ERROR: memory allocation failure! - data_double(nx, ny)'
              !Read data
              read(lin,*,end=10,err=99) data_double
        close(lin)
        return
 10     write(*,'(a,a)') 'ERROR, END OF FILE: ', trim(ascfl)
        close(lin)
        return
 99     write(*,'(a,a)') 'ERROR READING FILE: ',trim(ascfl)
        close(lin)
        return
        return
    end subroutine readgrid

end module esri_io
