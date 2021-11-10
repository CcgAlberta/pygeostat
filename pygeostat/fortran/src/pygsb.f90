! gsb_interface.f90 - contains routines for interfacing with the GSB file format
!
! (c) 2016 Jared Deutsch

subroutine pyreadgsbheader(gsbfl, nvar, nx, ny, nz, nreal, errorvalue)
!-----------------------------------------------------------------------------
! Load in GSB header from a file
!   Parameters:
!     gsbfl (char*512) - GSB file name
!   Returns:
!     nvar (int) - number of variables in the file
!     nx (int) - dim(1) in GSB file
!     ny (int) - dim(2) in GSB file
!     nz (int) - dim(3) in GSB file
!     nreal (int) - number of realizations
!     errorvalue (int) - error value reported from gslib_binary
!
! .. codeauthor:: Jared Deutsch - 2016-02-19
!-----------------------------------------------------------------------------
  use gslib_binary
  implicit none

  character*512, intent(in) :: gsbfl
  integer, intent(out) :: nvar, nx, ny, nz, nreal, errorvalue

  type(gsb_info) :: gsb_in
  type(gsb_var), allocatable :: gsbv_in(:)
  integer :: lin

  ! Collect header information - and keep the errorcode
  errorvalue = read_bin_header( lin, gsbfl, gsb_in, gsbv_in )
  if(errorvalue == 0)then
    nx = gsb_in%dims(1)
    ny = 1
    nz = 1
    if(gsb_in%ndim > 1) ny = gsb_in%dims(2)
    if(gsb_in%ndim > 1) nz = gsb_in%dims(3)
    nvar = gsb_in%nvar
    nreal = gsb_in%nreal
  else
    print '(a,x,i0)', 'pygsbreadheader failed with code: ',errorvalue
  endif
  ! Cleanup
  close(lin)
  if(allocated(gsbv_in)) deallocate(gsbv_in)

end subroutine pyreadgsbheader


subroutine pyreadgsbdata(gsbfl, nvar, nxyz, ireal, nullval, errorvalue, datamat, vnames)
!-----------------------------------------------------------------------------
! Load in GSB data from a file
!   Parameters:
!     gsbfl (char*512) - GSB file name
!     nvar (int) - number of variables in the file
!     nxyz (int) - number of records (nx*ny*nz)
!     ireal (int) - realization number to read
!     nullval (dble) - the value to fill the matrix with
!   Returns:
!     datamat(nvar, nxyz) - matrix of data values
!     errorvalue (int) - error value reported from gslib_binary
!     vnames(50) (char*64) - variable names reported from gslib_binary
!
! .. codeauthor:: Jared Deutsch - 2016-02-19
!-----------------------------------------------------------------------------
  use gslib_binary
  implicit none

  character*512, intent(in) :: gsbfl
  integer, intent(in) :: nxyz, nvar, ireal
  real*8, intent(in) :: nullval
  integer, intent(out) :: errorvalue
  real*8, dimension(nvar, nxyz), intent(out) :: datamat
  character*64, dimension(nvar), intent(out) :: vnames

  integer :: i, lin
  integer, dimension(nvar) :: vcols
  type(gsb_info) :: gsb_in
  type(gsb_var), allocatable :: gsbv_in(:)

  ! Header information
  errorvalue = read_bin_header( lin, gsbfl, gsb_in, gsbv_in )

  ! Grab variable names
  do i=1,gsb_in%nvar
    vnames(i) = gsb_in%vnames(i)
  end do

  ! Check for errors
  datamat = nullval
  call set_gsb_null_value( nullval )
  if ((errorvalue .ne. 0) .or. (gsb_in%nreal < ireal)) then
    if (gsb_in%nreal < ireal) then
      errorvalue = -987
      write(*,*) 'ERROR: ireal > nreal',ireal,gsb_in%nreal
    end if
    return
  end if

  ! Reset the error value and seed the column numbers
  errorvalue = 0
  do i=1,nvar
    vcols(i) = i
  enddo

  ! Get the correct realization and make sure that keyouts are
  ! correctly read
  do i=1,ireal
    ! Read it in
    errorvalue = read_bin_data( lin, gsb_in, gsbv_in, vcols, ireal, 1 )
    if(errorvalue /= 0) exit
    ! Get the data matrix
    errorvalue = gsbv2array( gsb_in, gsbv_in, vcols, datamat )
    if(errorvalue /= 0) exit
  end do

  ! Cleanup
  close(lin)
  deallocate(gsbv_in)

end subroutine pyreadgsbdata


subroutine pywritegsbdata(gsbfl, datamat, vnames, nvar, nreal, nx, ny, nz, tmin, &
                          tmax, tvar, vkinds, errorvalue)
!-----------------------------------------------------------------------------
! Write out GSB data to a file
!   Parameters:
!     gsbfl (char*512) - output GSB file name
!     datamat(nvar, nxyz) - matrix of data values
!     vnames(50) (char*64) - variable names for gslib_binary
!     nvar (int) - number of variables in the file
!     nx (int) - number of cells in x direction
!     ny (int) - number of cells in y direction
!     nz (int) - number of cells in z direction
!     tmin (float) - lower trimming limit
!     tmax (float) - upper trimming limit
!     tvar (int) - uses keyout trimming with tmin/tmax based on variable
!                  tvar if tvar > 0
!     vkinds (nvar) - format numbers for each variable (1=i2,2=r4,3=r8)
!   Returns:
!     errorvalue (int) - error value reported from gslib_binary
!
! .. codeauthor:: Jared Deutsch - 2016-02-19
!-----------------------------------------------------------------------------
  use gslib_binary
  implicit none

  character*512, intent(in) :: gsbfl
  integer, intent(in) :: nx, ny, nz, nvar, nreal, tvar
  integer, dimension(nvar), intent(in) :: vkinds
  real*8, intent(in) :: tmin, tmax
  real*8, dimension(nvar, nreal*nx*ny*nz), intent(in) :: datamat
  character*64, dimension(nvar), intent(in) :: vnames
  integer, intent(out) :: errorvalue

  integer :: i, lout, istart, ifin, j
  integer, dimension(nvar) :: vcols
  type(gsb_info) :: gsb_out
  type(gsb_var), allocatable :: gsbv_out(:)

  ! Assign the required output GSB parameters
  gsb_out%header = 'pygeostat output data'
  gsb_out%nvar = nvar
  allocate(gsb_out%vnames(nvar), gsb_out%kinds(nvar))
  gsb_out%vnames = vnames
  do i=1,nvar
    !vkinds(i) = 3
    vcols(i) = i
  end do
  gsb_out%kinds = vkinds ! Variable kinds (1=i2,2=r4,3=r8)
  gsb_out%ndim = 3 ! Number of dimensions (not counting realizations)
  gsb_out%dims = [nx,ny,nz] ! Dimension lengths
  gsb_out%nreal = nreal ! Number of realizations
  if (tvar .gt. 0) then
    gsb_out%ikeyout = 1 ! Use trimming limits for keyout?
  else
    gsb_out%ikeyout = 0 ! Use trimming limits for keyout?
  end if
  gsb_out%ikeyout_real = 1 ! If keyout, use variable trimming per real
  gsb_out%tmin = tmin ! Minimum trimming limit
  gsb_out%tmax = tmax ! Maximum trimming limit
  gsb_out%tvar = tvar ! Variable for trimming

  ! Write the header information and prepare arrays
  ! Note that gsb_out and gsbv_out are allocated within the subroutine
  errorvalue = write_bin_header( lout, gsbfl, gsb_out, gsbv_out )

  ! Write out the data
  istart = 1
  do i = 1, nreal
    ifin = i * nx * ny * nz
    errorvalue =  array2gsbv( gsb_out, gsbv_out, vcols, datamat(:, istart:ifin) )
    errorvalue = write_bin_data( lout, gsb_out, gsbv_out, i )
    istart = ifin + 1
  enddo

  ! Cleanup
  close(lout)
  deallocate(gsbv_out)

end subroutine pywritegsbdata
