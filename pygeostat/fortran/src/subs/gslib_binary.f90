
!------------------------------------------------------------------------------------------------
! Module for reading and writing GSLIB-style binary data in compressed binary format (GSB)
!
! Dependencies: filehandling module
!
! Modified: Ryan M. Barnett, March 2018   (Version 5.0.1)
! - 5.0.0 introduced some minor issues with backwards compatability of old program
!   implementations. These are now corrected/handled.
!
! Modified: John G. Manchuk, December 2017   (Version 5.0.0)
! - implemented new open, read, write interfaces to simplify usage
!   and permit more generalized error checking.
!
! Modified: Ryan M. Barnett, November 2017   (Version 4.1.0)
! - increased the capped number of variables to 1000 from 50, initialization of vnames avoids
!   null characters
! Modified: Ryan M. Barnett, October 2017   (Version 4.004)
! - write in a single statement, using a temporary array to avoid stack overflow issues
! Modified: Ryan M. Barnett, September 2017
! - converted several subroutines to functions for improved error handling
! - added gsbv2array routines for converting to and from GSB objects
! Modified: Ryan M. Barnett, January 2016
! - placed the write within a loop to avoid large indexes that can lead to stack overflow
! Modified: Ryan M. Barnett, May 2015
! - added functionality, annotation, and cleaned up many routines
! Author: Ryan M. Barnett, 2014
!
! Users should run the demo-gslib_binary program to understand the implemenation of this module
! within a calling program.
!------------------------------------------------------------------------------------------------

module gslib_binary

    use filehandling, only : file_open, str2lower
    implicit none

!-------------------------------------------------------------------------------------------
! Public routines and variables
!-------------------------------------------------------------------------------------------
    ! Interface routines by JGM
    public :: opengsb, &            !Prep a gsb file for writing (calls write_bin_header and builds gsb type)
              gsbwrite, &           !Write data to a GSB file (calls write_bin_data)
              gsbread               !Read data from a GSB file (calls read_bin_data)

    public :: write_bin_data, &     ! Write data
              read_bin_data, &      ! Read data
              gsbv2array, &         ! Convert GSBV object to matrix array
              array2gsbv, &         ! Convert matrix array to GSBV object
              checkGSB, &           ! Check if a file has a GSB extension
              killgsb, &            ! Deallocate/nullify the gsb and gsbv arrays
              set_gsb_null_value    ! Set the null site fill value, default = -999.99d0


    ! Deprecated
    public :: write_bin_header, &   ! Write a binary header
              read_bin_header       ! Read a binary header

    public :: gsb_info, gsb_var

    ! Interface for reading binary data - one (vcol) or multiple (vcols) variables
    interface read_bin_data
        module procedure read_bin_data_vcol, &
                         read_bin_data_vcols
    end interface read_bin_data

    ! Interface for converting different kinds of gsbv data to a single matrix or vector array.
    ! May be useful for accumulating data for operations after reading in the variables within gsbv
    ! objects. The interface allows for a single subroutine call (gsbv2array)  to be used with
    ! output arrays of i) one-dimension (vector) or two dimensions (matrix), and ii) integer*2,
    ! real*4 or real*8 formats.
    interface gsbv2array
        module procedure gsbv2vector_i2, &
                         gsbv2vector_r4, &
                         gsbv2vector_r8, &
                         gsbv2matrix_i2, &
                         gsbv2matrix_r4, &
                         gsbv2matrix_r8
    end interface gsbv2array

    ! Interface for converting a single matrix or vector array into different kinds of gsbv data arrays.
    ! May be useful for distributing data into the required gsbv objects prior to writing out.
    ! The interface allows for a single subroutine call (array2gsbv)  to be used with
    ! input data arrays of i) one-dimension (vector) or two dimensions (matrix), and ii) integer*2,
    ! real*4 or real*8 formats.
    interface array2gsbv
        module procedure vector2gsbv_i2, &
                         vector2gsbv_r4, &
                         vector2gsbv_r8, &
                         matrix2gsbv_i2, &
                         matrix2gsbv_r4, &
                         matrix2gsbv_r8
    end interface array2gsbv

    private

    ! Value to fill keyed out sites on array reconstruction from GSB objects
    real*8, private :: ASCII_NULL_VALUE = -999.99d0

    ! Data array for a variable (3 available kinds), run length encoding if i2, etc.
    type gsb_var
        ! Users do not need to know anyting about these variables if using the gsbv2array/array2gsbv routines
        integer*4 :: nrle                       ! Length of the rle arrays
        integer*4, allocatable :: rle(:)        ! Run length encoding runs
        integer*4, allocatable :: rle_index(:)  ! Index of data_ value for each run
        integer*2, allocatable :: data_i2(:)    ! Integer type 2 storage
        real*4, allocatable :: data_r4(:)       ! Real kind=4 storage
        real*8, allocatable :: data_r8(:)       ! Real kind=8 storage
    end type gsb_var

    ! Stores header information, variable names, variable kinds, keyout array, etc
    type gsb_info
        integer :: fid                  ! File handle
        character*512 :: filename       ! Name of associated file
        character*512 :: header         ! Arbitrary text information about the file
        integer*2 :: nvar               ! Number of variables in the file
        character*64, allocatable :: vnames(:)    ! Names of the variables
        integer*2, allocatable :: kinds(:)        ! Kind code of the variables (1=i2,2=r4,3=r8) (previously called fmts)
        integer*2 :: ndim               ! Number of dimensions for all variables
        integer*4 :: dims(3)            ! Size of each dimension
        integer*4 :: nrec               ! Total number of records (equals multiplication of dimensions)
        integer*2 :: ikeyout            ! Is a keyout array used?
        real*4 :: tmin                  ! Minimum trimming limits to use for masking
        real*4 :: tmax                  ! Minimum trimming limits to use for masking
        integer*2 :: tvar               ! Variable to use for trimming/masking
        integer*4 :: nreal              ! Number of realizations
        integer*2 :: ikeyout_real       ! Is a unique keyout used for each realization?
        integer*4 :: nkeyout                    ! Length of the keyout array
        integer*4 :: nko_index                  ! Length of the keyout index
        logical,  allocatable :: mask(:)        ! TODO UNDEFINED???
        integer*4, allocatable :: keyout(:)     ! Keyout run length array
        integer*4, allocatable :: ko_index(:)   ! Index of non-null data_ values
        type(gsb_var), allocatable :: gsbv(:)  ! Data storage
    end type gsb_info

    ! Temporary storage arrays
    logical, allocatable, target, dimension(:) :: ltemp11
    integer*2, allocatable, dimension(:) :: i2temp11, i2temp12, i2temp13, i2temp14
    real*4, allocatable, dimension(:) :: r4temp11, r4temp12
    real*8, allocatable, dimension(:) :: r8temp11, r8temp12

    integer, parameter :: INT2B = 1, REAL4B = 2, REAL8B = 3

    contains

    !----------------------------------------------------------------------------------------
    ! Set (or reset, noargs) the null value to fill keyed out sites with
    !----------------------------------------------------------------------------------------
    subroutine set_gsb_null_value( nullval )
        real*8, intent(in), optional :: nullval
        if(present(nullval))then
            ASCII_NULL_VALUE = nullval
        else
            ASCII_NULL_VALUE = -999.99d0
        endif
    end subroutine set_gsb_null_value

    !----------------------------------------------------------------------------------------
    ! The following function tests if a file extension (in a character string) is .gsb%
    ! This may be useful for determining if a data file is a GSB format (rather than ASCII
    ! or conventional binary format) based on a character string from a parameter file.
    !----------------------------------------------------------------------------------------
    logical function checkgsb ( fname )
        character*(*), intent (in) :: fname
        character*3 :: test
        integer :: los
        checkgsb = .false.
        los = index(fname,'.',back=.true.)
        if(los < 1)then !No extension case
            return
        endif
        test = fname(los+1:los+3)
        call str2lower(test)
        if( test == 'gsb' )then
            checkgsb = .true.
        endif
    end function checkGSB

!-------------------------------------------------------------------------------------------
! Routines for writing a binary file
!-------------------------------------------------------------------------------------------

    !Open a GSB file and initialize the necessary GSB type - this wrapper was implemented
    !by JGM to facilitate more user friendly initialization/usage of gsb files and allow
    !some error checking.
    integer function opengsb ( gsb, fname, rw, nvar, vtypes, dims, header, vnames, nsim, &
                               trimkey, varkey, tmin, tmax, tvar ) result ( fid )
        type(gsb_info), intent (inout) :: gsb
        character*(*), intent (in) :: fname
        character*1, intent (in) :: rw !read 'r', or write 'w'
        integer, intent (in), optional :: nvar
        integer, intent (in), optional :: vtypes(:)
        integer, intent (in), optional :: dims(:)

        !Defaults assigned if optional variables not passed in
        character*(*), intent (in), optional :: header
        character*(*), intent (in), optional :: vnames(:)
        integer, intent (in), optional :: nsim, trimkey, varkey
        real, intent (in), optional :: tmin, tmax
        integer, intent (in), optional :: tvar

        character*16 :: sstr
        integer :: i

        fid = -1 !Indicates error condition - file was not opened

        gsb%filename = fname

        if(rw == 'w')then !Check for required inputs in this case
            if(.not. present(nvar) .or. &
               .not. present(vtypes) .or. &
               .not. present(dims))then
                write(*,'(/,a)') 'nvar, vtypes and dims are required for opengsb in write mode'
                return
            endif
        elseif(rw == 'r')then
            i = read_bin_header( fid, fname, gsb, gsb%gsbv )
            if(i /= 0)then
                fid = -1
            endif
            gsb%fid = fid
            return
        else
            write(*,'(/,a)') 'opengsb modes are r (read) or w (write)'
            return
        endif

        if(present(header))then
            gsb%header = trim(header)               !Header information
        else
            gsb%header = 'GSB OUTPUT'
        endif

        if(nvar < 1)then
            write(*,'(/,a)') 'opengsb requires nvar > 0'
            return
        endif
        gsb%nvar = nvar                           !Number of variables

        if(ubound(vtypes,1) > 1)then              !If 1 type specified, assumes all variables are same type
            if(ubound(vtypes,1) /= nvar)then
                write(*,'(/,a)') 'length(vtypes) /= nvar in opengsb'
                return
            endif
        endif

        if(present(vnames))then
            if(ubound(vnames,1) /= nvar)then
                write(*,'(/,a)') 'length(vnames) /= nvar in opengsb'
                return
            endif
        endif

        if(allocated(gsb%vnames)) deallocate(gsb%vnames)
        allocate(gsb%vnames(nvar),stat=i)
        if(i /= 0)then
            write(*,'(/,a)') 'allocation failure in opengsb (gsb%vnames)'
            return
        endif

        if(present(vnames))then
            do i = 1, nvar                          !Variable names
                gsb%vnames(i) = trim(vnames(i))
            enddo
        else
            do i = 1, nvar                          !default Variable names
                write(sstr,'(i5.5)') i
                gsb%vnames(i) = 'variable_' // sstr(1:5)
            enddo
        endif

        if(allocated(gsb%kinds)) deallocate(gsb%kinds)
        allocate(gsb%kinds(nvar),stat=i)
        if(i /= 0)then
            write(*,'(/,a)') 'allocation failure in opengsb (gsb.kinds)'
            deallocate(gsb%vnames,stat=i)
            return
        endif
        if(ubound(vtypes,1) == 1)then
            gsb%kinds = vtypes(1)                   !Set all formats to same
        else
            gsb%kinds = vtypes                      !Variable formats (1=i2,2=r4,3=r8)
        endif

        gsb%ndim = ubound(dims,1)                   !Number of dimensions (not counting realizations)

        gsb%dims(1:gsb%ndim) = dims                 !Dimension lengths

        if(present(nsim))then
            if(nsim < 1)then
                write(*,'(/,a)') 'opengsb requires nsim > 0'
                return
            endif
            gsb%nreal = nsim                        !Number of realizations
        else
            gsb%nreal = 1
        endif

        if(present(trimkey))then                    !Use trimming limits for keyout?
            gsb%ikeyout = trimkey
        else
            gsb%ikeyout = 0
        endif

        if(present(varkey))then                     !If keyout, use variable trimming per real
            gsb%ikeyout_real = varkey
        else
            gsb%ikeyout_real = 0
        endif

        if(present(tmin))then                       !Lower trimming limit
            gsb%tmin = tmin
        else
            gsb%tmin = -huge(gsb%tmin)
        endif

        if(present(tmax))then                       !Upper trimming limit
            gsb%tmax = tmax
        else
            gsb%tmax = -huge(gsb%tmax)
        endif

        if(present(tvar))then                       !Variable for trimming
            gsb%tvar = tvar
        else
            gsb%tvar = 1
        endif

        i = write_bin_header ( fid, fname, gsb, gsb%gsbv )
        if(i /= 0)then
            fid = -1
        endif
        gsb%fid = fid

    end function opengsb

    !Write data out to a GSB file.
    ! thedata is passed in as the highest type accepted by this module (R8)
    ! and casted to the required type for output based on gsb%kinds(:).
    integer function gsbwrite ( gsb, thedata, ireal ) result ( err )
        type(gsb_info), target, intent (inout) :: gsb
        real*8, intent (in) :: thedata(:,:) !nvar * n
        integer, intent (in) :: ireal
        integer :: i, j, k, n
        type(gsb_var), pointer :: pgsbv

        err = -1

        !Map data for output
        n = ubound(thedata,2)
        if(n /= gsb%nrec)then
            i = error_report ( 4, gsb )
            return
        endif

        do i = 1, gsb%nvar
            pgsbv => gsb%gsbv(i)
            select case ( gsb%kinds(i) )
            case ( INT2B )
                if(ubound(pgsbv%data_i2,1) < n)then !Memory check
                    deallocate(pgsbv%data_i2,stat=k)
                    allocate(pgsbv%data_i2(n),stat=k)
                    if(k /= 0)then
                        err = error_report ( 1, gsb )
                        return
                    endif
                endif
                do j = 1, n
                    pgsbv%data_i2(j) = nint(thedata(i,j))
                enddo
            case ( REAL4B )
                if(ubound(pgsbv%data_r4,1) < n)then !Memory check
                    deallocate(pgsbv%data_r4,stat=k)
                    allocate(pgsbv%data_r4(n),stat=k)
                    if(k /= 0)then
                        err = error_report ( 1, gsb )
                        return
                    endif
                endif
                do j = 1, n
                    pgsbv%data_r4(j) = real(thedata(i,j),4)
                enddo
            case ( REAL8B )
                if(ubound(pgsbv%data_r8,1) < n)then !Memory check
                    deallocate(pgsbv%data_r8,stat=k)
                    allocate(pgsbv%data_r8(n),stat=k)
                    if(k /= 0)then
                        err = error_report ( 1, gsb )
                        return
                    endif
                endif
                do j = 1, n
                    pgsbv%data_r8(j) = thedata(i,j)
                enddo
            end select
        enddo

        err = write_bin_data ( gsb%fid, gsb, gsb%gsbv, ireal )

    end function gsbwrite

    !Read data from a GSB file.
    ! Data is casted to R8 on output based on types specified in gsb%kinds(:)
    integer function gsbread ( gsb, thedata, ireal, ivar ) result ( err )
        type(gsb_info), target, intent (inout) :: gsb
        real*8, intent (out) :: thedata(:,:) !nvar * n
        integer, intent (in) :: ireal, ivar(:) !Realization number, variable number(s)
        integer :: iv, i, j, k, n
        type(gsb_var), pointer :: pgsbv

        err = -1

        !Checks
        n = ubound(thedata,2)
        if(n /= gsb%nrec)then
            i = error_report ( 4, gsb )
            return
        endif
        if(ubound(thedata,1) /= ubound(ivar,1))then
            write(*,'(a)') 'more columns were requested than are available in thedata array'
            return
        endif

        !Load the data - note freal is set to 1 because anything else makes no sense
        ! TODO - the original logic behind "freal" needs to be clarified.
        err = read_bin_data_vcols ( gsb%fid, gsb, gsb%gsbv, ivar, ireal, freal = 1 )
        if(err /= 0) return

        !Map data for output
        do i = 1, ubound(ivar,1)
            iv = ivar(i)
            pgsbv => gsb%gsbv(iv)
            select case ( gsb%kinds(iv) )
            case ( INT2B )
                do j = 1, n
                    thedata(i,j) = real(pgsbv%data_i2(j),8)
                enddo
            case ( REAL4B )
                do j = 1, n
                    thedata(i,j) = real(pgsbv%data_r4(j),8)
                enddo
            case ( REAL8B )
                do j = 1, n
                    thedata(i,j) = pgsbv%data_r8(j)
                enddo
            end select
        enddo
        err = 0
    end function gsbread

    integer function write_bin_header( fid, fname, gsb, gsbv ) result(errcode)
        ! Opens the binary output file and writes the header information contained in the
        ! gsb object. Also allocates the gsbv data arrays based on the kinds in gsb%
        character*(*), intent(in) :: fname
        type(gsb_info), intent(inout) :: gsb
        type(gsb_var), allocatable, intent(inout) :: gsbv(:)
        integer, intent(out) :: fid
        integer*4 :: j, i

        errcode = 0
        ! Users may forget to deallocate arrays from entry program - insure everything
        ! is ready for allocation
        ! This deallocates required info like vnames and kinds... ...commenting out - RMB
        !call killgsb ( gsb )

        ! Open the output file with stream access
        fid = file_open( fname, .true., 'stream' )
        if( fid < 0 )then
            errcode = error_report ( 2, gsb )
            return
        endif
        j = len_trim(gsb%header)
        write(fid,err=98) j
        write(fid,err=98) gsb%header(1:j) ! Header info (char)
        write(fid,err=98) gsb%nvar ! Number of variables
        write(fid,err=98) gsb%ikeyout ! Keyout array?
        write(fid,err=98) gsb%ndim
        write(fid,err=98) gsb%dims(1:gsb%ndim)
        write(fid,err=98) gsb%nreal
        if( gsb%ikeyout > 0 .and. gsb%nreal > 1 ) write(fid,err=98) gsb%ikeyout_real

        do i = 1, gsb%nvar ! Write information for each variable
            j = len_trim(gsb%vnames(i))
            write(fid,err=98) j
            write(fid,err=98) gsb%vnames(i)(1:j) ! Variable name (char)
            write(fid,err=98) gsb%kinds(i) ! Type code (1=i2,2=r4,3=r8)
        enddo

        gsb%nrec = 1
        do i = 1, gsb%ndim
            gsb%nrec = gsb%nrec*gsb%dims(i)
        enddo

        ! Allocate the required arrays
        errcode = 0
        ! For coding convenience, a previously used GSB (e.g., input)
        ! may be converted to output, meaning that these arrays
        ! may be allocated. Check and deallocate in case.
        if( gsb%ikeyout > 0)then
            if(allocated(gsb%mask)) deallocate(gsb%mask)
            if(allocated(gsb%keyout)) deallocate(gsb%keyout)
            if(allocated(gsb%ko_index)) deallocate(gsb%ko_index)
            allocate( gsb%mask(gsb%nrec), gsb%keyout(gsb%nrec), &
                      gsb%ko_index(gsb%nrec), stat = errcode )
        endif
        allocate( gsbv(gsb%nvar), stat = errcode )
        do i = 1,gsb%nvar
            if(allocated(gsbv(i)%rle)) deallocate(gsbv(i)%rle)
            if(allocated(gsbv(i)%rle_index)) deallocate(gsbv(i)%rle_index)
            if( gsb%kinds(i) == 1 )then
                if(allocated(gsbv(i)%data_i2)) deallocate(gsbv(i)%data_i2)
                allocate(gsbv(i)%data_i2(gsb%nrec),gsbv(i)%rle(gsb%nrec), &
                        gsbv(i)%rle_index(gsb%nrec), stat = errcode )
            elseif( gsb%kinds(i) == 2 )then
                if(allocated(gsbv(i)%data_r4)) deallocate(gsbv(i)%data_r4)
                allocate( gsbv(i)%data_r4(gsb%nrec), stat = errcode )
            else
                if(allocated(gsbv(i)%data_r8)) deallocate(gsbv(i)%data_r8)
                allocate( gsbv(i)%data_r8(gsb%nrec), stat = errcode )
            endif
        enddo
        if( errcode /= 0 )then
            errcode = error_report ( 1, gsb )
            return
        endif
        return
98      errcode = error_report ( 16, gsb )
    end function write_bin_header

    integer function write_bin_data ( fid, gsb, gsbv, ireal ) result(errcode)
        ! Write out the binary data. Must be called in an iterative manner for each
        ! realization (if writing multiple realizations).
        type(gsb_info), intent(inout) :: gsb
        type(gsb_var), intent(inout) :: gsbv(:)
        integer, intent(in) :: fid
        integer*4, intent(in) :: ireal
        integer*4 i, bytes_real, n
        integer*4, allocatable, dimension(:) :: bytes_var
        integer*2, allocatable, dimension(:) :: itemp
        real*4, allocatable, dimension(:) :: r4temp
        real*8, allocatable, dimension(:) :: r8temp
        integer*8 :: j
        integer :: test
        real :: start, finish

        errcode = 0
        if( size(gsbv) /= gsb%nvar )then
            errcode = error_report ( 7, gsb )
            return
        endif
        if( ireal > gsb%nreal )then
            errcode = error_report ( 8, gsb )
            return
        endif
        if( gsb%ikeyout > 0 .and. gsb%tvar == 0 )then
            !TODO - this should actually kill the program, force users to specify...
            write(*,'(/,/,a,/,a,/,/)') 'WARNING: tvar, tmin and tmax not set for gsb object', &
                                       'SETTING DEFAULTS: tvar = 1, tmin = -8, tmax = 1.0e21'
            gsb%tvar = 1
            gsb%tmin = -8 !TODO - why -8, should be -huge()
            gsb%tmax = 1.0e21
        endif

        !
        ! Prepare to write the data
        !

        bytes_real = 0 ! Tracks the size of this realization

        allocate( bytes_var(gsb%nvar) )
        ! If required, calculate the keyout array
        if( gsb%ikeyout > 0)then
            if( ireal == 1 .or. gsb%ikeyout_real  > 0 )then
                if( gsb%kinds(gsb%tvar) == 1 )then
                    gsb%mask = logical( gsbv(gsb%tvar)%data_i2 >= gsb%tmin .and. &
                                        gsbv(gsb%tvar)%data_i2 <= gsb%tmax )
                elseif( gsb%kinds(gsb%tvar) == 2 )then
                    gsb%mask = logical( gsbv(gsb%tvar)%data_r4 >= gsb%tmin .and. &
                                        gsbv(gsb%tvar)%data_r4 <= gsb%tmax )
                else
                    gsb%mask = logical( gsbv(gsb%tvar)%data_r8 >= gsb%tmin .and. &
                                        gsbv(gsb%tvar)%data_r8 <= gsb%tmax )
                endif
                ! Apply the mask to generate a keyout array
                ! and values with null values removed
                call calc_keyout( gsb%mask, gsb%nkeyout, gsb%keyout, gsb%nko_index, gsb%ko_index )
                ! Position = (nkeyout + keyout)*int4
                bytes_real = bytes_real + (1 + gsb%nkeyout)*4
            endif
        endif

        do i = 1,gsb%nvar
            if( gsb%kinds(i) > 1 )then
                ! Real values - no run length encoding (though
                ! potentially a keyout array)
                if( gsb%ikeyout > 0)then
                    n = gsb%nko_index
                else
                    n = gsb%nrec
                endif
                bytes_var(i) = n * (2**gsb%kinds(i))
            else
                if( gsb%ikeyout  > 0)then
                    call calc_rle_keyout( gsbv(i)%data_i2, gsb%keyout, &
                         gsb%nko_index, gsb%ko_index, n,gsbv(i)%rle, gsbv(i)%rle_index )

                else
                    call calc_rle( gsbv(i)%data_i2, n, gsbv(i)%rle, gsbv(i)%rle_index )
                endif
                gsbv(i)%nrle = n ! Store the length of the rle sequence
                ! Size = (nrle + rle + rle_vals)*int4
                bytes_var(i) = (1 + n)*4 + n*2
            endif
            bytes_real = bytes_real + bytes_var(i) + 4
        enddo

        write(fid,err=98) bytes_real
        if( gsb%ikeyout > 0)then
            if( ireal == 1 .or. gsb%ikeyout_real > 0)then
                write(fid,err=98) gsb%nkeyout
                write(fid,err=98) gsb%keyout(1:gsb%nkeyout)
            endif
        endif
        do i = 1,gsb%nvar
            write(fid,err=98) bytes_var(i)
            if( gsb%kinds(i) >  1 )then
                if( gsb%ikeyout < 1 )then
                    if( gsb%kinds(i) ==  2 )then
                        write(fid,err=98) gsbv(i)%data_r4(1:gsb%nrec)
                    else
                        write(fid,err=98) gsbv(i)%data_r8(1:gsb%nrec)
                    endif
                else
                    if( gsb%kinds(i) ==  2 )then
                        allocate( r4temp(gsb%nko_index), stat = test   )
                        if( test /= 0 ) go to 96
                        do j = 1, gsb%nko_index
                            r4temp(j) = gsbv(i)%data_r4( gsb%ko_index(j) )
                        enddo
                        write(fid,err=98) r4temp
                        deallocate(r4temp)
                    else
                        allocate( r8temp(gsb%nko_index), stat = test )
                        if( test /= 0 ) go to 96
                        do j = 1, gsb%nko_index
                            r8temp(j) = gsbv(i)%data_r8( gsb%ko_index(j) )
                        enddo
                        write(fid,err=98) r8temp
                        deallocate(r8temp)
                    endif
                endif
            else
                allocate( itemp(gsbv(i)%nrle), stat = test )
                if( test /= 0 ) go to 96
                do j = 1, gsbv(i)%nrle
                    itemp(j) = gsbv(i)%data_i2( gsbv(i)%rle_index(j) )
                enddo
                write(fid,err=98) gsbv(i)%nrle
                write(fid,err=98) gsbv(i)%rle( 1:gsbv(i)%nrle )
                write(fid,err=98) itemp
                deallocate(itemp)
            endif
        enddo

        return
96      errcode = error_report ( 1, gsb )
        return
98      errcode = error_report ( 16, gsb )
    end function write_bin_data

!-------------------------------------------------------------------------------------------
! Routines for reading a binary file
!-------------------------------------------------------------------------------------------

    integer function read_bin_header( fid, fname, gsb, gsbv ) result(errcode)
        ! Read in the binary header information. Populate the gsb object and allocate
        ! the gsbv data arrays
        character*(*), intent(in)  :: fname
        integer, intent(out) :: fid
        type(gsb_info), intent(inout)  :: gsb
        type(gsb_var), allocatable, intent(inout)  :: gsbv(:)
        integer*4 len, i, j

        errcode = 0
        ! Users may forget to deallocate arrays from entry program - insure everything
        ! is ready for allocation
        call killgsb( gsb ) !, gsbv )
        ! Open the input data file with stream access
        fid = file_open( fname, .false., 'stream' )
        if( fid < 0 )then
            errcode = error_report ( 3, gsb )
            return
        endif
        read(fid,err=97,end=97) len
        read(fid,err=97,end=97) gsb%header(1:len) ! Header info (char)
        read(fid,err=97,end=97) gsb%nvar ! Number of variables
        read(fid,err=97,end=97) gsb%ikeyout ! Keyout array?
        read(fid,err=97,end=97) gsb%ndim
        read(fid,err=97,end=97) gsb%dims(1:gsb%ndim) ! Dimensions
        read(fid,err=97,end=97) gsb%nreal
        ! Initialize in case it doesn't get read...
        gsb%ikeyout_real = 0
        if( gsb%nreal > 1 )then
            if( gsb%ikeyout > 0) read(fid,err=97,end=97) gsb%ikeyout_real
        endif
        allocate(gsb%vnames(gsb%nvar),gsb%kinds(gsb%nvar),stat=i)
        if(i /= 0)then
            errcode = error_report ( 1, gsb )
            return
        endif
        do i = 1,gsb%nvar
            gsb%vnames(i) = ' '
            read(fid,err=97,end=97) len
            read(fid,err=97,end=97) gsb%vnames(i)(1:len) ! Variable name (char)
            read(fid,err=97,end=97) gsb%kinds(i) ! Type code (1=i2,2=r4,3=r8)
        enddo
        gsb%nrec = 1
        do i = 1,gsb%ndim
            gsb%nrec = gsb%nrec*gsb%dims(i)
        enddo
        ! Allocate the required arrays
        errcode = 0
        if( gsb%ikeyout > 0 )then
            allocate( gsb%mask(gsb%nrec), gsb%keyout(gsb%nrec), &
                      gsb%ko_index(gsb%nrec), stat = errcode )
        endif
        allocate( gsbv(gsb%nvar), stat = errcode )
        do i = 1,gsb%nvar
            if( gsb%kinds(i) == 1 )then
               allocate(gsbv(i)%data_i2(gsb%nrec),gsbv(i)%rle(gsb%nrec),&
                        gsbv(i)%rle_index(gsb%nrec), stat = errcode )
            elseif( gsb%kinds(i) == 2 )then
                allocate( gsbv(i)%data_r4(gsb%nrec), stat = errcode )
            else
                allocate( gsbv(i)%data_r8(gsb%nrec), stat = errcode )
            endif
        enddo
        if( errcode /= 0 )then
            errcode = error_report ( 1, gsb )
            return
        endif
        return
97      errcode = error_report ( 15, gsb )
    end function read_bin_header


    integer function read_bin_data_vcols ( fid, gsb, gsbv, vcols, ireal, freal ) result(errcode)
        ! Read in binary data. The required variables are specified by vcols
        integer, intent(in) :: fid, ireal, freal, vcols(:)
        type(gsb_info), intent(inout) :: gsb
        type(gsb_var), intent(inout) :: gsbv(:)
        integer*4 :: bytes
        integer :: icol, i, nvari

        errcode = 0
        if( size(gsbv) /= gsb%nvar )then
            errcode = error_report ( 7, gsb )
            return
        endif
        if( ireal > gsb%nreal )then
            errcode = error_report ( 8, gsb )
            return
        endif
        if( freal > ireal )then
            errcode = error_report ( 9, gsb )
            return
        endif
        if( freal < 1 )then
            errcode = error_report ( 10, gsb )
            return
        endif
        if( maxval(vcols) > gsb%nvar )then
            errcode = error_report ( 6, gsb )
            return
        endif
        nvari = size(vcols)
        read(fid,err=97,end=97) bytes
        if( ireal < freal)then
            if( ireal == 1 .and. gsb%ikeyout > 0 .and. gsb%ikeyout_real < 1 )then
                ! Have to read and store the keyout
                errcode = read_bin_keyout( fid, gsb  )
                bytes = bytes - (gsb%nkeyout + 1)*4
            endif
            ! Skip the remainder of this realization
            errcode = skip_pos( fid, bytes )
            if( errcode /= 0 )then
                return
            endif
        elseif( ireal == 1 .and. gsb%ikeyout > 0)then
            errcode = read_bin_keyout( fid, gsb )
        elseif( gsb%ikeyout_real > 0 )then
            errcode = read_bin_keyout( fid, gsb )
        endif

        icol = 1
        do i = 1,gsb%nvar
            read(fid, err = 97,end=97) bytes
            if( icol > nvari )then
                errcode = skip_pos( fid, bytes )
                if( errcode /= 0 )then
                    return
                endif
                cycle
            elseif( i /= vcols(icol) )then
                errcode = skip_pos( fid, bytes )
                if( errcode /= 0 )then
                    return
                endif
                cycle
            endif

            ! Read in this variable
            errcode = read_bin_var( fid, gsb, gsbv, vcols(icol)  )
            icol = icol + 1
        enddo
        return
97      errcode = error_report ( 15, gsb )
    end function read_bin_data_vcols

    integer function read_bin_data_vcol( fid, gsb, gsbv, vcol, ireal, freal ) result(errcode)
        ! Read in binary data. The required variable is specified by vcol
        integer, intent(in) :: fid, ireal, freal, vcol
        type(gsb_info), intent(inout) :: gsb
        type(gsb_var), intent(inout) :: gsbv(:)
        integer*4 :: bytes
        integer :: icol, i, nvari

        errcode = 0
        if( size(gsbv) /= gsb%nvar )then
            errcode = error_report ( 7, gsb )
            return
        endif
        if( ireal > gsb%nreal )then
            errcode = error_report ( 8, gsb )
            return
        endif
        if( freal > ireal )then
            errcode = error_report ( 9, gsb )
            return
        endif
        if( freal < 1 )then
            errcode = error_report ( 10, gsb )
            return
        endif
        if( vcol > gsb%nvar )then
            errcode = error_report ( 6, gsb )
            return
        endif
        read(fid,err=97,end=97) bytes
        if( ireal < freal)then
            if( ireal == 1 .and. gsb%ikeyout > 0 .and. gsb%ikeyout_real < 1 )then
                ! Have to read and store the keyout
                errcode = read_bin_keyout( fid, gsb  )
                bytes = bytes - (gsb%nkeyout + 1)*4
            endif
            ! Skip the remainder of this realization
            errcode = skip_pos( fid, bytes )
            if( errcode /= 0 )then
                return
            endif
        elseif( ireal == 1 .and. gsb%ikeyout > 0)then
            errcode = read_bin_keyout( fid, gsb )
        elseif( gsb%ikeyout_real > 0 )then
            errcode = read_bin_keyout( fid, gsb )
        endif
        do i = 1,gsb%nvar
            read(fid, err = 97,end=97) bytes
            if( i /= vcol )then

                errcode = skip_pos( fid, bytes )
                if( errcode /= 0 )then
                    errcode = 0
                    return
                endif
                cycle
            endif
            ! Read in this variable
            errcode = read_bin_var( fid, gsb, gsbv, vcol  )
        enddo
        return
97      errcode = error_report ( 15, gsb )
    end function read_bin_data_vcol


!-------------------------------------------------------------------------------------------
! Subroutines for converting a gsbv array to a vector or matrix array
!-------------------------------------------------------------------------------------------

    integer function gsbv2vector_i2( gsb, gsbv, vcol, data_  ) result(errcode)
        ! Outputs integer*2 data_ (one dimension)
        type(gsb_info), intent(inout) :: gsb
        type(gsb_var), intent(in) :: gsbv(:)
        integer, intent(in) :: vcol
        integer*2, intent(out) :: data_(:)
        integer i, nvari
        errcode = 0
        if( size(data_,1) /= gsb%nrec )then
            errcode = error_report ( 4, gsb )
            return
        elseif( vcol > gsb%nvar )then
            errcode = error_report ( 6, gsb )
            return
        endif
        if( gsb%kinds(vcol) == 1 )then
            data_ = gsbv(vcol)%data_i2
        elseif( gsb%kinds(vcol) == 2 )then
            data_ = nint(gsbv(vcol)%data_r4)
        else
            data_ = nint(gsbv(vcol)%data_r8)
        endif
    end function gsbv2vector_i2

    integer function gsbv2vector_r4( gsb, gsbv, vcol, data_  ) result(errcode)
        ! Outputs real*4 data_ (one dimension)
        type(gsb_info), intent(inout) :: gsb
        type(gsb_var), intent(in) :: gsbv(:)
        integer, intent(in) :: vcol
        real*4, intent(out) :: data_(:)
        integer i, nvari
        errcode = 0
        if( size(data_,1) /= gsb%nrec )then
            errcode = error_report ( 4, gsb )
            return
        elseif( vcol > gsb%nvar )then
            errcode = error_report ( 6, gsb )
            return
        endif
        if( gsb%kinds(vcol) == 1 )then
            data_ = real(gsbv(vcol)%data_i2)
        elseif( gsb%kinds(vcol) == 2 )then
            data_ = gsbv(vcol)%data_r4
        else
            data_ = gsbv(vcol)%data_r8
        endif
    end function gsbv2vector_r4

    integer function gsbv2vector_r8( gsb, gsbv, vcol, data_  ) result(errcode)
        ! Outputs real*8 data_ (one dimension)
        type(gsb_info), intent(inout) :: gsb
        type(gsb_var), intent(in) :: gsbv(:)
        integer, intent(in) :: vcol
        real*8, intent(out) :: data_(:)
        integer i, nvari
        errcode = 0
        if( size(data_,1) /= gsb%nrec )then
            errcode = error_report ( 4, gsb )
            return
        elseif( vcol > gsb%nvar )then
            errcode = error_report ( 6, gsb )
            return
        endif
        if( gsb%kinds(vcol) == 1 )then
            data_ = dble(gsbv(vcol)%data_i2)
        elseif( gsb%kinds(vcol) == 2 )then
            data_ = dble(gsbv(vcol)%data_r4)
        else
            data_ = gsbv(vcol)%data_r8
        endif
    end function gsbv2vector_r8

    integer function gsbv2matrix_i2( gsb, gsbv, vcols, data_  ) result(errcode)
        ! Outputs integer*2 data_ (two dimensions)
        type(gsb_info), intent(inout) :: gsb
        type(gsb_var), intent(in) :: gsbv(:)
        integer, intent(in) :: vcols(:)
        integer*2, intent(out) :: data_(:,:)
        integer i, nvari
        errcode = 0
        nvari = size(vcols)
        if( size(data_,2) /= gsb%nrec )then
            errcode = error_report ( 4, gsb )
            return
        elseif( size(data_,1) /= nvari )then
            errcode = error_report ( 5, gsb )
            return
        elseif( maxval(vcols) > gsb%nvar )then
            errcode = error_report ( 6, gsb )
            return
        endif
        do i = 1, nvari
            if( gsb%kinds(vcols(i)) == 1 )then
                data_(i,:) = gsbv(vcols(i))%data_i2
            elseif( gsb%kinds(vcols(i)) == 2 )then
                data_(i,:) = nint(gsbv(vcols(i))%data_r4)
            else
                data_(i,:) = nint(gsbv(vcols(i))%data_r8)
            endif
        enddo
    end function gsbv2matrix_i2

    integer function gsbv2matrix_r4( gsb, gsbv, vcols, data_  ) result(errcode)
        ! Outputs real*4 data_ (two dimensions)
        type(gsb_info), intent(inout) :: gsb
        type(gsb_var), intent(in) :: gsbv(:)
        integer, intent(in) :: vcols(:)
        real*4, intent(out) :: data_(:,:)
        integer i, nvari
        errcode = 0
        nvari = size(vcols)
        if( size(data_,2) /= gsb%nrec )then
            errcode = error_report ( 4, gsb )
            return
        elseif( size(data_,1) /= nvari )then
            errcode = error_report ( 5, gsb )
            return
        elseif( maxval(vcols) > gsb%nvar )then
            errcode = error_report ( 6, gsb )
        endif
        do i = 1, nvari
            if( gsb%kinds(vcols(i)) == 1 )then
                data_(i,:) = real(gsbv(vcols(i))%data_i2)
            elseif( gsb%kinds(vcols(i)) == 2 )then
                data_(i,:) = gsbv(vcols(i))%data_r4
            else
                data_(i,:) = gsbv(vcols(i))%data_r8
            endif
        enddo
    end function gsbv2matrix_r4

    integer function gsbv2matrix_r8( gsb, gsbv, vcols, data_  ) result(errcode)
        ! Outputs real*8 data_ (two dimensions)
        type(gsb_info), intent(inout) :: gsb
        type(gsb_var), intent(in) :: gsbv(:)
        integer, intent(in) :: vcols(:)
        real*8, intent(out) :: data_(:,:)
        integer i, nvari
        errcode = 0
        nvari = size(vcols)
        if( size(data_,2) /= gsb%nrec )then
            errcode = error_report ( 4, gsb )
            return
        elseif( size(data_,1) /= nvari )then
            errcode = error_report ( 5, gsb )
            return
        elseif( maxval(vcols) > gsb%nvar )then
            errcode = error_report ( 6, gsb )
        endif
        do i = 1, nvari
            if( gsb%kinds(vcols(i)) == 1 )then
                data_(i,:) = dble(gsbv(vcols(i))%data_i2)
            elseif( gsb%kinds(vcols(i)) == 2 )then
                data_(i,:) = dble(gsbv(vcols(i))%data_r4)
            else
                data_(i,:) = gsbv(vcols(i))%data_r8
            endif
        enddo
    end function gsbv2matrix_r8

!-------------------------------------------------------------------------------------------
! Subroutines for converting a vector or matrix array into potentially
! different gsbv kinds (preparing to write out to binary).
! R. Barnett - 2015
!-------------------------------------------------------------------------------------------

    integer function vector2gsbv_i2( gsb, gsbv, vcol, data_  ) result(errcode)
        ! Input integer*2 data_ (one dimension)
        type(gsb_info), intent(inout) :: gsb
        integer, intent(in) :: vcol
        integer*2, intent(in) :: data_(:)
        type(gsb_var), intent(inout) :: gsbv(:)
        integer i
        errcode = 0
        if( size(data_) /= gsb%nrec )then
            errcode = error_report ( 4, gsb )
            return
        elseif( vcol > gsb%nvar )then
            errcode = error_report ( 6, gsb )
            return
        endif
        if( gsb%kinds(vcol) == 1 )then
            gsbv(vcol)%data_i2 = data_
        elseif( gsb%kinds(vcol) == 2 )then
            gsbv(vcol)%data_r4 = real(data_)
        else
            gsbv(vcol)%data_r8 = dble(data_)
        endif
    end function vector2gsbv_i2

    integer function vector2gsbv_r4( gsb, gsbv, vcol, data_  ) result(errcode)
        ! Input real*4 data_ (one dimension)
        type(gsb_info), intent(inout) :: gsb
        integer, intent(in) :: vcol
        real*4, intent(in) :: data_(:)
        type(gsb_var), intent(inout) :: gsbv(:)
        integer i
        errcode = 0
        if( size(data_) /= gsb%nrec )then
            errcode = error_report ( 4, gsb )
            return
        elseif( vcol > gsb%nvar )then
            errcode = error_report ( 6, gsb )
            return
        endif
        if( gsb%kinds(vcol) == 1 )then
            gsbv(vcol)%data_i2 = nint(data_)
        elseif( gsb%kinds(vcol) == 2 )then
            gsbv(vcol)%data_r4 = data_
        else
            gsbv(vcol)%data_r8 = data_
        endif
    end function vector2gsbv_r4

    integer function vector2gsbv_r8( gsb, gsbv, vcol, data_  ) result(errcode)
        ! Input real*8 data_ (one dimension)
        type(gsb_info), intent(inout) :: gsb
        integer, intent(in) :: vcol
        real*8, intent(in) :: data_(:)
        type(gsb_var), intent(inout) :: gsbv(:)
        integer i
        errcode = 0
        if( size(data_) /= gsb%nrec )then
            errcode = error_report ( 4, gsb )
            return
        elseif( vcol > gsb%nvar )then
            errcode = error_report ( 6, gsb )
            return
        endif
        if( gsb%kinds(vcol) == 1 )then
            gsbv(vcol)%data_i2 = nint(data_)
        elseif( gsb%kinds(vcol) == 2 )then
            gsbv(vcol)%data_r4 = data_
        else
            gsbv(vcol)%data_r8 = data_
        endif
    end function vector2gsbv_r8

    integer function matrix2gsbv_i2( gsb, gsbv, vcols, data_  ) result(errcode)
        ! Input integer*2 data_ (two dimensions)
        type(gsb_info), intent(inout) :: gsb
        integer, intent(in) :: vcols(:)
        integer*2, intent(in) :: data_(:,:)
        type(gsb_var), intent(inout) :: gsbv(:)
        integer i, nvari
        nvari = size(vcols)
        errcode = 0
        if( size(data_,2) /= gsb%nrec )then
            errcode = error_report ( 4, gsb )
            return
        elseif( size(data_,1) /= nvari )then
            errcode = error_report ( 5, gsb )
            return
        elseif( maxval(vcols) > gsb%nvar )then
            errcode = error_report ( 6, gsb )
            return
        endif
        do i = 1, nvari
            if( gsb%kinds(vcols(i)) == 1 )then
                gsbv(vcols(i))%data_i2 = data_(i,:)
            elseif( gsb%kinds(vcols(i)) == 2 )then
                gsbv(vcols(i))%data_r4 = real(data_(i,:))
            else
                gsbv(vcols(i))%data_r8 = dble(data_(i,:))
            endif
        enddo
    end function matrix2gsbv_i2

    integer function matrix2gsbv_r4( gsb, gsbv, vcols, data_  ) result(errcode)
        ! Input real*4 data_ (two dimensions)
        type(gsb_info), intent(inout) :: gsb
        integer, intent(in) :: vcols(:)
        real*4, intent(in) :: data_(:,:)
        type(gsb_var), intent(inout) :: gsbv(:)
        integer i, nvari
        errcode = 0
        nvari = size(vcols)
        if( size(data_,2) /= gsb%nrec )then
            errcode = error_report ( 4, gsb )
            return
        elseif( size(data_,1) /= nvari )then
            errcode = error_report ( 5, gsb )
            return
        elseif( maxval(vcols) > gsb%nvar )then
            errcode = error_report ( 6, gsb )
            return
        endif
        do i = 1, nvari
            if( gsb%kinds(vcols(i)) == 1 )then
                gsbv(vcols(i))%data_i2 = nint(data_(i,:))
            elseif( gsb%kinds(vcols(i)) == 2 )then
                gsbv(vcols(i))%data_r4 = data_(i,:)
            else
                gsbv(vcols(i))%data_r8 = data_(i,:)
            endif
        enddo
    end function matrix2gsbv_r4

    integer function matrix2gsbv_r8( gsb, gsbv, vcols, data_  ) result(errcode)
        ! Input real*8 data_ (two dimensions)
        type(gsb_info), intent(inout) :: gsb
        integer, intent(in) :: vcols(:)
        real*8, intent(in) :: data_(:,:)
        type(gsb_var), intent(inout) :: gsbv(:)
        integer i, nvari
        errcode = 0
        nvari = size(vcols)
        if( size(data_,2) /= gsb%nrec )then
            errcode = error_report ( 4, gsb )
            return
        elseif( size(data_,1) /= nvari )then
            errcode = error_report ( 5, gsb )
            return
        elseif( maxval(vcols) > gsb%nvar )then
            errcode = error_report ( 6, gsb )
            return
        endif
        do i = 1, nvari
            if( gsb%kinds(vcols(i)) == 1 )then
                gsbv(vcols(i))%data_i2 = nint(data_(i,:))
            elseif( gsb%kinds(vcols(i)) == 2 )then
                gsbv(vcols(i))%data_r4 = data_(i,:)
            else
                gsbv(vcols(i))%data_r8 = data_(i,:)
            endif
        enddo
    end function matrix2gsbv_r8

!----------------------------------------------------------------------------
! Routines for compressing and uncompressing a runlength encoded keyout array
!----------------------------------------------------------------------------

    subroutine calc_keyout( mask, nkeyout, keyout, nko_index, ko_index)
        ! Calculate a masking keyout vector based on whether
        ! the associated data value is trimmed (0) or not (1)
        logical, intent(in) :: mask(:)
        integer*4, intent(out) :: keyout(:), ko_index(:)
        integer*4 :: nrec, nkeyout, nko_index, irec, nrun

        ! The first value indicates whether the start of the sequence
        ! begins with a null value (0) or not (1). Null values are
        ! not actually stored/written/read.
        nrec = ubound(mask,1)
        irec = 0
        nkeyout = 1
        if( mask(1) )then
            keyout(nkeyout) = 1
        else
            keyout(nkeyout) = 0
            ! The main do loop must be initialized on a
            ! non-null value so move through this run of
            ! null values
            irec = irec + 1
            nkeyout = nkeyout + 1
            nrun = 1
            if( irec+1  > nrec )then
                keyout(nkeyout) = nrun
                write(*,'(a)') 'Warning - all data are null values!'
            else
                do while( .not. mask(irec+1) )
                    irec = irec + 1
                    nrun = nrun + 1
                    if( irec == nrec ) exit
                enddo
                keyout(nkeyout) = nrun
            endif
        endif
        nko_index = 0
        do while( irec < nrec )
            ! Run of non-null values
            irec = irec + 1
            nkeyout = nkeyout + 1
            nko_index = nko_index + 1
            ko_index(nko_index) = irec
            nrun = 1
            if( irec+1  > nrec )then
                keyout(nkeyout) = nrun
                exit
            else
                do while( mask(irec+1) )
                    irec = irec + 1
                    nrun = nrun + 1
                    nko_index = nko_index + 1
                    ko_index(nko_index) = irec ! Stores the index to the original data file
                    if( irec == nrec ) exit
                enddo
                keyout(nkeyout) = nrun
            endif

            if( irec+1  > nrec ) exit

            ! Run of null values
            irec = irec + 1
            nkeyout = nkeyout + 1
            nrun = 1
            if( irec+1  > nrec )then
                keyout(nkeyout) = nrun
            else
                do while( .not. mask(irec+1) )
                    irec = irec + 1
                    nrun = nrun + 1
                    if( irec == nrec ) exit
                enddo
                keyout(nkeyout) = nrun
            endif
         enddo

    end subroutine calc_keyout

    subroutine calc_nko_index( nkeyout, keyout, nko_index )
        ! Calculate a masking keyout vector based on whether
        ! the associated data value is trimmed (0) or not (1)
        integer*4, intent(in) :: nkeyout
        integer*4, intent(in) :: keyout(:)
        integer*4, intent(out) :: nko_index
        integer*4 :: ikeyout
        ! The first value indicates whether the start of the sequence
        ! begins with a null value (0) or not (1). Null values are
        ! not actually stored/written/read.

        ikeyout = 1
        if( keyout(ikeyout) == 0 )then
            ! First values are null, so initialize
            ! as such to facilitate subsequent loop
            ikeyout = ikeyout + 1
        endif
        nko_index = 0
        do while( ikeyout < nkeyout )
            ikeyout = ikeyout + 1
            nko_index = nko_index + keyout(ikeyout)
            if( ikeyout == nkeyout ) exit
            ikeyout = ikeyout + 1
        enddo
    end subroutine calc_nko_index

    subroutine calc_keyout_inv_i2( valsin, nkeyout, keyout, valsout )
        ! Calculate a masking keyout vector based on whether
        ! the associated data value is trimmed (0) or not (1)
        integer*2, intent(in) :: valsin(:)
        integer*4, intent(in) :: nkeyout
        integer*4, intent(in) :: keyout(:)
        integer*2, intent(out) :: valsout(:)
        integer*4 :: nrec, irec, nko_index, ikeyout, iko_index

        ! The first value indicates whether the start of the sequence
        ! begins with a null value (0) or not (1). Null values are
        ! not actually stored/written/read.
        ikeyout = 1
        irec = 0
        if( keyout(ikeyout) == 0 )then
            ! First values are null, so initialize
            ! as such to facilitate subsequent loop
            ikeyout = ikeyout + 1
            nrec = irec + keyout(ikeyout)
            valsout(irec+1:nrec) = int(ASCII_NULL_VALUE)
            irec = nrec
        endif
        iko_index = 0
        do while( ikeyout < nkeyout )
            ikeyout = ikeyout + 1
            nko_index = iko_index + keyout(ikeyout)
            nrec = irec + keyout(ikeyout)
            valsout(irec+1:nrec) = valsin(iko_index+1:nko_index)
            iko_index = nko_index
            irec = nrec
            if( ikeyout == nkeyout ) exit
            ikeyout = ikeyout + 1
            nrec = irec + keyout(ikeyout)
            valsout(irec+1:nrec) = int(ASCII_NULL_VALUE)
            irec = nrec
        enddo
    end subroutine calc_keyout_inv_i2

    subroutine calc_keyout_inv_r4( valsin, nkeyout, keyout, valsout )
        ! Calculate a masking keyout vector based on whether
        ! the associated data value is trimmed (0) or not (1)
        real*4, intent(in) :: valsin(:)
        integer*4, intent(in) :: nkeyout
        integer*4, intent(in) :: keyout(:)
        real*4, intent(out) :: valsout(:)
        integer*4 :: nrec, irec, nko_index, ikeyout, iko_index

        ! The first value indicates whether the start of the sequence
        ! begins with a null value (0) or not (1). Null values are
        ! not actually stored/written/read.
        ikeyout = 1
        irec = 0
        if( keyout(ikeyout) == 0 )then
            ! First values are null, so initialize
            ! as such to facilitate subsequent loop
            ikeyout = ikeyout + 1
            nrec = irec + keyout(ikeyout)
            valsout(irec+1:nrec) = real(ASCII_NULL_VALUE, 4)
            irec = nrec
        endif
        iko_index = 0
        do while( ikeyout < nkeyout )
            ikeyout = ikeyout + 1
            nko_index = iko_index + keyout(ikeyout)
            nrec = irec + keyout(ikeyout)
            valsout(irec+1:nrec) = valsin(iko_index+1:nko_index)
            iko_index = nko_index
            irec = nrec
            if( ikeyout == nkeyout ) exit
            ikeyout = ikeyout + 1
            nrec = irec + keyout(ikeyout)
            valsout(irec+1:nrec) = real(ASCII_NULL_VALUE, 4)
            irec = nrec
        enddo
    end subroutine calc_keyout_inv_r4

    subroutine calc_keyout_inv_r8( valsin, nkeyout, keyout, valsout )
        ! Calculate a masking keyout vector based on whether
        ! the associated data value is trimmed (0) or not (1)
        real*8, intent(in) :: valsin(:)
        integer*4, intent(in) :: nkeyout
        integer*4, intent(in) :: keyout(:)
        real*8, intent(out) :: valsout(:)
        integer*4 :: nrec, irec, nko_index, ikeyout, iko_index

        ! The first value indicates whether the start of the sequence
        ! begins with a null value (0) or not (1). Null values are
        ! not actually stored/written/read.
        ikeyout = 1
        irec = 0
        if( keyout(ikeyout) == 0 )then
            ! First values are null, so initialize
            ! as such to facilitate subsequent loop
            ikeyout = ikeyout + 1
            nrec = irec + keyout(ikeyout)
            valsout(irec+1:nrec) = ASCII_NULL_VALUE
            irec = nrec
        endif
        iko_index = 0
        do while( ikeyout < nkeyout )
            ikeyout = ikeyout + 1
            nko_index = iko_index + keyout(ikeyout)
            nrec = irec + keyout(ikeyout)
            valsout(irec+1:nrec) = valsin(iko_index+1:nko_index)
            iko_index = nko_index
            irec = nrec
            if( ikeyout == nkeyout ) exit
            ikeyout = ikeyout + 1
            nrec = irec + keyout(ikeyout)
            valsout(irec+1:nrec) = ASCII_NULL_VALUE
            irec = nrec
        enddo
    end subroutine calc_keyout_inv_r8

!---------------------------------------------------------------------
! Routines for compressing and uncompressing a runlength encoded array
!---------------------------------------------------------------------

    subroutine calc_rle_keyout( vals, keyout, nko_index, ko_index,&
                                               &  nrle, rle, rle_index )
        ! Calculate the run length encoding (RLE) of an input
        ! integer array. The RLE may be performed as a second level
        ! of compression beyond the keyout (which masks out null values)
        integer*2, intent(in) :: vals(:)
        integer*4, intent(in) :: keyout(:), ko_index(:)
        integer*4, intent(out) :: rle(:), rle_index(:)
        integer*4 :: iko_index, nko_index, nrle, nrun, ikeyout

        ! Initialize keyout
        if( keyout(1) > 0 )then
            ! Non-null value in first record so initialize here
            ikeyout = 1
        else
            ! Null value in first record so initialize in second
            ikeyout = 2
        endif

        ! Determine the modified run length encoding
        nrle = 0; iko_index = 0
        do while( iko_index < nko_index )
            ! Query the current keyout run of non-null values
            ikeyout = ikeyout + 1
            iko_index = iko_index + 1
            nrle = nrle + 1
            rle_index(nrle) = ko_index(iko_index)
            nrun = 1
            if( iko_index== nko_index )then
                ! Record the run of this integer
                rle(nrle) = nrun
            else
                do while( vals(ko_index(iko_index)) == vals(ko_index(iko_index+1)) )
                    iko_index = iko_index + 1
                    nrun = nrun + 1
                    if( iko_index == nko_index ) exit
                enddo
                rle(nrle) = nrun
                ! Move forward through the non-null values
                ikeyout = ikeyout + 1
            endif
         enddo
    end subroutine calc_rle_keyout

    subroutine calc_rle( vals, nrle, rle, rle_index)
        ! Calculate the run length encoding (RLE) of an input
        ! integer array.
        integer*2, intent(in) :: vals(:)
        integer*4, intent(out) :: rle(:), rle_index(:)
        integer*4 :: irec, nrec, nrle, nrun

        ! Determine the modified run length encoding
        irec = 0 ! Index  of the original records
        nrec = ubound(vals,1)
        nrle = 0 ! Index of the number of of runs
        do while( irec < nrec )
            irec = irec + 1
            nrle = nrle + 1
            ! Record the index of this integer
            rle_index(nrle) = irec
            nrun = 1 ! Length of the run
            if( irec == nrec )then
                ! Record the run of this integer
                rle(nrle) = nrun
            else
                do while( vals(irec) == vals(irec+1) )
                    irec = irec + 1
                    nrun = nrun + 1
                    if( irec == nrec ) exit
                enddo
                rle(nrle) = nrun
            endif
         enddo
    end subroutine calc_rle

    subroutine calc_rle_inv( valsin, nrle, rle, valsout )
        ! Calculate the inverse of run length encoding (RLE),
        ! yielding the uncompressed values.
        integer*2, intent(in) :: valsin(:)
        integer*4, intent(in) :: rle(:)
        integer*2, intent(out) :: valsout(:)
        integer*4 :: irec, nrec, nrle, irle
        irec = 0 ! Index of the uncompressed data array
        do irle = 1,nrle ! Index of the compressed data array
            nrec = irec + rle(irle)
            valsout(irec+1:nrec) = valsin(irle)
            irec = nrec
        enddo
    end subroutine calc_rle_inv

!--------------------------------------------------------------------------------
! Routines for determining the position in a binary file, skipping to a position,
! reading a keyout, reading variables, etc.
!--------------------------------------------------------------------------------


    integer function skip_pos( fid, bytes ) result(errcode)
        integer, intent(in) :: fid
        integer*4, intent(in) :: bytes
        integer*8 :: mypos
        integer*1 :: junk
        errcode = 0
        read(fid, err = 97,end=97) junk
        inquire( fid, POS = mypos)
        mypos = mypos + bytes - 1
        !RMB - end isn't necessarily an error if the end of the file
        ! in the case of the GSB format
        read( fid,  POS = mypos, end = 96, err = 97 )
96      return
97      errcode = 15
    end function skip_pos

    integer function read_bin_keyout( fid, gsb ) result(errcode)
        ! Read in the keyout information
        integer, intent(in) :: fid
        type(gsb_info), intent(inout) :: gsb
        errcode = 0
        read(fid,err=97,end=97) gsb%nkeyout
        read(fid,err=97,end=97) gsb%keyout(1:gsb%nkeyout)
        call calc_nko_index( gsb%nkeyout, gsb%keyout, gsb%nko_index )
        return
97      errcode = error_report ( 15, gsb )
    endfunction read_bin_keyout

    integer function read_bin_var( fid, gsb, gsbv, ivar  ) result(errcode)
        ! Read in a binary variable
        integer, intent(in) :: fid
        integer, intent (in) :: ivar
        type(gsb_info), intent(inout) :: gsb
        type(gsb_var), intent(inout) :: gsbv(:)
        integer*4 n, i
        errcode = 0
        if( gsb%kinds(ivar) >  1 )then
            if( gsb%ikeyout < 1 )then
                if( gsb%kinds(ivar) == 2 )then
                    read(fid,err=97,end=97) gsbv(ivar)%data_r4
                else
                    read(fid,err=97,end=97) gsbv(ivar)%data_r8
                endif
            else
                if( gsb%kinds(ivar) == 2 )then
                    allocate( r4temp11(gsb%nko_index) )
                    read(fid,err=97,end=97) r4temp11 ! Temporary storage of input compressed keyout data
                    call calc_keyout_inv_r4( r4temp11, gsb%nkeyout, gsb%keyout, gsbv(ivar)%data_r4 )
                    deallocate(r4temp11)
                 else
                    allocate( r8temp11(gsb%nko_index) )
                    read(fid,err=97,end=97) r8temp11 ! Temporary storage of input compressed keyout data
                    call calc_keyout_inv_r8( r8temp11, gsb%nkeyout, gsb%keyout, gsbv(ivar)%data_r8 )
                    deallocate(r8temp11)
                 endif
            endif
        else
            read(fid,err=97,end=97) gsbv(ivar)%nrle
            read(fid,err=97,end=97) gsbv(ivar)%rle( 1:gsbv(ivar)%nrle )
            n = sum( gsbv(ivar)%rle( 1:gsbv(ivar)%nrle ) )
            allocate( i2temp11(gsbv(ivar)%nrle), i2temp12(n) )
            ! i2temp11 is temporary storage of input compressed rle data
            ! i2temp12 is temporary storage of output uncompressed rle data
            read(fid,err=97,end=97) i2temp11
            call calc_rle_inv( i2temp11, gsbv(ivar)%nrle, gsbv(ivar)%rle, i2temp12 )
            deallocate(i2temp11)
            if( gsb%ikeyout < 1 )then
                gsbv(ivar)%data_i2 = i2temp12
            else
                call calc_keyout_inv_i2( i2temp12, gsb%nkeyout, gsb%keyout, gsbv(ivar)%data_i2 )
            endif
            deallocate(i2temp12)
        endif
        return
97      errcode = error_report ( 15, gsb )
    end function read_bin_var

!-----------------------------------------------------------------------------
! Other miscelanneous functions and subroutines
!-----------------------------------------------------------------------------

    subroutine killgsb ( gsb )
        ! Subroutine for deallocating gsb and gsbv arrays
        type(gsb_info), target, intent(inout) :: gsb
        type(gsb_var), pointer :: pgsbv(:)
        integer i

        pgsbv => gsb%gsbv

        if( allocated(gsb%mask) ) deallocate(gsb%mask)
        if( allocated(gsb%keyout) ) deallocate(gsb%keyout)
        if( allocated(gsb%ko_index) ) deallocate(gsb%ko_index)
        if( allocated(gsb%gsbv) )then
            do i = 1,gsb%nvar
                if( allocated(pgsbv(i)%data_i2) )then
                    deallocate( pgsbv(i)%data_i2, pgsbv(i)%rle)
                elseif( allocated( pgsbv(i)%data_r4  ) )then
                    deallocate( pgsbv(i)%data_r4  )
                elseif( allocated( pgsbv(i)%data_r8  ) )then
                    deallocate( pgsbv(i)%data_r8 )
                endif
            enddo
        endif
        if( allocated(gsb%vnames) ) deallocate(gsb%vnames)
        if( allocated(gsb%kinds) ) deallocate(gsb%kinds)
        pgsbv => null()
        deallocate(gsb%gsbv,stat=i)
     endsubroutine killgsb

    integer function error_report ( test, gsb ) result(errcode)
        ! Function for reporting common errors to the user.
        ! Only call in even of an error as it deletes all memory associated with gsb
        integer, intent(in) :: test
        type(gsb_info), intent (inout) :: gsb
        call killgsb ( gsb )
        select case (test)
            case(1)
                write(*,'(/,a,/)') 'ERROR: allocation failed for GSB array'
            case(2)
                write(*,'(/,a,/)') 'ERROR: could not open output GSB file'
            case(3)
                write(*,'(/,a,/)') 'ERROR: could not open input GSB file'
            case(4)
                write(*,'(/,a,/,a)') 'ERROR: second dimension of data_ array ',&
                                     '       should match gsb%nrec'
            case(5)
                write(*,'(/,a,/,a)') 'ERROR: first dimension of data_ array ',&
                                     '       should match length of vcols'
            case(6)
                write(*,'(/,a,/)') 'ERROR: value in vcols exceeds gsb%nvar'
            case(7)
                write(*,'(/,a,/)') 'ERROR: length of gsbv must match gsb%nvar'
            case(8)
                write(*,'(/,a,/)') 'ERROR: ireal greater than gsb%nreal'
            case(9)
                write(*,'(/,a,/)') 'ERROR: ireal must be greater than or equal to freal'
            case(10)
                write(*,'(/,a,/)') 'ERROR: freal must be greater than 1'
            case(15)
                write(*,'(/,a,/,a)') 'ERROR: could not read from GSB file'
            case(16)
                write(*,'(/,a,/,a)') 'ERROR: could not write to GSB file'
        end select
        errcode = test
     end function error_report

end module gslib_binary


