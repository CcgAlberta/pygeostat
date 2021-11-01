! public
module fortran_read
! ----------------------------------------------------------------------
! Pygeostat fortan module for reading a file
! ----------------------------------------------------------------------
    implicit none
    real*8, allocatable, dimension(:, :), public :: data_double
    real, allocatable, dimension(:, :), public   :: data_single
    integer, public :: datalength, num_var, real_read

contains
    subroutine read_gslib_grid(filename, gridsize, nr, real_to_read, double)
    ! -------------------------------------------------------------------
    !
    ! Simple fortran subroutine for reading a file for pygeostat import
    !   - USE For gridded Data
    !
    !
    ! author: Ryan Martin & Tyler Acorn
    ! -------------------------------------------------------------------
        character*(*), intent(in)   :: filename
        character*512               :: str
        integer                     :: lin, ir, err, num_var, i
        integer                     :: gridsize, nr, real_to_read
        !real*8, allocatable         :: tmp(:)
        logical                     :: testfl, double

        ! -------------------------------------------------
        !standard gslib file IO to get the length of the data_ file  and the number
        !of variables to use to import using read_gslibfl below
        lin=1
        ! -------------------------------------------------
        ! Open File and read header
        ! -------------------------------------------------
        open( lin, file=filename, status='OLD' )
        read( lin, *, end=10, err=99 )
        read( lin, *, end=10, err=99 ) num_var
        do i = 1, num_var
            read(lin,*,end=10,err=99)
        enddo

        ! -------------------------------------------------
        ! Check Memory
        ! -------------------------------------------------
        if( real_to_read == 0 ) then ! read in all realizations
            datalength = gridsize * nr
        else
            datalength = gridsize ! only read in specified realization
        end if

        if( allocated(data_double) ) deallocate(data_double)
        if( allocated(data_single) ) deallocate(data_single)

        if( double .eqv. .true.) then
            allocate( data_double(num_var, datalength), stat=err )
            if (err .ne. 0 ) stop ' ERROR: memory allocation failure! - data_double(num_var, datalength)'
        else
            allocate( data_single(num_var, datalength), stat=err )
            if (err .ne. 0 ) stop ' ERROR: memory allocation failure! - data_single(num_var, datalength)'
        end if

        ! -------------------------------------------------
        ! Read Data
        ! -------------------------------------------------
        if( real_to_read == 0 ) then ! read in all realizations
            if( double .eqv. .true.) read(lin,*,end=10,err=99) data_double
            if( double .eqv. .false.) read(lin,*,end=10,err=99) data_single
        else
            do ir = 1, real_to_read ! read in only specified realization
                if( double .eqv. .true.) read(lin,*,end=10,err=99) data_double
                if( double .eqv. .false.) read(lin,*,end=10,err=99) data_single
                real_read = ir
            end do
        end if

        close(lin)
        return
 10     write(*,'(a,a)') 'ERROR, END OF FILE: ', trim(filename)
        close(lin)
        return
 99     write(*,'(a,a)') 'ERROR READING FILE: ',trim(filename)
        close(lin)
        return

    end subroutine read_gslib_grid

    subroutine read_gslib_point(filename, double)
    ! -------------------------------------------------------------------
    !
    ! Simple fortran subroutine for reading a file for pygeostat import
    !   - USE For size unkown Data I.E. point data_
    !
    ! author: Ryan Martin & Tyler Acorn
    ! -------------------------------------------------------------------
        character*(*), intent(in)   :: filename
        character*512               :: str
        integer                     :: lin,i, err, test
        logical                     :: double

        !standard gslib file IO to get the length of the data_ file  and the number
        !of variables to use to import using read_gslibfl below
        lin=1; datalength=0; num_var=0

        ! -------------------------------------------------
        ! Open File get data_ size
        ! -------------------------------------------------
        open(lin,file=trim(filename),status='OLD')
        read(lin,*,err=99)
        read(lin,*,err=99) num_var
        do i=1,num_var
           read(lin,*)
        end do
        do
            read( lin, *, end=9, err=99)
            datalength = datalength + 1
        end do
 9      rewind(lin)

        ! -------------------------------------------------
        ! Read header
        ! -------------------------------------------------
        do i = 1, num_var+2
            read(lin, *) str
        end do

        ! -------------------------------------------------
        ! Check Memory
        ! -------------------------------------------------
        if( allocated(data_double) ) deallocate(data_double)
        if( allocated(data_single) ) deallocate(data_single)

        if( double .eqv. .true.) then
            allocate( data_double(num_var, datalength), stat=err )
            if (err .ne. 0 ) stop ' ERROR: memory allocation failure! - data_double(num_var, datalength)'
        else
            allocate( data_single(num_var, datalength), stat=err )
            if (err .ne. 0 ) stop ' ERROR: memory allocation failure! - data_single(num_var, datalength)'
        end if

        ! -------------------------------------------------
        ! Read Data
        ! -------------------------------------------------
        if( double .eqv. .true.) read(lin,*,end=10,err=99) data_double
        if( double .eqv. .false.) read(lin,*,end=10,err=99) data_single

        close(lin)

        return
 10     write(*,'(a,a)') 'ERROR, END OF FILE: ', trim(filename)
        close(lin)
        return
 99     write(*,'(a,a,a,i2)') 'ERROR READING FILE: ',trim(filename),' error:', err
        close(lin)
        return
    end subroutine read_gslib_point

    subroutine write_gslib_r8(filename, datarray, num_var, num_dat, basefmt)
        ! -------------------------------------------------------------------
        !
        ! Simple fortran subroutine for buffered write to file - real*8 array
        !
        ! NOTE: subroutine assumes the header, nvar and columns are previously
        !       written to this file in python.
        !       For some reason this subroutine is wicked slow with GNU fortran
        !
        ! author: Ryan Martin
        ! -------------------------------------------------------------------
        character*(*), intent(in)   :: filename, basefmt
        integer, intent(in)         :: num_var, num_dat
        real*8, intent(in)          :: datarray(num_var, num_dat)
        character*100000            :: fmtst
        integer                     :: lout = 1, i

        !open the file for output
        open(lout,file=trim(filename),status='UNKNOWN')
        call build_formatstring( num_var, num_dat, datarray, basefmt, fmtst )
        ! ! write the format statement ! old one
        ! write(fmtst, '(a, i0, 2a)' ) '(', num_var, trim(basefmt), ')'
        ! loop past the pre-written header
        do i = 1, num_var + 2
            read(lout,*)
        enddo
        ! write the data and close
        write(lout, trim(fmtst)) datarray
        close(lout)
    end subroutine write_gslib_r8

    subroutine write_gslib_i4(filename, datarray, num_var, num_dat, fmstring)
        ! -------------------------------------------------------------------
        !
        ! Simple fortran subroutine for buffered write to file - integer data array
        !
        ! NOTE: subroutine assumes the header, nvar and columns are previously
        !       written to this file in python.
        !       For some reason this subroutine is wicked slow with GNU fortran
        !
        ! author: Ryan Martin
        ! -------------------------------------------------------------------
        character*(*), intent(in)   :: filename, fmstring
        integer, intent(in)         :: num_var, num_dat
        integer, intent(in)         :: datarray(num_var, num_dat)
        character*512               :: fmtst
        integer                     :: lout = 1, i

        !open the file for output
        open(lout,file=trim(filename),status='UNKNOWN')
        ! write the format statement
        write(fmtst, '(a, i0, 2a)' ) '(', num_var, trim(fmstring), ')'
        ! loop past the pre-written header
        do i = 1, num_var + 2
            read(lout,*)
        enddo
        ! write the data and close
        write(lout, fmtst) datarray
        close(lout)
    end subroutine write_gslib_i4

    subroutine build_formatstring( nvar, ndata, datarray, basefmt, formatstring )
    !----------------------------------------------------------------------------------------------
    ! Generalized function to build a format string for the case when a single real*8 array
    ! contains integers and space can be saved by writing out 'f0.0' for those integer values
    ! Parameters:
    !   nvar (int) :: the number of variables, first dimension of datarray
    !   ndata (int) :: the number of data, second dimension of datarray
    !   datarray (real*8, dim=2) :: nvar x ndata sized array
    !   basefmt (char) :: the base format string to consider should default 'f0.6'
    !
    ! .. codeauthor:: Ryan Martin - 10-17-2017
    !----------------------------------------------------------------------------------------------
        integer, intent(in) :: nvar, ndata
        real*8, intent(in) :: datarray(nvar, ndata)
        character*(*), intent(in) :: basefmt
        character*100000, intent(out) :: formatstring
        real*8 :: tempmod(ndata)
        integer :: ivar
        tempmod = 0.0d0
        formatstring = '('
        ! loop through variables and figure out what format is appropriate
        do ivar = 1, nvar
            if(len(trim(formatstring)) > 999950)then
                print *, " Exceeded max string length! " // &
                         " recompile fgslib_io with increased string format size "
                stop
            endif
            tempmod = mod( abs(datarray(ivar, :)), 1.0d0 )
            if(maxval(tempmod) == 0)then
                formatstring = trim(formatstring) // 'f0.0'  ! add an integer
            else
                formatstring = trim(formatstring) // trim(basefmt)  ! add the base fmt
            endif
            if(ivar /= nvar)then
                formatstring = trim(formatstring) // ',x,'
            endif
        enddo
        formatstring = trim(formatstring) // ')'
    end subroutine build_formatstring

    ! subroutine write_gslib_strfmt(filename, datarray, num_var, num_dat, fmt_array)
    !     ! -------------------------------------------------------------------
    !     !
    !     ! Simple fortran subroutine for writing a file from the given data
    !     ! array using the string formatting style to remove unneccessary
    !     ! whitespace
    !     !
    !     ! NOTE: subroutine assumes the header, nvar and columns are previously
    !     !       written to this file in python.
    !     !       Will be the slower write option but conserves space
    !     !
    !     ! author: Tyler Acorn
    !     ! -------------------------------------------------------------------
    !
    !     character*512, intent(in)           :: filename
    !     character*1000000, intent(inout)    :: str_write
    !     character*7 , intent(in)            :: fmt_array(num_var)
    !     integer, intent(in)                 :: num_var, num_dat
    !     real*8, intent(in)                  :: datarray(num_var, num_dat)
    !     character*512                       :: fmtst
    !     integer                             :: lout = 1, i_test, f_test
    !     integer                             :: var, dat
    !     integer                             :: vartype(num_var)
    !
    !     ! ----------- Open the File -------------------
    !     open(lout,file=trim(filename),status='UNKNOWN')
    !
    !     ! ---------------------------------------------
    !     ! Scan format string to see what
    !     ! is integer and what is float
    !     ! ---------------------------------------------
    !
    !     do var = 1, num_var
    !         ! Scan format to see if it is integer, or float
    !         i_test = SCAN(fmt_array(var). "iI")
    !         f_test = SCAN(fmt_array(var), "fgFG")
    !         if ( f_test .gt. 0 ) then
    !             vartype(var) = 0
    !         elseif ( i_test .gt. 0 ) then
    !             vartype(var) = 1
    !         else
    !             stop 'ERROR: Format string is not of type I, F, or G'
    !         end if
    !
    !     end do
    !
    !     ! ---------------------------------------------
    !     ! Loop through the data, format as string and
    !     ! write out the results
    !     ! ---------------------------------------------
    !
    !     do dat = 1, num_dat
    !         do var = 1, num_var
    !             str_len = 1
    !             str_write = ''
    !             if ( vartype(var) .eq. 0 ) call format_real(str_write, data_array(var, dat), fmt_array(var), str_len)
    !             if ( vartype(var) .eq. 1 ) call format_int(str_write, data_array(var, dat), fmt_array(var), str_len)
    !         end do
    !         write(lout, '(a)') str_write
    !     end do
    !
    !     return
    !
    ! end subroutine write_gslib_strfmt
    !
    ! subroutine format_real(str_write, r_variable, fmtr, str_len)
    !     implicit none
    !     character*1000000, intent(inout)    :: str_write
    !     character*15                        :: str_variable
    !     character*7, intent(in)             :: fmtr
    !     integer                             :: var_len
    !     integer, intent(inout)              :: str_len
    !     real, intent(in)                    :: r_variable
    !
    !     write(str_variable, fmtr) r_variable
    !     str_variable = adjustl( str_variable )
    !     var_len = len_trim( str_variable )
    !     str_write( str_len : str_len+var_len ) = str_variable( 1:var_len )
    !     str_len = str_len + var_len + 1
    !
    !     return
    !
    ! end subroutine format_real
    !
    ! subroutine format_int(str_write, i_variable, fmti, str_len)
    !     implicit none
    !     character*1000000, intent(inout)    :: str_write
    !     character*15                        :: str_variable
    !     character*7, intent(in)             :: fmti
    !     integer                             :: var_len
    !     integer, intent(in)                 :: i_variable
    !     integer, intent(inout)              :: str_len
    !
    !     write(str_variable, fmti) i_variable
    !     str_variable = adjustl( str_variable )
    !     var_len = len_trim( str_variable )
    !     str_write( str_len : str_len+var_len ) = str_variable( 1:var_len )
    !     str_len = str_len + var_len + 1
    !
    !     return
    !
    ! end subroutine format_int

    subroutine deallocate_arrays()
    ! -------------------------------------------------------------------
    !
    ! Subroutine to deallocate the module variables, since the memory is
    ! copied in python to the dataframe
    !
    ! author: Ryan Martin
    ! -------------------------------------------------------------------
    if( allocated(data_single) ) deallocate(data_single)
    if( allocated(data_double) ) deallocate(data_double)

    end subroutine deallocate_arrays

end module fortran_read
