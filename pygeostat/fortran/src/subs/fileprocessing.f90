! Module for general file processing operations
! Author: John Manchuk, 2010
!
! NOTES:
!   1. The file_open function has been written to return a valid
!      file handle so hard coding a handle value (ex. lin = 1, lout = 2)
!      is not required.  The function determines a valid handle and
!      returns it to the calling program.
!   2. Reading through a GSLIB style text file is straightforward using
!      the "datafilesize" and "readfileheader" routines.  The file does
!      not have to be pre-opened to call these.  The following is
!      recommended to read in a GSLIB file:
!           call datafilesize ( fid, fname, ncol, nrow )
!           ! allocate a real array ncol by nrow, ex.
!           allocate(filedata(ncol,nrow))
!           ncol = readfileheader ( fid, fname )
!           read(fid,*) filedata
!
! Latest version removes the datafilesize_stream functionality and writefile
! functionality, which requires the use of ifport... therefore making the module Intel-dependent
module filehandling

    implicit none
    
    public :: getparfile            !get the parameter file name from the user
    
    public :: file_exists, &        !check if a file exists
              file_open, &          !open a file and return the handle (integer function)
              readfileheader, &     !read through a gslib file header and return number of columns (integer function)
              datafilesize, &       !determine the number of rows and columns in a file
              checkname, &
              getsubstring, &
              sscanf, &
              getvarname
    
    public :: str2upper, str2lower
    
    ! public :: datafilesize_stream   !Beta to get the file size using steaming input (faster than datafilesize)
    !                                 !Not sure of portability
    private
    
    contains
    
    !Determine if a file exists
    logical function file_exists ( fname )
        character*(*), intent (in) :: fname
        logical :: does_exist
        file_exists = .false.
        inquire ( file = fname, exist = does_exist )
        if(does_exist) file_exists = .true.
    end function file_exists
    
    !Open a file
    ! This returns the handle value
    ! A handle of -1 indicates an error
    integer function file_open ( fname, create, access_ ) result ( fid )
        character*(*), intent (in) :: fname
        logical, intent (in) :: create
        character*(*), optional :: access_
        logical :: isopen
        integer :: fidt,ivar
        
        !get next available link
        fidt = 0
        isopen = .true.
        do while ( isopen )
            fidt = fidt + 1
            if(fidt == 6) fidt = fidt + 1 !skip command window handle
            inquire( fidt, OPENED = isopen )
        enddo
        
        fid = -1
        if(create)then
            fid = fidt
            if(present(access_))then
                open ( fid, file = fname, status = 'replace', access = access_, iostat=ivar )
            else
                open ( fid, file = fname, status = 'replace', iostat=ivar )
            endif
            if(ivar /= 0)then
                write(*,'(a,a)') 'ERROR OPENING ',trim(adjustl(fname))
                write(*,'(a)')   '  THE FILE MAY BE OPEN OR READ-ONLY'
                stop
            endif
        else
            if( .not. file_exists ( fname ) ) return
            fid = fidt
            if(present(access_))then
                open ( fid, file = fname, status = 'old', access = access_ )
            else
                open ( fid, file = fname, status = 'old' )
            endif
        endif
        return
    end function file_open
    
    !Get a parameter file from the user
    subroutine getparfile ( defaultname, parfilename, makepar )
        character*(*), intent (in) :: defaultname       !Default name of the parameter file
        character*(*), intent (inout) :: parfilename    !User input name
        character*512 :: str
        integer :: dlen
        logical :: does_exist
        external :: makepar
        
        dlen = len_trim(defaultname)
        
        str = ' '
        call getarg(1,str)
        if(str(1:1) == ' ')then
            write(*,*) 'Which parameter file do you want to use?'
            read (*,'(a)') str
        endif
        if(str(1:1) == ' ') str(1:dlen) = defaultname
        inquire ( file = str, exist = does_exist)
        if(does_exist)then
            parfilename = str
            return
        else
            write(*,*) 'ERROR - the parameter file does not exist,'
            write(*,*) '        check for the file and try again  '
            write(*,*)
            if(str(1:dlen) == defaultname) then
                write(*,*) '        creating a blank parameter file'
                call makepar ( )
                write(*,*)
            endif
            stop
        endif
    end subroutine getparfile
    
    !Read through a GSLIB file header
    integer function readfileheader ( fid, fname, vnames ) result ( ncol )
        integer, intent (inout) :: fid
        character*(*), intent (in) :: fname
        character*(*), intent (out), allocatable, optional :: vnames(:) ! RMB Allocatable array
        logical :: isopen, retnames
        character*512 :: str
        integer :: i, test, slen
        
        inquire(fid,OPENED = isopen)
        if(.not.isopen)then
            fid = file_open ( fname, .false. )
            if(fid < 0) return
        else
            rewind(fid)
        endif
        
        read(fid,'(a)',end=10,err=11) str
        read(fid,*,end=10,err=11) ncol
        ! RMB - 
        retnames = .false.
        if( present(vnames) )then
            retnames = .true.
            if( allocated(vnames) ) deallocate(vnames)
            allocate(vnames(ncol),stat=test)
            if( test /= 0 )then
                write(*,'(/,a,/)') 'ERROR ALLOCATING VNAMES ARRAY'
                stop
            endif
            do i = 1,ncol
                vnames(i) = '' !Weird bugs can occur with spaces if you don't initialize
            enddo
        endif
        do i = 1, ncol
            read(fid,'(a)',end=10,err=11) str
            if( retnames )then
                str = adjustl(str)
                slen = len_trim(str)
                vnames(i)(1:slen) = str(1:slen)
            endif
        enddo
        
        return
        
    10  write(*,'(a,a)') 'ERROR, END OF FILE: ',trim(fname)
        stop
    11  write(*,'(a,a)') 'ERROR READING FILE: ',trim(fname)
        stop
    end function readfileheader
    
    !Extract a substring from a list of variable names delimited by ::
    character*512 function getvarname ( varnames, j ) result ( vname )
        character*(*), intent (in) :: varnames
        integer, intent (in) :: j
        integer :: i, i0, k
        
        i0 = 1
        do i = 1, j
            k = index(varnames(i0:),'::')
            if(k == 0)then
                k = i0
                k = len_trim(varnames) + 1
                exit
            endif
            k = k + i0 - 1
            if(i < j ) i0 = k + 2
        enddo
        vname = ' '
        vname = varnames(i0:k-1)
    end function getvarname
        
    
    !Determine the size of a data file with the typical row-by-row method
    ! TODO - this should position the file back for loading the data after the call
    subroutine datafilesize ( fid, fname, ncol, nrow )
        integer, intent (inout) :: fid !file handle
        character*(*), intent (in) :: fname
        integer, intent (out) :: ncol
        integer*8, intent (out) :: nrow
        logical :: isopen
        real :: rec
        
        inquire(fid,OPENED = isopen)
        if(.not.isopen)then
            fid = file_open ( fname, .false. )
            if(fid < 0) return
        else
            rewind(fid)
        endif
                
        ncol = readfileheader ( fid, fname )
        nrow = 0
        
        do
            read(fid,*,end=9,err=11) rec
            nrow = nrow + 1
        enddo
    9   rewind(fid)
        
        return
    10  write(*,'(a,a)') 'ERROR: END OF FILE: ',trim(fname)
        stop
    11  write(*,'(a,a)') 'ERROR READING FILE: ',trim(fname)
        stop
    end subroutine datafilesize
    
    ! This subroutine takes the character string "str" of length "len" and
    ! removes all leading blanks and blanks out all characters after the
    ! first blank found in the string (leading blanks are removed first).
    subroutine checkname(temp)
        character*(*) :: temp
        integer       :: i

        ! Remove leading blanks:
        temp = adjustl(temp)
        
        ! find first two blanks and blank out remaining characters:
        i = index(temp,'  ')  ; if (i > 0) temp(i:)=' '

        ! Look for "-fi"
        i = index(temp,' -fi') ; if (i > 0) temp(i:)=' '

        ! Look for "\fi"
        i = index(temp,'\fi') ; if (i > 0) temp(i:)=' '
        return
    end subroutine checkname
    
    !Extract a particular substring from a longer string based on a
    !delimiter and index
    subroutine getsubstring ( str, k, delim, j0, j1 )
        character*(*), intent (in) :: str
        integer, intent (in) :: k
        character*(*), intent (in) :: delim
        integer, intent (out) :: j0, j1
        integer :: i, j, nd, nc
        
        nd = len_trim ( delim )
        nc = len_trim ( str )
        i = 0
        j = 1
        do while ( j > 0 )
            j1 = index ( str(j:nc), delim )
            i = i + 1
            if(i == k)then
                j0 = j
                if(j1 < 1)then
                    j1 = nc
                else
                    j1 = j0 + j1 - 2
                endif
                return
            endif
            j0 = j
            j = j0 + j1 + nd - 1
        enddo
        j0 = -1 !Not found
        j1 = -1
    end subroutine getsubstring
    
    
    !Scan a string for floating point numbers
    integer function sscanf ( s, f ) result ( n )
        character*(*), intent (in) :: s         !String to process
        real*8, optional, intent (out) :: f(:)  !Optional storage to pass back numbers
        integer :: i, k
        real*8 :: a
        
        k = len_trim ( s )
        n = 0
        i = 0
        SSL : do while ( i < k )
            i = i + 1
            if(s(i:i) /= ' ')then
                read(s(i:),*,err=50,end=49) a
                n = n + 1
                if(present(f)) f(n) = a
                do while ( s(i:i) /= ' ' )
                    i = i + 1
                    if(i > k) goto 50
                enddo
            endif
        enddo SSL
49      n = n + 1
        return
50      return  
    end function sscanf
    
    subroutine str2upper ( str )
        character*(*), intent (inout) :: str
        integer :: charcode, i
        do i = 1, len(str)
            charcode = iachar(str(i:i))
            if(charcode >= iachar('a') .and. charcode <= iachar('z') ) then
                str(i:i) = achar(iachar(str(i:i)) - 32)
            endif
        enddo
    end subroutine str2upper
    
    subroutine str2lower ( str )
        character*(*), intent (inout) :: str
        integer :: charcode, i
        do i = 1, len(str)
            charcode = iachar(str(i:i))
            if(charcode >= iachar('A') .and. charcode <= iachar('Z') ) then
                str(i:i) = achar(iachar(str(i:i)) + 32)
            endif
        enddo
    end subroutine str2lower

end module filehandling