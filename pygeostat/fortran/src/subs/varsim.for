c
c Module to declare dynamic arrays in multiple subroutines:
c
      module geostat
      
      use filehandling
      use gslib_binary
      
      real,allocatable         :: vr(:),sumgam(:,:,:),sumnp(:,:,:),
     +               sumhm(:,:,:),sumtm(:,:,:),sumhv(:,:,:),
     +               sumtv(:,:,:)

      real*8,allocatable      :: sills(:),gam(:),hm(:),tm(:),hv(:),
     +                           tv(:),np(:), data_(:,:),dataread(:,:),
     +                           vals(:,:) 
      integer*2, allocatable, dimension(:) :: datalith
      integer,parameter        :: lout=789,lout2=891
      integer,allocatable      :: ixd(:),iyd(:),izd(:),ivtail(:),
     +                            ivhead(:),ivtype(:),lith(:,:,:)
      character*12,allocatable :: names(:)
      character*64, allocatable, dimension(:) :: vnames
      real EPSLON,VERSION,xsiz,ysiz,zsiz,tmin,tmax
      integer  nlag,nx,ny,nz,nxy,nxyz,ndir,isill,nvarg,test,ireal,endfi
      integer  ncol, nrow, linlith, ivrlith, vcols(100), ldat, nreal
      character outfl*512,strdatafl*40,outfl2*512,lithfl*512,str1*512,
     +          parfl1*512
      logical :: isopen
      
      type(gsb_info) :: gsb, gsb_lith  ! General information
      type(gsb_var), allocatable :: gsbv(:), gsbv_lith(:)  
      
      end module
c
c
c
      module varsim_wrap
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C                                                                      %
C Copyright (C) 1996, The Board of Trustees of the Leland Stanford     %
C Junior University.  All rights reserved.                             %
C                                                                      %
C The programs in GSLIB are distributed in the hope that they will be  %
C useful, but WITHOUT ANY WARRANTY.  No author or distributor accepts  %
C responsibility to anyone for the consequences of using them or for   %
C whether they serve any particular purpose or work at all, unless he  %
C says so in writing.  Everyone is granted permission to copy, modify  %
C and redistribute the programs in GSLIB, but only under the condition %
C that this notice and the above copyright notice remain intact.       %
C                                                                      %
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c-----------------------------------------------------------------------
c
c                Variogram of Data on a Regular Grid
c                ***********************************
c
c This is a template driver program for GSLIB's "gam" subroutine.  The
c input data is ordered rowwise (x cycles fastest, then y, then z) in a
c GEOEAS format file.  The User's Guide contains more details.
c
c The program is executed with no command line arguments.  The user
c will be prompted for the name of a parameter file.  The parameter
c file is described in the documentation (see the example gam.par)
c
c
c The output file will contain each directional variogram ordered by
c direction and then variogram (the directions cycle fastest then the
c variogram number).  For each variogram there will be a one line
c description and then "nlag" lines with the following:
c
c        a) lag number (increasing from 1 to nlag)
c        b) separation distance
c        c) the "variogram" value
c        d) the number of pairs for the lag
c        e) the mean of the data contributing to the tail
c        f) the mean of the data contributing to the head
c
c Modified to calculate variograms for all realizations of a simulation
c
c AUTHOR: Clayton V. Deutsch                             DATE: 1989-1999
c MODIFIED: Chad Neufeld                                 DATE: 2002
c MODIFIED: Julian Ortiz                                 DATE: 2003
c MODIFIED: Jared Deutsch                                DATE: 2013
c     Fixed standardization of categorical indicator variograms
c MODIFIED: Jared Deutsch                                DATE: 2014
c     Converted to write out new Var program file format.
c MODIFIED: Jared Deutsch                                DATE: 2015
c     Fixed indexing of variograms to be 1,2,3,4...
c MODIFIED: Jared Deutsch                                DATE: 2015
c     Swapped order of loops to match with varcalc for ease of use
c MODIFIED: Ryan Barnett                                DATE: May 2015
c     Added GSB (binary) file detection and read capability 
c MODIFIED: Ryan Barnett                                DATE: Jan 2016
c     Fixed a couple of bugs relating to trimming of indicator variables 
c     and indicator variogram standardization
c-----------------------------------------------------------------------
      use geostat
      contains
      subroutine runvarsim(parfl)
        character*(*), intent(in) :: parfl

      EPSLON  = 1.0e-20
      VERSION = 1.500
c
c Note VERSION number:
c
      write(*,9999) VERSION
 9999 format(/' VarSim (conversion of gamsim_ave) Version: ',f5.3/)
c
c Read the Parameter File:
c
      parfl1 = trim(adjustl(parfl))

      ireal=0
      nreal = 10000
      do
      ireal=ireal+1
      
      if( ireal > nreal )then
        endfi = 1
        exit
      endif
      
      call readparm

c
c Call gam to compute the required variograms:
c
      write(*,*) 'working on realization number',ireal
      call gamma
c
c Write Results:
c
      call writeoutnew
      end do
      call writeoutnew
c
c Finished:
c
!      close(ldat)
!      close(lout)
c      write(*,*) 'Number of variograms = 'ireal-1
      write(*,9998) VERSION
 9998 format(/' VarSim Version: ',f5.3,
     +       ' Finished'/)
      !stop
      inquire(lout, OPENED = isopen)
      if( isopen ) then
        close(lout)
      endif
      inquire(lout2, OPENED = isopen)
      if( isopen ) then
        close(lout2)
      endif
      end subroutine
      end module
 
 
      subroutine readparm
c-----------------------------------------------------------------------
c
c                  Initialization and Read Parameters
c                  **********************************
c
c The input parameters and data are read in from their files. Some quick
c error checking is performed and the statistics of all the variables
c being considered are written to standard output.
c
c
c
c-----------------------------------------------------------------------
      use geostat
      use filehandling
      use gslib_binary
      parameter(MV=500)
      real      var(MV),cut(500)
      real*8    avg(MV),ssq(MV)
      integer   ivar(MV),num(MV),ivc(500),indflag(500)
      character datafl*512,str*512
      logical   testfl, usegsb, usegsb_lith
      data      lin/123/,ncut/0/
      save
     
      
c
c Declare dynamic arrays:
c
      real, allocatable :: vrmin(:),vrmax(:)
c
c Get the name of the parameter file - try the default name if no input:
c
      ! inquire( linlith, OPENED = isopen )
      ! if(ireal.eq.1) then
      ! do i=1,512
      !       str(i:i) = ' '
      ! end do
      ! call getarg(1,str)
      ! if(str(1:1).eq.' ')then
      !       write(*,*) 'Which parameter file do you want to use?'
      !       read (*,'(a)') str
      ! end if
      ! if(str(1:1).eq.' ') str(1:20) = 'varsim.par              '
      ! inquire(file=str,exist=testfl)
      ! if(.not.testfl) then
      !       write(*,*) 'ERROR - the parameter file does not exist,'
      !       write(*,*) '        check for the file and try again  '
      !       write(*,*)
      !       if(str(1:20).eq.'varsim.par              ') then
      !             write(*,*) '        creating a blank parameter file'
      !             call makepar
      !             write(*,*)
      !       end if
      !       stop
      ! endif
      !    strdatafl(1:40)=str(1:40)
      ! else
      !    str(1:40)=strdatafl(1:40)
      ! end if
      open(lin,file=parfl1,status='OLD')
c
c Find Start of Parameters:
c
 1    read(lin,'(a4)',end=98) str(1:4)
      if(str(1:4).ne.'STAR') go to 1
c
c Read Input Parameters:
c
      read(lin,'(a512)',err=98) lithfl
      call chknam(lithfl,512)
      if(ireal.eq.1)then
        usegsb_lith = checkGSB( lithfl )
        if( usegsb_lith )then
          write(*,*) ' GSB-style lithology file = ',lithfl(1:40)
        else
          write(*,*) ' GSLIB-style lithology file = ',lithfl(1:40)
        endif
      endif

      read(lin,*,err=98) ivrlith,ilith
      if(ireal.eq.1) 
     +write(*,*) ' lithology column, code = ',ivrlith,ilith

      read(lin,'(a512)',err=98) datafl
      call chknam(datafl,512)
      if(ireal.eq.1)then
        usegsb = checkGSB( datafl )
        if( usegsb )then
          write(*,*) ' GSB-style data file = ',datafl(1:40)
        else
          write(*,*) ' GSLIB-style data file = ',datafl(1:40)
        endif
      endif

      read(lin,*,err=98) nvar
      if(ireal.eq.1) 
     +write(*,*) ' number of variables = ',nvar
      backspace lin

      read(lin,*,err=98) j,(ivar(i),i=1,nvar)
      if(ireal.eq.1) 
     +write(*,*) ' columns = ',(ivar(i),i=1,nvar)

      read(lin,*,err=98) tmin,tmax
      if(ireal.eq.1) 
     +write(*,*) ' trimming limits = ',tmin,tmax

      read(lin,'(a512)',err=98) outfl
      call chknam(outfl,512)
      if(ireal.eq.1) 
     +write(*,*) ' output file = ',outfl(1:40)

      read(lin,'(a512)',err=98) outfl2
      call chknam(outfl2,512)
      if(ireal.eq.1) 
     +write(*,*) ' output file for average = ',outfl2(1:40)

      isim=1

      read(lin,*,err=98) nx,xmn,xsiz
      if(ireal.eq.1) 
     +write(*,*) ' nx,xmn,xsiz = ',nx,xmn,xsiz

      read(lin,*,err=98) ny,ymn,ysiz
      if(ireal.eq.1) 
     +write(*,*) ' ny,ymn,ysiz = ',ny,ymn,ysiz

      read(lin,*,err=98) nz,zmn,zsiz
      if(ireal.eq.1) 
     +write(*,*) ' nz,zmn,zsiz = ',nz,zmn,zsiz
     
      read(lin,*,err=98) nreal
      if(ireal.eq.1) 
     +write(*,*) ' nreal = ',nreal

      nxy  = nx * ny
      nxyz = nx * ny * nz
      if(nx.lt.1)         stop 'nx must be at least 1: check parameters'
      if(ny.lt.1)         stop 'ny must be at least 1: check parameters'
      if(nz.lt.1)         stop 'nz must be at least 1: check parameters'

      read(lin,*,err=98) ndir,nlag
      if(ireal.eq.1) 
     +write(*,*) ' ndir,nlag = ',ndir,nlag
      if(ndir.lt.1)       stop 'ndir is too small: check parameters'
      if(nlag.lt.1)       stop 'nlag is too small: check parameters'

      if(ireal.eq.1) then
      allocate (ixd(ndir),stat = test)
      if (test.ne.0) then
            write(*,*) 'Error: Allocation failed due to ',
     +                 'insufficient memory!', test
            stop
      end if  
      
      allocate (iyd(ndir),stat = test)
      if (test.ne.0) then
            write(*,*) 'Error: Allocation failed due to ',
     +                 'insufficient memory!', test
            stop
      end if  
     
      allocate (izd(ndir),stat = test)
      if (test.ne.0) then
            write(*,*) 'Error: Allocation failed due to ',
     +                 'insufficient memory!', test
            stop
      end if
      end if

      do i=1,ndir
            read(lin,*,err=98) ixd(i),iyd(i),izd(i)
            if(ireal.eq.1) 
     +write(*,*) ' direction = ',ixd(i),iyd(i),izd(i)
      end do
      read(lin,*,err=98) isill
      if(ireal.eq.1) 
     +write(*,*) ' flag to standardize sills = ',isill

      read(lin,*,err=98) nvarg
      if(ireal.eq.1) 
     +write(*,*) ' number of variograms = ',nvarg
      if(nvarg.lt.1)      stop 'nvarg is too small: check parameters'
c
c Allocate the needed memory:
c
      if(ireal.eq.1) then
      allocate (names(nvar+nvarg),stat = test)
      if (test.ne.0) then
            write(*,*) 'Error: Allocation failed due to ',
     +                 'insufficient memory!', test
            stop
      end if
      mxdlv = ndir*nlag*nvarg
      allocate (gam(mxdlv),stat = test) 
      if (test.ne.0) then
            write(*,*) 'Error: Allocation failed due to ',
     +                 'insufficient memory!', test
            stop
      end if      
      allocate (hm(mxdlv),stat = test)
      if (test.ne.0) then
            write(*,*) 'Error: Allocation failed due to ',
     +                 'insufficient memory!', test
            stop
      end if      
      allocate (tm(mxdlv),stat = test)
      if (test.ne.0) then
            write(*,*) 'Error: Allocation failed due to ',
     +                 'insufficient memory!', test
            stop
      end if      
      allocate (hv(mxdlv),stat = test)
      if (test.ne.0) then
            write(*,*) 'Error: Allocation failed due to ',
     +                 'insufficient memory!', test
            stop
      end if      
      allocate (tv(mxdlv),stat = test)  
      if (test.ne.0) then
            write(*,*) 'Error: Allocation failed due to ',
     +                 'insufficient memory!', test
            stop
      end if      
      allocate (np(mxdlv),stat = test)  
      if (test.ne.0) then
            write(*,*) 'Error: Allocation failed due to ',
     +                 'insufficient memory!', test
            stop
      end if      
      allocate (ivtail(nvarg+2),stat = test)
      if (test.ne.0) then
            write(*,*) 'Error: Allocation failed due to ',
     +                 'insufficient memory!', test
            stop
      end if      
      allocate (ivhead(nvarg+2),stat = test)
      if (test.ne.0) then
            write(*,*) 'Error: Allocation failed due to ',
     +                 'insufficient memory!', test
            stop
      end if      
      allocate (ivtype(nvarg+2),stat = test)
      if (test.ne.0) then
            write(*,*) 'Error: Allocation failed due to ',
     +                 'insufficient memory!', test
            stop
      end if      
      allocate (sumgam(nvarg,ndir,nlag),stat = test)
      if (test.ne.0) then
            write(*,*) 'Error: Allocation failed due to ',
     +                 'insufficient memory!', test
            stop
      end if      
      allocate (sumnp(nvarg,ndir,nlag),stat = test)
      if (test.ne.0) then
            write(*,*) 'Error: Allocation failed due to ',
     +                 'insufficient memory!', test
            stop
      end if      
      allocate (sumhm(nvarg,ndir,nlag),stat = test)
      if (test.ne.0) then
            write(*,*) 'Error: Allocation failed due to ',
     +                 'insufficient memory!', test
            stop
      end if      
      allocate (sumtm(nvarg,ndir,nlag),stat = test)
      if (test.ne.0) then
            write(*,*) 'Error: Allocation failed due to ',
     +                 'insufficient memory!', test
            stop
      end if      
      allocate (sumhv(nvarg,ndir,nlag),stat = test)
      if (test.ne.0) then
            write(*,*) 'Error: Allocation failed due to ',
     +                 'insufficient memory!', test
            stop
      end if      
      allocate (sumtv(nvarg,ndir,nlag),stat = test)
      if (test.ne.0) then
            write(*,*) 'Error: Allocation failed due to ',
     +                 'insufficient memory!', test
            stop
      end if      
      sumgam=0
      sumnp=0
      sumhm=0
      sumtm=0
      sumhv=0
      sumtv=0
      end if
c
      ncut = 0
      do i=1,nvarg
            read(lin,*,err=98) ivtail(i),ivhead(i),ivtype(i)
            if(ireal.eq.1) 
     +write(*,*) ' tail,head,type = ',
     +                   ivtail(i),ivhead(i),ivtype(i)
            if(ivtype(i).eq.9.or.ivtype(i).eq.10) then
                   ncut = ncut + 1
                   if(tmin.gt.0.0)stop'tmin interferes with indicators!'
                   if(tmax.le.1.0)stop'tmax interferes with indicators!'
                   backspace lin
                   read(lin,*,err=98) ii,jj,kk,cut(ncut)
                   if(ivtype(i).eq.9)  indflag(ncut) = 1
                   if(ivtype(i).eq.10) indflag(ncut) = 0
                   ivc(ncut) = ivtail(i)
                   ivtail(i) = nvar + ncut
                   ivhead(i) = nvar + ncut
                   write(names(nvar+ncut),140) ncut
 140               format('Indicator ',i2)
                   if(ireal.eq.1) 
     +write(*,*) ' indicator threshold:  ',cut(ncut)
            endif
      end do
      if(ireal.eq.1) 
     +write(*,*)
      close(lin)
c
c Determine the maximum size needed and allocate the needed memory:
c
      MAXVAR = nvar + ncut
      if(ireal.eq.1) then
      allocate (vr(nxyz*MAXVAR),stat = test)
      if (test.ne.0) then
            write(*,*) 'Error: Allocation failed due to ',
     +                 'insufficient memory!', test
            stop
      end if
      allocate (sills(MAXVAR),stat = test)
      if (test.ne.0) then
            write(*,*) 'Error: Allocation failed due to ',
     +                 'insufficient memory!', test
            stop
      end if
      end if

      allocate (vrmin(MAXVAR),stat = test)
      if (test.ne.0) then
            write(*,*) 'Error: Allocation failed due to ',
     +                 'insufficient memory!', test
            stop
      end if
c
      allocate (vrmax(MAXVAR),stat = test)
      if (test.ne.0) then
            write(*,*) 'Error: Allocation failed due to ',
     +                 'insufficient memory!', test
            stop
      end if

c
c Check lithology information
c
      if(ireal == 1) then
        allocate (lith(nx,ny,nz),stat = test)
            if (test.ne.0) then
                  write(*,*) 'Error: Allocation failed due to ',
     +                       'insufficient memory!', test
                  stop
            end if
        if(ivrlith.gt.0) then
c
c Read in the lithology information
c
            inquire(file=lithfl,exist=testfl)
            if(.not.testfl) then
                write(*,*) 'ERROR - the lithology file does not exist,'
                write(*,*) '        check for the file and try again  '
                stop
            endif
c
c The lithology file exists so open the file and read in the header
c information. Initialize the storage that will be used to summarize
c the data found in the file:
c
            allocate( datalith(nxyz), stat = test)
            if(test /= 0) go to 96 
            if( .not. usegsb_lith )then
                ncol = readfileheader( linlith, lithfl )
                if(ivrlith > ncol) then
                  write(*,*) ' ivr > the number in data file!'
                  stop
                end if
              allocate( dataread(ncol,nxyz), stat = test)
              if(test /= 0) go to 96 
            else ! if( .not. usegsb_lith )then
                test = read_bin_header( linlith, lithfl, gsb_lith, 
     +                                                     gsbv_lith )
            endif ! if( .not. usegsb_lith )then
            
        else ! if(ivrlith.gt.0) then
            lith=1
        end if ! if(ivrlith.gt.0) then
      end if !if(ireal.eq.1) then
      
      if( ivrlith > 0 )then
        if( .not. usegsb_lith )then
          read(linlith,*) dataread
          datalith = nint(dataread(ivrlith,:))
        else
          test = read_bin_data( linlith, gsb_lith, gsbv_lith,
     +                                                 ivrlith,ireal,1)
          test = gsbv2array( gsb_lith, gsbv_lith, ivrlith, datalith)
        endif
        lith=0
        i = 0
        j = 0
          do iz=1,nz
             do iy=1,ny
                do ix=1,nx
                 i = i + 1
                 if( datalith(i) == ilith )then
                     lith(ix,iy,iz) = 1
                     j = j + 1
                  endif
                end do
             end do
          end do
        write(str,*) j
        write(str1,*) ireal
      write(*,'(a,/,a)') '  found '//trim(adjustl(str))//' lith nodes',
     +                   '  for realization '//trim(adjustl(str1))
      endif
c
c Check to make sure the data file exists, then either read in the
c data or write an error message and stop:
c
      if(ireal == 1) then
      
        inquire(file=datafl,exist=testfl)
            if(.not.testfl) then
            write(*,*) 'ERROR data file ',datafl,' does not exist!'
            stop
        endif
c
c The data file exists so open the file and read in the header
c information. Initialize the storage that will be used to summarize
c the data found in the file:
c
        if( .not. usegsb )then
            ncol = readfileheader( ldat, datafl, vnames )
            do i=1,ncol
                do iv=1,nvar
                  j=ivar(iv)
                  if(i.eq.j) names(iv) = vnames(i)(1:12)
                end do
            enddo
            allocate( data_(ncol,nxyz), stat = test)
            if(test /= 0) go to 96
         else
            test = read_bin_header( ldat, datafl, gsb, gsbv )
            do i=1,gsb%nvar
                do iv=1,nvar
                  j=ivar(iv)
                  if(i.eq.j) names(iv) = gsb%vnames(i)(1:12)
                end do
            enddo
         endif
        allocate( vals(nvar,nxyz), stat = test)
        if(test /= 0) go to 96
      end if ! if(ireal.eq.1) then


      if( .not. usegsb )then
        read(ldat,*,err=99) data_
        vals = data_(ivar(1:nvar),:)
      else
        test = read_bin_data( ldat, gsb, gsbv, ivar(1:nvar), ireal, 1)
        test = gsbv2array( gsb, gsbv, ivar(1:nvar), vals  )
      endif
c
c Read the regular grid information row wise (x cycles fastest):
c
      num=0
      avg=0.0
      ssq=0.0
      i = 0
      do iz=1,nz
        do iy=1,ny
            do ix=1,nx
                i = i + 1
                do iv=1,nvar
                    index=ix+(iy-1)*nx+(iz-1)*nxy+(iv-1)*nxyz
                    vr(index) = vals(iv,i)
                    if(vr(index) >= tmin .and. vr(index) < tmax .and.
     +                           lith(ix,iy,iz) == 1) then
                            num(iv) = num(iv) + 1
                            avg(iv) = avg(iv) + dble(vr(index))
                            ssq(iv) = ssq(iv)
     +                                      + dble(vr(index)*vr(index))
                    else
                      vr(index) = tmin - 1;
                    end if
               end do
            end do
         end do
      end do
c
c Compute the averages and variances as an error check for the user:
c
      do iv=1,nvar
            sills(iv) = -999.
            if(num(iv).gt.0) then
                  avg(iv) = avg(iv) / dble(num(iv))
                  ssq(iv) =(ssq(iv) / dble(num(iv))) - avg(iv) * avg(iv)
                  sills(iv) = ssq(iv)
                  if(ireal.eq.1) 
     +            write(*,*) 'Variable number ',iv
                  if(ireal.eq.1) 
     +            write(*,*) '  Number   = ',num(iv)
                  if(ireal.eq.1) 
     +            write(*,*) '  Average  = ',real(avg(iv))
                  if(ireal.eq.1) 
     +            write(*,*) '  Variance = ',real(ssq(iv))
            endif
      end do
c
c Construct Indicator Variables if necessary:
c
      do ic=1,ncut
            iv   = ivc(ic)
            jv   = nvar + ic
            ptot = 0.0
            p1   = 0.0
            do ix=1,nx
            do iy=1,ny
            do iz=1,nz
                  index = ix+(iy-1)*nx+(iz-1)*nxy+(iv-1)*nxyz
                  jndex = ix+(iy-1)*nx+(iz-1)*nxy+(jv-1)*nxyz
                  if(vr(index).lt.tmin.or.vr(index).ge.tmax) then
                        !vr(jndex) = tmin - EPSLON 
                        vr(jndex) = tmin - 1  ! EPSLON causing issues later
                  else
                        if(indflag(ic).eq.1) then
                              if(vr(index).lt.cut(ic)) then
                                    vr(jndex) = 0.0
                              else
                                    vr(jndex) = 1.0
                              end if
                              p1   = p1   + vr(index)
                              ptot = ptot + 1.0
                        else
                          vr(jndex) = 0.0
                          if(int(vr(index)+0.5).eq.int(cut(ic)+0.5))then
                            vr(jndex) = 1.0
                            p1   = p1   + vr(jndex) ! This wasn't previously within an if statement
                          endif
                          ptot = ptot + 1.0
                        end if
                  end if
            end do
            end do
            end do
            p1        = p1 / max(ptot,1.0)
            sills(jv) = dble (p1*(1.0-p1))
      end do
c
c Establish minimums and maximums:
c
      do i=1,MAXVAR
            vrmin(i) =  1.0e21
            vrmax(i) = -1.0e21
      end do
      do ix=1,nx
       do iy=1,ny
        do iz=1,nz
         do iv=1,nvar+ncut
          index = ix+(iy-1)*nx+(iz-1)*nxy+(iv-1)*nxyz
          if(vr(index).ge.tmin.and.vr(index).lt.tmax) then
            if(vr(index).lt.vrmin(iv)) vrmin(iv) = vr(index)
            if(vr(index).gt.vrmax(iv)) vrmax(iv) = vr(index)
          end if
         end do 
        end do 
       end do 
      end do 
c
c Check on the variogams that were requested:
c
      call check(vrmin,vrmax)
      deallocate (vrmin,stat = test)
      deallocate (vrmax,stat = test)
c
c Return:
c
      return
c
c Error in an Input File Somewhere:
c
 96   stop 'ERROR allocating array'
 98   stop 'ERROR in parameter file!'
 99   stop 'ERROR in data file!'

      end subroutine



      subroutine gamma
c-----------------------------------------------------------------------
c
c                Variogram of Data on a Regular Grid
c                ***********************************
c
c This subroutine computes any of eight different measures of spatial
c continuity for regular spaced 3-D data.  Missing values are allowed
c and the grid need not be cubic.
c
c
c
c INPUT VARIABLES:
c
c   nlag             Maximum number of lags to be calculated
c   nx               Number of units in x (number of columns)
c   ny               Number of units in y (number of lines)
c   nz               Number of units in z (number of levels)
c   ndir             Number of directions to consider
c   ixd(ndir)        X (column) indicator of direction - number of grid
c                      columns that must be shifted to move from a node
c                      on the grid to the next nearest node on the grid
c                      which lies on the directional vector
c   iyd(ndir)        Y (line) indicator of direction - similar to ixd,
c                      number of grid lines that must be shifted to
c                      nearest node which lies on the directional vector
c   izd(ndir)        Z (level) indicator of direction - similar to ixd,
c                      number of grid levels that must be shifted to
c                      nearest node of directional vector
c   nv               The number of variables
c   vr(nx*ny*nz*nv)  Array of data
c   tmin,tmax        Trimming limits
c   isill            1=attempt to standardize, 0=do not
c   sills            the sills (variances) to standardize with
c   nvarg            Number of variograms to compute
c   ivtail(nvarg)    Variable for the tail of the variogram
c   ivhead(nvarg)    Variable for the head of the variogram
c   ivtype(nvarg)    Type of variogram to compute:
c                      1. semivariogram
c                      2. cross-semivariogram
c                      3. covariance
c                      4. correlogram
c                      5. general relative semivariogram
c                      6. pairwise relative semivariogram
c                      7. semivariogram of logarithms
c                      8. madogram
c                      9. indicator semivariogram: an indicator variable
c                         is constructed in the main program.
c
c OUTPUT VARIABLES:  The following arrays are ordered by direction,
c                    then variogram, and finally lag, i.e.,
c                      iloc = (id-1)*nvarg*nlag+(iv-1)*nlag+il
c
c   np()             Number of pairs
c   gam()            Semivariogram, covariance, correlogram,... value
c   hm()             Mean of the tail data
c   tm()             Mean of the head data
c   hv()             Variance of the tail data
c   tv()             Variance of the head data
c
c
c
c Original:  A.G. Journel                                           1978
c Revisions: B.E. Buxton                                       Apr. 1982
c-----------------------------------------------------------------------
      use geostat
c
c Initialize the summation arrays for each direction, variogram, and lag
c
      nxyz = nx*ny*nz
      nsiz = ndir*nvarg*nlag
      do i=1,nsiz
            np(i)  = 0.
            gam(i) = 0.0
            hm(i)  = 0.0
            tm(i)  = 0.0
            hv(i)  = 0.0
            tv(i)  = 0.0
      end do
c
c First fix the location of a seed point on the grid (ix,iy,iz):
c
      do ix=1,nx
      do iy=1,ny
      do iz=1,nz
c
c For the fixed seed point, loop through all directions:
c
            do id=1,ndir
              ixinc = ixd(id)
              iyinc = iyd(id)
              izinc = izd(id)
              ix1   = ix
              iy1   = iy
              iz1   = iz
c
c For this direction, loop over all the lags:
c
              do il=1,nlag
c
c Check to be sure that the point being considered is still in the
c grid - if not, then finished with this direction:
c
                ix1 = ix1 + ixinc
                if(ix1.lt.1.or.ix1.gt.nx) go to 3
                iy1 = iy1 + iyinc
                if(iy1.lt.1.or.iy1.gt.ny) go to 3
                iz1 = iz1 + izinc
                if(iz1.lt.1.or.iz1.gt.nz) go to 3
c
c For this direction and lag, loop over all variograms:
c
                do iv=1,nvarg
                  it = ivtype(iv)
c
c Get the head value, skip this value if missing:
c
                  i     = ivhead(iv)
                  index = ix+(iy-1)*nx+(iz-1)*nxy+(i-1)*nxyz
                  vrt   = vr(index)
                  if(vrt.lt.tmin.or.vrt.ge.tmax) go to 5
c
c Get the tail value, skip this value if missing:
c
                  i     = ivtail(iv)
                  index = ix1+(iy1-1)*nx+(iz1-1)*nxy+(i-1)*nxyz
                  vrh   = vr(index)
                  if(vrh.lt.tmin.or.vrh.ge.tmax) go to 5
c
c Need increment for the cross semivariogram only:
c
                  if(it.eq.2) then
                        i     = ivtail(iv)
                        index = ix+(iy-1)*nx+(iz-1)*nxy+(i-1)*nxyz
                        vrhpr = vr(index)
                        if(vrhpr.lt.tmin.or.vrhpr.ge.tmax) go to 5
                        i     = ivhead(iv)
                        index = ix1+(iy1-1)*nx+(iz1-1)*nxy+(i-1)*nxyz
                        vrtpr = vr(index)
                        if(vrtpr.lt.tmin.or.vrtpr.ge.tmax) go to 5
                  endif
c
c We have an acceptable pair, therefore accumulate all the statistics
c that are required for the variogram:
c
                  i      = (id-1)*nvarg*nlag+(iv-1)*nlag+il
                  np(i)  = np(i) + 1.
                  tm(i)  = tm(i) + dble(vrt)
                  hm(i)  = hm(i) + dble(vrh)
c
c Choose the correct variogram type and keep relevant sums:
c
                  if(it.eq.1.or.it.ge.9) then
                    gam(i) = gam(i) + dble((vrh-vrt)*(vrh-vrt))
                  else if(it.eq.2) then
                    gam(i) = gam(i) + dble((vrhpr-vrh)*(vrt-vrtpr))
                  else if(abs(it).eq.3) then
                    gam(i) = gam(i) +  dble(vrh*vrt)
                  else if(it.eq.4) then
                    gam(i) = gam(i) +  dble(vrh*vrt)
                    hv(i)  = hv(i)  +  dble(vrh*vrh)
                    tv(i)  = tv(i)  +  dble(vrt*vrt)
                  else if(it.eq.5) then
                    gam(i) = gam(i) + dble((vrh-vrt)*(vrh-vrt))
                  else if(it.eq.6) then
                    if((vrt+vrh).lt.EPSLON) then
                        np(i)  = np(i) - 1.
                        tm(i)  = tm(i) - dble(vrt)
                        hm(i)  = hm(i) - dble(vrh)
                    else
                        tempvar= 2.0*(vrt-vrh)/(vrt+vrh)
                        gam(i) = gam(i) + dble(tempvar*tempvar)
                    endif
                  else if(it.eq.7) then
                    if(vrt.lt.EPSLON.or.vrh.lt.EPSLON) then
                        np(i)  = np(i) - 1.
                        tm(i)  = tm(i) - dble(vrt)
                        hm(i)  = hm(i) - dble(vrh)
                    else
                        tempvar= alog(vrt)-alog(vrh)
                        gam(i) = gam(i) + dble(tempvar*tempvar)
                    endif
                  else if(it.eq.8) then
                    gam(i) = gam(i) + dble(abs(vrt-vrh))
                  endif
 5              continue
                end do
 4            continue
              end do
 3          continue
            end do
      end do
      end do
      end do
c
c Get average values for gam, hm, tm, hv, and tv, then compute
c the correct "variogram" measure:
c
      do id=1,ndir
      do iv=1,nvarg
      do il=1,nlag
            i = (id-1)*nvarg*nlag+(iv-1)*nlag+il
            if(np(i).eq.0.) go to 6
            rnum   = np(i)
            gam(i) = gam(i) / dble(rnum)
            hm(i)  = hm(i)  / dble(rnum)
            tm(i)  = tm(i)  / dble(rnum)
            hv(i)  = hv(i)  / dble(rnum)
            tv(i)  = tv(i)  / dble(rnum)
            it     = ivtype(iv)
c
c Attempt to standardize:
c

            if(isill.eq.1) then
                  if(ivtail(iv).eq.ivhead(iv)) then
                        iii = ivtail(iv)
                        if((it.eq.1.or.it.ge.9).and.sills(iii).gt.0.0)
     +                    gam(i) = gam(i) / sills(iii)
                  end if
            end if
c
c 1. report the semivariogram rather than variogram
c 2. report the cross-semivariogram rather than variogram
c 3. the covariance requires "centering"
c 4. the correlogram requires centering and normalizing
c 5. general relative requires division by lag mean
c 6. report the semi(pairwise relative variogram)
c 7. report the semi(log variogram)
c 8. report the semi(madogram)
c
            if(it.eq.1.or.it.eq.2) then
                  gam(i) = 0.5 * gam(i)
            else if(abs(it).eq.3) then
                  gam(i) = gam(i) - hm(i)*tm(i)
                  if(it.lt.0) then
                        if(sills(ivtail(iv)).lt.0.0.or.
     +                     sills(ivhead(iv)).lt.0.0) then
                              gam(i) = -999.0
                        else
                              variance = ( sqrt(sills(ivtail(iv)))
     +                                 *   sqrt(sills(ivhead(iv))) )
                              gam(i) = variance - gam(i)
                        end if
                  end if
            else if(it.eq.4) then
                  hv(i)  = hv(i)-hm(i)*hm(i)
                  if(hv(i).le.0.0) hv(i) = 0.0
                  hv(i)  = sqrt(hv(i))
                  tv(i)  = tv(i)-tm(i)*tm(i)
                  if(tv(i).le.0.0) tv(i) = 0.0
                  tv(i)  = sqrt(tv(i))
                  if((hv(i)*tv(i)).lt.EPSLON) then
                        gam(i) = 0.0
                  else
                        gam(i) =(gam(i)-hm(i)*tm(i))/(hv(i)*tv(i))
                  endif
c
c Square "hv" and "tv" so that we return the variance:
c
                  hv(i)  = hv(i)*hv(i)
                  tv(i)  = tv(i)*tv(i)
            else if(it.eq.5) then
                  htave  = 0.5*(hm(i)+tm(i))
                  htave  = htave   *   htave
                  if(htave.lt.EPSLON) then
                        gam(i) = 0.0
                  else
                        gam(i) = gam(i)/dble(htave)
                  endif
            else if(it.ge.6) then
                  gam(i) = 0.5 * gam(i)
            endif
 6    continue
      end do
      end do
      end do
      return
      end subroutine
 
 
 
      subroutine writeout
c-----------------------------------------------------------------------
c
c                  Write Out the Results of GAM
c                  ****************************
c
c An output file will be written which contains each directional
c variogram ordered by direction and then variogram (the directions
c cycle fastest then the variogram number).  For each variogram there
c will be a one line description and then "nlag" lines with:
c
c        a) lag number (increasing from 1 to nlag)
c        b) separation distance
c        c) the "variogram" value
c        d) the number of pairs for the lag
c        e) the mean of the data contributing to the tail
c        f) the mean of the data contributing to the head
c        g) IF the correlogram - variance of tail values
c        h) IF the correlogram - variance of head values
c
c
c-----------------------------------------------------------------------
      use geostat

      character title*100,titleavg*100
c      data      lout/1/
c
c Loop over all the variograms that have been computed:
c
      if(ireal.eq.1) then
      open(lout,file=outfl,status='UNKNOWN')
      open(lout2,file=outfl2,status='UNKNOWN')
      endif
c
      if(endfi.ne.1) then
      do iv=1,nvarg
c
c Construct a title that reflects the variogram type and the variables
c that were used to calculate the variogram:
c
      it = abs(ivtype(iv))
      if(it.eq. 1) title(1:24) = 'Semivariogram           '
      if(it.eq. 2) title(1:24) = 'Cross Semivariogram     '
      if(it.eq. 3) title(1:24) = 'Covariance              '
      if(it.eq. 4) title(1:24) = 'Correlogram             '
      if(it.eq. 5) title(1:24) = 'General Relative        '
      if(it.eq. 6) title(1:24) = 'Pairwise Relative       '
      if(it.eq. 7) title(1:24) = 'Variogram of Logarithms '
      if(it.eq. 8) title(1:24) = 'Semimadogram            '
      if(it.eq. 9) title(1:24) = 'Indicator 1/2 Variogram '
      if(it.eq.10) title(1:24) = 'Indicator 1/2 Variogram '
      write(title(25:79),100) names(ivtail(iv)),names(ivhead(iv)),ireal
 100  format('tail:',a12,' head:',a12,'realization:'i4)
c
c Loop over all the directions (note the direction in the title):
c
      do id=1,ndir
            write(title(80:92),101) id
 101        format('direction ',i2)
            write(lout,'(a)') title(1:92)
c
c Compute the unit lag distance along the directional vector:
c
            dis = sqrt( max(((ixd(id)*xsiz)**2 + (iyd(id)*ysiz)**2 +
     +                  (izd(id)*zsiz)**2),0.0) )
c
c Write out all the lags:
c
            do il=1,nlag
                  i = (id-1)*nvarg*nlag+(iv-1)*nlag+il
                  disl = real(il)*dis
                  nump = int(np(i))
                  if(it.eq.4) then
                        write(lout,102) il,disl,gam(i),nump,
     +                                  hm(i),tm(i),hv(i),tv(i)

c
c Compute average variogram statistics
c            
                   sumgam(iv,id,il)=sumgam(iv,id,il)+gam(i)
                   sumnp(iv,id,il)=sumnp(iv,id,il)+nump
                   sumhm(iv,id,il)=sumhm(iv,id,il)+hm(i)
                   sumtm(iv,id,il)=sumtm(iv,id,il)+tm(i)
                   sumhv(iv,id,il)=sumhv(iv,id,il)+hv(i)
                   sumtv(iv,id,il)=sumtv(iv,id,il)+tv(i)
                  else
                        write(lout,102) il,disl,gam(i),nump,
     +                                  hm(i),tm(i)
                   sumgam(iv,id,il)=sumgam(iv,id,il)+gam(i)
                   sumnp(iv,id,il)=sumnp(iv,id,il)+nump
                   sumhm(iv,id,il)=sumhm(iv,id,il)+hm(i)
                   sumtm(iv,id,il)=sumtm(iv,id,il)+tm(i)
                  endif
 102              format(1x,i3,1x,f12.3,1x,f12.5,1x,i8,4(1x,f14.5))
            end do
      end do
      end do
      end if
c
c Write out averages
c
      if(endfi.eq.1) then 
         nsim=ireal-1
      do iv=1,nvarg
c
c Construct a title that reflects the variogram type and the variables
c that were used to calculate the variogram:
c
      it = abs(ivtype(iv))
      if(it.eq. 1) titleavg(1:24) = 'Semivariogram           '
      if(it.eq. 2) titleavg(1:24) = 'Cross Semivariogram     '
      if(it.eq. 3) titleavg(1:24) = 'Covariance              '
      if(it.eq. 4) titleavg(1:24) = 'Correlogram             '
      if(it.eq. 5) titleavg(1:24) = 'General Relative        '
      if(it.eq. 6) titleavg(1:24) = 'Pairwise Relative       '
      if(it.eq. 7) titleavg(1:24) = 'Variogram of Logarithms '
      if(it.eq. 8) titleavg(1:24) = 'Semimadogram            '
      if(it.eq. 9) titleavg(1:24) = 'Indicator 1/2 Variogram '
      if(it.eq.10) titleavg(1:24) = 'Indicator 1/2 Variogram '
      write(titleavg(25:76),103) names(ivtail(iv)),names(ivhead(iv))
 103  format('tail:',a12,' head:',a12)
c
c Loop over all the directions (note the direction in the title):
c
      do id=1,ndir
            write(titleavg(76:88),104) id
 104        format('direction ',i2)
            write(lout2,'(a)') titleavg(1:88)
c
c Compute the unit lag distance along the directional vector:
c
            dis = sqrt( max(((ixd(id)*xsiz)**2 + (iyd(id)*ysiz)**2 +
     +                  (izd(id)*zsiz)**2),0.0) )
c
c Write out all the lags:
c
            do il=1,nlag
                  i = (id-1)*nvarg*nlag+(iv-1)*nlag+il
                  disl = real(il)*dis
                  nump = int(np(i))
                  if(it.eq.4) then
                        write(lout2,105) il,disl,sumgam(iv,id,il)/nsim,
     +                  int(sumnp(iv,id,il)/nsim),sumhm(iv,id,il)/nsim,
     +                  sumtm(iv,id,il)/nsim,sumhv(iv,id,il)/nsim,
     +                  sumtv(iv,id,il)/nsim
                  else
                        write(lout2,105) il,disl,sumgam(iv,id,il)/nsim,
     +                  int(sumnp(iv,id,il)/nsim),sumhm(iv,id,il)/nsim,
     +                  sumtm(iv,id,il)/nsim
                  endif
 105              format(1x,i3,1x,f12.3,1x,f12.5,1x,i8,4(1x,f14.5))
            end do
      end do
      end do
      end if
c      close(lout)
      return
      end subroutine
 
      subroutine writeoutnew
c-----------------------------------------------------------------------
c
c         Write Out the Results of GAM - using the new format
c         ***************************************************
c
c An output file will be written which contains each directional
c variogram ordered by direction and then variogram (the directions
c cycle fastest then the variogram number).  For each variogram there
c will be a one line description and then "nlag" lines with:
c
c        a) lag number (increasing from 1 to nlag)
c        b) separation distance
c        c) the "variogram" value
c        d) the number of pairs for the lag
c        e) the mean of the data contributing to the tail
c        f) the mean of the data contributing to the head
c        g) IF the correlogram - variance of tail values
c        h) IF the correlogram - variance of head values
c
c
c-----------------------------------------------------------------------
      use geostat
      real, parameter :: RAD2DEG = 57.2957795

      integer :: varidx
      real :: azm,dip

      character title*100,titleavg*100
c      data      lout/1/
c
c Loop over all the variograms that have been computed:
c
      if(ireal.eq.1) then
        open(lout,file=outfl,status='UNKNOWN')
        write(lout,'(a)') 'Gridded Variograms'
        write(lout,'(a)') '8'
        write(lout,'(a)') 'Variogram Index'
        write(lout,'(a)') 'Lag Distance'
        write(lout,'(a)') 'Number of Pairs'
        write(lout,'(a)') 'Variogram Value'
        write(lout,'(a)') 'Variogram Number'
        write(lout,'(a)') 'Calculation Azimuth'
        write(lout,'(a)') 'Calculation Dip'
        write(lout,'(a)') 'Variogram Type'

        open(lout2,file=outfl2,status='UNKNOWN')
        write(lout2,'(a)') 'Average Gridded Variograms'
        write(lout2,'(a)') '8'
        write(lout2,'(a)') 'Variogram Index'
        write(lout2,'(a)') 'Lag Distance'
        write(lout2,'(a)') 'Number of Pairs'
        write(lout2,'(a)') 'Variogram Value'
        write(lout2,'(a)') 'Variogram Number'
        write(lout2,'(a)') 'Calculation Azimuth'
        write(lout2,'(a)') 'Calculation Dip'
        write(lout2,'(a)') 'Variogram Type'
      endif
c
      if(endfi.ne.1) then


      varidx = 0

      do id=1,ndir
      
cJLD - swapped direction and variable places
          
      do iv=1,nvarg
c
c Construct a title that reflects the variogram type and the variables
c that were used to calculate the variogram:
c
      it = abs(ivtype(iv))
      if(it.eq. 1) title(1:24) = 'Semivariogram           '
      if(it.eq. 2) title(1:24) = 'Cross Semivariogram     '
      if(it.eq. 3) title(1:24) = 'Covariance              '
      if(it.eq. 4) title(1:24) = 'Correlogram             '
      if(it.eq. 5) title(1:24) = 'General Relative        '
      if(it.eq. 6) title(1:24) = 'Pairwise Relative       '
      if(it.eq. 7) title(1:24) = 'Variogram of Logarithms '
      if(it.eq. 8) title(1:24) = 'Semimadogram            '
      if(it.eq. 9) title(1:24) = 'Indicator 1/2 Variogram '
      if(it.eq.10) title(1:24) = 'Indicator 1/2 Variogram '
      write(title(25:79),100) names(ivtail(iv)),names(ivhead(iv)),ireal
 100  format('tail:',a12,' head:',a12,'realization:'i4)
c
c Loop over all the directions (note the direction in the title):
c
         
          
          
        write(title(80:92),101) id
 101    format('direction ',i2)
c        write(lout,'(a)') title(1:92)
c
c Compute the unit lag distance along the directional vector:
c
        dis = sqrt( max(((ixd(id)*xsiz)**2 + (iyd(id)*ysiz)**2 +
     +              (izd(id)*zsiz)**2),0.0) )
        if (abs(iyd(id)*ysiz) .gt. 0.0) then
          azm = atan((ixd(id)*xsiz)/(iyd(id)*ysiz))*RAD2DEG
        else
          if (ixd(id) > 0)then
            azm = 90.0
          else
            azm = 270.0
          end if
        end if
        if (sqrt((ixd(id)*xsiz)**2 + (iyd(id)*ysiz)**2) .gt. 0.0) then
          dip = atan((izd(id)*zsiz)/
     +          sqrt((ixd(id)*xsiz)**2 + (iyd(id)*ysiz)**2))*RAD2DEG
        else
          dip = 90.0
        end if
c
c Write out all the lags:
c
        do il=1,nlag
          varidx = varidx + 1
          i = (id-1)*nvarg*nlag+(iv-1)*nlag+il
          disl = real(il)*dis
          nump = int(np(i))
          if(it.eq.4) then
c            write(lout,102) il,disl,gam(i),nump,
c     +                                  hm(i),tm(i),hv(i),tv(i)

c
c Compute average variogram statistics
c            
            sumgam(iv,id,il)=sumgam(iv,id,il)+gam(i)
            sumnp(iv,id,il)=sumnp(iv,id,il)+nump
            sumhm(iv,id,il)=sumhm(iv,id,il)+hm(i)
            sumtm(iv,id,il)=sumtm(iv,id,il)+tm(i)
            sumhv(iv,id,il)=sumhv(iv,id,il)+hv(i)
            sumtv(iv,id,il)=sumtv(iv,id,il)+tv(i)
          else
c            write(lout,102) il,disl,gam(i),nump,
c     +                          hm(i),tm(i)
            sumgam(iv,id,il)=sumgam(iv,id,il)+gam(i)
            sumnp(iv,id,il)=sumnp(iv,id,il)+nump
            sumhm(iv,id,il)=sumhm(iv,id,il)+hm(i)
            sumtm(iv,id,il)=sumtm(iv,id,il)+tm(i)
          endif
          write(lout,'(i5,f16.5,i12,f16.5,i6,2f12.4,i4)') 
     +                       (ireal-1)*(nvarg*ndir)+(id-1)*nvarg+iv,
     +                       disl,
     +                       nump,gam(i),(id-1)*nvarg+iv,azm,dip,it
 102      format(1x,i3,1x,f12.3,1x,f12.5,1x,i8,4(1x,f14.5))
        end do
      end do
      end do
      end if
c
c Write out averages
c
      if(endfi.eq.1) then 
         nsim=ireal-1
         
c
c Loop over all the directions (note the direction in the title):
c
      do id=1,ndir
         
         
cJLD Swapped positions         
         
      do iv=1,nvarg
          
c
c Construct a title that reflects the variogram type and the variables
c that were used to calculate the variogram:
c
      it = abs(ivtype(iv))
      if(it.eq. 1) titleavg(1:24) = 'Semivariogram           '
      if(it.eq. 2) titleavg(1:24) = 'Cross Semivariogram     '
      if(it.eq. 3) titleavg(1:24) = 'Covariance              '
      if(it.eq. 4) titleavg(1:24) = 'Correlogram             '
      if(it.eq. 5) titleavg(1:24) = 'General Relative        '
      if(it.eq. 6) titleavg(1:24) = 'Pairwise Relative       '
      if(it.eq. 7) titleavg(1:24) = 'Variogram of Logarithms '
      if(it.eq. 8) titleavg(1:24) = 'Semimadogram            '
      if(it.eq. 9) titleavg(1:24) = 'Indicator 1/2 Variogram '
      if(it.eq.10) titleavg(1:24) = 'Indicator 1/2 Variogram '
      write(titleavg(25:76),103) names(ivtail(iv)),names(ivhead(iv))
 103  format('tail:',a12,' head:',a12)

            write(titleavg(76:88),104) id
 104        format('direction ',i2)
c            write(lout2,'(a)') titleavg(1:88)
c
c Compute the unit lag distance along the directional vector:
c
            dis = sqrt( max(((ixd(id)*xsiz)**2 + (iyd(id)*ysiz)**2 +
     +                  (izd(id)*zsiz)**2),0.0) )
        if (abs(iyd(id)*ysiz) .gt. 0.0) then
          azm = atan((ixd(id)*xsiz)/(iyd(id)*ysiz))*RAD2DEG
        else
          azm = 90.0
        end if
        if (sqrt((ixd(id)*xsiz)**2 + (iyd(id)*ysiz)**2) .gt. 0.0) then
          dip = atan((izd(id)*zsiz)/
     +          sqrt((ixd(id)*xsiz)**2 + (iyd(id)*ysiz)**2))*RAD2DEG
        else
          dip = 90.0
        end if
c
c Write out all the lags:
c
            do il=1,nlag
                  i = (id-1)*nvarg*nlag+(iv-1)*nlag+il
                  disl = real(il)*dis
                  nump = int(np(i))
c                  if(it.eq.4) then
c                        write(lout2,105) il,disl,sumgam(iv,id,il)/nsim,
c     +                  int(sumnp(iv,id,il)/nsim),sumhm(iv,id,il)/nsim,
c     +                  sumtm(iv,id,il)/nsim,sumhv(iv,id,il)/nsim,
c     +                  sumtv(iv,id,il)/nsim
c                  else
c                        write(lout2,105) il,disl,sumgam(iv,id,il)/nsim,
c     +                  int(sumnp(iv,id,il)/nsim),sumhm(iv,id,il)/nsim,
c     +                  sumtm(iv,id,il)/nsim
c                  endif
            write(lout2,'(i4,f16.5,i12,f16.5,i6,2f12.4,i4)') 
     +                       (id-1)*nvarg+iv,
     +                       disl,
     +                       int(sumnp(iv,id,il)/nsim),
     +                       sumgam(iv,id,il)/nsim,
     +                       (id-1)*nvarg+iv,azm,dip,it
 105              format(1x,i3,1x,f12.3,1x,f12.5,1x,i8,4(1x,f14.5))
            end do
      end do
      end do
      end if
c      close(lout)
      return
      end subroutine
 
 
      subroutine check(vrmin,vrmax)
c-----------------------------------------------------------------------
c
c                Error Check and Note Variogram Types
c                ************************************
c
c Go through each variogram type and note the type to the screen and
c report any possible errors.
c
c
c
c
c
c-----------------------------------------------------------------------
      use       geostat
      real      vrmin(*),vrmax(*)
      character title*100
c
c Loop over all the variograms to be computed:
c
      do iv=1,nvarg
c
c Note the variogram type and the variables being used:
c
c
      it = abs(ivtype(iv))
      if(it.eq. 1) title(1:24) = 'Semivariogram          :'
      if(it.eq. 2) title(1:24) = 'Cross Semivariogram    :'
      if(it.eq. 3) title(1:24) = 'Covariance             :'
      if(it.eq. 4) title(1:24) = 'Correlogram            :'
      if(it.eq. 5) title(1:24) = 'General Relative       :'
      if(it.eq. 6) title(1:24) = 'Pairwise Relative      :'
      if(it.eq. 7) title(1:24) = 'Variogram of Logarithms:'
      if(it.eq. 8) title(1:24) = 'Semimadogram           :'
      if(it.eq. 9) title(1:24) = 'Indicator 1/2 Variogram:'
      if(it.eq.10) title(1:24) = 'Indicator 1/2 Variogram:'
      write(title(25:64),100) names(ivtail(iv)),names(ivhead(iv))
 100  format('  tail=',a12,' head=',a12)
      if(ireal.eq.1) 
     +write(*,101) iv,title(1:64)
      if(ireal.eq.1) 
     +write(*,*)
 101  format(' Variogram ',i2,1x,a64)
c
c Check for possible errors or inconsistencies:
c
      if(it.eq.2) then
            if(ivtail(iv).eq.ivhead(iv)) write(*,201)
 201        format('  WARNING: cross variogram with the same variable!')
      else if(it.eq.5) then
            if(ivtail(iv).ne.ivhead(iv)) write(*,501)
            if(vrmin(ivtail(iv)).lt.0.0.and.vrmax(ivtail(iv)).gt.0.0)
     +            write(*,502)
            if(vrmin(ivhead(iv)).lt.0.0.and.vrmax(ivhead(iv)).gt.0.0)
     +            write(*,502)
 501        format('  WARNING: cross general relative variogram are',
     +             ' difficult to interpret!')
 502        format('  WARNING: there are both positive and negative',
     +             ' values - lag mean could be zero!')
      else if(it.eq.6) then
            if(ivtail(iv).ne.ivhead(iv)) write(*,601)
            if(vrmin(ivtail(iv)).lt.0.0.and.vrmax(ivtail(iv)).gt.0.0)
     +            write(*,602)
            if(vrmin(ivhead(iv)).lt.0.0.and.vrmax(ivhead(iv)).gt.0.0)
     +            write(*,602)
 601        format('  WARNING: cross pairwise relative variogram are',
     +             ' difficult to interpret!')
 602        format('  WARNING: there are both positive and negative',
     +             ' values - pair means could be zero!')
      else if(it.eq.7) then
            if(ivtail(iv).ne.ivhead(iv)) write(*,701)
            if(vrmin(ivtail(iv)).lt.0.0.or.vrmin(ivhead(iv)).lt.0.0)
     +      write(*,702)
 701        format('  WARNING: cross logarithmic variograms may be',
     +             ' difficult to interpret!')
 702        format('  WARNING: there are zero or negative',
     +             ' values - logarithm undefined!')
      else if(it.eq.8) then
            if(ivtail(iv).ne.ivhead(iv)) write(*,901)
 901        format('  WARNING: cross madograms may be difficult to',
     +             ' interpret!')
      endif
c
c Loop over all variograms:
c
      end do
      return
      end subroutine



      subroutine makepar
c-----------------------------------------------------------------------
c
c                      Write a Parameter File
c                      **********************
c
c
c
c-----------------------------------------------------------------------
      lun = 99
      open(lun,file='varsim.par',
     +     status='UNKNOWN')
      write(lun,10)
 10   format('                 Parameters for VarSim',/,
     +       '                 *********************',/,/,
     +       'START OF PARAMETERS:')
         write(lun,05)
 05   format('../data/lithology.dat        ',
     +       '-file with lithology information')
      write(lun,06)
 06   format('1   7                        ',
     +       '-   lithology column (0=not used), code')
      write(lun,11)
 11   format('../data/true.dat      ',
     +       '-file with data')
      write(lun,12)
 12   format('2   1   2             ',
     +       '-   number of variables, column numbers')
      write(lun,13)
 13   format('-1.0e21     1.0e21    ',
     +       '-   trimming limits')
      write(lun,14)
 14   format('varsim_reals.out      ',
     +       '-output file for variograms of realizations')
      write(lun,15)
 15   format('varsim_avg.out        ',
     +       '-output file for average variogram')
      write(lun,16)
 16   format('50   0.5   1.0        ',
     +       '-nx, xmn, xsiz')
      write(lun,17)
 17   format('50   0.5   1.0        ',
     +       '-ny, ymn, ysiz')
      write(lun,18)
 18   format(' 1   0.5   1.0        ',
     +       '-nz, zmn, zsiz')
      write(lun,180)
 180  format(' 100                  ',
     +       '-number of realizations')
      write(lun,19)
 19   format('2  10                 ',
     +       '-number of directions, number of lags')
      write(lun,20)
 20   format(' 1  0  0              ',
     +       '-ixd(1),iyd(1),izd(1)')
      write(lun,21)
 21   format(' 0  1  0              ',
     +       '-ixd(2),iyd(2),izd(2)')
      write(lun,22)
 22   format('1                     ',
     +       '-standardize sill? (0=no, 1=yes)')
      write(lun,23)
 23   format('5                     ',
     +       '-number of variograms')
      write(lun,24)
 24   format('1   1   1             ',
     +       '-tail variable, head variable, variogram type')
      write(lun,25)
 25   format('1   1   3             ',
     +       '-tail variable, head variable, variogram type')
      write(lun,26)
 26   format('2   2   1             ',
     +       '-tail variable, head variable, variogram type')
      write(lun,27)
 27   format('2   2   3             ',
     +       '-tail variable, head variable, variogram type')
      write(lun,28)
 28   format('1   1   9  2.5        ',
     +       '-tail variable, head variable, variogram type')
      write(lun,40)
 40   format(//,'type 1 = traditional semivariogram',/,
     +          '     2 = traditional cross semivariogram',/,
     +          '     3 = covariance',/,
     +          '     4 = correlogram',/,
     +          '     5 = general relative semivariogram',/,
     +          '     6 = pairwise relative semivariogram',/,
     +          '     7 = semivariogram of logarithms',/,
     +          '     8 = semimadogram',/,
     +          '     9 = indicator semivariogram - continuous',/,
     +          '     10= indicator semivariogram - categorical')
      write(lun,41)
 41   format(//,'VarSim is a conversion of gamsim_ave to support',/,
     +         'the file format used by varcalc, varmodel and varplot.')
      close(lun)
      return
      end subroutine



      subroutine chknam(str,len)
c-----------------------------------------------------------------------
c
c                   Check for a Valid File Name
c                   ***************************
c
c This subroutine takes the character string "str" of length "len" and
c removes all leading blanks and blanks out all characters after the
c first blank found in the string (leading blanks are removed first).
c
c
c
c-----------------------------------------------------------------------
      character(len=*), intent(inout) :: str
      integer itrim
c
c Remove leading blanks:
      str=adjustl(str)
c
c find first two blanks and blank out remaining characters:
      itrim=index(str,'   ')
      if (itrim > 0) str(itrim:)=' '
c
c Look for "-fi"
      itrim=index(str,'-fi')
      if (itrim > 0) str(itrim:)=' '
c
c Look for "\fi"
      itrim=index(str,'\fi')
      if (itrim > 0) str(itrim:)=' '
c
c Return with modified file name:
      return
      end subroutine
