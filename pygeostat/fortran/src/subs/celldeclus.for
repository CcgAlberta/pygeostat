       program main
c-----------------------------------------------------------------------
c
c           CellDeclus: An Updated Cell Declustering Program
c           ************************************************
c
c - saves the grid index, but not a grid of indices
c - randomizes the cell grid origin (not systematic)
c
c
c
c AUTHOR: Clayton V. Deutsch                             DATE: 1989-1999
c UPDATED:                                                          2015
c-----------------------------------------------------------------------
      parameter (MAXLEN=512,VERSION=4.000,
     +           KORDEI=12,MAXOP1=KORDEI+1,MAXINT=2**30)

      real      var(512)
      integer   test,nocc
      character datafl*512,sumfl*512,outfl*512,str*512,strlin*512
      logical   testfl,inflag
      real, allocatable    :: x(:),y(:),z(:),wt(:),vr(:),wtopt(:)
      integer, allocatable :: cellid(:),numincell(:),dincell(:)
c
c For random number generator:
c
      real*8  acorni
      common /iaco/ ixv(MAXOP1)
      data   ixv/MAXOP1*0.0/

      lin  = 1
      lout = 2
      lsum = 3
c
c Note VERSION number:
c
      write(*,9999) VERSION
 9999 format(/' CELLDECLUS Version: ',f5.3/)
c
c Get the name of the parameter file - try the default name if no input:
c
      do i=1,512
            str(i:i) = ' '
      end do
      call getarg(1,str)
      if(str(1:1).eq.' ')then
            write(*,*) 'Which parameter file do you want to use?'
            read (*,'(a)') str
      end if
      if(str(1:1).eq.' ') str(1:20) = 'CellDeclus.par      '
      inquire(file=str,exist=testfl)
      if(.not.testfl) then
            write(*,*) 'ERROR - the parameter file does not exist,'
            write(*,*) '        check for the file and try again  '
            write(*,*)
            if(str(1:20).eq.'CellDeclus.par      ') then
                  write(*,*) '        creating a blank parameter file'
                  call makepar
                  write(*,*)
            end if
            stop
      endif
      open(lin,file=str,status='OLD')
c
c Find Start of Parameters:
c
 1    read(lin,'(a4)',end=97) str(1:4)
      if(str(1:4).ne.'STAR') go to 1
c
c Read Input Parameters:
c
      read(lin,'(a512)',err=97) datafl
      call chknam(datafl,512)
      write(*,*) ' data file = ',datafl(1:40)

      read(lin,*,err=97) ix,iy,iz,ivr
      write(*,*) ' columns = ',ix,iy,iz,ivr

      read(lin,*,err=97) tmin,tmax
      write(*,*) ' tmin,tmax = ',tmin,tmax

      read(lin,'(a512)',err=97) sumfl
      call chknam(sumfl,512)
      write(*,*) ' summary file = ',sumfl(1:40)

      read(lin,'(a512)',err=97) outfl 
      call chknam(outfl,512)
      write(*,*) ' output file = ',outfl(1:40)

      read(lin,*,err=97) anisy,anisz
      write(*,*) ' anisotropy = ',anisy,anisz

      read(lin,*,err=97) ncell,cmin,cmax
      if(ncell.eq.1) cmax = cmin
      write(*,*) ' ncell min max = ',ncell,cmin,cmax

      read(lin,*,err=97) minmax,tarsize
      write(*,*) ' minmax flag = ',minmax,tarsize

      write(*,*)
      close(lin)

      norig  = 500
      ixv(1) = 69069
      do i=1,10000
            zz = real(acorni(idum))
      end do
c
c Make sure that we have a data file:
c
      inquire(file=datafl,exist=testfl)
      if(.not.testfl) then
            write(*,*) 'ERROR: ',datafl,' does not exist'
            stop
      endif
c
c Open up the input and output files:
c
      open(lin,file=datafl,status='OLD')
      open(lsum,file=sumfl,status='UNKNOWN')
      open(lout,file=outfl,status='UNKNOWN')
c
c Read the header off the data file, find MAXDAT, prepare output files:
c
      read(lin,*,err=98)
      read(lin,*,err=98) nvari
      do i=1,nvari
            read(lin,*,err=98)
      end do
      maxdat = 0
 20   read(lin,*,end=40,err=98)(var(j),j=1,nvari)
      if(var(ivr).lt.tmin.or.var(ivr).ge.tmax)go to 20
      maxdat = maxdat + 1
      go to 20
 40   continue
c
c Allocate the needed memory:
c
      allocate (x(maxdat),stat = test)
      if (test.ne.0) then
            write(*,*) 'Error: Allocation failed', test
            stop
      end if
      allocate (y(maxdat),stat = test)
      if (test.ne.0) then
            write(*,*) 'Error: Allocation failed', test
            stop
      end if
      allocate (z(maxdat),stat = test)
      if (test.ne.0) then
            write(*,*) 'Error: Allocation failed', test
            stop
      end if
      allocate (wt(maxdat),stat = test)
      if (test.ne.0) then
            write(*,*) 'Error: Allocation failed', test
            stop
      end if
      allocate (vr(maxdat),stat = test)
      if (test.ne.0) then
            write(*,*) 'Error: Allocation failed', test
            stop
      end if
      allocate (wtopt(maxdat),stat = test)
      if (test.ne.0) then
            write(*,*) 'Error: Allocation failed', test
            stop
      end if
      allocate (cellid(maxdat),stat = test)
      if (test.ne.0) then
            write(*,*) 'Error: Allocation failed', test
            stop
      end if
      allocate (numincell(maxdat),stat = test)
      if (test.ne.0) then
            write(*,*) 'Error: Allocation failed', test
            stop
      end if
      allocate (dincell(maxdat),stat = test)
      if (test.ne.0) then
            write(*,*) 'Error: Allocation failed', test
            stop
      end if
      rewind(lin)
      read(lin,'(a)') str
      write(lout,'(a)') trim(str)
      write(lsum,'(a)') trim(str)
      read(lin,*,err=98) nvari
      write(lout,'(i3)') nvari+1
      do i=1,nvari
            read(lin,'(a)') str
            write(lout,'(a)') trim(str)
      end do
      write(lout,101)
 101  format('Cell Declustering Weight')
      write(lsum,102)
 102  format('2',/,'Cell Size',/,'Declustered Mean')
c
c Now, read in the actual data:
c
      nt = 0
      nd = 0
      vrmin =  1.0e21
      vrmax = -1.0e21
      xmin  =  1.0e21
      ymin  =  1.0e21
      zmin  =  1.0e21
      xmax  = -1.0e21
      ymax  = -1.0e21
      zmax  = -1.0e21
      vrav  =  0.0
 3    read(lin,*,end=4,err=98) (var(j),j=1,nvari)
      if(var(ivr).lt.tmin.or.var(ivr).ge.tmax) then
            nt = nt + 1
            go to 3
      endif
c
c Found a valid data:
c
      nd     = nd + 1
      vr(nd) = var(ivr)
      vrav   = vrav + vr(nd)
      if(vr(nd).lt.vrmin) vrmin = vr(nd)
      if(vr(nd).gt.vrmax) vrmax = vr(nd)
      wtopt(nd) = 1.0
      x(nd)     = 0.0
      y(nd)     = 0.0
      z(nd)     = 0.0
      if(ix.ge.1) x(nd) = var(ix)
      if(iy.ge.1) y(nd) = var(iy)
      if(iz.ge.1) z(nd) = var(iz)
      if(x(nd).lt.xmin) xmin=x(nd)
      if(x(nd).gt.xmax) xmax=x(nd)
      if(y(nd).lt.ymin) ymin=y(nd)
      if(y(nd).gt.ymax) ymax=y(nd)
      if(z(nd).lt.zmin) zmin=z(nd)
      if(z(nd).gt.zmax) zmax=z(nd)
      go to 3
 4    close(lin)
      if(nd.le.1) then
            write(*,*) ' ERROR: too few data'
            stop
      endif
      vrav = vrav / real(nd)
c
c Write out some statistics:
c
      write(*,900) nd,vrav,vrmin,vrmax,(xmax-xmin),(ymax-ymin),
     +             (zmax-zmin)
 900  format(/' There are ',i8,' data with:',/,
     +        '   mean value            = ',g15.8,/,
     +        '   minimum and maximum   = ',2g15.8,/,
     +        '   size of data vol in X = ',g15.8,/,
     +        '   size of data vol in Y = ',g15.8,/,
     +        '   size of data vol in Z = ',g15.8,/)
c
c Make sure cmin is not too small:
c
      if(cmin.lt.(xmax-xmin)*0.001) cmin = (xmax-xmin)*0.001
      if(cmin.lt.(ymax-ymin)*0.001) cmin = (ymax-ymin)*0.001
      if(cmin.lt.(zmax-zmin)*0.001) cmin = (zmax-zmin)*0.001
c
c Initialize the "best" weight values:
c
      vrop = vrav
      best = 0.0
      write(lsum,300)best,vrop
 300  format(g15.8,2x,g15.8)
      diffcs = 1.0e21
c
c Define a "lower" origin to use for the cell sizes:
c
      xo1 = xmin - 0.01*(xmax-xmin)
      yo1 = ymin - 0.01*(ymax-ymin)
      zo1 = zmin - 0.01*(zmax-zmin)
c
c Define some cell size parameters:
c
      xinc = (cmax-cmin) / real(ncell)
      yinc = anisy * xinc
      zinc = anisz * xinc
      xcs  =  cmin        - xinc
      ycs  = (cmin*anisy) - yinc
      zcs  = (cmin*anisz) - zinc
c
c MAIN LOOP over cell sizes:
c
      do lp=1,ncell+1
            xcs = xcs + xinc
            ycs = ycs + yinc
            zcs = zcs + zinc
            wt  = 0.0
c
c Determine the maximum number of grid cells in the network:
c
            ncellx = int((xmax-(xo1-xcs))/xcs)+2
            ncelly = int((ymax-(yo1-ycs))/ycs)+2
            ncellz = int((zmax-(zo1-zcs))/zcs)+2
c
c Loop over all the origin offsets selected:
c
            do kp=1,norig
                  xo = xo1 - real(acorni(idum))*xcs
                  yo = yo1 - real(acorni(idum))*ycs
                  zo = zo1 - real(acorni(idum))*zcs
c
c Initialize the cumulative weight indicators:
c
                  nocc = 0
                  cellid = 0
                  numincell = 0
c
c Determine which cell each datum is in:
c
                  do i=1,nd
                        call getindx(ncellx,xo,xcs,x(i),icellx,inflag)
                        if(.not.inflag) then
                              write(*,*) ncellx,xo,xcs,x(i),icellx
                              stop
                        end if      
                        call getindx(ncelly,yo,ycs,y(i),icelly,inflag)
                        if(.not.inflag) then
                              write(*,*) ncelly,yo,ycs,y(i),icelly
                              stop
                        end if      
                        call getindx(ncellz,zo,zcs,z(i),icellz,inflag)
                        if(.not.inflag) then
                              write(*,*) ncellz,zo,zcs,z(i),icellz
                              stop
                        end if      
                        icell =  icellx + 
     +                          (icelly-1)*ncellx +
     +                          (icellz-1)*ncelly*ncellx
                        if(nocc.ge.1) then
                              do iocc=1,nocc
                                    if(icell.eq.cellid(iocc)) then
                                          dincell(i) = iocc
                                          numincell(iocc) = 
     +                                    numincell(iocc)+1
                                          go to 10
                                    end if      
                              end do      
                        end if
                        nocc = nocc + 1
                        cellid(nocc)    = icell
                        dincell(i)      = nocc
                        numincell(nocc) = 1
 10                     continue
                  end do
c
c The weight assigned to each datum is inversely proportional to the
c number of data in the cell:
c
                  do i=1,nd
                        iocc = dincell(i)
                        wt(i) = wt(i) + 1.0/(nocc*numincell(iocc))
                  end do
c
c End loop over all random origins:
c
            end do
c
c Compute the weighted average for this cell size:
c
            sumw  = 0.0
            sumwg = 0.0
            do i=1,nd
                  sumw  = sumw  + wt(i)
                  sumwg = sumwg + wt(i)*vr(i)
            end do
            vrcr  = sumwg / sumw
            write(lsum,300) xcs,vrcr
c
c see if this weighting is optimal:
c
            if((minmax.eq.-1.and.vrcr.lt.vrop).or.
     +         (minmax.eq. 1.and.vrcr.gt.vrop).or.(ncell.eq.1)) then
                  best = xcs
                  vrop = vrcr
                  do i=1,nd
                        wtopt(i) = wt(i)
                  end do
            endif

            if(minmax.eq.0.and.abs(xcs-tarsize).lt.diffcs) then
                  diffcs = abs(xcs-tarsize)
                  vrop   = vrcr
                  do i=1,nd
                        wtopt(i) = wt(i)
                  end do
            end if
c
c END MAIN LOOP over all cell sizes:
c
      end do
      close(lsum)
c
c Get the optimal weights:
c
      sumw = 0.0
      do i=1,nd
            sumw = sumw + wtopt(i)
      end do
      wtmin = 99999.
      wtmax =-99999.
      facto = real(nd) / sumw
      do i = 1,nd
            wtopt(i) = wtopt(i) * facto
            if(wtopt(i).lt.wtmin) wtmin = wtopt(i)
            if(wtopt(i).gt.wtmax) wtmax = wtopt(i)
      end do
c
c Read the header off the data file:
c
      open(lin,file=datafl,status='OLD')
      read(lin,'()',err=98)
      read(lin,*,err=98) nvari
      do i=1,nvari
            read(lin,'()',err=98)
      end do
c
c Read through input file appending the declustering weight to the
c  output file:
c
      nd2 = 0
 5    read(lin,*,end=6) (var(i),i=1,nvari)
      if(var(ivr).lt.tmin.or.var(ivr).ge.tmax) then
            thewt = -999.0
      else
            nd2   = nd2 + 1
            thewt = wtopt(nd2)
      end if
c
c Write out the results:
c
      backspace lin
      read(lin,'(a)') strlin
      call strlen(strlin,MAXLEN,lostr)
      write(lout,'(a,1x,g15.8)') strlin(1:lostr),thewt
      go to 5
 6    continue
      if(nd.ne.nd2) then
            write(*,*)
            write(*,*) 'ERROR in data somewhere - changed during run!'
            write(*,*)
      end if
c
c Some Debugging Information:
c
      write(*,901) vrop,wtmin,wtmax,1.0
 901  format('   declustered mean      = ',g15.8,/,
     +       '   min and max weight    = ',2g15.8,/,
     +       '   equal weighting       = ',g15.8,/)
c
c Finished:
c
      close(lin)
      close(lout)
      close(lsum)
      write(*,9998) VERSION
 9998 format(/' CELLDECLUS Version: ',f5.3, ' Finished'/)
      stop
 97   stop 'ERROR in parameter file'
 98   stop 'ERROR in data file'
      end



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
      open(lun,file='CellDeclus.par',status='UNKNOWN')
      write(lun,10)
 10   format('                  Parameters for CELLDECLUS',/,
     +       '                  *************************',/,/,
     +       'START OF PARAMETERS:')

      write(lun,11)
 11   format('data.dat                    ',
     +       '-file with data')
      write(lun,12)
 12   format('1   2   0   3               ',
     +       '-  columns for X, Y, Z, and variable')
      write(lun,13)
 13   format('-1.0e21     1.0e21          ',
     +       '-  trimming limits')
      write(lun,14)
 14   format('CellDeclusSumm.dat          ',
     +       '-file for summary output')
      write(lun,15)
 15   format('CellDeclus.dat              ',
     +       '-file for output with data and weights')
      write(lun,16)
 16   format('1.0   1.0                   ',
     +       '-Y and Z cell anisotropy (Ysize=size*Yanis)')
      write(lun,17)
 17   format('24  1.0  25.0               ',
     +       '-number of cell sizes, min size, max size')
      write(lun,18)
 18   format('-1  10.0                    ',
     +       '-cell size to keep: -1 = minimum')
      write(lun,19)
 19   format('                            ',
     +       '                     0 = specified')
      write(lun,20)
 20   format('                            ',
     +       '                    +1 = maximum')

      close(lun)
      return
      end



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
      end



      subroutine getindx(n,min,siz,loc,index,inflag)
c-----------------------------------------------------------------------
c
c     Gets the coordinate index location of a point within a grid
c     ***********************************************************
c
c
c n       number of "nodes" or "cells" in this coordinate direction
c min     origin at the center of the first cell
c siz     size of the cells
c loc     location of the point being considered
c index   output index within [1,n]
c inflag  true if the location is actually in the grid (false otherwise
c         e.g., if the location is outside then index will be set to
c         nearest boundary
c
c
c
c-----------------------------------------------------------------------
      integer   n,index
      real      min,siz,loc
      logical   inflag
c
c Compute the index of "loc":
c
      index = int( (loc-min)/siz + 1.5 )
c
c Check to see if in or out:
c
      if(index.lt.1) then
            index  = 1
            inflag = .false.
      else if(index.gt.n) then
            index  = n
            inflag = .false.
      else
            inflag = .true.
      end if
c
c Return to calling program:
c
      return
      end



      double precision function acorni(idum)
c-----------------------------------------------------------------------
c
c Fortran implementation of ACORN random number generator of order less
c than or equal to 12 (higher orders can be obtained by increasing the
c parameter value MAXORD).
c
c
c NOTES: 1. The variable idum is a dummy variable. The common block
c           IACO is used to transfer data into the function.
c
c        2. Before the first call to ACORN the common block IACO must
c           be initialised by the user, as follows. The values of
c           variables in the common block must not subsequently be
c           changed by the user.
c
c             KORDEI - order of generator required ( must be =< MAXORD)
c
c             MAXINT - modulus for generator, must be chosen small
c                      enough that 2*MAXINT does not overflow
c
c             ixv(1) - seed for random number generator
c                      require 0 < ixv(1) < MAXINT
c
c             (ixv(I+1),I=1,KORDEI)
c                    - KORDEI initial values for generator
c                      require 0 =< ixv(I+1) < MAXINT
c
c        3. After initialisation, each call to ACORN generates a single
c           random number between 0 and 1.
c
c        4. An example of suitable values for parameters is
c
c             KORDEI   = 10
c             MAXINT   = 2**30
c             ixv(1)   = an odd integer in the (approximate) range 
c                        (0.001 * MAXINT) to (0.999 * MAXINT)
c             ixv(I+1) = 0, I=1,KORDEI
c
c
c
c Author: R.S.Wikramaratna,                           Date: October 1990
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      parameter (KORDEI=12,MAXOP1=KORDEI+1,MAXINT=2**30)
      common/iaco/ ixv(MAXOP1)
      do i=1,KORDEI
            ixv(i+1)=(ixv(i+1)+ixv(i))
            if(ixv(i+1).ge.MAXINT) ixv(i+1)=ixv(i+1)-MAXINT
      end do
      acorni=dble(ixv(KORDEI+1))/MAXINT
      return
      end



      subroutine strlen(str,MAXLEN,lostr)
c-----------------------------------------------------------------------
c
c      Determine the length of the string minus trailing blanks
c
c
c
c-----------------------------------------------------------------------
      character(MAXLEN) :: str
      lostr = MAXLEN
      do i=1,MAXLEN
            j = MAXLEN - i + 1
            if(str(j:j).ne.' ') return
            lostr = lostr - 1
      end do
      return
      end
