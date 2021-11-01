c These routines for DGM are essentially the exact routine found in
c the original histscale code for calculating the discrete Gaussian model.
c Intents for f2py were added as comments, but the routines were otherwise
c left unchanged.
c
c Compile with f2py -c -m --fcompiler=gnu95 dgm dgm.for
c
c These subroutines are used under the license distributed with GSLIB:
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C                                                                      %
C Copyright (C) 2003, Statios Software and Services Incorporated.  All %
C rights reserved.                                                     %
C                                                                      %
C This program has been modified from the one distributed in 1996 (see %
C below).  This version is also distributed in the hope that it will   %
C be useful, but WITHOUT ANY WARRANTY. Compiled programs based on this %
C code may be redistributed without restriction; however, this code is %
C for one developer only. Each developer or user of this source code   %
C must purchase a separate copy from Statios.                          %
C                                                                      %
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C                                                                      %
C Copyright (C) 1996, The Board of Trustees of the Leland Stanford    %
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

      subroutine dgm(nd,vr1,vr2,coef,xm,xv,zcerr,var,accerr,np,fi,zc,h)
c-----------------------------------------------------------------------
c  
c                       Discrete Gaussian Model 
c                       ***********************
c
c-----------------------------------------------------------------------   
c np is the number of Hermite polynomials
c accerr is the tolerance level
      parameter(PI=3.1415926535,EXP=2.718281828)
      real    accerr
      integer np,nd
      real    vr1(nd),vr2(nd),L,coef,xm,xv,zcerr,var
      real*8  cump
      integer test,p
c  Extra arrays
c      real, allocatable :: y(:,:),g(:),zc(:),h(:,:),fi(:)

      real y(nd,3),g(nd),zc(nd),h(np+2,nd),fi(np+2)

cf2py intent(in) nd,vr1,xm,xv,accerr,np,coef
cf2py intent(out) vr2,zcerr,var,fi,zc,h

c      allocate (y(nd,3),stat= test)
c      if(test.ne.0) then
c            write(*,*) 'Error: Allocation failed!', test
c            go to 321
c      end if
c      allocate (g(nd),stat= test)
c      if(test.ne.0) then
c            write(*,*) 'Error: Allocation failed!', test
c            go to 321
c      end if
c      allocate (zc(nd),stat= test)
c      if(test.ne.0) then
c            write(*,*) 'Error: Allocation failed!', test
c            go to 321
c      end if
c      allocate (h(np+2,nd),stat= test)
c      if(test.ne.0) then
c            write(*,*) 'Error: Allocation failed!', test
c            go to 321
c      end if
c      allocate (fi(np+2),stat= test)
c      if(test.ne.0) then
c            write(*,*) 'Error: Allocation failed!', test
c            go to 321
c      end if
    
c Get the distribution and normal score transform:
      do i=1,nd
            cump  = dble(real(i) / real(nd+1))
            call gauinv(cump,vrrg,ierr)
            y(i,1)  = vrrg
      end do

c The Hermite polynomials:
      do p=1,np+1
            fi(p) = 0.
      end do
      do i=1,nd
            H(1,i)=1
            H(2,i)=-y(i,1)
      end do
      do p=1,np-1
            do i=1,nd
                  H(p+2,i)=-y(i,1)*H(p+1,i)/(sqrt(real(p+1)))-
     +                     sqrt(real(p))/sqrt(real(p)+1)*H(p,i)
            end do
      end do

c The probability for the transformed values:
      do i=1,nd
            g(i)=(1/sqrt(2*PI))*(EXP**(-(y(i,1)**2)/2))
      end do

c Calculate the Hermite coefficient depending on the order np
      do i=1,nd
            fi(1)=fi(1)+vr1(i)/nd
      end do
      do p=2,np+1
            do i=2,nd
                  fi(p)=fi(p)+(vr1(i-1)-vr1(i))*
     +                   (1/sqrt(real(p-1)))*H(p-1,i)*g(i)
            end do
      end do

c Check the anamorphosis by calculating the point scale values and 
c the variance - not doing anything with these values, but we could:
      zc    = 0.0
      zcerr = 0.0
      do i=1,nd
            do p=1,np+1
                  zc(i)=zc(i)+fi(p)*H(p,i)
            end do
            zcerr = zcerr + (vr1(i)-zc(i))**2/nd
      end do
      var = 0.0
      do p=2,np+1
            var=var+fi(p)**2
      end do
      var = var - xv
    
c Calculate the change of support coefficient such that it minimizes 
c the difference between the sum of fi(p)^2*r^(2*p) and the variance 
c for block support obtained using gammabar
      varblock = xv*coef
      a   = 0.0
      b   = 1.0
97    L=b-a
      rm=(a+b)/2
      r1=a+(L/4)
      r2=b-(L/4)
      varr1=0
      varr2=0
      varrm=0
      do p=2,np+1
            varr1=varr1+(fi(p)**2)*(r1**(2*(p-1)))  
            varr2=varr2+(fi(p)**2)*(r2**(2*(p-1)))  
            varrm=varrm+(fi(p)**2)*(rm**(2*(p-1)))  
      end do
      difr1=abs(varr1-varblock)
      difr2=abs(varr2-varblock)
      difrm=abs(varrm-varblock)
      if(difrm.lt.accerr) then
            goto 98
      else
      if(difr1.lt.difrm) then
            b=rm
            goto 97
      else
      if(difr1.ge.difrm) then
          if(difr2.lt.difrm) then
                a=rm
                goto 97
          else
          if(difr2.ge.difrm) then
                a=r1
                b=r2
                goto 97
          else 
                write(*,*) 'ERROR'
          end if
          end if
      end if
      end if
      end if
98    continue
	  write(*,*) 'rm =',rm
    
c The anamorphosis at block support:
      do i=1,nd
            vr2(i) = 0.0
      end do
      do i=1,nd
            do p=1,np+1
                  vr2(i)=vr2(i)+fi(p)*H(p,i)*rm**(p-1)
            end do
      end do
  
c Finished:
321   continue
c    deallocate(y,stat=test)
c    deallocate(g,stat=test)
c    deallocate(zc,stat=test)
c    deallocate(h,stat=test)
c    deallocate(fi,stat=test)
      return
      end subroutine dgm

      subroutine gauinv(p,xp,ierr)
c-----------------------------------------------------------------------
c
c Computes the inverse of the standard normal cumulative distribution
c f unction with a numerical approximation from : Statistical Computing,
c by W.J. Kennedy, Jr. and James E. Gentle, 1980, p. 95.
c
c INPUT/OUTPUT:
c
c   p      = double precision cumulative probability value: dble(psingle)
c   xp   = G^-1 (p) in single precision
c   ierr = 1 - then error situation (p out of range), 0 - OK
c
c-----------------------------------------------------------------------
      real*8 p0,p1,p2,p3,p4,q0,q1,q2,q3,q4,y,pp,lim,p
      save   p0,p1,p2,p3,p4,q0,q1,q2,q3,q4,lim

cf2py intent(in) p,ierr
cf2py intent(out) xp

      ! Coefficients of approximation:
      data lim/1.0e-10/
      data p0/-0.322232431088/,p1/-1.0/,p2/-0.342242088547/,
     + p3/-0.0204231210245/,p4/-0.0000453642210148/
      data q0/0.0993484626060/,q1/0.588581570495/,q2/0.531103462366/,
     + q3/0.103537752850/,q4/0.0038560700634/

      ! Check for an error situation:
      ierr = 1
      if(p.lt.lim) then
              xp = -1.0e10
              return
      end if
      if(p.gt.(1.0-lim)) then
              xp =  1.0e10
              return
      end if
      ierr = 0        

      ! Get k for an error situation:
      pp   = p
      if(p.gt.0.5) pp = 1 - pp
      xp   = 0.0
      if(p.eq.0.5) return

      ! Approximate the function:
      y  = dsqrt(dlog(1.0/(pp*pp)))
      xp = real( y + ((((y*p4+p3)*y+p2)*y+p1)*y+p0) / 
     + ((((y*q4+q3)*y+q2)*y+q1)*y+q0) )
      if(real(p).eq.real(pp)) xp = -xp
      
      ! Return with G^-1(p):
      return
      end subroutine gauinv
