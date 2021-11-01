c public
c These routines for covabar are essentially the exact routine found in
c the original GSLIB code for calculating covariances. Intents for f2py
c were added as comments, but the routines were otherwise left unchanged
c except for an increase in the number of digits for pi

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

      subroutine cova3(x1,y1,z1,x2,y2,z2,ivarg,nst,MAXNST,c0,it,cc,aa,
     +                 irot,MAXROT,rotmat,cmax,cova)
c-----------------------------------------------------------------------
c
c                    Covariance Between Two Points
c                    *****************************
c
c This subroutine calculated the covariance associated with a variogram
c model specified by a nugget effect and nested varigoram structures.
c The anisotropy definition can be different for each nested structure.
c
c
c
c INPUT VARIABLES:
c
c   x1,y1,z1         coordinates of first point
c   x2,y2,z2         coordinates of second point
c   nst(ivarg)       number of nested structures (maximum of 4)
c   ivarg            variogram number (set to 1 unless doing cokriging
c                       or indicator kriging)
c   MAXNST           size of variogram parameter arrays
c   c0(ivarg)        isotropic nugget constant
c   it(i)            type of each nested structure:
c                      1. spherical model of range a;
c                      2. exponential model of parameter a;
c                           i.e. practical range is 3a
c                      3. gaussian model of parameter a;
c                           i.e. practical range is a*sqrt(3)
c                      4. power model of power a (a must be gt. 0  and
c                           lt. 2).  if linear model, a=1,c=slope.
c   cc(i)            multiplicative factor of each nested structure.
c                      (sill-c0) for spherical, exponential,and gaussian
c                      slope for linear model.
c   aa(i)            parameter "a" of each nested structure.
c   irot             index of the rotation matrix for the first nested
c                    structure (the second nested structure will use
c                    irot+1, the third irot+2, and so on)
c   MAXROT           size of rotation matrix arrays
c   rotmat           rotation matrices
c
c
c OUTPUT VARIABLES:
c
c   cmax             maximum covariance
c   cova             covariance between (x1,y1,z1) and (x2,y2,z2)
c
c
c
c EXTERNAL REFERENCES: sqdist    computes anisotropic squared distance
c                      rotmat    computes rotation matrix for distance
c NOTE : RM - 2019 - pygeostat cova3 modified for full double precision
c-----------------------------------------------------------------------
      parameter(PI=3.14159265359,PMX=1.e10,EPSLON=1.e-10)
      integer   nst,it(*),ivarg,MAXNST
      real*8    c0,cc(*),aa(*),x1,y1,z1,x2,y2,z2,cmax,cova
      real*8    rotmat(MAXROT,3,3),hsqd,sqdist
cf2py intent(in) x1,y1,z1,x2,y2,z2,ivarg,nst,MAXNST,c0,it,cc,aa
cf2py intent(in) irot,MAXROT,rotmat,cmax
cf2py intent(out) cova

c
c Calculate the maximum covariance value (used for zero distances and
c for power model covariance):
c
      istart = 1 + (ivarg-1)*MAXNST
      cmax   = c0
      do is=1,nst
            ist = istart + is - 1
            if(it(ist).eq.4) then
                  cmax = cmax + PMX
            else
                  cmax = cmax + cc(ist)
            endif
      end do
c
c Check for "zero" distance, return with cmax if so:
c
      hsqd = sqdist(x1,y1,z1,x2,y2,z2,irot,MAXROT,rotmat)
      if(real(hsqd).lt.EPSLON) then
            cova = cmax
            return
      endif
c
c Loop over all the structures:
c
      cova = 0.0
      do is=1,nst
            ist = istart + is - 1
c
c Compute the appropriate distance:
c
            if(ist.ne.1) then
                  ir = min((irot+is-1),MAXROT)
                  hsqd=sqdist(x1,y1,z1,x2,y2,z2,ir,MAXROT,rotmat)
            end if
            h = real(dsqrt(hsqd))
c
c Spherical Variogram Model?
c
            if(it(ist).eq.1) then
                  hr = h/aa(ist)
                  if(hr.lt.1.) cova=cova+cc(ist)*(1.-hr*(1.5-.5*hr*hr))
c
c Exponential Variogram Model?
c
            else if(it(ist).eq.2) then
                  cova = cova + cc(ist)*exp(-3.0*h/aa(ist))
c
c Gaussian Variogram Model?
c
            else if(it(ist).eq.3) then
                  cova = cova + cc(ist)*exp(-(3.0*h/aa(ist))
     +                                      *(h/aa(ist)))
c
c Power Variogram Model?
c
            else if(it(ist).eq.4) then
                  cova = cova + cmax - cc(ist)*(h**aa(ist))
c
c Hole Effect Model?
c
            else if(it(ist).eq.5) then
c                 d = 10.0 * aa(ist)
c                 cova = cova + cc(ist)*exp(-3.0*h/d)*cos(h/aa(ist)*PI)
                  cova = cova + cc(ist)*cos(h/aa(ist)*PI)
            endif
      end do
c
c Finished:
c
      return
      end


      subroutine pwcova3(locations,ivarg,nst,MAXNST,c0,it,cc,aa,
     +                   irot,MAXROT,rotmat,cova,ndata,ndim)
c-----------------------------------------------------------------------
c
c                Pairwise Covariance between locations
c                *************************************
c
c compute the pairwise covariance between locations in the passed array
c
c Ryan Martin March 2019
c-----------------------------------------------------------------------
      real*8, intent(in) :: locations(ndim, ndata),
     +                      c0,cc(nst),aa(nst),
     +                      rotmat(MAXROT,nst,nst)
      integer, intent(in) :: ivarg,nst,MAXNST,it(nst),
     +                       irot,MAXROT,ndata,ndim
      real*8, intent(out) :: cova(ndata, ndata)
      real*8 :: c, x1, y1, z1, x2, y2, z2, cmax
      integer :: i, j
c     fill the diagonal
      call cova3(0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,ivarg,
     +           nst,MAXNST,c0,it,cc,aa,irot,MAXROT,rotmat,cmax,c)
      do i = 1, ndata
            cova(i, i) = cmax
      enddo
c     fill the rest of the matrix
      do i = 1, ndata
            x1 = locations(1, i)
            y1 = locations(2, i)
            z1 = locations(3, i)
            do j = i + 1, ndata
                  x2 = locations(1, j)
                  y2 = locations(2, j)
                  z2 = locations(3, j)
                  call cova3(x1,y1,z1,x2,y2,z2,ivarg,nst,MAXNST,c0,it,
     +                       cc,aa,irot,MAXROT,rotmat,cmax,c)
                  cova(i, j) = c
                  cova(j, i) = c
            enddo
      enddo
      end subroutine pwcova3



      real*8 function sqdist(x1,y1,z1,x2,y2,z2,ind,MAXROT,rotmat)
c-----------------------------------------------------------------------
c
c    Squared Anisotropic Distance Calculation Given Matrix Indicator
c    ***************************************************************
c
c This routine calculates the anisotropic distance between two points
c  given the coordinates of each point and a definition of the
c  anisotropy.
c
c
c INPUT VARIABLES:
c
c   x1,y1,z1         Coordinates of first point
c   x2,y2,z2         Coordinates of second point
c   ind              The rotation matrix to use
c   MAXROT           The maximum number of rotation matrices dimensioned
c   rotmat           The rotation matrices
c
c
c
c OUTPUT VARIABLES:
c
c   sqdist           The squared distance accounting for the anisotropy
c                      and the rotation of coordinates (if any).
c
c
c NO EXTERNAL REFERENCES
c
c
c Author: C. Deutsch                                Date: September 1989
c-----------------------------------------------------------------------
      integer ind,MAXROT
      real*8 rotmat(MAXROT,3,3),cont,dx,dy,dz
      real*8 x1,y1,z1,x2,y2,z2
c
c Compute component distance vectors and the squared distance:
c
      dx = dble(x1 - x2)
      dy = dble(y1 - y2)
      dz = dble(z1 - z2)
      sqdist = 0.0
      do i=1,3
            cont   = rotmat(ind,i,1) * dx
     +             + rotmat(ind,i,2) * dy
     +             + rotmat(ind,i,3) * dz
            sqdist = sqdist + cont * cont
      end do
      return
      end


      subroutine setrot(ang1,ang2,ang3,anis1,anis2,ind,MAXROT,rotmat)
c-----------------------------------------------------------------------
c
c              Sets up an Anisotropic Rotation Matrix
c              **************************************
c
c Sets up the matrix to transform cartesian coordinates to coordinates
c accounting for angles and anisotropy (see manual for a detailed
c definition):
c
c
c INPUT PARAMETERS:
c
c   ang1             Azimuth angle for principal direction
c   ang2             Dip angle for principal direction
c   ang3             Third rotation angle
c   anis1            First anisotropy ratio
c   anis2            Second anisotropy ratio
c   ind              matrix indicator to initialize
c   MAXROT           maximum number of rotation matrices dimensioned
c   rotmat           rotation matrices
c
c
c NO EXTERNAL REFERENCES
c
c
c Author: C. Deutsch                                Date: September 1989
c-----------------------------------------------------------------------
      parameter(DEG2RAD=3.14159265359/180.0,EPSLON=1.e-20)
      integer   ind,MAXROT
      real*8    rotmat(MAXROT,3,3),afac1,afac2,sina,sinb,sint,
     +          cosa,cosb,cost,ang1,ang2,ang3,anis1,anis2

cf2py intent(in) ang1,ang2,ang3,anis1,anis2,ind,MAXROT
cf2py intent(out) rotmat

c
c Converts the input angles to three angles which make more
c  mathematical sense:
c
c         alpha   angle between the major axis of anisotropy and the
c                 E-W axis. Note: Counter clockwise is positive.
c         beta    angle between major axis and the horizontal plane.
c                 (The dip of the ellipsoid measured positive down)
c         theta   Angle of rotation of minor axis about the major axis
c                 of the ellipsoid.
c
      if(ang1.ge.0.0.and.ang1.lt.270.0) then
            alpha = (90.0   - ang1) * DEG2RAD
      else
            alpha = (450.0  - ang1) * DEG2RAD
      endif
      beta  = -1.0 * ang2 * DEG2RAD
      theta =        ang3 * DEG2RAD
c
c Get the required sines and cosines:
c
      sina  = dble(sin(alpha))
      sinb  = dble(sin(beta))
      sint  = dble(sin(theta))
      cosa  = dble(cos(alpha))
      cosb  = dble(cos(beta))
      cost  = dble(cos(theta))
c
c Construct the rotation matrix in the required memory:
c
      afac1 = 1.0 / dble(max(anis1,EPSLON))
      afac2 = 1.0 / dble(max(anis2,EPSLON))
      rotmat(ind,1,1) =       (cosb * cosa)
      rotmat(ind,1,2) =       (cosb * sina)
      rotmat(ind,1,3) =       (-sinb)
      rotmat(ind,2,1) = afac1*(-cost*sina + sint*sinb*cosa)
      rotmat(ind,2,2) = afac1*(cost*cosa + sint*sinb*sina)
      rotmat(ind,2,3) = afac1*( sint * cosb)
      rotmat(ind,3,1) = afac2*(sint*sina + cost*sinb*cosa)
      rotmat(ind,3,2) = afac2*(-sint*cosa + cost*sinb*sina)
      rotmat(ind,3,3) = afac2*(cost * cosb)
c
c Return to calling program:
c
      return
      end
