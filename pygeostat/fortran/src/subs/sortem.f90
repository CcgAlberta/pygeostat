module sortem
  public :: modsortem,dblemodsortem
contains
      subroutine modsortem(a,n,iperm,b,c,d,e,f,g,h)
!-----------------------------------------------------------------------
!
! WARNING: this is a modified version of sortem - the array length
!          is specified - this is compatible with slicing
!
!                      Quickersort Subroutine
!                      **********************
!
! This is a subroutine for sorting a real array in ascending order. This
! is a Fortran translation of algorithm 271, quickersort, by R.S. Scowen
! in collected algorithms of the ACM.
!
! The method used is that of continually splitting the array into parts
! such that all elements of one part are less than all elements of the
! other, with a third part in the middle consisting of one element.  An
! element with value t is chosen arbitrarily (here we choose the middle
! element). i and j give the lower and upper limits of the segment being
! split.  After the split a value q will have been found such that 
! a(q)=t and a(l)<=t<=a(m) for all i<=l<q<m<=j.  The program then
! performs operations on the two segments (i,q-1) and (q+1,j) as follows
! The smaller segment is split and the position of the larger segment is
! stored in the lt and ut arrays.  If the segment to be split contains
! two or fewer elements, it is sorted and another segment is obtained
! from the lt and ut arrays.  When no more segments remain, the array
! is completely sorted.
!
!
! INPUT PARAMETERS:
!
!   ib,ie        start and end index of the array to be sorteda
!   a            array, a portion of which has to be sorted.
!   iperm        0 no other array is permuted.
!                1 array b is permuted according to array a
!                2 arrays b,c are permuted.
!                3 arrays b,c,d are permuted.
!                4 arrays b,c,d,e are permuted.
!                5 arrays b,c,d,e,f are permuted.
!                6 arrays b,c,d,e,f,g are permuted.
!                7 arrays b,c,d,e,f,g,h are permuted.
!               >7 no other array is permuted.
!
!   b,c,d,e,f,g,h  arrays to be permuted according to array a.
!
! OUTPUT PARAMETERS:
!
!    a      = the array, a portion of which has been sorted.
!
!    b,c,d,e,f,g,h  =arrays permuted according to array a (see iperm)
!
! NO EXTERNAL ROUTINES REQUIRED:
!
!-----------------------------------------------------------------------
      !dimension a(*),b(*),c(*),d(*),e(*),f(*),g(*),h(*)
      integer, intent(in) :: n,iperm
      real, dimension(n), intent(inout) :: a
      real, dimension(n), optional, intent(inout) :: b,c,d,e,f,g,h
!
! The dimensions for lt and ut have to be at least log (base 2) n
!
      integer   lt(64),ut(64),i,j,k,m,p,q,ib,ie

      ib = 1
      ie = n
!
! Initialize:
!
      j     = ie
      m     = 1
      i     = ib
      iring = iperm+1
      if (iperm.gt.7) iring=1
!
! If this segment has more than two elements  we split it
!
 10   if (j-i-1) 100,90,15
!
! p is the position of an arbitrary element in the segment we choose the
! middle element. Under certain circumstances it may be advantageous
! to choose p at random.
!
 15   p    = (j+i)/2
      ta   = a(p)
      a(p) = a(i)
      go to (21,19,18,17,16,161,162,163),iring
 163     th   = h(p)
         h(p) = h(i)
 162     tg   = g(p)
         g(p) = g(i)
 161     tf   = f(p)
         f(p) = f(i)
 16      te   = e(p)
         e(p) = e(i)
 17      td   = d(p)
         d(p) = d(i)
 18      tc   = c(p)
         c(p) = c(i)
 19      tb   = b(p)
         b(p) = b(i)
 21   continue
!
! Start at the beginning of the segment, search for k such that a(k)>t
!
      q = j
      k = i
 20   k = k+1
      if(k.gt.q)     go to 60
      if(a(k).le.ta) go to 20
!
! Such an element has now been found now search for a q such that a(q)<t
! starting at the end of the segment.
!
 30   continue
      if(a(q).lt.ta) go to 40
      q = q-1
      if(q.gt.k)     go to 30
      go to 50
!
! a(q) has now been found. we interchange a(q) and a(k)
!
 40   xa   = a(k)
      a(k) = a(q)
      a(q) = xa
      go to (45,44,43,42,41,411,412,413),iring
 413     xh   = h(k)
         h(k) = h(q)
         h(q) = xh
 412     xg   = g(k)
         g(k) = g(q)
         g(q) = xg
 411     xf   = f(k)
         f(k) = f(q)
         f(q) = xf
 41      xe   = e(k)
         e(k) = e(q)
         e(q) = xe
 42      xd   = d(k)
         d(k) = d(q)
         d(q) = xd
 43      xc   = c(k)
         c(k) = c(q)
         c(q) = xc
 44      xb   = b(k)
         b(k) = b(q)
         b(q) = xb
 45   continue
!
! Update q and search for another pair to interchange:
!
      q = q-1
      go to 20
 50   q = k-1
 60   continue
!
! The upwards search has now met the downwards search:
!
      a(i)=a(q)
      a(q)=ta
      go to (65,64,63,62,61,611,612,613),iring
 613     h(i) = h(q)
         h(q) = th
 612     g(i) = g(q)
         g(q) = tg
 611     f(i) = f(q)
         f(q) = tf
 61      e(i) = e(q)
         e(q) = te
 62      d(i) = d(q)
         d(q) = td
 63      c(i) = c(q)
         c(q) = tc
 64      b(i) = b(q)
         b(q) = tb
 65   continue
!
! The segment is now divided in three parts: (i,q-1),(q),(q+1,j)
! store the position of the largest segment in lt and ut
!
      if (2*q.le.i+j) go to 70
      lt(m) = i
      ut(m) = q-1
      i = q+1
      go to 80
 70   lt(m) = q+1
      ut(m) = j
      j = q-1
!
! Update m and split the new smaller segment
!
 80   m = m+1
      go to 10
!
! We arrive here if the segment has  two elements we test to see if
! the segment is properly ordered if not, we perform an interchange
!
 90   continue
      if (a(i).le.a(j)) go to 100
      xa=a(i)
      a(i)=a(j)
      a(j)=xa
      go to (95,94,93,92,91,911,912,913),iring
 913     xh   = h(i)
         h(i) = h(j)
         h(j) = xh
 912     xg   = g(i)
         g(i) = g(j)
         g(j) = xg
 911     xf   = f(i)
         f(i) = f(j)
         f(j) = xf
   91    xe   = e(i)
         e(i) = e(j)
         e(j) = xe
   92    xd   = d(i)
         d(i) = d(j)
         d(j) = xd
   93    xc   = c(i)
         c(i) = c(j)
         c(j) = xc
   94    xb   = b(i)
         b(i) = b(j)
         b(j) = xb
   95 continue
!
! If lt and ut contain more segments to be sorted repeat process:
!
 100  m = m-1
      if (m.le.0) go to 110
      i = lt(m)
      j = ut(m)
      go to 10
 110  continue
      return
      end subroutine modsortem

     subroutine dblemodsortem(a,n,iperm,b,c,d,e,f,g,h)
!-----------------------------------------------------------------------
!
! WARNING: this is a modified version of sortem - the array length
!          is specified - this is compatible with slicing
!
!                      Quickersort Subroutine
!                      **********************
!
! This is a subroutine for sorting a real array in ascending order. This
! is a Fortran translation of algorithm 271, quickersort, by R.S. Scowen
! in collected algorithms of the ACM.
!
! The method used is that of continually splitting the array into parts
! such that all elements of one part are less than all elements of the
! other, with a third part in the middle consisting of one element.  An
! element with value t is chosen arbitrarily (here we choose the middle
! element). i and j give the lower and upper limits of the segment being
! split.  After the split a value q will have been found such that 
! a(q)=t and a(l)<=t<=a(m) for all i<=l<q<m<=j.  The program then
! performs operations on the two segments (i,q-1) and (q+1,j) as follows
! The smaller segment is split and the position of the larger segment is
! stored in the lt and ut arrays.  If the segment to be split contains
! two or fewer elements, it is sorted and another segment is obtained
! from the lt and ut arrays.  When no more segments remain, the array
! is completely sorted.
!
!
! INPUT PARAMETERS:
!
!   ib,ie        start and end index of the array to be sorteda
!   a            array, a portion of which has to be sorted.
!   iperm        0 no other array is permuted.
!                1 array b is permuted according to array a
!                2 arrays b,c are permuted.
!                3 arrays b,c,d are permuted.
!                4 arrays b,c,d,e are permuted.
!                5 arrays b,c,d,e,f are permuted.
!                6 arrays b,c,d,e,f,g are permuted.
!                7 arrays b,c,d,e,f,g,h are permuted.
!               >7 no other array is permuted.
!
!   b,c,d,e,f,g,h  arrays to be permuted according to array a.
!
! OUTPUT PARAMETERS:
!
!    a      = the array, a portion of which has been sorted.
!
!    b,c,d,e,f,g,h  =arrays permuted according to array a (see iperm)
!
! NO EXTERNAL ROUTINES REQUIRED:
!
!-----------------------------------------------------------------------
      !dimension a(*),b(*),c(*),d(*),e(*),f(*),g(*),h(*)
      integer, intent(in) :: n,iperm
      real(kind=8), dimension(n), intent(inout) :: a
      real(kind=8), dimension(n), optional, intent(inout) :: b,c,d,e,f,g,h
      real(kind=8) th,tg,tf,te,td,tc,tb,ta,xh,xg,xf,xe,xd,xc,xb,xa
!
! The dimensions for lt and ut have to be at least log (base 2) n
!
      integer   lt(64),ut(64),i,j,k,m,p,q,ib,ie

      ib = 1
      ie = n
!
! Initialize:
!
      j     = ie
      m     = 1
      i     = ib
      iring = iperm+1
      if (iperm.gt.7) iring=1
!
! If this segment has more than two elements  we split it
!
 10   if (j-i-1) 100,90,15
!
! p is the position of an arbitrary element in the segment we choose the
! middle element. Under certain circumstances it may be advantageous
! to choose p at random.
!
 15   p    = (j+i)/2
      ta   = a(p)
      a(p) = a(i)
      go to (21,19,18,17,16,161,162,163),iring
 163     th   = h(p)
         h(p) = h(i)
 162     tg   = g(p)
         g(p) = g(i)
 161     tf   = f(p)
         f(p) = f(i)
 16      te   = e(p)
         e(p) = e(i)
 17      td   = d(p)
         d(p) = d(i)
 18      tc   = c(p)
         c(p) = c(i)
 19      tb   = b(p)
         b(p) = b(i)
 21   continue
!
! Start at the beginning of the segment, search for k such that a(k)>t
!
      q = j
      k = i
 20   k = k+1
      if(k.gt.q)     go to 60
      if(a(k).le.ta) go to 20
!
! Such an element has now been found now search for a q such that a(q)<t
! starting at the end of the segment.
!
 30   continue
      if(a(q).lt.ta) go to 40
      q = q-1
      if(q.gt.k)     go to 30
      go to 50
!
! a(q) has now been found. we interchange a(q) and a(k)
!
 40   xa   = a(k)
      a(k) = a(q)
      a(q) = xa
      go to (45,44,43,42,41,411,412,413),iring
 413     xh   = h(k)
         h(k) = h(q)
         h(q) = xh
 412     xg   = g(k)
         g(k) = g(q)
         g(q) = xg
 411     xf   = f(k)
         f(k) = f(q)
         f(q) = xf
 41      xe   = e(k)
         e(k) = e(q)
         e(q) = xe
 42      xd   = d(k)
         d(k) = d(q)
         d(q) = xd
 43      xc   = c(k)
         c(k) = c(q)
         c(q) = xc
 44      xb   = b(k)
         b(k) = b(q)
         b(q) = xb
 45   continue
!
! Update q and search for another pair to interchange:
!
      q = q-1
      go to 20
 50   q = k-1
 60   continue
!
! The upwards search has now met the downwards search:
!
      a(i)=a(q)
      a(q)=ta
      go to (65,64,63,62,61,611,612,613),iring
 613     h(i) = h(q)
         h(q) = th
 612     g(i) = g(q)
         g(q) = tg
 611     f(i) = f(q)
         f(q) = tf
 61      e(i) = e(q)
         e(q) = te
 62      d(i) = d(q)
         d(q) = td
 63      c(i) = c(q)
         c(q) = tc
 64      b(i) = b(q)
         b(q) = tb
 65   continue
!
! The segment is now divided in three parts: (i,q-1),(q),(q+1,j)
! store the position of the largest segment in lt and ut
!
      if (2*q.le.i+j) go to 70
      lt(m) = i
      ut(m) = q-1
      i = q+1
      go to 80
 70   lt(m) = q+1
      ut(m) = j
      j = q-1
!
! Update m and split the new smaller segment
!
 80   m = m+1
      go to 10
!
! We arrive here if the segment has  two elements we test to see if
! the segment is properly ordered if not, we perform an interchange
!
 90   continue
      if (a(i).le.a(j)) go to 100
      xa=a(i)
      a(i)=a(j)
      a(j)=xa
      go to (95,94,93,92,91,911,912,913),iring
 913     xh   = h(i)
         h(i) = h(j)
         h(j) = xh
 912     xg   = g(i)
         g(i) = g(j)
         g(j) = xg
 911     xf   = f(i)
         f(i) = f(j)
         f(j) = xf
   91    xe   = e(i)
         e(i) = e(j)
         e(j) = xe
   92    xd   = d(i)
         d(i) = d(j)
         d(j) = xd
   93    xc   = c(i)
         c(i) = c(j)
         c(j) = xc
   94    xb   = b(i)
         b(i) = b(j)
         b(j) = xb
   95 continue
!
! If lt and ut contain more segments to be sorted repeat process:
!
 100  m = m-1
      if (m.le.0) go to 110
      i = lt(m)
      j = ut(m)
      go to 10
 110  continue
      return
      end subroutine dblemodsortem
end module sortem
