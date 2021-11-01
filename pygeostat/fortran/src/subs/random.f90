! public
!Module with types and functions for generating random numbers
!in gslib type programs

!No actual defined types are present

! Copyright (C) 1997 Makoto Matsumoto and Takuji Nishimura.
! When you use this, send an email to: matumoto@math.keio.ac.jp
! with an appropriate reference to your work.

!***********************************************************************
! Fortran translation by Hiroshi Takano.  Jan. 13, 1999.

!   genrand()      -> double precision function grnd()
!   sgenrand(seed) -> subroutine sgrnd(seed)
!                     integer seed

! This program uses the following standard intrinsics.
!   ishft(i,n): If n > 0, shifts bits in i by n positions to left.
!               If n < 0, shifts bits in i by n positions to right.
!   iand (i,j): Performs logical AND on corresponding bits of i and j.
!   ior  (i,j): Performs inclusive OR on corresponding bits of i and j.
!   ieor (i,j): Performs exclusive OR on corresponding bits of i and j.

!***********************************************************************

!WARNING - DO NOT CHANGE THE INTEGER TYPE, OTHERWISE BIT-SHIFT WILL NOT
!          OPERATE PROPERLY AND RANDOM NUMBERS WILL BE WRONG.

module random
    implicit none
    PUBLIC :: init_genrand, & !Initialize the generator with a seed
              grnd            !draw a random number

    !PRIVATE :: tshftu, tshfts, tshftt, tshftl

    integer, private, parameter :: n = 624, n1 = n+1, m = 397, mata = -1727483681 ! constant vector a
    integer, private, parameter :: umask = -2147483647 - 1                        ! most significant w-r bits
    integer, private, parameter :: lmask =  2147483647                            ! least significant r bits
    integer, private, parameter :: tmaskb= -1658038656, tmaskc= -272236544        ! the array for the state vector
    integer, private, save      :: mt(0:n-1), mti = n1                            ! mti==N+1 means mt[N] is not initialized

contains

    subroutine init_genrand(seed)
    ! This initialization is based upon the multiplier given on p.106 of the
    ! 3rd edition of Knuth, The Art of Computer Programming Vol. 2.
    ! This version assumes that integer overflow does NOT cause a crash.
      integer, intent (in) :: seed
      integer :: latest

      mt(0) = seed
      latest = seed
      do mti = 1, n-1
        latest = IEOR( latest, ISHFT( latest, -30 ) )
        latest = latest * 1812433253 + mti
        mt(mti) = latest
      enddo

      return
    end subroutine init_genrand

    real*8 function grnd ( ) result ( fn_val )
        integer, SAVE :: mag01(0:1) = (/ 0, mata /) ! mag01(x) = x * MATA for x=0,1
        integer       :: kk, y

        ! These statement functions have been replaced with separate functions
        ! tshftu(y) = ISHFT(y,-11)
        ! tshfts(y) = ISHFT(y,7)
        ! tshftt(y) = ISHFT(y,15)
        ! tshftl(y) = ISHFT(y,-18)

        if( mti >= n )then ! generate N words at one time
            if( mti == n+1 )then ! if sgrnd() has not been called,
                call init_genrand ( 69069 )  ! a default initial seed is used
            endif
            do kk = 0, n - m - 1
                y = IOR(IAND(mt(kk),umask), IAND(mt(kk+1),lmask))
                mt(kk) = IEOR(IEOR(mt(kk+m), ISHFT(y,-1)),mag01(IAND(y,1)))
            enddo
            do kk = n-m, n-2
                y = IOR(IAND(mt(kk),umask), IAND(mt(kk+1),lmask))
                mt(kk) = IEOR(IEOR(mt(kk+(m-n)), ISHFT(y,-1)),mag01(IAND(y,1)))
            enddo
            y = IOR(IAND(mt(n-1),umask), IAND(mt(0),lmask))
            mt(n-1) = IEOR(IEOR(mt(m-1), ISHFT(y,-1)),mag01(IAND(y,1)))
            mti = 0
        endif

        y = mt(mti)
        mti = mti + 1
        y = IEOR(y, tshftu(y))
        y = IEOR(y, IAND(tshfts(y),tmaskb))
        y = IEOR(y, IAND(tshftt(y),tmaskc))
        y = IEOR(y, tshftl(y))

        if( y < 0 )then
            fn_val = (real(y,8) + 2.0D0**32) / (2.0D0**32 - 1.0D0)
        else
            fn_val = real(y,8) / (2.0D0**32 - 1.0D0)
        endif

        return
    end function grnd

    integer function tshftu ( y ) RESULT( fn_val )
      integer, intent(in) :: y
      fn_val = ISHFT(y,-11)
      return
    end function tshftu

    integer function tshfts(y) RESULT(fn_val)
      integer, intent(in) :: y
      fn_val = ISHFT(y,7)
      return
    end function tshfts

    integer function tshftt(y) RESULT(fn_val)
      integer, intent(in) :: y
      fn_val = ISHFT(y,15)
      return
    end function tshftt

    integer function tshftl(y) RESULT(fn_val)
      integer, intent(in) :: y
      fn_val = ISHFT(y,-18)
      return
    end function tshftl

end module random
