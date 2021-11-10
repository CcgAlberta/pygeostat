! public
! ACORN random number generation module.
! Author: John Manchuk, 2010
!
! function ACORNI originally from GSLIB (Deutsch, C.V., 1992)
!
! This was created so that multiple independent instatiations of the ACORN
! RNG could be implemented within the same program
!
module acorn_rng
    implicit none

    public :: acorninit, & !initialize an RNG with a seed (must be called first)
              acorni, &    !generate a random number
              acorniperm    !generate a random permutation

    type, public :: acorn_type
        integer :: KORDEI = 12
        integer :: MAXOP1 = 13 !KORDEI + 1
        real*8  :: MAXINT = 2**30
        integer :: ixv(13)
    end type acorn_type

    contains

    real*8 function acorni( acorn )
    !-----------------------------------------------------------------------
    ! Fortran implementation of ACORN random number generator of order less
    ! than or equal to 12 (higher orders can be obtained by increasing the
    ! parameter value MAXORD).
    ! NOTES: 1. The variable idum is a dummy variable. The common block
    !           IACO is used to transfer data into the function.
    !        2. Before the first call to ACORN the common block IACO must
    !           be initialised by the user, as follows. The values of
    !           variables in the common block must not subsequently be
    !           changed by the user.
    !
    !             KORDEI - order of generator required ( must be =< MAXORD)
    !             MAXINT - modulus for generator, must be chosen small
    !                      enough that 2*MAXINT does not overflow
    !             ixv(1) - seed for random number generator
    !                      require 0 < ixv(1) < MAXINT
    !             (ixv(I+1),I=1,KORDEI)
    !                    - KORDEI initial values for generator
    !                      require 0 =< ixv(I+1) < MAXINT
    !        3. After initialisation, each call to ACORN generates a single
    !           random number between 0 and 1.
    !        4. An example of suitable values for parameters is
    !             KORDEI   = 10
    !             MAXINT   = 2**30
    !             ixv(1)   = an odd integer in the (approximate) range
    !                        (0.001 * MAXINT) to (0.999 * MAXINT)
    !             ixv(I+1) = 0, I=1,KORDEI
    ! Author: R.S.Wikramaratna,                           Date: October 1990
    !-----------------------------------------------------------------------
        type(acorn_type), intent (inout) :: acorn
        integer :: i

        do i = 1, acorn%KORDEI
            acorn%ixv(i+1) = acorn%ixv(i+1) + acorn%ixv(i)
            if(acorn%ixv(i+1) >= acorn%MAXINT)then
                acorn%ixv(i+1) = acorn%ixv(i+1) - acorn%MAXINT
            endif
        end do

        acorni = real(acorn%ixv(acorn%MAXOP1), 8) / acorn%MAXINT
        return
    end function acorni

    !Initialize a random number generator
    subroutine acorninit(acorn,seed)
        type(acorn_type), intent (inout) :: acorn
        integer, intent (in) :: seed
        integer :: i
        real*8 :: p
        acorn%ixv = 0
        acorn%ixv(1) = seed
        do i = 1,seed
            p = acorni(acorn)
        enddo
        return
    end subroutine acorninit

    !Generate a random permutation
    subroutine acorniperm ( acorn, perm )
        type(acorn_type), intent (inout) :: acorn
        integer, intent (inout) :: perm(:)
        integer :: i, j, i0, i1, n

        n = ubound(perm,1)
        perm = 0
        do i = 1, n
            j = nint(acorni(acorn) * n)
            if(j < 1)then
                j = 1
            elseif(j > n)then
                j = n
            endif
            if(perm(j) /= 0)then
                i0 = j
                i1 = j
                do
                    if(i0 > 1)then
                        i0 = i0 - 1
                        if(perm(i0) < 1)then
                            perm(i0) = i
                            exit
                        endif
                    endif
                    if(i1 < n)then
                        i1 = i1 + 1
                        if(perm(i1) < 1)then
                            perm(i1) = i
                            exit
                        endif
                    endif
                enddo
            else
                perm(j) = i
            endif
        enddo

        return
    end subroutine acorniperm

end module acorn_rng
