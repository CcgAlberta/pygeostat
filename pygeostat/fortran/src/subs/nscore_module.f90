!
! Module for normal score transform related routines
! Author: John Manchuk, 2010
!
! Module Dependencies:
!   NORMALDIST - used to compute inverse normal distribution
!   QUICKSORT  - for sorting transform tables
!   ACORN_RNG  - random number generator for breaking ties
! 
module normalscore
    use normaldist, ONLY : snorm_inv, snorm_cdf, locate, powint
    use quicksort
    use acorn_rng
    implicit none
    
    !These are general functions.  All other functions in this module are
    !called from these depending on the inputs.
    public :: build_transtable, normal_score 
    
    type :: transtable !used for jagged arrays
        integer :: n
        real*8 :: sumw
        real*8 :: mint
        real*8, allocatable :: table(:,:)
    end type transtable
    public :: transtable
    

    !Machine precision
    real*8, parameter :: MACHP = sqrt(epsilon(real(1,8)))
    
    !Random number generator for breaking ties
    integer, parameter :: seed = 69069
    logical :: isacorninit = .false.
    type(acorn_type) :: rng
    
    !Work array
    real*8, allocatable :: rwork(:,:)
    
    contains
    
    !Compute the normal scores
    subroutine normal_score ( data_, colvr, colct, ct, numq, tmin, tmax, ptt, nsvr, test )
        real*8, pointer :: data_(:,:), nsvr(:)
        integer, intent (in) :: colvr, colct !variable and category columns
        integer, intent (in) :: ct !the category of interest
        integer, intent (in) :: numq
        real*8, intent (in) :: tmin, tmax !trimming limits
        type(transtable), target, intent(in) :: ptt
        logical, intent(out) :: test
        integer :: ncol, nrow, nq
        integer :: i, j, maxiter
        real*8, pointer :: pvr(:),pvrns(:)
        
        !Initialize RNG
        if(.not.isacorninit)then
            call acorninit(rng,seed)
            isacorninit = .true.
        endif
        
        !Choose the version based on inputs
        if( colct < 1 ) then
            ncol = ubound(data_,1)
            nrow = ubound(data_,2)
            nq = ptt%n
            pvr => ptt%table(1,:)
            pvrns => ptt%table(2,:)
            do i = 1, nrow
                if(data_(colvr,i) < tmin .or. data_(colvr,i) >= tmax) then
                    nsvr(i) = -999.
                else
                    nsvr(i) = data_(colvr,i) + acorni(rng)*MACHP
                    maxiter = 0; test = .true.
                    do while ( nsvr(i) < pvr(1) .or. nsvr(i) > pvr(nq) )
                        nsvr(i) = data_(colvr,i) + acorni(rng)*MACHP
                        maxiter = maxiter + 1
                        if(maxiter > 1E7)then 
                            test = .false.
                            return
                        endif
                    enddo
                    call locate(pvr,nq,1,nq,nsvr(i),j)
                    j = min(max(1,j),(nq-1))
                    nsvr(i) = powint(pvr(j),pvr(j+1),pvrns(j),pvrns(j+1),nsvr(i),1.D0)
                endif
            enddo
        else
            call normal_score_by_category ( data_, colvr, colct, ct, &
                numq, tmin, tmax, ptt, nsvr )
            test = .true.
        endif
        return
    end subroutine normal_score
    
    !Compute the normal scores of a category
    subroutine normal_score_by_category ( data_, colvr, colct, ct, numq, tmin, tmax, ptt, nsvr )
        real*8, pointer :: data_(:,:), nsvr(:)
        integer, intent (in) :: colvr, colct !variable and category columns
        integer, intent (in) :: ct !the category of interest
        real*8, intent (in) :: tmin, tmax !trimming limits
        type(transtable), target, intent(in) :: ptt
        integer :: ncol, nrow, numq
        integer :: i, j
        real*8, pointer :: pvr(:),pvrns(:)

        ncol = ubound(data_,1)
        nrow = ubound(data_,2)
        numq = ptt%n
        pvr => ptt%table(1,:)
        pvrns => ptt%table(2,:)
        
        do i = 1, nrow
            if(data_(colvr,i) < tmin .or. data_(colvr,i) >= tmax) then
                nsvr(i) = -999.
            else
                if(floor(data_(colct,i)+0.1) == ct)then
                    nsvr(i) = data_(colvr,i) + acorni(rng)*MACHP
                    call locate(pvr,numq,1,numq,nsvr(i),j)
                    j = min(max(1,j),(numq-1))
                    nsvr(i) = powint(pvr(j),pvr(j+1),pvrns(j),pvrns(j+1),nsvr(i),1.D0)
                endif
            endif
        enddo
        return
    end subroutine normal_score_by_category
    
    !build a transformation table
    logical function build_transtable ( data_, colvr, colwt, colct, ct, numq, tmin, tmax, ptt )
        real*8, pointer :: data_(:,:) !a block of data for the table, must be ncol by nrows
        integer, intent (in) :: colvr, colwt, colct !variable, weight, and category columns
        integer, intent (in) :: ct !the category of interest
        integer, intent (in) :: numq !number of quantiles
        real*8, intent (in) :: tmin, tmax !trimming limits
        type(transtable), intent(out) :: ptt !Transform table to be built
        integer :: ncol, nrow
        integer :: i, j, test
        integer :: numdata
        real*8 :: mindat
        
        !Some initial checks
        build_transtable = .false.
        ncol = ubound(data_,1)
        nrow = ubound(data_,2)
        if(colvr < 1)then
            write(*,'(a)') 'ERROR BUILDING TRANSTABLE, INVALID VARIABLE COLUMN'
            return
        endif
        if(colvr > ncol)then
            write(*,'(a)') 'ERROR BUILDING TRANSTABLE, REQUESTED VARIABLE COLUMN TOO LARGE'
            return
        endif
        if(colwt > ncol)then
            write(*,'(a)') 'ERROR BUILDING TRANSTABLE, REQUESTED WEIGHT COLUMN TOO LARGE'
            return
        endif
        if(colct > ncol)then
            write(*,'(a)') 'ERROR BUILDING TRANSTABLE, REQUESTED CATEGORY COLUMN TOO LARGE'
            return
        endif
        if(numq > 0 .and. numq < 3)then
            write(*,'(a)') 'ERROR BUILDING TRANSTABLE, NUMBER OF QUANTILES MUST BE > 2'
            return
        endif
        if(nrow < 2)then
            write(*,'(a)') 'ERROR BUILDING TRANSTABLE, NOT ENOUGH DATA'
            return
        endif
        if(allocated(ptt%table)) deallocate(ptt%table)
        
        call acorninit(rng,seed)
        numdata = 0
        
        !Choose the version of this routine depending on the inputs
        if(colwt < 1 .and. colct < 1)then
            !Count the data
            do i = 1, nrow
                if(data_(colvr,i) >= tmin .and. data_(colvr,i) < tmax)then
                    numdata = numdata + 1
                endif
            enddo
            if(numdata < 2) return
            
            !Fill in the transform table
            allocate(ptt%table(2,numdata), stat = test)
            if(test /= 0)then
                write(*,'(/,a,/)') 'ERROR BUILDING TRANSTABLE, ALLOCATION FAILED'
                return
            endif
            j = 0
            mindat = huge(mindat)
            do i = 1, nrow
                if(data_(colvr,i) >= tmin .and. data_(colvr,i) < tmax)then
                    j = j + 1
                    if(data_(colvr,i) < mindat) mindat = data_(colvr,i)
                    ptt%table(1,j) = data_(colvr,i) + acorni(rng)*MACHP
                    ptt%table(2,j) = 1.
                endif
            enddo
            ptt%mint = mindat
                
            !define size of the table
            ptt%n = numdata
            ptt%sumw = real(numdata,8)

            ! Compute nscores
            call nscore_transtable ( ptt, numq )
            
            build_transtable = .true.
        else
            if(colwt < 1)then !colct > 0
                build_transtable = build_transtable_by_category ( data_, &
                             colvr, colct, ct, numq, tmin, tmax, ptt )
            elseif(colct < 1)then !colwt > 0
                build_transtable = build_transtablew ( data_, &
                             colvr, colwt, numq, tmin, tmax, ptt )
            else !colwt and colct > 0
                build_transtable = build_transtablew_by_category ( data_, &
                             colvr, colwt, colct, ct, numq, tmin, tmax, ptt )
            endif
        endif
    
        return
    end function build_transtable
    
    !Build a weighted transform table
    ! this routine is called from build_transtable so error checking is not needed
    logical function build_transtablew ( data_, colvr, colwt, numq, tmin, tmax, ptt )
        real*8, pointer :: data_(:,:) !a block of data for the table, must be ncol by nrows
        integer, intent (in) :: colvr, colwt !variable and category columns
        integer, intent (in) :: numq !number of quantiles
        real*8, intent (in) :: tmin, tmax !trimming limits
        type(transtable), intent(inout) :: ptt !pointer to the transform table to be built
        integer :: ncol, nrow
        integer :: i, j, test
        integer :: numdata
        real*8 :: sumweight
        
        build_transtablew = .false.
        ncol = ubound(data_,1)
        nrow = ubound(data_,2)
        numdata = 0
        sumweight = 0.
        
        !Count the data
        do i = 1, nrow
            if(data_(colvr,i) >= tmin .and. data_(colvr,i) < tmax)then
                numdata = numdata + 1
            endif
        enddo
        if(numdata < 2) return
        
        !Fill in the transform table
        allocate(ptt%table(2,numdata), stat = test)
        if(test /= 0)then
            write(*,'(/,a,/)') 'ERROR BUILDING TRANSTABLE, ALLOCATION FAILED'
            return
        endif
        j = 0
        do i = 1, nrow
            if(data_(colvr,i) >= tmin .and. data_(colvr,i) < tmax )then
                j = j + 1
                ptt%table(1,j) = data_(colvr,i) + acorni(rng)*MACHP
                ptt%table(2,j) = data_(colwt,i)
                if(ptt%table(2,j) < MACHP )then ! Check if the weight is numerically small
                    ptt%table(2,j) = MACHP
                endif
                sumweight = sumweight + ptt%table(2,j)
            endif
        enddo
        
        !define size of the table
        ptt%n = numdata
        ptt%sumw = sumweight
        
        ! Compute nscores
        call nscore_transtable ( ptt, numq )
        
        build_transtablew = .true.
        return
    end function build_transtablew
    
    !Build a transform table by category
    ! this routine is called from build_transtable so error checking is not needed
    logical function build_transtable_by_category ( data_, colvr, colct, ct, numq, tmin, tmax, ptt )
        real*8, pointer :: data_(:,:) !a block of data for the table, must be ncol by nrows
        integer, intent (in) :: colvr, colct !variable and category columns
        integer, intent (in) :: ct !category of interest
        integer, intent (in) :: numq !number of quantiles
        real*8, intent (in) :: tmin, tmax !trimming limits
        type(transtable), intent(inout) :: ptt !pointer to the transform table to be built
        integer :: ncol, nrow
        integer :: i, j, test
        integer :: numdata
        
        build_transtable_by_category = .false.
        ncol = ubound(data_,1)
        nrow = ubound(data_,2)
        numdata = 0
        
        !Count the data
        do i = 1, nrow
            if(data_(colvr,i) >= tmin .and. data_(colvr,i) < tmax .and. &
               floor(data_(colct,i)+0.1) == ct)then
                numdata = numdata + 1
            endif
        enddo
        if(numdata < 2) return
        
        !Fill in the transform table
        allocate(ptt%table(2,numdata), stat = test)
        if(test /= 0)then
            write(*,'(/,a,/)') 'ERROR BUILDING TRANSTABLE, ALLOCATION FAILED'
            return
        endif
        j = 0
        do i = 1, nrow
            if(data_(colvr,i) >= tmin .and. data_(colvr,i) < tmax .and. &
               floor(data_(colct,i)+0.1) == ct)then
                j = j + 1
                ptt%table(1,j) = data_(colvr,i) + acorni(rng)*MACHP
                ptt%table(2,j) = 1.
            endif
        enddo
        
        !define size of the table
        ptt%n = numdata
        ptt%sumw = real(numdata,8)
        
        ! Compute nscores
        call nscore_transtable ( ptt, numq )
        
        build_transtable_by_category = .true.
        return
    end function build_transtable_by_category
    
    !Build a weighted transform table by category
    ! this routine is called from build_transtable so error checking is not needed
    logical function build_transtablew_by_category ( data_, colvr, colwt, colct, ct, numq, tmin, tmax, ptt )
        real*8, pointer :: data_(:,:) !a block of data for the table, must be ncol by nrows
        integer, intent (in) :: colvr, colwt, colct !variable and category columns
        integer, intent (in) :: ct !category of interest
        integer, intent (in) :: numq !number of quantiles
        real*8, intent (in) :: tmin, tmax !trimming limits
        type(transtable), intent(inout) :: ptt !pointer to the transform table to be built
        integer :: ncol, nrow
        integer :: i, j, test
        integer :: numdata
        real*8 :: sumweight
        
        build_transtablew_by_category = .false.
        ncol = ubound(data_,1)
        nrow = ubound(data_,2)
        numdata = 0
        sumweight = 0.
        
        !Count the data
        do i = 1, nrow
            if(data_(colvr,i) >= tmin .and. data_(colvr,i) < tmax .and. &
               floor(data_(colct,i)+0.1) == ct)then
                numdata = numdata + 1
            endif
        enddo
        if(numdata < 2) return
        
        !Fill in the transform table
        allocate(ptt%table(2,numdata), stat = test)
        if(test /= 0)then
            write(*,'(/,a,/)') 'ERROR BUILDING TRANSTABLE, ALLOCATION FAILED'
            return
        endif
        j = 0
        do i = 1, nrow
            if(data_(colvr,i) >= tmin .and. data_(colvr,i) < tmax .and. &
               floor(data_(colct,i)+0.1) == ct)then
                j = j + 1
                ptt%table(1,j) = data_(colvr,i) + acorni(rng)*MACHP
                ptt%table(2,j) = data_(colwt,i)
                if(ptt%table(2,j) < MACHP )then ! Check if the weight is numerically small
                    ptt%table(2,j) = MACHP
                endif
                sumweight = sumweight + ptt%table(2,j)
            endif
        enddo
        
        !define size of the table
        ptt%n = numdata
        ptt%sumw = sumweight
        
        ! Compute nscores
        call nscore_transtable ( ptt, numq )
        
        build_transtablew_by_category = .true.
        return
    end function build_transtablew_by_category
    
    ! Compute the normal scores for a transform table
    ! The table contains the number of data and total weight needed.
    ! If numq is passed in > 0, then the size of the table is
    ! cut back.
    subroutine nscore_transtable ( ptt, numq )
        type(transtable), target, intent(inout) :: ptt
        integer, intent (in) :: numq
        real*8 :: w, wf, cp0, cp1
        real*8, pointer :: z(:),q(:)
        integer :: i, j, n
        
        ! Sort data by value (this sort call sorts trtab by column 1)
        n = ptt%n
        i = qsort(ptt%table,1,n,(/1/))
        wf = 1. / ptt%sumw
        cp0 = 0.
        cp1 = 0.
        
        if(n > numq .and. numq > 0)then ! Resize for specified numq
            !Compute CDF
            do i = 1, n
                w =  wf * ptt%table(2,i)
                cp1 = cp1 + w
                ptt%table(2,i) = (cp1 + cp0) / 2.
                cp0 = cp1
            enddo
            
            !Generate a set of equally spaced normal values covering
            !the range of quantiles in the full size transform table
            !cp0 and cp1 are reused for miny and maxy and w as the diff
            allocate(rwork(2,numq))
            z => ptt%table(1,:)
            q => ptt%table(2,:)
            cp0 = snorm_inv(q(1))
            cp1 = snorm_inv(q(n))
            w = (cp1 - cp0) / real(numq-1,8)

            do i = 1, numq
                rwork(2,i) = snorm_cdf(cp0+real(i-1,8)*w)
                call locate(q,n,1,n,rwork(2,i),j)
                j = min(max(1,j),(n-1))
                rwork(1,i) = powint(q(j),q(j+1),z(j),z(j+1),rwork(2,i),1.D0)
                rwork(2,i) = snorm_inv ( rwork(2,i) )
            enddo
            rwork(1,numq) = ptt%table(1,n)
            
            deallocate(ptt%table)
            allocate(ptt%table(2,numq))
            ptt%n = numq
            ptt%table(1,:) = rwork(1,:)
            ptt%table(2,:) = rwork(2,:)
            
            !Add small number to max so there is room for breaking ties
            ptt%table(1,1) = ptt%mint
            ptt%table(1,numq) = ptt%table(1,numq) + MACHP
            
            nullify(z)
            nullify(q)
            deallocate(rwork)
        else
            ! Compute the cumulative probabilities and nscore values
            do i = 1, n
                w =  wf * ptt%table(2,i)
                cp1 = cp1 + w
                ptt%table(2,i) = snorm_inv ( (cp1 + cp0) / 2. )
                cp0 = cp1
            enddo
            
            !Add small number to max so there is room for breaking ties
            ptt%table(1,n) = ptt%table(1,n) + MACHP
            
        endif

        return
    end subroutine nscore_transtable
    
end module normalscore