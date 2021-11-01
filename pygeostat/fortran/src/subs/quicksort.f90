! Module for sorting various arrays
! Author: John Manchuk, 2010
! QUICKSORT modified from:
!   Algorithms and Data Structures in F and Fortran
!   Robin A. Vowels
!   Unicomp, 1998
module quicksort
    implicit none
    
    public :: qsort
    
    public :: insertion_sort
    
    interface qsort
        module procedure integer_sort_1d
        module procedure real_sort_1d
        module procedure double_sort_1d
        module procedure integer_sort_2d
        module procedure real_sort_2d
        module procedure double_sort_2d
    end interface qsort
    
    private
    
    !pointers for recursive sort routines
    integer, pointer :: pa(:,:) => null()
    real,    pointer :: par(:,:) => null()
    real*8,  pointer :: pad(:,:) => null()
    integer, pointer :: pt(:) => null()
    real,    pointer :: ptr(:) => null()
    real*8,  pointer :: ptd(:) => null()
    integer, pointer :: pi(:) => null()

contains

    !Sort a one-dimensional array of integers
    integer function integer_sort_1d ( a, left, right, ord, b ) result ( err )
        integer, intent (in out) :: a(:)
        integer, intent (in) :: left, right
        integer, optional, intent (out) :: ord(:)
        integer, optional, intent (inout) :: b(:)
        integer, allocatable :: indexes(:)
        integer, allocatable :: temp(:)
        integer :: i
        
        !Setup indexing array
        allocate(indexes(left:right), stat = err)
        if(err /= 0)then
            return
        endif
        
        do i = left, right
            indexes(i) = i
        enddo
        
        !Sort the array
        call quick_sort_int (a, indexes, left, right)
        
        !Provide the order?
        if(present(ord))then
            ord = indexes
        endif
        
        !Sort an additional array?
        if(present(b))then
            allocate(temp(left:right), stat = err)
            if(err /= 0)then
                return
            endif
            do i = left, right
                temp(i) = b(i)
            enddo
            do i = left, right
                b(i) = temp(indexes(i))
            enddo
            deallocate(temp)
        endif
        deallocate(indexes)
    end function integer_sort_1d
    
    !Sort a one-dimensional array of reals
    integer function real_sort_1d ( a, left, right, ord, b ) result ( err )
        real*4, intent (in out) :: a(:)
        integer, intent (in) :: left, right
        integer, optional, intent (out) :: ord(:)
        integer, allocatable :: indexes(:)
        real*4, optional, intent (inout) :: b(:)
        real*4, allocatable :: temp(:)
        integer :: i
        
        !Setup indexing array
        allocate(indexes(left:right), stat = err)
        if(err /= 0)then
            return
        endif
        
        do i = left, right
            indexes(i) = i
        enddo
        
        !Sort the array
        call quick_sort_real (a, indexes, left, right)
        
        !Provide the order?
        if(present(ord))then
            ord = indexes
        endif
        
        !Sort an additional array?
        if(present(b))then
            allocate(temp(left:right), stat = err)
            if(err /= 0)then
                return
            endif
            do i = left, right
                temp(i) = b(i)
            enddo
            do i = left, right
                b(i) = temp(indexes(i))
            enddo
            deallocate(temp)
        endif
        deallocate(indexes)
    end function real_sort_1d
    
    !Sort a one-dimensional array of doubles
    integer function double_sort_1d ( a, left, right, ord, b ) result ( err )
        real*8, intent (in out) :: a(:)
        integer, intent (in) :: left, right
        integer, optional, intent (out) :: ord(:)
        integer, allocatable :: indexes(:)
        real*8, optional, intent (inout) :: b(:)
        real*8, allocatable :: temp(:)
        integer :: i
        
        !Setup indexing array
        allocate(indexes(left:right), stat = err)
        if(err /= 0)then
            return
        endif
        
        do i = left, right
            indexes(i) = i
        enddo
        
        !Sort the array
        call quick_sort_double (a, indexes, left, right)
        
        !Provide the order?
        if(present(ord))then
            ord = indexes
        endif
        
        !Sort an additional array?
        if(present(b))then
            allocate(temp(left:right), stat = err)
            if(err /= 0)then
                return
            endif
            do i = left, right
                temp(i) = b(i)
            enddo
            do i = left, right
                b(i) = temp(indexes(i))
            enddo
            deallocate(temp)
        endif
        deallocate(indexes)
    end function double_sort_1d




    !Two dimensional arrays
    
    !Sort a two dimensional array of integers
    integer function integer_sort_2d ( a, left, right, col ) result ( err )
        integer, target, intent (in out) :: a(:,:)
        integer, intent (in) :: left, right
        integer, intent (in) :: col(:)
        integer, allocatable, target :: temp(:)
        integer, allocatable, target :: ids(:)
        integer :: i
        
        !Setup indexing arrays
        allocate(ids(left:right),temp(left:right), stat = err)
        if(err /= 0)then
            return
        endif
        
        do i = left, right
            ids(i) = i
        enddo
        
        pa => a
        pt => temp
        pi => ids
        
        call integer_sort_2d_ ( left, right, col, 1 )

        deallocate(IDS,TEMP)
    end function integer_sort_2d

    !Recursive sort for integer_sort_2d
    recursive subroutine integer_sort_2d_ ( left, right, cols, icol )
        integer, intent (in) :: left, right !left and right bounds
        integer, intent (in) :: cols(:) !columns being sorted
        integer, intent (in) :: icol !current column index
        integer :: i, j, i0, i1
        
        !Do nothing for 1 element
        if(right == left) return
        
        !Sort the current column and reorder the others
        do j = left, right
            pt(j) = pa(cols(icol),j)
        enddo
        
        call quick_sort_int (pt, pi, left, right)
        RLOOP : do i = 1, ubound(pa,1)
            !do not re-order a column that is already sorted
            do j = 1, icol-1
                if(i == cols(j)) cycle RLOOP
            enddo
            do j = left, right
                pt(j) = pa(i,j)
            enddo
            do j = left, right
                pa(i,j) = pt(pi(j))
            enddo
        enddo RLOOP
        
        !Compute bounds based on the most recently sorted column
        i0 = left
        j = cols(icol)
        do while (i0 < right)
            i1 = i0
            pi(i0) = i0
            do while ( pa(j,i1) == pa(j,i0) )
                i1 = i1 + 1
                if(i1 > right) exit
                pi(i1) = i1
            enddo
            i1 = i1 - 1
            
            !Determine sorted order of elements
            if(icol < ubound(cols,1))then
                call integer_sort_2d_ ( i0, i1, cols, icol + 1 )
            endif
            
            !Next subset
            i0 = i1 + 1
        enddo

        return
    end subroutine integer_sort_2d_

    !Sort a two dimensional array of reals
    integer function real_sort_2d ( a, left, right, col ) result ( err )
        real, target, intent (in out) :: a(:,:)
        integer, intent (in) :: left, right
        integer, intent (in) :: col(:)
        real, allocatable, target :: temp(:)
        integer, allocatable, target :: ids(:)
        integer :: i
        
        !Setup indexing arrays
        allocate(ids(left:right),temp(left:right), stat = err)
        if(err /= 0)then
            return
        endif
        
        do i = left, right
            ids(i) = i
        enddo
        
        par => a
        ptr => temp
        pi => ids
        
        call real_sort_2d_ ( left, right, col, 1 )

        deallocate(IDS,TEMP)
    end function real_sort_2d

    !Recursive sort for integer_sort_2d
    recursive subroutine real_sort_2d_ ( left, right, cols, icol )
        integer, intent (in) :: left, right !left and right bounds
        integer, intent (in) :: cols(:) !columns being sorted
        integer, intent (in) :: icol !current column index
        integer :: i, j, i0, i1
        
        !Do nothing for 1 element
        if(right == left) return
        
        !Sort the current column and reorder the others
        do j = left, right
            ptr(j) = par(cols(icol),j)
        enddo
        
        call quick_sort_real (ptr, pi, left, right)
        RLOOP : do i = 1, ubound(par,1)
            !do not re-order a column that is already sorted
            do j = 1, icol-1
                if(i == cols(j)) cycle RLOOP
            enddo
            do j = left, right
                ptr(j) = par(i,j)
            enddo
            do j = left, right
                par(i,j) = ptr(pi(j))
            enddo
        enddo RLOOP
        
        !Compute bounds based on the most recently sorted column
        i0 = left
        j = cols(icol)
        do while (i0 < right)
            i1 = i0
            pi(i0) = i0
            do while ( par(j,i1) == par(j,i0) )
                i1 = i1 + 1
                if(i1 > right) exit
                pi(i1) = i1
            enddo
            i1 = i1 - 1
            
            !Determine sorted order of elements
            if(icol < ubound(cols,1))then
                call real_sort_2d_ ( i0, i1, cols, icol + 1 )
            endif
            
            !Next subset
            i0 = i1 + 1
        enddo

        return
    end subroutine real_sort_2d_
    
    !Sort a two dimensional array of reals
    integer function double_sort_2d ( a, left, right, col ) result ( err )
        real*8, target, intent (in out) :: a(:,:)
        integer, intent (in) :: left, right
        integer, intent (in) :: col(:)
        real*8, allocatable, target :: temp(:)
        integer, allocatable, target :: ids(:)
        integer :: i
        save :: temp, ids
        
        err = 0
        
        !Setup indexing arrays
        if(allocated(ids))then
            if(lbound(ids,1) > left .or. ubound(ids,1) < right)then
                deallocate(ids,temp)
                allocate(ids(left:right),temp(left:right), stat = err)
                if(err /= 0)then
                    return
                endif
            endif
        else
            allocate(ids(left:right),temp(left:right), stat = err)
            if(err /= 0)then
                return
            endif
        endif
        
        do i = left, right
            ids(i) = i
        enddo
        
        pad => a
        ptd => temp
        pi => ids
        
        call double_sort_2d_ ( left, right, col, 1 )
    end function double_sort_2d

    !Recursive sort for integer_sort_2d
    recursive subroutine double_sort_2d_ ( left, right, cols, icol )
        integer, intent (in) :: left, right !left and right bounds
        integer, intent (in) :: cols(:) !columns being sorted
        integer, intent (in) :: icol !current column index
        integer :: i, j, i0, i1
        
        !Do nothing for 1 element
        if(right == left) return
        
        !Sort the current column and reorder the others
        do j = left, right
            ptd(j) = pad(cols(icol),j)
        enddo
        
        call quick_sort_double (ptd, pi, left, right)
        RLOOP : do i = 1, ubound(pad,1)
            !do not re-order a column that is already sorted
            do j = 1, icol-1
                if(i == cols(j)) cycle RLOOP
            enddo
            do j = left, right
                ptd(j) = pad(i,j)
            enddo
            do j = left, right
                pad(i,j) = ptd(pi(j))
            enddo
        enddo RLOOP
        
        !Compute bounds based on the most recently sorted column
        i0 = left
        j = cols(icol)
        do while (i0 < right)
            i1 = i0
            pi(i0) = i0
            do while ( pad(j,i1) == pad(j,i0) )
                i1 = i1 + 1
                if(i1 > right) exit
                pi(i1) = i1
            enddo
            i1 = i1 - 1
            
            !Determine sorted order of elements
            if(icol < ubound(cols,1))then
                call double_sort_2d_ ( i0, i1, cols, icol + 1 )
            endif
            
            !Next subset
            i0 = i1 + 1
        enddo

        return
    end subroutine double_sort_2d_
    

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    !! This is a professional version of Quicksort.
    recursive subroutine quick_sort_int ( A, IDS, LLeft, RRight )
       !! This subroutine implements an improved Quicksort, based on the algorithm by
       !! C. A. R. Hoare.
       !! INCOMING: A      = an array whose elements are to be sorted;
       !!           LLeft  = a pointer to the left-hand end of the partition to be sorted;
       !!           RRight = a pointer to the right-hand end of the partition to be sorted;
       !! OUTGOING: A      = an array whose elements A(LLeft) through A(RRight) are sorted.
       !!           IDS    = an array of indices holding original positions
       integer, intent (inout) :: A( : )
       integer, intent (in out) :: IDS( : )
       integer, intent (in) :: LLeft, RRight
       integer :: Mid, Left, Right
       integer :: Ref, Temp
       integer :: iRef,iTemp, L, R, La
       integer :: No_Swaps, J, Rml
       integer, parameter :: One = 1

       Left = LLeft
       Right = RRight
       if (Right <= Left) then                        !! The partition is empty or contains one element.
          return                                      !! There's no work to do.
       end  if

       !! Section to select the median of three elements, and to move it to the left.
       Mid = (Left + Right)/2
       if (A(Mid)  > A(Right))  then
          call SWAPI (A(Mid), A(Right), IDS(Mid), IDS(Right))
       endif

       if (Left+1 == Right) then                      !! There are 2 elements in the partition,
          return                                      !! & they are now in sort.
       endif
       if (A(Left) > A(Mid))  then
          call SWAPI (A(Left), A(Mid), IDS(Left), IDS(Mid))
       endif
       if (A(Mid)  > A(Right))  then
          call SWAPI (A(Mid), A(Right), IDS(Mid), IDS(Right))
       endif
       if (Left+ 2 == Right) then                     !! There are 3 elements in the partition,
          return                                      !! & they are now in sort.
       endif
       if (A(Mid) == A(Right)) then                   !! Some elements are equal!
          Ref = A(Left)                               !! Forces the left partition to omit equal elements.
         iRef = IDS(Left)
       else
          call SWAPI (A(Left), A(Mid), IDS(Left), IDS(Mid))
          Ref = A(Left)                               !! Select the Reference Element.
         iRef = IDS(Left)
       endif

       L = Left
       R = Right + 1

       !! Partition the elements into three groups.
       No_Swaps = 0
       do
          if (L >= R) then
             exit
          endif
          do L = L + 1, R - 1                         !! Scan from the left for an element
             if (A(L) > Ref) then                     !! larger than the Reference.
                exit
             endif
          enddo

          do                                          !! Scan from the right for an element
             R = R - 1
             if (A(R) <= Ref) then                    !! less than or equal to the Reference.
                exit
             endif
          enddo

          if (L < R) then                             !! Swap two elements that are in the wrong partitions.
             Temp = A(R) ; iTemp  = IDS(R)
             A(R) = A(L) ; IDS(R) = IDS(L)
             A(L) = Temp ; IDS(L) = iTemp
             No_Swaps = No_Swaps + 1                  !! Count each swap as we go.
          endif
       enddo
                                                      !! Partitioning is complete.
       if (Left < R) then                             !! Swap the Reference Element into its final position R in the array.
          A(Left) = A(R) ; IDS(Left) = IDS(R)
          A(R)    = Ref  ; IDS(R)    = iRef
       endif
       !! At this point, A(R) is in its correct position in the list.  Elements A(Left) to A(R-1)
       !! are less than or equal to A(R), and elements A(R+1) to A(Right) are greater then A(R).

       !! Section to find out why no swaps were performed.
       if (No_Swaps == 0) then                        !! Something funny happened: not one element was moved. Investigate cause.
          INCREASING: do
             !! Look for any pre-existing order.
             do J = Left, Right-1
                if (A(J) > A(J+1))  then
                   exit  INCREASING
                endif
             enddo
             return !! The elements are already in order.
          enddo INCREASING

       !! Section to take a strong hand when the maximum number of elements is swapped.
       !! It's possible that the elements were in decreasing order.  Check it.
       elseif (No_Swaps+1 == (Right-Left)/2) then
                                                      !! All possible pairs were swapped. Perhaps the elements were in reverse
       DECREASING: do                                 !! order?  Find out why.
             Rml = Right - Left
             if (iand(Rml, One) /= 0) then            !! A partition containing an even number f elements was disarranged during partitioning.
                if (Left < R-1)   then
                   call SWAPI (A(Left), A(R-1), IDS(Left), IDS(R-1))        !! Restore order.
                endif
             endif
             do J = Left, Right-1                     !! Check that the elements are sorted.
                if (A(J) > A(J+1)) then
                   exit DECREASING
                endif
             enddo
             return                                   !! The entire sub-list is sorted.
          enddo DECREASING
       endif

       do La = R-1, Left+1, -1                        !! Pass over any elements that are equal to Ref.
          if (A(La) /= Ref) then
             exit
          endif
       enddo
        !! At this point, elements A(La+1) through A(R) are equal.
        !! A(L) is in its correct position too, even
        !! if it is not equal to A(Mid)!
        !! But if Left=La already, the partition
        !! was just lop-sided.

       if (Left < La) then
          call QUICK_SORT_INT (A, IDS, Left, La)               !! Partition the left segment.
       endif
                                                      !! The element at R is in its correct position.
       if (R+1 < Right) then
          call QUICK_SORT_INT (A, IDS, R+1, Right)             !! Partition the right segment.
       endif

    end subroutine quick_sort_int
    
    recursive subroutine quick_sort_real ( a, ids, lleft, rright )
        real, intent (inout) :: a( : )
        integer, intent (in out) :: ids( : )
        integer, intent (in) :: lleft, rright
        integer :: mid, left, right
        real :: ref, temp
        integer :: iref, itemp, l, r, la
        integer :: no_swaps, j, rml
        integer, parameter :: one = 1

        if(rright <= lleft)then
            return
        endif
        left = lleft
        right = rright

        mid = (left + right)/2
        if(ptr(mid) > ptr(right))then
            call swapf (ptr(mid), ptr(right), pi(mid), pi(right))
        endif

        if(left+1 == right)then
            return
        endif
        if(ptr(left) > ptr(mid))then
            call swapf (ptr(left), ptr(mid), pi(left), pi(mid))
        endif
        if(ptr(mid)  > ptr(right))then
            call swapf (ptr(mid), ptr(right), pi(mid), pi(right))
        endif
        if(left+2 == right) then
            return
        endif
        if(ptr(mid) == ptr(right))then
            ref = ptr(left)
            iref = ids(left)
        else
            call swapf (ptr(left), ptr(mid), pi(left), pi(mid))
            ref = ptr(left)
            iref = pi(left)
        endif

        l = left
        r = right + 1

        no_swaps = 0
        do
            if(l >= r)then
                exit
            endif
            do l = l + 1, r - 1
                if(ptr(l) > ref)then
                    exit
                endif
            enddo

            do
                r = r - 1
                if(ptr(r) <= ref)then
                    exit
                endif
            enddo

            if(l < r)then 
                temp = ptr(r)
                ptr(r) = ptr(l)
                ptr(l) = temp
                itemp = pi(r)
                pi(r) = pi(l)
                pi(l) = itemp
                no_swaps = no_swaps + 1
            endif
        enddo

        if(left < r)then
            ptr(left) = ptr(r)
            ptr(r) = ref
            pi(left) = pi(r)
            pi(r) = iref
        endif

        if(no_swaps == 0)then
            increasing: do
                do j = left, right-1
                    if (ptr(j) > ptr(j+1))  then
                        exit  increasing
                    endif
                enddo
                return
            enddo increasing
        elseif (no_swaps+1 == (right-left)/2) then
            decreasing: do
                rml = right - left
                if(iand(rml, one) /= 0)then
                    if(left < r-1)then
                        call swapf (ptr(left), ptr(r-1), pi(left), pi(r-1))
                    endif
                endif
                do j = left, right-1
                    if(ptr(j) > ptr(j+1))then
                        exit decreasing
                    endif
                enddo
                return
            enddo decreasing
        endif

        do la = r-1, left+1, -1
            if(ptr(la) /= ref)then
                exit
            endif
        enddo

        if (left < la) then
            call quick_sort_real (a, ids, left, la)
        endif
        
        if (r+1 < right) then
            call quick_sort_real (a, ids, r+1, right)
        endif
        
    end subroutine quick_sort_real
    
    recursive subroutine quick_sort_double ( a, ids, lleft, rright )
        real*8, intent (inout) :: a( : )
        integer, intent (in out) :: ids( : )
        integer, intent (in) :: lleft, rright
        integer :: mid, left, right
        real*8 :: ref, temp
        integer :: iref, itemp, l, r, la
        integer :: no_swaps, j, rml
        integer, parameter :: one = 1

        if(rright <= lleft)then
            return
        endif
        left = lleft
        right = rright

        mid = (left + right)/2
        if(a(mid) > a(right))then
            call swapd (a(mid), a(right), ids(mid), ids(right))
        endif

        if(left+1 == right)then
            return
        endif
        if(a(left) > a(mid))then
            call swapd (a(left), a(mid), ids(left), ids(mid))
        endif
        if(a(mid)  > a(right))then
            call swapd (a(mid), a(right), ids(mid), ids(right))
        endif
        if(left+2 == right) then
            return
        endif
        if(a(mid) == a(right))then
            ref = a(left)
            iref = ids(left)
        else
            call swapd (a(left), a(mid), ids(left), ids(mid))
            ref = a(left)
            iref = ids(left)
        endif

        l = left
        r = right + 1

        no_swaps = 0
        do
            if(l >= r)then
                exit
            endif
            do l = l + 1, r - 1
                if(a(l) > ref)then
                    exit
                endif
            enddo

            do
                r = r - 1
                if(a(r) <= ref)then
                    exit
                endif
            enddo

            if(l < r)then 
                temp = a(r)
                a(r) = a(l)
                a(l) = temp
                itemp = ids(r)
                ids(r) = ids(l)
                ids(l) = itemp
                no_swaps = no_swaps + 1
            endif
        enddo

        if(left < r)then
            a(left) = a(r)
            a(r) = ref
            ids(left) = ids(r)
            ids(r) = iref
        endif

        if(no_swaps == 0)then
            increasing: do
                do j = left, right-1
                    if (a(j) > a(j+1))  then
                        exit  increasing
                    endif
                enddo
                return
            enddo increasing
        elseif (no_swaps+1 == (right-left)/2) then
            decreasing: do
                rml = right - left
                if(iand(rml, one) /= 0)then
                    if(left < r-1)then
                        call swapd (a(left), a(r-1), ids(left), ids(r-1))
                    endif
                endif
                do j = left, right-1
                    if(a(j) > a(j+1))then
                        exit decreasing
                    endif
                enddo
                return
            enddo decreasing
        endif

        do la = r-1, left+1, -1
            if(a(la) /= ref)then
                exit
            endif
        enddo

        if (left < la) then
            call quick_sort_double (a, ids, left, la)
        endif
        
        if (r+1 < right) then
            call quick_sort_double (a, ids, r+1, right)
        endif
        
    end subroutine quick_sort_double

    

    !! this subroutine swaps the element left_element with right_element.
    subroutine swapi ( left_element, right_element, left_id, right_id)
        integer, intent (in out) :: left_element, right_element
        integer, intent (in out) :: left_id, right_id
        integer :: temp, tempid
        temp = left_element
        left_element = right_element
        right_element = temp
        tempid = left_id
        left_id = right_id
        right_id = tempid
    end subroutine swapi

    subroutine swapf ( left_element, right_element, left_id, right_id)
        real, intent (in out) :: left_element, right_element
        integer, intent (in out) :: left_id, right_id
        real :: temp
        integer :: tempid
        temp = left_element
        left_element = right_element
        right_element = temp
        tempid = left_id
        left_id = right_id
        right_id = tempid
    end subroutine swapf

    subroutine swapd ( left_element, right_element, left_id, right_id)
        real*8, intent (in out) :: left_element, right_element
        integer, intent (in out) :: left_id, right_id
        real*8 :: temp
        integer :: tempid
        temp = left_element
        left_element = right_element
        right_element = temp
        tempid = left_id
        left_id = right_id
        right_id = tempid
    end subroutine swapd
    
    !
    ! Simple sort for small arrays (faster than quicksort)
    !
    
    !Insertion sort for 2D array sorted by 1 column
    subroutine insertion_sort ( a, left, right, col )
        real*8, intent (inout) :: a(:,:)
        integer, intent (in) :: left, right, col
        real*8 :: v(ubound(a,1)) 
        integer :: i, j
        do i = left + 1, right
            v = a(:,i)
            j = i - 1
            do while ( j >= left )
                if( a(col,j) > v(col) )then
                    a(:,j + 1) = a(:,j)
                    j = j - 1
                else
                    exit
                endif
            enddo
            a(:,j + 1) = v
        enddo
        return
    end subroutine insertion_sort

end module quicksort