module nscore 
use normalscore
implicit none

! Module Variables
logical, public :: init = .false.
type(transtable), target, allocatable, private :: ttb(:,:)
    contains
subroutine buildtranstables(nd, nvar, nquantile, ncat, catcol, catcodes, wtcol, ltrim, utrim, &
                            indat, test)
    !-----------------------------------------------------------------------------------------
    ! Subroutine to get the normal scores of the input vector contained in indat
    ! array
    !
    ! Parameters:
    !   nd : number of data considered
    !   nvar : number of variables considered
    !   nquantile : number of quantiles to discretize the transform table 
    !   ncat : number of rocktypes to transform according to
    !   catcol : location in the indat table containing the categories, 0 if none 
    !   catcodes: numcat long array of category codes
    !   wtcol : location in the indat table containing the weights, 0 if no weights
    !   ltrim : lower trimming limit
    !   utrim : upper trimming limit
    !   indat : [nvar+(1+1) x nd) array of data to normal score ( weights and categories 
    !           appended to this table if used)
    ! Returns: 
    !   test : false if something went wrong!
    !
    ! Ryan Martin : python-callable wrapper for nscore_module.f90 by John Manchuk
    !-----------------------------------------------------------------------------------------
    integer, intent(in) :: nd, nvar, nquantile, ncat, wtcol, catcol, catcodes(ncat)
    real*8, intent(in)  :: ltrim, utrim
    real*8, target, intent(in) :: indat(nvar,nd)
    logical, intent(out) :: test
    real*8, pointer :: dat_pt(:,:) => null()
    real*8, pointer :: ns_slice(:) => null()
    integer  :: i, j, ivr, ict, wt_col, cat_code, num_variable

    ! setup the module variables
    init = .false.

    ! Clear the old storage, if it exists:
    if(allocated(ttb))then 
        do i = 1,ubound(ttb,1)
            do j = 1,ubound(ttb,2)
                if(allocated(ttb(i,j)%table)) deallocate(ttb(i,j)%table)
            enddo
        enddo
        deallocate(ttb)
    endif
    ! Figure out how many extra columns added to the array of data: 
    if(catcol > 0 .and. wtcol > 0)then 
        num_variable = nvar - 2
    elseif(catcol > 0 .or. wtcol > 0)then 
        num_variable = nvar - 1
    else
        num_variable = nvar
    endif
    allocate(ttb(num_variable,ncat))

    dat_pt => indat
    test = .true.

    do ivr = 1,num_variable          
        do ict = 1,ncat
            cat_code = catcodes(ict)      
            test = build_transtable( dat_pt, ivr, wtcol, catcol, cat_code, &
                                     nquantile, ltrim, utrim, ttb(ivr,ict)     )
            if(.not.test) return ! the logical returned false, something went wrong
        enddo
    enddo

    ! set the module variable to true
    init = .true.
    nullify(dat_pt)
end subroutine buildtranstables

subroutine normalscoredata(nd, nvar, nquantile, ncat, catcol, catcodes, ltrim, utrim, indat, &
                           nsdat, test)
    !--------------------------------------------------------------------------------
    ! Subroutine to get the normal scores of the input vector contained in indat
    ! array
    !
    ! Parameters:
    !   nd : number of data considered
    !   nvar : number of variables considered
    !   nquantile : number of quantiles to discretize the transform table 
    !   ncat : number of rocktypes to transform according to
    !   catcol : location in the indat table containing the categories, 0 if none 
    !   catcodes: numcat long array of category codes
    !   ltrim : lower trimming limit
    !   utrim : upper trimming limit
    !   indat : (nvar x nd) array of data to normal score
    ! Returns: 
    !   nsdat : (nvar x nd) array of normal scored data
    !   test  : returns false if the thing fails where despiking probably matters!
    !
    ! Ryan Martin : python / simplified wrapper for nscore_module.f90 by John Manchuk
    !--------------------------------------------------------------------------------
    integer, intent(in) :: nd, nvar, nquantile, ncat, catcol, catcodes(ncat)
    real*8, intent(in)  :: ltrim, utrim
    real*8, target, intent(in) :: indat(nvar,nd)
    real*8, target, intent(out) :: nsdat(nvar,nd)
    logical, intent(out) :: test
    real*8, pointer :: ns_slice(:) => null()
    real*8, pointer :: dat_pt(:,:) => null()
    integer  :: ict, i, num_variable, ivr

    if(.not. init)then
        write(*,*) 'Please initialize this module by calling buildtranstables()'
        return
    endif
    ! Figure out how many extra columns added to the array of data: 
    if(catcol > 0)then 
        num_variable = nvar - 1
    else
        num_variable = nvar
    endif

    dat_pt => indat
    nsdat = 0.D+00

    do ivr = 1,num_variable
        do ict = 1,ncat
            ns_slice => nsdat(ivr,:)
            call normal_score( dat_pt, ivr, catcol, catcodes(ict), nquantile, &
                               ltrim, utrim, ttb(ivr,ict), ns_slice,  test   )
            if(.not.test) return ! the logical returned false, something went wrong
        enddo
    enddo

    nullify(dat_pt,ns_slice)

end subroutine normalscoredata

subroutine gettranstable(varnum, category, nquantile, outtrans)
    !--------------------------------------------------------------------------------
    ! Subroutine to get the transform table for the requested variable 
    !
    ! Parameters:
    !   varnum : the 1-indexed variable number from the original data array to get the 
    !            the trans
    !   category : The category from which to pull the data..
    !   nquantile : number of quantiles to discretize the transform table 
    ! Returns: 
    !   outtrans : (2 x nquantile) array of normal scored data
    !
    ! Ryan Martin : python / simplified wrapper for nscore_module.f90 by John Manchuk
    !--------------------------------------------------------------------------------

    integer, intent(in) :: varnum, nquantile, category
    real*8,  intent(out):: outtrans(2, nquantile)

    if(varnum < 1 .or. varnum > ubound(ttb,1))then 
        write(*,*) 'Invalid variable inquiry'
        return
    endif

    if(.not. init)then
        write(*,*) 'Please initialize this module by calling buildtranstables()'
        return
    endif

    outtrans = ttb(varnum,category)%table

end subroutine gettranstable

end module nscore