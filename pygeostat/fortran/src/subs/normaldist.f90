!Module with types and functions for using normal distributions
!in gslib type programs
! Author: John Manchuk, 2010
! References:   Michael Wichura,
!               The Percentage Points of the Normal Distribution,
!               Algorithm AS 241,
!               Applied Statistics,
!               Volume 37, Number 3, pages 477-484, 1988.
! Some routines were written by John Burkardt
!
module normaldist
    implicit none
    
    !Normal distribution functions
    PUBLIC ::  norm_pdf, & !calculate the non-standard normal density
              snorm_pdf, & !calculate the standard normal density
               norm_cdf, & !calculate the non-standard normal cumulative probability
              snorm_cdf, & !calculate the standard normal cumulative probability
              snorm_inv    !calculate the inverse standard normal cdf
              
    PUBLIC :: locate, powint
    
    PRIVATE
    
    real*8, parameter :: PI = 3.141592653589793238462643383279502884197
    real*8, parameter :: SCALE = 1. / ( 2. * PI )
    
contains

    !***********************************!
    !   NORMAL DISTRIBUTION FUNCTIONS   !
    !***********************************!

    !Calculate the non-standard probability density for a value
    real*8 function norm_pdf ( x, m, v ) result ( y )
        real*8, intent (in) :: x, m, v !value, mean, variance
        y = snorm_pdf ( (x - m) / dsqrt( v ) )
    end function norm_pdf
    
    !Calculate the standard normal probability density for a value
    real*8 function snorm_pdf ( x ) result ( y )
        real*8, intent (in) :: x
        y = dsqrt ( SCALE ) * exp ( -0.5 * x**2 )
    end function snorm_pdf
    
    !Calculate the non-standard cumulative probability for a value
    real*8 function norm_cdf ( x, m, v ) result ( cdf )
        real*8, intent (in) :: x, m, v !value, mean, variance
        cdf = snorm_cdf ( (x - m) / dsqrt( v ) )
    end function norm_cdf
    
    !Calculate the standard cumulative probability for a value
    real*8 function snorm_cdf ( x ) result ( cdf )
        real*8, intent (in) :: x
        real*8, parameter :: a1 = 0.398942280444D+00, a2 = 0.399903438504D+00, a3 = 5.75885480458D+00, &
	                         a4 = 29.8213557808D+00,  a5 = 2.62433121679D+00,  a6 = 48.6959930692D+00, &
	                         a7 = 5.92885724438D+00
	    real*8, parameter :: b0 = 0.398942280385D+00, b1 = 3.8052D-08,         b2 = 1.00000615302D+00, &
	                         b3 = 3.98064794D-04,     b4 = 1.98615381364D+00,  b5 = 0.151679116635D+00, &
	                         b6 = 5.29330324926D+00,  b7 = 4.8385912808D+00,   b8 = 15.1508972451D+00, &
	                         b9 = 0.742380924027D+00, b10 = 30.789933034D+00,  b11 = 3.99019417011D+00
	    real*8 :: q, y

        !  |X| <= 1.28.
	    if( abs( x ) <= 1.28 )then
	        y = 0.5 * x**2
	        q = 0.5 - abs ( x ) * (a1 - a2 * y / (y + a3 - a4 / (y + a5 + a6 / (y + a7 ))))

        !  1.28 < |X| <= 12.7
	    elseif( abs( x ) <= 12.7 )then
	        y = 0.5 * x**2
	        q = exp( -y ) * b0 / ( abs ( x ) - b1                 &
                + b2  / (abs( x ) + b3 + b4 / (abs( x ) - b5   &
                + b6  / (abs( x ) + b7 - b8 / (abs( x ) + b9   &
                + b10 / (abs( x ) + b11))))))
    
        !  12.7 < |X|
	    else
	        q = 0.0
	    end if
    
        !  Take account of negative X.
        if(x < 0.0)then
            cdf = q
        else
            cdf = 1.0 - q
        end if
    end function snorm_cdf
    
    !Invert the normal cdf at a specified quantile
    real*8 function snorm_inv ( p ) result ( x )
        !  Discussion: The result is accurate to about 1 part in 10**16.
        !  Modified: 27 December 2004
        !  Author: Michael Wichura. FORTRAN90 version by John Burkardt
        !  Reference:
        !    Michael Wichura,
        !    The Percentage Points of the Normal Distribution,
        !    Algorithm AS 241,
        !    Applied Statistics,
        !    Volume 37, Number 3, pages 477-484, 1988.
        real*8, intent (in) :: p
        real*8, parameter, dimension ( 8 ) :: a = (/ 3.3871328727963666080D+00, 1.3314166789178437745D+02, &
            1.9715909503065514427D+03, 1.3731693765509461125D+04, 4.5921953931549871457D+04, &
            6.7265770927008700853D+04, 3.3430575583588128105D+04, 2.5090809287301226727D+03 /)
        real*8, parameter, dimension ( 8 ) :: b = (/ 1.0000000000000000000D+00, 4.2313330701600911252D+01, &
            6.8718700749205790830D+02, 5.3941960214247511077D+03, 2.1213794301586595867D+04, &
            3.9307895800092710610D+04, 2.8729085735721942674D+04, 5.2264952788528545610D+03 /)
        real*8, parameter, dimension ( 8 ) :: c = (/ 1.42343711074968357734D+00, 4.63033784615654529590D+00, &
            5.76949722146069140550D+00, 3.64784832476320460504D+00, 1.27045825245236838258D+00, &
            2.41780725177450611770D-01, 2.27238449892691845833D-02, 7.74545014278341407640D-04 /)
        real*8, parameter, dimension ( 8 ) :: d = (/ 1.00000000000000000000D+00, 2.05319162663775882187D+00, &
            1.67638483018380384940D+00, 6.89767334985100004550D-01, 1.48103976427480074590D-01, &
            1.51986665636164571966D-02, 5.47593808499534494600D-04, 1.05075007164441684324D-09 /)
        real*8, parameter, dimension ( 8 ) :: e = (/ 6.65790464350110377720D+00, 5.46378491116411436990D+00, &
            1.78482653991729133580D+00, 2.96560571828504891230D-01, 2.65321895265761230930D-02, &
            1.24266094738807843860D-03, 2.71155556874348757815D-05, 2.01033439929228813265D-07 /)
        real*8, parameter, dimension ( 8 ) :: f = (/ 1.00000000000000000000D+00, 5.99832206555887937690D-01, &
            1.36929880922735805310D-01, 1.48753612908506148525D-02, 7.86869131145613259100D-04, &
            1.84631831751005468180D-05, 1.42151175831644588870D-07, 2.04426310338993978564D-15 /)
        real*8, parameter :: split1 = 0.425
        real*8, parameter :: split2 = 5.0
        real*8, parameter :: const1 = 0.180625
        real*8, parameter :: const2 = 1.6
        real*8 :: q, r

        !Is p valid ( 0 <= p <= 1 )
        if ( p <= 0.0 ) then
            x = -huge ( p ) ; return
        endif

        if ( p >= 1.0 ) then
            x = huge ( p ) ; return
        endif

        q = p - 0.5

        if( abs ( q ) <= split1 )then
            r = const1 - q * q
            x = q * dpoly_value ( 8, a, r ) / dpoly_value ( 8, b, r )
        else
            if( q < 0.0 )then
                r = p
            else
                r = 1.0 - p
            endif

            if( r <= 0.0 )then
                x = -1.0
                STOP
            endif

            r = dsqrt ( -dlog ( r ) )

            if( r <= split2 )then
                r = r - const2
                x = dpoly_value ( 8, c, r ) / dpoly_value ( 8, d, r )
            else
                r = r - split2
                x = dpoly_value ( 8, e, r ) / dpoly_value ( 8, f, r )
            endif

            if ( q < 0.0 ) then
                x = -x
            endif
        endif
        
        contains
        real*8 function dpoly_value ( n, a, x ) result ( p )
            !Given N and A, the form of the polynomial is:
            !  p(x) = a(1) + a(2) * x + ... + a(n-1) * x^(n-2) + a(n) * x^(n-1)
            !  Modified: 13 August 2004
            !  Author: John Burkardt
            integer, intent (in) :: n
            real*8,    intent (in) :: x
            real*8, dimension ( n ), intent (in) :: a
            integer :: i
            p = 0.0
            do i = n, 1, -1
                p = p * x + a(i)
            enddo
        end function dpoly_value
    end function snorm_inv

    !Power interpolate the value of y between (xlow,ylow) and (xhigh,yhigh)
    !for a value of x and a power pow.
    real*8 function powint ( xlow, xhigh, ylow, yhigh, xval, pow)
        real*8, parameter   :: EPSLON = 1.0e-20
        real*8, intent (in) :: xlow, xhigh, ylow, yhigh, xval, pow

        if((xhigh - xlow) < EPSLON) then
            powint = (yhigh + ylow) / 2.0
        else
            powint = ylow + (yhigh-ylow) * (((xval-xlow)/(xhigh-xlow))**pow)
        end if
    end function powint


    subroutine locate ( xx, n, is, ie, x, j )
    !-----------------------------------------------------------------------
    ! Given an array "xx" of length "n", and given a value "x", this routine
    ! returns a value "j" such that "x" is between xx(j) and xx(j+1).  xx
    ! must be monotonic, either increasing or decreasing.  j=is-1 or j=ie is
    ! returned to indicate that x is out of range.
    ! Modified from GSLIB version (Deutsch, C.V., 1992)
    !-----------------------------------------------------------------------
        integer, intent (in)    :: n,is,ie
        integer, intent (inout) :: j
        real*8,  pointer    :: xx(:)
        real*8,  intent (inout) :: x
        integer :: jl, ju, jm

        ! Initialize lower and upper methods:
        jl = lbound(xx, 1)
        if(is < jl .or. ie > n) &
            STOP 'ERROR IN LOCATE ROUTINE: LOWER BOUNDS EXCEEDED'
        
        jl = is - 1
        ju = ie
        if(xx(n) <= x) then
            j = ie ; return
        end if

        ! If we are not done then compute a midpoint:
        do while (ju - jl > 1)
            jm = (ju + jl) / 2
            if((xx(ie) > xx(is)).eqv.(x > xx(jm))) then
                jl = jm ; else
                ju = jm ; endif
        enddo
        j = jl
    end subroutine locate
 
end module normaldist