module varsim_wrapper
    use varsim_wrap
    implicit none
contains
    subroutine callvarsim(parfl)
        character*(*), intent(in) :: parfl
        call runvarsim(parfl)
    end subroutine
end module varsim_wrapper
    
