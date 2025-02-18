!  solve.f90 
!
!  FUNCTIONS/SUBROUTINES exported from solve.dll:
!  solve      - AX = B
!
subroutine solve ( A, X, B, N, M, INFO )
    implicit none
  ! Variables
    integer, intent (in) :: N, M !Rank and Solutions

    real*8, dimension(N,N), intent (inout) :: A
    real*8, dimension(N,M), intent (inout) :: X
    real*8, dimension(N,M), intent (in   ) :: B
    integer, dimension(N) :: IPIV
    integer, intent (out) :: INFO !return state of the program

    integer :: LDA, LDB !Leading dimension of A and B
    integer :: NRHS !number of righthandsides

    ! Body of solve
    LDA = N
    LDB = N
    NRHS = M
    X = B
    CALL DGESV( N, NRHS, A, LDA, IPIV, X, LDB, INFO ) 
end subroutine solve


subroutine invert ( A, N )
    implicit none
  ! Variables
    integer, intent (in) :: N !Rank and Solutions

    real*8, dimension(N,N), intent (inout) :: A
    integer, dimension(N) :: IPIV

    integer :: LDA !Leading dimension of A and B
    integer :: INFO  !return state of the program
    integer :: LWORK !length of array work
  
    real*8, allocatable, dimension(:) :: WORK

    ! Body of solve
    LDA = N
    LWORK = 3*N-1     !Set size of WORK array
    allocate(WORK(LWORK))
    WORK = 0.0        !Zero the work array
    INFO = 0
    CALL DGETRF( N, N, A, LDA, IPIV, INFO )
    CALL DGETRI( N, A, LDA, IPIV, WORK, LWORK, INFO )
    return
end subroutine invert