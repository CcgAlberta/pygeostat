!  eig.f90 
!
!  FUNCTIONS/SUBROUTINES exported from eig.dll:
!  eig      - subroutine 
!
subroutine eig ( A, E, V, N )
    implicit none
  ! Variables
  integer, intent (in) :: N
  
  real*8, dimension(N,N), intent (in) :: A !Input matrix
  real*8, dimension( N ), intent (out) :: E  !Eigenvalues
  real*8, dimension(N,N), intent (out) :: V  !Eigenvectors
!!!  
  character*1, parameter :: JOBS = 'V' !Always calculate eigenvalues and vectors
  character*1, parameter :: UPLO = 'U' !Always store the upper triangle of A
  
  integer :: LDA !leading dimension of A ( = N )
  integer :: LWORK !length of array work
  integer :: INFO  !return state of the program
  
  real*8, allocatable, dimension(:) :: WORK

  ! Body of eig
  LDA = N           !Set LDA
  E = 0.0           !Zero the eigenvalue array (W in dsyev.f)
  LWORK = 3*N-1     !Set size of WORK array
  allocate(WORK(LWORK))
  WORK = 0.0        !Zero the work array
  INFO = 0
  
  !Call the eigenvalue decomposition routine
  V = A
  call DSYEV( JOBS, UPLO, N, V, LDA, E, WORK, LWORK, INFO )
  deallocate(WORK)
end subroutine eig
