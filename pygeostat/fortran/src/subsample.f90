! public
      subroutine subsample(datafl,col,ncell,nsub,nreal,rseed,subsamp)
!-----------------------------------------------------------------------
!
!                  Sub-sample Realizations of a Dataset
!                  ************************************
!
!  Description
!  -----------
!  When datasets become too large to efficiently use, a sub-sample may
!  be extracted and used in lieu. A Fisher-Yates shuffle is implemented.
!  This subroutine was designed with the intent of wrapping it for use
!  in python. Motivated from Jared L. Deutsch's use of Fisher-Yates
!  shuffle in gslib program histpltsim.
!
!  Assumes that the data file being read is in GSLIB format. Outputs an
!  1D array in the case that only one realization is sub-sampled, or a
!  2D array in the case where multiple realizations are sub-sampled.
!  Each realizations sub-sample is within a unique column.
!
!  Compile using the command:
!  f2py -c -m --compiler=cygwin subsample subsample.f90 random.f90
!
!  Parameters
!  ----------
!  datafl:  A single input datafile with the variable being sampled
!  col:     Column containing data to sub-sample
!  ncell:   The number of cells within a single realization, can also be
!           interpreted as the number of data to read, then sub-sample
!  nsub:    Number of sub-samples to output
!  nreal:   Number of realizations to sub-sample
!  rseed:   A seed value to pass to the random number generator
!  subsamp: Output array with a realization in each column
!
!
!  Contributors
!  ------------
!  Author  - Warren E. Black                    DATE: September 22, 2015
!
!  Revised - Name                               DATE:
!            Description of revisions.
!
!  (c) 2015, Warren E. Black
!-----------------------------------------------------------------------

      use random

      implicit none

      integer, intent(in)   :: col, ncell, nsub, nreal, rseed
      character, intent(in) :: datafl*512
      integer, parameter    :: dp=kind(0.d0)
      real(dp), intent(out) :: subsamp(nsub,nreal)
      real(dp)              :: s_temp, temp(col)
      real(dp), allocatable :: dat(:)
      integer               :: s, test, i, j, k, nvar

      !Initialize the random number generator
      call init_genrand(rseed)

      !Read over the header
      open(20,file=datafl,action='read',status='old')
      read(20,*)
      read(20,*) nvar
      do i=1, nvar
        read(20,*)
      end do

      !Loop over all realizations
      do i=1,nreal
        !Allocate the required arrays
        allocate(dat(ncell),stat = test)
        if (test.ne.0) stop 'Error:  Allocation failed'
        !Read in the data
        do j=1,ncell
          read(20,*) (temp(k),k=1,col)
          dat(j)=temp(col)
        end do
        !Shuffle the required number of subsamples
        do k=ncell, (ncell-nsub+1),-1
            !Draw a random integer
            s = max(1,int(real(grnd())*real(k)))
            !Swap the required subsample
            s_temp = dat(k)
            dat(k) = dat(s)
            dat(s) = s_temp
            !Save the swapped value to the output array
            j=ncell-k+1
            subsamp(j,i)=dat(k)
        end do
        !Dump the saved data
        deallocate(dat)
      end do

      !Finish up
      close(20)
      return
      end subroutine subsample
