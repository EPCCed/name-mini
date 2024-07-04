module pa_aos0_gpu

#include "assertion.h"

  use assertion
  use constants
  use maths
  use pa_aos0
  use storage_aos0
  use source
  use util_omp_lib

  implicit none
  private

  ! Particle AOS storage management intended for device use

  integer, parameter :: NAME_PA_BLOCKSZ = 128

  type, public, extends(pa_aos0_t) :: pa_aos0_gpu_t
    private
    integer (int32)               :: nmaxLocal        ! particles allocated
    integer (int32)               :: nBlocks          ! of size NAME_PA_BLOCKSZ
    integer (int32)               :: lastBlock        ! last active block
    integer (int32)               :: lastActive       ! last active particle
    integer (int32), allocatable  :: nActiveBlock(:)  ! per block active count
    integer (int64)               :: nReleased        ! single source
  contains
    procedure :: printState        => pa_aos0_gpu_printState
    procedure :: releaseFromSource => pa_aos0_gpu_releaseFromSource
    procedure :: removeInactive    => pa_aos0_gpu_removeInactive
    procedure :: lastParticle      => pa_aos0_gpu_lastParticle
    procedure :: nameScheduleSize  => pa_aos0_gpu_nameScheduleSize
    procedure :: destroy           => pa_aos0_gpu_release
  end type pa_aos0_gpu_t

  public :: pa_aos0_gpu_create_pointer
  public :: pa_aos0_gpu_create
  public :: pa_aos0_gpu_release

  public :: pa_aos0_gpu_releaseFromSource
  public :: pa_aos0_gpu_printState
  public :: pa_aos0_gpu_lastParticle
  public :: pa_aos0_gpu_nameScheduleSize

contains

  !---------------------------------------------------------------------------

  function pa_aos0_gpu_create_pointer(maxLocal) result(pa)

    ! Return a pointer to a new object

    integer,               intent(in) :: maxLocal
    class (pa_aos0_gpu_t), pointer    :: pa

    allocate(pa)
    call pa_aos0_gpu_create(pa, maxLocal)

  end function pa_aos0_gpu_create_pointer

  !---------------------------------------------------------------------------

  subroutine pa_aos0_gpu_create(pa, maxLocal)

    ! Instantiate an object

    type (pa_aos0_gpu_t), intent(out) :: pa
    integer (int32),      intent(in)  :: maxLocal

    ! Storage: allocate a whole number of blocks
    ! Furthermore, for the purposes of sorting, allocate 2^n blocks

    pa%nBlocks   = mathsNearestPowerOfTwo(1 + maxLocal/NAME_PA_BLOCKSZ)
    pa%nmaxLocal = pa%nBlocks*NAME_PA_BLOCKSZ

    allocate(pa%particles(pa%nmaxLocal))
    allocate(pa%nActiveBlock(pa%nBlocks))

    ! Management
    pa%lastBlock  = 0
    pa%lastActive = 0
    pa%nReleased  = 0

    !$omp target enter data map(always, to: pa)
    !$omp target enter data map(alloc: pa%particles)
    !$omp target enter data map(alloc: pa%nActiveBlock)

    ! Ensure all particles are inactive ...
    block
      integer :: iP

      !$omp target teams num_teams(pa%nBlocks)
      !$omp distribute parallel do num_threads(NAME_PA_BLOCKSZ)
      do iP = 1, pa%nmaxLocal
        pa%particles(iP)%active = .false.
      end do
      !$omp end distribute parallel do
      !$omp end target teams
    end block

  end subroutine pa_aos0_gpu_create

  !---------------------------------------------------------------------------

  subroutine pa_aos0_gpu_release(self)

    ! Release resources

    class (pa_aos0_gpu_t), intent(inout) :: self

    !$omp target exit data map(delete: self%nActiveBlock)
    !$omp target exit data map(delete: self%particles)
    !$omp target exit data map(delete: self)

    deallocate(self%nActiveBlock)
    deallocate(self%particles)

    self%lastActive = 0
    self%nReleased  = 0
    self%lastBlock  = 0
    self%nBlocks    = 0
    self%nmaxLocal  = 0

  end subroutine pa_aos0_gpu_release

  !---------------------------------------------------------------------------

  subroutine pa_aos0_gpu_releaseFromSource(self, nRelease, source)

    ! Release nRelease particles from source

    class (pa_aos0_gpu_t), intent(inout) :: self
    integer,               intent(in)    :: nRelease
    type (source_t),       intent(in)    :: source

    integer (int32) :: iR
    integer (int64) :: iUP
    integer (int32) :: iP

    ! If nRelease << NAME_PA_BLOCKSZ, blocks can get used up rather
    ! quickly. Better to offset by just lastParticle?

    ! Check we have enough capacity to store nRelease new particles

    if (self%lastParticle() + nRelease > self%nmaxLocal) then
      print *, "Attempt to release too many particles"
      print *, "nRelease:                            ", nRelease
      print *, "lastParticle()                       ", self%lastParticle()
      print *, "nmaxLocal                            ", self%nmaxLocal
      print *, "nBlocks                              ", self%nBlocks
      print *, "Block size                           ", NAME_PA_BLOCKSZ
      stop     "Fatal"
    end if

    !$omp target teams map(to: nRelease) thread_limit(NAME_PA_BLOCKSZ)
    !$omp distribute parallel do private(iR, iUP, iP)
    do iR = 1, nRelease

      iUP = source%iUP0 + self%nReleased + iR
      iP  = self%lastActive + iR
      self%particles(iP) = storage_aos0_initialiseFromSource(iUP, source)

    end do
    !$omp end distribute parallel do
    !$omp end target teams

    self%lastActive = self%lastActive + nRelease
    !$omp target update to(self%lastActive)

    self%lastBlock = (NAME_PA_BLOCKSZ + self%lastActive - 1)/NAME_PA_BLOCKSZ
    !$omp target update to(self%lastBlock)

    ! Total particles released for one source for all time

    self%nReleased = self%nReleased + nRelease
    !$omp target update to(self%nReleased)


  end subroutine pa_aos0_gpu_releaseFromSource

  !---------------------------------------------------------------------------

  subroutine pa_aos0_gpu_removeInactive(self)

    class (pa_aos0_gpu_t), intent(inout) :: self

    ! A two-stage process.
    ! 1. Intra-block sort for each block
    ! 2. A merge between blocks

    call pa_aos0_gpu_update_nactiveblock(self)
    call pa_aos0_gpu_sort_ordered(self)
    call pa_aos0_gpu_update_nactiveblock(self)
    call pa_aos0_gpu_merge(self)

    ! could be streamlined ...
    ! could print some information ...

    self%lastActive = self%nActive()
    !$omp target update to(self%lastActive)

    self%lastBlock = (NAME_PA_BLOCKSZ + self%lastActive - 1)/NAME_PA_BLOCKSZ
    !$omp target update to(self%lastBlock)

  end subroutine pa_aos0_gpu_removeInactive

  !----------------------------------------------------------------------------

  subroutine pa_aos0_gpu_update_nactiveblock(self)

    ! Set number of active entries per block

    class (pa_aos0_gpu_t), intent(inout) :: self

    integer :: bid
    integer :: na, iL

    !$omp target teams num_teams(self%nBlocks) private(bid, iL, na)
    bid = omp_get_team_num()
    na = 0
    do iL = 1, NAME_PA_BLOCKSZ
      if (self%particles(bid*NAME_PA_BLOCKSZ + iL)%active) na = na + 1
    end do
    self%nActiveBlock(1 + bid) = na
    !$omp end target teams

  end subroutine pa_aos0_gpu_update_nactiveblock

  !---------------------------------------------------------------------------

  subroutine pa_aos0_gpu_sort_ordered(self)

    ! Block-level sort (preserves order). Serial. Slow.
    ! A truly parallel block-level sort is pending implementation.

    class (pa_aos0_gpu_t), intent(inout) :: self

    integer :: bid, iPLast, iL, iP, ioffset

    !$omp target teams num_teams(self%nBlocks) thread_limit(1) &
    !$omp              private(bid, ioffset, iL)

    bid     = omp_get_team_num()
    ioffset = bid*NAME_PA_BLOCKSZ

    do iL = 1, NAME_PA_BLOCKSZ
      if (self%particles(ioffset + iL)%active) cycle
      ! Find next active and move to iL
      do iPLast = iL + 1, NAME_PA_BLOCKSZ
         if (.not. self%particles(ioffset + iPlast)%active) cycle
         self%particles(ioffset + iL) = self%particles(ioffset + iPLast)
         self%particles(ioffset + iPLast)%active = .false.
         exit
      end do
    end do

    !$omp end target teams

  end subroutine pa_aos0_gpu_sort_ordered

  !---------------------------------------------------------------------------

  subroutine pa_aos0_gpu_merge(self)

    ! Heirarchical merge for arb. nPart

    class (pa_aos0_gpu_t), intent(inout) :: self

    integer :: nteams
    integer :: istr
    integer :: iteam, iact, irank
    integer :: isrc0, idst0
    integer :: iblkdst, iblksrc
    integer :: nSpaces, nCopy
    integer :: iL

    assert(mathsIsPowerOfTwo(self%nBlocks) .eqv. .true.)

    ! Number of teams and "stage" or "stride" counter
    nteams  = self%nBlocks / 2
    istr    = 1

    if (nteams < 1) return

    do
      !$omp target teams num_teams(nteams) thread_limit(NAME_PA_BLOCKSZ)

      iteam = omp_get_team_num()
      ! acting team = team/istr     0, 1, ...
      ! acting block size istr*nBlocksize
      ! acting rank
      iact   = iteam/istr
      irank  = mod(iteam, istr)

      ! Index into nActiveBlock(0:nBlocks-1)
      iblkdst   = (2*iact    )*istr
      iblksrc   = (2*iact + 1)*istr

      ! spaces in the destination block
      ! number to copy from source block ...
      nSpaces = istr*NAME_PA_BLOCKSZ - self%nActiveBlock(1 + iblkdst)
      nCopy   = min(nSpaces, self%nActiveBlock(1 + iblksrc))

      ! Offsets into active(1:nPart)
      idst0  = iblkdst*NAME_PA_BLOCKSZ + self%nActiveBlock(1 + iblkdst)
      isrc0  = iblksrc*NAME_PA_BLOCKSZ &
             + max(0, self%nActiveBlock(1 + iblksrc) - nCopy)

      ! Parallel region to merge using threads...

      !$omp parallel private(iL)
      iL = 1 + irank*omp_get_num_threads() + omp_get_thread_num()
      if (iL <= nCopy) then
        assert(.not. self%particles(idst0 + iL)%active)
        assert(      self%particles(isrc0 + iL)%active)
        self%particles(idst0 + iL) = self%particles(isrc0 + iL)
        self%particles(isrc0 + iL)%active = .false.
      end if
      !$omp end parallel

      !$omp end target teams


      ! Teams/blocks must synchronise before update to nActiveBlock(:)
      ! There could be a more handy way to do this which does not
      ! involve a separate kernel...

      !$omp target teams num_teams(nteams) thread_limit(1)

      iteam = omp_get_team_num()
      iact   = iteam/istr
      irank  = mod(iteam, istr)

      ! Index into nActiveBlock(0:nBlocks-1)
      iblkdst   = (2*iact    )*istr
      iblksrc   = (2*iact + 1)*istr

      ! Update active block count (only once, at acting root)
      if (irank == 0) then
        self%nActiveBlock(1 + iblkdst) = self%nActiveBlock(1 + iblkdst) &
                                       + self%nActiveBlock(1 + iblksrc)
      end if
      !$omp end target teams

      istr = istr*2
      if (istr >= self%nBlocks) exit

    end do

  end subroutine pa_aos0_gpu_merge

  !---------------------------------------------------------------------------

  subroutine pa_aos0_gpu_printState(self)

    ! Produce a report on the current state of the management

    class (pa_aos0_gpu_t), intent(in) :: self

    integer :: iP
    integer :: nactive

    ! There should be no active particles beyond the last block

    nactive = 0

    !$omp target enter data map(always, to: nactive)

    !$omp target teams reduction(+: nactive)
    !$omp distribute parallel do private(iP) reduction(+: nactive)
    do iP = 1, self%lastBlock*NAME_PA_BLOCKSZ
      if (self%particles(iP)%active) nactive = nactive + 1
    end do
    !$omp end distribute parallel do
    !$omp end target teams

    !$omp target exit data map(from: nactive)

    print *, "Manager: pa_aos0_gpu               "

    print *, "Manager: NAME_PA_BLOCSZ         ", NAME_PA_BLOCKSZ
    print *, "Manager: maxLocal               ", self%nmaxLocal
    print *, "Manager: nBlocks                ", self%nBlocks
    print *, "Manager: omp_get_num_devices()  ", omp_get_num_devices()

    print *, "Manager: nactive                ", nactive
    print *, "Manager: lastBlock              ", self%lastBlock
    print *, "Manager: lastActive             ", self%lastActive
    print *, "Manager: lastParticle()         ", self%lastParticle()

  end subroutine pa_aos0_gpu_printState

  !---------------------------------------------------------------------------

  pure function pa_aos0_gpu_lastParticle(self) result(last)

    class (pa_aos0_gpu_t), intent(in) :: self
    integer                        :: last

    last = self%lastBlock*NAME_PA_BLOCKSZ

  end function pa_aos0_gpu_lastParticle

  !---------------------------------------------------------------------------

  function pa_aos0_gpu_nameScheduleSize(self) result(nScheduleSize)

    ! The static schedule size is 1 for the GPU

    class (pa_aos0_gpu_t), intent(in) :: self
    integer                        :: nScheduleSize

    nScheduleSize = 1

  end function pa_aos0_gpu_nameScheduleSize

end module pa_aos0_gpu
