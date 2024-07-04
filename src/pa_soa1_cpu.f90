module pa_soa1_cpu

#include "assertion.h"

  use assertion
  use constants
  use pa_soa1
  use storage_soa1
  use source
  use util_omp_lib

  implicit none
  private

  ! See manager_soa0_cpu.f90 for a description
  ! There is a considerable amount of repeated code at the moment

  integer, parameter :: NAME_PA_BLOCKSZ = 256

  type, public, extends(pa_soa1_t) :: pa_soa1_cpu_t
    private
    integer (int32)               :: nThreads         ! omp_get_max_threads()
    integer (int32)               :: nBlocks          ! allocated
    integer (int32)               :: nmaxLocal        ! particles allocated
    integer (int32)               :: lastBlock        ! last active block
    integer (int32), allocatable  :: last(:)          ! per thread
  contains
    procedure :: printState        => pa_soa1_cpu_printState
    procedure :: releaseFromSource => pa_soa1_cpu_releaseFromSource
    procedure :: removeInactive    => pa_soa1_cpu_removeInactive
    procedure :: lastParticle      => pa_soa1_cpu_lastParticle
    procedure :: nameScheduleSize  => pa_soa1_cpu_nameScheduleSize
    procedure :: destroy           => pa_soa1_cpu_release
  end type pa_soa1_cpu_t

  public :: pa_soa1_cpu_create_pointer
  public :: pa_soa1_cpu_create
  public :: pa_soa1_cpu_release

  public :: pa_soa1_cpu_printState
  public :: pa_soa1_cpu_releaseFromSource
  public :: pa_soa1_cpu_removeInactive
  public :: pa_soa1_cpu_lastParticle
  public :: pa_soa1_cpu_nameScheduleSize

contains

  !---------------------------------------------------------------------------

  function pa_soa1_cpu_create_pointer(maxLocal) result(pa)

    ! Return pointer to a new object

    integer,               intent(in) :: maxLocal
    class (pa_soa1_cpu_t), pointer    :: pa

    allocate(pa)
    call pa_soa1_cpu_create(pa, maxLocal)

  end function pa_soa1_cpu_create_pointer

  !---------------------------------------------------------------------------

  subroutine pa_soa1_cpu_create(pa, maxLocal)

    ! Instantiate the manager object
    ! As this is cpu-only, there is no mapping to device memory

    type (pa_soa1_cpu_t), intent(out) :: pa
    integer (int32),      intent(in)  :: maxLocal

    integer :: nGroup

    ! allocate a whole number of large blocks worth ...

    pa%nThreads  = omp_get_max_threads()

    nGroup       = 1 + (maxLocal - 1)/(pa%nThreads*NAME_PA_BLOCKSZ)
    pa%nBlocks   = nGroup*pa%nThreads
    pa%nmaxLocal = pa%nBlocks*NAME_PA_BLOCKSZ
    assert(pa%nmaxLocal >= maxLocal)

    ! Storage (no GPU mapping)
    call storage_soa1_create(pa%particles, pa%nmaxLocal, map = .false.)

    ! Management
    pa%lastBlock = 0
    allocate(pa%last(0:pa%nThreads-1))

    block
      integer :: tid
      do tid = 0, pa%nThreads-1
        pa%last(tid) = 0
      end do
    end block

    ! Ensure all particles are inactive ...
    block
      integer :: iP

      !$omp target teams
      !$omp distribute parallel do schedule(static, NAME_PA_BLOCKSZ)
      do iP = 1, pa%nmaxLocal
        pa%particles%active(iP) = .false.
      end do
      !$omp end distribute parallel do
      !$omp end target teams
    end block

  end subroutine pa_soa1_cpu_create

  !---------------------------------------------------------------------------

  subroutine pa_soa1_cpu_release(self)

    ! Release resources (no unmapping from device memory)

    class (pa_soa1_cpu_t), intent(inout) :: self

    call storage_soa1_release(self%particles, unmap = .false.)
    deallocate(self%last)

    self%lastBlock = 0
    self%nmaxLocal = 0

  end subroutine pa_soa1_cpu_release

  !---------------------------------------------------------------------------

  subroutine pa_soa1_cpu_releaseFromSource(self, nRelease, source)

    ! Release nRelease particles from source

    class (pa_soa1_cpu_t), intent(inout) :: self
    integer,               intent(in)    :: nRelease
    type (source_t),       intent(in)    :: source

    integer (int32) :: iR
    integer (int64) :: iUP
    integer (int32) :: iP

    if (self%lastBlock*NAME_PA_BLOCKSZ + nRelease > self%nmaxLocal) then
      print *, "Attempt to release too many particles"
      print *, "nRelease:                            ", nRelease
      print *, "lastBlock                            ", self%lastBlock
      print *, "nmaxLocal                            ", self%nmaxLocal
      stop     "Fatal"
    end if

    ! Release particles
    ! We want to add roughly even numbers per thread, so default schedule

    !$omp parallel do private(iR, iUP, iP)
    do iR = 1, nRelease

      iUP = source%iUP0 + iR                 ! Revisit as iUP0 needs updating
      iP  = pa_soa1_cpu_next_index(self)
      call storage_soa1_initialiseFromSource(self%particles, iP, iUP, source)

    end do
    !$omp end parallel do

  end subroutine pa_soa1_cpu_releaseFromSource

  !---------------------------------------------------------------------------

  subroutine pa_soa1_cpu_removeInactive(self)

    class (pa_soa1_cpu_t), intent(inout) :: self

    ! Move active particles to the left, and inactive particles to the
    ! right in the array on each thread.

    ! This effectively clears out inactive particles, and we can add new
    ! particles at the first free position again without synchronisation.

    integer :: tid
    integer :: iplast
    integer :: iP
    integer :: lastblock

    ! For each thread, march forward from the front, and look for holes.
    ! At the same time, perform a backwards search from the end for actives.
    ! Fill the holes at the start with active particles from the end
    ! and mark the latter storage location inactive.

    lastblock = 0

    !$omp parallel private(tid, iplast, iP) reduction(max: lastblock)

    tid = omp_get_thread_num()

    !$omp do schedule (static, NAME_PA_BLOCKSZ)
    do iP = 1, self%lastParticle()
      iplast = pa_soa1_cpu_last_to_index(self, tid)
      if (iP >= iplast) exit
      if (self%particles%active(iP)) cycle
      ! Copy iplast to iP
      call storage_soa1_copy(self%particles, iP, iplast)
      self%particles%active(iplast) = .false.
      do
        self%last(tid) = self%last(tid) - 1
        iplast = pa_soa1_cpu_last_to_index(self, tid)
        if (self%particles%active(iplast)) exit
      end do
    end do
    !$omp end do

    lastblock = max(lastblock, (iplast - 1)/NAME_PA_BLOCKSZ)

    !$omp end parallel

    self%lastBlock = lastblock

  end subroutine pa_soa1_cpu_removeInactive

  !---------------------------------------------------------------------------

  subroutine pa_soa1_cpu_printState(self)

    ! Provide some information on the current state of the management

    class (pa_soa1_cpu_t), intent(in) :: self

    integer, allocatable :: npart(:)
    integer, allocatable :: nlast(:)
    integer              :: iP, tid

    allocate(npart(0:omp_get_max_threads()))
    allocate(nlast(0:omp_get_max_threads()))

    npart(:) = 0
    nlast(:) = 0

    !$omp parallel do private(iP, tid) schedule(static, NAME_PA_BLOCKSZ)
    do iP = 1, self%nmaxLocal
      tid = omp_get_thread_num()
      if (self%particles%active(iP)) then
         npart(tid) = npart(tid) + 1
         nlast(tid) = iP
      end if
    end do

    print *, "Manager: pa_soa1_cpu"
    print *, "Manager: NAME_PA_BLOCSZ: ", NAME_PA_BLOCKSZ
    print *, "Manager: nThreads:       ", self%nThreads
    print *, "Manager: nBlocks:        ", self%nBlocks
    print *, "Manager: maxLocal:       ", self%nmaxLocal
    print *, "Manager: total active    ", sum(npart(:))
    print *, "Manager: lastBlock       ", self%lastblock

    print *, ""
    print *, "Thread, active particles, last particle"
    print *, "------------------------------------------"
    do tid = 0, omp_get_max_threads() - 1
      write (*, '(i6,i16,i16)') tid, npart(tid), nlast(tid)
    end do

    deallocate(nlast)
    deallocate(npart)

  end subroutine pa_soa1_cpu_printState

  !---------------------------------------------------------------------------

  pure function pa_soa1_cpu_lastParticle(self) result(last)

    class (pa_soa1_cpu_t), intent(in) :: self
    integer                           :: last

    ! blocks count from zero ...
    last = (self%lastBlock + 1)*NAME_PA_BLOCKSZ

  end function pa_soa1_cpu_lastParticle

  !---------------------------------------------------------------------------

  function pa_soa1_cpu_nameScheduleSize(self) result(nScheduleSize)

    ! Is the block size for cpu threads ...

    class (pa_soa1_cpu_t), intent(in) :: self
    integer                        :: nScheduleSize

    nScheduleSize = NAME_PA_BLOCKSZ

  end function pa_soa1_cpu_nameScheduleSize

  !---------------------------------------------------------------------------
  !---------------------------------------------------------------------------

  function pa_soa1_cpu_next_index(self) result(iP)

    type (pa_soa1_cpu_t), intent(inout)  :: self
    integer (int32)                   :: iP

    integer :: tid
    integer :: iB
    integer :: iL

    assert(allocated(self%last))

    ! Add one particle on the current thread

    tid = omp_get_thread_num()

    self%last(tid) = self%last(tid) + 1

    iP = self%last(tid)

    ! Update the last block record if we have landed in a new block
    ! if iB > lastBlock, need to update e.g. "omp atomic update, compare",
    ! or ...

    iB = (iP - 1)/NAME_PA_BLOCKSZ
    iB = iB*self%nThreads + tid

    iL = mod(iP - 1, NAME_PA_BLOCKSZ)         ! 0 ... NAME_PA_BLOCKSZ - 1
    if (iL == 0) then
      !$omp atomic
      self%lastBlock = max(self%lastBlock, iB)
    end if

    ! The is the index in the particle array...

    iP = 1 + iL + iB*NAME_PA_BLOCKSZ

    assert(1 <= iP .and. iP <= self%nmaxLocal)

  end function pa_soa1_cpu_next_index

  !---------------------------------------------------------------------------

  pure function pa_soa1_cpu_last_to_index(self, tid) result(iP)

    ! The last(tid) entry is the last particle per thread; translate
    ! to the corresponding array index allowing for the block structure.

    ! A thread may have no particles, in which case iP = 0

    class (pa_soa1_cpu_t), intent(in) :: self
    integer,               intent(in) :: tid
    integer                           :: iP

    integer :: iL
    integer :: iBlock

    assert(0 <= tid .and. tid < self%nThreads)

    iL     = mod(self%last(tid) - 1, NAME_PA_BLOCKSZ)
    iBlock = pa_soa1_cpu_last_to_block(self, tid)
    iP     = 1 + iL + iBlock*NAME_PA_BLOCKSZ

    assert(0 <= iP .and. iP <= self%nmaxLocal)

  end function pa_soa1_cpu_last_to_index

  !----------------------------------------------------------------------------

  pure function pa_soa1_cpu_index_to_last(self, iP) result(last)

    ! For array index iP, what is the corresponding "last" value
    ! in the block structure. This only depends on the block size
    ! and the total number of threads.

    class (pa_soa1_cpu_t), intent(in) :: self
    integer,               intent(in) :: iP
    integer                           :: last

    assert(1 <= iP  .and. iP  <= self%nmaxLocal)

    block
      integer :: iBlock
      integer :: iL

      iBlock = (iP - 1) / (NAME_PA_BLOCKSZ*self%nThreads)
      iL     = mod(iP - 1, NAME_PA_BLOCKSZ)
      last   = 1 + iL + iBlock*NAME_PA_BLOCKSZ
    end block

    assert(1 <= last .and. last <= (self%nmaxLocal/self%nThreads))

  end function pa_soa1_cpu_index_to_last

  !---------------------------------------------------------------------------

  pure function pa_soa1_cpu_last_to_block(self, tid) result(iblock)

    ! What is the block index (in the global picture) of the last
    ! particle for thread tid?

    class (pa_soa1_cpu_t), intent(in) :: self
    integer,               intent(in) :: tid
    integer                        :: iblock

    assert(0 <= tid .and. tid < self%nThreads)

    iblock = ((self%last(tid) - 1)/NAME_PA_BLOCKSZ)*self%nThreads + tid

    assert(0 <= iblock .and. iblock < self%nBlocks)

  end function pa_soa1_cpu_last_to_block

  !---------------------------------------------------------------------------

  subroutine pa_soa1_cpu_update_last(self)

    ! Set the last(:) entries and lastBlock from the current
    ! particle state.

    class (pa_soa1_cpu_t), intent(inout) :: self

    integer :: iP
    integer :: lastBlock
    integer :: tid

    lastBlock = 0

    !$omp parallel private(iP, tid) reduction(max: lastBlock)
    tid = omp_get_thread_num()
    !$omp do schedule(static, NAME_PA_BLOCKSZ)
    do iP = 1, self%nmaxLocal
      if (.not. self%particles%active(iP)) cycle
      self%last(tid) = pa_soa1_cpu_index_to_last(self, iP)
      lastBlock      = pa_soa1_cpu_last_to_block(self, tid)
    end do
    !$omp end do
    !$omp end parallel

    self%lastBlock = lastBlock

  end subroutine pa_soa1_cpu_update_last

end module pa_soa1_cpu
