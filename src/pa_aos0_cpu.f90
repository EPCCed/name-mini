module pa_aos0_cpu

#include "assertion.h"

  use assertion
  use constants
  use util_omp_lib
  use pa_aos0
  use storage_aos0
  use source

  implicit none
  private

  ! Particle AOS storage managed by block intended for host use only

  ! Structure
  ! We will use a number of blocks each of size NAME_PA_BLOCKSZ
  ! with blocks indexed 0, 1, ..., nblocks - 1
  !
  ! We take account of the maximum number of threads at run
  ! time omp_get_max_threads() tid = 0, omp_get_max_threads() - 1.
  !
  ! We allocate space for nGroups*nThreads*NAME_PA_BLOCKSZ particles
  ! where nGroups is large enough to accommodate the total number of
  ! particles requested. For the purposes of sorting, we also want a
  ! power-of-two number of blocks.

  ! The management of memory between CPU threads is opaque.

  ! The object life cycle is:
  !
  !   call aos0_cpu_create(pa, maxParticles)
  !
  !   call aos0_cpu_release(pa)

  ! The block structure allows clients to use something of the form
  !
  ! !$omp parallel do schedule (static, NAME_PA_BLOCKSZ)
  ! do iP = 1, manager%lastParticle()
  !   ...

  ! One would have less complexity here by using, e.g.:
  !
  ! !$omp parallel
  ! do iP = manager%firstParticle(), manager%lastParticle()
  !
  ! which would simplify the block structure. NOT SURE HOW THAT WORKS
  ! in "target teams" PICTURE.

  ! The "last(:)" entry per thread is a little opaque. It is the
  ! number of the last active particle relevant for this thread
  ! where "number" would be the index if only thread blocks were
  ! considered.
  ! E.g., a last active particle at position 1 in the second block
  ! for a given thread would have last = NAME_PA_BLOCKSZ + 1

  integer, parameter :: NAME_PA_BLOCKSZ = 256

  type, public, extends(pa_aos0_t) :: pa_aos0_cpu_t
    private
    integer (int32)               :: nThreads         ! omp_get_max_threads()
    integer (int32)               :: nBlocks          ! allocated
    integer (int32)               :: nmaxLocal        ! particles allocated
    integer (int32)               :: lastBlock        ! last active block
    integer (int32), allocatable  :: last(:)          ! per thread
    integer (int64)               :: nReleased        ! single source
    integer (int32)               :: iRelease         ! episode counter
  contains
    procedure :: printState        => pa_aos0_cpu_printState
    procedure :: releaseFromSource => pa_aos0_cpu_releaseFromSource
    procedure :: removeInactive    => pa_aos0_cpu_removeInactive
    procedure :: lastParticle      => pa_aos0_cpu_lastParticle
    procedure :: nameScheduleSize  => pa_aos0_cpu_nameScheduleSize
    procedure :: destroy           => pa_aos0_cpu_release
  end type pa_aos0_cpu_t

  public :: pa_aos0_cpu_create_pointer
  public :: pa_aos0_cpu_create
  public :: pa_aos0_cpu_release

  public :: pa_aos0_cpu_printState
  public :: pa_aos0_cpu_releaseFromSource
  public :: pa_aos0_cpu_removeInactive
  public :: pa_aos0_cpu_lastParticle
  public :: pa_aos0_cpu_nameScheduleSize

  public :: pa_aos0_cpu_self_test

contains

  !---------------------------------------------------------------------------

  function pa_aos0_cpu_create_pointer(maxLocal) result(pa)

    ! Return a pointer to a new object

    integer,               intent(in) :: maxLocal
    class (pa_aos0_cpu_t), pointer    :: pa

    allocate(pa)
    call pa_aos0_cpu_create(pa, maxLocal)

  end function pa_aos0_cpu_create_pointer

  !---------------------------------------------------------------------------

  subroutine pa_aos0_cpu_create(pa, maxLocal)

    ! Instantiate the object

    type (pa_aos0_cpu_t), intent(out) :: pa
    integer (int32),      intent(in)  :: maxLocal

    integer :: nGroup

    ! allocate a whole number of large blocks worth ...

    pa%nThreads = omp_get_max_threads()

    nGroup = 1 + (maxLocal - 1)/(pa%nThreads*NAME_PA_BLOCKSZ)
    pa%nBlocks  = nGroup*pa%nThreads
    pa%nmaxLocal = pa%nBlocks*NAME_PA_BLOCKSZ
    assert(pa%nmaxLocal >= maxLocal)

    ! Note we have, by construction, at least one block per thread ...
    assert(pa%nBlocks >= pa%nThreads)

    ! Storage
    allocate(pa%particles(pa%nmaxLocal))

    ! Management
    pa%lastBlock = 0
    allocate(pa%last(0:pa%nThreads-1))

    block
      integer :: tid
      do tid = 0, pa%nThreads-1
        pa%last(tid) = 0
      end do
    end block

    ! May need to investigate how this looks for first touch ...
    block
      integer :: iP

      !$omp target teams
      !$omp distribute parallel do schedule(static, NAME_PA_BLOCKSZ)
      do iP = 1, pa%nmaxLocal
        pa%particles(iP)%active = .false.
      end do
      !$omp end distribute parallel do
      !$omp end target teams
    end block

    pa%iRelease  = 0
    pa%nReleased = 0

  end subroutine pa_aos0_cpu_create

  !---------------------------------------------------------------------------

  subroutine pa_aos0_cpu_release(self)

    ! Release resources

    class (pa_aos0_cpu_t), intent(inout) :: self

    deallocate(self%particles)
    deallocate(self%last)

    self%lastBlock = 0
    self%nmaxLocal = 0
    self%nReleased = 0

  end subroutine pa_aos0_cpu_release

  !---------------------------------------------------------------------------

  subroutine pa_aos0_cpu_releaseFromSource(self, nRelease, source)

    ! Release nRelease particles from source

    class (pa_aos0_cpu_t), intent(inout) :: self
    integer,               intent(in)    :: nRelease
    type (source_t),       intent(in)    :: source

    integer (int32) :: iR
    integer (int64) :: iUP
    integer (int32) :: iP

    ! Check we have enough capacity to store nRelease new particles

    ! The criterion here is: have all threads got enough capacity?
    ! Based on a static distribution of nRelease particles, we can
    ! work this out...

    block
      integer :: tid
      integer :: nRel ! per thread
      logical :: flag

      flag = .false.
      nRel = 1 + nRelease/self%nThreads

      do tid = 0, self%nThreads - 1
        if (self%last(tid) + nRel > self%nmaxLocal/self%nThreads) flag = .true.
      end do

      if (flag) then
        print *, "Attempt to release too many particles"
        print *, "Total released from source:          ", int(self%nReleased)
        print *, "Currently active:                    ", self%nActive()
        print *, "nRelease:                            ", nRelease
        print *, "nmaxLocal                            ", self%nmaxLocal
        stop     "Fatal"
      end if
    end block

    ! if nRelease <~ nThreads, one might get an unbalanced
    ! distribution of particles between threads, as new
    ! particles will tend to end up in low numbered threads.

    if (nRelease < 2*self%nThreads) then
       ! Expedient for small release balance
       block
         integer :: tid0
         tid0 = mod(self%iRelease, self%nThreads)
         !$omp parallel
         if (omp_get_thread_num() == tid0) then
            do iR = 1, nRelease
               iUP = source%iUP0 + self%nReleased + iR
               iP  = pa_aos0_cpu_next_index(self)
               self%particles(iP) = storage_aos0_initialiseFromSource(iUP, source)
            end do
         end if
         !$omp end parallel
       end block
    else
    !$omp parallel do private(iR, iUP, iP)
    do iR = 1, nRelease

      iUP = source%iUP0 + self%nReleased + iR
      iP  = pa_aos0_cpu_next_index(self)
      self%particles(iP) = storage_aos0_initialiseFromSource(iUP, source)

    end do
    !$omp end parallel do
    end if

    self%iRelease  = self%iRelease  + 1
    self%nReleased = self%nReleased + nRelease

  end subroutine pa_aos0_cpu_releaseFromSource

  !---------------------------------------------------------------------------

  subroutine pa_aos0_cpu_removeInactive(self)

    ! Move active particles to the left, and inactive particles to the
    ! right in the array on each thread.

    ! This effectively clears out inactive particles, and we can add new
    ! particles at the first free position again without synchronisation.

    class (pa_aos0_cpu_t), intent(inout) :: self

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
      iplast = pa_aos0_cpu_last_to_index(self, tid)
      if (iP >= iplast) exit
      if (self%particles(iP)%active) cycle
      ! Move. PENDING Repeat for masses if have masses
      self%particles(iP) = self%particles(iplast)
      self%particles(iplast)%active = .false.
      do
        self%last(tid) = self%last(tid) - 1
        iplast = pa_aos0_cpu_last_to_index(self, tid)
        if (self%particles(iplast)%active) exit
      end do
    end do
    !$omp end do

    lastblock = max(lastblock, (iplast - 1)/NAME_PA_BLOCKSZ)

    !$omp end parallel

    self%lastBlock = lastblock

  end subroutine pa_aos0_cpu_removeInactive

  !---------------------------------------------------------------------------

  subroutine pa_aos0_cpu_printState(self)

    ! Provide some information on the current state of the management

    class (pa_aos0_cpu_t), intent(in) :: self

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
      if (self%particles(iP)%active) then
         npart(tid) = npart(tid) + 1
         nlast(tid) = iP
      end if
    end do

    print *, ""
    print *, "Particle management information"
    print *, "Manager: pa_aos0_cpu"
    print *, "Manager: NAME_PA_BLOCSZ: ", NAME_PA_BLOCKSZ
    print *, "Manager: nThreads:       ", self%nThreads
    print *, "Manager: nBlocks:        ", self%nBlocks
    print *, "Manager: maxLocal:       ", self%nmaxLocal
    print *, "Manager: total released  ", int(self%nReleased)
    print *, "Manager: release episodes", self%iRelease
    print *, "Manager: total active    ", sum(npart(:))
    print *, "Manager: lastBlock       ", self%lastblock

    ! The "High water mark" may want a better name

    print *, ""
    print *, "Per thread information"
    print *, "Thread      active       HWM        global last"
    print *, "-----------------------------------------------"
    do tid = 0, omp_get_max_threads() - 1
      write (*, '(i6,i12,i12,i12)') tid, npart(tid), self%last(tid), nlast(tid)
    end do

    deallocate(nlast)
    deallocate(npart)

  end subroutine pa_aos0_cpu_printState

  !---------------------------------------------------------------------------

  pure function pa_aos0_cpu_lastParticle(self) result(last)

    class (pa_aos0_cpu_t), intent(in) :: self
    integer                        :: last

    ! blocks count from zero ...
    last = (self%lastBlock + 1)*NAME_PA_BLOCKSZ

  end function pa_aos0_cpu_lastParticle

  !---------------------------------------------------------------------------

  function pa_aos0_cpu_nameScheduleSize(self) result(nScheduleSize)

    class (pa_aos0_cpu_t), intent(in) :: self
    integer                        :: nScheduleSize

    nScheduleSize = NAME_PA_BLOCKSZ

  end function pa_aos0_cpu_nameScheduleSize

  !---------------------------------------------------------------------------
  !---------------------------------------------------------------------------

  function pa_aos0_cpu_next_index(self) result(iP)

    type (pa_aos0_cpu_t), intent(inout)  :: self
    integer (int32)                      :: iP

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

  end function pa_aos0_cpu_next_index

  !----------------------------------------------------------------------------

  pure function pa_aos0_cpu_last_to_index(self, tid) result(iP)

    ! The last(tid) entry is the last particle per thread; translate
    ! to the corresponding array index allowing for the block structure.

    ! A thread may have no particles, in which case iP = 0

    class (pa_aos0_cpu_t), intent(in) :: self
    integer,            intent(in) :: tid
    integer                        :: iP

    integer :: iL
    integer :: iBlock

    assert(0 <= tid .and. tid < self%nThreads)

    iL     = mod(self%last(tid) - 1, NAME_PA_BLOCKSZ)
    iBlock = pa_aos0_cpu_last_to_block(self, tid)
    iP     = 1 + iL + iBlock*NAME_PA_BLOCKSZ

    assert(0 <= iP .and. iP <= self%nmaxLocal)

  end function pa_aos0_cpu_last_to_index

  !----------------------------------------------------------------------------

  pure function pa_aos0_cpu_index_to_last(self, iP) result(last)

    ! For array index iP, what is the corresponding "last" value
    ! in the block structure. This only depends on the block size
    ! and the total number of threads.

    class (pa_aos0_cpu_t), intent(in) :: self
    integer,            intent(in) :: iP
    integer                        :: last

    assert(1 <= iP  .and. iP  <= self%nmaxLocal)

    block
      integer :: iBlock
      integer :: iL
      
      iBlock = (iP - 1) / (NAME_PA_BLOCKSZ*self%nThreads)
      iL     = mod(iP - 1, NAME_PA_BLOCKSZ)
      last   = 1 + iL + iBlock*NAME_PA_BLOCKSZ
    end block

    assert(1 <= last .and. last <= (self%nmaxLocal/self%nThreads))

  end function pa_aos0_cpu_index_to_last

  !---------------------------------------------------------------------------

  pure function pa_aos0_cpu_last_to_block(self, tid) result(iblock)

    ! What is the block index (in the global picture) of the last
    ! particle for thread tid?

    class (pa_aos0_cpu_t), intent(in) :: self
    integer,            intent(in) :: tid
    integer                        :: iblock

    assert(0 <= tid .and. tid < self%nThreads)

    iblock = ((self%last(tid) - 1)/NAME_PA_BLOCKSZ)*self%nThreads + tid

    assert(0 <= iblock .and. iblock < self%nBlocks)

  end function pa_aos0_cpu_last_to_block

  !---------------------------------------------------------------------------

  subroutine pa_aos0_cpu_update_last(self)

    ! Set the last(:) entries and lastBlock from the current
    ! particle state.

    class (pa_aos0_cpu_t), intent(inout) :: self

    integer :: iP
    integer :: lastBlock
    integer :: tid

    lastBlock = 0

    !$omp parallel private(iP, tid) reduction(max: lastBlock)
    tid = omp_get_thread_num()
    !$omp do schedule(static, NAME_PA_BLOCKSZ)
    do iP = 1, self%nmaxLocal
      if (.not. self%particles(iP)%active) cycle
      self%last(tid) = pa_aos0_cpu_index_to_last(self, iP)
      lastBlock      = pa_aos0_cpu_last_to_block(self, tid)
    end do
    !$omp end do
    !$omp end parallel

    self%lastBlock = lastBlock

  end subroutine pa_aos0_cpu_update_last

  !---------------------------------------------------------------------------
  !---------------------------------------------------------------------------

  subroutine pa_aos0_cpu_self_test()

    ! Test driver

    call pa_aos0_cpu_self_test_removeInactive()

  end subroutine pa_aos0_cpu_self_test

  !---------------------------------------------------------------------------

  subroutine pa_aos0_cpu_self_test_removeInactive()

    ! Test
    test_block_with_no_inactives: &
    block

      ! Add a certain number of particles with no gaps, remove, and check
      integer, parameter :: nmax = 2000
      integer, parameter :: nactive = nmax/2
      type (pa_aos0_cpu_t)  :: aos0
      integer            :: iP

      call pa_aos0_cpu_create(aos0, nmax)

      !$omp parallel private(iP)
      !$omp do schedule(static, NAME_PA_BLOCKSZ)
      do iP = 1, nmax
        aos0%particles(iP)%active = .false.
        if (ip <= nactive) aos0%particles(iP)%active = .true.
      end do
      !$omp end do
      !$omp end parallel

      call aos0%removeInactive()

      !$omp parallel private(iP)
      !$omp do schedule(static, NAME_PA_BLOCKSZ)
      do iP = 1, nmax
         if (ip <= nactive) then
           assert(aos0%particles(iP)%active .eqv. .true.)
         end if
         if (ip >  nactive) then
           assert(aos0%particles(iP)%active .eqv. .false.)
         end if
      end do
      !$omp end do
      !$omp end parallel

      call pa_aos0_cpu_release(aos0)

    end block test_block_with_no_inactives

    ! Test

    test_with_inactive_start_of_block: &
    block
      ! add particles uniformly, and then introduce some gaps at the
      ! start of the blocks.
      integer, parameter :: nmax = 4000
      integer            :: nactive = nmax/2
      type (pa_aos0_cpu_t)  :: aos0
      integer            :: iP
      integer            :: tid

      call pa_aos0_cpu_create(aos0, nmax)

      !$omp parallel default(shared) private(iP, tid)
      tid = omp_get_thread_num()
      !$omp do schedule(static, 256)
      do iP = 1, nmax
        associate (particle => aos0%particles(iP))
        particle%active = .false.
        if (ip <= nactive) then
           particle%active = .true.
           particle%ie     = tid
        end if
        end associate
      end do
      !$omp end do
      !$omp end parallel

      call pa_aos0_cpu_update_last(aos0)

      ! Inactivate 16 particles at the start of each block

      !$omp parallel private(iP)
      !$omp do schedule(static, NAME_PA_BLOCKSZ)
      do iP = 1, nmax
        associate(particle => aos0%particles(iP))
          if (mod(iP-1, NAME_PA_BLOCKSZ) < 16) particle%active = .false.
        end associate
      end do
      !$omp end do
      !$omp end parallel

      ! We've removed 16 particles per block (not per thread)
      ! Run the removeInactive()

      nactive = nactive - (1 + nactive/NAME_PA_BLOCKSZ)*16

      call aos0%removeInactive()

      ! check total number of active particles has not changed in the remove

      block
        integer :: nobserved

        nobserved = 0

        !$omp parallel do reduction(+: nobserved)
        do iP = 1, nmax
          if (aos0%particles(iP)%active) nobserved = nobserved + 1
        end do
        !$omp end parallel do

        assert(nobserved == nactive)
        assert(sum(aos0%last(:)) == nactive)
      end block

      ! check the right particles are active

      block
        integer :: iplast

        !$omp parallel private(ip, iplast, tid)
        tid = omp_get_thread_num()
        do iplast = 1, aos0%last(tid)
          iP = pa_aos0_cpu_last_to_index(aos0, tid)
          assert(aos0%particles(iP)%active .eqv. .true.)
          assert(aos0%particles(iP)%ie == tid)
        end do
        !$omp end parallel
      end block

      assert(aos0%lastBlock == nactive/NAME_PA_BLOCKSZ )

      call pa_aos0_cpu_release(aos0)

    end block test_with_inactive_start_of_block

    ! Test
    test_with_inactive_end_of_block: &
    block

      integer, parameter :: nmax = 4000
      integer            :: nactive = nmax
      type (pa_aos0_cpu_t)  :: aos0
      integer            :: iP
      integer            :: tid

      call pa_aos0_cpu_create(aos0, nmax)

      !$omp parallel default(shared) private(iP, tid)
      tid = omp_get_thread_num()
      !$omp do schedule(static, 256)
      do iP = 1, nactive
        associate (particle => aos0%particles(iP))
        particle%active = .false.
        if (mod(iP - 1, NAME_PA_BLOCKSZ) < 16) then
          particle%active = .true.
        end if
        end associate
      end do
      !$omp end do
      !$omp end parallel

      call pa_aos0_cpu_update_last(aos0)
      call aos0%removeInactive()

      ! We should now have (nBlocks/nThreads) x 16 particles for
      ! each thread contiguous at left-hand end of the relevant
      ! blocks

      nactive = (aos0%nBlocks/aos0%nThreads)*16

      !$omp parallel
      !$omp do schedule(static, NAME_PA_BLOCKSZ)
      do iP = 1, aos0%nmaxLocal
        if (pa_aos0_cpu_index_to_last(aos0, iP) <= nactive) then
          assert(aos0%particles(iP)%active .eqv. .true.)
        else
          assert(aos0%particles(iP)%active .eqv. .false.)
        end if
      end do
      !$omp end do
      !$omp end parallel

      call pa_aos0_cpu_release(aos0)

    end block test_with_inactive_end_of_block

  end subroutine pa_aos0_cpu_self_test_removeInactive

end module pa_aos0_cpu
