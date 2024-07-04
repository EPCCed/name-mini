module pa_soa2_gpu

#include "assertion.h"

  use assertion
  use constants
  use maths
  use pa_soa2
  use storage_soa2
  use source
  use util_omp_lib

  implicit none
  private

  ! Particle SOA2 storage management intended for device use

  integer, parameter :: NAME_PA_BLOCKSZ = 128

  type, public, extends(pa_soa2_t) :: pa_soa2_gpu_t
    private
    integer (int32)               :: nmaxLocal        ! particles allocated
    integer (int32)               :: nBlocks          ! size NAME_PA_BLOCKSZ
    integer (int32)               :: lastBlock        ! last active block
  contains
    procedure :: printState        => pa_soa2_gpu_printState
    procedure :: releaseFromSource => pa_soa2_gpu_releaseFromSource
    procedure :: removeInactive    => pa_soa2_gpu_removeInactive
    procedure :: lastParticle      => pa_soa2_gpu_lastParticle
    procedure :: nameScheduleSize  => pa_soa2_gpu_nameScheduleSize
    procedure :: destroy           => pa_soa2_gpu_release
  end type pa_soa2_gpu_t

  public :: pa_soa2_gpu_create_pointer
  public :: pa_soa2_gpu_create
  public :: pa_soa2_gpu_release

  public :: pa_soa2_gpu_releaseFromSource
  public :: pa_soa2_gpu_printState
  public :: pa_soa2_gpu_lastParticle
  public :: pa_soa2_gpu_nameScheduleSize

contains

  !---------------------------------------------------------------------------

  function pa_soa2_gpu_create_pointer(maxLocal) result(pa)

    ! Return a pointer to a new object; assume this is the top-level
    ! object and the device mapping should always take place.

    integer,               intent(in) :: maxLocal
    class (pa_soa2_gpu_t), pointer    :: pa

    allocate(pa)
    call pa_soa2_gpu_create(pa, maxLocal, map = .true.)

  end function pa_soa2_gpu_create_pointer

  !---------------------------------------------------------------------------

  subroutine pa_soa2_gpu_create(soa2, maxLocal, map)

    type (pa_soa2_gpu_t), intent(out) :: soa2
    integer (int32),      intent(in)  :: maxLocal
    logical, optional,    intent(in)  :: map

    logical :: maptarget

    maptarget = .false.
    if (present(map)) maptarget = map

    ! Storage: allocate a whole number of blocks
    ! Furthermore, allocate a power-of-two number of blocks

    soa2%nBlocks   = mathsNearestPowerOfTwo(1 + maxLocal/NAME_PA_BLOCKSZ)
    soa2%nmaxLocal = soa2%nBlocks*NAME_PA_BLOCKSZ
    call storage_soa2_create(soa2%particles, soa2%nmaxLocal, map = .false.)

    ! Management
    soa2%lastBlock = 0

    if (maptarget) then
      !$omp target enter data map(always, to: soa2)

      !$omp target enter data map(alloc: soa2%particles)

      !$omp target enter data map(alloc: soa2%particles%iUP)
      !$omp target enter data map(alloc: soa2%particles%iE)

      !$omp target enter data map(alloc: soa2%particles%x)
      !$omp target enter data map(alloc: soa2%particles%xOld)
      !$omp target enter data map(alloc: soa2%particles%t)
      !$omp target enter data map(alloc: soa2%particles%tOld)
      !$omp target enter data map(alloc: soa2%particles%t0)
      !$omp target enter data map(alloc: soa2%particles%needToSetVel)
      !$omp target enter data map(alloc: soa2%particles%active)
      !$omp target enter data map(alloc: soa2%particles%marked)

      !$omp target enter data map(alloc: soa2%particles%iHCoord)
      !$omp target enter data map(alloc: soa2%particles%iZCoord)
      !$omp target enter data map(alloc: soa2%particles%iHCoordOld)
      !$omp target enter data map(alloc: soa2%particles%iZCoordOld)
      !$omp target enter data map(alloc: soa2%particles%iSource)

      !$omp target enter data map(alloc: soa2%particles%h)
      !$omp target enter data map(alloc: soa2%particles%wSed)
      !$omp target enter data map(alloc: soa2%particles%diameter)
      !$omp target enter data map(alloc: soa2%particles%density)
      !$omp target enter data map(alloc: soa2%particles%particleShape)
      !$omp target enter data map(alloc: soa2%particles%elongation)
      !$omp target enter data map(alloc: soa2%particles%flatness)
      !$omp target enter data map(alloc: soa2%particles%shapeSchemeCode)

      !$omp target enter data map(alloc: soa2%particles%s)
    end if

    ! Make sure all particles are inactive
    block
      integer :: iP

      !$omp target teams
      !$omp distribute parallel do
      do iP = 1, soa2%nmaxLocal
        soa2%particles%active(iP) = .false.
      end do
      !$omp end distribute parallel do
      !$omp end target teams
    end block

  end subroutine pa_soa2_gpu_create

  !---------------------------------------------------------------------------

  subroutine pa_soa2_gpu_release(self)

    class (pa_soa2_gpu_t), intent(inout) :: self

    logical :: unmaptarget

    unmaptarget = .true.

    if (unmaptarget) then

      !$omp target exit data map(delete: self%particles%s)

      !$omp target exit data map(delete: self%particles%shapeSchemeCode)
      !$omp target exit data map(delete: self%particles%flatness)
      !$omp target exit data map(delete: self%particles%elongation)
      !$omp target exit data map(delete: self%particles%particleShape)
      !$omp target exit data map(delete: self%particles%density)
      !$omp target exit data map(delete: self%particles%diameter)
      !$omp target exit data map(delete: self%particles%wSed)

      !$omp target exit data map(delete: self%particles%h)
       
      !$omp target exit data map(delete: self%particles%iSource)
      !$omp target exit data map(delete: self%particles%iZCoordOld)
      !$omp target exit data map(delete: self%particles%iHCoordOld)
      !$omp target exit data map(delete: self%particles%iZCoord)
      !$omp target exit data map(delete: self%particles%iHCoord)

      !$omp target exit data map(delete: self%particles%marked)
      !$omp target exit data map(delete: self%particles%active)
      !$omp target exit data map(delete: self%particles%needToSetVel)
      !$omp target exit data map(delete: self%particles%t0)
      !$omp target exit data map(delete: self%particles%tOld)
      !$omp target exit data map(delete: self%particles%t)
      !$omp target exit data map(delete: self%particles%xOld)
      !$omp target exit data map(delete: self%particles%x)

      !$omp target exit data map(delete: self%particles%iE)
      !$omp target exit data map(delete: self%particles%iUP)

      !$omp target exit data map(delete: self%particles)
      !$omp target exit data map(delete: self)
    end if

    call storage_soa2_release(self%particles, unmap = .false.)

    self%lastBlock = 0
    self%nBlocks   = 0
    self%nmaxLocal = 0

  end subroutine pa_soa2_gpu_release

  !---------------------------------------------------------------------------

  subroutine pa_soa2_gpu_releaseFromSource(self, nRelease, source)

    ! Release nRelease particles from source

    class (pa_soa2_gpu_t), intent(inout) :: self
    integer,               intent(in)    :: nRelease
    type (source_t),       intent(in)    :: source

    integer (int32) :: iR
    integer (int64) :: iUP
    integer (int32) :: iP

    ! Check we have enough capacity to store nRelease new particles

    if (self%lastParticle() + nRelease > self%nmaxLocal) then
      print *, "Attempt to release too many particles"
      print *, "nRelease:                            ", nRelease
      print *, "lastParticle()                       ", self%lastParticle()
      print *, "nmaxLocal                            ", self%nmaxLocal
      stop     "Fatal"
    end if

    ! Contiguous release from first free position

    !$omp target teams
    !$omp distribute parallel do private(iR, iUP, iP)
    do iR = 1, nRelease

      iUP = source%iUP0 + iR                 ! Revisit as iUP0 needs updating
      iP = self%lastBlock*NAME_PA_BLOCKSZ + iR
      call storage_soa2_initialiseFromSource(self%particles, iP, iUP, source)

    end do
    !$omp end distribute parallel do
    !$omp end target teams

    block
      integer :: nblock
      ! Always move forward a whole number of blocks
      nblock = (NAME_PA_BLOCKSZ + nRelease - 1)/NAME_PA_BLOCKSZ
      self%lastBlock = self%lastBlock + nblock
      !$omp target update to(self%lastBlock)
    end block

  end subroutine pa_soa2_gpu_releaseFromSource

  !---------------------------------------------------------------------------

  subroutine pa_soa2_gpu_removeInactive(self)

    class (pa_soa2_gpu_t), intent(inout) :: self

    ! Kernel ...

  end subroutine pa_soa2_gpu_removeInactive

  !---------------------------------------------------------------------------

  subroutine pa_soa2_gpu_printState(self)

    ! Produce a report on the current state of the management

    class (pa_soa2_gpu_t), intent(in) :: self

    integer :: nactive

    nactive = self%nActive()

    print *, ""
    print *, "Manager: pa_soa2_gpu            "

    print *, "Manager: NAME_PA_BLOCSZ         ", NAME_PA_BLOCKSZ
    print *, "Manager: maxLocal               ", self%nmaxLocal
    print *, "Manager: nBlocks                ", self%nBlocks
    print *, "Manager: omp_get_num_devices()  ", omp_get_num_devices()

    print *, "Manager: nactive                ", nactive
    print *, "Manager: lastBlock              ", self%lastBlock
    print *, "Manager: lastActiveParticle     ", self%lastParticle()

  end subroutine pa_soa2_gpu_printState

  !---------------------------------------------------------------------------

  pure function pa_soa2_gpu_lastParticle(self) result(last)

    ! Whole number of blocks

    class (pa_soa2_gpu_t), intent(in) :: self
    integer                        :: last

    last = self%lastBlock*NAME_PA_BLOCKSZ

  end function pa_soa2_gpu_lastParticle

  !---------------------------------------------------------------------------

  function pa_soa2_gpu_nameScheduleSize(self) result(nScheduleSize)

    ! The static schedule size is 1 for GPU

    class (pa_soa2_gpu_t), intent(in) :: self
    integer                        :: nScheduleSize

    nScheduleSize = 1

  end function pa_soa2_gpu_nameScheduleSize

end module pa_soa2_gpu
