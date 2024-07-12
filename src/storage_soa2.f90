module storage_soa2

  ! SOA particle storage characterised by storage order
  ! favouring complete coalescing on the GPU. E.g., for
  ! for the random number state:
  !
  !   integer (int32) :: s(nParticles, 4)
  !
  ! All the allocatable compoenets are actually common
  ! between this and the soa1 implmemtation, but it is
  ! repeated here.

#include "assertion.h"

  use assertion
  use constants
  use maths
  use rtime
  use source

  implicit none
  private

  type, public :: soa2_t
    integer                      :: nParticles

    integer (int64), allocatable :: iUP(:)
    integer (int32), allocatable :: iE(:)
    real (std),      allocatable :: x(:, :)
    real (std),      allocatable :: xOld(:, :)

    type (sTime_t),  allocatable :: t(:)
    type (sTime_t),  allocatable :: tOld(:)
    type (sTime_t),  allocatable :: t0(:)
    logical,         allocatable :: needToSetVel(:)
    logical,         allocatable :: active(:)
    integer,         allocatable :: marked(:)

    integer,         allocatable :: iHCoord(:)
    integer,         allocatable :: iZCoord(:)
    integer,         allocatable :: iHCoordOld(:)
    integer,         allocatable :: iZCoordOld(:)
    integer,         allocatable :: iSource(:)

    real (std),      allocatable :: h(:)

    real (std),      allocatable :: wSed(:)
    real (std),      allocatable :: diameter(:)
    real (std),      allocatable :: density(:)
    real (std),      allocatable :: particleShape(:)
    real (std),      allocatable :: elongation(:)
    real (std),      allocatable :: flatness(:)
    integer (int32), allocatable :: shapeSchemeCode(:)

    integer (int32), allocatable :: s(:, :)
  end type soa2_t

  public :: storage_soa2_create
  public :: storage_soa2_release
  public :: storage_soa2_initialiseFromSource
  public :: storage_soa2_copy

contains

  !---------------------------------------------------------------------------

  subroutine storage_soa2_create(pa, nParticles, map)

    type (soa2_t),     intent(out) :: pa
    integer,           intent(in)  :: nParticles
    logical, optional, intent(in)  :: map

    logical :: maptarget

    maptarget = .false.
    if (present(map)) maptarget = map

    allocate(pa%iUP(nParticles))
    allocate(pa%iE(nParticles))

    allocate(pa%x(nParticles, 3))
    allocate(pa%xOld(nParticles, 3))

    allocate(pa%t(nParticles))
    allocate(pa%tOld(nParticles))
    allocate(pa%t0(nParticles))
    allocate(pa%needToSetVel(nParticles))
    allocate(pa%active(nParticles))
    allocate(pa%marked(nParticles))

    allocate(pa%iHCoord(nParticles))
    allocate(pa%iZCoord(nParticles))
    allocate(pa%iHCoordOld(nParticles))
    allocate(pa%iZCoordOld(nParticles))
    allocate(pa%iSource(nParticles))

    allocate(pa%h(nParticles))
    allocate(pa%wSed(nParticles))
    allocate(pa%diameter(nParticles))
    allocate(pa%density(nParticles))
    allocate(pa%particleShape(nParticles))
    allocate(pa%elongation(nParticles))
    allocate(pa%flatness(nParticles))
    allocate(pa%shapeSchemeCode(nParticles))

    allocate(pa%s(nParticles, 4))

    pa%nParticles = nParticles

    if (maptarget) then
      !$omp target enter data map(to:    pa)
      !$omp target enter data map(alloc: pa%iUP)
      !$omp target enter data map(alloc: pa%iE)

      !$omp target enter data map(alloc: pa%x)
      !$omp target enter data map(alloc: pa%xOld)
      !$omp target enter data map(alloc: pa%t)
      !$omp target enter data map(alloc: pa%tOld)
      !$omp target enter data map(alloc: pa%t0)
      !$omp target enter data map(alloc: pa%needToSetVel)
      !$omp target enter data map(alloc: pa%active)
      !$omp target enter data map(alloc: pa%marked)

      !$omp target enter data map(alloc: pa%iHCoord)
      !$omp target enter data map(alloc: pa%iZCoord)
      !$omp target enter data map(alloc: pa%iHCoordOld)
      !$omp target enter data map(alloc: pa%iZCoordOld)
      !$omp target enter data map(alloc: pa%iSource)

      !$omp target enter data map(alloc: pa%h)
      !$omp target enter data map(alloc: pa%wSed)
      !$omp target enter data map(alloc: pa%diameter)
      !$omp target enter data map(alloc: pa%density)
      !$omp target enter data map(alloc: pa%particleShape)
      !$omp target enter data map(alloc: pa%elongation)
      !$omp target enter data map(alloc: pa%flatness)
      !$omp target enter data map(alloc: pa%shapeSchemeCode)

      !$omp target enter data map(alloc: pa%s)
    end if

  end subroutine storage_soa2_create

  !---------------------------------------------------------------------------

  subroutine storage_soa2_release(pa, unmap)

    type (soa2_t),     intent(inout) :: pa
    logical, optional, intent(in)    :: unmap

    logical :: unmaptarget

    unmaptarget = .false.
    if (present(unmap)) unmaptarget = unmap

    if (unmaptarget) then
      !$omp target exit data map(delete: pa%s)

      !$omp target exit data map(delete: pa%shapeSchemeCode)
      !$omp target exit data map(delete: pa%flatness)
      !$omp target exit data map(delete: pa%elongation)
      !$omp target exit data map(delete: pa%particleShape)
      !$omp target exit data map(delete: pa%density)
      !$omp target exit data map(delete: pa%diameter)
      !$omp target exit data map(delete: pa%wSed)

      !$omp target exit data map(delete: pa%h)
       
      !$omp target exit data map(delete: pa%iSource)
      !$omp target exit data map(delete: pa%iZCoordOld)
      !$omp target exit data map(delete: pa%iHCoordOld)
      !$omp target exit data map(delete: pa%iZCoord)
      !$omp target exit data map(delete: pa%iHCoord)

      !$omp target exit data map(delete: pa%marked)
      !$omp target exit data map(delete: pa%active)
      !$omp target exit data map(delete: pa%needToSetVel)
      !$omp target exit data map(delete: pa%t0)
      !$omp target exit data map(delete: pa%tOld)
      !$omp target exit data map(delete: pa%t)
      !$omp target exit data map(delete: pa%xOld)
      !$omp target exit data map(delete: pa%x)

      !$omp target exit data map(delete: pa%iE)
      !$omp target exit data map(delete: pa%iUP)
      !$omp target exit data map(delete: pa)
    end if

    deallocate(pa%s)

    deallocate(pa%shapeSchemeCode)
    deallocate(pa%flatness)
    deallocate(pa%elongation)
    deallocate(pa%particleShape)
    deallocate(pa%density)
    deallocate(pa%diameter)
    deallocate(pa%wSed)

    deallocate(pa%h)
       
    deallocate(pa%iSource)
    deallocate(pa%iZCoordOld)
    deallocate(pa%iHCoordOld)
    deallocate(pa%iZCoord)
    deallocate(pa%iHCoord)

    deallocate(pa%marked)
    deallocate(pa%active)
    deallocate(pa%needToSetVel)
    deallocate(pa%t0)
    deallocate(pa%tOld)
    deallocate(pa%t)
    deallocate(pa%xOld)
    deallocate(pa%x)

    deallocate(pa%iE)
    deallocate(pa%iUP)

  end subroutine storage_soa2_release

  !---------------------------------------------------------------------------

  subroutine storage_soa2_initialiseFromSource(pa, iP, iUP, source)

    type (soa2_t),   intent(inout) :: pa
    integer (int32), intent(in)    :: iP
    integer (int64), intent(in)    :: iUP
    type (source_t), intent(in)    :: source

    !$omp declare target

    assert(0 < iP .and. iP <= pa%nParticles)

    pa%iUP(iP)           = iUP
    pa%x(iP, 1)          = source%x(1)
    pa%x(iP, 2)          = source%x(2)
    pa%x(iP, 3)          = source%x(3)
    pa%xold(iP, 1)       = source%x(1)
    pa%xold(iP, 2)       = source%x(2)
    pa%xold(iP, 3)       = source%x(3)

    pa%t(iP)             = sTimeZero()
    pa%tOld(iP)          = pa%t(iP)
    pa%t0(iP)            = pa%t(iP)

    pa%needToSetVel(iP)  = .false.
    pa%active(iP)        = .true.
    pa%marked(iP)        = 0
    pa%iHCoord(iP)       = 0
    pa%iZCoord(iP)       = 0
    pa%iHCoordOld(iP)    = 0
    pa%iZCoordOld(iP)    = 0
    pa%iSource(iP)       = 0

    pa%h(iP)             = 0.0
    pa%wSed(iP)          = 0.0
    pa%diameter(iP)      = 0.0
    pa%density(iP)       = 0.0
    pa%particleShape(iP) = 0.0
    pa%elongation(iP)    = 0.0
    pa%flatness(iP)      = 0.0
    pa%shapeSchemeCode(iP) = 0

    block
      integer (int32) :: s0(4)
      s0(1:4)     = mathsRNGStateAdvanceFrom(source%seed0, iUP)
      pa%s(iP, 1) = s0(1)
      pa%s(iP, 2) = s0(2)
      pa%s(iP, 3) = s0(3)
      pa%s(iP, 4) = s0(4)
    end block

  end subroutine storage_soa2_initialiseFromSource

  !----------------------------------------------------------------------------

  subroutine storage_soa2_copy(pa, idst, isrc)

    ! Copy particle data from isrc to idst.

    type (soa2_t),   intent(inout) :: pa
    integer (int32), intent(in)    :: idst
    integer (int32), intent(in)    :: isrc

    !$omp declare target

    assert(0 < idst .and. idst <= pa%nParticles)
    assert(0 < isrc .and. isrc <= pa%nParticles)

    pa%iUP(idst)             = pa%iUP(isrc)
    pa%x(idst, :)            = pa%x(isrc, :)
    pa%xold(idst, :)         = pa%xOld(isrc, :)

    pa%t(idst)               = pa%t(isrc)
    pa%tOld(idst)            = pa%tOld(isrc)
    pa%t0(idst)              = pa%t0(isrc)

    pa%needToSetVel(idst)    = pa%needToSetVel(isrc)
    pa%active(idst)          = pa%active(isrc)
    pa%marked(idst)          = pa%marked(isrc)
    pa%iHCoord(idst)         = pa%iHCoord(isrc)
    pa%iZCoord(idst)         = pa%iZCoord(isrc)
    pa%iHCoordOld(idst)      = pa%iHCoordOld(isrc)
    pa%iZCoordOld(idst)      = pa%iZCoordOld(isrc)
    pa%iSource(idst)         = pa%iSource(isrc)

    pa%h(idst)               = pa%h(isrc)
    pa%wSed(idst)            = pa%wSed(isrc)
    pa%diameter(idst)        = pa%diameter(isrc)
    pa%density(idst)         = pa%density(isrc)
    pa%particleShape(idst)   = pa%particleShape(isrc)
    pa%elongation(idst)      = pa%elongation(isrc)
    pa%flatness(idst)        = pa%flatness(isrc)
    pa%shapeSchemeCode(idst) = pa%shapeSchemeCode(isrc)

    pa%s(idst, 1:4)          = pa%s(isrc, 1:4)

  end subroutine storage_soa2_copy

end module storage_soa2
