module pa_aos0

  ! Extends pa_t to make concrete the storage
  ! This allows the update to be defined; management is still deferred.

  use constants
  use pa
  use storage_aos0
  use flow
  use maths

  implicit none
  private

  type, extends(pa_t), abstract, public :: pa_aos0_t
    type (aos0_t), allocatable :: particles(:)
    real (std),    allocatable :: masses(:)
  contains
    ! Defer
    procedure (if_x_releaseFromSource), pass, deferred :: releaseFromSource
    procedure (if_x_removeInactive),    pass, deferred :: removeInactive
    procedure (if_x_lastParticle),      pass, deferred :: lastParticle
    procedure (if_x_printState),        pass, deferred :: printState
    procedure (if_x_nameScheduleSize),  pass, deferred :: nameScheduleSize
    ! Define
    procedure :: nActive    => pa_aos0_nactive
    procedure :: updateStep => pa_aos0_update
    ! Defer
    procedure (if_x_destroy),           pass, deferred :: destroy
  end type pa_aos0_t


  abstract interface

    subroutine if_x_releaseFromSource(self, nRelease, source)
      import pa_aos0_t
      import source_t
      class (pa_aos0_t), intent(inout) :: self
      integer,           intent(in)    :: nRelease
      type (source_t),   intent(in)    :: source
    end subroutine if_x_releaseFromSource

    subroutine if_x_removeInactive(self)
      import pa_aos0_t
      class (pa_aos0_t), intent(inout) :: self
    end subroutine if_x_removeInactive

    pure function if_x_lastParticle(self) result(last)
      import pa_aos0_t
      class (pa_aos0_t), intent(in) :: self
      integer                       :: last
    end function if_x_lastParticle

    subroutine if_x_printState(self)
      import pa_aos0_t
      class (pa_aos0_t), intent(in) :: self
    end subroutine if_x_printState

    function if_x_nameScheduleSize(self) result(nScheduleSize)
      import pa_aos0_t
      class (pa_aos0_t), intent(in) :: self
      integer                       :: nScheduleSize
    end function if_x_nameScheduleSize

    subroutine if_x_destroy(self)
      import pa_aos0_t
      class (pa_aos0_t), intent(inout)  :: self
    end subroutine if_x_destroy

  end interface

  public :: pa_aos0_nactive
  public :: pa_aos0_update

contains

  !---------------------------------------------------------------------------

  function pa_aos0_nactive(self) result(nActiveParticles)

    ! Count the current number of active particles

    ! See note in pa_soa1_t analogue

    class (pa_aos0_t), intent(in) :: self
    integer                       :: nActiveParticles

    integer :: nactive
    integer :: nmaxLocal
    integer :: iP

    nactive   = 0
    nmaxLocal = size(self%particles(:))

    !$omp target enter data map(always, to: nmaxLocal, nactive)

    !$omp target teams reduction(+:nactive)
    !$omp distribute parallel do private(iP) reduction(+: nactive)
    do iP = 1, nmaxLocal
      if (self%particles(iP)%active) nactive = nactive + 1
    end do
    !$omp end distribute parallel do
    !$omp end target teams

    !$omp target exit data map(from: nactive)
    !$omp target exit data map(release: nmaxLocal)

    nActiveParticles = nactive

  end function pa_aos0_nactive

  !---------------------------------------------------------------------------

  subroutine pa_aos0_update(self, flow, dt)

    ! Update driver: update each particle one step of length dt seconds
    ! with background flow

    class (pa_aos0_t), intent(inout) :: self
    type (flow_t),     intent(in)    :: flow
    real (std),        intent(in)    :: dt

    integer :: lastActiveParticle
    integer :: nSchedule
    integer :: iP

    ! Update particle position for [t, t + dt]
    ! The flow should already be present on the target

    lastActiveParticle = self%lastParticle()
    nSchedule          = self%nameScheduleSize()

    !$omp target enter data map(always, to: lastActiveParticle, nSchedule)

    !$omp target teams
    !$omp distribute parallel do private(iP) schedule(static, nschedule)
    do iP = 1, lastActiveParticle
      if (self%particles(iP)%active .eqv. .false.) cycle
      call pa_aos0_update_particle(self%particles, iP, flow, dt)
    end do
    !$omp end distribute parallel do
    !$omp end target teams

    !$omp target exit data map(release: lastActiveParticle, nSchedule)

  end subroutine pa_aos0_update

  !----------------------------------------------------------------------------
  !----------------------------------------------------------------------------

  subroutine pa_aos0_update_particle(particles, iP, flow, dt)

    !$omp declare target

    ! Update for particles(iP)

    type (aos0_t), intent(inout) :: particles(:)
    integer,       intent(in)    :: iP
    type (flow_t), intent(in)    :: flow
    real (std),    intent(in)    :: dt

    real (std)                   :: rv(3)
    integer                      :: k

    call mathsRandomGaussian3(particles(iP)%s, rv)

    ! need to check the appropriate flow level for the current particle

    k = flowZToLevel(flow, particles(iP)%x(3))

    block
      real (std) :: dx(3)
      dx(1) = (flow%u(1,k)*dt + 2.0*sqrt(flow%k(1)*dt)*rv(1))
      dx(2) = (flow%u(2,k)*dt + 2.0*sqrt(flow%k(2)*dt)*rv(2))
      dx(3) = (flow%u(3,k)*dt + 2.0*sqrt(flow%k(3)*dt)*rv(3))
      particles(iP)%x(1) = particles(iP)%x(1) + dx(1)
      particles(iP)%x(2) = particles(iP)%x(2) + dx(2)
      particles(iP)%x(3) = particles(iP)%x(3) + dx(3)
    end block

    ! Check lower and upper boundaries ...
    if (particles(iP)%x(3) < flow%zb(1)) particles(iP)%active = .false.
    if (particles(iP)%x(3) > flow%zb(flow%nZ+1)) particles(iP)%active = .false.

    return

  end subroutine pa_aos0_update_particle

end module pa_aos0
