module pa_soa1

  ! Extends pa_t to make concrete storage
  ! This alos allows the update to be defined

  use constants
  use pa
  use storage_soa1
  use flow
  use maths

  implicit none
  private

  type, extends(pa_t), abstract, public :: pa_soa1_t
    type (soa1_t)           :: particles
    real (std), allocatable :: masses(:,:)
  contains
    ! Defer
    procedure (if_x_releaseFromSource), pass, deferred :: releaseFromSource
    procedure (if_x_removeInactive),    pass, deferred :: removeInactive
    procedure (if_x_lastParticle),      pass, deferred :: lastParticle
    procedure (if_x_printState),        pass, deferred :: printState
    procedure (if_x_nameScheduleSize),  pass, deferred :: nameScheduleSize
    ! Define
    procedure :: nActive    => pa_soa1_nactive
    procedure :: updateStep => pa_soa1_update
    ! Defer
    procedure (if_x_destroy),           pass, deferred ::  destroy
  end type pa_soa1_t

  abstract interface

    subroutine if_x_releaseFromSource(self, nRelease, source)
      import pa_soa1_t
      import source_t
      class (pa_soa1_t), intent(inout) :: self
      integer,           intent(in)    :: nRelease
      type (source_t),   intent(in)    :: source
    end subroutine if_x_releaseFromSource

    subroutine if_x_removeInactive(self)
      import pa_soa1_t
      class (pa_soa1_t), intent(inout) :: self
    end subroutine if_x_removeInactive

    pure function if_x_lastParticle(self) result(last)
      import pa_soa1_t
      class (pa_soa1_t), intent(in) :: self
      integer                  :: last
    end function if_x_lastParticle

    subroutine if_x_printState(self)
      import pa_soa1_t
      class (pa_soa1_t), intent(in) :: self
    end subroutine if_x_printState

    function if_x_nameScheduleSize(self) result(nScheduleSize)
      import pa_soa1_t
      class (pa_soa1_t), intent(in) :: self
      integer                  :: nScheduleSize
    end function if_x_nameScheduleSize

    subroutine if_x_destroy(self)
      import pa_soa1_t
      class (pa_soa1_t), intent(inout)  :: self
    end subroutine if_x_destroy

  end interface

  public :: pa_soa1_nactive
  public :: pa_soa1_update

contains

  !---------------------------------------------------------------------------

  function pa_soa1_nactive(self) result(nActiveParticles)

    ! Count the current number of active particles

    ! Note at nvidia SDK 24.5
    ! Don't try to perform the reduction on the result variable as
    ! this seems to provide the wrong answer in the device context.
    ! Use a local variable.
    
    class (pa_soa1_t), intent(in) :: self
    integer                       :: nActiveParticles

    integer :: nactive
    integer :: nmaxLocal
    integer :: iP

    nActive   = 0
    nmaxLocal        = size(self%particles%active(:))

    !$omp target enter data map(always, to: nmaxLocal, nActive)

    !$omp target teams reduction(+: nActive)
    !$omp distribute parallel do private(iP) reduction(+: nActive)
    do iP = 1, nmaxLocal
      if (self%particles%active(iP)) nActive = nActive + 1
    end do
    !$omp end distribute parallel do
    !$omp end target teams

    !$omp target exit data map(from: nActive)
    !$omp target exit data map(release: nmaxLocal)

    nActiveParticles = nactive

  end function pa_soa1_nactive

  !---------------------------------------------------------------------------

  subroutine pa_soa1_update(self, flow, dt)

    ! Update driver: update each particle one self of length dt seconds
    ! with background flow

    class (pa_soa1_t), intent(inout) :: self
    type  (flow_t),    intent(in)    :: flow
    real  (std),       intent(in)    :: dt

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
      if (self%particles%active(iP) .eqv. .false.) cycle
      call pa_soa1_update_particle(self%particles, iP, flow, dt)
    end do
    !$omp end distribute parallel do
    !$omp end target teams

    !$omp target exit data map(release: nSchedule, lastActiveParticle)

  end subroutine pa_soa1_update

  !---------------------------------------------------------------------------

  subroutine pa_soa1_update_particle(particles, iP, flow, dt)

    !$omp declare target

    type (soa1_t), intent(inout), target :: particles
    integer,       intent(in)            :: iP
    type (flow_t), intent(in)            :: flow
    real (std),    intent(in)            :: dt

    integer                         :: k
    real (std)                      :: rv(3)
    real (std), pointer, contiguous :: x(:)

    x => particles%x(1:3, iP)

    ! Best effort c. 40 registers per thread

    call mathsRandomGaussian3(particles%s(1:4, iP), rv)

    ! Determine the current level

    k = flowZToLevel(flow, x(3))

    x(1) = x(1) + (flow%u(1,k)*dt + 2.0*sqrt(flow%k(1)*dt)*rv(1))
    x(2) = x(2) + (flow%u(2,k)*dt + 2.0*sqrt(flow%k(2)*dt)*rv(2))
    x(3) = x(3) + (flow%u(3,k)*dt + 2.0*sqrt(flow%k(3)*dt)*rv(3))

    ! Check lower and upper boundaries ...
    if (x(3) < flow%zb(1)) particles%active(iP) = .false.
    if (x(3) > flow%zb(flow%nZ+1)) particles%active(iP) = .false.

  end subroutine pa_soa1_update_particle

end module pa_soa1
