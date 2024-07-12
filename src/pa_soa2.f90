module pa_soa2

  ! Extends pa_t to make concrete storage

  use constants
  use pa
  use storage_soa2
  use flow
  use maths

  implicit none
  private

  type, extends(pa_t), abstract, public :: pa_soa2_t
    type (soa2_t)           :: particles
    real (std), allocatable :: masses(:,:)
  contains
    ! Defer
    procedure (if_x_releaseFromSource), pass, deferred :: releaseFromSource
    procedure (if_x_removeInactive),    pass, deferred :: removeInactive
    procedure (if_x_lastParticle),      pass, deferred :: lastParticle
    procedure (if_x_printState),        pass, deferred :: printState
    procedure (if_x_nameScheduleSize),  pass, deferred :: nameScheduleSize
    ! Define
    procedure :: nActive    => pa_soa2_nactive
    procedure :: updateStep => pa_soa2_update
    ! Defer
    procedure (if_x_destroy),           pass, deferred :: destroy
  end type pa_soa2_t

  abstract interface

    subroutine if_x_releaseFromSource(self, nRelease, source)
      import pa_soa2_t
      import source_t
      class (pa_soa2_t), intent(inout) :: self
      integer,           intent(in)    :: nRelease
      type (source_t),   intent(in)    :: source
    end subroutine if_x_releaseFromSource

    subroutine if_x_removeInactive(self)
      import pa_soa2_t
      class (pa_soa2_t), intent(inout) :: self
    end subroutine if_x_removeInactive

    pure function if_x_lastParticle(self) result(last)
      import pa_soa2_t
      class (pa_soa2_t), intent(in) :: self
      integer                       :: last
    end function if_x_lastParticle

    subroutine if_x_printState(self)
      import pa_soa2_t
      class (pa_soa2_t), intent(in) :: self
    end subroutine if_x_printState

    function if_x_nameScheduleSize(self) result(nScheduleSize)
      import pa_soa2_t
      class (pa_soa2_t), intent(in) :: self
      integer                       :: nScheduleSize
    end function if_x_nameScheduleSize

    subroutine if_x_destroy(self)
      import pa_soa2_t
      class (pa_soa2_t), intent(inout)  :: self
    end subroutine if_x_destroy

  end interface

  public :: pa_soa2_nactive
  public :: pa_soa2_update

contains

  !---------------------------------------------------------------------------

  function pa_soa2_nactive(self) result(nActiveParticles)

    ! Count the current number of active particles

    ! See note in pa_soa1_t analogue

    class (pa_soa2_t), intent(in) :: self
    integer                       :: nActiveParticles

    integer :: nactive
    integer :: nmaxLocal
    integer :: iP

    nactive   = 0
    nmaxLocal = size(self%particles%active(:))

    !$omp target enter data map(always, to: nmaxLocal, nactive)

    !$omp target teams reduction(+: nactive)
    !$omp distribute parallel do private(iP) reduction(+: nactive)
    do iP = 1, nmaxLocal
      if (self%particles%active(iP)) nactive = nactive + 1
    end do
    !$omp end distribute parallel do
    !$omp end target teams

    !$omp target exit data map(from: nactive)
    !$omp target exit data map(release: nmaxLocal)

    nActiveParticles = nactive

  end function pa_soa2_nactive

  !---------------------------------------------------------------------------

  subroutine pa_soa2_update(self, flow, dt)

    ! Update driver

    class (pa_soa2_t), intent(inout) :: self
    type  (flow_t),    intent(in)    :: flow
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
      if (self%particles%active(iP) .eqv. .false.) cycle
      call pa_soa2_update_particle(self%particles, iP, flow, dt)
    end do
    !$omp end distribute parallel do
    !$omp end target teams

    !$omp target exit data map(release: lastActiveParticle, nSchedule)

  end subroutine pa_soa2_update

  !----------------------------------------------------------------------------

  subroutine pa_soa2_update_particle(particles, iP, flow, dt)

    !$omp declare target

    type (soa2_t), intent(inout) :: particles
    integer,       intent(in)    :: iP
    type (flow_t), intent(in)    :: flow
    real (std),    intent(in)    :: dt

    real (std) :: rv(3)
    integer    :: k

    if (.false.) then
      ! Don't do this. Performance disaster awaits ...
      call mathsRandomGaussian3(particles%s(iP, :), rv)
    else
      block
        ! We must form a contiguous argument ...
        ! The 1:4 must be unrolled explicitly
        integer (int32) :: s(4)
        s(1) = particles%s(iP, 1)
        s(2) = particles%s(iP, 2)
        s(3) = particles%s(iP, 3)
        s(4) = particles%s(iP, 4)
        call mathsRandomGaussian3(s, rv)
        particles%s(iP, 1) = s(1)
        particles%s(iP, 2) = s(2)
        particles%s(iP, 3) = s(3)
        particles%s(iP, 4) = s(4)
      end block
   end if

   ! Determine the current level
   k = flowZToLevel(flow, particles%x(iP, 3))

   particles%x(iP, 1) = particles%x(iP, 1) &
                      + (flow%u(1,k)*dt + 2.0*sqrt(flow%k(1)*dt)*rv(1))
   particles%x(iP, 2) = particles%x(iP, 2) &
                      + (flow%u(2,k)*dt + 2.0*sqrt(flow%k(2)*dt)*rv(2))
   particles%x(iP, 3) = particles%x(iP, 3) &
                      + (flow%u(3,k)*dt + 2.0*sqrt(flow%k(3)*dt)*rv(3))

    ! Check lower and upper boundaries ...
    if (particles%x(iP, 3) < flow%zb(1)) particles%active(iP) = .false.
    if (particles%x(iP, 3) > flow%zb(flow%nZ+1)) particles%active(iP) = .false.

  end subroutine pa_soa2_update_particle

end module pa_soa2
