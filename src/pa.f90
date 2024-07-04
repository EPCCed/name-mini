module pa

  ! Abstract particle array
  !
  ! Actual storage is deferred to one of the extending classes pa_*.f90
  ! Sub-module pa_submodule.f90 provides a factory method for an object
  ! of the abstract class.

  use constants
  use flow
  use options
  use source

  implicit none
  public

  type, abstract, public :: pa_t
  contains
    procedure (if_releaseFromSource), pass, deferred :: releaseFromSource
    procedure (if_removeInactive),    pass, deferred :: removeInactive
    procedure (if_lastParticle),      pass, deferred :: lastParticle
    procedure (if_printState),        pass, deferred :: printState
    procedure (if_nameScheduleSize),  pass, deferred :: nameScheduleSize
    procedure (if_nActive),           pass, deferred :: nActive
    procedure (if_updateStep),        pass, deferred :: updateStep
    procedure (if_destroy),           pass, deferred :: destroy
  end type pa_t

  abstract interface

    subroutine if_releaseFromSource(self, nRelease, source)
      import pa_t
      import source_t
      class (pa_t),    intent(inout) :: self
      integer,         intent(in)    :: nRelease
      type (source_t), intent(in)    :: source
    end subroutine if_releaseFromSource

    subroutine if_removeInactive(self)
      import pa_t
      class (pa_t), intent(inout) :: self
    end subroutine if_removeInactive

    pure function if_lastParticle(self) result(last)
      import pa_t
      class (pa_t), intent(in) :: self
      integer                  :: last
    end function if_lastParticle

    subroutine if_printState(self)
      import pa_t
      class (pa_t), intent(in) :: self
    end subroutine if_printState

    function if_nameScheduleSize(self) result(nScheduleSize)
      import pa_t
      class (pa_t), intent(in) :: self
      integer                  :: nScheduleSize
    end function if_nameScheduleSize

    function if_nActive(self) result(nActiveParticles)
      import pa_t
      class (pa_t), intent(in) :: self
      integer                  :: nActiveParticles
    end function if_nActive

    subroutine if_updateStep(self, flow, dt)
      import pa_t
      import flow_t
      import std
      class (pa_t),  intent(inout) :: self
      type (flow_t), intent(in)    :: flow
      real (std),    intent(in)    :: dt
    end subroutine if_updateStep

    subroutine if_destroy(self)
      import pa_t
      class (pa_t), intent(inout)  :: self
    end subroutine if_destroy

  end interface

  ! Factory method

  interface if_pa_factory
    module function pa_factory(nParticles, options) result(pa)
      integer,          intent(in) :: nParticles
      type (options_t), intent(in) :: options
      class (pa_t),     pointer    :: pa
    end function pa_factory
  end interface if_pa_factory

  public :: pa_factory

end module pa
