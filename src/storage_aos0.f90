module storage_aos0

  ! AOS particle storage following the standard approach
  ! excepting that the random number state is:
  !
  !   integer (int32) :: s(4)

#include "assertion.h"

  use assertion
  use constants
  use maths
  use rtime
  use source

  implicit none
  private

  type, public :: aos0_t
    integer (int64)  :: iUP                ! an index unique to the particle
    integer          :: iE                 ! Store particle extra index.
    real (std)       :: X(3)               ! Particle position.
    real (std)       :: XOld(3)            ! Position at previous time step.

    type (sTime_t)   :: T                  ! Travel time.
    type (sTime_t)   :: Told               ! Travel time at previous time step.
    type (sTime_t)   :: T0                 ! Release time.
    logical          :: NeedToSetVel       !
    logical          :: Active             ! Indicates particle is active.
    integer          :: Marked             ! 

    integer          :: iHCoord            !
    integer          :: iZCoord            ! vertical coordinate look-up index
    integer          :: iHCoordOld         ! Horizontal coords for xOld
    integer          :: iZCoordOld         !
    integer          :: iSource            ! Index of source of the particle
    real (std)       :: H                  ! Boundary layer depth (in metres).

    ! Particulates have ...
    real (std)        :: WSed              ! Gravitational settling velocity.
    real (std)        :: Diameter          ! Diameter
    real (std)        :: Density           ! Density
    real (std)        :: ParticleShape     ! Particle Shape
    real (std)        :: Elongation        ! Elongation
    real (std)        :: Flatness          ! Flatness
    integer           :: ShapeSchemeCode   ! Shape Scheme

    integer (int32)   :: s(4)              ! random state
  end type aos0_t

  public :: storage_aos0_stateZero
  public :: storage_aos0_initialiseFromSource

contains

  function storage_aos0_stateZero(iUP, seed0) result(particle)

    ! Initialise an return a particle with unique particle index iUP
    ! and random number generator (RNG) state seed0.

    ! All other state is "defaulted", including iE, which is expected
    ! to provide an indirection to default Extra_ information.

    ! The default active flag is .true.

    integer (int64),  intent(in) :: iUP
    integer (int32),  intent(in) :: seed0(4)
    type (aos0_t)                :: particle

    particle%iUP     = iUP
    particle%iE      = 0
    particle%x(:)    = 0.0
    particle%xold(:) = 0.0

    particle%t    = sTimeZero()
    particle%tOld = particle%t
    particle%t0   = particle%t

    particle%needToSetVel = .false.
    particle%active       = .true.
    particle%marked       = 0

    particle%iHCoord    = 0
    particle%iZCoord    = 0
    particle%iHCoordOld = 0
    particle%iZCoordOld = 0
    particle%iSource    = 0

    particle%h        = 0.0
    particle%wSed     = 0.0
    particle%diameter = 0.0
    particle%density  = 0.0

    particle%particleShape   = 0.0
    particle%shapeSchemeCode = 0

    particle%s(1:4) = seed0(1:4)

  end function storage_aos0_stateZero

  !---------------------------------------------------------------------------
  !---------------------------------------------------------------------------

  pure function storage_aos0_initialiseFromSource(iUP, source) result(particle)

    !$omp declare target

    integer (int64), intent(in) :: iUP
    type (source_t), intent(in) :: source
    type (aos0_t)               :: particle

    particle%iUP             = iUP
    particle%iE              = 0
    particle%x(1:3)          = source%x(1:3)
    particle%xold(1:3)       = source%x(1:3)

    particle%t               = sTimeZero()
    particle%tOld            = particle%t
    particle%t0              = particle%t

    particle%needToSetVel    = .false.
    particle%active          = .true.
    particle%marked          = 0

    particle%iHCoord         = 0
    particle%iZCoord         = 0
    particle%iHCoordOld      = 0
    particle%iZCoordOld      = 0
    particle%iSource         = 0

    particle%h               = 0.0
    particle%wSed            = 0.0
    particle%diameter        = 0.0
    particle%density         = 0.0

    particle%particleShape   = 0.0
    particle%shapeSchemeCode = 0

    particle%s(1:4)          = mathsRNGStateAdvanceFrom(source%seed0, iUP)

  end function storage_aos0_initialiseFromSource

end module storage_aos0
