module source

  ! A simple source type which can produce particles at a given rate.
  ! It also defines the number of species.

#include "assertion.h"

  use assertion
  use constants
  use options
  use rtime

  implicit none
  private

  type, public :: source_t
    integer                :: iD             ! source iD
    integer                :: nSpecies       ! for masses(nSpecies, nParticles)
    type (calendar_t)      :: cal            ! calendar
    type (sTime_t)         :: startTime      ! start producing particles
    type (sTime_t)         :: stopTime       ! stop
    type (sTimeInterval_t) :: tperiod        ! period (top-hat only)
    real (std)             :: rfraction      ! scale factor for release
    real (std)             :: x(3)           ! position
    real (std)             :: releaseRate    ! particles per second
    integer (int64)        :: iUP0           ! base iUP for this source
    integer (int32)        :: seed0(4)       ! base RNG seed for this source
  end type source_t

  public :: sourceInstantaneous
  public :: sourcePiecewiseConstant
  public :: sourceTopHat
  public :: sourceReleaseInterval

  integer, parameter, private :: SOURCE_TYPE_INSTANTANEOUS      = 1
  integer, parameter, private :: SOURCE_TYPE_PIECEWISE_CONSTANT = 2
  integer, parameter, private :: SOURCE_TYPE_TOP_HAT            = 3

contains

  !---------------------------------------------------------------------------
  !---------------------------------------------------------------------------

  function sourceInstantaneous(options) result(source)

    ! All particles are release at startTime == stopTime == 0.
    ! The maximum number of particles possible should be released
    ! in one go. (Specified elsewhere.)

    ! The source position is somewhat arbitrary. However, should not
    ! be too close to the ground.

    type (options_t), intent(in) :: options
    type (source_t)              :: source

    source%iD          = SOURCE_TYPE_INSTANTANEOUS
    source%nSpecies    = options%nSpecies
    source%cal         = calendar_t()         ! default
    source%startTime   = sTimeZero()
    source%stopTime    = sTimeZero()
    source%x(1)        = 0.0
    source%x(2)        = 0.0
    source%x(3)        = 300.0
    source%releaseRate = 1.0*options%maxLocal
    source%iUP0        = 0
    source%seed0       = (/ 11111111, 22222222, 33333333, 44444444 /)

  end function sourceInstantaneous

  !----------------------------------------------------------------------------

  function sourcePiecewiseConstant(options, t0, t1) result(source)

    type (options_t), intent(in) :: options
    type (sTime_t),   intent(in) :: t0
    type (sTime_t),   intent(in) :: t1
    type (source_t)              :: source

    source%iD          = SOURCE_TYPE_PIECEWISE_CONSTANT
    source%nSpecies    = options%nSpecies
    source%cal         = calendar_t()
    source%startTime   = t0
    source%stopTime    = t1
    source%x(1)        = 0.0
    source%x(2)        = 0.0
    source%x(3)        = 300.0
    source%releaseRate = options%maxLocal ! see below
    source%iUP0        = 0
    source%seed0       = (/ 11111111, 22222222, 33333333, 44444444 /)

    block
      ! Compute the extent of source activity in seconds, which goes
      ! into the release rate computation
      type (sTimeInterval_t) :: dt
      real (std) :: ss
      dt = sTimeSubtract(source%cal, t1, t0)
      ss = sTimeIntervalToReal(source%cal, dt, 1.0_std)
      assert(ss > 0.0)
      source%releaseRate = options%maxLocal / ss
    end block

  end function sourcePiecewiseConstant

  !----------------------------------------------------------------------------

  function sourceTopHat(options,  tperiod, rfraction) result(source)

    ! Source with top-hat profile in time of given period.
    ! Source will start at "on" at t = 0, and switch to "off"
    ! a tperiod/2, and back to "on" again at tperiod, and so on.

    ! For each "on" episode, the source should produce rfraction
    ! of the maximum number of particles specified by the options
    ! (0 < rfraction, and probably rfraction <= 1.0).

    ! The source will release indefinitely.

    type (options_t),       intent(in) :: options
    type (sTimeInterval_t), intent(in) :: tperiod
    real (std),             intent(in) :: rfraction
    type (source_t)                    :: source

    source%iD          = SOURCE_TYPE_TOP_HAT
    source%nSpecies    = options%nSpecies
    source%cal         = calendar_t()
    source%startTime   = sTimeZero()
    source%stopTime    = sTimeInfFuture(source%cal)
    source%tperiod     = tperiod
    source%rfraction   = rfraction
    source%x(1)        = 0.0
    source%x(2)        = 0.0
    source%x(3)        = 100.0            ! nearer the bottom here
    source%releaseRate = options%maxLocal ! see below
    source%iUP0        = 0
    source%seed0       = (/ 11111111, 22222222, 33333333, 44444444 /)

    block
      ! Compute the extent of source activity in seconds, which goes
      ! into the release rate computation
      real (std) :: ss
      ss = sTimeIntervalToReal(source%cal, tperiod, 1.0_std)
      ss = ss / 2.0 ! Release is over half the period
      assert(ss > 0.0)
      source%releaseRate = rfraction*options%maxLocal / ss
    end block

  end function sourceTopHat

  !----------------------------------------------------------------------------

  function sourceReleaseInterval(source, t0, deltat) result(nRelease)

    ! Compute the number of particles to be released between t0 and t0 + dt.
    ! This assumes all real times, and times are positive.

    type (source_t),        intent(in) :: source
    type (sTime_t),         intent(in) :: t0
    type (sTimeInterval_t), intent(in) :: deltat
    integer                            :: nRelease

    nRelease = 0

    if (source%iD == SOURCE_TYPE_INSTANTANEOUS) then
      if (sTimeEqual(t0, sTimeZero())) nRelease = nint(source%releaseRate)
    end if

    if (source%iD == SOURCE_TYPE_PIECEWISE_CONSTANT) then
      block
        real (std) :: tstart
        real (std) :: tstop
        real (std) :: tmin
        real (std) :: tmax
        real (std) :: dtsource
        real (std) :: dtrelease
        type (sTime_t)         :: t1
        type (sTimeInterval_t) :: sdt

        ! source active between [tstart, tstop] (seconds)
        sdt    = sTimeSubtract(source%cal, source%startTime, sTimeZero())
        tstart = sTimeIntervalToReal(source%cal, sdt, 1.0_std)

        sdt    = sTimeSubtract(source%cal, source%stopTime, sTimeZero())
        tstop  = sTimeIntervalToReal(source%cal, sdt, 1.0_std)

        ! work out overlap between this and time of interest [t0, t0 + dt]
        sdt    = sTimeSubtract(source%cal, t0, sTimeZero())
        tmin   = max(tstart, sTimeIntervalToReal(source%cal, sdt, 1.0_std))

        t1     = sTimeAdd(source%cal, t0, deltat)
        sdt    = sTimeSubtract(source%cal, t1, sTimeZero())
        tmax   = min(tstop,  sTimeIntervalToReal(source%cal, sdt, 1.0_std))

        dtsource  = tstop - tstart
        dtrelease = max(0.0_std, min(tmax - tmin, dtsource))

        nRelease  = nint(source%releaseRate*dtrelease)
      end block
    end if

    if (source%iD == SOURCE_TYPE_TOP_HAT) then
      ! Overlap between d + dt and nearest "on" episode
      ! We assume dt << tperiod so there's only a single "on" episode
      ! involved.

      ! This looks rather more painful than it should ...
      block
        type (sTimeInterval_t) :: sdt
        type (sTime_t)         :: t1
        real (std) :: tp
        real (std) :: ts, te
        real (std) :: tmin, tmax
        real (std) :: dtsource
        real (std) :: dtrelease
        integer    :: np

        ! Determine start (ts) and end (te) of "on" period relevant for t0
        tp  = sTimeIntervalToReal(source%cal, source%tperiod, 1.0_std)
        sdt = sTimeSubtract(source%cal, t0, sTimeZero())
        np  = sTimeIntervalDivide(sdt, source%tperiod)
        ts  = (1.0*np      )*tp
        te  = (1.0*np + 0.5)*tp
        dtsource = te - ts

        ! Overlap between [ts, te] and [t0, t0 + dt]

        sdt    = sTimeSubtract(source%cal, t0, sTimeZero())
        tmin   = max(ts, sTimeIntervalToReal(source%cal, sdt, 1.0_std))

        t1     = sTimeAdd(source%cal, t0, deltat)
        sdt    = sTimeSubtract(source%cal, t1, sTimeZero())
        tmax   = min(te,  sTimeIntervalToReal(source%cal, sdt, 1.0_std))

        dtrelease = max(0.0_std, min(tmax - tmin, dtsource))
        nRelease  = nint(source%releaseRate*dtrelease)
      end block
    end if

  end function sourceReleaseInterval

end module source
