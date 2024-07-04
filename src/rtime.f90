module rtime

  ! An experimental replacement short time type which
  !
  !  1. decouples date-like and interval-like time structures into
  !     separate sTime_t and sTimeInterval_t
  !  2. avoids run time failures and messages which would preclude
  !     device versions via !$omp target device_type(any)
  !  3. has a re-entrant version free of global state for each
  !     function (and a non re-entrant version which uses a global
  !     calendar_t object).
  !
  ! Note that not all functions require device versions at the
  ! present time.

#include "assertion.h"

  use assertion
  use constants

  implicit none
  private

  ! A reduced "global state" for which there are four possibilities
  ! No relative time frame is provided

  type, public :: calendar_t
    private
    logical :: isBackwards = .false.      ! Arrow of time
    logical :: isGregorian = .false.      ! Otherwise 360 day year
    logical :: initialised = .false.      ! Must be true for time operations
  end type calendar_t

  public :: calendar_initialise
  public :: calendar_finalise

  ! Class static calendar

  type (calendar_t), private :: calendar_
  !$omp declare target (calendar_)

  ! Short time. This really must be opaque, as the implementation
  !             should change.

  type, public :: sTime_t
    private
    integer (int64) :: fractions          ! 10^-7 sec since reference time
  end type sTime_t

  public :: sTimeZero
  public :: sTimeInfFuture
  public :: sTimeEqual
  public :: sTimeAdd
  public :: sTimeSubtract

  interface sTimeGreaterThanOrEqual
    module procedure sTimeGreaterThanOrEqual    ! Explicit calendar arg
    module procedure sTimeGreaterThanOrEqual_   ! Uses calendar_
  end interface sTimeGreaterThanOrEqual

  interface sTimeSubtract
    module procedure sTimeSubtract
    module procedure sTimeSubtract_
  end interface sTimeSubtract

  !---------------------------------------------------------------------------

  ! A separate short time interval.

  type, public :: sTimeInterval_t
    private
    integer (int64)   :: fractions
  end type sTimeInterval_t

  ! "Constructors"

  interface sTimeInterval_t
    module procedure sTimeIntervalreal32
    module procedure sTimeIntervalreal64
  end interface

  public :: sTimeIntervalToReal
  interface sTimeIntervalToReal
    module procedure sTimeIntervalToReal32
    module procedure sTimeIntervalToReal64
    module procedure sTimeIntervalToReal32_
    module procedure sTimeIntervalToReal64_
  end interface sTimeIntervalToReal

  public :: sTimeIntervalDivide

  ! We will use the equivalent of INT_MIN and INT MAX to represent
  ! +/- infinity. Note that minimum is -huge() - 1 (two's complement).
  ! These are labelled STIME_INF_PAST and STIME_INF_FUTURE.
  !
  ! To avoid flirting with integer overflow or underflow, STIME_MIN
  ! and STIME_MAX are limits for conversion from real64 values.

  integer (int64), parameter :: STIME_FRACS_PER_SEC = 10000000
  integer (int64), parameter :: STIME_INF_PAST      = -huge(0_int64) - 1
  integer (int64), parameter :: STIME_INF_FUTURE    = +huge(0_int64)

  real (real64),   parameter :: STIME_MAX           = 0.5*STIME_INF_FUTURE
  real (real64),   parameter :: STIME_MIN           = 0.5*STIME_INF_PAST

contains

  !---------------------------------------------------------------------------
  !---------------------------------------------------------------------------

  subroutine calendar_initialise(isBackwards, isGregorian)

    logical, intent(in), optional :: isBackwards
    logical, intent(in), optional :: isGregorian

    assert(calendar_%initialised .eqv. .false.)

    calendar_ = Calendar_t()
    if (present(isBackwards)) calendar_%isBackwards = isBackwards
    if (present(isGregorian)) calendar_%isGregorian = isGregorian

    calendar_%initialised = .true.
    
  end subroutine calendar_initialise

  !----------------------------------------------------------------------------
  !----------------------------------------------------------------------------

  subroutine calendar_finalise()

    calendar_ = Calendar_t()

    assert(calendar_%initialised .eqv. .false.)

  end subroutine calendar_finalise

  !----------------------------------------------------------------------------
  !----------------------------------------------------------------------------

  pure function calendars_same(cal1, cal2) result(same)

    type (Calendar_t), intent(in) :: cal1
    type (Calendar_t), intent(in) :: cal2
    logical                       :: same

    !$omp declare target

    same = .true.
    if (cal1%isBackwards .neqv. cal2%isBackwards) same = .false.
    if (cal1%isGregorian .neqv. cal2%isGregorian) same = .false.

  end function calendars_same

  !---------------------------------------------------------------------------
  !---------------------------------------------------------------------------

  pure function sTimeIsInfPast(calendar, st) result(infinitePast)

    type (calendar_t), intent(in) :: calendar
    type (sTime_t),    intent(in) :: st
    logical                       :: infinitePast

    infinitePast = (st%fractions == STIME_INF_PAST)
    if (calendar%isBackwards) then
      infinitePast = (st%fractions == STIME_INF_FUTURE)
    end if

  end function sTimeIsInfPast

  !---------------------------------------------------------------------------

  pure function sTimeIsInfFuture(calendar, st) result(infiniteFuture)

    type (calendar_t), intent(in) :: calendar
    type (sTime_t),    intent(in) :: st
    logical                       :: infiniteFuture

    infiniteFuture = (st%fractions == STIME_INF_FUTURE)
    if (calendar%isBackwards) then
      infiniteFuture = (st%fractions == STIME_INF_PAST)
    end if

  end function sTimeIsInfFuture

  !---------------------------------------------------------------------------

  pure function sTimeIsInf(st) result(infinite)

    ! Don't care about the calandar (same forwards or backwards)

    type (sTime_t), intent(in) :: st
    logical                    :: infinite

    infinite = (st%fractions == STIME_INF_PAST)   .or. &
               (st%fractions == STIME_INF_FUTURE)

  end function sTimeIsInf

  !---------------------------------------------------------------------------

  pure function sTimeInfPast(calendar) result(st)

    type (calendar_t), intent(in) :: calendar
    type (sTime_t)                :: st

    st%fractions = STIME_INF_PAST
    if (calendar%isBackwards) st%fractions = STIME_INF_FUTURE

  end function sTimeInfPast

  !---------------------------------------------------------------------------

  pure function sTimeInfPast_() result(st)

    type (sTime_t) :: st

    assert(calendar_%initialised)

    st = sTimeInfPast(calendar_)

  end function sTimeInfPast_

  !---------------------------------------------------------------------------

  pure function sTimeInfFuture(calendar) result(st)

    type (calendar_t), intent(in) :: calendar
    type (sTime_t)                :: st

    st%fractions = STIME_INF_FUTURE
    if (calendar%isBackwards) st%fractions = STIME_INF_PAST

  end function sTimeInfFuture

  !---------------------------------------------------------------------------

  pure function sTimeInfFuture_() result(st)

    type (sTime_t) :: st

    assert(calendar_%initialised)

    st = sTimeInfFuture(calendar_)

  end function sTimeInfFuture_

  !---------------------------------------------------------------------------

  pure function sTimeEqual(t1, t2) result(equal)

    ! Are the two times equal (don't care about the calendar here).

    type (sTime_t), intent(in) :: t1
    type (sTime_t), intent(in) :: t2
    logical                    :: equal

    !$omp declare target

    equal = (t1%fractions == t2%fractions)

  end function sTimeEqual

  !---------------------------------------------------------------------------

  pure function sTimeGreaterThanOrEqual(calendar, t1, t2) result(gte)

    ! Return true if t1 >= t2 in the calendar supplied

    type (calendar_t), intent(in) :: calendar
    type (sTime_t),    intent(in) :: t1
    type (sTime_t),    intent(in) :: t2
    logical                       :: gte

    !$omp declare target

    if (calendar%isBackwards) then
      ! Backwards ...
      gte = (t1%fractions <= t2%fractions)
    else
      ! As one expects.
      gte = (t1%fractions >= t2%fractions)
    end if

  end function sTimeGreaterThanOrEqual

  !---------------------------------------------------------------------------

  pure function sTimeGreaterThanOrEqual_(t1, t2) result(gte)

    ! Return true if t1 >= t2 using the static calendar_

    type (sTime_t), intent(in) :: t1
    type (sTime_t), intent(in) :: t2
    logical                    :: gte

    !$omp declare target

    assert(calendar_%initialised)

    gte = sTimeGreaterThanOrEqual(calendar_, t1, t2)

  end function sTimeGreaterThanOrEqual_

  !---------------------------------------------------------------------------

  pure function sTimeZero() result(t)

    type (sTime_t) :: t

    !$omp declare target

    t%fractions  = 0

  end function sTimeZero

  !---------------------------------------------------------------------------

  function sTimeAdd(cal, t, dt) result(tsum)

    ! tsum = t + dt (with various cases related to infinities ...)

    type (calendar_t),      intent(in)    :: cal
    type (sTime_t),         intent(in)    :: t
    type (sTimeInterval_t), intent(in)    :: dt
    type (sTime_t)                        :: tsum

    !$omp declare target device_type(any)

    ! Not +inf + -inf or -inf + +inf ("cancelling infinities")
    ! PENDING "programmer error" or "run time error"?
    ! PENDING calendar unused if no assertions

    assert(.not.(sTimeIsInfPast(cal,t)   .and. sTimeIntervalIsPlusInf(cal,dt)))
    assert(.not.(sTimeIsInfFuture(cal,t) .and. sTimeIntervalIsMinusInf(cal,dt)))

    tsum = sTime_t(0)

    if (sTimeIsInf(t)) then
      tsum%fractions = t%fractions
    else if (sTimeIntervalIsInf(dt)) then
      tsum%fractions = dt%fractions
    else
      ! finite case
      tsum%fractions = t%fractions + dt%fractions
    end if

  end function sTimeAdd

  !---------------------------------------------------------------------------
  !---------------------------------------------------------------------------

  pure function sTimeSubtract(cal, t1, t2) result(tdifference)

    ! tdifference = t1 minus t2
    ! +inf - +inf or -inf - -inf are not allowed ("cancelling infinities")

    type (calendar_t), intent(in) :: cal
    type (sTime_t),    intent(in) :: t1
    type (sTime_t),    intent(in) :: t2
    type (sTimeInterval_t)        :: tdifference

    ! The calendar only affects the identification of +/- inf
    assert(.not. (sTimeIsInfFuture(cal, t1) .and. sTimeIsInfFuture(cal, t2)))
    assert(.not. (sTimeIsInfPast(cal, t1)   .and. sTimeIsInfPast(cal, t2)))

    tDifference = sTimeInterval_t(0)

    if (sTimeIsInf(t1)) then
      tdifference%fractions = t1%fractions
    else if (sTimeIsInf(t2)) then
      tdifference% fractions = t2%fractions
    else
      tdifference%fractions = t1%fractions - t2%fractions
    end if

  end function sTimeSubtract

  !---------------------------------------------------------------------------

  pure function sTimeSubtract_(t1, t2) result(tdifference)

    type (sTime_t),    intent(in) :: t1
    type (sTime_t),    intent(in) :: t2
    type (sTimeInterval_t)        :: tdifference

    assert(calendar_%initialised)

    tdifference = sTimeSubtract(calendar_, t1, t2)

  end function sTimeSubtract_

  !---------------------------------------------------------------------------
  !---------------------------------------------------------------------------

  pure function sTimeIntervalMinusInf(calendar) result(sti)

    type (calendar_t), intent(in) :: calendar
    type (sTimeInterval_t)        :: sti

    sti%fractions = STIME_INF_PAST
    if (calendar%isBackwards) sti%fractions = STIME_INF_FUTURE

  end function sTimeIntervalMinusInf

  !---------------------------------------------------------------------------

  pure function sTimeIntervalPlusInf(calendar) result(sti)

    type (calendar_t), intent(in) :: calendar
    type (sTimeInterval_t)        :: sti

    sti%fractions = STIME_INF_FUTURE
    if (calendar%isBackwards) sti%fractions = STIME_INF_PAST

  end function sTimeIntervalPlusInf

  !---------------------------------------------------------------------------

  pure function sTimeIntervalIsMinusInf(calendar, sti) result(minusInf)

    type (calendar_t),      intent(in) :: calendar
    type (sTimeInterval_t), intent(in) :: sti
    logical                            :: minusInf

    minusInf = (sti%fractions == STIME_INF_PAST)
    if (calendar%isBackwards) then
      minusInf = (sti%fractions == STIME_INF_FUTURE)
    end if

  end function sTimeIntervalIsMinusInf

  !---------------------------------------------------------------------------

  pure function sTimeIntervalIsPlusInf(calendar, sti) result(plusInf)

    type (calendar_t),      intent(in) :: calendar
    type (sTimeInterval_t), intent(in) :: sti
    logical                            :: plusInf

    plusInf = (sti%fractions == STIME_INF_FUTURE)
    if (calendar%isBackwards) then
      plusInf = (sti%fractions == STIME_INF_PAST)
    end if

  end function sTimeIntervalIsPlusInf

  !---------------------------------------------------------------------------

  pure function sTimeIntervalIsInf(sti) result(infinite)

    ! Same in forwards or backwards calendar

    type (sTimeInterval_t), intent(in) :: sti
    logical                            :: infinite

    infinite = (sti%fractions == STIME_INF_PAST)   .or. &
               (sti%fractions == STIME_INF_FUTURE)

  end function sTimeIntervalIsInf

  !---------------------------------------------------------------------------

  pure function sTimeIntervalreal32(calendar, seconds) result(sti)

    type (calendar_t), intent(in) :: calendar
    real (real32),     intent(in) :: seconds
    type (sTimeInterval_t)        :: sti

    real (real64) :: s64

    s64 = real(seconds, real64)
    sti = sTimeIntervalreal64(calendar, s64)

  end function sTimeIntervalreal32

  !---------------------------------------------------------------------------

  pure function sTimeIntervalreal64(calendar, seconds) result(sti)

    ! Convert seconds to an sTimeInterval_t in the current calendar

    type (calendar_t), intent(in) :: calendar
    real (real64),     intent(in) :: seconds
    type (sTimeInterval_t)        :: sti

    real (real64) :: s

    s = seconds*STIME_FRACS_PER_SEC

    sti%fractions = int(s, int64)
    if (calendar%isBackwards) sti%fractions = -sti%fractions

    if (s > STIME_MAX) sti = sTimeIntervalPlusInf(calendar)
    if (s < STIME_MIN) sti = sTimeIntervalMinusInf(calendar)

  end function sTimeIntervalreal64

  !---------------------------------------------------------------------------

  pure function sTimeIntervalToReal32(cal, stime, rkind) result(t)

    type (calendar_t),      intent(in) :: cal
    type (sTimeInterval_t), intent(in) :: stime
    real (real32),          intent(in) :: rkind
    real (real32)                      :: t

    assert(.not. sTimeIntervalIsInf(stime))

    t = real(stime%fractions, real32) / real(STIME_FRACS_PER_SEC, real32)

    if (cal%isBackwards) t = -t

  end function sTimeIntervalToReal32

  !---------------------------------------------------------------------------

  pure function sTimeIntervalToReal64(cal, stime, rkind) result(t)

    ! Return time interval in seconds

    type (calendar_t),      intent(in) :: cal
    type (sTimeInterval_t), intent(in) :: stime
    real (real64),          intent(in) :: rkind
    real (real64)                      :: t

    assert(.not. sTimeIntervalIsInf(stime))

    t = real(stime%fractions, real64) / real(STIME_FRACS_PER_SEC, real64)

    if (cal%isBackwards) t = -t

  end function sTimeIntervalToReal64

  !----------------------------------------------------------------------------

  pure function sTimeIntervalToReal32_(stime, rkind) result(t)

    type (sTimeInterval_t), intent(in) :: stime
    real (real32),          intent(in) :: rkind
    real (real32)                      :: t

    assert(calendar_%initialised)

    t = sTimeIntervalToReal32(calendar_, stime, rkind)

  end function sTimeIntervalToReal32_

  !----------------------------------------------------------------------------

  pure function sTimeIntervalToReal64_(stime, rkind) result(t)

    type (sTimeInterval_t), intent(in) :: stime
    real (real64),          intent(in) :: rkind
    real (real64)                      :: t

    assert(calendar_%initialised)

    t = sTimeIntervalToReal64(calendar_, stime, rkind)

  end function sTimeIntervalToReal64_

  !----------------------------------------------------------------------------

  pure function sTimeIntervalDivide(dt1, dt2) result(iratio)

    ! Divide one interval by another (integer division)

    type (sTimeInterval_t), intent(in) :: dt1
    type (sTimeInterval_t), intent(in) :: dt2
    integer (int32)                    :: iratio

    assert(.not. sTimeIntervalIsInf(dt1))
    assert(.not. sTimeIntervalIsInf(dt2))
    assert(dt2%fractions /= 0)

    iratio = int(dt1%fractions/dt2%fractions, int32)

  end function sTimeIntervalDivide

end module rtime
