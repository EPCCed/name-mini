module maths

  ! Provides versions of uniform an Gaussian deistribted random
  ! numbers based on L'Ecuyer and Andres (1997) LCG.

  ! For the GPU version we have reverted to "plain old data" in
  ! the RNG state "integer (int32) :: s(4)" The mulitpliers and
  ! the moduli are inlined explicitly.

  ! nb. the qualifier "device_type(any)" in "!$omp declare target"
  ! is causing problems for some compilers at the moment and so is
  ! ommitted. All could be host/devive functions.

  ! (c) 2024 The University of Edinburgh

#include "assertion.h"

  use assertion
  use constants

  implicit none
  private

  public :: mathsRNGStateAdvanceFrom

  ! Uniform RNG interface
  public :: mathsRandomUniform
  public :: mathsRandomUniform64

  interface mathsRandomUniform
    module procedure mathsRandomUniform64
  end interface mathsRandomUniform

  ! Gaussian RNG interface

  public :: mathsRandomGaussian2_64
  public :: mathsRandomGaussian2_32
  public :: mathsRandomGaussian1_64
  public :: mathsRandomGaussian1_32

  public :: mathsRandomGaussian3
  public :: mathsRandomGaussian3_64
  public :: mathsRandomGaussian3_32

  interface mathsRandomGaussian3
    module procedure mathsRandomGaussian3_32
    module procedure mathsRandomGaussian3_64
  end interface mathsRandomGaussian3

  public :: mathsIsPowerOfTwo
  public :: mathsNearestPowerOfTwo

contains

  !---------------------------------------------------------------------------
  !---------------------------------------------------------------------------

  subroutine mathsRandomUniform64(s, u01)

    ! L'Ecuyer and Andres (1997).
    ! Return a uniformly distributed variate 0.0 <= u01 < 1.0.

    !$omp declare target

    integer (int32), intent(inout) :: s(4)
    real (real64),   intent(out)   :: u01

    integer (int64) :: s64

    assert(1 <= s(1) .and. s(1) < 2147483647)
    assert(1 <= s(2) .and. s(2) < 2147483543)
    assert(1 <= s(3) .and. s(3) < 2147483423)
    assert(1 <= s(4) .and. s(4) < 2147483323)

    s64 = s(1)
    s(1) = int(mod(45991 * s64, 2147483647_int64), int32)
    u01  = (1.0_real64/2147483647)*s(1)

    s64 = s(2)
    s(2) = int(mod(207707 * s64, 2147483543_int64), int32)
    u01 = u01 - (1.0_real64/2147483543)*s(2)
    if (u01 <  0.0_real64) u01 = u01 + 1.0_real64

    s64 = s(3)
    s(3) = int(mod(138556 * s64, 2147483423_int64), int32)
    u01 = u01 + (1.0_real64/2147483423)*s(3)
    if (u01 >= 1.0_real64) u01 = u01 - 1.0_real64

    s64 = s(4)
    s(4) = int(mod(49689 * s64, 2147483323_int64), int32)
    u01 = u01 - (1.0_real64/2147483323)*s(4)
    if (u01 <  0.0_real64) u01 = u01 + 1.0_real64

    assert(0.0 <= u01 .and. u01 < 1.0_real64)

  end subroutine mathsRandomUniform64

  !---------------------------------------------------------------------------
  !---------------------------------------------------------------------------

  subroutine mathsRandomGaussian3_64(s, rv)

    ! Convenience for three variates in one go

    !$omp declare target

    integer (int32), intent(inout) :: s(4)
    real (real64),   intent(out)   :: rv(3)

    real (real64) :: r(2)

    call mathsRandomGaussian2_64(s, r)
    rv(1:2) = r(1:2)
    call mathsRandomGaussian2_64(s, r)
    rv(3) = r(1) ! 1 variate discarded

  end subroutine mathsRandomGaussian3_64

  !---------------------------------------------------------------------------

  subroutine mathsRandomGaussian3_32(s, rv)

    ! Convenience for three variates in one go

    !$omp declare target

    integer (int32), intent(inout) :: s(4)
    real (real32),   intent(out)   :: rv(3)

    real (real64) :: r(2)

    call mathsRandomGaussian2_64(s, r)
    rv(1:2) = r(1:2)
    call mathsRandomGaussian2_64(s, r)
    rv(3) = r(1) ! 1 variate discarded

  end subroutine mathsRandomGaussian3_32

  !---------------------------------------------------------------------------
  !---------------------------------------------------------------------------

  subroutine mathsRandomGaussian2_64(s, rv)

    ! Return two random variates from Gaussian distribution with zero mean
    ! and unit variance.

    !$omp declare target

    integer (int32), intent(inout) :: s(4)
    real (real64),   intent(out)   :: rv(2)

    real (real64) :: r1, r2
    real (real64) :: rsq
    real (real64) :: factor

    ! Box-Mueller

    do
      call mathsRandomUniform64(s, r1)
      call mathsRandomUniform64(s, r2)

      r1  = 2.0*r1 - 1.0
      r2  = 2.0*r2 - 1.0
      rsq = r1*r1 + r2*r2
      if (rsq < 1.0 .and. rsq /= 0.0) exit
    end do

    factor = sqrt(-2.0 * log(rsq) / rsq)
    rv(1) = r1 * factor
    rv(2) = r2 * factor

  end subroutine mathsRandomGaussian2_64

  !----------------------------------------------------------------------------
  !----------------------------------------------------------------------------

  subroutine mathsRandomGaussian1_32(s, rv)

    ! An overloaded version of mathsRandomGaussian1().

    !$omp declare target

    integer (int32),  intent(inout) :: s(4)
    real (real32),    intent(out)   :: rv

    real (real64) :: r1, r2
    real (real64) :: rsq
    real (real64) :: factor

    ! Box-Mueller

    do
      call mathsRandomUniform64(s, r1)
      call mathsRandomUniform64(s, r2)

      r1 = 2.0 * r1 - 1.0
      r2 = 2.0 * r2 - 1.0
      rsq = r1 * r1 + r2 * r2
      if (rsq < 1.0 .and. rsq /= 0.0) exit
    end do

    factor  = sqrt(-2.0 * log(rsq) / rsq)
    rv = real(r1 * factor, real32)

  end subroutine mathsRandomGaussian1_32

  !----------------------------------------------------------------------------
  !----------------------------------------------------------------------------

  subroutine mathsRandomGaussian1_64(s, rv)

    ! An overloaded version of mathsRandomGaussian1().

    !$omp declare target

    integer (int32),  intent(inout) :: s(4)
    real (real64),    intent(out)   :: rv

    real (real64) :: r1, r2
    real (real64) :: rsq

    ! Box-Mueller
    do
      call mathsRandomUniform64(s, r1)
      call mathsRandomUniform64(s, r2)

      r1 = 2.0*r1 - 1.0
      r2 = 2.0*r2 - 1.0
      rsq = r1*r1 + r2*r2
      if (rsq < 1.0 .and. rsq /= 0.0) exit
    end do

    rv = r1*sqrt(-2.0*log(rsq)/rsq)

  end subroutine mathsRandomGaussian1_64

  !----------------------------------------------------------------------------
  !----------------------------------------------------------------------------

  subroutine mathsRandomGaussian2_32(s, rv)

    ! An overloaded version of mathsRandomGaussian2().

    !$omp declare target

    integer (int32), intent(inout) :: s(4)
    real (real32),   intent(out)   :: rv(2)

    real (real64) :: r64(2)

    call mathsRandomGaussian2_64(s, r64)
    rv(:) = real(r64(:), real32)

  end subroutine mathsRandomGaussian2_32

  !---------------------------------------------------------------------------
  !---------------------------------------------------------------------------

  pure function mathsIsPowerOfTwo(n) result(isPoT)

    ! Determine whether the integer argument is a power of two.
    ! If n = 2^k, then the binary representations of n and (n-1)
    ! are
    !  n   ...00010000...
    !  n-1 ...00001111...

    integer, intent(in) :: n
    logical             :: isPoT

    assert(n > 0)

    isPoT = (iand(n, n-1) == 0)

  end function mathsIsPowerOfTwo

  !----------------------------------------------------------------------------
  !----------------------------------------------------------------------------

  pure function mathsLEcuyerMultiplier(j) result(a)

    ! Return the multipliers needed to advance the state by k places,
    ! that is, a^k, with k = 2^j.

    !$omp declare target

    integer (int32),  intent(in) :: j
    integer (int64)              :: a(4)

    integer :: i

    a(1) = 45991
    a(2) = 207707
    a(3) = 138556
    a(4) = 49689
    do i = 1, j
      a(1) = mod(a(1)*a(1), 2147483647_int64)
      a(2) = mod(a(2)*a(2), 2147483543_int64)
      a(3) = mod(a(3)*a(3), 2147483423_int64)
      a(4) = mod(a(4)*a(4), 2147483323_int64)
    end do

  end function mathsLEcuyerMultiplier

  !---------------------------------------------------------------------------
  !---------------------------------------------------------------------------

  pure function mathsLEcuyerAdvance(a, s0) result(snew)

    ! Advance the state provided by k = 2^j places and return the new state.
    ! The multipliers "a"  can be generated via a call to
    !     mathLEcuyerMulitiplier(j)
    ! and using such multipliers, we get: s_{i+k} = (multi * s_i) mod m.

    !$omp declare target

    integer (int64),  intent(in) :: a(4)
    integer (int32),  intent(in) :: s0(4)
    integer (int32)              :: snew(4)

    integer (int64) :: s

    s = s0(1)
    s = mod(a(1)*s, 2147483647_int64)
    snew(1) = int(s, int32)

    s = s0(2)
    s = mod(a(2)*s, 2147483543_int64)
    snew(2) = int(s, int32)

    s = s0(3)
    s = mod(a(3)*s, 2147483423_int64)
    snew(3) = int(s, int32)

    s = s0(4)
    s = mod(a(4)*s, 2147483323_int64)
    snew(4) = int(s, int32)

  end function mathsLEcuyerAdvance

  !---------------------------------------------------------------------------
  !---------------------------------------------------------------------------

  pure function mathsRNGStateAdvanceFrom(seed0, iUP) result(s)

    ! Advance seed0 an arbitrary number of places iUP (blocks of 2^31)
    ! and return the new state

    !$omp declare target

    integer (int32),  intent(in) :: seed0(4)
    integer (int64),  intent(in) :: iUP
    integer (int32)              :: s(4)

    assert(iUP > 0)

    s(1:4) = seed0(1:4)

    block
      integer (int64)  :: m
      integer          :: j
      integer (int64)  :: a(4)

      m = iUP
      j = 0
      do
        if (mod(m, 2_int64) == 1) then
          a = mathsLEcuyerMultiplier(31 + j)
          s = mathsLEcuyerAdvance(a, s)
        end if
        m = m/2
        j = j + 1
        if (m <= 0) exit
      end do
    end block

  end function mathsRNGStateAdvanceFrom

  !----------------------------------------------------------------------------

  pure function mathsNearestPowerOfTwo(n) result(k)

    ! For integer n, return k such that k = 2^m >= n; requires n > 0
    ! m = ceiling (log_10 n / log_10 2)

    integer, intent(in) :: n
    integer             :: k

    assert(n > 0)

    k = ceiling(log10(1.0*n) / log10(2.0))
    k = 2**k

    assert(k >= n)

  end function mathsNearestPowerOfTwo

end module maths
