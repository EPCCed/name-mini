
! This implementation of assertion_stop() is kept in a separate file to the
! explicit interface in assertion.f90. This prevents the compiler noticing that
! the interface is `pure` but the implementation is not.

subroutine assertion_stop(expr, file, line)

  use, intrinsic :: iso_fortran_env
  implicit none

  !$omp declare target

  ! Stop if expr is not .true.

  logical,             intent(in) :: expr
  character (len = *), intent(in) :: file
  integer,             intent(in) :: line

  if (expr) return

  ! Device notes ...
  ! no fmt available;
  ! no intrinsic trim() available;
  ! no string catenations via '//' are available;
  ! A bare list-directed write is all that is available ...

#ifdef __INTEL_LLVM_COMPILER
  ! Intel does not like strings
  ! error #5623: **Internal compiler error: internal abort**
  ! So, ...
  print *, "Assertion failed. [Line] ", line
#else
  write (error_unit, *) "Assertion failed. [File, Line]: ", file, line
#endif

  stop

end subroutine assertion_stop
