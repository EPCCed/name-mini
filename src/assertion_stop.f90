
! This implementation of assertion_stop() is kept in a separate file to the
! explicit interface in assertion.f90. This prevents the compiler noticing that
! the interface is `pure` but the implementation is not.

subroutine assertion_stop(expr, file, line)

  use, intrinsic :: iso_fortran_env
  implicit none

  !$omp declare target device_type(any)

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
  write (error_unit, *) "Assertion failed. [File, Line]: ", file, line

  stop

end subroutine assertion_stop
