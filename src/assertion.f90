module assertion

  ! A C-like assert() macro is provided via C-preprocessor directives
  ! in the file assertion.h. The actual implementation is via Fortran
  ! stop in the external procedure assertion_stop().
  !
  ! A wrapped version assertion() declared with a pure interface is
  ! provided to allow assertions to be placed in pure procedures.
  ! This is a sleight-of-hand which may generate a compiler warning
  ! (e.g., with Cray ftn).
  !
  ! The assertion_stop() subroutine implementation is deferred to a
  ! separate file to prevent the compiler noticing that the interface
  ! is `pure` but implementation is not. (assertion_stop.f90)

  implicit none
  public

contains

  !---------------------------------------------------------------------------
  !---------------------------------------------------------------------------

  pure subroutine assertion_impl(expr, file, line)

    ! Fail if expr is not .true.

    !$omp declare target

    logical,             intent(in) :: expr
    character (len = *), intent(in) :: file
    integer,             intent(in) :: line

    ! We have to simulate a pure subroutine.
    ! From an idea of Harold Anlauf.
    interface
      pure subroutine assertion_stop(expr, file, line)
        !$omp declare target
        logical,             intent(in) :: expr
        character (len = *), intent(in) :: file
        integer,             intent(in) :: line
      end subroutine assertion_stop
    end interface

    if (expr) return
    call assertion_stop(expr, file, line)

  end subroutine assertion_impl

end module assertion

