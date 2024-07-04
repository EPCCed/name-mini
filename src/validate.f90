module validate

  ! Utility to help to output quantities for validation

  use constants
  use util_omp_lib

  implicit none
  private

  public :: validate_filename
  public :: validate_write

contains

  function validate_filename(stub, nParticles) result(filename)

    ! E.g., stub-gpu-10-001000000.dat

    character (len = *), intent(in)  :: stub
    integer,             intent(in)  :: nParticles
    character (len = :), allocatable :: filename

    character (len = 3)              :: nthreads
    character (len = 9)              :: nparts

    write (nthreads, '(i3.3)') omp_get_max_threads()
    write (nparts,   '(i9.9)') nparticles

    if (omp_get_num_devices() > 0) then
      filename = trim(stub)//"-gpu-"//nthreads//"-"//nparts//".dat"
    else
      filename = trim(stub)//"-cpu-"//nthreads//"-"//nparts//".dat"
    end if

  end function validate_filename

  !---------------------------------------------------------------------------

  subroutine validate_write(nparticles, xdata, my_unit, ierr)

    integer,    intent(in)  :: nparticles
    real (std), intent(in)  :: xdata(nparticles, 3)
    integer,    intent(in)  :: my_unit
    integer,    intent(out) :: ierr

    character (len = *), parameter :: afmt32 = "(i8,3(1x,e14.7))"
    character (len = *), parameter :: afmt64 = "(i8,3(1x,e22.15))"

    character (len = *), parameter :: bfmt32 = "(i8,1x,a1,1x,b32.32)"
    character (len = *), parameter :: bfmt64 = "(i8,1x,a1,1x,b64.64)"

    character (len = :), allocatable :: afmt
    character (len = :), allocatable :: bfmt

    integer :: iP

    ierr = 0

    ! Format for either real32 or real64
    afmt = afmt32
    bfmt = bfmt32
    if (std == real64) then
      afmt = afmt64
      bfmt = bfmt64
    end if

    ! ascii first
    do iP = 1, nparticles
      write (my_unit, afmt, err = 999) iP, xdata(iP, :)
    end do

    ! then binary
    do iP = 1, nparticles
      write (my_unit, bfmt, err = 999) iP, "x", xdata(iP, 1)
      write (my_unit, bfmt, err = 999) iP, "y", xdata(iP, 2)
      write (my_unit, bfmt, err = 999) iP, "z", xdata(iP, 3)
    end do

    return

999 continue
    ierr = 999

  end subroutine validate_write
     
end module validate
