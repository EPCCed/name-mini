module util_omp_lib

#ifdef _OPENMP
  use omp_lib
  implicit none
  public
#else

  use iso_fortran_env
  implicit none
  public

contains

  function omp_get_max_threads() result(nthreads)
    integer :: nthreads
    nthreads = 1
  end function omp_get_max_threads

  function omp_get_num_devices() result(ndevice)
    integer :: ndevice
    ndevice = 0
  end function omp_get_num_devices

  function omp_get_thread_num() result(tid)
    integer :: tid
    tid = 0
  end function omp_get_thread_num

  function omp_get_num_threads() result(nthreads)
    integer :: nthreads
    nthreads = 1
  end function omp_get_num_threads

  function omp_get_initial_device() result(idevice)
    integer :: idevice
    idevice = -1
  end function omp_get_initial_device

  function omp_get_default_device() result(idevice)
    integer :: idevice
    idevice = 0
  end function omp_get_default_device

  function omp_get_team_num() result(iteam)
    integer :: iteam
    iteam = 0
  end function omp_get_team_num

  function omp_get_num_teams() result(nteams)
    integer :: nteams
    nteams = 1
  end function omp_get_num_teams

  function omp_get_wtime() result(t)
    real (real64) :: t
    call cpu_time(t)
  end function omp_get_wtime

  function omp_is_initial_device() result(initial)

    logical :: initial

    initial = .true.

  end function omp_is_initial_device

#endif

end module util_omp_lib
