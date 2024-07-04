program main

  use constants
  use benchmarks
  use options
  use util_omp_lib

  implicit none

  type (options_t) :: rt_options

  logical :: saturation = .true.
  integer :: ierr
  integer :: nlog

  character (len = *), parameter :: compiler_ = compiler_version()
  character (len = *), parameter :: options_  = compiler_options()

  write (*, '(a)')
  write (*, '(a)')    "Name mini benchmark"
  write (*, '(a,i4)') "omp_get_num_devices()     ", omp_get_num_devices()
  write (*, '(a,i4)') "omp_get_max_threads()     ", omp_get_max_threads()
  write (*, '(a,i4)') "omp_get_initial_device()  ", omp_get_initial_device()
  write (*, '(a,i4)') "omp_get_default_device()  ", omp_get_default_device()
  write (*, '(a)')

  write (*, '(a)')    "Compiler:                 ", compiler_
  write (*, '(a)')    "Compiler options:         ", options_
  write (*, '(a)')

  ! "Parse command line"

  rt_options%aos0       = .true.
  rt_options%soa1       = .true.
  rt_options%soa2       = .true.
  rt_options%gpu        = main_have_target()
  rt_options%ntimestep  = 6*24*7

  write (*, '(a)')       "Benchmark options         "
  write (*, '(a,l)')     "options.aos0              ", rt_options%aos0
  write (*, '(a,l)')     "options.soa1              ", rt_options%soa1
  write (*, '(a,l)')     "options.soa2              ", rt_options%soa2
  write (*, '(a,l)')     "options.gpu               ", rt_options%gpu
  write (*, '(a,i10)')   "options.maxLocal          ", rt_options%maxLocal
  write (*, '(a,i10)')   "options.nSpecies          ", rt_options%nSpecies
  write (*, '(a,i10)')   "options.ntimestep         ", rt_options%ntimestep
  write (*, '(a,f10.2)') "options.dt                ", rt_options%dt

  ! Benchmarks

  if (.false.) then

     print *, ""
     print *, "Validation benchmark for time-dependent source"
     print *, ""
    ! TOP HAT TEST
    rt_options%memory     = MEMORY_ORDER_AOS0
    rt_options%maxLocal   = 1000000
    rt_options%ntimestep  = 6*24*7
    write (*, '(a,i10)') "Number of particles ", rt_options%maxLocal
    call benchmark_time_dependent_source(rt_options, ierr)
    stop "END TOP HAT TEST"
  end if

  block
    ! Validation: Usual range is maxnlog 4, 7, but 4, 6 will is quicker
    ! Satauration: 2^25 is c. 32 million particles
    integer :: maxnval = 6
    integer :: maxsat  = 25
    logical :: ifail = .false.

    if (rt_options%aos0) then
      rt_options%memory = MEMORY_ORDER_AOS0
      ifail = main_validation_benchmark(rt_options, "instantaneous", maxnval)
      ifail = main_validation_benchmark(rt_options, "top-hat", maxnval)
      if (saturation) call main_saturation_benchmark(rt_options, maxsat)
    end if

    if (rt_options%soa1) then
      rt_options%memory = MEMORY_ORDER_SOA1
      ifail = main_validation_benchmark(rt_options, "instantaneous", maxnval)
      if (saturation) call main_saturation_benchmark(rt_options, maxsat)
    end if

    if (rt_options%soa2) then
      rt_options%memory = MEMORY_ORDER_SOA2
      ifail = main_validation_benchmark(rt_options, "instantaneous", maxnval)
      if (saturation) call main_saturation_benchmark(rt_options, maxsat)
    end if
  end block

  print *, ""
  print *, "Benchmarks completed normally"
  print *, ""

contains

  function main_have_target() result(havegpu)

    logical :: havegpu

    havegpu = .false.

    !$omp target teams num_teams(1) map(from: havegpu)
    havegpu = .not. omp_is_initial_device()
    !$omp end target teams

  end function main_have_target

  function main_validation_benchmark(rt_options, id, maxnlog) result(ifail)

    type (options_t),    intent(in) :: rt_options
    character (len = *), intent(in) :: id
    integer,             intent(in) :: maxnlog
    logical                         :: ifail

    ifail = .false.

    print *, ""
    print *, "** Validation benchmark for source: "//id
    print *, ""

    do nlog = 4, maxnlog

       block
         type (options_t) :: opts
         integer          :: ivalid

         opts = rt_options
         opts%maxLocal = 10**nlog

         write (*, '(a)')
         write (*, '(a,i10)') "Number of particles ", opts%maxLocal
         if (id == "instantaneous") then
           call benchmark_instantaneous_source(opts, ivalid)
         end if
         if (id == "top-hat") then
           call benchmark_time_dependent_source(opts, ivalid)
         end if
       end block

    end do

  end function main_validation_benchmark

  !---------------------------------------------------------------------------

  subroutine main_saturation_benchmark(rt_options, nmaxlog)

    ! 2^25 particles is the limit for 16 GB device memory

    type (options_t), intent(in) :: rt_options
    integer,          intent(in) :: nmaxlog

    type (options_t) :: opts
    integer          :: ivalid

    print *, ""
    print *, "** Saturation benchmark for source: instantaneous"
    print *, ""

    do nlog = 10, nmaxlog

      opts = rt_options
      opts%maxLocal = 2**nlog

      write (*, '(a)')
      write (*, '(a,i10)') "Number of particles ", opts%maxLocal
      call benchmark_instantaneous_source(opts, ivalid)

    end do

  end subroutine main_saturation_benchmark

end program main
