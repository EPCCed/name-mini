module benchmarks

  ! The actual benchmarks.
  !
  ! 1. Instantaneous source at time t = 0
  !    This allocates the number of particles required and does not
  !    require any memory "recycling".
  !
  ! 2. A simple time-dependent source
  !    An example using removeInactive() to provide recycling

  use constants
  use pa
  use flow
  use options
  use rtime
  use source
  use util_omp_lib
  use validate

  implicit none
  private

  public :: benchmark_instantaneous_source
  public :: benchmark_time_dependent_source

contains

  !---------------------------------------------------------------------------

  subroutine benchmark_instantaneous_source(options, ivalid)

    ! Instantaneous source at t = 0, no particle memory recycling

    type (options_t), intent(in)  :: options
    integer,          intent(out) :: ivalid

    real (real64)    :: t0, t1, t2, t3
    integer          :: nparticles

    type (source_t)  :: source
    type (flow_t)    :: flow
    real (std)       :: dt
    integer          :: nt, ntimestep
    integer          :: nRelease

    type (calendar_t)      :: cal
    type (sTime_t)         :: t
    type (sTimeInterval_t) :: deltat

    class (pa_t), pointer  :: ppa => null()

    ivalid = -1 ! No validation by default, 0 = fail, +ve = pass

    dt = options%dt
    ntimestep  = options%ntimestep
    nparticles = options%maxLocal

    write (*, '(a)') ""
    print *, "Real kind:            ", std
    print *, "Timestep:             ", dt
    print *, "Iterations            ", ntimestep
    print *, "Initialise particles: ", nParticles

    t0 = omp_get_wtime()

    ! initialise source
    ! initialise particle array and management
    ! initialise a flow

    source = sourceInstantaneous(options)

    ppa => pa_factory(nParticles, options)

    call flowSubMock4(flow)
    !call flowMockUM6R(flow)

    cal = source%cal
    t   = sTimeZero()
    deltat = sTimeInterval_t(cal, dt)

    print *, ""
    print *, "Start time step loop"
    print *, ""

    ! Fixed flow throughout ...
    !$omp target update to(flow%u)

    t1 = omp_get_wtime()

    do nt = 1, ntimestep

      ! Release from source for [t, t + dt]
      nRelease = sourceReleaseInterval(source, t, deltat)
      if (nRelease > 0) call ppa%releaseFromSource(nRelease, source)

      ! Update dt
      call ppa%updateStep(flow, dt)

      ! Update current time t -> t + dt
      t = sTimeAdd(cal, t, deltat)
    end do

    t2 = omp_get_wtime()

    t3 = omp_get_wtime()

    print *, "Times (sec)"
    print *, "Establish data        ", t1-t0
    print *, "Main loop             ", t2-t1
    print *, "Total                 ", t3-t0
    print *, "Particle it/s (total) ", 1.0*nParticles*ntimestep/(t3-t0)
    print *, "Particle it/s (loops) ", 1.0*nParticles*ntimestep/(t2-t1)
    print *, "Millions it/s (total) ", 0.001*0.001*nParticles*ntimestep/(t3-t0)
    print *, "Millions it/s (loops) ", 0.001*0.001*nParticles*ntimestep/(t2-t1)

    ! Check the total number of active particles in the end of simulation.
    ! This is just the total number of active particles.
    ! This is 4-level flow result...

    block
      integer :: ntotal

      ntotal = ppa%nActive()

      write (*, '(a)') ""
      write (*, '(a,i8)') "Particles still active: ", ntotal
      write (*, '(a)') ""

      select case (nParticles)
      case (10000)
         if (ntotal /= 7959) then
            write (*, '(a)') 'Validation FAIL: expected 7959 active'
            ivalid = 0
         else
            write (*, '(a)') 'Validation PASS: expected 7959 active'
            ivalid = nParticles
         end if
      case (100000)
         if (ntotal /= 79504) then
            write (*, '(a)') 'Validation FAIL: expected 79504 active'
            ivalid = 0
         else
            write (*, '(a)') 'Validation PASS: expected 79504 active'
            ivalid = nParticles
         end if
      case (1000000)
         if (ntotal /= 795821) then
            write (*, '(a)') 'Validation FAIL: expected 795821 active'
            ivalid = 0
         else
            write (*, '(a)') 'Validation PASS: expected 795821 active'
            ivalid = nParticles
         end if
      case (10000000)
         if (ntotal /= 7957394) then
            write (*, '(a)') 'Validation FAIL: expected 7957394 active'
            ivalid = 0
         else
            write (*, '(a)') 'Validation PASS: expected 7957394 actives'
            ivalid = nParticles
         end if
      end select

      ! File output
      if (options%fileOutput) then
#ifdef OLD
        block
        integer :: npart = 10
        real (std), allocatable :: x(:,:)
        character (len = :), allocatable :: filename
        integer :: my_unit
        integer :: ierr

        allocate(x(npart, 3))
        do iP = 1, npart
           x(iP, 1:3) = pa%particles(iP)%x(1:3)
        end do

        filename = validate_filename("instantaneous.dat", nparticles)
        open (newunit = my_unit, file = filename, status = "replace", &
             action = "write", form = "formatted")
        call validate_write(npart, x, my_unit, ierr)
        assert(ierr == 0)
        close (unit = my_unit, status = "keep")
      end block
#else
      stop "NO WRITE IMPLEMENTED"
#endif
   end if

    end block

    ! Report some detail on particle management
    call ppa%printState()

    ! Release resources
    ! Could manage the final deallocation of the pointer more elegantly
    ! - it's coming from the factory method

    call flowRelease(flow)
    call ppa%destroy()
    deallocate(ppa)

  end subroutine benchmark_instantaneous_source

  !---------------------------------------------------------------------------

  subroutine benchmark_time_dependent_source(options, ivalid)

    ! Draft

    type (options_t), intent(in)  :: options
    integer,          intent(out) :: ivalid

    type (source_t)        :: source
    integer                :: nparticles
    integer                :: ntimestep
    real (real64)          :: t0, t1, t2, t3

    type (calendar_t)      :: cal
    type (sTime_t)         :: t
    type (sTimeInterval_t) :: deltat

    class (pa_t), pointer  :: ppa

    type (flow_t)          :: flow

    ! A slightly contrived case with time dependent source which
    ! requires use of pa%removeInactive()

    ! We set up a source with period half of the total run, and
    ! ask it to produce 60% of the total capacity nmaxLocal in
    ! each of two episodes.

    ! At half-time, check the current nActive and run the clean-up.

    ntimestep  = options%ntimestep
    nparticles = options%maxLocal

    ivalid = -1

    ! Gotcha sTimeInterval(real) constructor returns rubiish???
    ! Why not compiler error?

    block
      type (sTimeInterval_t) :: period
      period = sTimeInterval_t(calendar_t(), 0.5*ntimestep*options%dt)
      source = sourceTopHat(options, period, rfraction = 0.6_std)
    end block

    ! Establish particle data

    ppa => pa_factory(nparticles, options)

    ! A flow

    call flowMockUM6R(flow)

    cal = source%cal
    t   = sTimeZero()
    deltat = sTimeInterval_t(cal, options%dt)

    print *, ""
    print *, "Start time step loop"
    print *, ""

    t0 = omp_get_wtime()

    block

      ! Time step loop

      integer :: nt
      integer :: nRelease

      do nt = 1, ntimestep

        nRelease = sourceReleaseInterval(source, t, deltat)
        if (nRelease > 0) call ppa%releaseFromSource(nRelease, source)

        ! Update

        call ppa%updateStep(flow, options%dt)

        if (nt == ntimestep/2) then
          t1 = omp_get_wtime()
          call ppa%removeInactive()
          t2 = omp_get_wtime()
          print *, ""
          print *, "Removed inactive particles at step ", nt
          print *, ""
       end if

        ! Update current time t -> t + dt
        t = sTimeAdd(cal, t, deltat)
      end do

    end block
    t3 = omp_get_wtime()

    ! Validation (serial result with no recycling but same problem)
    ! 10^4 gives  nactive == 5829
    ! 10^5 gives  nacive  == 57024
    ! 10^6 gives          == 569385

    block
      integer :: ntotal

      ntotal = ppa%nActive()

      write (*, '(a)') ""
      write (*, '(a,i8)') "Particles still active: ", ntotal
      write (*, '(a)') ""

      select case (nParticles)
      case (10000)
         if (ntotal /= 5829) then
            write (*, '(a)') 'Validation FAIL: expected 5829 active'
            ivalid = 0
         else
            write (*, '(a)') 'Validation PASS: expected 5829 active'
            ivalid = nParticles
         end if
      case (100000)
         if (ntotal /= 57024) then
            write (*, '(a)') 'Validation FAIL: expected 57024 active'
            ivalid = 0
         else
            write (*, '(a)') 'Validation PASS: expected 57024 active'
            ivalid = nParticles
         end if
      case (1000000)
         if (ntotal /= 569385) then
            write (*, '(a)') 'Validation FAIL: expected 569385 active'
            ivalid = 0
         else
            write (*, '(a)') 'Validation PASS: expected 569385 active'
            ivalid = nParticles
         end if
      end select
    end block

    ! Final details.
    call ppa%printState()

    print *, ""
    print *, "Times (sec)"
    print *, "Total per loop iter.  ", (t3 - t0)/ntimestep
    print *, "Time per sort         ", (t2 - t1)

    ! Clean up
    call flowRelease(flow)
    call ppa%destroy()
    deallocate(ppa)

  end subroutine benchmark_time_dependent_source

end module benchmarks
