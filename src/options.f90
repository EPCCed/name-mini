module options

  ! Options for the benchmark itself

  use constants
  implicit none
  private

  ! Memory model
  integer, parameter, public :: MEMORY_ORDER_AOS0 = 1
  integer, parameter, public :: MEMORY_ORDER_SOA1 = 2
  integer, parameter, public :: MEMORY_ORDER_SOA2 = 3

  ! start / stop time (whole run)
  ! initialise at start or use source? (related to source start/sop) 

  type, public :: options_t
    logical               :: aos0 = .true.         ! run aos0
    logical               :: soa1 = .true.         ! run soa1
    logical               :: soa2 = .false.        ! run soa2
    logical               :: gpu = .false.         ! Run on GPU
    logical               :: fileOutput = .false.  ! sample output to file
    integer               :: maxLocal              ! particles
    integer               :: nSpecies = 1          ! masses
    integer               :: ntimestep             ! 6*24*7
    real (std)            :: dt = 600.0            ! time step (s)
    integer               :: memory = -1           ! current memory model
  end type options_t
  
end module options
