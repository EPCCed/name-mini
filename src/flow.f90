module flow

  ! A fake flow object to provide a 1-d vertical profile

  ! Note: the flow object is not part of any composition, so device
  ! data mapping always takes place here. Such mappping must take
  ! place in a subroutine (rather than a function) to retain the
  ! device references.

  use constants

  implicit none
  private

  type, public :: flow_t
    integer                        :: nz         ! A vertical profile only
    real (std)                     :: k(3)       ! Eddy diffusivity.
    real (std), allocatable        :: u(:,:)     ! Mean velocity.
    real (std), allocatable        :: zb(:)      ! boundaries zb
  end type flow_t

  public :: flowSubMock4
  public :: flowMockUM6R
  public :: flowRelease
  public :: flowZToLevel

contains

  !---------------------------------------------------------------------------
  !---------------------------------------------------------------------------

  subroutine flowSubMock4(flow)

    ! A flow with 4 levels

    type (flow_t), intent(out) :: flow

    integer :: k

    flow%nz = 4

    allocate(flow%u(3, flow%nz))
    allocate(flow%zb(flow%nz + 1))

    flow%k(1) = 20.0_std
    flow%k(2) = 20.0_std
    flow%k(3) = 1.00_std

    do k = 1 , flow%nz
      flow%u(1, k) = 2.5*k
      flow%u(2, k) = 2.5*k
      flow%u(3, k) = 0.01
    end do

    flow%zb(1) = 0.00000_std
    flow%zb(2) = 1000.00_std
    flow%zb(3) = 2000.00_std
    flow%zb(4) = 5000.00_std
    flow%zb(5) = 10000.0_std

    !$omp target enter data map(always, to: flow)
    !$omp target enter data map(always, to: flow%u)
    !$omp target enter data map(always, to: flow%zb)

  end subroutine flowSubMock4

  !---------------------------------------------------------------------------
  !---------------------------------------------------------------------------

  subroutine flowMockUM6R(flow)

    type (flow_t), intent(out) :: flow

    ! These are the w level from MetDefnUM6R using nz = 32
    ! and the levels are height above sea level

    integer, parameter                   :: nZ = 32
    real (std), dimension(nZ), parameter :: z = &
         (/  0.0000, 20.0018, 80.0011, 179.998, 319.998, &
            500.002, 719.997, 980.002, 1280.00, 1620.00, &
            2000.00, 2420.00, 2880.00, 3380.00, 3920.00, &
            4500.00, 5120.00, 5780.00, 6480.00, 7220.00, &
            8000.00, 8820.00, 9680.00, 10580.0, 11520.0, &
            12500.0, 13520.0, 14580.0, 15680.0, 16820.0, &
            18000.0, 19180.0 /)

    integer :: k

    flow%nZ = nZ - 1

    allocate(flow%u(3, flow%nZ))
    allocate(flow%zb(flow%nZ + 1))

    do k = 1, nZ
       flow%zb(k) = z(k)
    end do

    ! Wind profile is a fiction: just really for the boundary layer.
    ! The vertical velocity is just "up"

    block
      real (std) :: z              ! mid-level height
      real (std) :: uref = 20.0    ! a reference speed
      real (std) :: zref = 100.0   ! a reference height
      real (std) :: alpha = 0.143  ! empircal constant

      do k = 1, flow%nZ
        z = 0.5*(flow%zb(k) + flow%zb(k+1))
        flow%u(1, k) = uref*(z/zref)**alpha
        flow%u(2, k) = 2.0
        flow%u(3, k) = 0.01
      end do
    end block

    ! Again, just some indicative numbers ...

    flow%k(1) = 20.0
    flow%k(2) = 20.0
    flow%k(3) = 1.0

    !$omp target enter data map(always, to: flow)
    !$omp target enter data map(always, to: flow%u)
    !$omp target enter data map(always, to: flow%zb)

  end subroutine flowMockUM6R

  !---------------------------------------------------------------------------
  !---------------------------------------------------------------------------

  subroutine flowRelease(flow)

    type (flow_t), intent(inout) :: flow

    !$omp target exit data map(delete: flow%zb)
    !$omp target exit data map(delete: flow%u)
    !$omp target exit data map(delete: flow)

    deallocate(flow%zb)
    deallocate(flow%u)

    flow%nz = 0
    flow%k(:) = 0.0

  end subroutine flowRelease

  !---------------------------------------------------------------------------
  !---------------------------------------------------------------------------

  function flowZToLevel(flow, z) result(k)

    ! Return integer level k containing vertical position z

    type (flow_t), intent(in) :: flow
    real (std),    intent(in) :: z
    integer                   :: k

    !$omp declare target

    integer :: iZ, kmin, kmax

    kmax = flow%nZ
    kmin = 1
    k = 0
    do
      if (kmin > kmax) exit

      k = kmin + (kmax - kmin)/2

      if (flow%zb(k) > z) then
        kmax = k-1
      else if (flow%zb(k+1) <= z) then
        kmin = k+1
      else
        exit
      end if
    end do

  end function flowZToLevel

end module flow
