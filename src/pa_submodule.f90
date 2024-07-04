submodule (pa) pa_submodule

  use pa_aos0_cpu
  use pa_aos0_gpu
  use pa_soa1_cpu
  use pa_soa1_gpu
  use pa_soa2_cpu
  use pa_soa2_gpu

  implicit none

contains

  module function pa_factory(nParticles, options) result(pa)

    integer,          intent(in) :: nParticles
    type (options_t), intent(in) :: options

    class (pa_t),     pointer    :: pa

    pa => null()

    select case (options%memory)
    case (MEMORY_ORDER_AOS0)
      if (      options%gpu) pa => pa_aos0_gpu_create_pointer(nParticles)
      if (.not. options%gpu) pa => pa_aos0_cpu_create_pointer(nParticles)
    case (MEMORY_ORDER_SOA1)
      if (      options%gpu) pa => pa_soa1_gpu_create_pointer(nParticles)
      if (.not. options%gpu) pa => pa_soa1_cpu_create_pointer(nParticles)
    case (MEMORY_ORDER_SOA2)
      if (      options%gpu) pa => pa_soa2_gpu_create_pointer(nParticles)
      if (.not. options%gpu) pa => pa_soa2_cpu_create_pointer(nParticles)
    case default
      print *, "pa_factory: invalid memory option ", options%memory
    end select

    if (.not. associated(pa)) stop "Factory method failed "

  end function pa_factory

end submodule pa_submodule
