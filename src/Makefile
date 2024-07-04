
# The default position is to build for "omp", that is build a
# single executable which can run on either cpu or gpu
# depending on run time choice.

default:
	$(MAKE) --no-print-directory omp

include ./current.mk

OBJS    = assertion_stop.o assertion.o util_omp_lib.o constants.o maths.o \
          options.o validate.o rtime.o source.o flow.o \
          pa.o \
	  storage_aos0.o storage_soa1.o storage_soa2.o \
          pa_aos0.o pa_soa1.o pa_soa2.o \
	  pa_aos0_cpu.o pa_aos0_gpu.o \
	  pa_soa1_cpu.o pa_soa1_gpu.o \
          pa_soa2_cpu.o pa_soa2_gpu.o \
          pa_submodule.o benchmarks.o

both:
	make --no-print-directory cpu
	make --no-print-directory omp

cpu:
	$(MAKE) "ARCH_OPTS=$(FC_OPTS_CPU)" cpu-t

omp:
	$(MAKE) "ARCH_OPTS=$(FC_OPTS_GPU)" omp-t

cpu-t:	$(OBJS)
	$(FC) $(OBJS) $(LD_OPTS_CPU) -o mini-cpu main.f90

omp-t:	$(OBJS)
	$(FC) $(OBJS) $(FC_OPTS) $(LD_OPTS_GPU) -o mini-omp main.f90


.SUFFIXES:
.SUFFIXES: .f90 .o

.f90.o:	
	$(FC) $(FC_OPTS) $(ARCH_OPTS) -c $?

.phony:	clean
clean:
	$(RM) -f *.i *.mod *.smod *.o mini-cpu mini-gpu
