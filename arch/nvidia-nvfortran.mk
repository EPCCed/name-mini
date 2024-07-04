# nvidia nvfortran
# last SDK 24.5

FC           = nvfortran
FC_OPTS      = -Mpreprocess -g -Ofast -DNDEBUG
FC_OPTS_CPU  = -mp
FC_OPTS_GPU  = -mp=gpu #-gpu=ptxinfo #-gpu=sm_70
LD_OPTS_CPU  = -mp
LD_OPTS_GPU  = -mp=gpu

# Comments:
# May need to maintain -O2 (not -O3) to keep binary reproduceability
# and that may be necssary but not sufficient

# e.g., -gpu=ptxinfo for ptx assembler information

# e.g., -gpu=sm_70 for specfic sm arch but can change result in terms of
#       number of registers etc cf. the same without it (bug?). Avoid.

# e.g., -gpu=maxregcount:32 but probably degrades performance significantly

# e.g., Do not try -Mipa=inline ("deprecated")

# e.g., can try -Minline but creates additional compilation problems without
#       apparent benefit

# For -G use "-gpu=debug"
