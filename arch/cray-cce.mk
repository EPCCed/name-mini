
# Cray CCE
# Last cce/16.0.1

# E.g.,
# For host        module load craype-accel-host
# For device      module load craype-accel-amd-gfx90a

FC           = ftn -g
FC_OPTS      = -e Z -O3 -hipa1 -M 878
FC_OPTS_CPU  = -fopenmp
FC_OPTS_GPU  = -fopenmp
LD_OPTS_GPU  = -fopenmp
LD_OPTS_GPU  = -fopenmp

# nb.

# Fortran warning 878 is about indirect inclusion of modules,
# which seems impossible to avoid, and seems OOT. So suppress with "-M 878".

# For serial compilation, one may wish to suppress 1077 concerning
# ignored directives.
