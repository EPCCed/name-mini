# name-mini

The code here is a simple proxy application for NAME, the
Met Office dispersion model. The code is written in Fortran
and uses OpenMP offload directives to admit execution on
both CPU and GPU back ends.

## Build

Sets of compiler options to build the code are available in
the `arch` directory. One must be made available for the
`Makefile` in the `src` directory. E.g.,
```
$ cd src
$ ln -s ../arch/nvidia-nvfortran.mk current.mk
$ make
```
The will produce an executable `mini-omp`.

The executable performs a self-validation which will report `PASS`
or `FAIL` at run time.

## Portability

The following limited set of platforms/compilers have been investigated ("ok" indicates successful compilation and execution):

| Compiler/Platform | CPU | GPU              | Notes                  |
|-------------------|-----|------------------|------------------------|
| Cray/AMD          | ok  | does not compile | Cray CCE 16.0          |
| `ifx`/Intel       | ok  | N/A              | `ifx` 2024.1           |
| nvfortran/NVIDIA  | ok  | ok               | SDK 2024.5 (Intel CPU) |


