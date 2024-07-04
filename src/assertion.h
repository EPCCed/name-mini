
! A C-like assert() pre-processor macro.

! The conditional here is "ExtraChecks" to be in line with the
! rest of the code (rather than the C standard NDEBUG).

! Note. Stringifying the assertion expression is not very portable
! so have not attempted to include it. The file and line are what
! really matter.

! Clients need to both include this .h file and "use"
! the associated module in assertion.f90

! Assertions must be restricted to a single source code line
! to avoid unexpected compilation problems.
! Do not try to use the continuation marker & to break lines.
! Do not attempt Fortran executable statements in which the absence of
! the macro will cause a compilation failure, e.g.,
! "if (.true.) assert(condition)"

#ifdef NDEBUG
#define assert(expr)
#else
#define assert(expr) call assertion_impl((expr), __FILE__, __LINE__)
#endif
