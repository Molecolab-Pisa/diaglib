# diaglib

## Authors and acknowledgment  
Ivan Gianni', Tommaso Nottoli, Riccardo Alessandro, Federica Pes, and Filippo Lipparini  
MoLECoLab Pisa  
Department of Chemistry and Industrial Chemistry  
University of Pisa  
Via G. Moruzzi 13, I-56124, Pisa, Italy  
Pisa, november 2022

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7680658.svg)](https://doi.org/10.5281/zenodo.7680658)

## License
diaglib is licensed under the Mozilla Public license 2.0

## Description
diaglib - a fortran library of matrix-free iterative algorithms to
compute a few eigenvalues and eigenvectors of large matrices, with
C and python interfaces.

the available drivers are

1) davidson-liu, for symmetric standard and generalized problems
   (`dgl_davidson_driver`)

2) locally optimal block preconditioned conjugate gradient, for symmetric
   standard and generalized problems (`dgl_lobpcg_driver`)

3) non-symmetric davidson, for the right and/or left eigenvectors of
   non-symmetric matrices (`dgl_davidson_nosym_driver`)

4) swapped metric-orthogonal generalized davidson (smo-gd), for linear
   response problems (`dgl_smogd_driver`)

all the drivers require user-provided routines to apply the matrix
(and, if needed, the metric) and a suitable preconditioner to a set of
vectors. such routines have the following interface:

  subroutine matvec(n,m,x,ax)
  subroutine precnd(n,m,shift,x,ax)

where n,m are integer(dgl_int) and x(n,m) and ax(n,m) are real(dgl_real)
arrays. a real(dgl_real) scalar shift is also passed to precnd: minus the lowest
non-converged approximate eigenvalue, to build shift-and-invert preconditioners
such as (D - lambda)^-1, as in davidson's method, or zero. the choice is made
with the optional dgl_precnd_shift argument (precnd_shift in C and python):
the default is to shift for the davidson drivers, and not to shift for lobpcg,
which works best with a positive definite, well conditioned preconditioner
(e.g., an approximation of the inverse of the matrix).

the optional shift argument of the drivers (not available for smo-gd) is only
added to the eigenvalues when they are printed, e.g., to print total energies
when the matrix does not include a constant energy contribution: it does not
change the computation nor the returned eigenvalues.

diaglib keeps no global state: a driver can be called from inside the
user-provided routines of another driver call, and different threads can call
the drivers at the same time (the C and python interfaces need OpenMP for this),
provided that the user-provided routines allow it.

all implementations favor numerical stability over efficiency and are
targeted at applications in molecular quantum chemistry, such as in
(full) ci or augmented hessian calculations, where typically m << n.

## Dependencies
- a Fortran and a C compiler, CMake (3.15 or newer; 3.22 or newer to select
  a BLAS/LAPACK library with 64-bit integers automatically)
- BLAS and LAPACK
- optionally, OpenMP (used for wall-clock timings, and to call the C interface
  from different threads at the same time)
- for the python interface, python 3.9 or newer with numpy

## Building and installing
    cmake -S . -B build [options]
    cmake --build build
    ctest --test-dir build        # or: cd build && ctest
    cmake --install build --prefix <prefix>

main options:
- `-DDGL_INTEGER_KIND=8`: use 64-bit integers (default: 4, i.e. 32-bit). A
  BLAS/LAPACK library with 64-bit integers is required: on Debian and Ubuntu, it
  is provided by `libopenblas64-dev`, while a 32-bit library will be rejected by
  CMake with `Could NOT find BLAS`.
- `-DBLA_VENDOR=<vendor>`: choose the BLAS/LAPACK library, e.g. `OpenBLAS`,
  `Intel10_64lp` (MKL, 32-bit integers), `Intel10_64ilp` (MKL, 64-bit integers).
  By default, the first library found by CMake is used.
- `-DBUILD_PYTHON=ON`: install the python interface (`pyDiaglib.py`).
- `-DBUILD_TESTING=OFF`: do not build the tests.
- `-DDGL_STRESS_TESTS=ON`: also build the stress tests (`ctest -L stress`), which run every
  driver many times from random guesses and compare the results with dense LAPACK.
- `-DDGL_STRICT_WARNINGS=ON`: compile the libraries with extra warnings, treated as errors.
- `-DDGL_NATIVE_ARCH=ON`: optimize for the build machine (not portable).
- `-DCMAKE_BUILD_TYPE=Debug`: build with run-time checks (default: `Release`).

## Usage
- Fortran: `use dgl_interface`, which provides the drivers, the kinds
  `dgl_int` and `dgl_real` and the error codes. `dgl_interface.mod` is the only
  module file installed: it has to be used with the same compiler (and version)
  that built the library. The installed CMake package can be
  used with `find_package(diaglib)` and the target `diaglib::diaglib`.
- C and C++: `#include "diaglib.h"`, which provides the drivers, the integer
  type `dgl_int` and the error codes, and link `libdiaglib_c` (CMake target
  `diaglib::diaglib_c`).
- python: `import pyDiaglib`, see the documentation in `pyDiaglib.py`. The
  path of `libdiaglib_c` can be passed explicitly or set in the
  `DIAGLIB_C_LIBRARY` environment variable. Alternatively, the module can be
  installed with `pip install src/python_interface`.

errors (invalid input, not enough memory, failures of LAPACK or of the
orthogonalizations) are reported through the optional `dgl_info` argument of the
Fortran drivers (if it is not present, the program is stopped), through the
`info` argument of the C functions, and as `RuntimeError` in python.
not converging within the maximum number of iterations is not an error: the
drivers then return `ok = .false.` and the latest approximations.

the tests in `test/` show how to use the drivers from Fortran, C and python.
`test/test_consumers` is a separate CMake project that uses an installed DiagLib the way a user
would (`find_package(diaglib)` from Fortran, C and C++, and the installed python module):

    cmake -S test/test_consumers -B build_consumers -DCMAKE_PREFIX_PATH=<prefix>
    cmake --build build_consumers
    ctest --test-dir build_consumers

all of this is run by the CI pipeline. Every push is built with gfortran and OpenBLAS, with
both integer kinds, and runs the tests, the strict warnings, the consumers, the stress
tests and the documentation. The jobs that use MKL and the Intel compilers download about
1.3 GB from the Intel oneAPI repository, and therefore only run on `main`, on tags, on
demand, and once a week.
The runner image is pinned, so that a result only changes when the repository does, with one
job following the newest image and its default compiler to warn when the two diverge.

## Documentation
the reference documentation is generated from the `!!` comments in the sources with
[FORD](https://forddocs.readthedocs.io), configured in `manual/diaglib.md`, together with
the guides in `manual/pages`:

    pip install ford
    ford manual/diaglib.md

the result is written to `doc/` (not tracked) and starts at `doc/index.html`. The CI
pipeline builds it on every push and fails if a `[[link]]` in it no longer resolves.
