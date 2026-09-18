---
title: Usage
---

# Usage

Every driver needs the same three things: the sizes of the problem, the routines that
apply the matrix and the preconditioner, and arrays that carry a guess in and the
eigenpairs out. Everything else is optional and has a default.

## The user-supplied routines

The interfaces are fixed, and are declared in [[dgl_interface]] as [[dgl_matvec]],
[[dgl_precnd]] and [[dgl_smogd_precnd]]:

```fortran
subroutine matvec(n, m, x, y)   ! y = A x, for m vectors at a time
subroutine precnd(n, m, shift, x, y)
```

`n` and `m` are `integer(dgl_int)`, `x(n, m)` and `y(n, m)` are `real(dgl_real)`. The
same interface is used for the metric-vector product. `shift` is minus the lowest
non-converged approximate eigenvalue, so that shift-and-invert preconditioners such as
\( (D - \lambda)^{-1} \) can be built, or zero: the choice is the `dgl_precnd_shift`
argument, true by default for the Davidson drivers and false for LOBPCG.

The routines are called back from inside the driver, once per block of vectors. They
may themselves call a DiagLib driver: DiagLib keeps no global state, so calls nest and
several threads may run their own calls at the same time.

## From Fortran

```fortran
use dgl_interface
integer(dgl_int) :: info
real(dgl_real) :: eig(n_max), evec(n, n_max)
logical :: ok

evec = 0.0_dgl_real             ! or a guess, which is completed if incomplete
call dgl_davidson_driver(n, n_targ, n_max, matvec, precnd, eig, evec, ok, &
                         dgl_tol=1.0e-8_dgl_real, dgl_info=info)
```

`dgl_interface.mod` is the only module file installed, and has to come from the same
compiler and version that built the library. With the installed CMake package:

```cmake
find_package(diaglib REQUIRED)
target_link_libraries(my_program PRIVATE diaglib::diaglib)
```

A generalized problem is selected by passing the optional `metvec` argument, which is a
**procedure pointer** rather than a procedure:

```fortran
procedure(dgl_matvec), pointer :: metvec_p => null()
metvec_p => my_metric_product
call dgl_davidson_driver(n, n_targ, n_max, matvec, precnd, eig, evec, ok, &
                         metvec=metvec_p, dgl_info=info)
```

## Optional arguments

| Argument | Meaning | Default |
| --- | --- | --- |
| `dgl_tol` | convergence threshold on the residual norms | \(10^{-7}\) |
| `dgl_max_iter` | maximum number of iterations | 100 |
| `dgl_dav_iter` | iterations before a Davidson restart | 25 |
| `dgl_memory`, `dgl_memory_unit` | memory DiagLib may allocate | 80, `MB` |
| `dgl_precnd_shift` | pass the shift to `precnd` | `.true.`, except LOBPCG |
| `dgl_shift` | constant added to the eigenvalues **only when printed** | 0 |
| `dgl_verbose` | print a convergence table at every iteration | `.false.` |
| `metvec` | metric-vector product, selects the generalized problem | absent |
| `dgl_info` | error status | see below |

A root is converged when the RMS norm of its residual is below `dgl_tol` and its largest
component is below ten times that.

`dgl_shift` is a printing convenience, for total energies when the matrix does not
include a constant contribution. It changes neither the iterations nor the returned
eigenvalues.

## Errors and non-convergence

They are different things. Not converging within `dgl_max_iter` iterations is **not** an
error: the driver returns `ok = .false.` together with the current approximations, which
are usually a good guess to restart from. Errors set `dgl_info` to one of the negative
codes declared in [[dgl_interface]]:

| Code | Meaning |
| --- | --- |
| `dgl_success` | no error |
| `dgl_err_input` | invalid input, or a user routine returned NaN or Inf |
| `dgl_err_memory` | allocation failure, or the memory limit was exceeded |
| `dgl_err_lapack` | a LAPACK routine failed |
| `dgl_err_ortho` | an orthogonalization failed, both the default route and the fallback |
| `dgl_err_mismatch` | the left and right runs of the non-symmetric driver disagree |

If `dgl_info` is not passed, DiagLib prints the error and stops the program, so an
unchecked failure cannot pass unnoticed.

## From C and C++

```c
#include "diaglib.h"

dgl_int info;
bool ok;
dgl_davidson_driver(n, n_targ, n_max, matvec, precnd, NULL, eig, evec, &ok, &info,
                    false, 1e-8, 100, 25, 0.0, true, 80, "MB");
```

`dgl_int` follows the integer kind of the build, and `dgl_integer_kind()` returns the
size in bytes of the integers of the library that was actually linked, so a mismatch
between header and library can be detected at run time. Link `diaglib::diaglib_c`. From
C++, the callbacks must be declared `extern "C"`, and an exception must not be allowed
to propagate through the Fortran frames.

## From python

```python
import numpy as np
import pyDiaglib as dgl

calc = dgl.diaglib(lib_path, n, n_targ, n_max, tol=1e-8)
eig = np.zeros(n_max)
evec = np.zeros((n, n_max), order="F")     # or a guess
ok = calc.dgl_davidson_driver(eig, evec, matvec, precnd)   # eig, evec filled in place
```

One `diaglib` object holds the options and can run any driver: `dgl_davidson_driver`,
`dgl_lobpcg_driver`, `dgl_davidson_nosym_driver` and `dgl_smogd_driver`, the first two
taking an optional `metvec` for the generalized problem. `calc.integer_kind` reports the
integer size of the library that was loaded.

The shared library is located from the argument, from the `DIAGLIB_C_LIBRARY`
environment variable, or through the loader. Errors are raised as `RuntimeError`. The
module can be installed with `pip install src/python_interface`, or through the CMake
option `-DBUILD_PYTHON=ON`.

## Examples

`test/` contains working programs for every interface, and `test/test_consumers` is a
small separate CMake project that uses an **installed** DiagLib from Fortran, C, C++ and
python, which is the shortest path to a working template.
