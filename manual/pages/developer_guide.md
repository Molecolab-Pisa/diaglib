---
title: Developer Guide
---

# Developer Guide

## How the sources are organized

`src/fortran` holds the library. The public interface is the module [[dgl_interface]],
which declares the kinds, the error codes, the interfaces of the user-supplied routines
and the interfaces of the four drivers, with their documentation. Each driver is a
**submodule** of it, so that `dgl_interface.mod` is the only module file installed and
everything else stays internal:

| File | Contents |
| --- | --- |
| `dgl_interface.f90` | public module: kinds, error codes, abstract interfaces, driver interfaces |
| `davidson_driver.f90` | submodule `dgl_davidson` |
| `lobpcg_driver.f90` | submodule `dgl_lobpcg` |
| `davidson_nosym_driver.f90` | submodule `dgl_davidson_nosym` |
| `smogd_driver.f90` | submodule `dgl_smogd` |
| `dgl_orthogonalizations.f90` | orthogonalization primitives |
| `dgl_global_utils.f90` | context, memory, errors, checks, timings |
| `dgl_minor_utils.f90` | guess checking, printing |
| `dgl_lapack.f90` | explicit interfaces for the BLAS and LAPACK routines used |

`src/c_interface` holds one `bind(C)` wrapper module per driver plus `diaglib.h`, and
`src/python_interface` the `ctypes` module. `test/` holds the test suites, and
`test/test_consumers` a separate project that consumes an installed DiagLib.

## Conventions the drivers follow

- **No global state.** Everything a call needs lives in a `dgl_context`, allocated on
  entry and passed down. That is what makes the drivers re-entrant and thread-safe, and
  it is easy to break: a local variable with an initializer implies `save` in Fortran,
  so it must not be used for anything that belongs to a call.
- **Repeated declarations.** The submodule procedures restate their arguments instead of
  using `module procedure`. This is deliberate: some gfortran versions lose the
  interfaces of dummy procedures inside a `module procedure` body, and then no longer
  check the calls to the user-supplied routines.
- **Errors do not stop the program.** A routine that detects an error records it in
  `ctx%status` and returns; the driver cleans up and reports it through `dgl_info`. Only
  a missing `dgl_info` turns an error into a stop.
- **Integers are `integer(ip)` everywhere**, including every argument passed to BLAS and
  LAPACK, which the explicit interfaces in `dgl_lapack.f90` enforce at compile time.
- **Optional procedure arguments are pointers.** A procedure cannot be passed directly
  to an optional dummy procedure, so the metric-vector product is a
  `procedure(dgl_matvec), pointer`.

## The orthogonalization primitives

[[orthonormalize_vs_x]] is the routine the drivers use: it orthogonalizes a block of new
vectors against the current space and among themselves. It takes the fast route first,
the Cholesky-based [[ortho_cd]], and measures the norm of the inverse factor it applies.
That norm is the growth factor of the transformation, and it is the signal that the
vectors are numerically linearly dependent. [[gs_fallback]], a per-vector Gram-Schmidt
repeated twice, is used only when Cholesky fails, when the growth factor exceeds
\(10^{10}\), or when the default route has not converged in three iterations. Changing
any of this is exactly the kind of change the stress tests exist for.

## Adding a driver

1. Declare its interface in [[dgl_interface]], with the documentation of every argument,
   and an abstract interface for any new kind of user-supplied routine.
2. Implement it as a new submodule, restating the arguments.
3. Allocate through the `mallocate` interface so the memory is accounted for, and check
   the results of user routines with [[dgl_check_finite]] before they reach LAPACK.
4. Add it to the Fortran, C and python test suites, and to the stress tests.

## Tests and the pipeline

`ctest` runs the Fortran, C, C++ and python suites. The Fortran suite compares against a
reference file produced by the `dgl_reference` executable, which diagonalizes the same
problems with LAPACK, so that executable has to be run before `dgl_test` — `ctest`
handles the order. The stress tests, enabled with `-DDGL_STRESS_TESTS=ON` and selected
with `ctest -L stress`, run every driver many times from random guesses against dense
references.

The CI pipeline builds and tests every push with gfortran and OpenBLAS, in both integer
kinds, with warnings as errors, plus the consumer project and the stress tests; the jobs
that use MKL and the Intel compilers run on `main`, on tags, on demand and weekly.

## Building this documentation

The pages are generated from the `!!` comments in the sources by
[FORD](https://forddocs.readthedocs.io), configured in `manual/diaglib.md`:

```console
pip install ford
ford manual/diaglib.md
```

The result lands in `doc/`, which is not tracked. Two settings in the project file are
worth knowing about. `warn: true` makes FORD list undocumented entities, which is a
backlog of about a thousand internal variables rather than a failure, so the `docs` job
of the pipeline fails only on a `[[link]]` that no longer resolves. And `fpp_extensions:
f90` together with `macro: DGL_INT_KIND_4` is needed because FORD does not preprocess
lower-case `.f90` files by default: without it, both branches of the `#ifdef` that picks
the integer kind end up in the documentation, and `dgl_int` is documented twice.

The project file is a markdown file whose front matter holds the settings. Note that a
`#` comment inside that front matter silently discards the keys that follow it.
