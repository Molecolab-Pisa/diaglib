---
title: Manual
---

# Welcome to the DiagLib manual

DiagLib solves eigenvalue problems without ever forming the matrix: the caller provides
a routine that applies it to a block of vectors, a preconditioner, and, for generalized
problems, the metric. It is written in Fortran 2008, with C and python interfaces, and
is used in production for linear-response CASSCF, full CI and Hartree-Fock Hessians, at
dimensions of a few million and for as many as 500 roots.

### Features

DiagLib features 4 algorithms:

- Davidson-Liu, for symmetric standard and generalized eigenvalue problems
  ([[dgl_davidson_driver]])
- LOBPCG, for symmetric standard and generalized eigenvalue problems: comparable
  performance with constant memory, as it keeps three blocks of vectors rather than a
  growing subspace ([[dgl_lobpcg_driver]])
- non-symmetric Davidson, for the right and/or left eigenvectors of a non-symmetric
  matrix ([[dgl_davidson_nosym_driver]])
- Swapped Metric Orthogonal Generalized Davidson (SMO-GD), which solves the
  [linear response equations in CASSCF](https://doi.org/10.1021/acs.jpca.5c03618)
  ([[dgl_smogd_driver]])

All of them rely on the same orthogonalization primitives ([[ortho_vs_x]],
[[ortho_cd]]): a Cholesky-based orthonormalization that measures the norm of the
transformation it applies, and a Gram-Schmidt fallback that is used when that
measurement says the vectors are numerically linearly dependent. They are internal
routines rather than part of the public interface, and are documented for developers.

### Where to look

- [Usage](./usage.html): how to call DiagLib from Fortran, C, C++ and python, the
  optional arguments and their defaults, and how errors are reported.
- [Developer guide](./developer_guide.html): how the sources are organized, the
  conventions the drivers follow, and how to build the documentation and run the tests.
