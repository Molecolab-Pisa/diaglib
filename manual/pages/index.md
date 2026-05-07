---
title: Manual
---

# Welcome to DiagLib Manual 1.0

DiagLib is a collection of matrix-free iterative eigensolvers. It is wrtitten in Fortran (compliant with Fortran2003 standard). It is compiled a shared library by the name `libdiaglib.so`. A secondary library, necessary to use DiagLib from C, is also produced (`libdiaglib_c.so`). Finally also a python interface (`pyDiaglib`) is provided and installed in an on-the-fly virtual enviroment as a `pip` package.

### Features
The current release of DiagLib features 4 algorithms:

- Davidson algorithm for the solution Symmetric, Standard and Generalized, eigenvalue problems
- LOBPCG algorithm for the solution Symmetric, Standard and Generalized, eigenvalue problems (similarly performant but memory limited in comparison to Davidson)
- Davidson algorithm for the solution of Non-Symmetrix eigenvalue problems
- Swapped Metric Orthogonal Generalized Davidson (SMOGD) algorithm (solves the [Linear Respons equations in CASSCF](https://doi.org/10.1021/acs.jpca.5c03618))

These algorithms make use of very stable primitives for the orthogonalization and orthonormalization of set of vectors ([[ortho_vs_x]],[[ortho_cd]]). These are also provided as part of DiagLib.
