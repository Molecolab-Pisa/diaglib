---
project: DiagLib
version: 2.0
project_github: https://github.com/Molecolab-Pisa/diaglib
summary: ### A collection of matrix-free iterative eigensolvers
author: I. Giannì, T. Nottoli, R. Alessandro, D. Cianchino, L. Lapi, F. Pes, A. Levitt, F. Lipparini
author_description: MolecoTheory subgroup of MolecoLab and collaborators
email: filippo.lipparini@unipi.it
src_dir: ../src
output_dir: ../doc
page_dir: ./pages
media_dir: ./media
fpp_extensions: f90
macro: DGL_INT_KIND_4
graph: false
sort: type-alpha
max_frontpage_items: 4
warn: true
---

![DiagLib logo](|media|/Diaglib_logo.png)

DiagLib is a collection of matrix-free iterative eigensolvers. It computes a few
eigenvalues and eigenvectors of matrices that are too large to store, given only
routines that apply the matrix, the metric and a preconditioner to a set of vectors.

It is written in Fortran (Fortran 2008: the drivers are submodules of the public
interface) and is built as the shared library `libdiaglib.so`. A second library,
`libdiaglib_c.so`, provides the C interface, which is also used by the python module
`pyDiaglib`.

These pages are generated from the source code. The public Fortran interface is the
module [[dgl_interface]]: it declares the kinds, the error codes, the interfaces of
the user-supplied routines and the interfaces of the four drivers, and it is the only
module that is installed. Everything else is internal, and is documented here for
developers. The Manual has the guides.

The sources are preprocessed with `DGL_INT_KIND_4` defined, so these pages describe the
default build, with 32-bit integers.
