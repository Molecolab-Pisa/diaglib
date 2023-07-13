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
diaglib is licensed under the LGPL v2.1 license

## Description
diaglib - a fortran library of matrix-free iterative algorithm to
compute a few eigenvalues and eigenvectors of large matrices.

diaglib provides an implementation of two matrix-free algorithms to
compute a few eigenvalues and eigenvectors of a large, possibly sparse
matrix.

the available algorithms are

1) locally optimal block preconditioned conjugate gradient 

2) davidson-liu

both algorithms require two user-provided routines to apply the matrx
and a suitable preconditioner to a set of vectors.
such routines have the following interface:

  subroutine matvec(n,m,x,ax)  
  subroutine precnd(n,m,shift,x,ax)

where n,m are integers and x(n,m) and ax(n,m) are double precision
arrays.
as using the first eigenvalue in a shift-and-invert spirit is very 
common, a double precision scalar shift is also passed to precnd.

both implementations favor numerical stability over efficiency and are
targeted at applications in molecular quantum chemistry, such as in
(full) ci or augmented hessian calculations, where typically m << n.

## Dependencies
BLAS and LAPACK

## Example use
A simple-minded driver is provided, so that the user can compile the
library and test it on toy matrices.

