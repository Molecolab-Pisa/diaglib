---
title: Notes
---
This page is not intended to be included in the final manual, it gathers some notes for this manual itself.

Important points

- DiagLib is built by default a a shared library. It is to be interfaced with Fortran code by using the module [[dgl_interface]], which exposes all the solvers contained in DiagLib plus the orthogonalization primitives.
- Some tests are built together with the library. A reference can be run to solve the exact problems with lapack procedure (`dgl_reference` executable). This produces a file that will be read and compared with the results from the DiagLib solvers using the `dgl_test` executable

- All solvers accept external procedures that have to be provided by the user. I currently do not know wheter its better like this or to have all procedures passed as pointers.
- Some solvers accept optional procedures, like the metric-vector product in [[davidson_driver]]. In this case, one can not directly pass the procedure. Insted one must pass a **pointer** associated to the desired procedure.

- There is one module that contains the interfaces for the external procedures that have to be passed.
- New procedures for new drivers should be included to grant that correct data types are passed. It is mandatory for optional procedures like metric-vetor products