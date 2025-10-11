#ifndef DIAGLIB_H
#define DIAGLIB_H

#include <stdbool.h>

void davidson_driver_c(
    bool verbose,
    int n,
    int n_targ,
    int n_max,
    int max_iter,
    int max_dav,
    double tol,
    double shift,
    void (*matvec)(int*, int*, double*, double*),
    void (*precnd)(int*, int*, double*, double*, double*),
    double* eig,
    double* evec,
    bool* ok
);

void lobpcg_driver_c(
    bool verbose,
    bool gen_eig,
    int n,
    int n_targ,
    int n_max,
    int max_iter,
    double tol,
    double shift,
    void (*matvec)(int*, int*, double*, double*),
    void (*precnd)(int*, int*, double*, double*, double*),
    void (*bvec)(int*, int*, double*, double*),
    double* eig,
    double* evec,
    bool* ok
);


#endif
