#ifndef DIAGLIB_H
#define DIAGLIB_H

#include <stdbool.h>

void davidson_driver_c(
#ifdef DGL_INT_KIND_4
    int n,
    int n_targ,
    int n_max,
#elif DGL_INT_KIND_8
    long int n,
    long int n_targ,
    long int n_max,
#endif
    void (*matvec)(int*, int*, double*, double*),
    void (*precnd)(int*, int*, double*, double*, double*),
    void (*metvec)(int*, int*, double*, double*),
    double* eig,
    double* evec,
    bool* ok,
    bool verbose,
    double tol,
#ifdef DGL_INT_KIND_4
    int max_iter,
    int dav_iter,
#elif DGL_INT_KIND_8
    int long max_iter,
    int long dav_iter,
#endif
    double shift,
#ifdef DGL_INT_KIND_4
    int memory,
#elif DGL_INT_KIND_8
    int long memory,
#endif
    const char* memory_unit
);

void nonsym_driver_c(
    bool verbose,
    int n,
    int n_targ,
    int n_max,
    int max_iter,
    double tol,
    int max_dav,
    double shift,
    void (*matvec_r)(int*, int*, double*, double*),
    void (*matvec_l)(int*, int*, double*, double*),
    void (*precnd)(int*, int*, double*, double*, double*),
    double* eig,
    double* evec_r,
    double* evec_l,
    int side,
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

void smogd_driver_c(
    bool verbose,
    int n,
    int n2,
    int n_targ,
    int n_max,
    int max_iter,
    double tol,
    int maxdav,
    void (*apbmul)(int*, int*, double*, double*),
    void (*ambmul)(int*, int*, double*, double*),
    void (*spdmul)(int*, int*, double*, double*),
    void (*smdmul)(int*, int*, double*, double*),
    void (*lrprec)(int*, int*, double*, double*, double*, double*, double*),
    double* eig,
    double* evec,
    bool* ok
);

#endif
