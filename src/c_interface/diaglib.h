#ifndef DIAGLIB_H
#define DIAGLIB_H

#include <stdbool.h>

extern void dgl_davidson_driver_c(
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

    
void dgl_lobpcg_driver_c(
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
#elif DGL_INT_KIND_8
    int long max_iter,
#endif
    double shift,
#ifdef DGL_INT_KIND_4
    int memory,
#elif DGL_INT_KIND_8
    int long memory,
#endif
    const char* memory_unit
);

void dgl_davidson_nosym_driver_c(
#ifdef DGL_INT_KIND_4
    int n,
    int n_targ,
    int n_max,
#elif DGL_INT_KIND_8
    long int n,
    long int n_targ,
    long int n_max,
#endif
    void (*matvec_r)(int*, int*, double*, double*),
    void (*metvec_l)(int*, int*, double*, double*),
    void (*precnd)(int*, int*, double*, double*, double*),
    const char* side,
    double* eig,
    double* evec_1,
    double* evec_2,
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

void dgl_smogd_driver_c(
#ifdef DGL_INT_KIND_4
    int n2,
    int n_targ,
    int n_max,
#elif DGL_INT_KIND_8
    long int n2,
    long int n_targ,
    long int n_max,
#endif
    void (*apbmul)(int*, int*, double*, double*),
    void (*ambmul)(int*, int*, double*, double*),
    void (*spdmul)(int*, int*, double*, double*),
    void (*smdmul)(int*, int*, double*, double*),
    void (*lrprec)(int*, int*, double*, double*, double*, double*, double*),
    double* eig,
    double* evec,
    bool* ok,
    bool verbose,
    double tol,
#ifdef DGL_INT_KIND_4
    int max_iter,
    int dav_iter,
    int memory,
#elif DGL_INT_KIND_8
    int long max_iter,
    int long dav_iter,
    int long memory,
#endif
    const char* memory_unit
);

#endif
