#ifndef DIAGLIB_H
#define DIAGLIB_H

#include <stdbool.h>
#include <stdint.h>
#include "diaglib_config.h"

/* Integer type matching the integers of the DiagLib build (see diaglib_config.h) */
#if DGL_INT_KIND == 8
typedef int64_t dgl_int;
#elif DGL_INT_KIND == 4
typedef int32_t dgl_int;
#else
#error "diaglib.h: unsupported DGL_INT_KIND, expected 4 or 8"
#endif

#ifdef __cplusplus
extern "C" {
#endif

/* Size in bytes of the integers of the DiagLib library actually linked (4 or 8).
 * It can be compared with DGL_INT_KIND to detect a mismatch between header and library. */
extern int dgl_integer_kind(void);

/*
 * Error codes returned in the info argument of the drivers. On error, ok is false
 * and eig/evec do not contain meaningful results. Not converging within max_iter
 * is not an error: info is DGL_SUCCESS and ok is false.
 */
#define DGL_SUCCESS       0  /* no error */
#define DGL_ERR_INPUT    -1  /* invalid input arguments */
#define DGL_ERR_MEMORY   -2  /* allocation failure or memory limit exceeded */
#define DGL_ERR_LAPACK   -3  /* a LAPACK routine failed */
#define DGL_ERR_ORTHO    -4  /* an orthogonalization procedure failed */
#define DGL_ERR_MISMATCH -5  /* left and right eigenvalues do not match (non-symmetric driver) */

extern void dgl_davidson_driver(
    dgl_int n,
    dgl_int n_targ,
    dgl_int n_max,
    void (*matvec)(dgl_int*, dgl_int*, double*, double*),
    void (*precnd)(dgl_int*, dgl_int*, double*, double*, double*),
    void (*metvec)(dgl_int*, dgl_int*, double*, double*),
    double* eig,
    double* evec,
    bool* ok,
    dgl_int* info,
    bool verbose,
    double tol,
    dgl_int max_iter,
    dgl_int dav_iter,
    double shift,       /* only added to the printed eigenvalues */
    bool precnd_shift,  /* shift passed to precnd: minus the lowest non-converged eigenvalue if true,
                           zero if false (the Fortran default is true) */
    dgl_int memory,
    const char* memory_unit
);

extern void dgl_lobpcg_driver(
    dgl_int n,
    dgl_int n_targ,
    dgl_int n_max,
    void (*matvec)(dgl_int*, dgl_int*, double*, double*),
    void (*precnd)(dgl_int*, dgl_int*, double*, double*, double*),
    void (*metvec)(dgl_int*, dgl_int*, double*, double*),
    double* eig,
    double* evec,
    bool* ok,
    dgl_int* info,
    bool verbose,
    double tol,
    dgl_int max_iter,
    double shift,       /* only added to the printed eigenvalues */
    bool precnd_shift,  /* shift passed to precnd: minus the lowest non-converged eigenvalue if true,
                           zero if false (the Fortran default is false: LOBPCG works best with a
                           positive definite, well conditioned preconditioner) */
    dgl_int memory,
    const char* memory_unit
);

extern void dgl_davidson_nosym_driver(
    dgl_int n,
    dgl_int n_targ,
    dgl_int n_max,
    void (*matvec_r)(dgl_int*, dgl_int*, double*, double*),
    void (*matvec_l)(dgl_int*, dgl_int*, double*, double*),
    void (*precnd)(dgl_int*, dgl_int*, double*, double*, double*),
    const char* side,
    double* eig,
    double* evec_1,
    double* evec_2,     /* only used if side is "LR", may be NULL otherwise. if converged, the first
                           n_targ left (evec_2) and right (evec_1) vectors are biorthonormal */
    bool* ok,
    dgl_int* info,
    bool verbose,
    double tol,
    dgl_int max_iter,
    dgl_int dav_iter,
    double shift,       /* only added to the printed eigenvalues */
    bool precnd_shift,  /* shift passed to precnd: minus the lowest non-converged eigenvalue if true,
                           zero if false (the Fortran default is true) */
    dgl_int memory,
    const char* memory_unit
);

extern void dgl_smogd_driver(
    dgl_int n2,
    dgl_int n_targ,
    dgl_int n_max,
    void (*apbmul)(dgl_int*, dgl_int*, double*, double*),
    void (*ambmul)(dgl_int*, dgl_int*, double*, double*),
    void (*spdmul)(dgl_int*, dgl_int*, double*, double*),
    void (*smdmul)(dgl_int*, dgl_int*, double*, double*),
    void (*lrprec)(dgl_int*, dgl_int*, double*, double*, double*, double*, double*),
    double* eig,
    double* evec,
    bool* ok,
    dgl_int* info,
    bool verbose,
    double tol,
    dgl_int max_iter,
    dgl_int dav_iter,
    dgl_int memory,
    const char* memory_unit
);

#ifdef __cplusplus
}
#endif

#endif
