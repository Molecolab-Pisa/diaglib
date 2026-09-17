#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "diaglib.h"

//
// simple minded test driver for diaglib called from c.
//
// list of functions:
// ==================
//
// matvec_c(dgl_int* n, dgl_int* m, double*x, double* ax)
//
//   applies the symmetric matrix 
//
//     a_{ij} = \delta_{ij} (i+1) + (1-\delta_{ij})/(i+j) 
//
//   to m vectors of length n (in x, assumed column major
//   as in fortran)
//
//
// precnd_c(dgl_int* n, dgl_int* m, double* shift, double*x, double* ax)
//
//   applies the preconditioner
//
//     p_{ij} = \delta_{ij} (a_{ii} + shift)
//
//   to m vectors of length n (in x, assumed column major
//   as in fortran), where a is the matrix from matvec_c.
//
// matvec_r_c(dgl_int* n, dgl_int* m, double* x, double* y)
// matvec_l_c(dgl_int* n, dgl_int* m, double* x, double* y)
//
//   apply the non-symmetric matrix
//
//     m = d^{-1} a d (matvec_r_c)
//
//   or its transpose 
//
//     m^t = d a d^{-1} (matvec_l_c)
//
//   to m vectors of length n (in x, assumed column major
//   as in fortran), where a is the matrix from matvec_c,
//   and d is a diagonal matrix defined as 
//
//     d_{ii} = 1/i
//    
// metvec_c(dgl_int* n, dgl_int* m, double* x, double* bx)
//
//   applies a positive definite matrix b (here, just the identity)
//   to m vectors of length n (in x, assumed column major
//   as in fortran)
// 
// apbmul_c(dgl_int* n, dgl_int* m, double* x, double* ax)
// ambmul_c(dgl_int* n, dgl_int* m, double* x, double* ax)
//
//   apply the SPD matrices (a+b) and (a-b), respectively, where
//
//     (a + b)_{ij} = \delta_{ij} (i+5) + (1-\delta_{ij})/(i+j)
//     (a - b)_{ij} = \delta_{ij} (i+2) + 0.2 * (1-\delta_{ij})/(i+j)
//
//   to m vectors of length n (in x, assumed column major
//   as in fortran). There matrices have the correct structure and 
//   signature for linear response problems and can be used to test
//   the smogd routine.
//
// spdmul_c(dgl_int* n, dgl_int* m, double* x, double* ax)
// smdmul_c(dgl_int* n, dgl_int* m, double* x, double* ax)
//
//   apply the matrices (\sigma+\delta) and (\sigma-\delta), respectively, where
//
//     \sigma_{ij} = \delta_{ij}
//     \delta_{ij} = 0,    i = j;
//                   0.05, i > j;
//                  -0.05, i < j.
//
//   to m vectors of length n (in x, assumed column major
//   as in fortran). 
//   note that sigma is SPD, while delta is antisymmetric, as they are
//   in linear response problems. these routines can thus be used to
//   test smogd. 
// 
// lrprec_c(dgl_int* n, dgl_int* m, double* fac, double* xp, double* xm,  double* yp, double* ym)
//
//   applies the smogd diagonal preconditioner to m vectors xp and xm (note that the two
//   quantities are used to define both yp and ym), which are assumed to be column major
//   as in fortran. 
//   this routine implements the same definitions of a and sigma as in the previous 
//   routines, and can be used to test smogd. 
//
void matvec_c(dgl_int* n, dgl_int* m, double* x, double* ax) {
  int N = *n;
  int M = *m;

  for (int j = 0; j < M; ++j) {
    for (int i = 0; i < N; ++i) {
      double sum = 0.0;
      for (int k = 0; k < N; ++k) {
        double A_ik = (i == k) ? (i + 2.0) : (1.0 / (i + k + 2.0)); // i,k zero-based
        sum += A_ik * x[k + j * N];
      }
      ax[i + j * N] = sum;
    }
  }
}

void precnd_c(dgl_int* n, dgl_int* m, double* shift, double* r, double* z) {
  int N = *n;
  int M = *m;

  for (int j = 0; j < M; ++j) {
    for (int i = 0; i < N; ++i) {
      double diag = (i + 2.0) + (*shift); // A_ii + shift
      z[i + j * N] = r[i + j * N] / diag;
    }
  }
}

void matvec_r_c(dgl_int* n, dgl_int* m, double* x, double* y) {
    int N = *n, M = *m;

    for (int k = 0; k < M; ++k)
        for (int i = 0; i < N; ++i) {
            double sum = 0.0;
            for (int j = 0; j < N; ++j) {
                double xval = x[j + N * k];
                sum += (j == i)
                    ? (i + 2.0) * xval
                    : ((i + 1.0) / (j + 1.0)) * xval / (i + j + 2.0);
            }
            y[i + N * k] = sum;
        }
}

void matvec_l_c(dgl_int* n, dgl_int* m, double* x, double* y) {
    int N = *n, M = *m;

    for (int k = 0; k < M; ++k)
        for (int i = 0; i < N; ++i) {
            double sum = 0.0;
            for (int j = 0; j < N; ++j) {
                double xval = x[j + N * k];
                sum += (j == i)
                    ? (i + 2.0) * xval
                    : ((j + 1.0) / (i + 1.0)) * xval / (i + j + 2.0);
            }
            y[i + N * k] = sum;
        }
}

void metvec_c(dgl_int* n, dgl_int* m, double* x, double* bx) {
  int N = *n;
  int M = *m;
  for (int j = 0; j < M; ++j)
    for (int i = 0; i < N; ++i)
      bx[i + j * N] = x[i + j * N];  
}

void apbmul_c(dgl_int* n, dgl_int* m, double* x, double* ax) {
  int N = *n;
  int M = *m;
  for (int j = 0; j < M; ++j) {
    for (int i = 0; i < N; ++i) {
      double sum = 0.0;
      for (int k = 0; k < N; ++k) {
        double val = (i == k) ? (6.0 + i) : 1.0 / (i + k + 2.0);
        sum += val * x[k + j * N];
      }
      ax[i + j * N] = sum;
    }
  }
}

void ambmul_c(dgl_int* n, dgl_int* m, double* x, double* ax) {
  int N = *n;
  int M = *m;
  for (int j = 0; j < M; ++j) {
    for (int i = 0; i < N; ++i) {
      double sum = 0.0;
      for (int k = 0; k < N; ++k) {
        double val = (i == k) ? (3.0 + i) : 0.2 / (i + k + 2.0);
        sum += val * x[k + j * N];
      }
      ax[i + j * N] = sum;
    }
  }
}

void spdmul_c(dgl_int* n, dgl_int* m, double* x, double* y) {
    int N = *n, M = *m;

    for (int k = 0; k < M; ++k)
        for (int i = 0; i < N; ++i) {
            double sum = 0.0;
            for (int j = 0; j < N; ++j) {
                double xval = x[j + k * N];
                sum += (i == j)     ? xval
                     : (i > j)      ? 0.05 * xval
                                   : -0.05 * xval;
            }
            y[i + k * N] = sum;
        }
}

void smdmul_c(dgl_int* n, dgl_int* m, double* x, double* y) {
    int N = *n, M = *m;

    for (int k = 0; k < M; ++k)
        for (int i = 0; i < N; ++i) {
            double sum = 0.0;
            for (int j = 0; j < N; ++j) {
                double xval = x[j + k * N];
                sum += (i == j)     ? xval
                     : (i > j)      ? -0.05 * xval
                                   : 0.05 * xval;
            }
            y[i + k * N] = sum;
        }
}

void lrprec_c(dgl_int* n, dgl_int* m, double* fac, double* xp, double* xm, 
              double* yp, double* ym) {
  int N = *n;
  int M = *m;
  double f = *fac;
  for (int j = 0; j < M; ++j) {
    for (int i = 0; i < N; ++i) {
      double val = 1.0 / ( f * f * (i+8.0) * (i+8.0) - 1.0 );
      double fp  =  val * (f * (i + 8.0 ) * xp[i + j * N] + xm[i + j * N]);
      double fm  =  val * (f * (i + 8.0 ) * xm[i + j * N] + xp[i + j * N]);
      yp[i + j * N] = fp;
      ym[i + j * N] = fm;
    }
  }
}

//
// ==========================================================================
// tests: every driver is run from a simple guess (unit vectors) and from a
// random guess (zero vectors, so that diaglib generates the guess). the
// converged eigenvalues are compared with reference values computed with
// lapack (see test/test_fortran/reference.f90), and the true residual of
// every returned eigenvector is computed with the matvecs above.
// ==========================================================================
//
#define N 700
#define N_TARG 5
#define N_MAX 10
#define TOL 1e-10
#define EIG_THRESH 1e-8
#define RES_THRESH 1e-7

// matrices applied by matvec_c (and, similar to it, by matvec_r_c / matvec_l_c)
static const double ref_sym[N_TARG] = {1.8693982134576630, 3.0004763487609458, 4.0177129623236798,
                                       5.0168124845115560, 6.0135237881685546};
// linear response problem defined by apbmul_c, ambmul_c, spdmul_c, smdmul_c
static const double ref_lr[N_TARG] = {3.9540942156164007, 4.7820108108607249, 5.6867559690540226,
                                      6.5222631332836727, 7.1025644429415040};

static int n_tests = 0, n_failed = 0;

typedef void (*matvec_t)(dgl_int*, dgl_int*, double*, double*);

// relative residual ||A x - lambda x|| / ||x|| of one vector
static double residual(matvec_t apply, dgl_int n, double* x, double lambda) {
  dgl_int one = 1;
  double* ax = malloc(n * sizeof(double));
  apply(&n, &one, x, ax);
  double nr = 0.0, nx = 0.0;
  for (dgl_int i = 0; i < n; ++i) {
    double r = ax[i] - lambda * x[i];
    nr += r * r;
    nx += x[i] * x[i];
  }
  free(ax);
  return sqrt(nr / nx);
}

// relative residual ||E z - omega M z|| / ||z|| of the linear response problem, with
// E = [A B; B A], M = [S D; -D -S], z = (y, x), and A+B, A-B, S+D, S-D given by the matvecs
static double lr_residual(dgl_int n, double* z, double omega) {
  dgl_int one = 1;
  double* p = malloc(n * sizeof(double)), *m = malloc(n * sizeof(double));
  double* ap = malloc(n * sizeof(double)), *am = malloc(n * sizeof(double));
  double* sp = malloc(n * sizeof(double)), *sm = malloc(n * sizeof(double));
  for (dgl_int i = 0; i < n; ++i) {
    p[i] = z[i] + z[n + i];
    m[i] = z[i] - z[n + i];
  }
  apbmul_c(&n, &one, p, ap);
  ambmul_c(&n, &one, m, am);
  spdmul_c(&n, &one, p, sp);
  smdmul_c(&n, &one, m, sm);
  double nr = 0.0, nz = 0.0;
  for (dgl_int i = 0; i < n; ++i) {
    double r1 = 0.5 * (ap[i] + am[i]) - omega * 0.5 * (sp[i] + sm[i]);
    double r2 = 0.5 * (ap[i] - am[i]) + omega * 0.5 * (sp[i] - sm[i]);
    nr += r1 * r1 + r2 * r2;
    nz += z[i] * z[i] + z[n + i] * z[n + i];
  }
  free(p); free(m); free(ap); free(am); free(sp); free(sm);
  return sqrt(nr / nz);
}

static void init_guess(dgl_int ld, double* evec, bool random_guess) {
  memset(evec, 0, ld * N_MAX * sizeof(double));
  if (!random_guess)
    for (dgl_int j = 0; j < N_MAX; ++j) evec[j + j * ld] = 1.0;
}

// record the result of a test
static void check(const char* label, bool random_guess, bool ok, dgl_int info, const double* eig,
                  const double* ref, double max_res) {
  double max_err = 0.0;
  for (int i = 0; i < N_TARG; ++i) max_err = fmax(max_err, fabs(eig[i] - ref[i]));
  bool pass = ok && info == DGL_SUCCESS && max_err < EIG_THRESH && max_res < RES_THRESH;
  n_tests++;
  if (!pass) n_failed++;
  printf("%-45s %s guess: ok=%d info=%ld max eigenvalue error=%.1e max residual=%.1e -> %s\n", label,
         random_guess ? "random" : "simple", ok, (long)info, max_err, max_res, pass ? "PASSED" : "FAILED");
}

static void check_error(const char* label, dgl_int info, dgl_int expected) {
  bool pass = info == expected;
  n_tests++;
  if (!pass) n_failed++;
  printf("%-45s info=%ld (expected %ld) -> %s\n", label, (long)info, (long)expected, pass ? "PASSED" : "FAILED");
}

static void test_davidson(bool generalized, bool random_guess) {
  double eig[N_MAX], evec[N * N_MAX], max_res = 0.0;
  bool ok = false;
  dgl_int info = -99;
  init_guess(N, evec, random_guess);
  dgl_davidson_driver(N, N_TARG, N_MAX, matvec_c, precnd_c, generalized ? metvec_c : NULL, eig, evec, &ok, &info,
                      false, TOL, 100, 20, 0.0, 1, "GB");
  // the metric is the identity: the residual is the same as for the standard problem
  for (int i = 0; i < N_TARG && ok; ++i) max_res = fmax(max_res, residual(matvec_c, N, evec + i * N, eig[i]));
  check(generalized ? "generalized Davidson" : "Davidson", random_guess, ok, info, eig, ref_sym, max_res);
}

static void test_lobpcg(bool generalized, bool random_guess) {
  double eig[N_MAX], evec[N * N_MAX], max_res = 0.0;
  bool ok = false;
  dgl_int info = -99;
  init_guess(N, evec, random_guess);
  dgl_lobpcg_driver(N, N_TARG, N_MAX, matvec_c, precnd_c, generalized ? metvec_c : NULL, eig, evec, &ok, &info,
                    false, TOL, 100, 0.0, 1, "GB");
  for (int i = 0; i < N_TARG && ok; ++i) max_res = fmax(max_res, residual(matvec_c, N, evec + i * N, eig[i]));
  check(generalized ? "generalized LOBPCG" : "LOBPCG", random_guess, ok, info, eig, ref_sym, max_res);
}

// the non-symmetric matrix of matvec_r_c is similar to the one of matvec_c:
// the eigenvalues are the same
static void test_nosym(const char* side, bool random_guess) {
  double eig[N_MAX], evec[N * N_MAX], evec_2[N * N_MAX], max_res = 0.0;
  bool ok = false, lr = strcmp(side, "LR") == 0;
  dgl_int info = -99;
  char label[64];
  init_guess(N, evec, random_guess);
  init_guess(N, evec_2, random_guess);
  dgl_davidson_nosym_driver(N, N_TARG, N_MAX, matvec_r_c, matvec_l_c, precnd_c, side, eig, evec,
                            lr ? evec_2 : NULL, &ok, &info, false, TOL, 100, 20, 0.0, 1, "GB");
  // in evec: right eigenvectors for "R" and "LR", left ones for "L"; in evec_2: left ones for "LR"
  for (int i = 0; i < N_TARG && ok; ++i) {
    matvec_t first = strcmp(side, "L") == 0 ? matvec_l_c : matvec_r_c;
    max_res = fmax(max_res, residual(first, N, evec + i * N, eig[i]));
    if (lr) max_res = fmax(max_res, residual(matvec_l_c, N, evec_2 + i * N, eig[i]));
  }
  // left and right eigenvectors have to be biorthonormal (the error is added to the residual)
  for (int i = 0; i < N_TARG && ok && lr; ++i)
    for (int j = 0; j < N_TARG; ++j) {
      double dot = 0.0;
      for (int k = 0; k < N; ++k) dot += evec_2[k + i * N] * evec[k + j * N];
      max_res = fmax(max_res, fabs(dot - (i == j ? 1.0 : 0.0)));
    }
  snprintf(label, sizeof label, "non-symmetric Davidson (%s)", side);
  check(label, random_guess, ok, info, eig, ref_sym, max_res);
}

static void test_smogd(bool random_guess) {
  double eig[N_MAX], evec[2 * N * N_MAX], max_res = 0.0;
  bool ok = false;
  dgl_int info = -99;
  init_guess(2 * N, evec, random_guess);
  dgl_smogd_driver(2 * N, N_TARG, N_MAX, apbmul_c, ambmul_c, spdmul_c, smdmul_c, lrprec_c, eig, evec, &ok, &info,
                   false, TOL, 100, 20, 1, "GB");
  for (int i = 0; i < N_TARG && ok; ++i) max_res = fmax(max_res, lr_residual(N, evec + i * 2 * N, eig[i]));
  check("SMO-GD", random_guess, ok, info, eig, ref_lr, max_res);
}

static void test_errors(void) {
  double eig[N_MAX], evec[N * N_MAX];
  bool ok;
  dgl_int info;
  init_guess(N, evec, false);
  dgl_davidson_driver(N, N_MAX + 1, N_MAX, matvec_c, precnd_c, NULL, eig, evec, &ok, &info, false, TOL, 100, 20, 0.0,
                      1, "GB");
  check_error("error: n_targ > n_max", info, DGL_ERR_INPUT);
  dgl_lobpcg_driver(N, N_TARG, N_MAX, NULL, precnd_c, NULL, eig, evec, &ok, &info, false, TOL, 100, 0.0, 1, "GB");
  check_error("error: NULL matvec", info, DGL_ERR_INPUT);
  dgl_davidson_nosym_driver(N, N_TARG, N_MAX, matvec_r_c, matvec_l_c, precnd_c, "LR", eig, evec, NULL, &ok, &info,
                            false, TOL, 100, 20, 0.0, 1, "GB");
  check_error("error: side = LR and evec_2 = NULL", info, DGL_ERR_INPUT);
  dgl_davidson_driver(N, N_TARG, N_MAX, matvec_c, precnd_c, NULL, eig, evec, &ok, &info, false, TOL, 100, 20, 0.0,
                      1, "KB");
  check_error("error: not enough memory", info, DGL_ERR_MEMORY);
}

// DiagLib called from inside a callback: nested_matvec_c solves a smaller problem with
// a different matvec (after a call with invalid input), and then applies the matrix of
// matvec_c. The callbacks of the outer call have to be restored after every nested call.
static int nested_calls = 0, nested_failed = 0;

static void nested_matvec_c(dgl_int* n, dgl_int* m, double* x, double* ax) {
  enum { N_IN = 100, N_TARG_IN = 2, N_MAX_IN = 4 };
  double eig[N_MAX_IN], evec[N_IN * N_MAX_IN];
  bool ok = true;
  dgl_int info = -99;
  nested_calls++;
  memset(evec, 0, sizeof evec);
  for (int j = 0; j < N_MAX_IN; ++j) evec[j + j * N_IN] = 1.0;
  dgl_davidson_driver(N_IN, N_MAX_IN + 1, N_MAX_IN, apbmul_c, precnd_c, NULL, eig, evec, &ok, &info, false, 1e-8,
                      100, 20, 0.0, 1, "GB");
  bool passed = !ok && info == DGL_ERR_INPUT;
  dgl_davidson_driver(N_IN, N_TARG_IN, N_MAX_IN, apbmul_c, precnd_c, NULL, eig, evec, &ok, &info, false, 1e-8,
                      100, 20, 0.0, 1, "GB");
  passed = passed && ok && info == DGL_SUCCESS && residual(apbmul_c, N_IN, evec, eig[0]) < 1e-6;
  if (!passed) nested_failed++;
  matvec_c(n, m, x, ax);
}

static void test_nested(void) {
  double eig[N_MAX], evec[N * N_MAX], max_res = 0.0;
  bool ok = false;
  dgl_int info = -99;
  init_guess(N, evec, false);
  dgl_davidson_driver(N, N_TARG, N_MAX, nested_matvec_c, precnd_c, metvec_c, eig, evec, &ok, &info, false, TOL, 100,
                      20, 0.0, 1, "GB");
  for (int i = 0; i < N_TARG && ok; ++i) max_res = fmax(max_res, residual(matvec_c, N, evec + i * N, eig[i]));
  printf("nested calls: %d, failed: %d\n", nested_calls, nested_failed);
  if (nested_calls == 0 || nested_failed > 0) ok = false;
  check("generalized Davidson, nested calls", false, ok, info, eig, ref_sym, max_res);
}

#ifdef _OPENMP
// DiagLib called at the same time from different threads, with different callbacks:
// the even tasks solve the problem of matvec_c, the odd ones the one of twice that matrix.
static void matvec2_c(dgl_int* n, dgl_int* m, double* x, double* ax) {
  matvec_c(n, m, x, ax);
  for (dgl_int i = 0; i < *n * *m; ++i) ax[i] *= 2.0;
}

static void precnd2_c(dgl_int* n, dgl_int* m, double* shift, double* r, double* z) {
  double half_shift = 0.5 * *shift;
  precnd_c(n, m, &half_shift, r, z);
  for (dgl_int i = 0; i < *n * *m; ++i) z[i] *= 0.5;
}

static void test_threads(void) {
  enum { N_TASKS = 8 };
  bool ok[N_TASKS];
  dgl_int info[N_TASKS];
  double max_err[N_TASKS];
  #pragma omp parallel for num_threads(4) schedule(dynamic, 1)
  for (int task = 0; task < N_TASKS; ++task) {
    double eig[N_MAX], *evec = malloc(N * N_MAX * sizeof(double));
    double scale = task % 2 == 0 ? 1.0 : 2.0;
    matvec_t mv = task % 2 == 0 ? matvec_c : matvec2_c;
    void (*pc)(dgl_int*, dgl_int*, double*, double*, double*) = task % 2 == 0 ? precnd_c : precnd2_c;
    init_guess(N, evec, false);
    ok[task] = false;
    info[task] = -99;
    if (task % 4 < 2)
      dgl_davidson_driver(N, N_TARG, N_MAX, mv, pc, NULL, eig, evec, &ok[task], &info[task], false, TOL, 100, 20,
                          0.0, 1, "GB");
    else
      dgl_lobpcg_driver(N, N_TARG, N_MAX, mv, pc, NULL, eig, evec, &ok[task], &info[task], false, TOL, 100, 0.0, 1,
                        "GB");
    max_err[task] = 0.0;
    for (int i = 0; i < N_TARG; ++i) max_err[task] = fmax(max_err[task], fabs(eig[i] - scale * ref_sym[i]));
    free(evec);
  }
  for (int task = 0; task < N_TASKS; ++task) {
    bool pass = ok[task] && info[task] == DGL_SUCCESS && max_err[task] < EIG_THRESH;
    n_tests++;
    if (!pass) n_failed++;
    printf("%-45s task %d: ok=%d info=%ld max eigenvalue error=%.1e -> %s\n",
           task % 4 < 2 ? "Davidson, concurrent calls" : "LOBPCG, concurrent calls", task, ok[task], (long)info[task],
           max_err[task], pass ? "PASSED" : "FAILED");
  }
}
#endif

//
// main program: test davidson, non-symmetric davidson, lobpcg and smogd.
//
int main(void) {
  n_tests++;
  if (dgl_integer_kind() != DGL_INT_KIND) n_failed++;
  printf("integer kind of the library: %d, of the header: %d -> %s\n", dgl_integer_kind(), DGL_INT_KIND,
         dgl_integer_kind() == DGL_INT_KIND ? "PASSED" : "FAILED");

  for (int guess = 0; guess < 2; ++guess) {
    bool random_guess = guess == 1;
    test_davidson(false, random_guess);
    test_davidson(true, random_guess);
    test_lobpcg(false, random_guess);
    test_lobpcg(true, random_guess);
    test_nosym("R", random_guess);
    test_nosym("L", random_guess);
    test_nosym("LR", random_guess);
    test_smogd(random_guess);
  }
  test_errors();
  test_nested();
#ifdef _OPENMP
  test_threads();
#else
  printf("OpenMP not available: concurrent calls are not tested\n");
#endif

  printf("Summary: %d failed tests out of %d.\n", n_failed, n_tests);
  return n_failed == 0 ? 0 : 1;
}
