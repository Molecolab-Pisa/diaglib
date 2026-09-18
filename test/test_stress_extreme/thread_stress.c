// Heavy concurrency stress test for the C interface, exercising the "no global state" claim in
// the README: many OpenMP threads call the drivers at the same time, repeatedly, each on its own
// independently-scaled problem, for a wall-clock time budget. Builds on test/test_c/test_c.c's
// test_threads (which only runs 8 tasks / 4 threads once): here every available core gets its own
// thread, and each thread loops for the whole budget with a different problem every iteration
// (varied by thread id and iteration count), so that any accidental cross-talk between concurrent
// calls (shared module state, mis-saved/restored callback pointers, etc.) has a real chance to show
// up as a wrong eigenvalue rather than being masked by lucky scheduling.
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "diaglib.h"

#define N 300
#define N_TARG 4
#define N_MAX 8
#define TOL 1e-9

// matvec: diagonal-dominant symmetric matrix, scaled by `scale` (thread/iteration dependent);
// eigenvalues are known in closed form up to O(1/n) perturbation from the off-diagonal part,
// so instead of a closed form we just check the residual (true test of correctness here).
typedef struct { double scale; } params_t;

static void matvec_c(dgl_int* n, dgl_int* m, double* x, double* ax, void* ctx) {
  params_t* p = (params_t*)ctx;
  int N_ = *n, M_ = *m;
  for (int j = 0; j < M_; ++j) {
    for (int i = 0; i < N_; ++i) {
      double sum = 0.0;
      for (int k = 0; k < N_; ++k) {
        double A_ik = (i == k) ? p->scale * (i + 2.0) : p->scale * 0.1 / (i + k + 2.0);
        sum += A_ik * x[k + j * N_];
      }
      ax[i + j * N_] = sum;
    }
  }
}

static void precnd_c(dgl_int* n, dgl_int* m, double* shift, double* x, double* z, void* ctx) {
  params_t* p = (params_t*)ctx;
  int N_ = *n, M_ = *m;
  for (int j = 0; j < M_; ++j)
    for (int i = 0; i < N_; ++i) {
      double diag = p->scale * (i + 2.0) + (*shift);
      if (fabs(diag) < 1e-3) diag = (diag < 0 ? -1 : 1) * 1e-3;
      z[i + j * N_] = x[i + j * N_] / diag;
    }
}

// diaglib's C matvec/precnd signatures take no user context (see diaglib.h); route the per-call
// scale through a thread-local instead, mirroring how a real caller with per-thread state would
// have to do it if it cannot use closures/context pointers in C.
static double tls_scale = 1.0;
#ifdef _OPENMP
#pragma omp threadprivate(tls_scale)
#endif

static void matvec_tls(dgl_int* n, dgl_int* m, double* x, double* ax) {
  params_t p = { tls_scale };
  matvec_c(n, m, x, ax, &p);
}
static void precnd_tls(dgl_int* n, dgl_int* m, double* shift, double* x, double* z) {
  params_t p = { tls_scale };
  precnd_c(n, m, shift, x, z, &p);
}

static double residual(dgl_int n, double scale, double* x, double lambda) {
  dgl_int one = 1;
  double* ax = malloc(n * sizeof(double));
  params_t p = { scale };
  matvec_c(&n, &one, x, ax, &p);
  double nr = 0.0, nx = 0.0;
  for (dgl_int i = 0; i < n; ++i) {
    double r = ax[i] - lambda * x[i];
    nr += r * r;
    nx += x[i] * x[i];
  }
  free(ax);
  return sqrt(nr / nx);
}

int main(int argc, char** argv) {
  double minutes_budget = argc > 1 ? atof(argv[1]) : 10.0;
#ifdef _OPENMP
  int nthreads = argc > 2 ? atoi(argv[2]) : omp_get_max_threads();
#else
  int nthreads = 1;
  printf("OpenMP not available: this is a SERIAL run, not a concurrency stress test\n");
#endif
  printf("thread_stress: %d threads, %.1f minute budget, N=%d\n", nthreads, minutes_budget, N);

  long total_iters = 0, total_bugs = 0, total_noconv = 0;
  time_t t0 = time(NULL);

#ifdef _OPENMP
  #pragma omp parallel num_threads(nthreads) reduction(+:total_iters,total_bugs,total_noconv)
#endif
  {
#ifdef _OPENMP
    int tid = omp_get_thread_num();
#else
    int tid = 0;
#endif
    long iter = 0;
    double eig[N_MAX];
    double* evec = malloc(N * N_MAX * sizeof(double));
    while (difftime(time(NULL), t0) < minutes_budget * 60.0) {
      double scale = 0.3 + 4.7 * ((double)((tid * 7919 + (int)iter * 104729) % 1000)) / 1000.0;
      tls_scale = scale;
      memset(evec, 0, N * N_MAX * sizeof(double));
      bool ok = false;
      dgl_int info = -99;
      if (iter % 2 == 0)
        dgl_davidson_driver(N, N_TARG, N_MAX, matvec_tls, precnd_tls, NULL, eig, evec, &ok, &info,
                            false, TOL, 200, 20, 0.0, true, 1, "GB");
      else
        dgl_lobpcg_driver(N, N_TARG, N_MAX, matvec_tls, precnd_tls, NULL, eig, evec, &ok, &info,
                          false, TOL, 200, 0.0, false, 1, "GB");
      total_iters++;
      bool bug = false;
      if (info != DGL_SUCCESS) {
        printf("thread %d iter %ld: *** BUG unexpected info=%ld (scale=%.4f)\n", tid, iter, (long)info, scale);
        bug = true;
      } else if (ok) {
        double max_res = 0.0;
        bool has_nan = false;
        for (int i = 0; i < N_TARG; ++i) {
          if (eig[i] != eig[i]) has_nan = true;
          max_res = fmax(max_res, residual(N, scale, evec + i * N, eig[i]));
        }
        if (has_nan || max_res > 1e-5) {
          printf("thread %d iter %ld: *** BUG wrong result, max_res=%.2e (scale=%.4f)%s\n", tid, iter,
                 max_res, scale, has_nan ? " NaN present" : "");
          bug = true;
        }
      } else {
        total_noconv++;
      }
      if (bug) total_bugs++;
      iter++;
    }
    free(evec);
  }

  printf("\n==================== thread_stress summary ====================\n");
  printf("threads: %d, total iterations: %ld, bugs: %ld, non-converged: %ld\n", nthreads, total_iters,
         total_bugs, total_noconv);
  printf(total_bugs == 0 ? "RESULT: CLEAN\n" : "RESULT: BUGS FOUND\n");
  return total_bugs == 0 ? 0 : 1;
}
