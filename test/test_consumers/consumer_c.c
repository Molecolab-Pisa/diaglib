/* C program using an installed DiagLib. Built by the CMake project in this directory, and
 * also compilable by hand with only -I<prefix>/include -L<prefix>/lib -ldiaglib_c.
 * Symmetric problem with known eigenvalues: A = s * H diag(d) H, H = I - 2 v v^T.
 * Non-symmetric: M = P diag(d) P^-1, P = I + u w^T. Threads are POSIX threads (not OpenMP). */
#include <math.h>
#include <pthread.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "diaglib.h"

#define N 1000
#define N_TARG 4
#define N_MAX 8
static double v[N], u[N], w[N], d[N], wu;
static _Thread_local double scale = 1.0;   /* per-thread matrix */
static int failures = 0;

static double dot(const double *a, const double *b) { double s = 0; for (int i = 0; i < N; ++i) s += a[i] * b[i]; return s; }

static void setup(void) {
  double nv = 0;
  for (int i = 0; i < N; ++i) {
    v[i] = sin(i + 1.0); nv += v[i] * v[i];
    u[i] = 0.3 * cos(2.0 * (i + 1)) / sqrt(N);
    w[i] = 0.5 * sin(3.0 * (i + 1)) / sqrt(N);
    d[i] = 1.0 + i + 0.5 * pow(sin(i + 1.0), 2);
  }
  for (int i = 0; i < N; ++i) v[i] /= sqrt(nv);
  wu = dot(w, u);
}

static void mv_sym(dgl_int *n, dgl_int *m, double *x, double *y) {
  for (dgl_int j = 0; j < *m; ++j) {
    double *xj = x + j * *n, *yj = y + j * *n;
    double c = 2.0 * dot(v, xj);
    for (int i = 0; i < N; ++i) yj[i] = xj[i] - c * v[i];
    for (int i = 0; i < N; ++i) yj[i] *= d[i];
    c = 2.0 * dot(v, yj);
    for (int i = 0; i < N; ++i) yj[i] = scale * (yj[i] - c * v[i]);
  }
}
static void pc_shift(dgl_int *n, dgl_int *m, double *shift, double *x, double *y) {
  for (dgl_int j = 0; j < *m; ++j)
    for (int i = 0; i < N; ++i) {
      double den = scale * d[i] + *shift;
      if (fabs(den) < 1e-3) den = den < 0 ? -1e-3 : 1e-3;
      y[i + j * *n] = x[i + j * *n] / den;
    }
}
static void pc_spd(dgl_int *n, dgl_int *m, double *shift, double *x, double *y) {
  (void)shift;
  for (dgl_int j = 0; j < *m; ++j)
    for (int i = 0; i < N; ++i) y[i + j * *n] = x[i + j * *n] / (scale * d[i]);
}
static void mv_r(dgl_int *n, dgl_int *m, double *x, double *y) {
  for (dgl_int j = 0; j < *m; ++j) {
    double *xj = x + j * *n, *yj = y + j * *n, c = dot(w, xj) / (1.0 + wu);
    for (int i = 0; i < N; ++i) yj[i] = d[i] * (xj[i] - c * u[i]);
    c = dot(w, yj);
    for (int i = 0; i < N; ++i) yj[i] += c * u[i];
  }
}
static void mv_l(dgl_int *n, dgl_int *m, double *x, double *y) {
  for (dgl_int j = 0; j < *m; ++j) {
    double *xj = x + j * *n, *yj = y + j * *n, c = dot(u, xj);
    for (int i = 0; i < N; ++i) yj[i] = d[i] * (xj[i] + c * w[i]);
    c = dot(u, yj) / (1.0 + wu);
    for (int i = 0; i < N; ++i) yj[i] -= c * w[i];
  }
}

static void check(const char *label, bool cond) {
  if (!cond) failures++;
  printf("%-45s %s\n", label, cond ? "PASSED" : "FAILED");
}
static double eig_err(const double *eig) {
  double e = 0;   /* d is increasing: the lowest eigenvalues are scale * d[0..] */
  for (int k = 0; k < N_TARG; ++k) e = fmax(e, fabs(eig[k] - scale * d[k]));
  return e;
}
static double residual(void (*mv)(dgl_int *, dgl_int *, double *, double *), const double *evec, const double *eig) {
  double r = 0, y[N], x[N]; dgl_int n = N, one = 1;
  for (int k = 0; k < N_TARG; ++k) {
    memcpy(x, evec + k * N, sizeof x);
    mv(&n, &one, x, y);
    double nr = 0;
    for (int i = 0; i < N; ++i) nr += pow(y[i] - eig[k] * x[i], 2);
    r = fmax(r, sqrt(nr / dot(x, x)));
  }
  return r;
}

static int nested_calls = 0;
static void mv_nested(dgl_int *n, dgl_int *m, double *x, double *y) {
  double e[4], ev[N * 4] = {0};
  bool ok; dgl_int info;
  __sync_fetch_and_add(&nested_calls, 1);
  dgl_lobpcg_driver(N, 2, 4, mv_sym, pc_spd, NULL, e, ev, &ok, &info, false, 1e-8, 200, 0.0, false, 100, "MB");
  if (!ok || info != DGL_SUCCESS || fabs(e[0] - scale * d[0]) > 1e-8) { printf("nested call failed\n"); exit(2); }
  mv_sym(n, m, x, y);
}

typedef struct { int id; bool pass; double err; } task_t;
static void *run_task(void *arg) {
  task_t *t = arg;
  double eig[N_MAX], *evec = calloc(N * N_MAX, sizeof(double));
  bool ok = false; dgl_int info = -99;
  scale = 1.0 + t->id;
  if (t->id % 2 == 0)
    dgl_davidson_driver(N, N_TARG, N_MAX, mv_nested, pc_shift, NULL, eig, evec, &ok, &info, false, 1e-9, 200, 20, 0.0, true, 100, "MB");
  else
    dgl_lobpcg_driver(N, N_TARG, N_MAX, mv_sym, pc_spd, NULL, eig, evec, &ok, &info, false, 1e-9, 300, 0.0, false, 100, "MB");
  t->err = eig_err(eig);
  t->pass = ok && info == DGL_SUCCESS && t->err < 1e-8;
  free(evec);
  return NULL;
}

int main(void) {
  double eig[N_MAX], evec[N * N_MAX], evec2[N * N_MAX];
  bool ok; dgl_int info;
  setup();
  printf("header integer kind %d, library %d\n", DGL_INT_KIND, dgl_integer_kind());
  check("integer kind matches", DGL_INT_KIND == dgl_integer_kind());

  memset(evec, 0, sizeof evec);
  dgl_davidson_driver(N, N_TARG, N_MAX, mv_sym, pc_shift, NULL, eig, evec, &ok, &info, false, 1e-9, 100, 25, 0.0, true, 100, "MB");
  printf("  eig err %.2e residual %.2e\n", eig_err(eig), residual(mv_sym, evec, eig));
  check("Davidson", ok && info == 0 && eig_err(eig) < 1e-8 && residual(mv_sym, evec, eig) < 1e-6);

  memset(evec, 0, sizeof evec);
  dgl_lobpcg_driver(N, N_TARG, N_MAX, mv_sym, pc_spd, NULL, eig, evec, &ok, &info, false, 1e-9, 300, 0.0, false, 100, "MB");
  check("LOBPCG", ok && info == 0 && eig_err(eig) < 1e-8 && residual(mv_sym, evec, eig) < 1e-6);

  memset(evec, 0, sizeof evec); memset(evec2, 0, sizeof evec2);
  dgl_davidson_nosym_driver(N, N_TARG, N_MAX, mv_r, mv_l, pc_shift, "LR", eig, evec, evec2, &ok, &info, false, 1e-9, 100, 25, 0.0, true, 100, "MB");
  double bi = 0;
  for (int k = 0; k < N_TARG; ++k) for (int l = 0; l < N_TARG; ++l) bi = fmax(bi, fabs(dot(evec2 + k * N, evec + l * N) - (k == l)));
  printf("  eig err %.2e residuals %.2e %.2e biortho %.2e\n", eig_err(eig), residual(mv_r, evec, eig), residual(mv_l, evec2, eig), bi);
  check("non-symmetric LR", ok && info == 0 && eig_err(eig) < 1e-8 && residual(mv_r, evec, eig) < 1e-6 &&
                            residual(mv_l, evec2, eig) < 1e-6 && bi < 1e-10);

  memset(evec, 0, sizeof evec);
  dgl_davidson_nosym_driver(N, N_TARG, N_MAX, mv_r, mv_l, pc_shift, "L", eig, evec, NULL, &ok, &info, false, 1e-9, 100, 25, 0.0, true, 100, "MB");
  check("non-symmetric L (evec_2 = NULL)", ok && info == 0 && eig_err(eig) < 1e-8 && residual(mv_l, evec, eig) < 1e-6);

  /* errors */
  dgl_davidson_driver(N, N_TARG, N_MAX, NULL, pc_shift, NULL, eig, evec, &ok, &info, false, 1e-9, 100, 25, 0.0, true, 100, "MB");
  check("error: NULL matvec", !ok && info == DGL_ERR_INPUT);
  dgl_davidson_nosym_driver(N, N_TARG, N_MAX, mv_r, mv_l, pc_shift, "LR", eig, evec, NULL, &ok, &info, false, 1e-9, 100, 25, 0.0, true, 100, "MB");
  check("error: LR with NULL evec_2", !ok && info == DGL_ERR_INPUT);
  dgl_davidson_nosym_driver(N, N_TARG, N_MAX, mv_r, mv_l, pc_shift, "X", eig, evec, NULL, &ok, &info, false, 1e-9, 100, 25, 0.0, true, 100, "MB");
  check("error: invalid side", !ok && info == DGL_ERR_INPUT);
  dgl_davidson_driver(N, N_TARG, N_MAX, mv_sym, pc_shift, NULL, eig, evec, &ok, &info, false, -1.0, 100, 25, 0.0, true, 100, "MB");
  check("error: negative tol", !ok && info == DGL_ERR_INPUT);
  dgl_davidson_driver(N, N_TARG, N_MAX, mv_sym, pc_shift, NULL, eig, evec, &ok, &info, false, 1e-9, 100, 25, 0.0, true, 100, NULL);
  check("NULL memory unit (default MB)", ok && info == 0);
  memset(evec, 0, sizeof evec);
  dgl_davidson_driver(N, N_TARG, N_MAX, mv_sym, pc_shift, NULL, eig, evec, &ok, &info, false, 1e-12, 2, 25, 0.0, true, 100, "MB");
  check("not converged: ok false, info 0", !ok && info == 0);

  /* concurrent calls from POSIX threads, with nested calls in half of them */
  enum { N_THREADS = 6 };
  pthread_t th[N_THREADS]; task_t tasks[N_THREADS];
  for (int t = 0; t < N_THREADS; ++t) { tasks[t].id = t; pthread_create(&th[t], NULL, run_task, &tasks[t]); }
  bool all = true; double maxerr = 0;
  for (int t = 0; t < N_THREADS; ++t) { pthread_join(th[t], NULL); all = all && tasks[t].pass; maxerr = fmax(maxerr, tasks[t].err); }
  printf("  nested calls %d, max eigenvalue error %.2e\n", nested_calls, maxerr);
  check("6 concurrent pthreads (nested calls in 3)", all && nested_calls > 0);

  printf("C consumer: %d failures\n", failures);
  return failures != 0;
}
