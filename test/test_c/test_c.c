#include <stdio.h>
#include <stdbool.h>
#include "diaglib.h"

void matvec_c(int* n, int* m, double* x, double* ax) {
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

void precnd_c(int* n, int* m, double* shift, double* r, double* z) {
  int N = *n;
  int M = *m;

  for (int j = 0; j < M; ++j) {
    for (int i = 0; i < N; ++i) {
      double diag = (i + 1.0) + (*shift); // A_ii + shift
      z[i + j * N] = r[i + j * N] / diag;
    }
  }
}

void matvec_r_c(int* n, int* m, double* x, double* y) {
    int N = *n;
    int M = *m;
    int i, j, k;
    double fac;

    for (i = 0; i < N; ++i)
        for (k = 0; k < M; ++k)
            y[i + N * k] = 0.0;

    for (k = 0; k < M; ++k) {
        for (i = 0; i < N; ++i) {
            for (j = 0; j < N; ++j) {
                if (j == i) {
                    y[i + N * k] += (i + 1.0 + 1.0) * x[j + N * k];
                } else {
                    fac = (i + 1.0) / (j + 1.0);
                    y[i + N * k] += fac * x[j + N * k] / (i + j + 2.0);
                }
            }
        }
    }
}

void matvec_l_c(int* n, int* m, double* x, double* y) {
    int i, j, k;
    double fac;
    int N = *n;
    int M = *m;

    for (i = 0; i < N; ++i)
        for (k = 0; k < M; ++k)
            y[i + N * k] = 0.0;

    for (k = 0; k < M; ++k) {
        for (i = 0; i < N; ++i) {
            for (j = 0; j < N; ++j) {
                if (j == i) {
                    y[i + N * k] += (i + 1.0 + 1.0) * x[j + N * k];
                } else {
                    fac = (j + 1.0) / (i + 1.0);
                    y[i + N * k] += fac * x[j + N * k] / (i + j + 2.0);
                }
            }
        }
    }
}


void bvec_c(int* n, int* m, double* x, double* bx) {
  int N = *n;
  int M = *m;
  for (int j = 0; j < M; ++j)
    for (int i = 0; i < N; ++i)
      bx[i + j * N] = x[i + j * N];  // identità
}

void apbmul_c(int* n, int* m, double* x, double* ax) {
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

void ambmul_c(int* n, int* m, double* x, double* ax) {
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

void spdmul_c(int* n, int* m, double* x, double* y) {
  int N = *n;
  int M = *m;

  // Inizializza y a zero
  for (int j = 0; j < M; ++j) {
    for (int i = 0; i < N; ++i) {
      y[i + j * N] = 0.0;
    }
  }

  // Applica la logica della matrice SPD modificata
  for (int k = 0; k < M; ++k) {
    for (int j = 0; j < N; ++j) {
      for (int i = 0; i < N; ++i) {
        if (i == j) {
          y[i + k * N] += x[i + k * N];
        } else if (i > j) {
          y[i + k * N] += 0.05 * x[j + k * N];
        } else {
          y[i + k * N] -= 0.05 * x[j + k * N];
        }
      }
    }
  }
}

void smdmul_c(int* n, int* m, double* x, double* y) {
  int N = *n;
  int M = *m;

  // Inizializza y a zero
  for (int j = 0; j < M; ++j) {
    for (int i = 0; i < N; ++i) {
      y[i + j * N] = 0.0;
    }
  }

  // Applica la logica della matrice SPD modificata
  for (int k = 0; k < M; ++k) {
    for (int j = 0; j < N; ++j) {
      for (int i = 0; i < N; ++i) {
        if (i == j) {
          y[i + k * N] += x[i + k * N];
        } else if (i > j) {
          y[i + k * N] -= 0.05 * x[j + k * N];
        } else {
          y[i + k * N] += 0.05 * x[j + k * N];
        }
      }
    }
  }
}

void lrprec_c(int* n, int* m, double* fac, double* xp, double* xm, 
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
  

int main() {
  const int n = 1000, n_targ = 20, n_max = 25, max_iter = 100, max_dav = 20;
  const double tol = 1e-8, shift = 0.0;
  double eig[n_max];
  double evec[n * n_max], evec_l[n * n_max];
  bool ok;

//  printf("\nCalling DAVIDSON driver...\n");
//  for (int i = 0; i < n * n_max; ++i)
//    evec[i] = 0.0;
//  for (int j = 0; j < n_max; ++j)
//    for (int i = 0; i < n; ++i)
//      evec[i + j * n] = (j == i) ? 1.0 : 0.0;
//
//  davidson_driver_c(true, n, n_targ, n_max, max_iter, max_dav, tol, shift,
//                    matvec_c, precnd_c, eig, evec, &ok);
//
//  if (ok) {
//    printf("Davidson converged.\nEigenvalues:\n");
//    for (int i = 0; i < n_targ; ++i)
//      printf("  %.8f\n", eig[i]);
//  } else {
//    printf("Davidson failed to converge.\n");
//  }

  printf("\nCalling non-symmetric DAVIDSON driver...\n");
  for (int i = 0; i < n * n_max; ++i)
    evec[i] = 0.0;
  for (int j = 0; j < n_max; ++j)
    for (int i = 0; i < n; ++i)
      evec[i + j * n] = (j == i) ? 1.0 : 0.0;
  for (int i = 0; i < n * n_max; ++i)
    evec_l[i] = 0.0;
  for (int j = 0; j < n_max; ++j)
    for (int i = 0; i < n; ++i)
      evec_l[i + j * n] = (j == i) ? 1.0 : 0.0;

  nonsym_driver_c(true, n, n_targ, n_max, max_iter, tol, max_dav, shift,
                  matvec_c, matvec_c, precnd_c, eig, evec, evec_l, 
//                  matvec_r_c, matvec_l_c, precnd_c, eig, evec, evec_l, 
                  4, &ok);

  if (ok) {
    printf("Non-symmetric Davidson converged.\nEigenvalues:\n");
    for (int i = 0; i < n_targ; ++i)
      printf("  %.8f\n", eig[i]);
  } else {
    printf("Non-symmetric Davidson failed to converge.\n");
  }

  printf("\nCalling LOBPCG driver...\n");

  for (int i = 0; i < n * n_max; ++i)
    evec[i] = 0.0;
  for (int j = 0; j < n_max; ++j)
    for (int i = 0; i < n; ++i)
      evec[i + j * n] = (j == i) ? 1.0 : 0.0;

  ok = false;
  
  lobpcg_driver_c(
    true,         // verbose
    false,        // gen_eig = false → problema standard
    n, n_targ, n_max, max_iter, tol, shift,
    matvec_c, precnd_c, bvec_c,
    eig, evec, &ok
  );
  
  if (ok) {
    printf("LOBPCG converged.\nEigenvalues:\n");
    for (int i = 0; i < n_targ; ++i)
      printf("  %.8f\n", eig[i]);
  } else {
    printf("LOBPCG failed to converge.\n");
  }

printf("\nCalling CASLR-EFF driver...\n");

  int n2 = 2 * n;
  double* evec2 = malloc(sizeof(double) * n2 * n_max);
  for (int i = 0; i < n2 * n_max; ++i)
    evec2[i] = 0.0;
  for (int i = 0; i < n_max && i < n2; ++i)
    evec2[i + i * n2] = 1.0;
  
  ok = false;
  
  caslr_eff_driver_c(
    true, n, n2, n_targ, n_max, max_iter, tol, max_dav,
    apbmul_c, ambmul_c, spdmul_c, smdmul_c, lrprec_c,
    eig, evec2, &ok
  );
  
  if (ok) {
    printf("CASLR-EFF converged.\nEigenvalues:\n");
    for (int i = 0; i < n_targ; ++i)
      printf("  %.8f\n", eig[i]);
  } else {
    printf("CASLR-EFF failed to converge.\n");
  }

free(evec2);

  return 0;
}
