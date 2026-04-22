#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include "diaglib.h"

//
// simple minded test driver for diaglib called from c.
//
// list of functions:
// ==================
//
// matvec_c(int* n, int*m, double*x, double* ax)
//
//   applies the symmetric matrix 
//
//     a_{ij} = \delta_{ij} (i+1) + (1-\delta_{ij})/(i+j) 
//
//   to m vectors of length n (in x, assumed column major
//   as in fortran)
//
//
// precnd_c(int* n, int*m, double* shift, double*x, double* ax)
//
//   applies the preconditioner
//
//     p_{ij} = \delta_{ij} (a_{ii} + shift)
//
//   to m vectors of length n (in x, assumed column major
//   as in fortran), where a is the matrix from matvec_c.
//
// matvec_r_c(int* n, int* m, double* x, double* y)
// matvec_l_c(int* n, int* m, double* x, double* y)
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
// metvec_c(int* n, int* m, double* x, double* bx)
//
//   applies a positive definite matrix b (here, just the identity)
//   to m vectors of length n (in x, assumed column major
//   as in fortran)
// 
// apbmul_c(int* n, int* m, double* x, double* ax)
// ambmul_c(int* n, int* m, double* x, double* ax)
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
// spdmul_c(int* n, int* m, double* x, double* ax)
// smdmul_c(int* n, int* m, double* x, double* ax)
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
// lrprec_c(int* n, int* m, double* fac, double* xp, double* xm,  double* yp, double* ym)
//
//   applies the smogd diagonal preconditioner to m vectors xp and xm (note that the two
//   quantities are used to define both yp and ym), which are assumed to be column major
//   as in fortran. 
//   this routine implements the same definitions of a and sigma as in the previous 
//   routines, and can be used to test smogd. 
//
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

void matvec_l_c(int* n, int* m, double* x, double* y) {
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

void metvec_c(int* n, int* m, double* x, double* bx) {
  int N = *n;
  int M = *m;
  for (int j = 0; j < M; ++j)
    for (int i = 0; i < N; ++i)
      bx[i + j * N] = x[i + j * N];  
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

void smdmul_c(int* n, int* m, double* x, double* y) {
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

//
// small function to output the results to a file:
//
void write_results_c(const char* label, int n, int n_targ, double* eig, double* evec, const char* filename) {
  FILE* f = fopen(filename, "a");
  if (!f) return;

  fprintf(f, "  %s results\n\n", label);
  fprintf(f, "  Eigenvalues:\n");
  for (int i = 0; i < n_targ; ++i)
    fprintf(f, "  %5d%14.6f\n", i + 1, eig[i]);

  for (int i = 0; i < n_targ; ++i) {
    fprintf(f, "   Eigenvector %3d:\n", i + 1);
    for (int j = 0; j < n; ++j)
      fprintf(f, "  %5d%14.6f\n", j + 1, evec[j + i * n]);
  }

  fprintf(f, "\n");
  fclose(f);
}

void write_results_c_1(const char* label, int n, int n_targ, double* eig, const char* filename) {
  FILE* f = fopen(filename, "a");
  if (!f) return;

  fprintf(f, "  %s results\n\n", label);
  fprintf(f, "  Eigenvalues:\n");
  for (int i = 0; i < n_targ; ++i)
    fprintf(f, "  %5d%14.6f\n", i + 1, eig[i]);

  fprintf(f, "\n");
  fclose(f);
}

void write_results_c_2(const char* label, int n, int n_targ, double* eig, double* evec, 
                       double* evec_l, const char* filename) {
  FILE* f = fopen(filename, "a");
  if (!f) return;

  fprintf(f, "  %s results\n\n", label);
  fprintf(f, "  Eigenvalues:\n");
  for (int i = 0; i < n_targ; ++i)
    fprintf(f, "  %5d%14.6f\n", i + 1, eig[i]);

  for (int i = 0; i < n_targ; ++i) {
    fprintf(f, "  right Eigenvector %3d:\n", i + 1);
    for (int j = 0; j < n; ++j)
      fprintf(f, "  %5d%12.4f\n", j + 1, evec[j + i * n]);
  }
  for (int i = 0; i < n_targ; ++i) {
    fprintf(f, "  left Eigenvector %3d:\n", i + 1);
    for (int j = 0; j < n; ++j)
      fprintf(f, "  %5d%12.4f\n", j + 1, evec_l[j + i * n]);
  }

  fprintf(f, "\n");
  fclose(f);
}

void fix_phase(int n, int n_targ, double* evec) {
  for (int i = 0; i < n_targ; ++i) {
    int idx = i * n;  // indice del primo elemento della colonna i
    if (evec[idx] < 0.0) {
      for (int j = 0; j < n; ++j)
        evec[j + i * n] = -evec[j + i * n];
    }
  }
}


void test_davidson(){
#ifdef DGL_INT_KIND_4
  const int n = 500, n_targ = 5, n_max = 10, max_iter = 100, dav_iter = 20;
  const int memory = 1;
#elif DGL_INT_KIND_8
  const long int n = 500, n_targ = 5, n_max = 10, max_iter = 100, dav_iter = 20;
  const long int memory = 1;
#endif
  const double tol = 1e-10, shift = 0.0;
  double eig[n_max];
  double evec[n * n_max];
  bool ok;
  bool verbose = false;
  const char* memory_unit = "GB";
//
  printf("\nCalling DAVIDSON driver...\n");
  for (int i = 0; i < n * n_max; ++i){
    evec[i] = 0.0;
  }
  for (int j = 0; j < n_max; ++j){
    for (int i = 0; i < n; ++i){
      evec[i + j * n] = (j == i) ? 1.0 : 0.0;
    }
  }
  dgl_davidson_driver_c(n, n_targ, n_max, matvec_c, precnd_c, NULL, eig, evec, &ok,
                   verbose, tol, max_iter, dav_iter, shift, memory, memory_unit);
  
  fix_phase(n,n_targ,evec);

  if (ok) {
    printf("Davidson converged.\n");
    write_results_c("Davidson", n, n_targ, eig, evec, "output_c.txt");
  } else {
    printf("Davidson failed to converge.\n");
  }

}

void test_davidson_generalized(){
#ifdef DGL_INT_KIND_4
  const int n = 500, n_targ = 5, n_max = 10, max_iter = 100, dav_iter = 20;
  const int memory = 1;
#elif DGL_INT_KIND_8
  const long int n = 500, n_targ = 5, n_max = 10, max_iter = 100, dav_iter = 20;
  const long int memory = 1;
#endif
  const double tol = 1e-10, shift = 0.0;
  double eig[n_max];
  double evec[n * n_max];
  bool ok;
  bool verbose = false;
  const char* memory_unit = "GB";
//
  printf("\nCalling GENERALIZED DAVIDSON driver...\n");
  for (int i = 0; i < n * n_max; ++i){
    evec[i] = 0.0;
  }
  for (int j = 0; j < n_max; ++j){
    for (int i = 0; i < n; ++i){
      evec[i + j * n] = (j == i) ? 1.0 : 0.0;
    }
  }
  dgl_davidson_driver_c(n, n_targ, n_max, matvec_c, precnd_c, metvec_c, eig, evec, &ok,
                   verbose, tol, max_iter, dav_iter, shift, memory, memory_unit);
  
  fix_phase(n,n_targ,evec);

  if (ok) {
    printf("Davidson converged.\n");
    write_results_c("Davidson", n, n_targ, eig, evec, "output_c.txt");
  } else {
    printf("Davidson failed to converge.\n");
  }

}

void test_lobpcg(){
  #ifdef DGL_INT_KIND_4
  const int n = 500, n_targ = 5, n_max = 10, max_iter = 100;
  const int memory = 1;
#elif DGL_INT_KIND_8
  const long int n = 500, n_targ = 5, n_max = 10, max_iter = 100;
  const long int memory = 1;
#endif
  const double tol = 1e-10, shift = 0.0;
  double eig[n_max];
  double evec[n * n_max];
  bool ok;
  bool verbose = false;
  const char* memory_unit = "GB";
  
  printf("\nCalling LOBPCG driver...\n");
  
  for (int i = 0; i < n * n_max; ++i)
    evec[i] = 0.0;
  for (int j = 0; j < n_max; ++j)
    for (int i = 0; i < n; ++i)
      evec[i + j * n] = (j == i) ? 1.0 : 0.0;
  
  ok = false;
  
  dgl_lobpcg_driver_c(n, n_targ, n_max, matvec_c, precnd_c, NULL, eig, evec, &ok,
                   verbose, tol, max_iter, shift, memory, memory_unit);

  fix_phase(n,n_targ,evec);
  
  if (ok) {
    printf("LOBPCG converged.\n");
    write_results_c("LOBPCG", n, n_targ, eig, evec, "output_c.txt");
  } else {
    printf("LOBPCG failed to converge.\n");
  }

}

void test_lobpcg_generalized(){
  #ifdef DGL_INT_KIND_4
  const int n = 500, n_targ = 5, n_max = 10, max_iter = 100;
  const int memory = 1;
#elif DGL_INT_KIND_8
  const long int n = 500, n_targ = 5, n_max = 10, max_iter = 100;
  const long int memory = 1;
#endif
  const double tol = 1e-10, shift = 0.0;
  double eig[n_max];
  double evec[n * n_max];
  bool ok;
  bool verbose = false;
  const char* memory_unit = "GB";
  
  printf("\nCalling GENERALIZED LOBPCG driver...\n");
  
  for (int i = 0; i < n * n_max; ++i)
    evec[i] = 0.0;
  for (int j = 0; j < n_max; ++j)
    for (int i = 0; i < n; ++i)
      evec[i + j * n] = (j == i) ? 1.0 : 0.0;
  
  ok = false;
  
  dgl_lobpcg_driver_c(n, n_targ, n_max, matvec_c, precnd_c, metvec_c, eig, evec, &ok,
                   verbose, tol, max_iter, shift, memory, memory_unit);

  fix_phase(n,n_targ,evec);
  
  if (ok) {
    printf("LOBPCG converged.\n");
    write_results_c("LOBPCG", n, n_targ, eig, evec, "output_c.txt");
  } else {
    printf("LOBPCG failed to converge.\n");
  }

}

void test_davidson_nosym_davidson(){
  #ifdef DGL_INT_KIND_4
  const int n = 500, n_targ = 5, n_max = 10, max_iter = 100, dav_iter = 20;
  const int memory = 1;
#elif DGL_INT_KIND_8
  const long int n = 500, n_targ = 5, n_max = 10, max_iter = 100, dav_iter = 20;
  const long int memory = 1;
#endif
  const double tol = 1e-10, shift = 0.0;
  double eig[n_max];
  double evec[n * n_max], evec_l[n * n_max];
  bool ok;
  bool verbose = false;
  const char* memory_unit = "GB";
  const char* side = "LR";

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
  
  dgl_davidson_nosym_driver_c(n, n_targ, n_max, matvec_r_c, matvec_l_c, precnd_c, side,
                              eig, evec, evec_l, &ok,
                              verbose, tol, max_iter, dav_iter, shift, memory, memory_unit);
  fix_phase(n,n_targ,evec);
  fix_phase(n,n_targ,evec_l);
  
  if (ok) {
    printf("Non-symmetric Davidson converged.\n");
    write_results_c_2("Non-Symmetric Davidson", n, n_targ, eig, evec, evec_l, "output_c.txt");
  } else {
    printf("Non-symmetric Davidson failed to converge.\n");
  }

}

void test_smogd(){
  #ifdef DGL_INT_KIND_4
  const int n = 500, n_targ = 5, n_max = 10, max_iter = 100, dav_iter = 20;
  const int memory = 1;
  int n2 = 2 * n;
#elif DGL_INT_KIND_8
  const long int n = 500, n_targ = 5, n_max = 10, max_iter = 100, dav_iter = 20;
  const long int memory = 1;
  long int n2 = 2 * n;
#endif
  const double tol = 1e-10;
  double eig[n_max];
  double evec2 [n2 * n_max];
  
  bool ok;
  bool verbose = true;
  const char* memory_unit = "GB";

  printf("\nCalling SMOGD driver...\n");
  for (int i = 0; i < n2 * n_max; ++i)
    evec2[i] = 0.0;
  for (int i = 0; i < n_max && i < n2; ++i)
    evec2[i + i * n2] = 1.0;
  
  ok = false;
  
  dgl_smogd_driver_c(n2, n_targ, n_max, apbmul_c, ambmul_c, spdmul_c, smdmul_c, 
                 lrprec_c, eig, evec2, &ok,
                 verbose, tol, max_iter, dav_iter, memory, memory_unit);
  
  if (ok) {
    printf("SMOGD converged.\n");
    write_results_c_1("SMOGD", n2, n_targ, eig, "output_c.txt");
  } else {
    printf("SMOGD failed to converge.\n");
  }

}
//
// main program: test davidson, non-symmetric davidson, lobpcg and smogd.
//
int main() {
//
// get rid of the output file if it's already present.
//
  remove("output_c.txt");

  test_davidson();
  test_davidson_generalized();
  test_lobpcg();
  test_lobpcg_generalized();
  test_davidson_nosym_davidson();
  test_smogd();

  return 0;
}
