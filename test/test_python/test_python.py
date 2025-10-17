import ctypes
import numpy as np

# Carica la libreria
lib = ctypes.CDLL('./../../lib/libdiaglib.so')

# Definisci i tipi delle callback
MATVEC = ctypes.CFUNCTYPE(None,
    ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double)
)

MATVEC2 = MATVEC

PRECND = ctypes.CFUNCTYPE(None,
    ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double)
)

BVEC = ctypes.CFUNCTYPE(None,
    ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double)
)

LRPREC = ctypes.CFUNCTYPE(None,
    ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double)
)

# Definisci la funzione Fortran
lib.davidson_driver_c.argtypes = [
    ctypes.c_bool, ctypes.c_int, ctypes.c_int, ctypes.c_int,
    ctypes.c_int, ctypes.c_int, ctypes.c_double, ctypes.c_double,
    MATVEC, PRECND,
    ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_bool)
]

lib.lobpcg_driver_c.argtypes = [
    ctypes.c_bool, ctypes.c_bool,
    ctypes.c_int, ctypes.c_int, ctypes.c_int, ctypes.c_int,
    ctypes.c_double, ctypes.c_double,
    MATVEC, PRECND, BVEC,
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_bool)
]

lib.nonsym_driver_c.argtypes = [
    ctypes.c_bool, ctypes.c_int, ctypes.c_int, ctypes.c_int, ctypes.c_int,
    ctypes.c_double, ctypes.c_int, ctypes.c_double,
    MATVEC2, MATVEC2, PRECND,
    ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double),
    ctypes.c_int, ctypes.POINTER(ctypes.c_bool)
]

lib.smogd_driver_c.argtypes = [
    ctypes.c_bool, ctypes.c_int, ctypes.c_int, ctypes.c_int, ctypes.c_int, ctypes.c_int,
    ctypes.c_double, ctypes.c_int,
    MATVEC2, MATVEC2, MATVEC2, MATVEC2, LRPREC,
    ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_bool)
]

# Matvec: A_{ii} = i+1, A_{ij} = 1/(i+j)
def matvec(n_ptr, m_ptr, x_ptr, ax_ptr):
    n = n_ptr[0]
    m = m_ptr[0]
    x = np.ctypeslib.as_array(x_ptr, shape=(n*m,))
    ax = np.ctypeslib.as_array(ax_ptr, shape=(n*m,))
    for j in range(m):
        for i in range(n):
            ax[i + j*n] = sum(
                ((i + 2.0) if i == k else 1.0 / (i + k + 2.0)) * x[k + j*n]
                for k in range(n)
            )

# Precondizionatore: solo diagonale
def precnd(n_ptr, m_ptr, shift_ptr, r_ptr, z_ptr):
    n = n_ptr[0]
    m = m_ptr[0]
    shift = shift_ptr[0]
    r = np.ctypeslib.as_array(r_ptr, shape=(n*m,))
    z = np.ctypeslib.as_array(z_ptr, shape=(n*m,))
    for j in range(m):
        for i in range(n):
            diag = (i + 2.0) + shift
            z[i + j*n] = r[i + j*n] / diag

def bvec(n_ptr, m_ptr, x_ptr, bx_ptr):
    n = n_ptr[0]
    m = m_ptr[0]
    x = np.ctypeslib.as_array(x_ptr, shape=(n*m,))
    bx = np.ctypeslib.as_array(bx_ptr, shape=(n*m,))
    for j in range(m):
        for i in range(n):
            bx[i + j*n] = x[i + j*n]

def matvec_r(n_ptr, m_ptr, x_ptr, y_ptr):
    n, m = n_ptr[0], m_ptr[0]
    x = np.ctypeslib.as_array(x_ptr, shape=(n*m,))
    y = np.ctypeslib.as_array(y_ptr, shape=(n*m,))
    for k in range(m):
        for i in range(n):
            y[i + k*n] = sum(
                ((i + 2.0) if j == i else ((i + 1.0)/(j + 1.0)) / (i + j + 2.0)) * x[j + k*n]
                for j in range(n)
            )

def matvec_l(n_ptr, m_ptr, x_ptr, y_ptr):
    n, m = n_ptr[0], m_ptr[0]
    x = np.ctypeslib.as_array(x_ptr, shape=(n*m,))
    y = np.ctypeslib.as_array(y_ptr, shape=(n*m,))
    for k in range(m):
        for i in range(n):
            y[i + k*n] = sum(
                ((i + 2.0) if j == i else ((j + 1.0)/(i + 1.0)) / (i + j + 2.0)) * x[j + k*n]
                for j in range(n)
            )

def apbmul(n_ptr, m_ptr, x_ptr, ax_ptr):
    n, m = n_ptr[0], m_ptr[0]
    x = np.ctypeslib.as_array(x_ptr, shape=(n*m,))
    ax = np.ctypeslib.as_array(ax_ptr, shape=(n*m,))
    for j in range(m):
        for i in range(n):
            ax[i + j*n] = sum(
                ((6.0 + i) if i == k else 1.0 / (i + k + 2.0)) * x[k + j*n]
                for k in range(n)
            )

def ambmul(n_ptr, m_ptr, x_ptr, ax_ptr):
    n, m = n_ptr[0], m_ptr[0]
    x = np.ctypeslib.as_array(x_ptr, shape=(n*m,))
    ax = np.ctypeslib.as_array(ax_ptr, shape=(n*m,))
    for j in range(m):
        for i in range(n):
            ax[i + j*n] = sum(
                ((3.0 + i) if i == k else 0.2 / (i + k + 2.0)) * x[k + j*n]
                for k in range(n)
            )

def spdmul(n_ptr, m_ptr, x_ptr, y_ptr):
    n, m = n_ptr[0], m_ptr[0]
    x = np.ctypeslib.as_array(x_ptr, shape=(n*m,))
    y = np.ctypeslib.as_array(y_ptr, shape=(n*m,))
    for k in range(m):
        for i in range(n):
            y[i + k*n] = sum(
                x[j + k*n] if i == j else (0.05 * x[j + k*n] if i > j else -0.05 * x[j + k*n])
                for j in range(n)
            )

def smdmul(n_ptr, m_ptr, x_ptr, y_ptr):
    n, m = n_ptr[0], m_ptr[0]
    x = np.ctypeslib.as_array(x_ptr, shape=(n*m,))
    y = np.ctypeslib.as_array(y_ptr, shape=(n*m,))
    for k in range(m):
        for i in range(n):
            y[i + k*n] = sum(
                x[j + k*n] if i == j else (-0.05 * x[j + k*n] if i > j else 0.05 * x[j + k*n])
                for j in range(n)
            )

def lrprec(n_ptr, m_ptr, fac_ptr, xp_ptr, xm_ptr, yp_ptr, ym_ptr):
    n, m = n_ptr[0], m_ptr[0]
    f = fac_ptr[0]
    xp = np.ctypeslib.as_array(xp_ptr, shape=(n*m,))
    xm = np.ctypeslib.as_array(xm_ptr, shape=(n*m,))
    yp = np.ctypeslib.as_array(yp_ptr, shape=(n*m,))
    ym = np.ctypeslib.as_array(ym_ptr, shape=(n*m,))
    for j in range(m):
        for i in range(n):
            denom = f*f*(i+8.0)*(i+8.0) - 1.0
            yp[i + j*n] = (f*(i+8.0)*xp[i + j*n] + xm[i + j*n]) / denom
            ym[i + j*n] = (f*(i+8.0)*xm[i + j*n] + xp[i + j*n]) / denom

# Parametri
n = 100
n_targ = 2
n_max = 4
max_iter = 100
max_dav = 20
tol = 1e-8
shift = 0.0

eig = np.zeros(n_max, dtype=np.float64)
evec = np.zeros((n, n_max), dtype=np.float64, order='F')
for i in range(min(n, n_max)):
    evec[i, i] = 1.0

ok = ctypes.c_bool(False)

# Chiamata
lib.davidson_driver_c(
    True, n, n_targ, n_max, max_iter, max_dav, tol, shift,
    MATVEC(matvec), PRECND(precnd),
    eig.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
    evec.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
    ctypes.byref(ok)
)

# Risultati
if ok.value:
    print("Davidson converged.")
    print("Eigenvalues:", eig[:n_targ])
else:
    print("Davidson failed.")



evec = np.zeros((n, n_max), dtype=np.float64, order='F')
for i in range(min(n, n_max)):
    evec[i, i] = 1.0
ok = ctypes.c_bool(False)

lib.lobpcg_driver_c(
    True, False,  # verbose, gen_eig
    n, n_targ, n_max, max_iter, tol, shift,
    MATVEC(matvec), PRECND(precnd), BVEC(bvec),
    eig.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
    evec.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
    ctypes.byref(ok)
)

if ok.value:
    print("LOBPCG converged.")
    print("Eigenvalues:", eig[:n_targ])
else:
    print("LOBPCG failed.")

evec_l = np.zeros((n, n_max), dtype=np.float64, order='F')
for i in range(min(n, n_max)):
    evec[i, i] = 1.0
    evec_l[i, i] = 1.0
ok = ctypes.c_bool(False)

lib.nonsym_driver_c(
    True, n, n_targ, n_max, max_iter, tol, max_dav, shift,
    MATVEC2(matvec_r), MATVEC2(matvec_l), PRECND(precnd),
    eig.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
    evec.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
    evec_l.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
    4, ctypes.byref(ok)
)

if ok.value:
    print("Non-symmetric Davidson converged.")
    print("Eigenvalues:", eig[:n_targ])
else:
    print("Non-symmetric Davidson failed.")

# === Chiamata smogd_driver_c ===

n2 = 2 * n
evec2 = np.zeros((n2, n_max), dtype=np.float64, order='F')
for i in range(min(n2, n_max)):
    evec2[i, i] = 1.0
ok = ctypes.c_bool(False)

lib.smogd_driver_c(
    True, n, n2, n_targ, n_max, max_iter, tol, max_dav,
    MATVEC2(apbmul), MATVEC2(ambmul), MATVEC2(spdmul), MATVEC2(smdmul), LRPREC(lrprec),
    eig.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
    evec2.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
    ctypes.byref(ok)
)

if ok.value:
    print("SMOGD converged.")
    print("Eigenvalues:", eig[:n_targ])
else:
    print("SMOGD failed.")
