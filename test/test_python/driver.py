import ctypes
import numpy as np

# Carica la libreria
lib = ctypes.CDLL('./../../lib/libdriver.so')

# Definisci i tipi delle callback
MATVEC = ctypes.CFUNCTYPE(None,
    ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double)
)

PRECND = ctypes.CFUNCTYPE(None,
    ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double)
)

BVEC = ctypes.CFUNCTYPE(None,
    ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int),
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

# Parametri
n = 1000
n_targ = 20
n_max = 25
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

