import numpy as np
import pyDiaglib as dgl

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

def metvec(n_ptr, m_ptr, x_ptr, bx_ptr):
    n = n_ptr[0]
    m = m_ptr[0]
    x = np.ctypeslib.as_array(x_ptr, shape=(n*m,))
    bx = np.ctypeslib.as_array(bx_ptr, shape=(n*m,))
    for j in range(m):
        for i in range(n):
            bx[i + j*n] = sum(
                ( 1.0 if i == k else 1.0 / (i + k + 2.0)) * x[k + j*n]
                for k in range(n)
            )

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


def reset_eigs(n, n_max):
    eig = np.zeros(n_max, dtype=np.float64)
    evec = np.zeros((n, n_max), dtype=np.float64, order='F')
    for i in range(min(n, n_max)):
        evec[i, i] = 1.0
    return eig, evec

def glance_results(ok,eig,evec,string):
    if ok.value:
        print(f"\n{string} converged.")
        print("Eigenvalues:", eig[:n_targ])
        print("Eigenvectors:", evec[:5,:n_targ])
    else:
        print(f"\n{string} failed.")

if __name__ == "__main__":

    n = 500
    n_targ = 5
    n_max = 10

    mycalc = dgl.diaglib("/home/i-gianni/software/diaglib/build/src/c_interface/libdiaglib_c.so", 4,
                n, n_targ, n_max)
    

    eig, evec = reset_eigs(n,n_max)
    ok = mycalc.dgl_davidson_driver(eig, evec, matvec, precnd)
    #glance_results(ok,eig,evec,"Davidson")

    eig, evec = reset_eigs(n,n_max)
    ok = mycalc.dgl_lobpcg_driver(eig, evec, matvec, precnd)
    #glance_results(ok,eig,evec,"LOBPCG")
    
    eig, evec = reset_eigs(n,n_max)
    ok = mycalc.dgl_davidson_driver(eig, evec, matvec, precnd, metvec=metvec)
    #glance_results(ok,eig,evec,"Generalized Davidson")

    eig, evec = reset_eigs(n,n_max)
    ok = mycalc.dgl_lobpcg_driver(eig, evec, matvec, precnd, metvec=metvec)
    #glance_results(ok,eig,evec,"Generalized LOBPCG")

    eig, evec = reset_eigs(n,n_max)
    eig, evec1 = reset_eigs(n,n_max)
    ok = mycalc.dgl_davidson_nosym_driver(eig, evec, evec1, "LR", matvec_r, matvec_l, precnd)
    #glance_results(ok,eig,evec,"Non-symmetric Davidson")

    eig, evec = reset_eigs(2*n,n_max)
    ok = mycalc.dgl_smogd_driver(eig, evec, apbmul, ambmul, spdmul, smdmul, lrprec)
    #glance_results(ok,eig,evec,"SMOGD")
