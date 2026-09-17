"""
Test of the python interface of DiagLib.

Every driver is run from a simple guess (unit vectors) and from a random guess (zero vectors,
so that DiagLib generates the guess): the converged eigenvalues are compared with the ones
computed by numpy, and the true residuals of the eigenvectors are computed.
The handling of invalid arrays, of DiagLib errors and of exceptions raised in the callbacks
is also checked.

usage: test_python.py [path of libdiaglib_c]
"""
import sys

import numpy as np

import pyDiaglib as dgl

N, N_TARG, N_MAX = 300, 4, 8
TOL = 1e-10
EIG_THRESH, RES_THRESH = 1e-8, 1e-7

n_tests = 0
n_failed = 0

# ------------------------------------------------------------------ test matrices

idx = np.arange(N)
I, K = np.meshgrid(idx, idx, indexing="ij")
# symmetric matrix and metric
A = 1.0 / (I + K + 2.0)
A[idx, idx] = idx + 2.0
S = 1.0 / (I + K + 2.0)
S[idx, idx] = 1.0
# non-symmetric matrix, similar to A
AR = ((I + 1.0) / (K + 1.0)) / (I + K + 2.0)
AR[idx, idx] = idx + 2.0
# linear response matrices
APB = 1.0 / (I + K + 2.0)
APB[idx, idx] = 6.0 + idx
AMB = 0.2 / (I + K + 2.0)
AMB[idx, idx] = 3.0 + idx
SPD = np.where(I > K, 0.05, np.where(I < K, -0.05, 1.0))
SMD = SPD.T.copy()


def as_matrix(ptr, n, m):
    return np.ctypeslib.as_array(ptr, shape=(n * m,)).reshape((n, m), order="F")


def make_matvec(matrix):
    def matvec(n_ptr, m_ptr, x_ptr, y_ptr):
        n, m = n_ptr[0], m_ptr[0]
        as_matrix(y_ptr, n, m)[:, :] = matrix @ as_matrix(x_ptr, n, m)
    return matvec


def precnd(n_ptr, m_ptr, shift_ptr, r_ptr, z_ptr):
    n, m = n_ptr[0], m_ptr[0]
    diag = (idx + 2.0 + shift_ptr[0])[:, None]
    as_matrix(z_ptr, n, m)[:, :] = as_matrix(r_ptr, n, m) / diag


def lrprec(n_ptr, m_ptr, fac_ptr, xp_ptr, xm_ptr, yp_ptr, ym_ptr):
    n, m = n_ptr[0], m_ptr[0]
    f = fac_ptr[0]
    d = (idx + 8.0)[:, None]
    xp, xm = as_matrix(xp_ptr, n, m), as_matrix(xm_ptr, n, m)
    denom = f * f * d * d - 1.0
    as_matrix(yp_ptr, n, m)[:, :] = (f * d * xp + xm) / denom
    as_matrix(ym_ptr, n, m)[:, :] = (f * d * xm + xp) / denom


# ------------------------------------------------------------------ references

ref_sym = np.linalg.eigvalsh(A)[:N_TARG]
L = np.linalg.cholesky(S)
Linv = np.linalg.inv(L)
ref_gen = np.linalg.eigvalsh(Linv @ A @ Linv.T)[:N_TARG]
# E z = omega M z, with E = [A B; B A], M = [S D; -D -S]
AA, BB = (APB + AMB) / 2, (APB - AMB) / 2
SS, DD = (SPD + SMD) / 2, (SPD - SMD) / 2
E = np.block([[AA, BB], [BB, AA]])
M = np.block([[SS, DD], [-DD, -SS]])
w = np.linalg.eigvals(np.linalg.solve(M, E))
ref_lr = np.sort(w.real[(abs(w.imag) < 1e-8) & (w.real > 0)])[:N_TARG]


# ------------------------------------------------------------------ checks

def record(label, passed, details=""):
    global n_tests, n_failed
    n_tests += 1
    if not passed:
        n_failed += 1
    print(f"{label:60s} {details} -> {'PASSED' if passed else 'FAILED'}")


def check(label, ok, eig, ref, residuals):
    err = np.max(np.abs(eig[:N_TARG] - ref))
    res = max(residuals)
    record(label, ok and err < EIG_THRESH and res < RES_THRESH,
           f"ok={ok} max eigenvalue error={err:.1e} max residual={res:.1e}")


def residual(matrix, x, lam, metric=None):
    bx = x if metric is None else metric @ x
    return np.linalg.norm(matrix @ x - lam * bx) / np.linalg.norm(x)


def guess(ld, random_guess):
    evec = np.zeros((ld, N_MAX), order="F")
    if not random_guess:
        evec[:N_MAX, :N_MAX] = np.eye(N_MAX)
    return np.zeros(N_MAX), evec


def expect_exception(label, exc_type, func):
    try:
        func()
    except exc_type as exc:
        record(label, True, f"{exc_type.__name__}: {exc}")
    except Exception as exc:
        record(label, False, f"unexpected {type(exc).__name__}: {exc}")
    else:
        record(label, False, "no exception")


if __name__ == "__main__":

    lib_path = sys.argv[1] if len(sys.argv) > 1 else None
    calc = dgl.diaglib(lib_path, N, N_TARG, N_MAX, tol=TOL, max_iter=200)
    print(f"DiagLib integer kind: {calc.integer_kind}")

    for random_guess in (False, True):
        kind = "random guess" if random_guess else "simple guess"

        eig, evec = guess(N, random_guess)
        ok = calc.dgl_davidson_driver(eig, evec, make_matvec(A), precnd)
        check(f"Davidson, {kind}", ok, eig, ref_sym, [residual(A, evec[:, i], eig[i]) for i in range(N_TARG)])

        eig, evec = guess(N, random_guess)
        ok = calc.dgl_davidson_driver(eig, evec, make_matvec(A), precnd, metvec=make_matvec(S))
        check(f"generalized Davidson, {kind}", ok, eig, ref_gen,
              [residual(A, evec[:, i], eig[i], S) for i in range(N_TARG)])

        eig, evec = guess(N, random_guess)
        ok = calc.dgl_lobpcg_driver(eig, evec, make_matvec(A), precnd)
        check(f"LOBPCG, {kind}", ok, eig, ref_sym, [residual(A, evec[:, i], eig[i]) for i in range(N_TARG)])

        eig, evec = guess(N, random_guess)
        ok = calc.dgl_lobpcg_driver(eig, evec, make_matvec(A), precnd, metvec=make_matvec(S))
        check(f"generalized LOBPCG, {kind}", ok, eig, ref_gen,
              [residual(A, evec[:, i], eig[i], S) for i in range(N_TARG)])

        for side in ("R", "L", "LR"):
            eig, evec1 = guess(N, random_guess)
            _, evec2 = guess(N, random_guess)
            ok = calc.dgl_davidson_nosym_driver(eig, evec1, evec2 if side == "LR" else None, side,
                                                make_matvec(AR), make_matvec(AR.T), precnd)
            first = AR.T if side == "L" else AR
            res = [residual(first, evec1[:, i], eig[i]) for i in range(N_TARG)]
            if side == "LR":
                res += [residual(AR.T, evec2[:, i], eig[i]) for i in range(N_TARG)]
                # biorthonormality error, added to the residuals
                res.append(np.max(np.abs(evec2[:, :N_TARG].T @ evec1[:, :N_TARG] - np.eye(N_TARG))))
            check(f"non-symmetric Davidson ({side}), {kind}", ok, eig, ref_sym, res)

        eig, evec = guess(2 * N, random_guess)
        ok = calc.dgl_smogd_driver(eig, evec, make_matvec(APB), make_matvec(AMB),
                                   make_matvec(SPD), make_matvec(SMD), lrprec)
        check(f"SMO-GD, {kind}", ok, eig, ref_lr, [residual(E, evec[:, i], eig[i], M) for i in range(N_TARG)])

    # invalid arrays
    eig, evec = guess(N, False)
    expect_exception("error: evec with the wrong dtype", TypeError,
                     lambda: calc.dgl_davidson_driver(eig, evec.astype(np.float32), make_matvec(A), precnd))
    expect_exception("error: C-ordered evec", ValueError,
                     lambda: calc.dgl_davidson_driver(eig, np.ascontiguousarray(evec), make_matvec(A), precnd))
    expect_exception("error: evec with the wrong shape", ValueError,
                     lambda: calc.dgl_davidson_driver(eig, evec[:, :N_TARG], make_matvec(A), precnd))
    expect_exception("error: side = LR without evec2", ValueError,
                     lambda: calc.dgl_davidson_nosym_driver(eig, evec, None, "LR", make_matvec(AR),
                                                            make_matvec(AR.T), precnd))

    # errors reported by DiagLib
    bad = dgl.diaglib(lib_path, N, N_MAX + 1, N_MAX, tol=TOL)
    expect_exception("error: n_targ > n_max", RuntimeError,
                     lambda: bad.dgl_davidson_driver(eig, evec, make_matvec(A), precnd))

    # exception raised in a callback
    calls = {"n": 0}

    def failing_matvec(n_ptr, m_ptr, x_ptr, y_ptr):
        calls["n"] += 1
        if calls["n"] == 3:
            raise ValueError("error in the user matvec")
        make_matvec(A)(n_ptr, m_ptr, x_ptr, y_ptr)

    eig, evec = guess(N, False)
    expect_exception("error: exception raised in a callback", RuntimeError,
                     lambda: calc.dgl_davidson_driver(eig, evec, failing_matvec, precnd))
    # the interface is still usable afterwards
    eig, evec = guess(N, False)
    ok = calc.dgl_davidson_driver(eig, evec, make_matvec(A), precnd)
    check("Davidson after a failed call", ok, eig, ref_sym, [residual(A, evec[:, i], eig[i]) for i in range(N_TARG)])

    print(f"Summary: {n_failed} failed tests out of {n_tests}.")
    sys.exit(1 if n_failed else 0)
