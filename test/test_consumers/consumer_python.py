"""Python program using an installed pyDiaglib (found through PYTHONPATH, with the library through
DIAGLIB_C_LIBRARY or its path as the first argument). The test problems have known eigenvalues
and are applied matrix-free with numpy."""
import sys
import threading

import numpy as np
import pyDiaglib as dgl

lib_path = sys.argv[1] if len(sys.argv) > 1 else None
N, N_TARG, N_MAX = 1000, 4, 8
i = np.arange(1, N + 1, dtype=float)
v = np.sin(i); v /= np.linalg.norm(v)
u = 0.3 * np.cos(2 * i) / np.sqrt(N)
w = 0.5 * np.sin(3 * i) / np.sqrt(N)
d = 1.0 + (i - 1) + 0.5 * np.sin(i) ** 2
s = 1.0 + 0.2 * np.cos(i) ** 2
a = 2.0 + i
b = 1.0 + 0.5 * i
failures = 0


def mat(ptr, n, m):
    return np.ctypeslib.as_array(ptr, shape=(n * m,)).reshape((n, m), order="F")


def hdh(diag, x):
    y = x - 2 * np.outer(v, v @ x)
    y = diag[:, None] * y
    return y - 2 * np.outer(v, v @ y)


def make_mv(diag, scale=1.0):
    def mv(n, m, x, y):
        mat(y, n[0], m[0])[:, :] = scale * hdh(diag, mat(x, n[0], m[0]))
    return mv


def make_pc(diag, spd=False):
    def pc(n, m, shift, x, y):
        den = diag + (0.0 if spd else shift[0])
        den = np.where(np.abs(den) < 1e-3, np.copysign(1e-3, den), den)
        mat(y, n[0], m[0])[:, :] = mat(x, n[0], m[0]) / den[:, None]
    return pc


def mv_r(n, m, x, y):
    X = mat(x, n[0], m[0])
    Y = X - np.outer(u, w @ X) / (1 + w @ u)
    Y = d[:, None] * Y
    mat(y, n[0], m[0])[:, :] = Y + np.outer(u, w @ Y)


def mv_l(n, m, x, y):
    X = mat(x, n[0], m[0])
    Y = d[:, None] * (X + np.outer(w, u @ X))
    mat(y, n[0], m[0])[:, :] = Y - np.outer(w, u @ Y) / (1 + w @ u)


def lrprec(n, m, fac, xp, xm, yp, ym):
    f = fac[0]
    den = f * f * a * b - 1.0
    den = np.where(np.abs(den) < 1e-3, 1e-3, den)[:, None]
    XP, XM = mat(xp, n[0], m[0]), mat(xm, n[0], m[0])
    mat(yp, n[0], m[0])[:, :] = (f * b[:, None] * XP + XM) / den
    mat(ym, n[0], m[0])[:, :] = (f * a[:, None] * XM + XP) / den


def check(label, cond, extra=""):
    global failures
    failures += not cond
    print(f"{label:45s} {'PASSED' if cond else 'FAILED'} {extra}")


def lowest(x, k=N_TARG):
    return np.sort(x)[:k]


def zeros(rows=N):
    return np.zeros(N_MAX), np.zeros((rows, N_MAX), order="F")


calc = dgl.diaglib(lib_path, N, N_TARG, N_MAX, tol=1e-9, max_iter=300)
print("integer kind of the library:", calc.integer_kind)

eig, evec = zeros()
ok = calc.dgl_davidson_driver(eig, evec, make_mv(d), make_pc(d))
err = np.abs(eig[:N_TARG] - lowest(d)).max()
res = np.linalg.norm(hdh(d, evec[:, :N_TARG]) - evec[:, :N_TARG] * eig[:N_TARG], axis=0).max()
check("Davidson", ok and err < 1e-8 and res < 1e-6, f"err {err:.1e} res {res:.1e}")

eig, evec = zeros()
ok = calc.dgl_lobpcg_driver(eig, evec, make_mv(d), make_pc(d / s, spd=True), metvec=make_mv(s))
err = np.abs(eig[:N_TARG] - lowest(d / s)).max()
check("generalized LOBPCG", ok and err < 1e-8, f"err {err:.1e}")

eig, evec1 = zeros(); _, evec2 = zeros()
ok = calc.dgl_davidson_nosym_driver(eig, evec1, evec2, "LR", mv_r, mv_l, make_pc(d))
err = np.abs(eig[:N_TARG] - lowest(d)).max()
bi = np.abs(evec2[:, :N_TARG].T @ evec1[:, :N_TARG] - np.eye(N_TARG)).max()
check("non-symmetric LR", ok and err < 1e-8 and bi < 1e-10, f"err {err:.1e} biortho {bi:.1e}")

eig, evec = zeros(2 * N)
ident = lambda n, m, x, y: mat(y, n[0], m[0]).__setitem__(slice(None), mat(x, n[0], m[0]))
ok = calc.dgl_smogd_driver(eig, evec, make_mv(a), make_mv(b), ident, ident, lrprec)
err = np.abs(eig[:N_TARG] - lowest(np.sqrt(a * b))).max()
check("SMO-GD", ok and err < 1e-8, f"err {err:.1e}")

# errors
eig, evec = zeros()
for label, exc, call in [
    ("error: C-ordered evec", ValueError, lambda: calc.dgl_davidson_driver(eig, np.zeros((N, N_MAX)), make_mv(d), make_pc(d))),
    ("error: int32 eig", TypeError, lambda: calc.dgl_davidson_driver(eig.astype(np.int32), evec, make_mv(d), make_pc(d))),
    ("error: invalid side", ValueError, lambda: calc.dgl_davidson_nosym_driver(eig, evec, None, "RL", mv_r, mv_l, make_pc(d))),
    ("error: DiagLib input error", RuntimeError,
     lambda: dgl.diaglib(lib_path, N, N_TARG, N_MAX, tol=-1.0).dgl_davidson_driver(eig, evec, make_mv(d), make_pc(d))),
    ("error: exception in callback", RuntimeError,
     lambda: calc.dgl_davidson_driver(eig, evec, lambda *args: 1 / 0, make_pc(d))),
]:
    try:
        call()
        check(label, False, "no exception")
    except exc as e:
        check(label, True, f"{type(e).__name__}: {e}")

# nested call, from a callback, on the same pyDiaglib object and on another one
inner = dgl.diaglib(lib_path, N, 2, 4, tol=1e-8, max_iter=400)
count = {"same": 0}


def mv_nested(n, m, x, y):
    for obj in (calc, inner):
        e2, v2 = np.zeros(obj.n_max), np.zeros((N, obj.n_max), order="F")
        okn = obj.dgl_davidson_driver(e2, v2, make_mv(d), make_pc(d))
        if not okn or abs(e2[0] - d.min()) > 1e-8:
            raise RuntimeError(f"nested call failed: {okn=} {obj is calc=} {e2[0] - d.min()=}")
    count["same"] += 1
    make_mv(d)(n, m, x, y)


eig, evec = zeros()
try:
    ok = calc.dgl_davidson_driver(eig, evec, mv_nested, make_pc(d))
    check("nested calls (same and other object)", ok and np.abs(eig[:N_TARG] - lowest(d)).max() < 1e-8 and count["same"] > 0)
except RuntimeError as e:
    check("nested calls (same and other object)", False, repr(e.__cause__))

# the object is still usable after all that
eig, evec = zeros()
check("usable afterwards", calc.dgl_davidson_driver(eig, evec, make_mv(d), make_pc(d)))

# concurrent calls from python threads (ctypes releases the GIL during the call)
results = {}


def task(t):
    obj = dgl.diaglib(lib_path, N, N_TARG, N_MAX, tol=1e-9)
    e, ev = np.zeros(N_MAX), np.zeros((N, N_MAX), order="F")
    sc = 1.0 + t
    ok = obj.dgl_davidson_driver(e, ev, make_mv(d, sc), make_pc(sc * d))
    results[t] = ok and np.abs(e[:N_TARG] - sc * lowest(d)).max() < 1e-8


threads = [threading.Thread(target=task, args=(t,)) for t in range(6)]
for th in threads:
    th.start()
for th in threads:
    th.join()
check("6 concurrent python threads", len(results) == 6 and all(results.values()))

print(f"python consumer: {failures} failures")
raise SystemExit(failures != 0)
