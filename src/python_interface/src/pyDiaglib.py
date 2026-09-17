"""
pyDiaglib: python interface to DiagLib, based on ctypes and numpy.

The drivers work in place on numpy arrays, which must be of type float64 and Fortran-ordered
(e.g., np.zeros((n, n_max), order='F')). Eigenvectors are stored by column.
If the eigenvector array is zero, DiagLib starts from a random guess.

The user-supplied callbacks receive ctypes pointers, e.g. for a matrix-vector product:

    def matvec(n_ptr, m_ptr, x_ptr, y_ptr):
        n, m = n_ptr[0], m_ptr[0]
        x = np.ctypeslib.as_array(x_ptr, shape=(n * m,)).reshape((n, m), order='F')
        y = np.ctypeslib.as_array(y_ptr, shape=(n * m,)).reshape((n, m), order='F')
        y[:, :] = A @ x

Errors reported by DiagLib, and exceptions raised in the callbacks, are raised as RuntimeError.
"""
import ctypes
import ctypes.util
import os

import numpy as np

# Error codes returned by DiagLib (see diaglib.h)
DGL_ERRORS = {
    -1: "invalid input arguments",
    -2: "allocation failure or memory limit exceeded",
    -3: "a LAPACK routine failed",
    -4: "an orthogonalization procedure failed",
    -5: "left and right eigenvalues do not match",
}


def _load_library(lib_path):
    """Load libdiaglib_c: from lib_path, from the DIAGLIB_C_LIBRARY environment variable,
    or from the standard library search path."""
    if lib_path is None:
        lib_path = os.environ.get("DIAGLIB_C_LIBRARY") or ctypes.util.find_library("diaglib_c")
    if lib_path is None:
        raise OSError("Could not find the DiagLib C library (libdiaglib_c): pass its path, "
                      "or set the DIAGLIB_C_LIBRARY environment variable")
    try:
        return ctypes.CDLL(lib_path)
    except OSError as exc:
        raise OSError(f"Could not load the DiagLib C library from {lib_path}: {exc}") from exc


class diaglib:
    """
    Interface to the DiagLib drivers.

    lib_path: path of libdiaglib_c (optional, see _load_library)
    n:        size of the problem. For SMO-GD, size of the A, B, S, D blocks: the eigenvectors
              then have 2*n components.
    n_targ:   number of required eigenpairs
    n_max:    size of the arrays of eigenvalues and eigenvectors (n_max >= n_targ)
    The other arguments are the options of the drivers, with the same defaults as in Fortran.
    shift is only added to the eigenvalues when they are printed (verbose=True).
    """

    def __init__(self, lib_path, n, n_targ, n_max,
                 verbose=False,
                 max_iter=100,
                 dav_iter=25,
                 tol=1e-7,
                 shift=0.0,
                 memory=80,
                 memory_unit="MB"):
        self.lib = _load_library(lib_path)

        # The size of the integers is the one of the library that has been loaded
        self.lib.dgl_integer_kind.argtypes = []
        self.lib.dgl_integer_kind.restype = ctypes.c_int
        kind = self.lib.dgl_integer_kind()
        if kind == 4:
            cInt = ctypes.c_int32
        elif kind == 8:
            cInt = ctypes.c_int64
        else:
            raise RuntimeError(f"Unexpected size of the DiagLib integers: {kind}")
        self.integer_kind = kind
        self.__cInt = cInt
        cDouble = ctypes.c_double
        cBool = ctypes.c_bool
        cCharP = ctypes.c_char_p
        pDouble = ctypes.POINTER(cDouble)

        self.n = n
        self.n_targ = n_targ
        self.n_max = n_max
        self.verbose = verbose
        self.max_iter = max_iter
        self.dav_iter = dav_iter
        self.tol = tol
        self.shift = shift
        self.memory = memory
        self.memory_unit = memory_unit

        MATVEC = ctypes.CFUNCTYPE(None,
            ctypes.POINTER(cInt), ctypes.POINTER(cInt),
            pDouble, pDouble
        )

        PRECND = ctypes.CFUNCTYPE(None,
            ctypes.POINTER(cInt), ctypes.POINTER(cInt),
            pDouble, pDouble, pDouble
        )

        LRPREC = ctypes.CFUNCTYPE(None,
            ctypes.POINTER(cInt), ctypes.POINTER(cInt), pDouble,
            pDouble, pDouble, pDouble, pDouble
        )

        self.__MATVEC = MATVEC
        self.__PRECND = PRECND
        self.__LRPREC = LRPREC
        self.__NONE = ctypes.cast(None, MATVEC)
        self.__callback_error = None

        self.lib.dgl_davidson_driver.argtypes = [
            cInt, cInt, cInt, MATVEC, PRECND, MATVEC,
            pDouble, pDouble, ctypes.POINTER(cBool), ctypes.POINTER(cInt),
            cBool, cDouble, cInt, cInt, cDouble, cInt, cCharP
        ]
        self.lib.dgl_davidson_driver.restype = None

        self.lib.dgl_lobpcg_driver.argtypes = [
            cInt, cInt, cInt, MATVEC, PRECND, MATVEC,
            pDouble, pDouble, ctypes.POINTER(cBool), ctypes.POINTER(cInt),
            cBool, cDouble, cInt, cDouble, cInt, cCharP
        ]
        self.lib.dgl_lobpcg_driver.restype = None

        self.lib.dgl_davidson_nosym_driver.argtypes = [
            cInt, cInt, cInt, MATVEC, MATVEC, PRECND, cCharP,
            pDouble, pDouble, pDouble, ctypes.POINTER(cBool), ctypes.POINTER(cInt),
            cBool, cDouble, cInt, cInt, cDouble, cInt, cCharP
        ]
        self.lib.dgl_davidson_nosym_driver.restype = None

        self.lib.dgl_smogd_driver.argtypes = [
            cInt, cInt, cInt, MATVEC, MATVEC, MATVEC, MATVEC, LRPREC,
            pDouble, pDouble, ctypes.POINTER(cBool), ctypes.POINTER(cInt),
            cBool, cDouble, cInt, cInt, cInt, cCharP
        ]
        self.lib.dgl_smogd_driver.restype = None

    # ------------------------------------------------------------------ helpers

    @staticmethod
    def __array(a, shape, name):
        """Check that a is an array DiagLib can work on in place, and return a pointer to its data"""
        if not isinstance(a, np.ndarray):
            raise TypeError(f"{name} must be a numpy array")
        if a.dtype != np.float64:
            raise TypeError(f"{name} must have dtype float64, not {a.dtype}")
        if a.shape != shape:
            raise ValueError(f"{name} must have shape {shape}, not {a.shape}")
        if not a.flags.f_contiguous:
            raise ValueError(f"{name} must be a Fortran-ordered, contiguous array "
                             f"(e.g., np.zeros({shape}, order='F'))")
        if not a.flags.writeable:
            raise ValueError(f"{name} must be writeable")
        return a.ctypes.data_as(ctypes.POINTER(ctypes.c_double))

    def __callback(self, func, ctype, outputs):
        """
        Wrap a user callback: an exception cannot propagate through the Fortran code, so it is
        stored, the outputs (the arguments in positions outputs) are filled with NaNs, and the
        exception is raised again when the driver returns.
        """
        def wrapper(*args):
            if self.__callback_error is None:
                try:
                    func(*args)
                    return
                except BaseException as exc:
                    self.__callback_error = exc
            size = args[0][0] * args[1][0]
            for k in outputs:
                np.ctypeslib.as_array(args[k], shape=(size,))[:] = np.nan
        return ctype(wrapper)

    def __matvec(self, func):
        return self.__callback(func, self.__MATVEC, (3,))

    def __precnd(self, func):
        return self.__callback(func, self.__PRECND, (4,))

    def __lrprec(self, func):
        return self.__callback(func, self.__LRPREC, (5, 6))

    def __check(self, info):
        if self.__callback_error is not None:
            exc, self.__callback_error = self.__callback_error, None
            raise RuntimeError("An exception was raised in a callback of DiagLib") from exc
        if info.value != 0:
            raise RuntimeError(f"DiagLib error {info.value}: {DGL_ERRORS.get(info.value, 'unknown error')}")

    def __memory_unit(self):
        return ctypes.c_char_p(bytes(self.memory_unit, "utf-8"))

    # ------------------------------------------------------------------ drivers

    def dgl_davidson_driver(self, eig, evec, matvec, precnd, metvec=None):
        """Davidson-Liu, standard or (if metvec is given) generalized symmetric problem.
        eig: (n_max,), evec: (n, n_max). Returns True if converged."""
        p_eig = self.__array(eig, (self.n_max,), "eig")
        p_evec = self.__array(evec, (self.n, self.n_max), "evec")
        c_matvec = self.__matvec(matvec)
        c_precnd = self.__precnd(precnd)
        c_metvec = self.__NONE if metvec is None else self.__matvec(metvec)
        ok = ctypes.c_bool(False)
        info = self.__cInt(0)

        self.lib.dgl_davidson_driver(
            self.n, self.n_targ, self.n_max,
            c_matvec, c_precnd, c_metvec,
            p_eig, p_evec,
            ctypes.byref(ok), ctypes.byref(info),
            self.verbose, self.tol, self.max_iter, self.dav_iter, self.shift,
            self.memory, self.__memory_unit()
        )
        self.__check(info)
        return bool(ok.value)

    def dgl_lobpcg_driver(self, eig, evec, matvec, precnd, metvec=None):
        """LOBPCG, standard or (if metvec is given) generalized symmetric problem.
        eig: (n_max,), evec: (n, n_max). Returns True if converged."""
        p_eig = self.__array(eig, (self.n_max,), "eig")
        p_evec = self.__array(evec, (self.n, self.n_max), "evec")
        c_matvec = self.__matvec(matvec)
        c_precnd = self.__precnd(precnd)
        c_metvec = self.__NONE if metvec is None else self.__matvec(metvec)
        ok = ctypes.c_bool(False)
        info = self.__cInt(0)

        self.lib.dgl_lobpcg_driver(
            self.n, self.n_targ, self.n_max,
            c_matvec, c_precnd, c_metvec,
            p_eig, p_evec,
            ctypes.byref(ok), ctypes.byref(info),
            self.verbose, self.tol, self.max_iter, self.shift,
            self.memory, self.__memory_unit()
        )
        self.__check(info)
        return bool(ok.value)

    def dgl_davidson_nosym_driver(self, eig, evec1, evec2, side, matvec_r, matvec_l, precnd):
        """Non-symmetric Davidson. side: "R", "L" or "LR".
        eig: (n_max,); evec1: (n, n_max), right eigenvectors for "R" and "LR", left ones for "L";
        evec2: (n, n_max), left eigenvectors, only used (and required) for "LR", may be None otherwise.
        For "LR", if converged, the first n_targ vectors are biorthonormal: evec2[:, i] @ evec1[:, j] = delta_ij.
        Returns True if converged."""
        if side not in ("R", "L", "LR"):
            raise ValueError(f"side must be 'R', 'L' or 'LR', not {side!r}")
        p_eig = self.__array(eig, (self.n_max,), "eig")
        p_evec1 = self.__array(evec1, (self.n, self.n_max), "evec1")
        if side == "LR":
            if evec2 is None:
                raise ValueError("evec2 is required for side = 'LR'")
            p_evec2 = self.__array(evec2, (self.n, self.n_max), "evec2")
        else:
            p_evec2 = None
        c_matvec_r = self.__matvec(matvec_r)
        c_matvec_l = self.__matvec(matvec_l)
        c_precnd = self.__precnd(precnd)
        ok = ctypes.c_bool(False)
        info = self.__cInt(0)

        self.lib.dgl_davidson_nosym_driver(
            self.n, self.n_targ, self.n_max,
            c_matvec_r, c_matvec_l, c_precnd,
            ctypes.c_char_p(bytes(side, "utf-8")),
            p_eig, p_evec1, p_evec2,
            ctypes.byref(ok), ctypes.byref(info),
            self.verbose, self.tol, self.max_iter, self.dav_iter, self.shift,
            self.memory, self.__memory_unit()
        )
        self.__check(info)
        return bool(ok.value)

    def dgl_smogd_driver(self, eig, evec, apbmul, ambmul, spdmul, smdmul, lrprec):
        """SMO-GD for linear response problems of size 2*n.
        eig: (n_max,), evec: (2*n, n_max). Returns True if converged."""
        p_eig = self.__array(eig, (self.n_max,), "eig")
        p_evec = self.__array(evec, (2 * self.n, self.n_max), "evec")
        c_apbmul = self.__matvec(apbmul)
        c_ambmul = self.__matvec(ambmul)
        c_spdmul = self.__matvec(spdmul)
        c_smdmul = self.__matvec(smdmul)
        c_lrprec = self.__lrprec(lrprec)
        ok = ctypes.c_bool(False)
        info = self.__cInt(0)

        self.lib.dgl_smogd_driver(
            2 * self.n, self.n_targ, self.n_max,
            c_apbmul, c_ambmul, c_spdmul, c_smdmul, c_lrprec,
            p_eig, p_evec,
            ctypes.byref(ok), ctypes.byref(info),
            self.verbose, self.tol, self.max_iter, self.dav_iter,
            self.memory, self.__memory_unit()
        )
        self.__check(info)
        return bool(ok.value)
