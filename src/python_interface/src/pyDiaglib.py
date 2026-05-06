import sys
import ctypes

class diaglib:

    def __init__(self,libPath,kind,
                n, n_targ, n_max,
                verbose = False,
                max_iter = 100,
                dav_iter = 20,
                tol = 1e-8,
                shift = 0.0,
                memory = 80,
                memory_unit = "MB"):
        # Load the library
        try:
            self.lib = ctypes.CDLL(libPath)
        except:
            print(f"\nCould not find diaglib at {libPath}")
            sys.exit(1)

        # Set correct kind for integers
        if kind == 4:
            cInt = ctypes.c_int
        elif kind == 8:
            cInt = ctypes.c_long
        cDouble = ctypes.c_double
        cBool = ctypes.c_bool
        cCharP = ctypes.c_char_p

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
            ctypes.POINTER(cDouble), ctypes.POINTER(cDouble)
        )

        PRECND = ctypes.CFUNCTYPE(None,
            ctypes.POINTER(cInt), ctypes.POINTER(cInt),
            ctypes.POINTER(cDouble),
            ctypes.POINTER(cDouble), ctypes.POINTER(cDouble)
        )

        METVEC = ctypes.CFUNCTYPE(None,
            ctypes.POINTER(cInt), ctypes.POINTER(cInt),
            ctypes.POINTER(cDouble), ctypes.POINTER(cDouble)
        )

        LRPREC = ctypes.CFUNCTYPE(None,
            ctypes.POINTER(cInt), ctypes.POINTER(cInt), ctypes.POINTER(cDouble),
            ctypes.POINTER(cDouble), ctypes.POINTER(cDouble),
            ctypes.POINTER(cDouble), ctypes.POINTER(cDouble)
        )

        NONE = ctypes.CFUNCTYPE(None,
            ctypes.POINTER(cInt), ctypes.POINTER(cInt),
            ctypes.POINTER(cDouble), ctypes.POINTER(cDouble)
        )

        self.__MATVEC = MATVEC
        self.__METVEC = METVEC
        self.__PRECND = PRECND
        self.__LRPREC = LRPREC
        self.__NONE = ctypes.cast(None,NONE)

        # Definisci la funzione Fortran
        self.lib.dgl_davidson_driver.argtypes = [
            cInt, cInt, cInt, MATVEC, PRECND, METVEC,
            ctypes.POINTER(cDouble), ctypes.POINTER(cDouble), ctypes.POINTER(cBool),
            cBool, cDouble, cInt, cInt, cDouble, cInt, cCharP
        ]

        self.lib.dgl_lobpcg_driver.argtypes = [
            cInt, cInt, cInt, MATVEC, PRECND, METVEC,
            ctypes.POINTER(cDouble), ctypes.POINTER(cDouble), ctypes.POINTER(cBool),
            cBool, cDouble, cInt, cDouble, cInt, cCharP
        ]

        self.lib.dgl_davidson_nosym_driver.argtypes = [
            cInt, cInt, cInt, MATVEC, MATVEC, PRECND, cCharP,
            ctypes.POINTER(cDouble), ctypes.POINTER(cDouble), ctypes.POINTER(cDouble), ctypes.POINTER(cBool),
            cBool, cDouble, cInt, cInt, cDouble, cInt, cCharP
        ]

        self.lib.dgl_smogd_driver.argtypes = [
            cInt, cInt, cInt, MATVEC, MATVEC, MATVEC, MATVEC, LRPREC,
            ctypes.POINTER(cDouble), ctypes.POINTER(cDouble), ctypes.POINTER(cBool),
            cBool, cDouble, cInt, cInt, cInt, cCharP
        ]

        self.lib.dgl_ortho_cd.argtypes = [
            cInt, cInt, ctypes.POINTER(cDouble), ctypes.POINTER(cDouble), ctypes.POINTER(cBool)
        ]

        self.lib.dgl_ortho_vs_x.argtypes = [
            cInt, cInt, cInt, ctypes.POINTER(cDouble), ctypes.POINTER(cDouble)
        ]

        self.lib.dgl_b_ortho_vs_x.argtypes = [
            cInt, cInt, cInt, ctypes.POINTER(cDouble), ctypes.POINTER(cDouble), ctypes.POINTER(cDouble)
        ]

    def dgl_davidson_driver(self, eig, evec, matvec, precnd, metvec=None):
        matvec = self.__MATVEC(matvec)
        precnd = self.__PRECND(precnd)
        if metvec is None:
            metvec = self.__NONE
        else:
            metvec = self.__METVEC(metvec)
        ok = ctypes.c_bool(False)

        self.lib.dgl_davidson_driver(
            self.n, self.n_targ, self.n_max,
            matvec, precnd, metvec,
            eig.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            evec.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            ctypes.byref(ok),
            self.verbose,
            self.tol,
            self.max_iter,
            self.dav_iter,
            self.shift,
            self.memory,
            ctypes.c_char_p(bytes(self.memory_unit,"utf-8"))
        )
        return ok


    def dgl_lobpcg_driver(self, eig, evec, matvec, precnd, metvec=None):
        matvec = self.__MATVEC(matvec)
        precnd = self.__PRECND(precnd)
        if metvec is None:
            metvec = self.__NONE
        else:
            metvec = self.__METVEC(metvec)
        ok = ctypes.c_bool(False)

        self.lib.dgl_lobpcg_driver(
            self.n, self.n_targ, self.n_max,
            matvec, precnd, metvec,
            eig.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            evec.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            ctypes.byref(ok),
            self.verbose,
            self.tol,
            self.max_iter,
            self.shift,
            self.memory,
            ctypes.c_char_p(bytes(self.memory_unit,"utf-8"))
        )
        return ok
    

    def dgl_davidson_nosym_driver(self, eig, evec1, evec2, side, matvec_r, matvec_l, precnd):
        matvec_r = self.__MATVEC(matvec_r)
        matvec_l = self.__MATVEC(matvec_l)
        precnd = self.__PRECND(precnd)
        ok = ctypes.c_bool(False)

        self.lib.dgl_davidson_nosym_driver(
            self.n, self.n_targ, self.n_max,
            matvec_r, matvec_l, precnd,
            ctypes.c_char_p(bytes(side,"utf-8")),
            eig.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            evec1.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            evec2.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            ctypes.byref(ok),
            self.verbose,
            self.tol,
            self.max_iter,
            self.dav_iter,
            self.shift,
            self.memory,
            ctypes.c_char_p(bytes(self.memory_unit,"utf-8"))
        )
        return ok

    def dgl_smogd_driver(self,eig,evec,apbmul,ambmul,spdmul,smdmul,lrprec):
        apbmul = self.__MATVEC(apbmul)
        ambmul = self.__MATVEC(ambmul)
        spdmul = self.__MATVEC(spdmul)
        smdmul = self.__MATVEC(smdmul)
        lrprec = self.__LRPREC(lrprec)
        ok = ctypes.c_bool(False)

        self.lib.dgl_smogd_driver(
            2*self.n, self.n_targ, self.n_max,
            apbmul, ambmul, spdmul, smdmul, lrprec,
            eig.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            evec.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            ctypes.byref(ok),
            self.verbose,
            self.tol,
            self.max_iter,
            self.dav_iter,
            self.memory,
            ctypes.c_char_p(bytes(self.memory_unit,"utf-8"))
        )
        return ok

    def dgl_ortho_cd(self,n,m,u):

        ok = ctypes.c_bool(False)
        growth = ctypes.c_double(0.)
        self.lib.dgl_ortho_cd(n, m, 
                              u.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                              growth, ok)

        return growth, ok
    
    def dgl_ortho_vs_x(self,n, m, k, x, u):

        self.lib.dgl_ortho_vs_x(n, m, k,
                                x.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                                u.ctypes.data_as(ctypes.POINTER(ctypes.c_double)))

    def dgl_b_ortho_vs_x(self,n, m, k, x, bx, u):

        self.lib.dgl_b_ortho_vs_x(n, m, k,
                                  x.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                                  bx.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
                                  u.ctypes.data_as(ctypes.POINTER(ctypes.c_double)))