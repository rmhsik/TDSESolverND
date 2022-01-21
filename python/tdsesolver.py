import ctypes as ct
import sys

path = '../lib/libtdsesolver.so'

lib = ct.cdll.LoadLibrary(path)


class TDSESolver:
    def __init__(self,param_class):
        lib.TDSESolver_new.argtypes = [ct.c_void_p]
        lib.TDSESolver_new.restype = ct.c_void_p
        lib.TDSESolver_ipropagate.argtypes = [ct.c_void_p]
        lib.TDSESolver_ipropagate.restype = ct.c_void_p

        lib.TDSESolver_propagate.argtypes = [ct.c_void_p]
        lib.TDSESolver_propagate.restype = ct.c_void_p

        self._obj = lib.TDSESolver_new(param_class)

class Parameters:
    def __init__(self):
        lib.Parameters_new.argtypes = []
        lib.Parameters_new.restype  = ct.c_void_p
        lib.Parameters_print.argtypes = [ct.c_void_p]
        lib.Parameters_print.restype = ct.c_void_p
        lib.Parameters_n_threads.argtypes = [ct.c_void_p, ct.c_int]
        lib.Parameters_n_threads.restype = ct.c_void_p
        lib.Parameters_init_wf.argtypes = [ct.c_void_p, ct.c_int]
        lib.Parameters_init_wf.restype = ct.c_void_p
        lib.Parameters_geometry.argtypes = [ct.c_void_p, ct.c_int]
        lib.Parameters_geometry.restype = ct.c_void_p
        lib.Parameters_ni.argtypes = [ct.c_void_p, ct.c_int]
        lib.Parameters_ni.restype = ct.c_void_p
        lib.Parameters_imin.argtypes = [ct.c_void_p, ct.c_double]
        lib.Parameters_imin.restype = ct.c_void_p
        lib.Parameters_imax.argtypes = [ct.c_void_p, ct.c_double]
        lib.Parameters_imax.restype = ct.c_void_p
        lib.Parameters_nk.argtypes = [ct.c_void_p, ct.c_int]
        lib.Parameters_nk.restype = ct.c_void_p
        lib.Parameters_kmin.argtypes = [ct.c_void_p, ct.c_double]
        lib.Parameters_kmin.restype = ct.c_void_p
        lib.Parameters_kmax.argtypes = [ct.c_void_p, ct.c_double]
        lib.Parameters_kmax.restype = ct.c_void_p
        lib.Parameters_w0.argtypes = [ct.c_void_p, ct.c_double]
        lib.Parameters_w0.restype = ct.c_void_p
        lib.Parameters_period.argtypes = [ct.c_void_p, ct.c_double]
        lib.Parameters_period.restype = ct.c_void_p
        lib.Parameters_tmax_ev.argtypes = [ct.c_void_p, ct.c_double]
        lib.Parameters_tmax_ev.restype = ct.c_void_p
        lib.Parameters_tmax_sim.argtypes = [ct.c_void_p, ct.c_double]
        lib.Parameters_tmax_sim.restype = ct.c_void_p
        lib.Parameters_dt.argtypes = [ct.c_void_p, ct.c_double]
        lib.Parameters_dt.restype = ct.c_void_p
        lib.Parameters_dt_ITP.argtypes = [ct.c_void_p, ct.c_double]
        lib.Parameters_dt_ITP.restype = ct.c_void_p
        lib.Parameters_nt.argtypes = [ct.c_void_p, ct.c_int]
        lib.Parameters_nt.restype = ct.c_void_p
        lib.Parameters_nt_ITP.argtypes = [ct.c_void_p, ct.c_int]
        lib.Parameters_nt_ITP.restype = ct.c_void_p
        lib.Parameters_nt_diag.argtypes = [ct.c_void_p, ct.c_int]
        lib.Parameters_nt_diag.restype = ct.c_void_p
        lib.Parameters_env.argtypes = [ct.c_void_p, ct.c_int]
        lib.Parameters_env.restype = ct.c_void_p
        lib.Parameters_w0Ei.argtypes = [ct.c_void_p, ct.c_double]
        lib.Parameters_w0Ei.restype = ct.c_void_p
        lib.Parameters_w0Ek.argtypes = [ct.c_void_p, ct.c_double]
        lib.Parameters_w0Ek.restype = ct.c_void_p
        lib.Parameters_w0Bi.argtypes = [ct.c_void_p, ct.c_double]
        lib.Parameters_w0Bi.restype = ct.c_void_p
        lib.Parameters_w0Bk.argtypes = [ct.c_void_p, ct.c_double]
        lib.Parameters_w0Bk.restype = ct.c_void_p
        lib.Parameters_E0i.argtypes = [ct.c_void_p, ct.c_double]
        lib.Parameters_E0i.restype = ct.c_void_p
        lib.Parameters_E0k.argtypes = [ct.c_void_p, ct.c_double]
        lib.Parameters_E0k.restype = ct.c_void_p
        lib.Parameters_B0i.argtypes = [ct.c_void_p, ct.c_double]
        lib.Parameters_B0i.restype = ct.c_void_p
        lib.Parameters_B0k.argtypes = [ct.c_void_p, ct.c_double]
        lib.Parameters_B0k.restype = ct.c_void_p
        lib.Parameters_phiEi.argtypes = [ct.c_void_p, ct.c_double]
        lib.Parameters_phiEi.restype = ct.c_void_p
        lib.Parameters_phiEk.argtypes = [ct.c_void_p, ct.c_double]
        lib.Parameters_phiEk.restype = ct.c_void_p
        lib.Parameters_phiBi.argtypes = [ct.c_void_p, ct.c_double]
        lib.Parameters_phiBi.restype = ct.c_void_p
        lib.Parameters_phiBk.argtypes = [ct.c_void_p, ct.c_double]
        lib.Parameters_phiBk.restype = ct.c_void_p
        lib.Parameters_acc_path.argtypes = [ct.c_void_p, ct.c_char_p]
        lib.Parameters_acc_path.restype = ct.c_void_p
        lib.Parameters_dip_path.argtypes = [ct.c_void_p, ct.c_char_p]
        lib.Parameters_dip_path.restype = ct.c_void_p
        self._obj = lib.Parameters_new()

    def print(self):
        lib.Parameters_print(self._obj)
    def set(self,param):
        self.n_threads(param["n_threads"])
        self.init_wf(param["init_wf"])
        self.geometry(param["geometry"])
        self.ni(param["ni"])
        self.imin(param["imin"])
        self.imax(param["imax"])
        self.nk(param["nk"])
        self.kmin(param["kmin"])
        self.kmax(param["kmax"])
        self.w0(param["w0"])
        self.period(param["period"])
        self.tmax_ev(param["tmax_ev"])
        self.tmax_sim(param["tmax_sim"])
        self.dt(param["dt"])
        self.dt_ITP(param["dt_ITP"])
        self.nt(param["nt"])
        self.nt_ITP(param["nt_ITP"])
        self.nt_diag(param["nt_diag"])
        self.env(param["env"])
        self.w0Ei(param["w0Ei"])
        self.w0Ek(param["w0Ek"])
        self.w0Bi(param["w0Bi"])
        self.w0Bk(param["w0Bk"])
        self.E0i(param["E0i"])
        self.E0k(param["E0k"])
        self.B0i(param["B0i"])
        self.B0k(param["B0k"])
        self.phiEi(param["phiEi"])
        self.phiEk(param["phiEk"])
        self.phiBi(param["phiBi"])
        self.phiBk(param["phiBk"])
        self.acc_path(param["acc_path"])
        self.dip_path(param["dip_path"])
    
    def n_threads(self, val):
        lib.Parameters_n_threads(self._obj, val)
    def init_wf(self, val):
        lib.Parameters_init_wf(self._obj, val)
    def geometry(self, val):
        lib.Parameters_geometry(self._obj, val)
    def ni(self, val):
        lib.Parameters_ni(self._obj, val)
    def imin(self, val):
        lib.Parameters_imin(self._obj, val)
    def imax(self, val):
        lib.Parameters_imax(self._obj, val)
    def nk(self, val):
        lib.Parameters_nk(self._obj, val)
    def kmin(self, val):
        lib.Parameters_kmin(self._obj, val)
    def kmax(self, val):
        lib.Parameters_kmax(self._obj, val)
    def w0(self, val):
        lib.Parameters_w0(self._obj, val)
    def period(self, val):
        lib.Parameters_period(self._obj, val)
    def tmax_ev(self, val):
        lib.Parameters_tmax_ev(self._obj, val)
    def tmax_sim(self, val):
        lib.Parameters_tmax_sim(self._obj, val)
    def dt(self, val):
        lib.Parameters_dt(self._obj, val)
    def dt_ITP(self, val):
        lib.Parameters_dt_ITP(self._obj, val)
    def nt(self, val):
        lib.Parameters_nt(self._obj, val)
    def nt_ITP(self, val):
        lib.Parameters_nt_ITP(self._obj, val)
    def nt_diag(self, val):
        lib.Parameters_nt_diag(self._obj, val)
    def env(self, val):
        lib.Parameters_env(self._obj, val)
    def w0Ei(self, val):
        lib.Parameters_w0Ei(self._obj, val)
    def w0Ek(self, val):
        lib.Parameters_w0Ek(self._obj, val)
    def w0Bi(self, val):
        lib.Parameters_w0Bi(self._obj, val)
    def w0Bk(self, val):
        lib.Parameters_w0Bk(self._obj, val)
    def E0i(self, val):
        lib.Parameters_E0i(self._obj, val)
    def E0k(self, val):
        lib.Parameters_E0k(self._obj, val)
    def B0i(self, val):
        lib.Parameters_B0i(self._obj, val)
    def B0k(self, val):
        lib.Parameters_B0k(self._obj, val)
    def phiEi(self, val):
        lib.Parameters_phiEi(self._obj, val)
    def phiEk(self, val):
        lib.Parameters_phiEk(self._obj, val)
    def phiBi(self, val):
        lib.Parameters_phiBi(self._obj, val)
    def phiBk(self, val):
        lib.Parameters_phiBk(self._obj, val)
    def acc_path(self, val):
        lib.Parameters_acc_path(self._obj, val.encode("utf-8"))
    def dip_path(self, val):
        lib.Parameters_dip_path(self._obj, val.encode("utf-8"))

