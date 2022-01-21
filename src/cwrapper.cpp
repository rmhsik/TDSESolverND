#include "parameters.h"

extern "C"{
    Parameters *Parameters_new(){ return new Parameters;}
    void Parameters_print(Parameters *p){p->print(); }
    void Parameters_n_threads(Parameters *p, int val){p->n_threads = val;}
    void Parameters_init_wf(Parameters *p, int val){p->init_wf = val;}
    void Parameters_geometry(Parameters *p, int val){p->geometry = val;}
    void Parameters_ni(Parameters *p, int val){p->ni = val;}
    void Parameters_imin(Parameters *p, double val){p->imin = val;}
    void Parameters_imax(Parameters *p, double val){p->imax = val;}
    void Parameters_nk(Parameters *p, int val){p->nk = val;}
    void Parameters_kmin(Parameters *p, double val){p->kmin = val;}
    void Parameters_kmax(Parameters *p, double val){p->kmax = val;}
    void Parameters_w0(Parameters *p, double val){p->w0 = val;}
    void Parameters_period(Parameters *p, double val){p->period = val;}
    void Parameters_tmax_ev(Parameters *p, double val){p->tmax_ev = val;}
    void Parameters_tmax_sim(Parameters *p, double val){p->tmax_sim = val;}
    void Parameters_dt(Parameters *p, double val){p->dt = val;}
    void Parameters_dt_ITP(Parameters *p, double val){p->dt_ITP = val;}
    void Parameters_nt(Parameters *p, int val){p->nt = val;}
    void Parameters_nt_ITP(Parameters *p, int val){p->nt_ITP = val;}
    void Parameters_nt_diag(Parameters *p, int val){p->nt_diag = val;}
    void Parameters_env(Parameters *p, int val){p->env = val;}
    void Parameters_w0Ei(Parameters *p, double val){p->w0Ei = val;}
    void Parameters_w0Ek(Parameters *p, double val){p->w0Ek = val;}
    void Parameters_w0Bi(Parameters *p, double val){p->w0Bi = val;}
    void Parameters_w0Bk(Parameters *p, double val){p->w0Bk = val;}
    void Parameters_E0i(Parameters *p, double val){p->E0i = val;}
    void Parameters_E0k(Parameters *p, double val){p->E0k = val;}
    void Parameters_B0i(Parameters *p, double val){p->B0i = val;}
    void Parameters_B0k(Parameters *p, double val){p->B0k = val;}
    void Parameters_phiEi(Parameters *p, double val){p->phiEi = val;}
    void Parameters_phiEk(Parameters *p, double val){p->phiEk = val;}
    void Parameters_phiBi(Parameters *p, double val){p->phiBi = val;}
    void Parameters_phiBk(Parameters *p, double val){p->phiEk = val;}
    void Parameters_acc_path(Parameters *p, char *val){p->set_acc_path(val);}
    void Parameters_dip_path(Parameters *p, char *val){p->set_dip_path(val);}
}
