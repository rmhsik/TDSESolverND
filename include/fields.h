#ifndef FIELDS_H
#define FIELDS_H
#include "parameters.h"

#define SIN2 0
#define TRAP 1
class Field{
private:
    Parameters *_param;
    int _nt;
    bool _flag = false;
    double _dt, _tmax, _w, _kw, _amp, _phi;
    double *_t;
    double **_field;
    double **_vecpot;
    double *_i,*_k,_di,_dk;
    double env_sin2(double ti);
    double env_trap(double ti);
    double plane_wave(double i, double k, double t);
    
public:
    Field();
    Field(double amp, double w, double kw, double phi, int env, double tmax, double *t, int nt);
    void set_params(Parameters *param);
    void set_geometry(double *i, double *k, const double di, const double dk);
    void calc_field(double t);
    void calc_pot();
    double get(int i, int k);
    double** get();
    double operator()(int i, int k);
    ~Field();
};

#endif
