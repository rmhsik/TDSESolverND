#include <cmath>
#include "debug.h"
#include "diagnostics.h"
#include "utils.h"

Diagnostics::Diagnostics(){

}

Diagnostics::Diagnostics(const int n_probes, std::string defs){
    _n_probes = n_probes;
    //_def = def;

    char delim= '\n';
    int counter = 0;
    std::stringstream ss(defs);
    std::string item;
    while(getline(ss,item,delim)){
        _def_vec.push_back(item);
        counter++;
    }
    if(counter != n_probes){debug0("[Diagnostics] n_probes != number of definitions");}
}

void Diagnostics::_set_parameters(Parameters *param){
    _param = param;
    
}

void Diagnostics::_set_geometry(double *i, double *k, double *t, const double di, const double dk){
    _i = i; _k = k; _di = di; _dk = dk, _t = t;
   _ni = _param->ni;
   _nk = _param->nk;
}

void Diagnostics::_set_wf(WF *wf){
    _wf = wf;
}

void Diagnostics::_set_ham(Hamiltonian *ham){
    _ham = ham;
}

void Diagnostics::_set_tempmask(){
    _tempmask = new cdouble[_param->nt];
    double tmax_mask; 
    if (_param->env==SIN2){
        tmax_mask = _param->tmax_ev;
    }
    else{
        tmax_mask = _param->tmax_ev + _param->period*2.0;
    }

    for(int i=0; i<_param->nt;i++){
        if(_t[i]<tmax_mask){
            _tempmask[i] = 1.0;
        }
        else{
            _tempmask[i] = exp(-pow(_t[i]-tmax_mask,2)/100.0);
        }
    }
    //path = "results/accmask.dat";
    //write_array(_accmask,_param->nt,path);


}

int Diagnostics::_get_type(std::string _def){
    char delim = ',';
    std::stringstream ss(_def);
    std::string item;
    std::vector<std::string> parsed;
    int type;
    while(getline(ss,item,delim)){
        parsed.push_back(item);
    }
    if(parsed[0] == "acc_i") 
        type = ACC_I;
    else if(parsed[0] == "acc_k")
        type = ACC_K;
    else if(parsed[0] == "dip_i")
        type = DIP_I;
    else if(parsed[0] == "dip_k")
        type = DIP_K;
    else if(parsed[0] == "pop")
        type = POP;
    return type;
}

void Diagnostics::_create_probes(){
    ProbeFactory factory;
    for(int i=0; i<_n_probes;i++){
        _probe_vec.push_back(factory.create(_param->geometry,_def_vec[i]));    
        _probe_vec[i]->set_wf(_wf);
        _probe_vec[i]->set_ham(_ham);
        _probe_vec[i]->set_param(_param);
        _probe_vec[i]->set_geometry(_i,_k,_di,_dk);
        _probe_vec[i]->set_tempmask(_tempmask);
    } 
}

void Diagnostics::run_diagnostic(const int idx){
    for(int i=0; i<_n_probes;i++){
        _probe_vec[i]->calc(idx);
    } 
}


