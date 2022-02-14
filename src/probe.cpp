#include <cmath>
#include "debug.h"
#include "diagnostics.h"
#include "utils.h"
#include "probe.h"

Probe::Probe(){

}

Probe::Probe(std::string def){
    _def = def;
    _parse_def();
    switch(_type){
        case ACC_I:
            _calc = &Probe::_acc_i;
            break;
        case ACC_K:
            _calc = &Probe::_acc_k;
            break;
        case DIP_I:
            _calc = &Probe::_dip_i;
            break;
        case DIP_K:
            _calc = &Probe::_dip_k;
            break;
        case POP:
            _calc = &Probe::_pop;
            break;
        case DENS:
            _calc = &Probe::_dens;
            break;
        case WFSNAP:
            _calc = &Probe::_wf_snap;
    } 
}

void Probe::set_wf(WF *wf){
    _wf = wf;
}

void Probe::set_ham(Hamiltonian *ham){
    _ham = ham;
}

void Probe::set_param(Parameters *param){
    _param = param;
}

void Probe::set_geometry(double* i, double *j, double *k, const double di, const double dj, const double dk){
    _i = i; _k = k; _ni = _param->ni; _nk =_param-> nk; 
    _di = di; _dk = dk; 
    _nt_diag = _param->nt_diag; 
    _nt = _param->nt;
    if (_type != DENS && _type != WFSNAP){
        std::cout<<"here \n";
        _data = new cdouble[_nt];
    }
    else{ 
        std::cout<<"here 2\n";
        _data = new cdouble[_nt/_nt_diag*_ni*_nk];
    }
}

void Probe::set_tempmask(cdouble* tempmask){
    _tempmask = tempmask;
}

void Probe::write_probe(){
    if(_type == DENS || _type == WFSNAP)
        write_array(_data,_nt/_nt_diag*_ni*_nk, _data_path);
    else
        write_array(_data,_nt, _data_path); 
}

void Probe::_parse_def(){
    char delim = ',';
    std::stringstream ss(_def);
    std::string item;
    std::vector<std::string> parsed;
    while(getline(ss,item,delim)){
        parsed.push_back(item);
    }
    if(parsed[0] == "acc_i") 
        _type = ACC_I;
    else if(parsed[0] == "acc_k")
        _type = ACC_K;
    else if(parsed[0] == "dip_i")
        _type = DIP_I;
    else if(parsed[0] == "dip_k")
        _type = DIP_K;
    else if(parsed[0] == "pop")
        _type = POP;
    else if(parsed[0] == "dens")
        _type = DENS;
    else if(parsed[0] == "wf")
        _type = WFSNAP;

    if (_type == ACC_I || _type == ACC_K || _type == DENS || _type == WFSNAP){
        _data_path = parsed[1];
        _int_imin = 0.0;
        _int_imax = 0.0;
        _int_kmin = 0.0;
        _int_kmax = 0.0;
    }
    else{
        _int_imin = std::stod(parsed[1]);
        _int_imax = std::stod(parsed[2]);
        _int_kmin = std::stod(parsed[3]);
        _int_kmax = std::stod(parsed[4]);
        _data_path = parsed[5];
    }

    std::cout<<"type: "<<_type<<std::endl;   
    std::cout<<"int_imin: "<<_int_imin<<std::endl;   
    std::cout<<"int_imax: "<<_int_imax<<std::endl;   
    std::cout<<"int_kmin: "<<_int_kmin<<std::endl;   
    std::cout<<"int_kmax: "<<_int_kmax<<std::endl;   
    std::cout<<"data_path: "<<_data_path<<std::endl;   
}

void Probe::calc(const int idx){
    (this->*(this->_calc))(idx);
}

Probe::~Probe(){
    delete[] _data;
}

