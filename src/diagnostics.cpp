#include <cmath>
#include "diagnostics.h"
#include "utils.h"

Diagnostics::Diagnostics(){

}

Diagnostics::Diagnostics(std::string def){
    _def = def;
   }

void Diagnostics::set_parameters(Parameters *param){
    _param = param;
    _data = new cdouble[_param->nt];

}

void Diagnostics::set_geometry(double *i, double *k, const double di, const double dk){
    _i = i; _k = k; _di = di; _dk = dk;
   _ni = _param->ni;
   _nk = _param->nk;
}

cdouble* Diagnostics::get_data(){
    return _data;
}

void Diagnostics::write_diagnostics(){
    write_array(_data, _param->nt, _data_path);
}

void Diagnostics::_parse_def(){
    char delim = ",";
    std::stringstream ss(_def);
    std::string item;
    std::vector<std::string> parsed;
    while(getline(ss,item,delim)){
        parsed.push_back(item);
    }
    _
    
}
