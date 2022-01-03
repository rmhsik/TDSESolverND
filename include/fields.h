#ifndef FIELDS_H
#define FIELDS_H

#define SIN2 0
#define TRAP 1
    class Field{
        private:
            int _nt;
            bool _flag = false;
            double _dt, _tmax, _w;
            double *_t;
            double *_field;
            double env_sin2(double ti);
            double env_trap(double ti);
            
        public:
            Field();
            Field(double amp, double w, double phi, int env, double tmax, double *t, int nt);
            void calc_pot();
            double get(int i);
            double* get();
            double operator[](int i);
    };

#endif
