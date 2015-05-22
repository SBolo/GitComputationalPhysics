#include <iostream>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

using namespace std;

#define SIMPSON_FUNCTION (0)
#define TRAPEZOID_FUNCTION (1)

class Integrate {
    public:
       Integrate();
       ~Integrate();
       double operator() (double (*)(double, double *), double *, int);
       double domain_start;
       double domain_end;

    private:
       double simpson(double (*)(double, double *), double *);
       double trapezoid(double (*)(double, double *), double *);
       double s_step;
       double t_step;

};

Integrate::Integrate() {
    s_step = 1.e-5;
    t_step = 1.e-4;
}
//------------------------------------------------------------------------------

Integrate::~Integrate() {
}
//------------------------------------------------------------------------------

double Integrate::operator() (double (*f)(double, double *), double *params, int method) {

    double r;

    if(method == 0) {
        r = simpson(f, params);
    } else if(method == 1) {
        r = trapezoid(f, params);
    }

    return r;
}
//------------------------------------------------------------------------------

double Integrate::simpson(double (*f)(double, double *), double *params) {
    double k, res = 0, domain_dimension = (domain_end - domain_start);
    int dim = (int)(domain_dimension / s_step);

    k = f(domain_start,params) + 4. * f(domain_start + s_step,params) + f(domain_start + s_step * (dim - 1), params) +
                    4. * f(domain_start + s_step * (dim - 2), params);

    for(int i = 2; i < dim-2; i+=2) {
        res += 4. * f(i * s_step + domain_start, params) + 2. * f((i+1) * s_step + domain_start, params);
    }

    return (k + res) * (s_step/3.);

}
//------------------------------------------------------------------------------

double Integrate::trapezoid(double (*f)(double, double *), double *params) {
    double k, res = 0, domain_dimension = (domain_end - domain_start);
    int dim = (int)(domain_dimension / t_step);

    for( int i = 0; i < dim; i++) {
        res += ( f(i * t_step + domain_start, params) + f( (i+1) * t_step + domain_start, params) ) * t_step/2.;
    }

    return res;
}
//------------------------------------------------------------------------------
