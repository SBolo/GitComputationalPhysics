#include <stdio.h>
#include <stdlib.h>

#include "a_lib.c"


double f(double x, void *params) {
    
    double *p = (double *)params;
    double T = p[0];
    double Tc = p[1];
    double h = p[2];
    
    return ( x - tanh( Tc * x / T + h/T ) );
}

int main() {
    
    double p[2];
    double t, s, dt = 1e-2, x; //indice della temperatura e media dello spin
    
    gsl_function F;
    
    F.function = f;
    F.params = p;
    
    p[1] = 10.; //temperatura critica
    p[2] = 0.;
    
    for(t = dt; t < 10; t += dt) {
        
        p[0] = t;
        s = zero(&F, 0.0001, 1.5);
        x = s;
        
        
        printf("%g \t %g \t %g\n", t, x, -p[2]/(p[1]*t));
    }

    
    return 0;
}





