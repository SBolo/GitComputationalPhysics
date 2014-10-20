#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/*
 PARAMETERS
 */
#define STEP 1.e-3

/*
 * dv/dt = -g/v
 * dx/dt = v
 *
 * d( x| v) = (v | -g/v)
 */


//________________________________________________EQUATIONS__________________________________________________
double f(double *fnz, double t, double *params) { //EQUATION FOR x (x sarà fnz[0])
    double g = params[0];
    double mu = params[1];
    double k = params[2];

    return -g/fnz[1];// - mu * fnz[1] - k * fnz[0];
}

double g(double *fnz, double t, double *params) { //EQUATION FOR MASS (la massa sarà fnz[1])
    return fnz[1];
}
//___________________________________________________END_____________________________________________________


//_________________________________________________RG4_ALGORITHM_________________________________________________
void refresh(double *fnz, double *x_temp, double l, double k, double d) { //refresh of the parameters in Runge-Kutta
	x_temp[0] = fnz[0] + STEP * k/d;
	x_temp[1] = fnz[1] + STEP * l/d;
} 

void II_ORD_RG4_STEP(double *fnz, double t, double *params) {
	
    double k1, k2, k3, k4, l1, l2, l3, l4;
	double x_temp[2];
	
	x_temp[0] = fnz[0];
    x_temp[1] = fnz[1];
	
    k1 = f(x_temp, t, params);
    l1 = g(x_temp, t, params);

        refresh(fnz, x_temp, l1, k1, 2.);
	
    k2 = f(x_temp, t, params);
    l2 = g(x_temp, t, params);

        refresh(fnz, x_temp, l2, k2, 2.);
	
    k3 = f(x_temp, t, params);
    l3 = g(x_temp, t, params);

        refresh(fnz, x_temp, l3, k3, 1.);
	
    k4 = f(x_temp, t, params);
    l4 = g(x_temp, t, params);
	
	fnz[0] += STEP * (1./6.) * (k1 + 2.*k2 + 2.*k3 + k4);
	fnz[1] += STEP * (1./6.) * (l1 + 2.*l2 + 2.*l3 + l4);
}
//_____________________________________________________END_RG4__________________________________________________


int main() {
	
    double fnz[2], params[3], t; //fnz[0] is the density, fnz[1] is the mass & t is the radius (generic parameter t in RG4)
	FILE *pt;

    params[0] = 2.;
    params[1] = 1.;
    params[2] = 3.;
    


        //initial conditions
        t = 0.; //initial time
        fnz[0] = 1.; //initial density - converted in x
        fnz[1] = 3.; //initial mass - no divergence for m(r = 0) = 0

          do {
                    t += STEP; //radius increment
                    II_ORD_RG4_STEP(fnz,t,params); //step of RG4 at radius t
                    printf("%g \t %g\n",t, fnz[1]);

                    if(fnz[1] < 1e-3) {
                        exit(0);
                    }
        
            } while ( t < 10 ); //stop integrating when you reach the 0 value of x (which is the minimum for the density)

   
	
	
    return 0;	
}
	 
