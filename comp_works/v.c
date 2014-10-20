#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/*
 PARAMETERS
 */
#define STEP 1.e-3


//________________________________________________EQUATIONS__________________________________________________
double f(double *fnz, double g) { //EQUATION FOR x (x sarà fnz[0])
    return -g/fnz[0];
}
//___________________________________________________END_____________________________________________________


//_________________________________________________RG4_ALGORITHM_________________________________________________
void refresh(double *fnz, double *x_temp, double k, double d) { //refresh of the parameters in Runge-Kutta
	x_temp[0] = fnz[0] + STEP * k/d;
} 

void II_ORD_RG4_STEP(double *fnz, double t, double g) {
	
    double k1, k2, k3, k4, l1, l2, l3, l4;
	double x_temp[2];
	
	x_temp[0] = fnz[0];
	
    k1 = f(x_temp, g);

        refresh(fnz, x_temp, k1, 2.);
	
    k2 = f(x_temp, g);

        refresh(fnz, x_temp, k2, 2.);
	
    k3 = f(x_temp, g);

        refresh(fnz, x_temp, k3, 1.);
	
    k4 = f(x_temp, g);
	
	fnz[0] += STEP * (1./6.) * (k1 + 2.*k2 + 2.*k3 + k4);
}
//_____________________________________________________END_RG4__________________________________________________


int main() {
	
    double fnz[1], t, g; //fnz[0] is the density, fnz[1] is the mass & t is the radius (generic parameter t in RG4)
	FILE *pt;

    g = 5.;

        //initial conditions
        t = 0.; //initial time
        fnz[0] = 3; //initial velocity

          do {
                    t += STEP; //radius increment
                    II_ORD_RG4_STEP(fnz,t,g); //step of RG4 at radius t
                    printf("%g \t %g\n",t, fnz[0]);
                    if (fnz[0] < 1e-3) {
                        exit(0);
                    }
        
            } while ( t < 10 ); //stop integrating when you reach the 0 value of x (which is the minimum for the density)

    return 0;
}
	 
