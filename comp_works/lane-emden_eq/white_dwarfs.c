#include <stdio.h>
#include <stdlib.h>
#include <math.h>

//#define N 3./2. //non relativistic white dwarf
#define N 3. //relativistic white dwarf
#define STEP 1e-6


double f(double *fnz, double x) {
	return fnz[1] / (x * x);
}

double g(double *fnz, double x) {
	return - pow(fnz[0], N) * x * x;
}

void refresh(double *fnz, double *x_temp, double l, double k, double d) {
	x_temp[0] = fnz[0] + STEP * k/d;
	x_temp[1] = fnz[1] + STEP * l/d;
} 

void II_ORD_RG4_STEP(double *fnz, double t) {
	double k1, k2, k3, k4, l1, l2, l3, l4;
	double x_temp[2];
	
	x_temp[0] = fnz[0];
    x_temp[1] = fnz[1];
	
	k1 = f(x_temp, t);
	l1 = g(x_temp, t);

        refresh(fnz, x_temp, l1, k1, 2.);
	
	k2 = f(x_temp, t);
	l2 = g(x_temp, t);

        refresh(fnz, x_temp, l2, k2, 2.);
	
	k3 = f(x_temp, t);
	l3 = g(x_temp, t);

        refresh(fnz, x_temp, l3, k3, 1.);
	
	k4 = f(x_temp, t);
	l4 = g(x_temp, t);
	
	fnz[0] += STEP * (1./6.) * (k1 + 2.*k2 + 2.*k3 + k4);
	fnz[1] += STEP * (1./6.) * (l1 + 2.*l2 + 2.*l3 + l4);
}


int main() {
	
	double fnz[2], t;
	
	FILE *pt;
    
    //initial conditions
	fnz[0] = 1;
	fnz[1] = 0;
    
    t = STEP; //nella funzione divido per t^2, quindi parto da dt

	//pt = fopen("white_dwarf","w");
        do {
        
                II_ORD_RG4_STEP(fnz,t);
                fprintf(pt,"%g \t %g \t %g\n",t, fnz[0], fabs(fnz[1]));
                t += STEP;
        
        } while (fnz[0] > 1e-7); //dt è  la precisione con cui mi aspetto di stimare lo 0
    
        printf("Valori interessanti: \n \t THETA(XI_1) = %g, DERIVATA*(XI_1)^2 = %g\n", t, fabs(fnz[1]));
			
	fclose(pt);
    
    //system("gnuplot wd-gp.gp");
	
	
    return 0;	
}
	 