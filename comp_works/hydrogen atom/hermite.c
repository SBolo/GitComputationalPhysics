#include <stdio.h>
#include <math.h>

#define STEP 0.001
#define GSL(F,x) (*((F)->function))(x,(F)->params)


typedef struct {
	double (*function)(double x, void *params);
	void *params;
} gsl_function;


double d2(gsl_function *f, double x, int n) {
    return ( GSL(f,x - STEP) + GSL(f,x + STEP) - 2. * GSL(f,x) ) / STEP*STEP ;
}

double simpson_integrate(double *f, int N) {
	double h, k, res = 0;
	int i = 0;
	
    k = f[0] + 4 * f[1] + f[N-1] + 4 * f[N-2];
    
	for(i = 2; i < N-2; i+=2) {
		res += 4 * f[i] + 2 * f[i+1];
	}
	
	return (k + res) * (STEP/3);
	
}

double polynomial(double x, int n) {
    
    if(n == 0) return 1;
        else if(n == 1) return 2.*x;
            else return ( 2.* x * polynomial(x,n-1) - 2.* (n-1) * polynomial(x,n-2) ); //Polinomi di Hermite

}


double basis(double x, void *params) {
    
    int *p = (int *)params;
    int n = p[0];
    //int dim = (int)(RMAX/STEP);
    //double tmp[dim];
    //double integ;
    
    /*
    int i = x/STEP, k = 0;
    
    for(k = 0; k < RMAX; k ++) {
        tmp[k] = polynomial(x,n) * exp(-x*x/2.);
    }
    
    integ = simpson_integrate(tmp,dim);
    
    for(k = 0; k < RMAX; k ++) {
        tmp[k] /= integ;
    }
    
    return tmp[i];*/
    return polynomial(x,n) * exp(-x*x/2.);
    
}

double hamiltonian(gsl_function *H, double x, int n) {
    
    
}

int main() {
    
    double x;
    int n[1];
    
    n[0] = 1;
    
    gsl_function B, V; //gsl_function per la base normalizzata e per il potenziale
    B.function = basis;
    B.params = n;
    
    //for(n = 0; n < 1000; n++) {
        for( x = 0; x < 10; x += STEP ) {
            printf("%g \t %g\n", x, GSL(&B,x));
        }
    //}
    
    return 0;
}
