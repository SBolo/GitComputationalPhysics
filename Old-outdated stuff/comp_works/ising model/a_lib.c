#include <math.h>

#define R 0.61803399
#define C (1.0-R) 
#define SHFT2(a,b,c) (a)=(b);(b)=(c); 
#define SHFT3(a,b,c,d) (a)=(b);(b)=(c);(c)=(d); 
#define GSL(F,x) (*((F)->function))(x,(F)->params)


typedef struct {
	double (*function)(double x, void *params);
	void *params;
} gsl_function;

double hyp_tan(double x) {
    
    double cos = exp(x) + exp(-x);
    double sin = exp(x) - exp(-x);
    
    return sin/cos;
}

double simpson_integrate(gsl_function *fn, double a, double b) {
//a,b estremi di integrazione
	double h = 0, somma = 0, k = 0, res = 0, risultato;
	int i = 0, t = 0, N = 1e5; //N = numero di campionamenti
	double f[N]; 
	
	h = (b-a)/(N); //passo di integrazione
	
				for(i=0; i<N; i++) {
					f[i] = GSL(fn,a+i*h); //campionamento della funzione
				}
				
			k = f[0] + 4 * f[1] + f[N-1] + 4 * f[N-2];

				for(i=2; i<N-2; i+=2) {
					res += 4 * f[i] + 2 * f[i+1];
				}	
	
			risultato = (k + res) * (h/3); //risultato dell'integrale
		
	return risultato;	
}

double golden_min(gsl_function *f, double ax, double bx, double cx, double tol)  { 
	double f1,f2,x0,x1,x2,x3;
	x0=ax; 
	x3=cx; 

	if (fabs(cx-bx) > fabs(bx-ax)) { 
		x1=bx; 
		x2=bx+C*(cx-bx); 
	} 
	else { 
		x2=bx; 
		x1=bx-C*(bx-ax); 
	}
	 
	f1=GSL(f,x1); 
	f2=GSL(f,x2); 
	
		while (fabs(x3-x0) > tol*(fabs(x1)+fabs(x2))) { 
			
			if (f2 < f1) { 
				SHFT3(x0,x1,x2,R*x1+C*x3)
				SHFT2(f1,f2,GSL(f,x2))
			} 
			else {  
				SHFT3(x3,x2,x1,R*x2+C*x0) 
				SHFT2(f2,f1,GSL(f,x1))  
			} 
		}
		  
	if (f1 < f2) { 
		return x1; 
	} 
	else { 
		return x2; 
	} 
} 

double golden_max(gsl_function *f, double ax, double bx, double cx, double tol)  { 
	double f1,f2,x0,x1,x2,x3;
	x0=ax; 
	x3=cx; 

	if (fabs(cx-bx) > fabs(bx-ax)) { 
		x1=bx; 
		x2=bx+C*(cx-bx); 
	} 
	else { 
		x2=bx; 
		x1=bx-C*(bx-ax); 
	}
	 
	f1=-GSL(f,x1); 
	f2=-GSL(f,x2); 
	
		while (fabs(x3-x0) > tol*(fabs(x1)+fabs(x2))) { 
			
			if (f2 < f1) { 
				SHFT3(x0,x1,x2,R*x1+C*x3)
				SHFT2(f1,f2,-GSL(f,x2))
			} 
			else {  
				SHFT3(x3,x2,x1,R*x2+C*x0) 
				SHFT2(f2,f1,-GSL(f,x1))  
			} 
		}
		  
	if (f1 < f2) { 
		return x1; 
	} 
	else { 
		return x2; 
	} 
} 

double zero(gsl_function *fn, double min, double max) {
// metodo della bisezione
	double m = 0, segno;
	
	do
		{

		segno = GSL(fn,min)*GSL(fn,max);
		
		if(segno > 0) {
			printf("#Non esiste uno zero in questo intervallo\n");
			break;
		}
		else {
			m = (min + max)/2;
			
			if(GSL(fn,m)*GSL(fn,min)<0) {
				max=m;
			}
			else {
				min=m;
			}
		}
		
	} while (max-min>0.00001);
		
		return m;	
}

void time(double t, double t_max) {
	printf("\t Progress: %g %% \r",100*t/t_max);
}
