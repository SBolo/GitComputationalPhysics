#include <stdio.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_eigen.h>
#include <math.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_integration.h>
#include <armadillo>

using namespace arma;

//_____________________________________________________________GSL_MATRIX_MACROS________________________________________________

#define GSL(F,x) (*((F)->function))(x,(F)->params)
#define STEP (9.e-3) //step definitivo che posso accettare
#define RMAX (20.) //dimensione mesh

//________________________________________________________________END_MACROS____________________________________________________

typedef struct {
    vec rho;
    gsl_integration_glfixed_table * w;
} schrodinger;

double integrate(vec g, double j, double h) { //array, primo estremo (primo indice), secondo estremo (secondo indice)

    int i;
    double risultato = 0.;

    for( i = j; i < h-1; i++) {
        risultato += (g(i) + g(i+1)) * STEP/2.;
    }

    return risultato;
}

double simpson(vec f, int n, int m) { //normalizzazione
    double h, k, res = 0;
    int i = 0;

    k = f(n) + 4 * f(n+1) + f(m-1) + 4 * f(m-2);

        for(i = n+2; i < m-2; i+=2) {
            res += 4 * f(i) + 2 * f(i+1);
        }

    return (k + res) * (STEP/3);

}



//OK
double hartree(double x, void *params) {

    schrodinger *D = (schrodinger *)params;
    int i;
    double ext, pot1, pot2;

    pot2 = simpson(D->rho, 0, (int)(x/STEP));

    return pot2;
}


int main() {

    vec rho(RMAX/STEP+1);
    schrodinger D;
    gsl_function F;

    F.function = hartree;
    F.params = &D;

    gsl_integration_glfixed_table * w = gsl_integration_glfixed_table_alloc(500);

    D.rho = rho;
    D.w = w;

   for(double x = 0; x < RMAX; x+=STEP) {
        D.rho((int)(x/STEP)) = x*exp(-x*x);
    }

   double result_int = integrate(D.rho, 0, (int)(5./STEP));
   //double result_gsl = gsl_integration_glfixed (&F, 0., 5, w);
   double result_simpson = simpson(D.rho, 0, (int)(5./STEP));

   cout << "Mia routine = " << result_int << endl;
   cout << "Hartree = " << hartree(5,&D) << endl;

    gsl_integration_glfixed_table_free(D.w);

    return 0;
}
