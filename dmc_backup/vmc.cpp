#include <iostream>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <gsl/gsl_rng.h>
#include <time.h>
#include "integrator.cpp"

using namespace std;

//#define DEBUG

//------------------------------------------MACROS-----------------------------------------
#define MDIM (10) //mesh dimension
#define WALKERS (400) //number of walkers
#define DTAU (1.e-4) //imaginary time step
#define TIME (10.) //propagation time
#define TMIN (2.) //minimum time after the energy sampling begins
#define NPARAMS (1) //number of variational parameters

#define A_MIN (0.3) //valore minimo per il parametro variazionale
#define A_MAX (0.8) //valore massimo per il parametro variazionale
#define A_STEP (0.1) //step per il parametro variazionale
#define EDIM (int)(A_MAX/A_STEP - A_MIN/A_STEP) //dimensione per l'array che conterrà le varie energie
//----------------------------------------END-MACROS---------------------------------------

class RandomNumbers {
    public:
        RandomNumbers(); //constructor
        ~RandomNumbers(); //destructor
        double operator() (); //operator overloading

    private:
        const gsl_rng_type *T; //holds static information about each type of generator
        gsl_rng *r; //describes an instance of a generator created from a given gsl_rng_type
};


class Vmc {
    public:

        Vmc(double (*)(double, double *), double (*)(double, double *), int, int, double); //constructor
        ~Vmc(); //destructor

        double Et; //trial energy
        double *prms; //array of variational parameters
        void thermalize(int); //termalization process: generates the population and thermalizes up to a time t_time
        void diffuse(); //starts the diffusion process
        double get_d(); //get the value of d
        double get_n(); //get the value of n
        void randomTest();

    private:
        RandomNumbers& rng; //we need the rand class to generate the random numbers
        bool thermalized; //tells me if the system is thermalized or not
        int walkers_num; //number of walkers
        int mesh_dim; //dimension of the mesh
        double dtau; // imaginary time step
        double n; //numerator: contains the sum w_i*EL
        double d; //denominator: contains the sum w_i
        double *population; //contains walker population
        double *weight; //contains associated weights
        double (*local_energy)(double, double *); //pointer to the local energy function
        double (*pseudo_force)(double, double *); //pointer to the pseudo force function
        double gauss_random_generator(double); //gaussian distributed random generator

};

RandomNumbers::RandomNumbers() {
#ifdef DEBUG //if the define DEBUG is present, compile this, so that the seed is always the same
    srand(0);
#else //otherwise, compile with the seed given by the present time
    srand(time(NULL)); //seed
#endif
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);
    gsl_rng_set(r, rand());
}
//------------------------------------------------------------------------------

RandomNumbers::~RandomNumbers() {
        gsl_rng_free(r);
}
//------------------------------------------------------------------------------

double RandomNumbers::operator() () {
    //calling rng() I get a random number
    return gsl_rng_uniform_pos(r);
}
//------------------------------------------------------------------------------

Vmc::Vmc(double (*l)(double, double *), double (*p)(double, double*), int wdim, int mdim, double step) : rng(*(new RandomNumbers())) {
    walkers_num = wdim;
    mesh_dim = mdim;
    dtau = step;
    population = new double[walkers_num];
    weight = new double[walkers_num];
    thermalized = false;
    local_energy = l;
    pseudo_force = p;

}
//------------------------------------------------------------------------------

Vmc::~Vmc() {
    delete population;
    delete weight;
    delete &rng;
}
//------------------------------------------------------------------------------

double Vmc::get_d() {
    return d;
}
//------------------------------------------------------------------------------

double Vmc::get_n() {
    return n;
}
//------------------------------------------------------------------------------

double local_energy(double x, double *alpha) { //local energy for the harmonic oscillator
    return alpha[0] + x * x * (1./2. - 2. * alpha[0] * alpha[0]);
}
//------------------------------------------------------------------------------

double pseudo_force(double x, double *alpha) { //pseudo force for the harmonic oscillator
    return - ( 2. * alpha[0] * x ) * DTAU;
}
//------------------------------------------------------------------------------

double Vmc::gauss_random_generator(double tau) { //given two uniformly distributed rng(), generates a gaussian distributed rng()
    return sqrt( - 2.*log(rng()) ) * cos(2.*M_PI*rng()) * sqrt(tau);
}
//------------------------------------------------------------------------------

void Vmc::thermalize(int t_time) { //thermalization process
    thermalized = false; //not thermalized

    for(int i = 0; i < walkers_num; i++ ) population[i] = ( .5 - rng() ) * MDIM; //generate a uniform population
    for(int k = 0; k <= t_time; k++) diffuse(); //diffuse for thermalization

    n = 0; //set parameters to be the good ones
    d = 0;

    thermalized = true; //now the system is thermalized and the system is ready for convergence
}
//------------------------------------------------------------------------------

void Vmc::diffuse() {

    double en_old, en_new;

    for( int i = 0; i < walkers_num; i++ ) {
        //compute the local energy for the present configuration
        en_old = local_energy(population[i], prms);

        //sum to each walker a random number gaussian distributed with the pseudo force component
        population[i] += ( gauss_random_generator(dtau) + pseudo_force(population[i],prms) );

        //re-compute the local energy
        en_new = local_energy(population[i],prms);

            if( thermalized ) { //if thermalization is achived, compute stuff

                //compute weight
                weight[i] = exp( - ( (en_new + en_old)/2. - Et ) * dtau );

                //compute energy terms, numerator and denominator
                n += weight[i] * local_energy(population[i],prms);
                d += weight[i];
            }
    }

}
//------------------------------------------------------------------------------

int main() {

    int k = 0, j = 0, therm, test;
    double sol = 100., prm;
    double *prms = new double[EDIM];
    FILE *pt;

    Vmc *V;
    V = new Vmc(&local_energy, &pseudo_force, WALKERS, MDIM, DTAU);
    V->prms = prms;

    pt = fopen("energy.txt", "w");

    printf("\n\t#----------------------------------");
    printf("\n\t#Alpha \t #Energy \t #Error\n");

    V->Et = 0.1;


    for( prms[0] = A_MIN; prms[0] <= A_MAX; prms[0] += A_STEP ) { //cycle over the values of the variational parameter

       V->thermalize((int)(TMIN/DTAU)); //thermalize

        do { //ready to badass cycle
            V->diffuse(); //diffuse
            k++;

                if( k%500 == 0 ) V->Et = V->get_n()/V->get_d(); //redefine trial energy every 500 steps

            } while( k < (TIME/DTAU) );

        k = 0;

        fflush(stdout);
        //print stuff
        printf("\t%g \t %.4f \t %.4f\n", prms[0], V->get_n()/V->get_d(),
                        fabs(prms[0] + 1./(8. * prms[0]) - prms[0]/2. - V->get_n()/V->get_d()) );
        fprintf(pt, "%g \t %g\n", prms[0], V->get_n()/V->get_d());

        if(sol > V->get_n()/V->get_d() ) {
            sol = V->get_n()/V->get_d(); //save the minimum value of the energy
            prm = prms[0];
        }

    }

    fclose(pt);
    delete prms;

    printf("\t#----------------------------------\n");
    printf("\n\t#SOLUTION: alfa = %g, energy = %g\n", prm, sol);
}
