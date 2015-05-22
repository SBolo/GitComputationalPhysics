#include <iostream>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <gsl/gsl_rng.h>
#include <time.h>

using namespace std;

//#define DEBUG //generates always to same sequence of random numbers in case debugging is needed

//------------------------------------------MACROS-----------------------------------------
#define MDIM (10) //mesh dimension
#define WALKERS (1000) //number of walkers
#define NFIX (2000) //mean number of walkers I want
#define DTAU (1.e-3) //imaginary time step
#define TIME (50.) //propagation time
#define TMIN (2.) //minimum time after the energy sampling begins
#define NPARAMS (1) //number of variational parameters
#define DIFFUSION (1) //to avoid thermalization process
#define VARIATIONAL (0) //to thermalize

#define A_MIN (0.3) //valore minimo per il parametro variazionale
#define A_MAX (0.4) //valore massimo per il parametro variazionale
#define THREADNUM (2) //numero di thread per il calcolo variazionale
/*
 * !!!!! N.B.:
 */
#define A_STEP (0.05) //step per il parametro variazionale
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


class MonteCarlo {
    public:

        MonteCarlo(double (*)(double, double *), double (*)(double, double *), int, int, double); //constructor
        ~MonteCarlo(); //destructor

        double Et; //trial energy
        double *prms; //array of variational parameters
        int mean;
        void thermalize(int, int); //termalization process: generates the population and thermalizes up to a time t_time
        void diffuseNoBranch(); //starts the diffusion process for the Variational Monte Carlo
        void diffuseBranch(); //starts the diffusione process for the Diffusion Monte Carlo
        double get_d(); //gets the value of d
        double get_n(); //gets the value of n
        double get_nErr(); //gets the value of nErr
        int get_wn(); //returns the number of walkers
        void mean_walkers(); //gets the mean value of walkers

    private:
        RandomNumbers& rng; //we need the rand class to generate the random numbers
        bool thermalized; //tells me if the system is thermalized or not
        int walkers_num; //number of walkers
        int mesh_dim; //dimension of the mesh
        double dtau; // imaginary time step
        double n; //numerator: contains the sum w_i*EL
        double d; //denominator: contains the sum w_i
        double *walkers; //struct which contains the population of walkers for each coordinate
        double *population;
        double *weight; //contains associated weights
        int *molteplicity; //contains associated molteplicity
        double nErr; //numerator for the calculation of the error
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

MonteCarlo::MonteCarlo(double (*l)(double, double *), double (*p)(double, double*),
                       int wdim, int mdim, double step) : rng(*(new RandomNumbers())) {
    walkers_num = wdim;
    mesh_dim = mdim;
    dtau = step;
    walkers = new double[walkers_num];
    population = new double[walkers_num];
    weight = new double[walkers_num];
    molteplicity = new int[walkers_num];
    thermalized = false;
    local_energy = l;
    pseudo_force = p;

}
//------------------------------------------------------------------------------

MonteCarlo::~MonteCarlo() {
    delete population;
    delete weight;
    delete &rng;
}
//------------------------------------------------------------------------------

double MonteCarlo::get_d() {
    return d;
}
//------------------------------------------------------------------------------

double MonteCarlo::get_n() {
    return n;
}
//------------------------------------------------------------------------------

double MonteCarlo::get_nErr() {
    return nErr;
}
//------------------------------------------------------------------------------

int MonteCarlo::get_wn() {
    return walkers_num;
}
//------------------------------------------------------------------------------

double local_energy(double x, double *var_param) { //local energy for the harmonic oscillator
    double a = var_param[0];

    return a + x * x * (1./2. - 2. * a * a);
}
//------------------------------------------------------------------------------

double pseudo_force(double x, double *var_param) { //pseudo force for the harmonic oscillator
    double a = var_param[0];

    return - ( 2. * a * x ) * DTAU;
}
//------------------------------------------------------------------------------

double MonteCarlo::gauss_random_generator(double tau) { //given two uniformly distributed rng(), generates a gaussian distributed rng()
    return sqrt( - 2.*log(rng()) ) * cos(2.*M_PI*rng()) * sqrt(tau);
}
//------------------------------------------------------------------------------

void MonteCarlo::thermalize(int method, int t_time) { //thermalization process

    if(method == VARIATIONAL) {

        thermalized = false; //not thermalized

        for(int i = 0; i < walkers_num; i++ ) population[i] = ( .5 - rng() ) * MDIM; //generates a uniform population
        for(int k = 0; k <= t_time; k++) diffuseNoBranch(); //diffuses for thermalization
    }

    n = 0; //set parameters to be the good ones
    d = 0;

    thermalized = true; //now the system is thermalized and the system is ready for convergence
}
//------------------------------------------------------------------------------

void MonteCarlo::diffuseNoBranch() { //variational MC

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

                //compute energy terms, numerator, denominator and error numerator
                n += weight[i] * local_energy(population[i],prms);
                d += weight[i];
            }
    }

}
//------------------------------------------------------------------------------

void MonteCarlo::diffuseBranch(bool thermalized) {

    double en_old, en_new, w_m;
    int M = 0, k = 0;

    for( int i = 0; i < walkers_num; i++ ) {

        //compute the local energy for the present configuration
        en_old = local_energy(population[i], prms);

        //sum to each walker a random number gaussian distributed with the pseudo force component
        population[i] += ( gauss_random_generator(dtau) + pseudo_force(population[i],prms) );

        //re-compute the local energy
        en_new = local_energy(population[i],prms);

        //compute weight
        //weight[i] = exp( - ( (en_new + en_old)/2. - Et ) * dtau );
        weight[i] = exp( - ( (en_new + en_old)/2. - (Et + 1./dtau * log(2000./walkers_num)) ) * dtau );

        molteplicity[i] = (int)( weight[i]  + rng() );
        //printf("weight = %g, exp(1/tau) = %g, walk/walk = %g \n", weight[i], exp(1./dtau), ((1.*WALKERS)/walkers_num) );

        M += molteplicity[i];

        if( thermalized == true ) {
            n += weight[i] * local_energy(population[i],prms);
            d += weight[i];
            nErr += weight[i] * local_energy(population[i],prms) * local_energy(population[i],prms);
        }
    }

    printf("\t#(-) Number of walkers: %d\r", M);
    fflush(stdout);

    double *tmp = new double[M];

    for( int i = 0; i < walkers_num; i++ ) {
        if( molteplicity[i] > 0 ) {
            for(int j = 0; j < molteplicity[i]; j++ ) {
                tmp[k] = population[i];
                k++;
            }
        } else if( molteplicity[i] == 0 ) {
            continue;
        }
    }

    walkers_num = M;

    delete population;
    delete molteplicity;
    delete weight;

    population = new double[walkers_num];
    molteplicity = new int[walkers_num];
    weight = new double[walkers_num];

    for(int i = 0; i < walkers_num; i++ ) {
        population[i] = tmp[i];
    }
}
//------------------------------------------------------------------------------

void MonteCarlo::mean_walkers() {
    mean += walkers_num;
}

int main() {

        int k = 0, j = 0, therm, test;
        double sol = 100., prm, mean_walkers, squareE = 0., Esquare = 0., blocks = 0.;
        double *prms = new double[EDIM];
        FILE *pt;

        MonteCarlo *M;
        M = new MonteCarlo(&local_energy, &pseudo_force, WALKERS, MDIM, DTAU);
        M->prms = prms;

        pt = fopen("energy.txt", "w");

        printf("\n\t#START VARIATIONAL MONTE CARLO");
        printf("\n\t#----------------------------------");
        printf("\n\t#Alpha \t #Energy \t #Error\n");

        M->Et = 0.1;


        for( prms[0] = A_MIN; prms[0] <= A_MAX; prms[0] += A_STEP ) { //cycle over the values of the variational parameter

           M->thermalize(VARIATIONAL, (int)(TMIN/DTAU)); //thermalize

            do { //ready to badass cycle
                M->diffuseNoBranch(); //diffuse
                k++;

                    if( k%200 == 0 ) M->Et = M->get_n()/M->get_d(); //redefine trial energy every 200 steps

            } while( k < (TIME/DTAU) );

            k = 0;

            fflush(stdout);
            //print stuff
            printf("\t%g \t %.4f \t %.4f\n", prms[0], M->get_n()/M->get_d(),
                            fabs(prms[0] + 1./(8. * prms[0]) - prms[0]/2. - M->get_n()/M->get_d()) );
            fprintf(pt, "%g \t %g\n", prms[0], M->get_n()/M->get_d());

            if(sol > M->get_n()/M->get_d() ) {
                sol = M->get_n()/M->get_d(); //save the minimum value of the energy
                prm = prms[0];
            }

        }


        printf("\t#----------------------------------\n");
        printf("\n\t#SOLUTION\n");
        printf("\t#----------------------------------\n");
        printf("\t#Alpha = %g\n", prm);
        printf("\t#Energy = %g\n",sol);
        printf("\t#----------------------------------\n");

        printf("\n\t#START DIFFUSION MONTE CARLO\n");
        printf("\t#----------------------------------\n");

        M->prms[0] = prm;
        M->Et = sol;
        k = 0;

        M->thermalize(VARIATIONAL, (int)(TMIN/DTAU)); //thermalize

         do { //ready to badass cycle
             M->diffuseNoBranch(); //diffuse
             k++;

                 if( k%200 == 0 ) M->Et = M->get_n()/M->get_d(); //redefine trial energy every 200 steps

         } while( k < (TIME/DTAU) );

         k = 0; //I use the resulting population of the variational MC as starting population for DMC

        M->thermalize(DIFFUSION, 0);
        M->mean = 0;

        pt = fopen("walkers.txt", "w");

        do { //ready to badass cycle
            fprintf(pt,"%d\t%d\n",k,M->get_wn());
            M->diffuseBranch(); //diffuse
            k++;

                if( k%100 == 0 ) {
                    M->Et = M->get_n()/M->get_d(); //redefine trial energy every 100 steps
                    squareE += M->get_nErr()/M->get_d();
                    Esquare += M->Et * M->Et;
                    blocks++;
                }

            M->mean_walkers();

        } while( k < (TIME/DTAU) );

        fclose(pt);

        printf("\t#Energy: %.3f +- %.3f [nat u]\n",
               M->get_n()/M->get_d(), sqrt(squareE - Esquare)/blocks);
        printf("\t#Mean number of walkers = %d\n", (int)(1.*(M->mean)/k));
        printf("\t#----------------------------------\n");


}
