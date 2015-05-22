#include <iostream>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <gsl/gsl_rng.h>
#include <time.h>

using namespace std;

//#define DEBUG

//------------------------------------------MACROS-----------------------------------------
#define N (1) //number of particles
#define MDIM (10) //mesh dimension
#define WALKERS (1000) //number of walkers
#define DTAU (5.e-5) //imaginary time step
#define TIME (10.) //propagation time
#define TMIN (2.) //minimum time after the energy sampling begins
#define NPARAMS (1) //number of variational parameters
#define DIFFUSION (1) //to avoid thermalization process
#define VARIATIONAL (0) //to thermalize

#define A_MIN (0.3) //valore minimo per il parametro variazionale
#define A_MAX (0.55) //valore massimo per il parametro variazionale
#define A_STEP (0.05) //step per il parametro variazionale
#define EDIM (int)(A_MAX/A_STEP - A_MIN/A_STEP) //dimensione per l'array che conterrà le varie energie
//----------------------------------------END-MACROS---------------------------------------

typedef struct {
    double *x; //contains walker population for x coordinate
    double *y; //contains walker population for y coordinate
    double *z; // containes walker population for z coordinate
} population; //structure which contains population, weight, molteplicity for each dimension
//most general way because I can handle any dimension


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
        MonteCarlo(double (*)(double *, double *), double (*)(double *, double *), int, int, int, double); //constructor
        ~MonteCarlo(); //destructor

        double Et; //trial energy
        double *prms; //array of variational parameters
        int mean;
        void thermalize(int, int, int); //termalization process: generates the population and thermalizes up to a time t_time
        void diffuseNoBranch(); //starts the diffusion process for the Variational Monte Carlo
        void diffuseBranch(); //starts the diffusione process for the Diffusion Monte Carlo
        double get_d(); //get the value of d
        double get_n(); //get the value of n
        double get_nErr(); //get the value of nErr
        void mean_walkers();

    private:
        RandomNumbers& rng; //we need the rand class to generate the random numbers
        bool thermalized; //tells me if the system is thermalized or not
        int particles_num; //number of particles in the system
        int walkers_num; //number of walkers
        int mesh_dim; //dimension of the mesh
        double dtau; // imaginary time step
        double n; //numerator: contains the sum w_i*EL
        double d; //denominator: contains the sum w_i
        double nErr; //numerator for the calculation of the error
        population *particles; //struct which contains the population of walkers for each coordinate
        double *weight; //contains associated weights
        int *molteplicity; //contains associated molteplicity
        double (*local_energy)(double *, double *); //pointer to the local energy function
        double (*pseudo_force)(double *, double *); //pointer to the pseudo force function
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

MonteCarlo::MonteCarlo(double (*l)(double*, double *), double (*p)(double*, double*),
                       int pdim, int wdim, int mdim, double step) : rng(*(new RandomNumbers())) {
    particles_num = pdim;
    walkers_num = wdim;
    mesh_dim = mdim;
    dtau = step;
    particles = new population[particles_num];
    double *x = new double[walkers_num];
    double *y = new double[walkers_num];
    double *z = new double[walkers_num];
    weight = new double[walkers_num];
    particles->x = x;
    particles->y = y;
    particles->z = z;
    molteplicity = new int[walkers_num];
    thermalized = false;
    local_energy = l;
    pseudo_force = p;

}
//------------------------------------------------------------------------------

MonteCarlo::~MonteCarlo() {
    delete particles;
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

double local_energy(double *r, double *var_param) { //local energy for the harmonic oscillator
    double x = r[0];
    double a = var_param[0];

    return a + x * x * (1./2. - 2. * a * a);
}
//------------------------------------------------------------------------------

double pseudo_force(double *r, double *var_param) { //pseudo force for the harmonic oscillator
    double x = r[0];
    double a = var_param[0];

    return - ( 2. * a * x ) * DTAU;
}
//------------------------------------------------------------------------------

double MonteCarlo::gauss_random_generator(double tau) { //given two uniformly distributed rng(), generates a gaussian distributed rng()
    return sqrt( - 2.*log(rng()) ) * cos(2.*M_PI*rng()) * sqrt(tau);
}
//------------------------------------------------------------------------------

void MonteCarlo::thermalize(int method, int t_time, int dimension) { //thermalization process

    if(method == VARIATIONAL) {

        thermalized = false; //not thermalized

        for(int j = 0; j < particles_num; j++) {
            for(int i = 0; i < walkers_num; i++ ) {
                particles[j].x[i] = ( .5 - rng() ) * MDIM; //generate a uniform population for x

                if( dimension == 2 ) { //if the dimensionality is 2, fill x and y coordinates
                    particles[j].y[i] = ( .5 - rng() ) * MDIM; //generate a uniform population for y
                } else if( dimension == 3 ) { //if the dimensionality is 3, fill all the coordinates
                    particles[j].y[i] = ( .5 - rng() ) * MDIM; //generate a uniform population for y
                    particles[j].z[i] = ( .5 - rng() ) * MDIM; //generate a uniform population for z
                }
            }
        }

        for(int k = 0; k <= t_time; k++) diffuseNoBranch(); //diffuse for thermalization
    }

    n = 0; //set parameters to be the good ones
    d = 0;

    thermalized = true; //now the system is thermalized and the system is ready for convergence
}
//------------------------------------------------------------------------------

void MonteCarlo::diffuseNoBranch() {

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

/*
void MonteCarlo::diffuseBranch() {

    double en_old, en_new;
    int M = 0, k = 0;
    double *tmp = new double[M];

    for( int i = 0; i < walkers_num; i++ ) {

        //compute the local energy for the present configuration
        en_old = local_energy(population[i], prms);

        //sum to each walker a random number gaussian distributed with the pseudo force component
        population[i] += ( gauss_random_generator(dtau) + pseudo_force(population[i],prms) );

        //re-compute the local energy
        en_new = local_energy(population[i],prms);

        //compute weight
        weight[i] = exp( - ( (en_new + en_old)/2. - Et ) * dtau );
        molteplicity[i] = (int)( floor( weight[i] + rng() ) );

        M += molteplicity[i];

        n += weight[i] * local_energy(population[i],prms);
        d += weight[i];
        nErr += weight[i] * local_energy(population[i],prms) * local_energy(population[i],prms);
    }

    printf("\t#(-) Number of walkers: %d\r", M);
    fflush(stdout);

    for( int i = 0; i < walkers_num; i++ ) {
        for(int j = 0; j < molteplicity[i]; j++ ) {
            tmp[k] = population[i];
            k++;
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
} */

int main() {

    int k = 0, j = 0, therm, test;
    double sol = 100., prm, mean_walkers, squareE = 0., Esquare = 0., blocks = 0.;
    double *prms = new double[EDIM];
    FILE *pt;

    MonteCarlo *M;
    M = new MonteCarlo(&local_energy, &pseudo_force, N, WALKERS, MDIM, DTAU);
    M->prms = prms;

    M->thermalize(VARIATIONAL, 0, 3);


    return 0;
}

//int main() {

//        int k = 0, j = 0, therm, test;
//        double sol = 100., prm, mean_walkers, squareE = 0., Esquare = 0., blocks = 0.;
//        double *prms = new double[EDIM];
//        FILE *pt;

//        MonteCarlo *M;
//        M = new MonteCarlo(&local_energy, &pseudo_force, WALKERS, MDIM, DTAU);
//        M->prms = prms;

//        pt = fopen("energy.txt", "w");

//        printf("\n\t#START VARIATIONAL MONTE CARLO");
//        printf("\n\t#----------------------------------");
//        printf("\n\t#Alpha \t #Energy \t #Error\n");

//        M->Et = 0.1;


//        for( prms[0] = A_MIN; prms[0] <= A_MAX; prms[0] += A_STEP ) { //cycle over the values of the variational parameter

//           M->thermalize(VARIATIONAL, (int)(TMIN/DTAU)); //thermalize

//            do { //ready to badass cycle
//                M->diffuseNoBranch(); //diffuse
//                k++;

//                    if( k%200 == 0 ) M->Et = M->get_n()/M->get_d(); //redefine trial energy every 500 steps

//            } while( k < (TIME/DTAU) );

//            k = 0;

//            fflush(stdout);
//            //print stuff
//            printf("\t%g \t %.4f \t %.4f\n", prms[0], M->get_n()/M->get_d(),
//                            fabs(prms[0] + 1./(8. * prms[0]) - prms[0]/2. - M->get_n()/M->get_d()) );
//            fprintf(pt, "%g \t %g\n", prms[0], M->get_n()/M->get_d());

//            if(sol > M->get_n()/M->get_d() ) {
//                sol = M->get_n()/M->get_d(); //save the minimum value of the energy
//                prm = prms[0];
//            }

//        }

//        fclose(pt);
//        //delete prms;

//        printf("\t#----------------------------------\n");
//        printf("\n\t#SOLUTION\n");
//        printf("\t#----------------------------------\n");
//        printf("\t#Alpha = %g\n", prm);
//        printf("\t#Energy = %g\n",sol);
//        printf("\t#----------------------------------\n");

//        printf("\n\t#START DIFFUSION MONTE CARLO\n");
//        printf("\t#----------------------------------\n");

//        prms[0] = prm;
//        M->Et = sol;
//        k = 0;

//        M->thermalize(DIFFUSION, 0);
//        M->mean = 0;

//        do { //ready to badass cycle
//            M->diffuseBranch(); //diffuse
//            k++;

//                if( k%100 == 0 ) {
//                    M->Et = M->get_n()/M->get_d(); //redefine trial energy every 500 steps
//                    squareE += M->get_nErr()/M->get_d();
//                    Esquare += M->Et * M->Et;
//                    blocks++;
//                }

//            M->mean_walkers();

//        } while( k < (TIME/DTAU) );

//        printf("\t#Energy: %.3f +- %.3f [nat u]\n",
//               M->get_n()/M->get_d(), sqrt(squareE - Esquare)/blocks);
//        printf("\t#Mean number of walkers = %d\n", (int)(1.*(M->mean)/k));
//        printf("\t#----------------------------------\n");


//}
