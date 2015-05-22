#include <iostream>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <gsl/gsl_rng.h>
#include <time.h>
#include <list>

using namespace std;

//#define DEBUG //generates always to same sequence of random numbers in case debugging is needed

//------------------------------------------MACROS-----------------------------------------
#define MDIM (10.) //mesh dimension
#define WALKERS (1000) //we generate a starting configuration for the DMC with WALKERS walkers
#define PARTICLE_NUM (10) //number of particles
#define NPARAMS (2) //number of variational parameters
#define STEPa (0.05) //step on the first variational parameter
#define STEPb (0.05) //step on the second variational parameter
#define STARTLEN (1.) //starting length for metropolis
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


class Walker {
    public:
        Walker(int); //constructor: needs the particles number
        Walker(const Walker &new_walker); //copy constructor
        ~Walker(); //destructor
        Walker& operator= (const Walker &obj); //= operator overloading

        double *x; //array of dimension particle_num: x coordinate
        double *y; //array of dimension particle_num: y coordinate
        double *z; //array of dimension particle_num: z coordinate
        double weight; //statistical weight associated with the walker
        int molteplicity; //molteplicity associated with the walker
        double en_old; //remebers the energy of the previous configuration
};

class MetropolisMC {
    public:
        MetropolisMC(int, int, double (*)(double *, double *, double *, double *, int),
                     double (*)(double *, double *, double *, double *, double *, double *, double *, int)); //constructor
        ~MetropolisMC(); //destructor
        double *prms; //array of variational parameters
        double err; //error on the variational esitmation
        double energy; //solution provided by metropolis algorithm
        void thermalize(); //termalization process: generates the population and thermalizes up to a time t_time
        void metropolis_step(); //one step of the metropolis algorithm
        void metropolis(int, FILE*, int); //metropolis algorithm

    private:
        RandomNumbers rng; //we need the rand class to generate the random numbers
        Walker rw; //random walker for the present population
        Walker rw_old; //remembers the old population
        bool thermalized; //tells me if the system is thermalized or not
        int particles_num; //number of particles in the system
        int mesh_dim; //dimension of the mesh
        int t_loops; //number of thermalization loops we want to do
        int p_loops; //number of production loops we want to do
        int m_step; //counts the number of matropolis steps
        int collection_steps; //every time you reach collection_steps, save the configuration
        double len; //typical dimension of the displacement
        int accepted; //number of accepted configurations
        double a_rate; //acceptance rate
        double n; //accumulator for the energy
        double nErr; //accumulator for the error on the energy
        double (*local_energy)(double *, double *, double *, double *, int); //pointer to the local energy function
        double (*acceptance_prob)(double *, double *, double *, double *, double *, double *, double *, int); //acceptance probability
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

Walker::Walker(int particle_num) {
    x = new double[particle_num];
    y = new double[particle_num];
    z = new double[particle_num];
    en_old = 0;
}
//------------------------------------------------------------------------------

Walker::Walker(const Walker &new_walker) { //copy constructor
// I'm building a way to teach the compiler how to copy correctly a class of kind Walker

    x = new double[PARTICLE_NUM]; //allocate the pointer members
    y = new double[PARTICLE_NUM];
    z = new double[PARTICLE_NUM];

    for( int i = 0; i < PARTICLE_NUM; i++ ) { //iterate over the number of particles
        x[i] = new_walker.x[i]; //and copy member by member
        y[i] = new_walker.y[i];
        z[i] = new_walker.z[i];
    }

    molteplicity = 1; //then set m = 1
    weight = 0.; //and statistical weight 0
    en_old = 0.;
}
//------------------------------------------------------------------------------

Walker& Walker::operator= (const Walker &obj) {

    if (this == &obj)
        return *this;

    // do the copy
    for( int i = 0; i < PARTICLE_NUM; i++ ) {
        x[i] = obj.x[i];
        y[i] = obj.y[i];
        z[i] = obj.z[i];
    }

    molteplicity = 1;
    weight = 0.;
    en_old = 0.;

    // return the existing object
    return *this;
}
//------------------------------------------------------------------------------

Walker::~Walker() {
    delete[] x;
    delete[] y;
    delete[] z;
}
//------------------------------------------------------------------------------

MetropolisMC::MetropolisMC(int mdim, int pdim, double (*l)(double *, double *, double *, double *, int),
                           double (*denst)(double *, double *, double *, double *, double *, double *, double *, int))
                            : rng(), rw(pdim), rw_old(pdim) {

    thermalized = false; //tells me if the system is thermalized or not
    particles_num = pdim; //number of particles in the system
    mesh_dim = mdim; //dimension of the mesh
    t_loops = 50000; //number of thermalization loops we want to do
    p_loops = 50000; //number of production loops we want to do
    collection_steps = 1000; //if k%collection_steps == 0 we collect interesting values
    m_step = 0; //counts the number of matropolis steps
    len = STARTLEN; //typical dimension of the displacement
    accepted = 0; //number of accepted configurations
    a_rate = 0; //acceptance rate

    local_energy = l; //local energy is incremented at each step in n
    acceptance_prob = denst; //defines the probability to accept the change in configuration
}
//------------------------------------------------------------------------------

double denst(double *x, double *y, double *z, double *x_old, double *y_old, double *z_old, double *prms, int particles_num) {
    double sum = 0;
    double r2;
    double a = prms[0];

    for( int i = 0; i < particles_num; i++ ) {
        r2 = ( x[i]*x[i] - x_old[i]*x_old[i] ) + ( y[i]*y[i] - y_old[i]*y_old[i] ) + ( z[i]*z[i] - z_old[i]*z_old[i] );
        sum += r2;
    }

    return exp(- 2. * a * sum);
}
//------------------------------------------------------------------------------

double local_energy(double *Rx, double *Ry, double *Rz, double *var_param, int particle_num) {
    //local energy for n 1d uncoupled harmonic oscillators
    double a = var_param[0];
    double sum = 0.;
    double r;

    for( int i = 0; i < particle_num; i++ ) {
        r = Rx[i]*Rx[i] + Ry[i]*Ry[i] + Rz[i]*Rz[i] ;
        sum += 3.*a - 2. * a * a * r + 1./2.*r; //nabla^2f/f+V_HO
    }

    return sum;
}
//------------------------------------------------------------------------------

void MetropolisMC::thermalize() { //thermalization process
    /*
     * we prepare the particles in a random position between -MDIM/2 and MDIM/2
     * we also save the value for the old configuration
     */
    for(int j = 0; j < particles_num; j++ ) {
        rw.x[j] = rw_old.x[j] = ( .5 - rng() ) * sqrt(prms[0]/particles_num);
        rw.y[j] = rw_old.y[j] = ( .5 - rng() ) * sqrt(prms[0]/particles_num);
        rw.z[j] = rw_old.z[j] = ( .5 - rng() ) * sqrt(prms[0]/particles_num);
    }

    rw.en_old = local_energy(rw.x, rw.y, rw.z, prms, particles_num);

    for(int k = 0; k <= t_loops; k++) {
        metropolis_step(); //do t_loops steps for thermalization: no energy is saved

        /* adaptive displancement length
         * the starting length len is not necessarly the good one
         * during thermalization we try to keep the accpetance rate to 0.5, modifying len
         */
        if( k < t_loops ) { //ONLY during thermalization
            if( k%collection_steps == 0 ) { //every collection_steps do the check
                if(a_rate < 0.45 ) { //acceptance rate is too small
                    len *= 0.9; //probably len is too big: reduce it!
                } else if(a_rate > 0.55) { //acceptance rate is too big
                    len *= 1.1; //probably len is too small: increment it!
                }
            }
        }
    }

    //set parameters to be the good ones
    n = 0.;
    nErr = 0.;
    accepted = 0;
    m_step = 0;
    len = STARTLEN;
    a_rate = 0.;
    thermalized = true; //now the system is thermalized and the system is ready for convergence
}
//------------------------------------------------------------------------------

void MetropolisMC::metropolis_step() { //one step for variational MC

    double en_new, r1, r2, r3, psi_new, p, alpha, rn;

    /*
     * loop over all the particles, and change the position of only one of them at the time
     * this is the most efficient way, because less configurations are rejected
     * change the position of the particle
     * calculate acceptance probability
     * draw a uniform rn and choose if the new configuration is accepted or not
     */

    for( int j = 0; j < particles_num; j++ ) {

        m_step++; //increment the number of steps

        //len is modified iteratively during the thermalization process
        r1 = 2. * len * ( rng() - 0.5 );
        r2 = 2. * len * ( rng() - 0.5 );
        r3 = 2. * len * ( rng() - 0.5 );

        //move the particle
        rw.x[j] += r1;
        rw.y[j] += r2;
        rw.z[j] += r3;

        //calculate acceptance probability
        alpha = acceptance_prob(rw.x, rw.y, rw.z, rw_old.x, rw_old.y, rw_old.z, prms, particles_num);

        if( rng() < alpha ) { //the configuration is accepted

            en_new = local_energy(rw.x, rw.y, rw.z, prms, particles_num); //compute the new energy
            rw.en_old = en_new; //save the new energy as the old one
            accepted++; //increment the accepted movements

            //if the configuration has been accepted I have to refresh the old config with the new one
            rw_old.x[j] = rw.x[j];
            rw_old.y[j] = rw.y[j];
            rw_old.z[j] = rw.z[j];

        } else { //configuration rejected

            //back to the old configuration!
            rw.x[j] -= r1;
            rw.y[j] -= r2;
            rw.z[j] -= r3;
            en_new = rw.en_old; //en_new is just the old one, no need to compute it again
        }

        a_rate = (1.*accepted)/(1.*m_step); //calculate the acceptance rate

        //re-compute the local energy

        if( thermalized ) { //if thermalization is achived, compute stuff
            n += en_new; //increment local energy
            nErr += en_new * en_new; //increment error
        }

    }
}
//------------------------------------------------------------------------------

void MetropolisMC::metropolis(int exit_print, FILE *walkers, int walkers_num) {
    //cycle needed to obtain an energy from a fixen variational parameter

    int k = 0;
    double rel_err;
    FILE *pt;
    int config_sample = int( (double)(p_loops)/(double)(walkers_num) );

    thermalize();

    pt = fopen("var_energy.txt","w");
    walkers = fopen("config.txt","w");

     do {
        metropolis_step(); //diffusion process for VMC: no branching done
        fprintf(pt,"%d\t%g\n",k,n/m_step);

        if( k%config_sample == 0 ) {
            for( int i = 0; i < particles_num; i++ ) {
                fprintf(walkers, "%g\t%g\t%g\n", rw.x[i], rw.y[i], rw.z[i]);
            }
        }

        if( k%collection_steps == 0 ) {
            fflush(stdout);
            printf("\t#Cycles: %d/%d\r", k, p_loops);
        }

        k++; //increment the counter

     } while( k < p_loops ); //redo until the total diffusion time has been reached

    fclose(pt);
    fclose(walkers);


    //print interesting stuff
    if( exit_print == 1 ) {
        err = sqrt((nErr/(m_step) - ((n/m_step)/(m_step))*((n/m_step)/(m_step)))/(m_step));
        rel_err = err/(n/m_step) * 100;
        fflush(stdout);
        printf("\r\t#-------------------------------------------------\n");
        printf("\t%g \t %g \t %.3f \t %.3f (%.2f%%)\n", prms[0], prms[1], (n/m_step), err, rel_err );
    } else {
        printf("\r\t#-------------------------------------------------\n");
    }

}
//------------------------------------------------------------------------------


int main() {

    int k = 0, j = 0, therm, test;
    double sol = 10000., prm1, prm2, err;
    double *prms = new double[NPARAMS];
    FILE *pt;
    FILE *walkers;
    FILE *config;

    MetropolisMC *M;

    M = new MetropolisMC(MDIM, PARTICLE_NUM, &local_energy, &denst);
    M->prms = prms;

    printf("\n\t#START VARIATIONAL MONTE CARLO");
    printf("\n\t#-------------------------------------------------");
    printf("\n\t#Alpha \t#Beta \t#Energy    #Error on E\n");


    M->prms[1] = 0.;
    for( M->prms[0] = 0.6; M->prms[0] <= 0.6; M->prms[0] += 0.1 ) {
        //for( M->prms[1] = LAMBDA; M->prms[1] <= LAMBDA; M->prms[1] += STEPb ) {
            //cycle over the values of the variational parameter

            M->metropolis(1, walkers, WALKERS);

            if( sol > M->energy + M->err ) {
                sol = M->energy; //save the minimum value of the energy
                prm1 = prms[0];
                prm2 = prms[1];
                err = M->err;
            }
        //}
    }

    printf("\t#-----------------------------------------\n");
    printf("\n\t#SOLUTION\n");
    printf("\t#-----------------------------------------\n");
    printf("\t#Alpha = %g +- %g\n", prm1, STEPa);
    printf("\t#Beta = %g +- %g\n", prm2, STEPb);
    printf("\t#Energy = %.3f +- %.3f\n",sol, err);

    printf("\t#-----------------------------------------\n");

    system("open config.txt");

    return 0;

}
