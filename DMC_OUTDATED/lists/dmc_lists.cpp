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
#define MDIM (20) //mesh dimension
#define WALKERS (500) //number of walkers
#define NFIX (500.) //mean number of walkers I want at convergence
#define PARTICLE_NUM (2) //number of particles
#define LAMBDA (0.25) //width of the harmic well
#define BETA (358.93) //constant in front of LJ potential
#define CONVERSION (0.297) //conversion factor for energy
#define NPARAMS (2) //number of variational parameters
#define DTAU (5.e-5) //imaginary time step
#define TIMEV (30.) //propagation time for VMC
#define TMINV (0.) //minimum time after the energy sampling begins
#define TIMED (50.) //propagation time for DMC
#define DIFFUSION (1) //to avoid thermalization process
#define VARIATIONAL (0) //to thermalize
#define STEPa (0.05) //step on the first variational parameter
#define STEPb (0.05) //step on the second variational parameter
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
};


class MonteCarlo {
    public:
        MonteCarlo(int, int, int, double); //constructor

        double Et; //trial energy
        double *prms; //array of variational parameters
        int mean; //mean number of walkers
        double err; //error on the variational esitmation
        void thermalize(int, int); //termalization process: generates the population and thermalizes up to a time t_time
        void diffuse_no_branch_cycle(int); //cycle to obtain an energy from the Variational Monte Carlo
        void diffuse_branch_cycle(); //cycle to obtain an energy from the Diffusion Monte Carlo
        double get_d(); //gets the value of d
        double get_n(); //gets the value of n
        double get_nErr(); //gets the value of nErr
        void mean_walkers(); //gets the mean value of walkers

    private:
        RandomNumbers rng; //we need the rand class to generate the random numbers
        list<Walker> population; //list of walkers
        bool thermalized; //tells me if the system is thermalized or not
        int particles_num;
        int walkers_num; //number of walkers
        int mesh_dim; //dimension of the mesh
        double dtau; // imaginary time step
        double n; //numerator: contains the sum w_i*EL
        double d; //denominator: contains the sum w_i
        double nErr; //numerator for the calculation of the error
        double local_energy(double *, double *, double *); //pointer to the local energy function
        double pseudo_force(double *, double *, double *, int, int); //pointer to the pseudo force function
        double gauss_random_generator(double); //gaussian distributed random generator
        void diffuseBranch(); //one step of the diffusione process for the Diffusion Monte Carlo
        void diffuseNoBranch(); //one step of the diffusion process for the Variational Monte Carlo
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

MonteCarlo::MonteCarlo(int wdim, int mdim, int pdim, double step) : rng() {
    walkers_num = wdim;
    particles_num = pdim;
    mesh_dim = mdim;
    dtau = step;
    thermalized = false;
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

double MonteCarlo::local_energy(double *Rx, double *Ry, double *Rz) {
    double a = prms[0];
    double b = prms[1];
    double b5 = b * b * b * b * b, l = LAMBDA * LAMBDA * LAMBDA * LAMBDA * LAMBDA * LAMBDA;
    double sum = 0., tot = 0.;
    double r, r_ij, r_ij7, r_ij6;

    //psi_T = f*g dove f = psi_HO e g = psi_Jastrow = psi_J
    //nabla^2psi_T/psi_T = nabla^2f/f + nabla^2g/g + 2 nablaf/f*nablag/g

    for( int i = 0; i < particles_num; i++ ) {
        r = sqrt( Rx[i]*Rx[i] + Ry[i]*Ry[i] + Rz[i]*Rz[i] );
        sum += ( 3. * a + r * r * (1./2. - 2. * a * a) ); //nabla^2f/f+V_HO

        for( int j = 0; j < particles_num; j++ ) {
            if( j != i ) {
                r_ij = sqrt( (Rx[j] - Rx[i]) * (Rx[j] - Rx[i]) + (Ry[j] - Ry[i]) * (Ry[j] - Ry[i]) + (Rz[j] - Rz[i]) * (Rz[j] - Rz[i]) );
                r_ij6 = r_ij * r_ij * r_ij * r_ij * r_ij * r_ij;
                r_ij7 = r_ij6 * r_ij;

                sum += - (25./8.) * b5 * b5 / (r_ij6 * r_ij6) + 5. * b5 / r_ij7 ; //nabla^2g/g
                sum -= (5./2.) * b5 * a * ( Rx[i] * (Rx[j] - Rx[i]) + Ry[i] * (Ry[j] - Ry[i]) + Rx[i] * (Rz[j] - Rz[i]) ) / r_ij7;
                 //nablaf/f*nablag/g
            }
        }
    }

    for( int i = 0; i < particles_num; i++ ) { //V_LJ
        for( int j = 0; j < i; j++ ) {
            r_ij = sqrt( (Rx[i] - Rx[j]) * (Rx[i] - Rx[j]) + (Ry[i] - Ry[j]) * (Ry[i] - Ry[j]) + (Rz[i] - Rz[j]) * (Rz[i] - Rz[j]) );
            r_ij6 = r_ij * r_ij * r_ij * r_ij * r_ij * r_ij;
            sum += BETA * ( l * l / (r_ij6 * r_ij6) - l / r_ij6 );
        }
    }

    return sum;
}
//------------------------------------------------------------------------------

double MonteCarlo::pseudo_force(double *Rx, double *Ry, double *Rz, int direction, int particle) {
    double a = prms[0];
    double b = prms[1];
    double b5 = b * b * b * b * b, sum = 0., r_ij, r_ij7;

    if( direction == 0 ) {
        for( int i = 0; i < particles_num; i++ ) {
            if( i != particle ) {
                r_ij = sqrt( (Rx[particle] - Rx[i]) * (Rx[particle] - Rx[i]) + (Ry[particle] - Ry[i]) *
                        (Ry[particle] - Ry[i]) + (Rz[particle] - Rz[i]) * (Rz[particle] - Rz[i]) );
                r_ij7 = r_ij * r_ij * r_ij * r_ij * r_ij * r_ij * r_ij;
                sum -= (5./2.) * b5 * (Rx[i] - Rx[particle]) / r_ij7;
            }
        }
        return ( - 2. * a * Rx[particle] + sum ) * DTAU;

    } else if (direction == 1 ) {
        for( int i = 0; i < particles_num; i++ ) {
            if( i != particle ) {
                r_ij = sqrt( (Rx[particle] - Rx[i]) * (Rx[particle] - Rx[i]) + (Ry[particle] - Ry[i]) *
                        (Ry[particle] - Ry[i]) + (Rz[particle] - Rz[i]) * (Rz[particle] - Rz[i]) );
                r_ij7 = r_ij * r_ij * r_ij * r_ij * r_ij * r_ij * r_ij;
                sum -= (5./2.) * b5 * (Ry[i] - Ry[particle]) / r_ij7;
            }
        }
        return ( - 2. * a * Ry[particle] + sum ) * DTAU;

    } else if( direction == 2 ) {
        for( int i = 0; i < particles_num; i++ ) {
            if( i != particle ) {
                r_ij = sqrt( (Rx[particle] - Rx[i]) * (Rx[particle] - Rx[i]) + (Ry[particle] - Ry[i]) *
                        (Ry[particle] - Ry[i]) + (Rz[particle] - Rz[i]) * (Rz[particle] - Rz[i]) );
                r_ij7 = r_ij * r_ij * r_ij * r_ij * r_ij * r_ij * r_ij;
                sum -= (5./2.) * b5 * (Rz[i] - Rz[particle]) / r_ij7;
            }
        }
        return ( - 2. * a * Rz[particle] + sum ) * DTAU;

    } else {
        printf("Not a direction");
        exit(0);
    }
}
//------------------------------------------------------------------------------

double MonteCarlo::gauss_random_generator(double tau) { //given two uniformly distributed rng(), generates a gaussian distributed rng()
    return sqrt( - 2.*log(rng()) ) * cos(2.*M_PI*rng()) * sqrt(tau);
}
//------------------------------------------------------------------------------

void MonteCarlo::thermalize(int method, int t_time) { //thermalization process

    if(method == VARIATIONAL) {

        //list<Walker>::iterator it; //define an iterator over a list of walkers
        thermalized = false; //not thermalized

        //now we are filling the list pushing back
        for( int i = 0; i < walkers_num; i++ ) { //runs over the desired number of starting walkers
            Walker new_walker(particles_num); //define a new class walker

            for(int j = 0; j < particles_num; j++ ) {
                new_walker.x[j] = ( .5 - rng() ) * MDIM; //uniform population between -5 and 5
                new_walker.y[j] = ( .5 - rng() ) * MDIM;
                new_walker.z[j] = ( .5 - rng() ) * MDIM;
            }

            population.push_back(new_walker); //push the walker at the new starting position of the filled list
        }

        for(int k = 0; k <= t_time; k++) diffuseNoBranch(); //diffuses for thermalization
    }

    n = 0.; //set parameters to be the good ones
    d = 0.;
    nErr = 0.;

    thermalized = true; //now the system is thermalized and the system is ready for convergence
}
//------------------------------------------------------------------------------

void MonteCarlo::diffuseNoBranch() { //variational MC

    double en_old, en_new;
    list<Walker>::iterator it; //define an iterator over a list of walkers

    for ( it = population.begin(); it != population.end(); it++ ) { //runs over all the walkers

        //compute the local energy for the present configuration
        en_old = local_energy(it->x, it->y, it->z);

        //sum to each walker a random number gaussian distributed with the pseudo force component
        for( int j = 0; j < particles_num; j++ ) {
            it->x[j] += ( gauss_random_generator(dtau) + pseudo_force(it->x, it->y, it->z, 0, j) );
            it->y[j] += ( gauss_random_generator(dtau) + pseudo_force(it->x, it->y, it->z, 1, j) );
            it->z[j] += ( gauss_random_generator(dtau) + pseudo_force(it->x, it->y, it->z, 2, j) );
        }

        //re-compute the local energy
        en_new = local_energy(it->x, it->y, it->z);

        if( thermalized ) { //if thermalization is achived, compute stuff

            //compute weight
            it->weight = exp( - ( (en_new + en_old)/2. - Et ) * dtau );

            //compute energy terms, numerator, denominator and error numerator
            n += it->weight * local_energy(it->x, it->y, it->z);
            d += it->weight;
            nErr += it->weight * local_energy(it->x, it->y, it->z) * local_energy(it->x, it->y, it->z);
        }
    }
}
//------------------------------------------------------------------------------

void MonteCarlo::diffuseBranch() { //diffusion MC

    double en_old, en_new;
    int M = 0, k = 0, dead_counter = 0, spawning_counter = 0, eff_weight;
    list<Walker>::iterator it; //define an iterator over a list of walkers

    for ( it = population.begin(); it != population.end(); it++ ) { //runs over all the walkers

        //compute the local energy for the present configuration
        en_old = local_energy(it->x, it->y, it->z);

        //sum to each walker a random number gaussian distributed with the pseudo force component
        for( int j = 0; j < particles_num; j++ ) {
            it->x[j] += ( gauss_random_generator(dtau) + pseudo_force(it->x, it->y, it->z, 0, j) );
            it->y[j] += ( gauss_random_generator(dtau) + pseudo_force(it->x, it->y, it->z, 1, j) );
            it->z[j] += ( gauss_random_generator(dtau) + pseudo_force(it->x, it->y, it->z, 2, j) );
        }

        //re-compute the local energy
        en_new = local_energy(it->x, it->y, it->z);

        //compute weight
        it->weight = exp( - ( (en_new + en_old)/2. - Et ) * dtau );
        eff_weight = exp( - ( (en_new + en_old)/2. - (Et + 1./dtau * log(NFIX/walkers_num)) ) * dtau );

        //compute molteplicity
        it->molteplicity = (int)( eff_weight + rng() );

        //compute energy terms, numerator, denominator and error numerator
        n += it->weight * local_energy(it->x, it->y, it->z);
        d += it->weight;
        nErr += it->weight * local_energy(it->x, it->y, it->z) * local_energy(it->x, it->y, it->z);
    }

    //list<Walker>::iterator end = population.end();

    for( it = population.begin(); it != population.end(); ) {
        if( it->molteplicity == 0 ) { //if molteplicity is 0 delete the walker
            it = population.erase(it); //using function erase, which returns the iterator pointing to this position
        } else {
            ++it;
        }
    }

    //end = population.end();

    for( it = population.begin(); it != population.end(); it++) {
       if(it->molteplicity > 1) { //if molteplicity is biggen than 1
            for( int j = 0; j < it->molteplicity; j++ ) { //copy molteplicity-times the walker
                population.push_back(*it); //push the walker at the new starting position of the filled list
                //*it is the dereferenced walker, which IS the class in the list
            }
        }
    }

    walkers_num = population.size(); //give me the number of walkers

    //printf("\tWalkers number: %d, \n",walkers_num);//dead = %d, spawned = %d\n",walkers_num, dead_counter, spawning_counter);
    fflush(stdout);
}
//------------------------------------------------------------------------------

void MonteCarlo::mean_walkers() {
    mean += walkers_num;
}
//------------------------------------------------------------------------------

void MonteCarlo::diffuse_no_branch_cycle(int exit_print) { //cycle needed to obtain an energy from a fixen variational parameter

    int k = 0, blocks = 0;
    double squareE = 0., Esquare = 0., rel_err;

    //start thermalization process
    //generate a new homogeneous population and set summation parameters to 0
    thermalize( VARIATIONAL, (int)(TMINV/DTAU) );

     do {
        diffuseNoBranch(); //diffusion process for VMC: no branching done
        if( k%1000 == 0 ) {
            Et = n/d; //redefine trial energy every 200 steps as the ratio n/d
            squareE += nErr/d;
            Esquare += Et*Et;
            blocks++;
            fflush(stdout);
            printf("\t#Cycles: %d/%d\r", k, (int)(TIMEV/DTAU));
        }
        k++; //increment the counter
     } while( k < (TIMEV/DTAU) ); //redo until the total diffusion time has been reached


    //print interesting stuff
    if( exit_print == 1 ) {
        err = sqrt(squareE - Esquare)/blocks;
        rel_err = err/(n/d) * 100;
        fflush(stdout);
        printf("\r\t#-----------------------------------------\n");
        printf("\t%g \t %g \t %.3f \t %.3f (%.2f%%)\n", prms[0], prms[1], n/d, err, rel_err );
    } else {
        printf("\r\t#-----------------------------------------\n");
    }
}
//------------------------------------------------------------------------------

void MonteCarlo::diffuse_branch_cycle() {

    int k = 0;
    double squareE = 0., Esquare = 0., blocks = 0.;
    FILE* pt;
    FILE* pt2;
    FILE *pt3;

    pt = fopen("walkers.txt", "w");
    pt2 = fopen("energy.txt","w");
    pt3 = fopen("energy_error.txt","w");

    thermalize(DIFFUSION, 0); //needed to reset n and d to 0
    mean = 0;

    do {
        fprintf(pt,"%d\t%d\n", k, walkers_num);
        diffuseBranch(); //diffuse

        if( k%1000 == 0 ) {
            Et = n/d; //redefine trial energy every 100 steps
            squareE += nErr/d;
            Esquare += Et*Et;
            blocks++;
            fprintf(pt2,"%d\t%g\n",k,n/d);
            printf("\t#(-) Cycles: %d/%d\r", k, (int)(TIMED/DTAU));
            fflush(stdout);
        }

        if( k%1000 == 0 ) {
            fprintf(pt3,"%d\t%g\t%g\n", k, n/d, sqrt(squareE - Esquare)/blocks);
        }

        k++;

        mean_walkers();
    } while( k < (TIMED/DTAU) );

    fclose(pt);
    fclose(pt2);
    fclose(pt3);

    printf("\t#Energy: %.3f +- %.3f [nat u]\n", n/d, sqrt(squareE - Esquare)/blocks);
    printf("\t#Mean number of walkers = %d\n", (int)(1.*(mean)/k));
    printf("\t#----------------------------------\n");
}
//------------------------------------------------------------------------------


int main() {

    int k = 0, j = 0, therm, test;
    double sol = 100., prm1, prm2;
    double *prms = new double[NPARAMS];
    FILE *pt;

    MonteCarlo *M;
    M = new MonteCarlo(WALKERS, MDIM, PARTICLE_NUM, DTAU);
    M->prms = prms;

    printf("\n\t#START VARIATIONAL MONTE CARLO");
    printf("\n\t#-----------------------------------------");
    printf("\n\t#Alpha \t#Beta \t#Energy    #Error on E\n");


    for( M->prms[0] = 0.7; M->prms[0] <= 0.7; M->prms[0] += STEPa ) {
        for( M->prms[1] = LAMBDA; M->prms[1] <= LAMBDA; M->prms[1] += STEPb ) {
            //cycle over the values of the variational parameter

            M->Et = 1.5;
            M->diffuse_no_branch_cycle(1);

            if(sol > M->get_n()/M->get_d()) {
                sol = M->get_n()/M->get_d(); //save the minimum value of the energy
                prm1 = prms[0];
                prm2 = prms[1];
            }
        }
    }

    printf("\t#-----------------------------------------\n");
    printf("\n\t#SOLUTION\n");
    printf("\t#-----------------------------------------\n");
    printf("\t#Alpha = %g +- %g\n", prm1, STEPa);
    printf("\t#Beta = %g +- %g\n", prm2, STEPb);
    printf("\t#Energy = %.3f +- %.3f\n",sol, M->err);

    printf("\t#-----------------------------------------\n");

    printf("\n\t#GENERATING POPULATION FOR DMC...\n");
    printf("\r\t#-----------------------------------------\n");

    M->prms[0] = prm1;
    M->prms[1] = prm2;
    M->Et = sol;

    M->diffuse_no_branch_cycle(0);

    printf("\n\t#START DIFFUSION MONTE CARLO\n");
    printf("\r\t#-----------------------------------------\n");

    M->diffuse_branch_cycle();

   return 0;

}
