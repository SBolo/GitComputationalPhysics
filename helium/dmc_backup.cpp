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
#define MDIM (20.) //mesh dimension
#define WALKERS (200) //num3ber of walkers
#define NFIX (200.) //mean number of walkers I want at convergence
#define PARTICLE_NUM (1) //number of particles
#define LAMBDA (0.25) //width of the harmic well
#define BETA (358.93) //constant in front of LJ potential
#define CONVERSION (0.297) //conversion factor for energy
#define NPARAMS (2) //number of variational parameters
#define DTAU (5.e-4) //imaginary time step
#define TIMEV (10.) //propagation time for VMC
#define TMINV (50.) //minimum time after the energy sampling begins
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
        bool bmflag;
        double bmnext;
        double en_old; //remember the local energy calculated in the previous step
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

    bmflag = false;
    en_old = 0.;
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

double MonteCarlo::gauss_random_generator(double tau) { //given two uniformly distributed rng(), generates a gaussian distributed rng()

    if(bmflag) {
        bmflag = false;
        return bmnext;
    } else {
        double r = rng();
        double theta = rng();
        bmnext = sqrt( - 2.*log(r) ) * cos(2.*M_PI*theta) * sqrt(tau);

        bmflag = true;
        return sqrt( - 2.*log(r) ) * sin(2.*M_PI*theta) * sqrt(tau);
    }

}
//------------------------------------------------------------------------------

double MonteCarlo::local_energy(double *Rx, double *Ry, double *Rz) {
    //local energy for n 1d uncoupled harmonic oscillators
    double a = prms[0];
    double b = prms[1];

    double rx_ij, ry_ij, rz_ij, r2, r_ij, r_ij7, r_ij6, lj6;
    double b5 = b * b * b * b * b, l = LAMBDA * LAMBDA * LAMBDA * LAMBDA * LAMBDA * LAMBDA;
    double b10 = b5*b5;
    double sum1 = 0., sum2 = 0., sum3 = 0., sum4 = 0., sum5 = 0., tot = 0.;

    for( int i = 0; i < particles_num; i++ ) {
        r2 =  Rx[i]*Rx[i] + Ry[i]*Ry[i] + Rz[i]*Rz[i];
        sum1 += r2; //nabla^2f/f+V_HO

//        for( int j = 0; j < particles_num; j++ ) {
//            if( j != i ) {
//                rx_ij = Rx[j] - Rx[i];
//                ry_ij = Ry[j] - Ry[i];
//                rz_ij = Rz[j] - Rz[i];

//                r_ij = sqrt( rx_ij * rx_ij + ry_ij * ry_ij + rz_ij * rz_ij );
//                r_ij6 = r_ij * r_ij * r_ij;
//                r_ij6 *= r_ij6;
//                r_ij7 = r_ij6 * r_ij;
//                lj6 = l/r_ij6;

//                sum2 += 1. / (r_ij6 * r_ij6);
//                sum3 += 1. / r_ij7;
//                sum4 += ( Rx[i] * rx_ij + Ry[i] * ry_ij + Rz[i] * rz_ij ) / r_ij7;
//            }

//            if( j < i ) {
//                sum5 += lj6 * (- 1. + lj6);
//            }
//        }
    }

    tot = 3. * a * particles_num + (0.5 - 2. * a * a ) * sum1;
//            - 25./8. * b10 * sum2 + 5. * b5 * sum3
//            - (5./2.) * b5 * a * sum4 + BETA * sum5;

    return tot;
}
//------------------------------------------------------------------------------


double MonteCarlo::pseudo_force(double *Rx, double *Ry, double *Rz, int direction, int particle) {
    //pseudo force for the harmonic oscillator
    double a = prms[0];
    double b = prms[1];
    double b5 = b * b * b * b * b, sum = 0., r_ij, r_ij7, rx_ij, ry_ij, rz_ij;

    if( direction == 0 ) {
//        for( int i = 0; i < particles_num; i++ ) {
//            if( i != particle ) {
//                rx_ij = Rx[particle] - Rx[i];
//                ry_ij = Ry[particle] - Ry[i];
//                rz_ij = Rz[particle] - Rz[i];

//                r_ij = sqrt( rx_ij * rx_ij + ry_ij * ry_ij + rz_ij * rz_ij );
//                r_ij7 = r_ij * r_ij * r_ij * r_ij * r_ij * r_ij * r_ij;
//                sum += rx_ij / r_ij7;
//            }
//        }
        return ( - 2. * a * Rx[particle] + (5./2.) * b5 * sum ) * DTAU;

    } else if (direction == 1 ) {
//        for( int i = 0; i < particles_num; i++ ) {
//            if( i != particle ) {
//                rx_ij = Rx[particle] - Rx[i];
//                ry_ij = Ry[particle] - Ry[i];
//                rz_ij = Rz[particle] - Rz[i];

//                r_ij = sqrt( rx_ij * rx_ij + ry_ij * ry_ij + rz_ij * rz_ij );
//                r_ij7 = r_ij * r_ij * r_ij * r_ij * r_ij * r_ij * r_ij;
//                sum += ry_ij / r_ij7;
//            }
//        }
        return ( - 2. * a * Ry[particle] + (5./2.) * b5 * sum ) * DTAU;

    } else if( direction == 2 ) {
//        for( int i = 0; i < particles_num; i++ ) {
//            if( i != particle ) {
//                rx_ij = Rx[particle] - Rx[i];
//                ry_ij = Ry[particle] - Ry[i];
//                rz_ij = Rz[particle] - Rz[i];

//                r_ij = sqrt( rx_ij * rx_ij + ry_ij * ry_ij + rz_ij * rz_ij );
//                r_ij7 = r_ij * r_ij * r_ij * r_ij * r_ij * r_ij * r_ij;
//                sum += rz_ij / r_ij7;
//            }
//        }
        return ( - 2. * a * Rz[particle] + (5./2.) * b5 * sum ) * DTAU;

    } else {
        printf("Not a direction");
        exit(0);
    }
}
//------------------------------------------------------------------------------

void MonteCarlo::thermalize(int method, int t_time) { //thermalization process

    FILE *config;

    config = fopen("config.txt","r");

    if(method == VARIATIONAL) {

        //list<Walker>::iterator it; //define an iterator over a list of walkers
        thermalized = false; //not thermalized

        //now we are filling the list pushing back
        for( int i = 0; i < walkers_num; i++ ) { //runs over the desired number of starting walkers
            Walker new_walker(particles_num); //define a new class walker

        for( int j = 0; j < particles_num; j++ ) {
            fscanf(config,"%lf\t%lf\t%lf\n",&new_walker.x[j], &new_walker.y[j], &new_walker.z[j]);
        }

            population.push_back(new_walker); //push the walker at the new starting position of the filled list
        }

        //for(int k = 0; k <= t_time; k++) diffuseNoBranch(); //diffuses for thermalization
    }

    fclose(config);

    n = 0.; //set parameters to be the good ones
    d = 0.;
    nErr = 0.;

    thermalized = true; //now the system is thermalized and the system is ready for convergence
}
//------------------------------------------------------------------------------

void MonteCarlo::diffuseNoBranch() { //variational MC

    double en_new, pf_1, pf_2, pf_3;
    list<Walker>::iterator it; //define an iterator over a list of walkers

    for( it = population.begin(); it != population.end(); it++ ) {
    //sum to each walker a random number gaussian distributed with the pseudo force component
        for( int j = 0; j < particles_num; j++ ) {

            pf_1 = pseudo_force(it->x, it->y, it->z, 0, j);
            pf_2 = pseudo_force(it->x, it->y, it->z, 1, j);
            pf_3 = pseudo_force(it->x, it->y, it->z, 2, j);

            it->x[j] += ( gauss_random_generator(dtau) + pf_1 );
            it->y[j] += ( gauss_random_generator(dtau) + pf_2 );
            it->z[j] += ( gauss_random_generator(dtau) + pf_3 );
        }

        //re-compute the local energy

        if( thermalized ) { //if thermalization is achived, compute stuff

            en_new = local_energy(it->x, it->y, it->z);
            it->weight = exp( - ( (en_new + en_old)/2. - Et ) * dtau );
            en_old = en_new;

            //compute energy terms, numerator, denominator and error numerator
            n += it->weight * en_new;
            d += it->weight;
            nErr += it->weight * en_new * en_new;
        }
    }
}
//------------------------------------------------------------------------------

void MonteCarlo::diffuseBranch() { //diffusion MC

    double en_old, en_new, pf_1, pf_2, pf_3;
    int M = 0, k = 0, dead_counter = 0, spawning_counter = 0;
    list<Walker>::iterator it; //define an iterator over a list of walkers
    list<Walker>::iterator end = population.end();

    for ( it = population.begin(); it != end; ) { //runs over all the walkers

        //compute the local energy for the present configuration
        en_old = local_energy(it->x, it->y, it->z);

        //sum to each walker a random number gaussian distributed with the pseudo force component
        for( int j = 0; j < particles_num; j++ ) {

            pf_1 = pseudo_force(it->x, it->y, it->z, 0, j);
            pf_2 = pseudo_force(it->x, it->y, it->z, 1, j);
            pf_3 = pseudo_force(it->x, it->y, it->z, 2, j);

            it->x[j] += ( gauss_random_generator(dtau) + pf_1 );
            it->y[j] += ( gauss_random_generator(dtau) + pf_2 );
            it->z[j] += ( gauss_random_generator(dtau) + pf_3 );
        }

        //re-compute the local energy
        en_new = local_energy(it->x, it->y, it->z);

        //compute weight
        it->weight = exp( - ( (en_new + en_old)/2. - Et ) * dtau );
        it->molteplicity = (int)( it->weight * NFIX/walkers_num + rng() );

        if( it->molteplicity == 0 ) { //if molteplicity is 0 delete the walker
            it = population.erase(it); //using function erase, which returns the iterator pointing to this position
            continue;
        }

        for( int j = it->molteplicity; j > 1; j-- ) { //copy molteplicity-times the walker
            population.push_back(*it); //push the walker at the new starting position of the filled list
            //*it is the dereferenced walker, which IS the class in the list
        }

        //compute energy terms, numerator, denominator and error numerator
        n += it->weight * en_new;
        d += it->weight;
        nErr += it->weight * en_new * en_new;

        ++it;
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

    int k = 0;
    double meanE = 0., meanE2 = 0., rel_err;

    //start thermalization process
    //generate a new homogeneous population and set summation parameters to 0
    thermalize( VARIATIONAL, (int)(TMINV/DTAU) );

     do {
        diffuseNoBranch(); //diffusion process for VMC: no branching done
        if( k%100 == 0 ) {
            Et = n/d; //redefine trial energy every 200 steps as the ratio n/d
            meanE += Et;
            meanE2 += Et*Et;
            fflush(stdout);
            printf("\t#Cycles: %d/%d\r", k, (int)(TIMEV/DTAU));
        }
        k++; //increment the counter
     } while( k < (TIMEV/DTAU) ); //redo until the total diffusion time has been reached


    //print interesting stuff
    if( exit_print == 1 ) {
        //err = sqrt(squareE - Esquare)/blocks;
        meanE /= (TIMED/DTAU);
        meanE2 /= (TIMED/DTAU);
        err = sqrt(meanE2 - meanE);
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
    double meanE, meanE2, blocks = 0., err;
    FILE* pt;
    FILE* pt2;
    FILE *pt3;

    pt = fopen("walkers.txt", "w");
    pt2 = fopen("energy.txt","w");
    pt3 = fopen("energy_error.txt","w");

    thermalize( VARIATIONAL, (int)(TMINV/DTAU) ); //needed to reset n and d to 0
    mean = 0;

    do {
        fprintf(pt,"%d\t%d\n", k, walkers_num);
        diffuseBranch(); //diffuse

        if( k%1000 == 0 ) {
            Et = n/d; //redefine trial energy every 200 steps as the ratio n/d
            meanE += Et/((TIMED/DTAU));
            meanE2 += Et*Et/(TIMED/DTAU);
            fprintf(pt2,"%d\t%g\n",k,n/d);
            printf("\t#(-) Cycles: %d/%d\r", k, (int)(TIMED/DTAU));
            fflush(stdout);
        }

        if( k%1000 == 0 ) {
            err = sqrt(meanE2 - meanE);
            fprintf(pt3,"%d\t%g\t%g\n", k, n/d, err);
        }

        k++;

        mean_walkers();
    } while( k < (TIMED/DTAU) );

    fclose(pt);
    fclose(pt2);
    fclose(pt3);

    printf("\t#Energy: %.3f +- %.3f [nat u]\n", n/d, err);
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

//    printf("\n\t#START VARIATIONAL MONTE CARLO");
//    printf("\n\t#-----------------------------------------");
//    printf("\n\t#Alpha \t#Beta \t#Energy    #Error on E\n");


//    M->prms[1] = 0.;//0.3;
//    for( M->prms[0] = 0.6; M->prms[0] <= 0.6; M->prms[0] += STEPa ) {
//        //for( M->prms[1] = LAMBDA; M->prms[1] <= LAMBDA; M->prms[1] += STEPb ) {
//            //cycle over the values of the variational parameter

//            M->Et = 1.5;
//            M->diffuse_no_branch_cycle(1);

//            if( sol > M->get_n()/M->get_d() ) {
//                sol = M->get_n()/M->get_d(); //save the minimum value of the energy
//                prm1 = prms[0];
//                prm2 = prms[1];
//            }
//        //}
//    }

//    printf("\t#-----------------------------------------\n");
//    printf("\n\t#SOLUTION\n");
//    printf("\t#-----------------------------------------\n");
//    printf("\t#Alpha = %g +- %g\n", prm1, STEPa);
//    printf("\t#Beta = %g +- %g\n", prm2, STEPb);
//    printf("\t#Energy = %.3f +- %.3f\n",sol, M->err);

//    printf("\t#-----------------------------------------\n");

    printf("\n\t#GENERATING POPULATION FOR DMC...\n");
    printf("\r\t#-----------------------------------------\n");

    M->prms[0] = 0.6;
    M->prms[1] = 1.523;
    M->Et = sol;

   // M->diffuse_no_branch_cycle(0);

    printf("\n\t#START DIFFUSION MONTE CARLO\n");
    printf("\r\t#-----------------------------------------\n");

    M->diffuse_branch_cycle();

   return 0;

}
