#include <iostream>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <gsl/gsl_rng.h>
#include <time.h>

using namespace std;

//#define DEBUG //generates always to same sequence of random numbers in case debugging is needed

//------------------------------------------MACROS-----------------------------------------
#define MDIM (20) //mesh dimension
#define WALKERS (20) //number of walkers
#define NFIX (20.) //mean number of walkers I want at convergence
#define PARTICLE_NUM (2) //number of particles
#define LAMBDA (0.25) //width of the harmic well
#define BETA (358.93) //constant in front of LJ potential
#define CONVERSION (0.297) //conversion factor for energy
#define NPARAMS (2) //number of variational parameters
#define DTAU (1.e-5) //imaginary time step
#define TIMEV (50.) //propagation time for VMC
#define TMINV (20.) //minimum time after the energy sampling begins
#define TIMED (70.) //propagation time for DMC
#define DIFFUSION (1) //to avoid thermalization process
#define VARIATIONAL (0) //to thermalize
#define STEPa (0.05) //step on the first variational parameter
#define STEPb (0.05) //step on the second variational parameter
//----------------------------------------END-MACROS---------------------------------------


class RandomNumbers { //random number generator
    public:
        RandomNumbers(); //constructor
        ~RandomNumbers(); //destructor
        double operator() (); //operator overloading

    private:
        const gsl_rng_type *T; //holds static information about each type of generator
        gsl_rng *r; //describes an instance of a generator created from a given gsl_rng_type
};


class Walker { //single walker
    public:
        Walker(int); //constructor: needs the particles number
        ~Walker(); //destructor

        double *x; //array of dimension particle_num: x coordinate
        double *y; //array of dimension particle_num: y coordinate
        double *z; //array of dimension particle_num: z coordinate
        double weight; //statistical weight associated with the walker
        int molteplicity; //molteplicity associated with the walker
};


class MonteCarlo { //everything about Monte Carlo
    public:
        MonteCarlo(double (*)(double *, double *, double *, double *, int),
                   double (*)(double *, double *, double *, double *, int, int, int), int, int, int, double); //constructor
        ~MonteCarlo(); //destructor

        double Et; //trial energy
        double *prms; //array of variational parameters
        int mean; //mean number of walkers
        double err; //error on the variational calculation of the energy
        void diffuse_no_branch_cycle(FILE *, int); //cycle needed to obtain a variational estimation for energy
        void diffuse_branch_cycle(); //cycle needed to get the energy via DMC
        double get_d(); //gets the value of d
        double get_n(); //gets the value of n
        double get_nErr(); //gets the value of nErr
        int get_wn(); //returns the number of walkers
        void mean_walkers(); //gets the mean value of walkers
        static void move(double *, int, double, double);

    private:
        RandomNumbers rng; //we need the rand class to generate the random numbers
        Walker **population; //array of classes which contains the population of walkers for each coordinate
        //I need to define it as a double pointer because it's the fast way to have an array of classes
        bool thermalized; //tells me if the system is thermalized or not
        int walkers_num; //number of walkers
        int particles_num; //number of particles
        int mesh_dim; //dimension of the mesh
        double dtau; // imaginary time step
        double n; //numerator: contains the sum w_i*EL
        double d; //denominator: contains the sum w_i
        double nErr; //numerator for the calculation of the error
        double (*local_energy)(double *, double *, double *, double *, int); //pointer to the local energy function
        double (*pseudo_force)(double *, double *, double *, double *, int, int, int); //pointer to the pseudo force function
        double gauss_random_generator(double); //gaussian distributed random generator
        void thermalize(int, int); //termalization process: generates the population and thermalizes up to a time t_time
        void diffuseNoBranch(); //starts the diffusion process for the Variational Monte Carlo
        void diffuseBranch(); //starts the diffusione process for the Diffusion Monte Carlo
        double psi_trial(double *, double *, double *); //needed for VMC
};


RandomNumbers::RandomNumbers() {
#ifdef DEBUG //if the define DEBUG is present, compile this, so that the seed is always the same
    srand(0);
#else //otherwise, compile with the seed given by the present time
    srand(time(NULL)); //seed
#endif
    //set the various stuff for the random number generator
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);
    gsl_rng_set(r, rand());
}
//------------------------------------------------------------------------------

RandomNumbers::~RandomNumbers() { //free the array r
        gsl_rng_free(r);
}
//------------------------------------------------------------------------------

double RandomNumbers::operator() () { //calling rng() I get a random number
    return gsl_rng_uniform_pos(r);
}
//------------------------------------------------------------------------------

Walker::Walker(int pn) { //when called, Walker allocates three arrays of dimension particle_number
    x = new double[pn];
    y = new double[pn];
    z = new double[pn];
}
//------------------------------------------------------------------------------

Walker::~Walker() { //delete the allocated arrays
    delete[] x;
    delete[] y;
    delete[] z;
}
//------------------------------------------------------------------------------

MonteCarlo::MonteCarlo(double (*l)(double *, double *, double *, double *, int),
                       double (*p)(double *, double*, double *, double *, int, int, int),
                       int wdim, int mdim, int pdim, double step) : rng() {
    //give values and allocate things
    walkers_num = wdim;
    particles_num = pdim;
    //this is the way we allocate an array of classes
    population = new Walker *[walkers_num]; //start allocating the pointer
    for(int i = 0; i < walkers_num; i++) { //cycle over walkers_num
       population[i] = new Walker(particles_num); //allocate each class separately
    }
    mesh_dim = mdim;
    dtau = step;
    thermalized = false; //when allocated, walkers are not thermalized
    local_energy = l;
    pseudo_force = p;
}
//------------------------------------------------------------------------------

MonteCarlo::~MonteCarlo() {
    //deallocate each single class in the array
    for( int i = 0; i < walkers_num; i++ ) {
        delete population[i];
    }
    delete[] population; //delete the array
}
//------------------------------------------------------------------------------

double MonteCarlo::get_d() {
    return d; //d is private, get its value
}
//------------------------------------------------------------------------------

double MonteCarlo::get_n() {
    return n; //n is private, get its value
}
//------------------------------------------------------------------------------

double MonteCarlo::get_nErr() {
    return nErr; //nErr is private, get its value
}
//------------------------------------------------------------------------------

int MonteCarlo::get_wn() {
    return walkers_num; //walkers_num is private, get its value
}
//------------------------------------------------------------------------------

double MonteCarlo::gauss_random_generator(double tau) { //given two uniformly distributed rng(), generates a gaussian distributed rng()
    return sqrt( - 2.*log(rng()) ) * cos(2.*M_PI*rng()) * sqrt(tau);
}
//------------------------------------------------------------------------------

double local_energy(double *Rx, double *Ry, double *Rz, double *var_param, int particle_num) {
    //local energy for n 1d uncoupled harmonic oscillators
    double a = var_param[0];
    double b = var_param[1];
    double b5 = b * b * b * b * b, l = LAMBDA * LAMBDA * LAMBDA * LAMBDA * LAMBDA * LAMBDA;
    double sum = 0., tot = 0.;
    double r, r_ij, r_ij7, r_ij6;

    //psi_T = f*g dove f = psi_HO e g = psi_Jastrow = psi_J
    //nabla^2psi_T/psi_T = nabla^2f/f + nabla^2g/g + 2 nablaf/f*nablag/g

    for( int i = 0; i < particle_num; i++ ) {
        r = sqrt( Rx[i]*Rx[i] + Ry[i]*Ry[i] + Rz[i]*Rz[i] );
        sum += ( 3. * a + r * r * (1./2. - 2. * a * a) ); //nabla^2f/f+V_HO

        for( int j = 0; j < particle_num; j++ ) {
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

    for( int i = 0; i < particle_num; i++ ) { //V_LJ
        for( int j = 0; j < i; j++ ) {
            r_ij = sqrt( (Rx[i] - Rx[j]) * (Rx[i] - Rx[j]) + (Ry[i] - Ry[j]) * (Ry[i] - Ry[j]) + (Rz[i] - Rz[j]) * (Rz[i] - Rz[j]) );
            r_ij6 = r_ij * r_ij * r_ij * r_ij * r_ij * r_ij;
            sum += BETA * ( l * l / (r_ij6 * r_ij6) - l / r_ij6 );
        }
    }

    return sum;
}
//------------------------------------------------------------------------------

double pseudo_force(double *Rx, double *Ry, double *Rz, double *var_param, int particle, int direction, int particle_num) {
    //pseudo force for the harmonic oscillator
    double a = var_param[0];
    double b = var_param[1];
    double b5 = b * b * b * b * b, sum = 0., r_ij, r_ij7;

    if( direction == 0 ) {
        for( int i = 0; i < particle_num; i++ ) {
            if( i != particle ) {
                r_ij = sqrt( (Rx[particle] - Rx[i]) * (Rx[particle] - Rx[i]) + (Ry[particle] - Ry[i]) *
                        (Ry[particle] - Ry[i]) + (Rz[particle] - Rz[i]) * (Rz[particle] - Rz[i]) );
                r_ij7 = r_ij * r_ij * r_ij * r_ij * r_ij * r_ij * r_ij;
                sum -= (5./2.) * b5 * (Rx[i] - Rx[particle]) / r_ij7;
            }
        }
        return ( - 2. * a * Rx[particle] + sum ) * DTAU;

    } else if (direction == 1 ) {
        for( int i = 0; i < particle_num; i++ ) {
            if( i != particle ) {
                r_ij = sqrt( (Rx[particle] - Rx[i]) * (Rx[particle] - Rx[i]) + (Ry[particle] - Ry[i]) *
                        (Ry[particle] - Ry[i]) + (Rz[particle] - Rz[i]) * (Rz[particle] - Rz[i]) );
                r_ij7 = r_ij * r_ij * r_ij * r_ij * r_ij * r_ij * r_ij;
                sum -= (5./2.) * b5 * (Ry[i] - Ry[particle]) / r_ij7;
            }
        }
        return ( -  2. * a * Ry[particle] + sum ) * DTAU;

    } else if( direction == 2 ) {
        for( int i = 0; i < particle_num; i++ ) {
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

void MonteCarlo::thermalize(int method, int t_time) { //thermalization process

    FILE *pt;

    if(method == VARIATIONAL) { //we start building the population for the variational MC
        thermalized = false; //not thermalized, no diffusion has already been done

        for( int i = 0; i < walkers_num; i++ ) { //for each walker
            for( int j = 0; j < particles_num; j++ ) { //for each particle in the system
                population[i]->x[j] = ( .5 - rng() ) * MDIM; //generates a uniform population x [-MDIM/2,MDIM/2]
                population[i]->y[j] = ( .5 - rng() ) * MDIM; //generates a uniform population y [-MDIM/2,MDIM/2]
                population[i]->z[j] = ( .5 - rng() ) * MDIM; //generates a uniform population z [-MDIM/2,MDIM/2]
            }
        }

        for(int k = 0; k <= t_time; k++) diffuseNoBranch(); //diffuses for thermalization, for time t_time
    }

    //set parameters to be the good ones, avoiding summing over garbage
    n = 0.;
    d = 0.;
    nErr = 0.;

    thermalized = true; //now the system is thermalized and the system is ready for convergence
}
//------------------------------------------------------------------------------

void MonteCarlo::diffuseNoBranch() { //variational MC

    double en_old, en_new, pf_1, pf_2, pf_3, rn;

    for( int i = 0; i < walkers_num; i++ ) { //for every walker

        //compute the local energy for the present configuration
        en_old = local_energy(population[i]->x, population[i]->y, population[i]->z, prms, particles_num);

        for( int j = 0; j < particles_num; j++ ) {
            //sum to each walker a random number gaussian distributed with the pseudo force component
            pf_1 = pseudo_force(population[i]->x,population[i]->y,population[i]->z,prms,j,0,particles_num);
            pf_2 = pseudo_force(population[i]->x,population[i]->y,population[i]->z,prms,j,1,particles_num);
            pf_3 = pseudo_force(population[i]->x,population[i]->y,population[i]->z,prms,j,2,particles_num);

            population[i]->x[j] += ( gauss_random_generator(dtau) + pf_1 );
            population[i]->y[j] += ( gauss_random_generator(dtau) + pf_2 );
            population[i]->z[j] += ( gauss_random_generator(dtau) + pf_3 );
        }

        //re-compute the local energy
        en_new = local_energy(population[i]->x, population[i]->y, population[i]->z, prms, particles_num);

        //if thermalization is achived, compute stuff
        if( thermalized ) {

            //compute weight
            population[i]->weight = exp( - ( (en_new + en_old)/2. - Et ) * dtau );

            //compute energy terms, numerator, denominator and error numerator
            n += population[i]->weight * en_new;
            d += population[i]->weight;
            nErr += population[i]->weight * en_new * en_new;
        }
    }
}
//------------------------------------------------------------------------------

void MonteCarlo::diffuseBranch() {

    double en_old, en_new, w_m, pf_1, pf_2, pf_3, eff_weight, rn;
    int M = 0, h = 0;
    Walker **tmp;

    for( int i = 0; i < walkers_num; i++ ) {

        //compute the local energy for the present configuration
        en_old = local_energy(population[i]->x, population[i]->y, population[i]->z, prms, particles_num);

        for( int j = 0; j < particles_num; j++ ) {
            //sum to each walker a random number gaussian distributed with the pseudo force component
            pf_1 = pseudo_force(population[i]->x,population[i]->y,population[i]->z,prms,j,0,particles_num);
            pf_2 = pseudo_force(population[i]->x,population[i]->y,population[i]->z,prms,j,1,particles_num);
            pf_3 = pseudo_force(population[i]->x,population[i]->y,population[i]->z,prms,j,2,particles_num);

            population[i]->x[j] += ( gauss_random_generator(dtau) + pf_1 );
            population[i]->y[j] += ( gauss_random_generator(dtau) + pf_2 );
            population[i]->z[j] += ( gauss_random_generator(dtau) + pf_3 );
        }

        //re-compute the local energy
        en_new = local_energy(population[i]->x, population[i]->y, population[i]->z, prms, particles_num);

        //compute weight
        population[i]->weight = exp( - ( (en_new + en_old)/2. - Et ) * dtau );
        eff_weight = exp( - ( (en_new + en_old)/2. - (Et + 1./dtau * log(NFIX/walkers_num)) ) * dtau );

        //compute molteplicity
        population[i]->molteplicity = (int)( eff_weight + rng() );

        //compute the new dimension of walkers
        M += population[i]->molteplicity;

        //compute statistical stuff
        n += population[i]->weight * en_new;
        d += population[i]->weight;
        nErr += population[i]->weight * en_new * en_new;
    }

    //Branching time
    tmp = new Walker *[M];
    for(int i = 0; i < M; i++) {
       tmp[i] = new Walker(particles_num);
    }

    for( int i = 0; i < walkers_num; i++ ) {
        if( population[i]->molteplicity > 0 ) {
            for( int j = 0; j < population[i]->molteplicity; j++ ) {
                for( int k = 0; k < particles_num; k++ ) {
                    tmp[h]->x[k] = population[i]->x[k];
                    tmp[h]->y[k] = population[i]->y[k];
                    tmp[h]->z[k] = population[i]->z[k];
                }
                h++;
            }
        } else if( population[i]->molteplicity == 0 ) {
            continue;
        }

        delete population[i];
    }

    delete[] population;

    walkers_num = M;

    population = new Walker *[M];
    for(int i = 0; i < M; i++) {
       population[i] = new Walker(particles_num);
    }

    for( int i = 0; i < walkers_num; i++ ) {
        for( int j = 0; j < particles_num; j++ ) {
            population[i]->x[j] = tmp[i]->x[j];
            population[i]->y[j] = tmp[i]->y[j];
            population[i]->z[j] = tmp[i]->z[j];
        }
        population[i]->weight = 0.;
        population[i]->molteplicity = 1;

        delete tmp[i];
    }

    delete[] tmp;

}
//------------------------------------------------------------------------------

void MonteCarlo::mean_walkers() {
    mean += walkers_num;
}
//------------------------------------------------------------------------------

void MonteCarlo::diffuse_no_branch_cycle(FILE *pt, int file_print) { //cycle needed to obtain an energy from a fixen variational parameter

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
        }

        if( k%5000 == 0 ) {
            fflush(stdout);
            printf("\t#Cycles: %d/%d\r", k, (int)(TIMEV/DTAU));
        }
        k++; //increment the counter
     } while( k < (TIMEV/DTAU) ); //redo until the total diffusion time has been reached


    //print interesting stuff
    if( file_print == 1 ) {
        err = sqrt(squareE - Esquare)/blocks;
        rel_err = err/(n/d) * 100;
        printf("\r\t#-----------------------------------------\n");
        printf("\t%g \t %g \t %.4f \t %g \t %.2f %%\n", prms[0], prms[1], n/d, err, rel_err );
        fprintf(pt, "\t%g \t %g \t %.4f \t %g\n", prms[0], prms[1], n/d, sqrt(squareE - Esquare)/blocks );
    }
}
//------------------------------------------------------------------------------

void MonteCarlo::diffuse_branch_cycle() {

    int k = 0;
    double squareE = 0., Esquare = 0., blocks = 0.;
    FILE* pt;
    FILE *pt2;

    pt = fopen("walkers.txt", "w");
    pt2 = fopen("energy_error.txt","w");

    thermalize(DIFFUSION, 0); //needed to reset n and d to 0
    mean = 0;

    do {
        diffuseBranch(); //diffuse

        if( k%200 == 0 ) {
            fprintf(pt,"%d\t%d\n", k, walkers_num);
        }

        if( k%1000 == 0 ) {
            Et = n/d; //redefine trial energy every 100 steps
            squareE += nErr/d;
            Esquare += Et*Et;
            blocks++;
        }

        if( k%5000 == 0 ) {
            fflush(stdout);
            fprintf(pt2,"%d\t%g\t%g\n", k, n/d, sqrt(squareE - Esquare)/blocks);
            printf("\t#(-) Cycles: %d/%d\r", k, (int)(TIMED/DTAU));
        }

        k++;
        mean_walkers();

    } while( k < (TIMED/DTAU) );

    fclose(pt);
    fclose(pt2);

    printf("\t#Energy: %.3f +- %.3f [nat u]\n", n/d, sqrt(squareE - Esquare)/blocks);
    printf("\t#Mean number of walkers = %d\n", (int)(1.*(mean)/k));
    printf("\t#----------------------------------\n");
}
//------------------------------------------------------------------------------


int main() {

    int k = 0, j = 0, therm, test;
    double sol = 100., prm1, prm2, err;
    double *prms = new double[NPARAMS];
    FILE *pt;

    MonteCarlo *M;
    M = new MonteCarlo(&local_energy, &pseudo_force, WALKERS, MDIM, PARTICLE_NUM, DTAU);
    M->prms = prms;

    pt = fopen("variational_parameters.txt", "w");

    printf("\n\t#START VARIATIONAL MONTE CARLO");
    printf("\n\t#-----------------------------------------");
    printf("\n\t#Alpha \t#Beta \t#Energy \t#Error \t#Relative error\n");


    for( M->prms[0] = 0.1; M->prms[0] <= 1; M->prms[0] += STEPa ) {
        for( M->prms[1] = LAMBDA; M->prms[1] <= LAMBDA+0.3; M->prms[1] += STEPb ) {
            //cycle over the values of the variational parameter

            M->Et = 3.;
            M->diffuse_no_branch_cycle(pt, 1);

            if(sol > (M->get_n()/M->get_d() + M->err) ) {
                sol = M->get_n()/M->get_d(); //save the minimum value of the energy
                prm1 = prms[0];
                prm2 = prms[1];
                err = M->err;
            }
        }
    }

    fclose(pt);

    printf("\t#-----------------------------------------\n");
    printf("\n\t#SOLUTION\n");
    printf("\t#-----------------------------------------\n");
    printf("\t#Alpha = %g +- %g\n", prm1, STEPa);
    printf("\t#Beta = %g +- %g\n", prm2, STEPb);
    printf("\t#Energy = %g +- %g\n",sol, err);


    printf("\t#-----------------------------------------\n");

    printf("\n\t#GENERATING POPULATION FOR DMC...\n");
    printf("\r\t#-----------------------------------------\n");

    M->prms[0] = prm1;
    M->prms[1] = prm2;
    M->Et = sol;

    //M->diffuse_no_branch_cycle(pt, 0);

    printf("\n\t#START DIFFUSION MONTE CARLO\n");
    printf("\r\t#-----------------------------------------\n");

    M->diffuse_branch_cycle();

    return 0;
}
