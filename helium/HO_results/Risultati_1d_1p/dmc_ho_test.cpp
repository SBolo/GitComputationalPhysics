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
#define WALKERS (1200) //number of walkers
#define NFIX (1000.) //mean number of walkers I want at convergence
#define PARTICLE_NUM (1) //number of particles
#define LAMBDA (0.25) //width of the harmic well
#define BETA (358.93) //constant in front of LJ potential
#define CONVERSION (0.297) //conversion factor for energy
#define NPARAMS (1) //number of variational parameters
#define DTAU (5.e-4) //imaginary time step
#define TIMEV (10.) //propagation time for VMC
#define TMINV (20.) //minimum time after the energy sampling begins
#define TIMED (25.) //propagation time for DMC
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
    delete [] x;
    delete [] y;
    delete [] z;
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
        delete [] population[i]->x;
        delete [] population[i]->y;
        delete [] population[i]->z;

        delete [] population[i];
    }
    delete [] population; //delete the array
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

double local_energy(double *Rx, double *Ry, double *Rz, double *var_param, int particle_num) {
    //local energy for n 1d uncoupled harmonic oscillators
    double a = var_param[0];
    //double b = var_param[1];
    //double b5 = b * b * b * b * b, l = LAMBDA * LAMBDA * LAMBDA * LAMBDA * LAMBDA * LAMBDA;
    double sum = 0., tot = 0., x;
    double r, r_ij, r_ij7, r_ij6;

    //psi_T = f*g dove f = psi_HO e g = psi_Jastrow = psi_J
    //nabla^2psi_T/psi_T = nabla^2f/f + nabla^2g/g + 2 nablaf/f*nablag/g

    for( int i = 0; i < particle_num; i++ ) {
        x = Rx[i];// + Ry[i]*Ry[i] + Rz[i]*Rz[i] );
        sum += a + x * x * (1./2. - 2. * a * a); //nabla^2f/f+V_HO
    }

    return sum;
}
//------------------------------------------------------------------------------

double pseudo_force(double *Rx, double *Ry, double *Rz, double *var_param, int particle, int direction, int particle_num) {
    //pseudo force for the harmonic oscillator
    double a = var_param[0];
    //double b = var_param[1];
    //double b5 = b * b * b * b * b, sum = 0., r_ij, r_ij7;

    if( direction == 0 ) {
        return ( - 2. * a * Rx[particle] ) * DTAU;

    } else if (direction == 1 ) {
        return ( - 2. * a * Ry[particle] ) * DTAU;

    } else if( direction == 2 ) {
        return ( - 2. * a * Rz[particle]) * DTAU;

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

double MonteCarlo::psi_trial(double *Rx, double *Ry, double *Rz) {
    double sum = 0., r2;
    double a = prms[0];

    for( int i = 0; i < particles_num; i++ ) {
        r2 = Rx[i] * Rx[i];// + Ry[i] * Ry[i] + Rz[i] * Rz[i];
        sum -= a * r2
    }

    return exp(sum);
}

void MonteCarlo::thermalize(int method, int t_time) { //thermalization process

    FILE *pt;

    if(method == VARIATIONAL) { //we start building the population for the variational MC
        thermalized = false; //not thermalized, no diffusion has already been done

        for( int i = 0; i < walkers_num; i++ ) { //for each walker
            for( int j = 0; j < particles_num; j++ ) { //for each particle in the system
                population[i]->x[j] = ( .5 - rng() ) * MDIM; //generates a uniform population x [-MDIM/2,MDIM/2]
                //population[i]->y[j] = ( .5 - rng() ) * MDIM; //generates a uniform population y [-MDIM/2,MDIM/2]
                //population[i]->z[j] = ( .5 - rng() ) * MDIM; //generates a uniform population z [-MDIM/2,MDIM/2]
            }
        }

        for(int k = 0; k <= t_time; k++) diffuseNoBranch(); //diffuses for thermalization, for time t_time
    }

    //set parameters to be the good ones, avoiding summing over garbage
    n = 0;
    d = 0;
    nErr = 0;

    thermalized = true; //now the system is thermalized and the system is ready for convergence
}
//------------------------------------------------------------------------------

void MonteCarlo::diffuseNoBranch() { //variational MC

    double en_old, en_new, pf, rn1, rn2, psit_new, psit_old, alpha, diff_pos_x, squared_pf_new, double_product;

    for( int i = 0; i < walkers_num; i++ ) { //for every walker

        double tmp_x[particles_num], tmp_y[particles_num], tmp_z[particles_num];

        rn1 = gauss_random_generator(dtau);

        //compute the local energy for the present configuration
        //printf("%g, %g, %g\n", population[i]->x[3], population[i]->y[3], population[i]->z[3]);
        en_old = local_energy(population[i]->x, population[i]->y, population[i]->z, prms, particles_num);
        //printf("en_old = %g, prms1 = %g, prms2 = %g, particles = %d\n", en_old, prms[0], prms[1], particles_num);

        psit_old = psi_trial(population[i]->x, population[i]->y, population[i]->z)
                * psi_trial(population[i]->x, population[i]->y, population[i]->z);



        for( int j = 0; j < particles_num; j++ ) {

            //sum to each walker a random number gaussian distributed with the pseudo force component
            //pf = pseudo_force(population[i]->x,population[i]->y,population[i]->z,prms,j,0,particles_num);
            //population[i]->x[j] += ( gauss_random_generator(dtau) + pf );
            tmp_x[j] = population[i]->x[j] + gauss_random_generator(dtau);
            //tmp_y[j] = population[i]->y[j] + gauss_random_generator(dtau);
            //tmp_z[j] = population[i]->y[j] + gauss_random_generator(dtau);

            diff_pos_x = (tmp_x[j] - population[i]->x[j]) * (tmp_x[j] - population[i]->x[j]);
//            pf = pseudo_force(population[i]->x,population[i]->y,population[i]->z,prms,j,1,particles_num);
//            population[i]->y[j] += ( gauss_random_generator(dtau) + pf );
//            pf = pseudo_force(population[i]->x,population[i]->y,population[i]->z,prms,j,2,particles_num);
//            population[i]->z[j] += ( gauss_random_generator(dtau) + pf );

            squared_pf_old = pseudo_force(population[i]->x,population[i]->y,population[i]->z,prms,j,0,particles_num) *
                    pseudo_force(population[i]->x,population[i]->y,population[i]->z,prms,j,0,particles_num) +
                    pseudo_force(population[i]->x,population[i]->y,population[i]->z,prms,j,1,particles_num) *
                    pseudo_force(population[i]->x,population[i]->y,population[i]->z,prms,j,1,particles_num) +
                    pseudo_force(population[i]->x,population[i]->y,population[i]->z,prms,j,2,particles_num) *
                    pseudo_force(population[i]->x,population[i]->y,population[i]->z,prms,j,2,particles_num);

            squared_pf_new = pseudo_force(tmp_x,tmp_y,tmp_z,prms,j,0,particles_num) *
                    pseudo_force(tmp_x,tmp_y,tmp_z,prms,j,0,particles_num) +
                    pseudo_force(tmp_x,tmp_y,tmp_z,prms,j,1,particles_num) *
                    pseudo_force(tmp_x,tmp_y,tmp_z,prms,j,1,particles_num) +
                    pseudo_force(tmp_x,tmp_y,tmp_z,prms,j,2,particles_num) *
                    pseudo_force(tmp_x,tmp_y,tmp_z,prms,j,2,particles_num);

            doble_product = diff_pos_x *

        }

        psit_new = psi_trial(tmp_x, tmp_y, tmp_z) * psi_trial(tmp_x, tmp_y, tmp_z);

        alpha = (psit_new/psit_old) * exp()




        //re-compute the local energy
        en_new = local_energy(population[i]->x, population[i]->y, population[i]->z, prms, particles_num);
        //fprintf(pt,"%d\t%g\n",i, en_new);

        //printf("en_new = %g\n",en_new);

        //if thermalization is achived, compute stuff
        if( thermalized ) {

            //compute weight
            population[i]->weight = exp( - ( (en_new + en_old)/2. - Et ) * dtau );
            //printf("Peso = %g\n",population[i]->weight);

            //compute energy terms, numerator, denominator and error numerator
            n += population[i]->weight * local_energy(population[i]->x, population[i]->y, population[i]->z, prms, particles_num);
            d += population[i]->weight;
            nErr += population[i]->weight * local_energy(population[i]->x, population[i]->y, population[i]->z, prms, particles_num)
                    * local_energy(population[i]->x, population[i]->y, population[i]->z, prms, particles_num);
        }
    }
}
//------------------------------------------------------------------------------

void MonteCarlo::diffuseBranch() {

    double en_old, en_new, w_m, pf, eff_weight;
    int M = 0, h = 0;

    for( int i = 0; i < walkers_num; i++ ) {

        //compute the local energy for the present configuration
        en_old = local_energy(population[i]->x, population[i]->y, population[i]->z, prms, particles_num);

        for( int j = 0; j < particles_num; j++ ) {
            //sum to each walker a random number gaussian distributed with the pseudo force component
            pf = pseudo_force(population[i]->x,population[i]->y,population[i]->z,prms,j,0,particles_num);
            population[i]->x[j] += ( gauss_random_generator(dtau) + pf );
//            pf = pseudo_force(population[i]->x,population[i]->y,population[i]->z,prms,j,1,particles_num);
//            population[i]->y[j] += ( gauss_random_generator(dtau) + pf );
//            pf = pseudo_force(population[i]->x,population[i]->y,population[i]->z,prms,j,2,particles_num);
//            population[i]->z[j] += ( gauss_random_generator(dtau) + pf );
        }

        //re-compute the local energy
        en_new = local_energy(population[i]->x, population[i]->y, population[i]->z, prms, particles_num);


        //compute weight
        population[i]->weight = exp( - ( (en_new + en_old)/2. - Et ) * dtau );;
        eff_weight = exp( - ( (en_new + en_old)/2. - (Et + 1./dtau * log(NFIX/walkers_num)) ) * dtau );


        //compute molteplicity
        population[i]->molteplicity = (int)( eff_weight + rng() );

        //compute the new dimension of walkers
        M += population[i]->molteplicity;

        //compute statistical stuff
        n += population[i]->weight * local_energy(population[i]->x, population[i]->y, population[i]->z, prms, particles_num);
        d += population[i]->weight;
        nErr += population[i]->weight * local_energy(population[i]->x, population[i]->y, population[i]->z, prms, particles_num)
                * local_energy(population[i]->x, population[i]->y, population[i]->z, prms, particles_num);
    }

    //Branching time
    Walker **tmp;
    tmp = new Walker *[M];
    for(int i = 0; i < M; i++) {
       tmp[i] = new Walker(particles_num);
    }

    for( int i = 0; i < walkers_num; i++ ) {
        if( population[i]->molteplicity > 0 ) {
            for( int j = 0; j < population[i]->molteplicity; j++ ) {
                for( int k = 0; k < particles_num; k++ ) {
                    tmp[h]->x[k] = population[i]->x[k];
                    //tmp[h]->y[k] = population[i]->y[k];
                    //tmp[h]->z[k] = population[i]->z[k];
                }
                h++;
            }
        } else if( population[i]->molteplicity == 0 ) {
            continue;
        }
    }

    delete [] population;

    walkers_num = M;

    population = new Walker *[M];
    for(int i = 0; i < M; i++) {
       population[i] = new Walker(particles_num);
    }

    for( int i = 0; i < walkers_num; i++ ) {
        for( int j = 0; j < particles_num; j++ ) {
            population[i]->x[j] = tmp[i]->x[j];
            //population[i]->y[j] = tmp[i]->y[j];
            //population[i]->z[j] = tmp[i]->z[j];
        }
        population[i]->weight = 0.;
        population[i]->molteplicity = 1;
    }

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
            fflush(stdout);
            printf("\t#Cycles: %d/%d\r", k, (int)(TIMEV/DTAU));
        }
        k++; //increment the counter
     } while( k < (TIMEV/DTAU) ); //redo until the total diffusion time has been reached


    //print interesting stuff
    if( file_print == 1 ) {
        err = sqrt(squareE - Esquare)/blocks;
        rel_err = err/(n/d) * 100;
        fflush(stdout);
        printf("\r\t#-----------------------------------------\n");
        printf("\t%g \t %.4f \t %g \t %.2f %%\n", prms[0], n/d, err, rel_err );
        fprintf(pt, "\t%g \t %.4f \t %g\n", prms[0], n/d, sqrt(squareE - Esquare)/blocks );
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
    double sol = 100., prm1, prm2, err;
    double *prms = new double[NPARAMS];
    FILE *pt;

    MonteCarlo *M;
    M = new MonteCarlo(&local_energy, &pseudo_force, WALKERS, MDIM, PARTICLE_NUM, DTAU);
    M->prms = prms;

    pt = fopen("variational_parameters.txt", "w");

    printf("\n\t#START VARIATIONAL MONTE CARLO");
    printf("\n\t#-----------------------------------------");
    printf("\n\t#Alpha \t#Energy \t#Error \t#Relative error\n");


    for( M->prms[0] = 0.1; M->prms[0] <= 1.; M->prms[0] += STEPa ) {
        //for( M->prms[1] = LAMBDA; M->prms[1] <= LAMBDA+0.2; M->prms[1] += STEPb ) {
            //cycle over the values of the variational parameter

            M->Et = 1.5;
            M->diffuse_no_branch_cycle(pt, 1);

            if(sol > M->get_n()/M->get_d()) {
                sol = M->get_n()/M->get_d(); //save the minimum value of the energy
                prm1 = prms[0];
                prm2 = prms[1];
                err = M->err;
            }
        //}
    }

    fclose(pt);

    printf("\t#-----------------------------------------\n");
    printf("\n\t#SOLUTION\n");
    printf("\t#-----------------------------------------\n");
    printf("\t#Alpha = %g +- %g\n", prm1, STEPa);
    //printf("\t#Beta = %g +- %g\n", prm2, STEPb);
    printf("\t#Energy = %g +- %g\n",sol, err);

/*
    printf("\t#-----------------------------------------\n");

    printf("\n\t#GENERATING POPULATION FOR DMC...\n");
    printf("\r\t#-----------------------------------------\n");

    M->prms[0] = prm1;
    //M->prms[1] = prm2;
    M->Et = sol;

    M->diffuse_no_branch_cycle(pt, 0);

    printf("\n\t#START DIFFUSION MONTE CARLO\n");
    printf("\r\t#-----------------------------------------\n");

    M->diffuse_branch_cycle();
    */

    return 0;
}
