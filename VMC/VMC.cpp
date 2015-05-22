#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <gsl/gsl_rng.h>
#include <time.h>
#include <string.h>
#include <armadillo>

using namespace arma;


//------------------------------------------MACROS-----------------------------------------
#define FN(F,x) (*((F)->function))(x,(F)->params) //macro per le funzioni
#define A (10) //dimensione della mesh
#define WALKERS (10) //numero di walkers
#define DTAU (1.e-4) //step temporale
#define TIME (20.) //tempo immaginario totale di propagazione
#define TMIN (5.) //tempo minimo dopo il quale campionare l'energia
#define A_MIN (0.1) //valore minimo per il parametro variazionale
#define A_MAX (1.6) //valore massimo per il parametro variazionale
#define A_STEP (0.1) //step per il parametro variazionale
#define EDIM (int)(A_MAX/A_STEP - A_MIN/A_STEP) //dimensione per l'array che conterrà le varie energie
//----------------------------------------END-MACROS---------------------------------------


//-----------------------------------------STRUCTS-----------------------------------------
typedef struct {
	double (*function)(double x, void *params);
	void *params;
} function;

typedef struct {
    double *population; //contiene la popolazione attuale
    double *weight; //array dei pesi da associare a ciascun punto
    double Et; //energia di trial, che andiamo a modificare dinamicamente
    double n; //numeratore: contiene l'accumulazione di wi*EL
    double d; //denominatore: contiene l'accumlazione di wi
    double alpha; //parametro variazionale
} diffusion;
//---------------------------------------END-STRUCTS----------------------------------------


//in input quello che serve per la gsl: genera un numero random distribuito tra 0 e 1
double uniform_random_generator(gsl_rng * r) {

    return gsl_rng_uniform(r);
}


/*
 * prende in pasto due numeri generati con distribuzione uniforme tra 0 e 1
 * e ritorna un numero con distribuzione gaussiana con varianza tau
 */
double gauss_random_generator(double x, double y, double tau) {
    
    return sqrt( - 2.*log(x) ) * cos(2.*M_PI*y) * sqrt(tau);
}


/*
 * generatore di popolazioni iniziali (starting population) uniformi di walkers tra -A/2 e A/2:
 * necessita della struct diffusion per salvare in population la popolazione
 */
void sp_generator(diffusion *D, gsl_rng *r) {
    
    for(int i = 0; i < WALKERS; i++ ) D->population[i] = ( .5 - uniform_random_generator(r) )*A;
}


/*
 * energia locale per il calcolo variazionale
 * può sempre essere calcolata analiticamente.
 * Nel caso in questione: psi_T = e^(-ax^2), potenziale = 1/2 x^2
 * E_L := H|psi_T> / |psi_T>
 */
double local_energy(double x, void *params) {
    diffusion *D = (diffusion *)params;
    double alpha = D->alpha;

    return alpha + x * x * (1./2. - 2. * alpha * alpha); //oscillatore armonico
}


/*
 * pseudoforza introdotta dall'importance sampling
 * PF := nabla|psi_T>/|psi_T>
 * nel caso in questione: psi_T = e^(-ax^2)
 * sempre calcolabile analiticamente
 */
double pseudo_force(double x, void *params) {
    diffusion *D = (diffusion *)params;
    double alpha = D->alpha;

    return - ( 2. * alpha * x ) * DTAU; //oscillatore armonico
}


/*
 * VUOLE IN PASTO:
 * energia ridotta, pseudoforza, struct diffusion, generatore di numeri random e condizione di termalizzazione
 * COSA FARE:
 * creare la nuova popolazione
 * propagare
 * calcolare i pesi
 * non branchare
 * le x sono contenute in population! NON vanno ricostruite come i*DELTA!
 * therm mi dice se cominciare a campionare per pescare l'energia (0 = non termalizzazione) o no (qualunque intero = termalizzazione)
 * IMPORTANCE SAMPLING:
 * alla propagazione degli walkers bisogna aggiungere un termine (Nabla psi_T/psi_T) * tau
 */
void diffuse(function *EL, function *PF, diffusion *D, gsl_rng *r, int therm) {

    double r1, r2, en_old, en_new;

    for( int i = 0; i < WALKERS; i++ ) {
        //calcola l'energia locale nella configurazione attuale
        en_old = FN(EL, D->population[i]);

            r1 = uniform_random_generator(r); //primo numero random per la distribuzione gaussiana
            r2 = uniform_random_generator(r); //secondo numero random per la distribuzione gaussiana

        //sommo a ciascun walker un numero random distribuito gaussiano ed il pezzo di importance sampling
        D->population[i] += ( gauss_random_generator(r1, r2, DTAU) + FN(PF, D->population[i]) );

        //ricalcola l'energia locale
        en_new = FN(EL, D->population[i]);

            if( therm == 0 ) { //se il ciclo NON è di termalizzazione calcola pesi e somma i termini per l'energia

                //calcola il peso
                D->weight[i] = exp( - ( (en_new + en_old)/2. - D->Et ) * DTAU );

                //calcola i termini per l'energia, numeratore (n) e denominatore (d)
                D->n += D->weight[i] * FN(EL, D->population[i]);
                D->d += D->weight[i];
            }
    }

}

//trova il minimo in un array
void find_min(double *energy, double *solution) {

    double sol = energy[0];
    int prm = 0;

    for(int j = 1; j < EDIM; j++) {
        if(sol > energy[j]) {
            sol = energy[j];
            prm = j;
        }
    }

    solution[0] = sol;
    solution[1] = prm;

}

//esegui tutti i passaggi per il calcolo variazionale
void VMC(function *EL, function *PF, diffusion *D, gsl_rng *r) {

    int k = 0, j = 0, prm, therm;
    double energy[EDIM], param[EDIM], solution[2];
    FILE *pt;

    pt = fopen("energy.txt", "w");

    printf("\t#Alpha \t Energy\n");

    for( D->alpha = A_MIN; D->alpha <= A_MAX; D->alpha += A_STEP ) {

        sp_generator(D,r); //generiamo una popolazione iniziale

        therm = 1; //cicli di termalizzazione
            do {
                diffuse(EL,PF,D,r,therm);
                k++;
            } while( k < (TMIN/DTAU) );

        //azzero le variabili che competono alla definizione dell'energia
            D->n = 0;
            D->d = 0;

        therm = 0; //cominciamo a campionare l'energia
            do {
                diffuse(EL,PF,D,r,therm); //diffondiamo
                k++;

                if( k%500 == 0 ) D->Et = D->n/D->d; //ridefiniamo l'energia di trial ogni 500 cicli

            } while( k < (TIME/DTAU) );

            k = 0;

        fflush(stdout);
        printf("\t%g \t %g\n", D->alpha, D->n/D->d);
        fprintf(pt, "%g \t %g\n", D->alpha, D->n/D->d);
        energy[j] = D->n/D->d;
        param[j] = D->alpha;
        j++;
    }

    fclose(pt);

    printf("\t------------------------\n");
    find_min(energy, solution);
    j = solution[1];
    printf("\n\t#SOLUTION: alfa = %g, energy = %g\n", param[j], solution[0]);

}


int main() {

    double population[WALKERS], weight[WALKERS];

    //--------------------------------------------------------------------------
    //Definiamo le struct di cui abbiamo bisogno: diffusion e funzioni
    diffusion D;
    function EL, PF;

    D.population = population;
    D.weight = weight;
    D.d = 0.;
    D.n = 0.;
    D.Et = 7.; //trial energy

    EL.function = local_energy;
    EL.params = &D;

    PF.function = pseudo_force;
    PF.params = &D;
    //--------------------------------------------------------------------------

    //--------------------------------------------------------------------------
    //Definiamo le quantità necessarie alla generazione di numeri random con GSL
    srand(time(NULL)); //seed per i numeri random
    const gsl_rng_type * T;
    gsl_rng * r;

    T = gsl_rng_default;
    r = gsl_rng_alloc (T);
    gsl_rng_set(r, rand()); //diverso seed
    //--------------------------------------------------------------------------

    printf("\t------------------------\n");
    printf("\t#VARIATIONAL MONTE CARLO\n");
    printf("\t------------------------\n");

        //svolgi il calcolo variazionale
        VMC(&EL, &PF, &D, r);

    gsl_rng_free (r);

	return 0;
}
