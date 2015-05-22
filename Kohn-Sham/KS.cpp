#include <stdio.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_eigen.h>
#include <math.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_integration.h>
#include <armadillo>

using namespace arma;

/*
 N.B. Prima dell'esecuzione del programma:
 -> settare numero di particelle (P_NUM)
 -> settare numero di shell energetiche (SHELL_NUM)
 -> eseguire
 */

//_____________________________________________________________GSL_MATRIX_MACROS________________________________________________

#define SOLVING(e,l,eval,deg) solve(energies,tmp,e,l,eval,deg, D)
#define STATE() system_energy(D)
#define ESUM() eigen_sum(D)

#define GSL(F,x) (*((F)->function))(x,(F)->params)
#define STEP (1.5e-2) //step di integrazione
#define RMAX (20.) //dimensione mesh
#define P_NUM (8.) //numero di elettroni nella sfera
#define SHELL_NUM (2) //numero di shell
#define RS_Li (4.) //raggio di Wiegner Seitz per cominciare
#define EPSILON (1.e-3) //precisione sulla somma degli autovalori
#define EDIM (6) //dimensione array energie

//________________________________________________________________END_MACROS____________________________________________________


//________________________________________________________________STRUCTURES____________________________________________________

typedef struct {
    /*Per la diagonalizzazione */
    mat matrix; //matrice quadrata: da diagonalizzare
    vec eigVal; //vettore GSL: conterrà gli autovalori
    mat eigVec; //matrice quadrata: contiene gli autovettori
    /* Parametri e altre cose */
    int dim; //dimensione della matrice
    vec psi; //funzione d'onda
    vec rho; //densità del sistema
    vec tmp; //per fare i conti se non voglio rovinare un array
    vec rho_tmp;
    vec angMom;
    vec mainQN;
    vec degDegree;
    double sphere_radius; //raggio della sfera carica
    double l; //momento angolare
} schrodinger;


//______________________________________________________________END_STRUCTURES____________________________________________________

//Simpson routine to normalize functions: old one, it works.
double simpson_normalize(vec psi, int N) {
    double h, k, res = 0;
    int i = 0;
    vec f(N);
	
        for(i = 0; i < N; i++) {
            f(i) = 4. * M_PI * psi(i) * pow(STEP*i, 2.); //solo per normalizzazione
        }
	
    k = f(0) + 4 * f(1) + f(N-1) + 4 * f(N-2);
    
        for(i = 2; i < N-2; i+=2) {
            res += 4 * f(i) + 2 * f(i+1);
        }
	
    return (k + res) * (STEP/3);
	
}

double integrate(vec g, double j, double h) { //array, primo estremo (primo indice), secondo estremo (secondo indice)
//integrale con trapezi: rapido e indolore
    int i;
    double risultato = 0.;

    for( i = j; i < h-1; i++) {
        risultato += (g(i) + g(i+1)) * STEP/2.;
    }

    return risultato;
}


double guess(double x, void *params) { //distribuzione di densità iniziale
    
    schrodinger * D = (schrodinger *)params;
    double radius = D->sphere_radius;


    double rho_plus = 3. / (4. * M_PI * pow(radius,3.));
    double mu = 1.; //deve essere piccolo: questo valore va approssimativamente bene
    
    double rho = rho_plus / (1. + exp(mu * (x - radius) ) );
    
    return 1./(1. + exp(mu * (x - radius))); //va integrata per essere normalizzata
    
}


double correlation_exchange(double r , schrodinger *D) {
    
	int i = (int) (r/STEP);
	double rho_r = D->rho[i];
    
	double r_s = pow((3.0 / (4.0 * M_PI * rho_r)), 1.0 / 3.0);
	double term = pow((3.0 * rho_r / M_PI), 1.0 / 3.0);
    
    return  - 0.44 / (7.8 + r_s) /*- 0.44 / 3.0 * r_s / ((7.8 + r_s) * (7.8 + r_s))*/ - term;
}


double coulomb(double x, schrodinger *D) { //porzione coulombiana del potenziale
    
    double radius = D->sphere_radius;
    double rho_plus = 3. / (4. * M_PI * pow(RS_Li,3.));
    double ext;
    
        if( x <= radius ) {
            ext = 2. * M_PI * rho_plus * ( (1./3.)*x*x - radius*radius );
        } else {
            ext = - 4./3. * M_PI * rho_plus * pow(radius,3.)/x;
        }

    return ext;
    
}


double hartree(double x, schrodinger *D) {
    
    int i;
    double ext, pot1, pot2;
    
        for(i = 0; i < D->dim; i++) {
            D->tmp(i) = D->rho(i) * (i * STEP);
        }
    
        pot1 = 4. * M_PI * integrate(D->tmp, (int)(x/STEP), (int)(RMAX/STEP));

        for(i = 0; i < D->dim; i++) {
            D->tmp(i) *= (i * STEP);
        }
    
        pot2 = 4. * M_PI * (1./(x+STEP)) * integrate(D->tmp, 2, (int)(x/STEP));

    return (pot1 + pot2);
}


double matrix_diagonal(double r, schrodinger *D) {
    //funzione per riempire la diagonale principale: potentiale + 1/(delta x)^2
        
        double l = (double)(D->l);
        double V;
        double CE = correlation_exchange(r, D);
        double H = hartree(r, D);
        double C = coulomb(r, D);

        V = CE + C + H + l * ( l + 1. ) / (2. * (r + STEP) * (r + STEP) );
    
    return ( V + 1./pow(STEP,2.) );
    
}


double diagonalization(int eval, schrodinger *D) { //scambio correlazione + coulomb + hartree
 
     //in input: gsl_function per il calcolo della diagonale, schrodinger per il settaggio della matrice e numero dell'autovalore desiderato.
     //L'autovettore è contenuto nella matrice all'interno della struttura. Lo prendo quando desidero, non serve farlo ora 
    
    D->matrix.fill(0.);
 
    for(int i = 0; i < (D->dim); i++) {
        for(int j = 0; j < (D->dim); j++) {
            
            if(i == j) D->matrix(i,j) = matrix_diagonal(i*STEP, D);//diagonale principale
            if(i == (j+1) || i == (j-1)) D->matrix(i, j) = - 1. / ( 2. * pow(STEP,2.) ); //diagonale superiore e inferiore
            
        }
    }
    
    eig_sym(D->eigVal, D->eigVec, D->matrix);
    
    //non normalizziamo, ma lo facciamo direttamente nella funzione self-consistence
    
    return D->eigVal(eval);
    
}

double system_energy(schrodinger *D) {
    
    double u_int, eps, eps_xc, rho, r_s;
    vec mean_uint(D->dim), rho2_de_xc(D->dim);
    
    for( int i = 0; i< D->dim; i++ ) {
        mean_uint(i) = hartree(i*STEP, D) * D->rho(i);
    }
    
    u_int = 1./2. * simpson_normalize(mean_uint, D->dim);
    
    for( int i = 0; i < D->dim; i++) {
        rho = D->rho(i);
        r_s = pow((3. / (4. * M_PI * rho)), 1. / 3.);
        rho2_de_xc[i] = -1./4. * pow( (3. / M_PI), 1./3. ) * pow( rho, 4./3. ) - 1./3. * ( 0.44 / pow((7.8 + r_s),2.) ) * pow( 3. / (4.*M_PI), 1./3.) * pow(rho, 2./3.);
    }
    
    eps_xc = simpson_normalize(rho2_de_xc, D->dim);
    
    
    return ( - u_int - eps_xc );
    
}

void solve(vec &energies, vec &tmp, int e, int l, int eval, int deg, schrodinger *D) {
//struttura, array per le energie, tmp per i valori temporanei, posto nell'array, momento angolare, autovalore, degenerazione
    
    double normalization;
    vec tmp_derivate(D->dim);
    int i;
    
    D->l = (1.*l);
     energies(e) = diagonalization(eval, D);
    
    for( i = 0; i < D->dim; i++ ) {
        D->psi(i) = pow( D->eigVec(i,eval) / ( (i+1)*STEP ), 2.);
    }
    
    normalization = simpson_normalize(D->psi, D->dim);
    
    for( i = 0; i < D->dim; i++ ) {
        D->psi(i) /= normalization;
        D->rho_tmp(i) += (1.*deg) * D->psi(i);
    }
}

void self_consistence(schrodinger *D, vec &energies, vec &tmp2, vec &tmp, double alpha, double &e_new) {

    for(int i = 0; i < D->dim; i++) {
        D->rho_tmp(i) = 0.;
        tmp2(i) = D->rho(i);
    }

    //solve(energies,tmp,e,l,eval,deg, D)

        for(int i = 0; i < SHELL_NUM; i++) {
            // arr l  n  deg
            SOLVING(i, D->angMom(i), D->mainQN(i), D->degDegree(i));
         }

    //salviamo il temporaneo dentro rho e ricalcoliamo l'energia dello stato
    for(int i = 0; i< D->dim; i++ ) {
        D->rho(i) = D->rho_tmp(i) * alpha + tmp2(i) * (1. - alpha);
    }

        e_new = 2.*energies(0) + 6.*energies(1) + 2.*energies(2) + 10.*energies(3) + 14.*energies(4) + 6.*energies(5); //40e
        //e_new = energies(0);

}


void SC_lopp(schrodinger *D, int e_dim) {
    
    double alpha, e_new, e_old, normalization, state;
    int i, k = 0, iteration = 0;
    vec tmp(D->dim), tmp2(D->dim), energies(e_dim);
    FILE *pt;
    char filen[20];
    int pn = (int)(P_NUM);
    
    sprintf(filen, "summary%de",pn);
    
    printf("\n");
    printf("----------\n");
    printf("SIMULATION\n");
    printf("----------\n");
    
    alpha = 0.05;
    
        for( i = 0; i < EDIM; i++ ) {
            energies(i) = 0.;
        }
    
        for(i = 0; i < D->dim; i++) {
            D->rho_tmp(i) = 0.;
            tmp2(i) = D->rho(i);
        }
    
    //per prima cosa risolvo l'equazione per il guess.
    self_consistence(D, energies, tmp2, tmp, alpha, e_new);
    
    //fino qui abbiamo calcolato la soluzione preliminare


            do {
                
                e_old = e_new; //salviamo in e_old la vecchia e_new
                
                self_consistence(D, energies, tmp2, tmp, alpha, e_new);
                
                iteration++;

                    printf("\n");
                    printf("\t# ITERATION = %d\r",iteration);
                    printf("\n");
                    printf("\t# e_old = %.3f | e_new = %.3f | (e_old - e_new) = %.4f\r",e_old, e_new,fabs(e_old - e_new));
                    printf("\n");
                    printf("\t# Single particle energies: -) 1s = %.3f\r",energies(0));
                    printf("\n");
                    printf("\t#                           -) 1p = %.3f\r", energies(1));
                    printf("\n");
                    printf("\t#                           -) 1d = %.3f\r", energies(2));
                    printf("\n");
                    printf("\t#                           -) 2s = %.3f\r", energies(3));
                    printf("\n");
                    printf("\t#                           -) 1f = %.3f\r", energies(4));
                    printf("\n");
                    printf("\t#                           -) 2p = %.3f\r", energies(5));
                    printf("\n");
                    printf("\t# State energy: E = %.3f\r",(e_new + STATE())/P_NUM);
                    printf("\n");
                    fflush(stdout);
            
            } while ( fabs( e_new - e_old ) > EPSILON);

    
    pt = fopen(filen,"w");
    //printf("\n");
    fprintf(pt,"-------------------------\n");
    fprintf(pt,"SUMMARY OF THE SIMULATION\n");
    fprintf(pt,"-------------------------\n");
            fprintf(pt,"#Number of particles: %g\n", P_NUM);
            fprintf(pt,"#Dimension of the matrix: %d\n", (int)(RMAX/STEP));
            fprintf(pt,"#Precision on the eigenvalues: %g\n", EPSILON);
            fprintf(pt,"#Correct sum of eigenvalues: %.4f\n", e_new);
            fprintf(pt,"#Correct single-level energies: 1s = %.3f H, %.3f eV\n                                1p = %.3f H, %.3f eV\n                                1d = %.3f H, %.3f eV\n                                2s = %.3f H, %.3f eV\n                                1f = %.3f H, %.3f eV\n                                2p = %.3f H, %.3f eV\n", energies[0], energies[0] * 27.21, energies[1], energies[1] * 27.21, energies[2], energies[2]*27.21, energies[3], energies[3]*27.21, energies[4], energies[4] * 27.21, energies[5], energies[5] * 27.21);
            fprintf(pt,"#Correct energy per particle: %.3f H\n", (e_new + STATE())/P_NUM);
            fprintf(pt,"#Correct energy per particle: %.2f eV\n", (e_new + STATE())/P_NUM * 27.21);
    //printf("\n");
    fclose(pt);
}

void printRho(schrodinger *D, char *filen) {

    FILE *pt;

    pt = fopen(filen,"w");

        for(int i = 0; i < D->dim; i++) {
            fprintf(pt, "%g \t %g \t \t \n", i*STEP, D->rho(i));
        }

    fclose(pt);

}

void printPot(schrodinger *D, char *filen) {

    FILE *pt;

    pt = fopen(filen,"w");

        for( int i = 0; i < D->dim; i++ ) {
           fprintf( pt, "%g \t %g \t %g \t %g\n", i*STEP, coulomb(i*STEP,D), hartree(i*STEP, D), correlation_exchange(i*STEP, D) );
           //potenziale autoconsistente
        }

    fclose(pt);

}


int main() {
    int i, j, eigen, sort, k = 0, pn;
    double valore, normalization, e, energy;
    char filen[20];
    
    pn = (int)(P_NUM);

    printf("\n");
    printf("--------------\n");
    printf("INITIAL CHECKS\n");
    printf("---------------\n");


    //preparo la struttura shrodinger
    schrodinger D;
    D.dim = RMAX/STEP;
    
    gsl_function RHO_GUESS;

    RHO_GUESS.function = guess;
    RHO_GUESS.params = &D;

    vec psi(D.dim), rho(D.dim), tmp(D.dim), rho_tmp(D.dim), rho_guess(D.dim), eigVal(D.dim);
    mat eigVec(D.dim, D.dim, fill::zeros), matrix(D.dim, D.dim, fill::zeros);

    vec angMom(6);
        angMom << 0 << 1 << 2 << 0 << 3 << 1;
    vec mainQN(6);
        mainQN << 0 << 0 << 0 << 1 << 0 << 1;
    vec degDegree(6);
        degDegree << 2 << 6 << 10 << 2 << 14 << 6;
    
    D.eigVal = eigVal;
    D.matrix = matrix;
    D.eigVec = eigVec;
    D.rho = rho;
    D.psi = psi;
    D.angMom = angMom;
    D.mainQN = mainQN;
    D.degDegree = degDegree;
    D.rho_tmp = rho_tmp;
    D.tmp = tmp;

    D.sphere_radius = RS_Li * pow(P_NUM,(1./3.));
    D.l = 0.;

    printf("\t-> Every vector correctly allocated\n");

    //normalizzazione del guess iniziale
        for(i=0; i< D.dim; i++) {
            D.rho(i) = GSL(&RHO_GUESS,i*STEP);
        }

            normalization = simpson_normalize(D.rho, D.dim);


        for (i=0; i < D.dim; i++) {
            D.rho(i) /= normalization;
            D.rho(i) *= P_NUM;
        }

    //fine normalizzazione

        printf("\t-> Guess density correctly normalized\n");
    
            sprintf(filen,"dens%de-guess",pn); //nomi dei file
            printRho(&D, filen);

        printf("\t-> Guess density correctly printed\n");
    

            sprintf(filen,"pot%de-guess",pn); //nomi dei file
            printPot(&D, filen);

        printf("\t-> Guess potentials correctly printed\n");

        printf("\t-> Starting the self-consistent iteration\n");
    
            //##########
            //SC_lopp(&D, EDIM);
            //#########
    
        printf("\n");
        printf("------------\n");
        printf("FINAL CHECKS\n");
        printf("-------------\n");

        sprintf(filen,"dens%de",pn); //nomi dei file
        printRho(&D, filen);

        printf("\t-> Final density correctly printed\n");

        sprintf(filen,"pot%de",pn); //nomi dei file
        printPot(&D, filen);

        printf("\t-> Final potentials correctly printed\n\n");


    return 0;
}
