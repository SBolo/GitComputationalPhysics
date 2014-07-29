#include <stdio.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_eigen.h>
#include <math.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_integration.h>
#include <armadillo>

using namespace arma;

//_____________________________________________________________GSL_MATRIX_MACROS________________________________________________

#define SOLVING(e,l,eval,deg) solve(energies,tmp,e,l,eval,deg, D)
#define STATE() system_energy(D)
#define ESUM() eigen_sum(D)

#define GSL(F,x) (*((F)->function))(x,(F)->params)
#define STEP (2ce-2)
#define RMAX (25.) //dimensione mesh
#define P_NUM (2.) //numero di elettroni nella sfera
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
    vec psi2;
    vec rho; //densità del sistema
    vec dpsi2; //derivata di psi al quadrato
    vec tmp; //per fare i conti se non voglio rovinare un array
    vec rho_tmp;
    double sphere_radius; //raggio della sfera carica
    double l; //momento angolare
} schrodinger;


//______________________________________________________________END_STRUCTURES____________________________________________________


double simpson_normalize(vec psi, int N) { //normalizzazione
	double h, k, res = 0;
	int i = 0;
    vec f(N);
	
        for(i = 0; i < N; i++) {
            f(i) = 4. * M_PI * psi(i) * pow(STEP*i, 2.);
        }
	
    k = f(0) + 4 * f(1) + f(N-1) + 4 * f(N-2);
    
        for(i = 2; i < N-2; i+=2) {
            res += 4 * f(i) + 2 * f(i+1);
        }
	
	return (k + res) * (STEP/3);
	
}

double integrate(vec g, double j, double h) { //array, primo estremo (primo indice), secondo estremo (secondo indice)

    int i;
    double risultato = 0.;
    
    for( i = j; i < h; i++) {
        risultato += g(i) * STEP;
    }
    
	return risultato;
}

double guess(double x, void *params) { //distribuzione di densità iniziale
    
    schrodinger * D = (schrodinger *)params;
    double radius = D->sphere_radius;
    double rho_plus = (3. /** P_NUM*/) / (4. * M_PI * pow(radius,3.));
    double mu = 1.; //deve essere piccolo: questo valore va approssimativamente bene
    
    double rho = rho_plus / (1. + exp(mu * (x - radius) ) );
    
    return rho; //va integrata per essere normalizzata
    
}


double correlation_exchange(double r , schrodinger *D){
    
    //schrodinger *D = (schrodinger *)params;
	int i = (int) (r/STEP);
	double rho_r = D->rho[i];
    
	double denom =  7.8 + pow((3.0 / (4.0 * M_PI * rho_r)), 1.0 / 3.0);
	double r_s = pow((3.0 / (4.0 * M_PI * rho_r)), 1.0 / 3.0);
	double term = pow((3.0 * rho_r / M_PI), 1.0 / 3.0);
    
	return  -0.44 / (7.8 + r_s) - 0.44 / 3.0 * r_s / ((7.8 + r_s) * (7.8 + r_s)) - term;
}


double coulomb(double x, schrodinger *D) { //porzione coulombiana del potenziale
    
    //schrodinger *D = (schrodinger *)params;
    double radius = D->sphere_radius;
    double rho_plus = (3. /** P_NUM*/) / (4. * M_PI * pow(RS_Li,3.));
    double ext;
    
        if( x <= radius ) {
            ext = 2. * M_PI * rho_plus * ( (1./3.)*x*x - radius*radius );
        } else {
            ext = - 4./3. * M_PI * rho_plus * pow(radius,3.)/x;
        }

    return ext;
    
}


double hartree(double x, schrodinger *D) {
    
    //schrodinger *D = (schrodinger *)params;
    int i;
    double ext, pot1, pot2;
    
        for(i = 0; i < D->dim; i++) {
            D->tmp(i) = D->rho(i) * (i * STEP);
        }
    
        pot1 = 4. * M_PI * integrate(D->tmp, x/STEP, RMAX/STEP);
    
        for(i = 0; i < D->dim; i++) {
            D->tmp(i) *= (i * STEP);
        }
    
        pot2 = 4. * M_PI * (1./(x+STEP)) * integrate(D->tmp, 0, x/STEP);
    
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


void self_consistence(schrodinger *D, int e_dim) {
    
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
    
    alpha = 0.04;
    
        for( i = 0; i < EDIM; i++ ) {
            energies(i) = 0.;
        }
    
        for(i = 0; i < D->dim; i++) {
            D->rho_tmp(i) = 0.;
            tmp2(i) = D->rho(i);
        }
    
    //per prima cosa risolvo l'equazione per il guess.
  
    //     arr l  n  deg
    SOLVING(0, 0, 0, 2); //1s
//    SOLVING(1, 1, 0, 6); //1p
//    SOLVING(2, 2, 0, 10); //1d
//    SOLVING(3, 0, 1, 2); //2s
//    SOLVING(4, 3, 0, 14); //1f
//    SOLVING(5, 1, 1, 6); //2p
    
    //in definitiva salviamo il temporaneo dentro rho e ricalcoliamo l'energia dello stato
        for(i = 0; i< D->dim; i++ ) {
            D->rho(i) = D->rho_tmp(i) * alpha + tmp2(i) * (1. - alpha);
        }

    e_new = 2.*energies(0) + 6.*energies(1) + 2.*energies(2) + 10.*energies(3) + 14.*energies(4) + 6.*energies(5); //40e
    
    //fino qui abbiamo calcolato la soluzione preliminare


            do {
                
                e_old = e_new; //salviamo in e_old la vecchia e_new
                
                for(i = 0; i < D->dim; i++) {
                    D->rho_tmp(i) = 0.;
                    tmp2(i) = D->rho(i);
                }
                
                //     arr l  n  deg
                SOLVING(0, 0, 0, 2); //1s
//                SOLVING(1, 1, 0, 6); //1p
//                SOLVING(2, 2, 0, 10); //1d
//                SOLVING(3, 0, 1, 2); //2s
//                SOLVING(4, 3, 0, 14); //1f
//                SOLVING(5, 1, 1, 6); //2p
                
                //in definitiva salviamo il temporaneo dentro rho e ricalcoliamo l'energia dello stato
                for( i = 0; i < D->dim; i++ ) {
                    D->rho(i) = tmp2(i) * (1. - alpha) + alpha * D->rho_tmp(i);
                }
                
                e_new = 2.*energies(0) + 6.*energies(1) + 2.*energies(2) + 10.*energies(3) + 14.*energies(4) + 6.*energies(5); //40e
                
                iteration++;
            
printf("# ITERATION = %d: e_old = %.4f, e_new = %.4f, diff = %.4f, 1s = %.3f, 1p = %.3f, 1d = %.3f, 2s = %.3f, 1f = %.3f, 2p = %.3f, state = %.3f\n", iteration, e_old, e_new, fabs(e_old - e_new), energies[0], energies[1], energies[2], energies[3], energies[4], energies[5], (e_new + STATE())/P_NUM);
                fflush(stdout);
            
            } while ( fabs( e_new - e_old ) > EPSILON);
    
    pt = fopen(filen,"w");
    //printf("\n");
    fprintf(pt,"-------------------------\n");
    fprintf(pt,"SUMMARY OF THE SIMULATION\n");
    fprintf(pt,"-------------------------\n");
            fprintf(pt,"#Number of particles: %g\n", P_NUM);
            fprintf(pt,"#Dimension of the matrix: %g\n", RMAX/STEP);
            fprintf(pt,"#Precision on the eigenvalues: %g\n", EPSILON);
            fprintf(pt,"#Correct sum of eigenvalues: %.4f\n", e_new);
            fprintf(pt,"#Correct single-level energies: 1s = %.3f H, %.3f eV\n                                1p = %.3f H, %.3f eV\n                                1d = %.3f H, %.3f eV\n                                2s = %.3f H, %.3f eV\n                                1f = %.3f H, %.3f eV\n                                2p = %.3f H, %.3f eV\n", energies[0], energies[0] * 27.21, energies[1], energies[1] * 27.21, energies[2], energies[2]*27.21, energies[3], energies[3]*27.21, energies[4], energies[4] * 27.21, energies[5], energies[5] * 27.21);
            fprintf(pt,"#Correct energy per particle: %.3f\n", (e_new + STATE())/P_NUM);
    //printf("\n");
    fclose(pt);
}



//N. B. : i è la riga, j la colonna
int main() {
    int i, j, eigen, sort, k = 0, pn;
    double valore, normalization, e, energy;
    FILE *pt;
    char filen[20];
    
    pn = (int)(P_NUM);
    
    //preparo la struttura shrodinger
    schrodinger D;
    D.dim = RMAX/STEP;
    
    gsl_function RHO_GUESS;
    RHO_GUESS.function = guess;
    
    RHO_GUESS.params = &D;

    printf("ciao!\n");

    vec psi(D.dim), psi2(D.dim), rho(D.dim), dpsi2(D.dim), tmp(D.dim), rho_tmp(D.dim), rho_guess(D.dim), eigVal(D.dim);
    mat eigVec(D.dim, D.dim, fill::zeros), matrix(D.dim, D.dim, fill::zeros);

    printf("ciao 2!\n");
    
    
//    double *psi = (double *)malloc((D.dim) * sizeof(double)); //funzione d'onda
//    double *psi2 = (double *)malloc((D.dim) * sizeof(double)); //funzione d'onda
//    double *rho = (double *)malloc((D.dim) * sizeof(double)); //densità del sistema
//    double *dpsi2 = (double *)malloc((D.dim) * sizeof(double)); //densità del sistema
//    double *tmp = (double *)malloc((D.dim) * sizeof(double)); //densità del sistema
//    double *rho_tmp = (double *)malloc((D.dim) * sizeof(double)); //densità del sistema
//    double *rho_guess = (double *)malloc((D.dim) * sizeof(double)); //densità del sistema
    
    D.eigVal = eigVal;
    D.matrix = matrix;
    D.eigVec = eigVec;
    D.rho = rho;
    D.psi = psi;
    D.psi2 = psi2;
    D.rho_tmp = rho_tmp;
    D.dpsi2 = dpsi2;
    D.tmp = tmp;
    D.sphere_radius = RS_Li * pow(P_NUM,(1./3.));
    
    D.l = 0.;

    printf("tutto definito!\n");
    
    //normalizzazione del guess iniziale
        for(i=0; i< D.dim; i++) {
            D.rho(i) = GSL(&RHO_GUESS,i*STEP);
        }
       printf("comicio a riempire\n");
            normalization = simpson_normalize(D.rho, D.dim);
            printf("Non ho sbagliato a integrare\n");

        for (i=0; i < D.dim; i++) {
            D.rho(i) /= normalization;
            D.rho(i) *= P_NUM;
        }
        printf("Quindi la normalizzazione funziona!\n");
    //fine normalizzazione
    
    sprintf(filen,"dens%de-guess",pn); //nomi dei file
    pt = fopen(filen,"w");
    
    for(i = 0; i < D.dim; i++) {
        fprintf(pt, "%g \t %g \t \t \n", i*STEP, D.rho(i));
    }

    printf("Ho stampato la densità\n");
    
    fclose(pt);
    
    sprintf(filen,"pot%de-guess",pn); //nomi dei file
    pt = fopen(filen,"w");
    
    for( i = 0; i < D.dim; i++ ) {
        fprintf( pt, "%g \t %g \t %g\n", i*STEP, hartree(i*STEP, &D), correlation_exchange(i*STEP, &D) ); //potenziale autoconsistente
    }
    printf("Ho pure stampato i potenziali\n");
    fclose(pt);
    
    printf("Se puoi leggere, allora l'errore è nel self-cons\n");
    
    //##########
    self_consistence(&D, EDIM);
    //#########
    
        sprintf(filen,"dens%de",pn); //nomi dei file
        pt = fopen(filen,"w");
    
            for(i = 0; i < D.dim; i++) {
                fprintf(pt, "%g \t %g \t \t \n", i*STEP, D.rho(i));
            }
    
        fclose(pt);

        sprintf(filen,"pot%de",pn); //nomi dei file
        pt = fopen(filen,"w");
    
            for( i = 0; i < D.dim; i++ ) {
                fprintf(pt, "%g \t %g \t %g\n", i*STEP, hartree(i*STEP, &D), correlation_exchange(i*STEP, &D)); //potenziale autoconsistente
            }
        fclose(pt);
   
    return 0;
}
