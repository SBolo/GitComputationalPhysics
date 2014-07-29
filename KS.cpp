#include <stdio.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_eigen.h>
#include <math.h>
#include <gsl/gsl_sf.h>

//_____________________________________________________________GSL_MATRIX_MACROS________________________________________________

#define MATRIX(n,m) gsl_matrix_alloc(n,m) 
#define GET_ELEMENT(m,i,j) gsl_matrix_get(m,i,j)
#define SET_ELEMENT(m,i,j,x) gsl_matrix_set(m,i,j,x)
#define GET_VECTOR(v,i) gsl_vector_get(v,i)
#define DIAG(x) matrix_diagonal(x,CE,C,H,D)
#define DIAG_M(x) matrix_diagonal(x,&CE,&C,&H,&D)
#define SOLVING(e,l,eval,deg) solve(D,CE,C,H,energies,tmp,e,l,eval,deg) 
#define STATE() system_energy(CE,C,H,D)
#define ESUM() eigen_sum(CE,C,H,D)

#define GSL(F,x) (*((F)->function))(x,(F)->params)
#define STEP (2.5e-2)
#define RMAX (25.) //dimensione mesh
#define P_NUM (2.) //numero di elettroni nella sfera
#define RS_Li (4.) //raggio di Wiegner Seitz per cominciare
#define EPSILON (1.e-3) //precisione sulla somma degli autovalori
#define EDIM (6) //dimensione array energie

//________________________________________________________________END_MACROS____________________________________________________


//________________________________________________________________STRUCTURES____________________________________________________

typedef struct {
	double (*function)(double x, void *params);
	void *params;
} gsl_function;


typedef struct {
    /*Per la diagonalizzazione */
    gsl_vector *eigenvalues; //vettore GSL: conterrà gli autovalori
    gsl_matrix *matrix; //matrice quadrata: da diagonalizzare
    gsl_matrix *eigenvectors; //matrice quadrata: contiene gli autovettori
    gsl_eigen_symmv_workspace *workspace; //spazio di memoria allocato per l'effettiva diagonalizzazione
    /* Parametri e altre cose */
    int dim; //dimensione della matrice
    double *psi; //funzione d'onda
    double *psi2;
    double *rho; //densità del sistema
    double *dpsi2; //derivata di psi al quadrato
    double *tmp; //per fare i conti se non voglio rovinare un array
    double *rho_tmp;
    double sphere_radius; //raggio della sfera carica
    double l; //momento angolare
} schrodinger;


//______________________________________________________________END_STRUCTURES____________________________________________________


double simpson_normalize(double *psi, int N) { //normalizzazione
	double h, k, res = 0;
	int i = 0;
	double f[N];
	
        for(i = 0; i < N; i++) {
            f[i] = 4. * M_PI * psi[i] * pow(STEP*i, 2.);
        }
	
    k = f[0] + 4 * f[1] + f[N-1] + 4 * f[N-2];
    
        for(i = 2; i < N-2; i+=2) {
            res += 4 * f[i] + 2 * f[i+1];
        }
	
	return (k + res) * (STEP/3);
	
}

double integrate(double *g, double j, double h) { //array, primo estremo (primo indice), secondo estremo (secondo indice)

    int i;
    double risultato = 0.;
    
    for( i = j; i <= h; i++) {
        risultato += g[i] * STEP;
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


double correlation_exchange(double r , void *params){
    
    schrodinger *D = (schrodinger *)params;
	int i = (int) (r/STEP);
	double rho_r = D->rho[i];
    
	double denom =  7.8 + pow((3.0 / (4.0 * M_PI * rho_r)), 1.0 / 3.0);
	double r_s = pow((3.0 / (4.0 * M_PI * rho_r)), 1.0 / 3.0);
	double term = pow((3.0 * rho_r / M_PI), 1.0 / 3.0);
    
	return  -0.44 / (7.8 + r_s) - 0.44 / 3.0 * r_s / ((7.8 + r_s) * (7.8 + r_s)) - term;
}


double coulomb(double x, void *params) { //porzione coulombiana del potenziale
    
    schrodinger *D = (schrodinger *)params;
    double radius = D->sphere_radius;
    double rho_plus = (3. /** P_NUM*/) / (4. * M_PI * pow(RS_Li,3.));
    int i;
    double ext, pot1, pot2;
    
        if( x <= radius ) {
            ext = 2. * M_PI * rho_plus * ( (1./3.)*x*x - radius*radius );
        } else {
            ext = - 4./3. * M_PI * rho_plus * pow(radius,3.)/x;
        }

    return ext;
    
}


double hartree(double x, void *params) {
    
    schrodinger *D = (schrodinger *)params;
    int i;
    double ext, pot1, pot2;
    
        for(i = 0; i < D->dim; i++) {
            D->tmp[i] = D->rho[i] * (i * STEP);
        }
    
        pot1 = 4. * M_PI * integrate(D->tmp, x/STEP, RMAX/STEP);
    
        for(i = 0; i < D->dim; i++) {
            D->tmp[i] *= (i * STEP);
        }
    
        pot2 = 4. * M_PI * (1./(x+STEP)) * integrate(D->tmp, 0, x/STEP);
    
    return (pot1 + pot2);
}


double matrix_diagonal(double x, gsl_function *CE, gsl_function *C, gsl_function *H, schrodinger *D) { //funzione per riempire la diagonale principale: potentiale + 1/(delta x)^2
        
        double l = D->l;
        double V;

        V = GSL(CE,x) + GSL(H,x) + GSL(C,x) + (1.*l) * ( (1.*l) + 1. ) / (2. * (x + STEP) * (x + STEP) );
    
    return ( V + 1./pow(STEP,2.) );
    
}


double diagonalization(gsl_function *CE, gsl_function *C, gsl_function *H, schrodinger *D, int eval) { //scambio correlazione + coulomb + hartree
 
     //in input: gsl_function per il calcolo della diagonale, schrodinger per il settaggio della matrice e numero dell'autovalore desiderato.
     //L'autovettore è contenuto nella matrice all'interno della struttura. Lo prendo quando desidero, non serve farlo ora 
    
    int i, j, eigen, sort, k = 0;
    double normalization, energy;
    
    gsl_matrix_set_zero(D->matrix);
 
    for(i = 0; i < (D->dim); i++) {
        for(j = 0; j < (D->dim); j++) {
            
            if(i == j) SET_ELEMENT( D->matrix, i, j, DIAG( i * STEP )); //diagonale principale
            if(i == (j+1) || i == (j-1)) SET_ELEMENT( D->matrix, i, j, - 1. / (2. * pow(STEP,2.) ) ); //diagonale superiore e inferiore
            
        }
    }
    
    
    gsl_eigen_symmv(D->matrix, D->eigenvalues, D->eigenvectors, D->workspace);
    gsl_eigen_symmv_sort(D->eigenvalues, D->eigenvectors, GSL_EIGEN_SORT_VAL_ASC); //ordina in ordine crescente
    
    //non normalizziamo, ma lo facciamo direttamente nella funzione self-consistence
    
    energy = GET_VECTOR(D->eigenvalues, eval); //trovo l'energia dello stato che mi interessa
    
    return energy;
    
}

double system_energy(gsl_function *CE, gsl_function *C, gsl_function *H, schrodinger *D) {
    
    double u_int, mean_uint[D->dim], rho2_de_xc[D->dim], eps, eps_xc, rho, r_s;
    int i;
    
    for( i = 0; i< D->dim; i++ ) {
        mean_uint[i] = GSL(H,i*STEP) * D->rho[i];
    }
    
    u_int = 1./2. * simpson_normalize(mean_uint, D->dim);
    
    for( i = 0; i < D->dim; i++) {
        rho = D->rho[i];
        r_s = pow((3. / (4. * M_PI * rho)), 1. / 3.);
        rho2_de_xc[i] = -1./4. * pow( (3. / M_PI), 1./3. ) * pow( rho, 4./3. ) - 1./3. * ( 0.44 / pow((7.8 + r_s),2.) ) * pow( 3. / (4.*M_PI), 1./3.) * pow(rho, 2./3.);
    }
    
    eps_xc = simpson_normalize(rho2_de_xc, D->dim);
    
    
    return ( - u_int - eps_xc );
    
}

void solve(schrodinger *D, gsl_function *CE, gsl_function *C, gsl_function *H, double *energies, double *tmp, int e, int l, int eval, int deg ) {
//struttura, array per le energie, tmp per i valori temporanei, posto nell'array, momento angolare, autovalore, degenerazione
    
    double normalization, tmp_derivate[D->dim];
    int i;
    
    D->l = (1.*l);
     energies[e] = diagonalization(CE,C,H,D,eval);
    
    for( i = 0; i < D->dim; i++ ) {
        D->psi[i] = pow( GET_ELEMENT(D->eigenvectors,i,eval) / ( (i+1)*STEP ), 2.);
    }
    
    normalization = simpson_normalize(D->psi, D->dim);
    
    for( i = 0; i < D->dim; i++ ) {
        D->psi[i] /= normalization;
        D->rho_tmp[i] += (1.*deg) * D->psi[i];
    }
}


void self_consistence(gsl_function *CE, gsl_function *C, gsl_function *H, schrodinger *D, int e_dim) {
    
    double alpha, e_new, e_old, normalization, state;
    int i, k = 0, iteration = 0;
    double tmp[D->dim], tmp2[D->dim], energies[e_dim];
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
            energies[i] = 0.;
        }
    
        for(i = 0; i < D->dim; i++) {
            D->rho_tmp[i] = 0.;
            tmp2[i] = D->rho[i];
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
            D->rho[i] = D->rho_tmp[i] * alpha + tmp2[i] * (1. - alpha);
        }

    e_new = 2.*energies[0] + 6.*energies[1] + 2.*energies[2] + 10.*energies[3] + 14.*energies[4] + 6.*energies[5]; //40e
    
    //fino qui abbiamo calcolato la soluzione preliminare


            do {
                
                e_old = e_new; //salviamo in e_old la vecchia e_new
                
                for(i = 0; i < D->dim; i++) {
                    D->rho_tmp[i] = 0.;
                    tmp2[i] = D->rho[i];
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
                    D->rho[i] = tmp2[i] * (1. - alpha) + alpha * D->rho_tmp[i];
                }
                
                e_new = 2.*energies[0] + 6.*energies[1] + 2.*energies[2] + 10.*energies[3] + 14.*energies[4] + 6.*energies[5]; //40e
                
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
    
    gsl_function RHO_GUESS, C, CE, H;
    RHO_GUESS.function = guess;
    C.function = coulomb;
    CE.function = correlation_exchange;
    H.function = hartree;
    
    RHO_GUESS.params = &D;
    C.params = &D;
    CE.params = &D;
    H.params = &D;
    
    
    double *psi = (double *)malloc((D.dim) * sizeof(double)); //funzione d'onda
    double *psi2 = (double *)malloc((D.dim) * sizeof(double)); //funzione d'onda
    double *rho = (double *)malloc((D.dim) * sizeof(double)); //densità del sistema
    double *dpsi2 = (double *)malloc((D.dim) * sizeof(double)); //densità del sistema
    double *tmp = (double *)malloc((D.dim) * sizeof(double)); //densità del sistema
    double *rho_tmp = (double *)malloc((D.dim) * sizeof(double)); //densità del sistema
    double *rho_guess = (double *)malloc((D.dim) * sizeof(double)); //densità del sistema
    
    D.eigenvalues = gsl_vector_alloc(D.dim);
    D.matrix = MATRIX(D.dim,D.dim);
    D.eigenvectors = MATRIX(D.dim,D.dim);
    D.workspace = gsl_eigen_symmv_alloc(D.dim);
    D.rho = rho;
    D.psi = psi;
    D.psi2 = psi2;
    D.rho_tmp = rho_tmp;
    D.dpsi2 = dpsi2;
    D.tmp = tmp;
    D.sphere_radius = RS_Li * pow(P_NUM,(1./3.));
    
    gsl_matrix_set_zero(D.matrix); //setto la matrice a 0
    D.l = 0.;
    
    //normalizzazione del guess iniziale
        for(i=0; i<D.dim; i++) {
            D.rho[i] = GSL(&RHO_GUESS,i*STEP);
        }
            normalization = simpson_normalize(D.rho, D.dim);

        for (i=0; i < D.dim; i++) {
            D.rho[i] /= normalization;
            D.rho[i] *= P_NUM;
        }
    //fine normalizzazione
    
    sprintf(filen,"dens%de-guess",pn); //nomi dei file
    pt = fopen(filen,"w");
    
    for(i = 0; i < D.dim; i++) {
        fprintf(pt, "%g \t %g \t \t \n", i*STEP, D.rho[i]);
    }
    
    fclose(pt);
    
    sprintf(filen,"pot%de-guess",pn); //nomi dei file
    pt = fopen(filen,"w");
    
    for( i = 0; i < D.dim; i++ ) {
        fprintf(pt, "%g \t %g \t %g\n", i*STEP, GSL(&H,i*STEP), GSL(&CE,i*STEP)); //potenziale autoconsistente
    }
    fclose(pt);
    
    
    //##########
    self_consistence(&CE, &C, &H, &D, EDIM);
    //#########
    
        sprintf(filen,"dens%de",pn); //nomi dei file
        pt = fopen(filen,"w");
    
            for(i = 0; i < D.dim; i++) {
                fprintf(pt, "%g \t %g \t \t \n", i*STEP, D.rho[i]);
            }
    
        fclose(pt);

        sprintf(filen,"pot%de",pn); //nomi dei file
        pt = fopen(filen,"w");
    
            for( i = 0; i < D.dim; i++ ) {
                fprintf(pt, "%g \t %g \t %g\n", i*STEP, GSL(&H,i*STEP), GSL(&CE,i*STEP)); //potenziale autoconsistente
            }
        fclose(pt);
   
    
    
//_______________________________________________________________FREE_VECTORS_____________________________________________________
    
    gsl_vector_free(D.eigenvalues);
    gsl_matrix_free(D.matrix);
    gsl_matrix_free(D.eigenvectors);
    gsl_eigen_symmv_free(D.workspace);
    
    free(psi);
    free(rho);
    free(tmp);
    free(dpsi2);
    free(psi2);
    free(rho_tmp);
    
//________________________________________________________________END_FREE________________________________________________________
    
    
    return 0;
}
