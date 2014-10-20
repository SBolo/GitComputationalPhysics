#include <stdio.h>
#include <math.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_laguerre.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_eigen.h>

//#define STEP 0.00001 //step sulla derivazione
#define STEP 0.001
#define STP 5e-4 //step sul plot
#define RMAX 10 //massima distanza
#define TOLERANCE 0.005
#define E_CONV 82.95 //coefficiente di conversione energia
#define GSL(F,x) (*((F)->function))(x,(F)->params)
#define L(n,a,x)  gsl_sf_laguerre_n(n,a,x)
#define FT(n) gsl_sf_fact(n) //vuole un unsigned int: fattoriale - ritorna double
#define SFT(n) gsl_sf_doublefact(n) //di nuovo unsigned: semifattoriale - ritorna double
#define HIJ(x) Hij(F,V,P,x,i,j) //funzione integranda
#define S(i,j,a,b) simpson(F,V,P,i,j,a,b) //simpson per Hij
#define MATRIX(n,m) gsl_matrix_alloc(n,m) //genera la matrice
#define GET_ELEMENT(m,i,j) gsl_matrix_get (m,i,j) //ottenere l'elemento ij-esimo della matrice
#define SET_ELEMENT(m,i,j,x) gsl_matrix_set(m,i,j,x) //definire l'elemento ij-esimo della matrice
#define GET_VECTOR(v,i) gsl_vector_get(v,i) //estrae l'elemento i-esimo dal vettore


typedef struct {
	double (*function)(double x, void *params);
	void *params;
} gsl_function;

typedef struct {
    int n; //numero quantico principale
    int l; //momento angolare del sistema
    int k; //momento angolare base
    int dim; //dimensione matrice
    double *norm;
    gsl_vector *eigenvalues; //vettore GSL: conterrà gli autovalori
    gsl_matrix *matrix; //matrice quadrata: da diagonalizzare
    gsl_matrix *eigenvectors; //matrice quadrata: contiene gli autovettori
    gsl_eigen_symmv_workspace *workspace; //spazio di memoria allocato per l'effettiva diagonalizzazione
} prms;

/*__________________________________________________FUNCTIONS_FOR_ANALYSIS________________________________________________*/
double d1(gsl_function *f, double x) { //derivata prima
    return (GSL(f,x - 2.*STEP) - 8. * GSL(f, x - STEP) + 8. * GSL(f, x + STEP) - GSL(f, x + 2. * STEP))/(12. * STEP);
}

double d2(gsl_function *f, double x) { //derivata seconda
    return (-GSL(f,x - 2.*STEP) + 16. * GSL(f, x - STEP) - 30. * GSL(f,x) + 16. * GSL(f, x + STEP) - GSL(f, x + 2. * STEP))/(12. * STEP * STEP);
}
/*______________________________________________________END_ANALYSIS____________________________________________________*/


double basis(double x, void *params) { //Base dell'oscillatore armonico 3D
    
    prms *P = (prms *)params;
    int n = P->n; //elemento da pescare
    double j = P->norm[n]; //elemento dell'array dei coefficienti di normalizzazione
    int l = P->k;
    
    return j * exp(-x*x/2.) * L(n, 1./2. + l, x*x) * pow(x, l);
}

double potential(double x, void *params) { //potenziale di Malfliet-Tjon
    prms *P = (prms *)params;
    int l = P->l;
    
    return (2.38/x) * ( 7.39 * exp(-3.110 * x) - 3.22 * exp(-1.555 * x) ) + l * (l + 1)/(2. * x * x);
}

double Hij(gsl_function *F, gsl_function *V, prms *P, double x, int i, int j) { //hamiltoniana
    
    P->n = j; //fisso il primo indice j per il calcolo delle derivate
    double dder = - (1./2.) * d2(F,x);
    double der = - (1./x) * d1(F,x);
    double pot = GSL(V,x) * GSL(F,x);
    
    P->n = i; //fisso il secondo indice i e rimoltiplico tutto
    
    return x * x * GSL(F,x) * (dder + der + pot);
    
}

double simpson(gsl_function *F, gsl_function *V, prms *P, int i, int j, double a, double b) {
    //integra in particolare l'hamiltoniana: non molto elegante, ma funziona
    
    //a,b estremi di integrazione
	double h = 0, somma = 0, m = 0, res = 0, risultato;
	int k = 0, t = 0, N = 1e3; //N = numero di campionamenti
	double f[N];
	
        h = (b-a)/(N); //passo di integrazione
	
    for(k = 0; k < N; k++) {
        f[k] = HIJ(a+k*h); //campionamento della funzione
    }
    
        m = f[0] + 4 * f[1] + f[N-1] + 4 * f[N-2];
    
    for(k = 2; k < N - 2; k += 2) {
        res += 4 * f[k] + 2 * f[k+1];
    }
	
        risultato = (m + res) * (h/3); //risultato dell'integrale
    
	return risultato;
}

double diagonalization(gsl_function *F, gsl_function *V, prms *P, double a, double b, int eval) { //routine di diagonalizzazione
    
    int h, k, eigen, sort;
    double norm, energy;
    
    gsl_matrix_set_zero(P->matrix);
    
    printf("\n");
    printf("\t#Riempio la matrice %d x %d\n", P->dim,P->dim);
    
    for(h = 0; h < P->dim; h++) {
        
        for(k = 0; k < P->dim; k++) {
            
            printf("\t%d posti riempiti su %d\r",(P->dim)*h + k, (P->dim)*(P->dim));
            fflush(stdout);
            
            SET_ELEMENT(P->matrix, h, k, S(h,k,a,b)); //inserisco gli elementi di matrice
            //printf("Elemento [%d,%d] = %g\n",h,k,S(h,k,a,b));
            
        }
    }
    
    printf("\t#Matrice riempita: procedo alla diagonalizzazione\n");
    
        eigen = gsl_eigen_symmv(P->matrix, P->eigenvalues, P->eigenvectors, P->workspace); //diagonalizzo
        sort = gsl_eigen_symmv_sort(P->eigenvalues, P->eigenvectors, GSL_EIGEN_SORT_VAL_ASC); //ordino gli autovalori in ordine crescente
    
    printf("\t#Diangoanlizzata e ordinata: ecco il risultato, con P.dim = %d: \n", P->dim);
    
    return GET_VECTOR(P->eigenvalues, eval); //ritorna l'autovalore desiderato
    
}

double print(prms *P, gsl_function *F, double x) { //per il print della funzione d'onda in modulo quadro
    int h;
    double f = 0.;
    
    for(h = 0; h < P->dim; h++) {
        P->n = h; 
        f += GET_ELEMENT(P->eigenvectors,h,0) * GSL(F,x);
    }
    
    return f*f;
}

void alloc(prms *P) {
    
    P->eigenvalues = gsl_vector_alloc(P->dim);
    P->matrix = MATRIX(P->dim, P->dim);
    P->eigenvectors = MATRIX(P->dim,P->dim);
    P->workspace = gsl_eigen_symmv_alloc(P->dim);
    
}

void dealloc(prms *P) {
    
    gsl_vector_free(P->eigenvalues);
    gsl_matrix_free(P->matrix);
    gsl_matrix_free(P->eigenvectors);
    gsl_eigen_symmv_free(P->workspace);
    
}



int main() {
    
    double x, integ, integ2, energy, energy1, energy2, energy3, df1, df2, df3, df4, df5, df6, norm[130];
    unsigned int i;
    int j, ex;
    double *psi; //funzione d'onda
    FILE *pt, *pt2;
    char filename[30], filename2[30];
    
    gsl_function F, V;
    prms P;
    
    F.function = basis;
    F.params = &P;
    V.function = potential;
    V.params = &P;
    P.norm = norm;
    P.l = 0;
    P.k = 0;
    
    sprintf(filename, "energy_lb%d_ls%d",P.k,P.l);
    sprintf(filename2, "wave_lb%d_ls%d",P.k,P.l);
    
    for( i = 0; i < 130; i++ ) { //costruisco un array con i coefficienti di normalizzazione per la base
        norm[i] = sqrt( sqrt( 1./(4. * M_PI) ) * pow(2.,i + 3 + 2. * P.k) * FT(i) * pow(1./2., P.k) / SFT(2 * i + 2. * P.k + 1) ); //coefficienti di normalizzazione della base
        //printf("%g\n", norm[i]);
    }
    
    printf("\n");
    printf("\t Risoluzione con l_base = %d e l_sis = %d\n",P.k,P.l);
    printf("\t Precisione: TOLERANCE = %g\n", TOLERANCE);
    
    pt = fopen(filename,"w");
    
    //set di condizioni con cui cominciare lo sviluppo
    P.dim = 1; //modificare qui per cambiare la dimensione della matrice
    //numeri qualunque per l'energia iniziale
    energy = 1.;
    energy1 = 1.1;
    energy2 = 1.2;
    energy3 = 1.3;
    ex = 0; //per cominciare il ciclo
    
    
    do {
        
        energy3 = energy2;
        energy2 = energy1;
        energy1 = energy*E_CONV;
        
        alloc(&P);
    
            energy = diagonalization(&F, &V, &P, STEP, 20., 0);
            printf("\t#ENERGIA: %.4f\n", energy*E_CONV);
        
        df1 = fabs(energy*E_CONV - energy1);
        df2 = fabs(energy*E_CONV - energy2);
        df3 = fabs(energy*E_CONV - energy3);
        df4 = fabs(energy1 - energy2);
        df5 = fabs(energy2 - energy3);
        df6 = fabs(energy1 - energy3);
        
            fprintf(pt,"%d \t %.4f \t %.4f \t %.4f \t %.4f \t %.4f \t %.4f \t %.4f\n", P.dim, energy*E_CONV, df1, df2, df3, df4, df5, df6);
    
            //printf("\t#ENERGY1: %.4f\n", energy1);
            //printf("\t#ENERGY2: %.4f\n", energy2);
            //printf("\t#ENERGY3: %.4f\n", energy3);
            //printf("\t#DIFFERENZE: %.4f, %.4f, %.4f, %.4f, %.4f, %.4f\n", df1, df2, df3, df4, df5, df6);
            //printf("\n");
     
        
            //if( df1 > TOLERANCE || df2 > TOLERANCE || df3 > TOLERANCE || df4 > TOLERANCE || df5 > TOLERANCE || df6 > TOLERANCE  ) {
        
        //if( (fabs(energy) - 2.22) < TOLERANCE) {
        
        if( P.dim <= 50) {
                //condizioni di uscita dal ciclo
                
                    ex = 1;
                
            } else { //se esco dal ciclo, allora plotto la corrispondente psi^2
                
                    ex = 0;
                    psi = (double *)malloc((RMAX/STP) * sizeof(double));
                
                        pt2 = fopen(filename2,"w");
                
                            for( j = 0; j < RMAX/STP; j++ ) {
                                psi[j] = print(&P,&F,j*STP);
                                fprintf(pt2,"%g \t %g\n", j*STP, psi[j]);
                            }
                
                        fclose(pt2);
                    free(psi);
            }
        
        (P.dim)++;
        dealloc(&P);
        
    } while ( ex == 1 );
        printf("Energia definitiva: %.4f\n", energy*E_CONV);
    
    
    fclose(pt);
    
    
    return 0;
}