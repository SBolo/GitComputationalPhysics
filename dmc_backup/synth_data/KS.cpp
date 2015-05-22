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

#define GSL(F,x) (*((F)->function))(x,(F)->params)

//________________________________________________________________END_MACROS____________________________________________________


//________________________________________________________________STRUCTURES____________________________________________________

typedef struct {
    /*Per la diagonalizzazione */
    double step;
    double rmax;
    mat matrix; //matrice quadrata: da diagonalizzare
    vec eigVal; //vettore GSL: conterrà gli autovalori
    mat eigVec; //matrice quadrata: contiene gli autovettori
    /* Parametri e altre cose */
    int dim; //dimensione della matrice
} schrodinger;


//______________________________________________________________END_STRUCTURES____________________________________________________

double matrix_diagonal(double r, schrodinger *D) {
    //funzione per riempire la diagonale principale: potentiale + 1/(delta x)^2

    return ( -1./(r+D->step) - 1./(1. + r * r * exp(r - 10)) );
}


void diagonalization(schrodinger *D) { //scambio correlazione + coulomb + hartree
 
     //in input: gsl_function per il calcolo della diagonale, schrodinger per il settaggio della matrice e numero dell'autovalore desiderato.
     //L'autovettore è contenuto nella matrice all'interno della struttura. Lo prendo quando desidero, non serve farlo ora 
    
    D->matrix.fill(0.);
 
    for(int i = 0; i < (D->dim); i++) {
        for(int j = 0; j < (D->dim); j++) {
            
            if(i == j) D->matrix(i,j) = matrix_diagonal(i*D->step, D) + 1./pow(D->step,2.);//diagonale principale
            if(i == (j+1) || i == (j-1)) D->matrix(i, j) = - 1. / ( 2. * pow(D->step,2.) ); //diagonale superiore e inferiore
            
        }
    }
    
    eig_sym(D->eigVal, D->eigVec, D->matrix);
    
}


int main() {

    //preparo la struttura shrodinger
    schrodinger D;
    FILE *pt;
    double c = -2.127;
    float energies[20];

    D.step = 0.07;
    D.rmax = 600;
    D.dim = D.rmax/D.step;

    vec eigVal(D.dim);
    mat eigVec(D.dim, D.dim, fill::zeros), matrix(D.dim, D.dim, fill::zeros);
    
    D.eigVal = eigVal;
    D.matrix = matrix;
    D.eigVec = eigVec;


    //diagonalization(&D);

    pt = fopen("b_energies", "rt");

   for(int i = 0; i < 20; i++) {
       fscanf(pt, "%f\n", &(energies[i]));
   }

   fclose(pt);

   pt = fopen("results", "w");

   fprintf(pt,"#Orbital \t BE \t Coulomb \t Delta \t |DE/E|\n");

    for(int i = 0; i < 20; i++) {
        fprintf(pt,"%d \t %.4f \t %.4f\t %.4f \t %.4f\n",i+1, -energies[i], 1. / (2. * (i+1) * (i+1)), 1. / (2. * (i+1) * (i+1)) +
                -c * 1./(sqrt(M_PI) * (i+1) * (i+1) * (i+1)), fabs( ( -energies[i] +
              1. / (2. * (i+1) * (i+1)) - c * 1./(sqrt(M_PI) * (i+1) * (i+1) * (i+1)) ) / energies[i] )  );
    }

    fclose(pt);



    return 0;
}
