#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>


//_______________________________________________BEGIN_MACROS_________________________________________________
#define FN(F,x) (*((F)->f))(x,(F)->params) //la mia macro per GSL function
#define XMAX 10.
//_______________________________________________END_MACROS___________________________________________________


//______________________________________________BEGIN_STRUCTS_________________________________________________
typedef struct { //struttura per una funzione
    double (*f)(double x, void *params);
	void *params;
} function;

typedef struct {
    double x0; //x_0 della gaussiana
    double x1; //x_1 dell'onda quadra
    double x2; //x_2 dell'onda quadra
} f_params;

typedef struct { //struttura per la risoluzione dell'equazione
    int time; //tempo della simulazione
    int dim; //A->dimensione degli array
    double xstep; //passo spaziale
    double tstep; //passo temporale
    double cf; //fattore di Courant
    double *solution; //array della soluzione
    double *tsol; //array tmp per la soluzione
    double *tsol_1; //secondo array tmp per la soluzione
} aveq;
//_______________________________________________END_STRUCTS___________________________________________________


double gaussian(double x, void *params) { //era initial_condition
    f_params *p = (f_params *)params;

    return exp( - (x - p->x0)*(x - p->x0) ) ;
}

double square_wave(double x, void *params) {
    f_params *p = (f_params *)params;

        if( x > p->x1 && x < p->x2 ) {
            return 1.;
        } else {
            return 0.;
        }
}

void fill_array_ic(aveq *A, function *G, function *S, int k) { //intero per decidere se utilizzare onda quadra o gaussiana

    int i;

        if( k == 1 ) { //condizione fissata per la gaussiana
            for(i = 0; i < A->dim; i++)
                    A->tsol[i] = A->solution[i] = A->tsol_1[i] = FN(G,i*A->xstep);

        } else if ( k == 0 ) { //condizione fissata per l'onda quadra
            for(i = 0; i < A->dim; i++)
                   A->tsol[i] = A->solution[i] = A->tsol_1[i] = FN(S,i*A->xstep);

        } else { //condizione a casaccio: errore
            printf("#Errore: inserire una condizione valida (0 o 1).\n");
            exit(0); //blocca esecuzione del programma??
        }

}

void FTCS(aveq *A, function *G, function *S, int k) { //struttura dei dati + intero per decidere se utilizzare onda quadra o gaussiana
    
    int i, t = 0, T = (int)( (A->time)/(A->tstep) );

    fill_array_ic(A, G, S, k); //riempiamo l'array con la condizione iniziale stabilita da k

        do {
                for( i = 1; i < A->dim-1; i++ ) {
                    A->solution[i] = A->tsol[i] - (A->cf/2.) * (A->tsol[i+1] - A->tsol[i-1]); //algoritmo
                }

            //condizioni al contorno periodiche
            A->solution[0] = A->tsol[0] - (A->cf/2.) * (A->tsol[1] - A->tsol[A->dim-1]);
            A->solution[A->dim-1] = A->tsol[A->dim-1] - (A->cf/2.) * (A->tsol[0] - A->tsol[A->dim-2]);
        
                for ( i = 0; i < A->dim; i++ ) {
                    A->tsol[i] = A->solution[i];
                }

            t++;

        } while(t <= T);
}


void LeFr(aveq *A, function *G, function *S, int k) {
    
    int i, t, T = (int)( (A->time)/(A->tstep) );

    fill_array_ic(A, G, S, k);//riempiamo l'array con la condizione iniziale stabilita da k

        for( t = 0; t <= T; t++ ) {
        
            if( t == 0) { } else if (t == 1) { } else {
            
                for( i = 1; i < A->dim-1; i++ ) {
                    A->solution[i] = A->tsol_1[i] - A->cf * (A->tsol[i+1] - A->tsol[i-1]); //algoritmo
                }
            
                //condizioni al bordo periodiche
                A->solution[0] = A->tsol_1[0] - A->cf * (A->tsol[1] - A->tsol[A->dim-1]);
                A->solution[A->dim-1] = A->tsol_1[A->dim-1] - A->cf * (A->tsol[0] - A->tsol[A->dim-2]);
            
                for ( i = 0; i < A->dim; i++ ) {
                    A->tsol_1[i] = A->tsol[i];
                    A->tsol[i] = A->solution[i];
                }
            }
        }
}


void LW(aveq *A, function *G, function *S, int k) {
    
    int i, t, T = (int)( (A->time)/(A->tstep) );
    
    fill_array_ic(A, G, S, k); //riempiamo l'array con la condizione iniziale stabilita da k
    
        for( t = 0; t <= T; t++ ) {
        
                for( i = 1; i < A->dim-1; i++ ) { //algoritmo
                    A->solution[i] = A->tsol[i] - (A->cf/2.) * (A->tsol[i+1] - A->tsol[i-1]) + A->cf * (A->cf/2.) * (A->tsol[i+1] - 2.*A->tsol[i] + A->tsol[i-1]);
                }
            //condizioni al contorno periodiche
            A->solution[0] = A->tsol[0] - (A->cf/2.) * (A->tsol[1] - A->tsol[A->dim-1]) + A->cf * (A->cf/2.) * (A->tsol[1] - 2.*A->tsol[0] + A->tsol[A->dim-1]);
            A->solution[A->dim-1] = A->tsol[A->dim-1] - (A->cf/2.) * (A->tsol[0] - A->tsol[A->dim-2]) + A->cf * (A->cf/2.) * ( A->tsol[0] - 2.* A->tsol[A->dim-1] + A->tsol[A->dim-1]);
        
                for ( i = 0; i < A->dim; i++ ) {
                    A->tsol[i] = A->solution[i];
                }
        }

}


void LF(aveq *A, function *G, function *S, int k) {
    
    int i, t = 0, T = (int)( (A->time)/(A->tstep) );
    
    fill_array_ic(A, G, S, k); //riempiamo l'array con la condizione iniziale stabilita da k
    
    
        do {
        
            if( t == 0 ) A->solution[i] = A->tsol[i];
            
            else {
        
                    for( i = 1; i < A->dim-1; i++ ) { //algoritmo
                        A->solution[i] = 1./2. * (A->tsol[i-1] + A->tsol[i+1]) - (A->cf/2.) * (A->tsol[i+1] - A->tsol[i-1]);
                    }

                //condizioni al bordo periodiche
                A->solution[0] = 1./2. * (A->tsol[A->dim-1] + A->tsol[1]) - (A->cf/2.) * (A->tsol[1] - A->tsol[A->dim-1]);
                A->solution[A->dim-1] = 1./2. * (A->tsol[A->dim-2] + A->tsol[0]) - (A->cf/2.) * (A->tsol[0] - A->tsol[A->dim-2]);
        
                    for ( i = 0; i < A->dim; i++ ) {
                        A->tsol[i] = A->solution[i];
                    }
            }
            
            t++;
            
        } while( t <= T );
    
}

/*
void L2norm(aveq *A, function *G, function *S, int k, int method) {

    double norm = 0.;
    int T = (long int)( 20./(A->tstep) ), t = 0, i;
    char filename[50];
    FILE *pt;

    double *solution = (double *)malloc((A->dim) * sizeof(double)); //soluzione
    double *tsol = (double *)malloc((A->dim) * sizeof(double)); //tmp soluzione
    double *tsol_1 = (double *)malloc((A->dim) * sizeof(double)); //tmp soluzione3

    //associa agli elementi della struct
    A->solution = solution;
    A->tsol = tsol;
    A->tsol_1 = tsol_1;

        //quale metodo di integrazione usiamo + creazione file all'interno della carrtella corrispondente
        //la cartella esiste sempre in quanto creata prima dell'esecuzione della norma L2 nel main
        //automaticamente settate condizioni iniziali dal parametro k
        switch(method) {
            case 0:
                if(k == 1) snprintf(filename, sizeof(filename), "FTCS_G/ftcs_l2norm_k%d", k);
                else snprintf(filename, sizeof(filename), "FTCS_S/ftcs_l2norm_k%d", k);
                pt = fopen(filename,"w");
                break;
            case 1:
                if(k == 1) snprintf(filename, sizeof(filename), "LEAPFROG_G/lpfr_l2norm_k%d", k);
                else snprintf(filename, sizeof(filename), "LEAPFROG_S/lpfr_l2norm_k%d",k);
                pt = fopen(filename,"w");
                break;
            case 2:
                if( k== 1) snprintf(filename, sizeof(filename), "LAXWENDROFF_G/lw_l2norm_k%d", k);
                else snprintf(filename, sizeof(filename), "LAXWENDROFF_S/lw_l2norm_k%d", k);
                pt = fopen(filename,"w");
                break;
            case 3:
                if( k== 1) snprintf(filename, sizeof(filename), "LAXFRIEDRICHS_G/lf_l2norm_k%d", k);
                else snprintf(filename, sizeof(filename), "LAXFRIEDRICHS_S/lf_l2norm_k%d", k);
                pt = fopen(filename,"w");
                break;

        }

        printf("Finito lo switch\n");

           // for(t = 0; t <= T; t++) {

        do {
                printf("Sono nel ciclo del t\n");

                A->time = t * A->tstep;
                printf("A->tstep = %g", A->tstep);
                printf("A->time = %d\n", A->time);
                printf("T = %d\n", T);

                switch(method) {
                    case 0:
                        FTCS(A, G, S, k);
                        break;
                    case 1:
                        LeFr(A, G, S, k);
                        break;
                    case 2:
                        LW(A, G, S, k);
                        break;
                    case 3:
                        LF(A, G, S, k);
                        break;
                }

                for( i = 0; i < A->dim; i++ ) {
                    norm += A->solution[i] * A->solution[i];
                }
                fprintf(pt,"%g \t %g\n", (t * A->tstep), sqrt((1./A->dim)*norm));
                norm = 0;

                t++;

      } while( t <= T);
   // }

    fclose(pt);
    free(solution);
    free(tsol);
    free(tsol_1);



}

void LF_wave(gsl_function *IC, double *solution, double time) {
    
    double r[A->dim], s[A->dim], tmpr[A->dim], tmps[A->dim], tsol[A->dim];
    int i, t = 0, T = (int)(time/TSTEP);
    
    for(i = 0; i < A->dim; i++ ) {
        tmps[i] = 0.;
        tmpr[i] = 2. * (i*XSTEP - 5.) * GSL(IC,i*XSTEP);
        
        tsol[i] = GSL(IC,i*XSTEP);
    }
    
    do {
        
        for( i = 1; i < A->dim-1; i++ ) {
            r[i] = 1./2. * (tmpr[i+1] + tmpr[i-1]) - CF/2. * ( tmps[i+1] - tmps[i-1] );
            s[i] = 1./2. * (tmps[i+1] + tmps[i-1]) - CF/2. * ( tmpr[i+1] - tmpr[i-1] );
            solution[i] = tsol[i] + TSTEP/2. * ( s[i] + tmps[i] );
        }
                
        r[0] = 1./2. * (tmpr[1] + tmpr[A->dim-1]) - CF/2. * ( tmps[1] - tmps[A->dim-1] );
        r[A->dim-1] = 1./2. * (tmpr[0] + tmpr[A->dim-2]) - CF/2. * ( tmps[0] - tmps[A->dim-2] );
        
        s[0] = 1./2. * (tmps[1] + tmps[A->dim-1]) - CF/2. * ( tmpr[1] - tmpr[A->dim-1] );
        s[A->dim-1] = 1./2. * (tmps[0] + tmps[A->dim-2]) - CF/2. * ( tmpr[0] - tmpr[A->dim-2] );
        
        solution[0] = tsol[0] + TSTEP/2. * ( s[0] + tmps[0] );
        solution[A->dim-1] = tsol[A->dim-1] + TSTEP/2. * ( s[A->dim-1] + tmps[A->dim-1] );
        
        for ( i = 0; i < A->dim; i++ ) {
            tsol[i] = solution[i];
            tmps[i] = s[i];
            tmpr[i] = r[i];
        }
        
        t++;
        
    } while( t <= T );
    
    
}

void LW_wave(gsl_function *IC, double *solution, double time) {
    
    double r[A->dim], s[A->dim], tmpr[A->dim], tmps[A->dim], tsol[A->dim];
    int i, t = 0, T = (int)(time/TSTEP);
    
    for(i = 0; i < A->dim; i++ ) {
        tmps[i] = 0.;
        tmpr[i] = 2. * (i*XSTEP - 5.) * GSL(IC,i*XSTEP);
       // tmpr[i] = ( GSL(IC, i*XSTEP + XSTEP) - GSL(IC,i*XSTEP) ) / XSTEP;
        
        tsol[i] = GSL(IC,i*XSTEP);
    }
    
    do {
        
        for( i = 1; i < A->dim-1; i++ ) {
            r[i] = tmpr[i] - CF/2. * ( tmps[i+1] - tmps[i-1] ) + CF * CF/2. * ( tmpr[i+1] - 2. * tmpr[i] + tmpr[i-1] );
            s[i] = tmps[i] - CF/2. * ( tmpr[i+1] - tmpr[i-1] ) + CF * CF/2. * ( tmps[i+1] - 2. * tmps[i] + tmps[i-1] );
            solution[i] = tsol[i] + TSTEP/2. * ( s[i] + tmps[i] );
        }
        
        r[0] = tmpr[0] - CF/2. * ( tmps[1] - tmps[A->dim-1] ) + CF * CF/2. * ( tmpr[1] - 2. * tmpr[0] + tmpr[A->dim-1] );
        r[A->dim-1] = tmpr[A->dim-1] - CF/2. * ( tmps[0] - tmps[A->dim-2] ) + CF * CF/2. * ( tmpr[0] - 2. * tmpr[A->dim-1] + tmpr[A->dim-2] );
        
        s[0] = tmps[0] - CF/2. * ( tmpr[1] - tmpr[A->dim-1] ) + CF * CF/2. * ( tmps[1] - 2. * tmps[0] + tmps[A->dim-1] );
        s[A->dim-1] = tmps[A->dim-1] - CF/2. * ( tmpr[0] - tmpr[A->dim-2] ) + CF * CF/2. * ( tmps[0] - 2. * tmps[A->dim-1] + tmps[A->dim-2] );
        
        //solution[0] = tsol[0] + TSTEP/2. * ( s[0] + tmps[0] );
        //solution[A->dim-1] = tsol[A->dim-1] + TSTEP/2. * ( s[A->dim-1] + tmps[A->dim-1] );
        
        for ( i = 0; i < A->dim; i++ ) {
            //tsol[i] = solution[i];
            tmps[i] = s[i];
            tmpr[i] = r[i];
        }
        
        t++;
        
    } while( t <= T );
    
    
}
*/

void simulation(aveq *A, function *G, function *S, int k, int method) { //l'intero method serve a designare il metodo di integrazione
    // 0 = FTCS
    // 1 = LeapFrog
    // 2 = LaxWendroff
    // 3 = LaxFriedrichs

    int i, t, it;
    char filename[50];
    FILE *pt;

    for(A->cf = 0.5; A->cf <= 1; A->cf += 0.5) { //ciclo su due diversi fattori di Courant

        for(A->time = 10; A->time <= 20; A->time += 10) { //ciclo su due diversi tempi di esecuzione

            A->xstep = 1e-1; //step di partenza

            for(it = 1; it <= 3; it++) { //ciclo sulle dimensioni dello step (e conseguentemente dell'array)

                A->tstep = A->cf * A->xstep; //definizione dello step temporale

                //-----------------------------definizione degli array
                A->dim = (int)(XMAX/A->xstep+1); //dimensione degli array

                //malloc non perchè grandi ma per evitare problemi di allocamento di memoria inutile nei cicli
                double *solution = (double *)malloc((A->dim) * sizeof(double)); //soluzione
                double *tsol = (double *)malloc((A->dim) * sizeof(double)); //tmp soluzione
                double *tsol_1 = (double *)malloc((A->dim) * sizeof(double)); //tmp soluzione3

                //associa agli elementi della struct
                A->solution = solution;
                A->tsol = tsol;
                A->tsol_1 = tsol_1;
                //-----------------------------

                printf("\t-> Simulating: cf = %g, execution time = %d, spatial step = %g\n", A->cf, A->time, A->xstep);

                //genera i file nelle cartelle che vengono create ed esegui in base al valore di method
                switch(method) {
                    case 0:
                        if(k == 1) snprintf(filename, sizeof(filename), "FTCS_G/ftcs_cf%g_xstep%g_T%d_k%d", A->cf, A->xstep, A->time,k);
                        else snprintf(filename, sizeof(filename), "FTCS_S/ftcs_cf%g_xstep%g_T%d_k%d", A->cf, A->xstep, A->time,k);
                        pt = fopen(filename,"w");
                        FTCS(A, G, S, k); //FTCS
                        break;
                    case 1:
                        if(k == 1) snprintf(filename, sizeof(filename), "LEAPFROG_G/lpfr_cf%g_xstep%g_T%d_k%d", A->cf, A->xstep, A->time,k);
                        else snprintf(filename, sizeof(filename), "LEAPFROG_S/lpfr_cf%g_xstep%g_T%d_k%d", A->cf, A->xstep, A->time,k);
                        pt = fopen(filename,"w");
                        LeFr(A, G, S, k); //Leapfrog
                        break;
                    case 2:
                        if( k== 1) snprintf(filename, sizeof(filename), "LAXWENDROFF_G/lw_cf%g_xstep%g_T%d_k%d", A->cf, A->xstep, A->time,k);
                        else snprintf(filename, sizeof(filename), "LAXWENDROFF_S/lw_cf%g_xstep%g_T%d_k%d", A->cf, A->xstep, A->time,k);
                        pt = fopen(filename,"w");
                        LW(A, G, S, k); //Lax Wendroff
                        break;
                    case 3:
                        if( k== 1) snprintf(filename, sizeof(filename), "LAXFRIEDRICHS_G/lf_cf%g_xstep%g_T%d_k%d", A->cf, A->xstep, A->time,k);
                        else snprintf(filename, sizeof(filename), "LAXFRIEDRICHS_S/lf_cf%g_xstep%g_T%d_k%d", A->cf, A->xstep, A->time,k);
                        pt = fopen(filename,"w");
                        LF(A, G, S, k); //Lax Friedrichs
                        break;

                }

                    //stampa finalmente i valori di solution su file
                    for(i = 0; i < A->dim; i++)
                        fprintf(pt, "%g \t %g\n", i*A->xstep, A->solution[i]);

                //libera tutto
                fclose(pt);
                free(solution);
                free(tsol);
                free(tsol_1);

                A->xstep /= 10.; //dividi per 10 lo step ogni volta che it viene incrementato di uno

           }
        }
     }


}

void create_directory(int method, int k) {

    switch(method) {
        case 0:
                printf("\t-> Integration method: FTCS\n");
                if(k == 1) mkdir("FTCS_G", S_IRWXU | S_IRWXG | S_IRWXO); //crea una cartella per la gaussiana
                else mkdir("FTCS_S", S_IRWXU | S_IRWXG | S_IRWXO); //crea una cartella per l'onda quadra
                break;
        case 1:
                printf("\t-> Integration method: Leapfrog\n");
                if(k == 1) mkdir("LEAPFROG_G", S_IRWXU | S_IRWXG | S_IRWXO);
                else mkdir("LEAPFROG_S", S_IRWXU | S_IRWXG | S_IRWXO);
                break;
        case 2:
                printf("\t-> Integration method: Lax-Wendroff\n");
                if(k == 1) mkdir("LAXWENDROFF_G", S_IRWXU | S_IRWXG | S_IRWXO);
                else mkdir("LAXWENDROFF_S", S_IRWXU | S_IRWXG | S_IRWXO);
                break;
        case 3:
                printf("\t-> Integration method: Lax-Friedrichs\n");
                if(k == 1) mkdir("LAXFRIEDRICHS_G", S_IRWXU | S_IRWXG | S_IRWXO);
                else mkdir("LAXFRIEDRICHS_S", S_IRWXU | S_IRWXG | S_IRWXO);
                break;
    }

}

int main() {

    int k, method;

    //----------------definisco gli elementi nella struct
    aveq A;
    f_params P;
    function S, G;

    G.f = gaussian;
    G.params = &P;
    
    S.f = square_wave;
    S.params = &P;

    P.x0 = 5.;
    P.x1 = 4.;
    P.x2 = 6.;
    //----------------



    k = 0; // 1 Gaussiana, 0 onda quadra
    method = 2; //metodo di integrazione (cfr. simulation)

    printf("\n");
    printf("\t----------\n");
    printf("\tSIMULATION\n");
    printf("\t----------\n");

    if(k == 1) {
        printf("\t-> Initial conditions: gaussian wave\n");
    } else {
        printf("\t-> Initial conditions: square wave\n");
    }

        create_directory(method, k);
        printf("\n");

        printf("\t-> Starting simulation for different parameters:\n");
        simulation(&A, &G, &S, k, method);
        //printf("\t-> Calculating L2norm\n");


    printf("\n");
    return 0;
}
