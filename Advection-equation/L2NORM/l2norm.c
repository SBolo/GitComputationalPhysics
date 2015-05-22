#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define GSL(F,x) (*((F)->function))(x,(F)->params)
#define XSTEP (1.e-2)
#define CF (0.5)
#define TSTEP (CF * XSTEP)
#define XMAX 10.
#define DIM (int)(XMAX/XSTEP+1)


typedef struct {
	double (*function)(double x, void *params);
	void *params;
} gsl_function;


double initial_condition(double x, void *params) {
    double *p = (double *)params;
    double x0 = p[0];
    
    return exp( - (x - x0)*(x - x0) ) ;
}

double square_wave(double x, void *params) {
    double *p = (double *)params;
    double x1 = p[0];
    double x2 = p[1];
    
    if( x > x1 && x < x2 ) {
        return 1.;
    } else {
        return 0.;
    }
    
}

void FTCS(gsl_function *IC, double *solution, double time) {
    
    int i, t, T = (int)(time/TSTEP);
    double tsol[DIM];
    
    
    for(i = 0; i < DIM; i++) tsol[i] = solution[i] = GSL(IC,i*XSTEP);

    
    for(t = 0; t <= T; t++ ) {
    
                for( i = 1; i < DIM-1; i++ ) {
                    solution[i] = tsol[i] - (CF/2.) * (tsol[i+1] - tsol[i-1]);
                }
        
            solution[0] = tsol[0] - (CF/2.) * (tsol[1] - tsol[DIM-1]);
            solution[DIM-1] = tsol[DIM-1] - (CF/2.) * (tsol[0] - tsol[DIM-2]);
        
                for ( i = 0; i < DIM; i++ ) {
                    tsol[i] = solution[i];
                }
            
        }
    
}

void LeFr(gsl_function *IC, double *solution, double time) {
    
    int i, t, T = (int)(time/TSTEP);
    double tsol[DIM], t_1sol[DIM]; //t1sol è al tempo t, t_2sol al tempo t-1
    
    
    for(i = 0; i < DIM; i++) t_1sol[i] = tsol[i] = solution[i] = GSL(IC,i*XSTEP);

    
    for( t = 0; t <= T; t++ ) {
        
        if( t == 0) { } else if (t == 1) { } else {
            
            for( i = 1; i < DIM-1; i++ ) {
                solution[i] = t_1sol[i] - CF * (tsol[i+1] - tsol[i-1]);
            }
            
            solution[0] = t_1sol[0] - CF * (tsol[1] - tsol[DIM-1]);
            solution[DIM-1] = t_1sol[DIM-1] - CF * (tsol[0] - tsol[DIM-2]);
            
            for ( i = 0; i < DIM; i++ ) {
                t_1sol[i] = tsol[i];
                tsol[i] = solution[i];
            }
        }
    }
    
    
}

void LW(gsl_function *IC, double *solution, double time) {
    
    int i, t, T = (int)(time/TSTEP);
    double tsol[DIM];
    
    
    for(i = 0; i < DIM; i++) {
        tsol[i] = GSL(IC,i*XSTEP);
        solution[i] = tsol[i];
    }
    
        for( t = 0; t <= T; t++ ) {
        
                for( i = 1; i < DIM-1; i++ ) {
                    solution[i] = tsol[i] - (CF/2.) * (tsol[i+1] - tsol[i-1]) + CF * (CF/2.) * (tsol[i+1] - 2.*tsol[i] + tsol[i-1]);
                }
        
            solution[0] = tsol[0] - (CF/2.) * (tsol[1] - tsol[DIM-1]) + CF * (CF/2.) * (tsol[1] - 2.*tsol[0] + tsol[DIM-1]);
            solution[DIM-1] = tsol[DIM-1] - (CF/2.) * (tsol[0] - tsol[DIM-2]) + CF * (CF/2.) * ( tsol[0] - 2.* tsol[DIM-1] + tsol[DIM-2]);
        
                for ( i = 0; i < DIM; i++ ) {
                    tsol[i] = solution[i];
                }
        }

}


void LF(gsl_function *IC, double *solution, double time) {
    
    int i, t = 0, T = (int)(time/TSTEP);
    double tsol[DIM];
    
    
    for(i = 0; i < DIM; i++) tsol[i] = GSL(IC,i*XSTEP);
    
    
        do {
        
            if( t == 0 ) solution[i] = tsol[i];
            
            else {
        
                    for( i = 1; i < DIM-1; i++ ) {
                        solution[i] = 1./2. * (tsol[i-1] + tsol[i+1]) - (CF/2.) * (tsol[i+1] - tsol[i-1]);
                    }
        
                solution[0] = 1./2. * (tsol[DIM-1] + tsol[1]) - (CF/2.) * (tsol[1] - tsol[DIM-1]);
                solution[DIM-1] = 1./2. * (tsol[DIM-2] + tsol[0]) - (CF/2.) * (tsol[0] - tsol[DIM-2]);
        
                    for ( i = 0; i < DIM; i++ ) {
                        tsol[i] = solution[i];
                    }
            }
            
            t++;
            
        } while( t <= T );
    
}


void LF_wave(gsl_function *IC, double *solution, double time) {
    
    double r[DIM], s[DIM], tmpr[DIM], tmps[DIM], tsol[DIM];
    int i, t = 0, T = (int)(time/TSTEP);
    
    for(i = 0; i < DIM; i++ ) {
        tmps[i] = 0.;
        tmpr[i] = 2. * (i*XSTEP - 5.) * GSL(IC,i*XSTEP);
        
        tsol[i] = GSL(IC,i*XSTEP);
    }
    
    do {
        
        for( i = 1; i < DIM-1; i++ ) {
            r[i] = 1./2. * (tmpr[i+1] + tmpr[i-1]) - CF/2. * ( tmps[i+1] - tmps[i-1] );
            s[i] = 1./2. * (tmps[i+1] + tmps[i-1]) - CF/2. * ( tmpr[i+1] - tmpr[i-1] );
            solution[i] = tsol[i] + TSTEP/2. * ( s[i] + tmps[i] );
        }
                
        r[0] = 1./2. * (tmpr[1] + tmpr[DIM-1]) - CF/2. * ( tmps[1] - tmps[DIM-1] );
        r[DIM-1] = 1./2. * (tmpr[0] + tmpr[DIM-2]) - CF/2. * ( tmps[0] - tmps[DIM-2] );
        
        s[0] = 1./2. * (tmps[1] + tmps[DIM-1]) - CF/2. * ( tmpr[1] - tmpr[DIM-1] );
        s[DIM-1] = 1./2. * (tmps[0] + tmps[DIM-2]) - CF/2. * ( tmpr[0] - tmpr[DIM-2] );
        
        solution[0] = tsol[0] + TSTEP/2. * ( s[0] + tmps[0] );
        solution[DIM-1] = tsol[DIM-1] + TSTEP/2. * ( s[DIM-1] + tmps[DIM-1] );
        
        for ( i = 0; i < DIM; i++ ) {
            tsol[i] = solution[i];
            tmps[i] = s[i];
            tmpr[i] = r[i];
        }
        
        t++;
        
    } while( t <= T );
    
    
}

void LW_wave(gsl_function *IC, double *solution, double time) {
    
    double r[DIM], s[DIM], tmpr[DIM], tmps[DIM], tsol[DIM];
    int i, t = 0, T = (int)(time/TSTEP);
    
    for(i = 0; i < DIM; i++ ) {
        tmps[i] = 0.;
        tmpr[i] = 2. * (i*XSTEP - 5.) * GSL(IC,i*XSTEP);
       // tmpr[i] = ( GSL(IC, i*XSTEP + XSTEP) - GSL(IC,i*XSTEP) ) / XSTEP;
        
        tsol[i] = GSL(IC,i*XSTEP);
    }
    
    do {
        
        for( i = 1; i < DIM-1; i++ ) {
            r[i] = tmpr[i] - CF/2. * ( tmps[i+1] - tmps[i-1] ) + CF * CF/2. * ( tmpr[i+1] - 2. * tmpr[i] + tmpr[i-1] );
            s[i] = tmps[i] - CF/2. * ( tmpr[i+1] - tmpr[i-1] ) + CF * CF/2. * ( tmps[i+1] - 2. * tmps[i] + tmps[i-1] );
            solution[i] = tsol[i] + TSTEP/2. * ( s[i] + tmps[i] );
        }
        
        r[0] = tmpr[0] - CF/2. * ( tmps[1] - tmps[DIM-1] ) + CF * CF/2. * ( tmpr[1] - 2. * tmpr[0] + tmpr[DIM-1] );
        r[DIM-1] = tmpr[DIM-1] - CF/2. * ( tmps[0] - tmps[DIM-2] ) + CF * CF/2. * ( tmpr[0] - 2. * tmpr[DIM-1] + tmpr[DIM-2] );
        
        s[0] = tmps[0] - CF/2. * ( tmpr[1] - tmpr[DIM-1] ) + CF * CF/2. * ( tmps[1] - 2. * tmps[0] + tmps[DIM-1] );
        s[DIM-1] = tmps[DIM-1] - CF/2. * ( tmpr[0] - tmpr[DIM-2] ) + CF * CF/2. * ( tmps[0] - 2. * tmps[DIM-1] + tmps[DIM-2] );
        
        solution[0] = tsol[0] + TSTEP/2. * ( s[0] + tmps[0] );
        solution[DIM-1] = tsol[DIM-1] + TSTEP/2. * ( s[DIM-1] + tmps[DIM-1] );
        
        for ( i = 0; i < DIM; i++ ) {
            tsol[i] = solution[i];
            tmps[i] = s[i];
            tmpr[i] = r[i];
        }
        
        t++;
        
    } while( t <= T );
    
    
}


void L2norm(gsl_function *IC, double *solution, double time) {
    
    double norm = 0.;
    int T = (int)(time/TSTEP), t, i;
    
    for(t = 1; t <= T; t++) {
        
        //FTCS(IC,solution,t*TSTEP);
        //LF(IC, solution, t*TSTEP);
        //LW(IC,solution,t*TSTEP);
        //LeFr(IC,solution,t*TSTEP);
        //LF_wave(IC,solution,t*TSTEP);
        LW_wave(IC,solution,t*TSTEP);
        
        for( i = 0; i < DIM; i++ ) {
            norm += solution[i] * solution[i];
        }
        
            printf("%g \t %g\n", (t*TSTEP), sqrt(1./(DIM)*norm));
         norm = 0;
    }
    
}


int main() {
    
    double solution[DIM], params[1], params2[2];
    int i;
    
    gsl_function IC, SW;
    IC.function = initial_condition;
    IC.params = params;
    
    SW.function = square_wave;
    SW.params = params2;
    
    params[0] = 5.;
    
    params2[0] = 4.;
    params2[1] = 6.;
        
    
   L2norm(&IC, solution, 20);

    
    return 0;
}


