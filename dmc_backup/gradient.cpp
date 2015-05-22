/********************************************************************************/
                //GRADIENT AND LAPLACIAN

void gradient_psi_test(configuration *params , phase_space *walker ,int i, int j, double *derive){
    double fl, fr, dr, r ,x,y,z;
    int k;

    dr  =params->delta_interval/5.0;
    x = (walker+i)->position_x[j];
    y = (walker+i)->position_y[j];
    z = (walker+i)->position_z[j];


    (walker+i)->position_x[j] = x -dr;
    fl =psi_test( params, walker, i  );
    (walker+i)->position_x[j] = x +dr;
    fr =psi_test( params, walker, i  );
    (walker+i)->position_x[j] =x;
    derive[0] = (0.5*(fr-fl)/dr);


    (walker+i)->position_y[j] = y -dr;
    fl =psi_test( params, walker, i  );
    (walker+i)->position_y[j] = y +dr;
    fr =psi_test( params, walker, i  );
    (walker+i)->position_y[j] =y;
    derive[1] = (0.5*(fr-fl)/dr);


    (walker+i)->position_z[j] = z -dr;
    fl =psi_test( params, walker, i  );
    (walker+i)->position_z[j] = z +dr;
    fr =psi_test( params, walker, i  );
    (walker+i)->position_z[j] =z;
    derive[2] = (0.5*(fr-fl)/dr);

    //printf("grad (%g , %g , %g )\n" , derive[0] , derive[1] , derive[2]);
}


double laplacian_psi_test(configuration *params , phase_space *walker , int i, int j){
    double fl[3], fr[3],f[3], ff[3] , dr ,x,y,z ;
    int k;

    dr = params->delta_interval/5.0; //delta r easy
    x = (walker+i)->position_x[j]; //walker[i]->position_x[j]
    y = (walker+i)->position_y[j];
    z = (walker+i)->position_z[j];


    (walker+i)->position_x[j] = x - dr;
    fl[0] = psi_test( params, walker, i  );
    (walker+i)->position_x[j] = x + dr; //questo è GIUSTO
    fr[0] =psi_test( params, walker, i  );
    (walker+i)->position_x[j] =x;
    f[0] = psi_test(params, walker , i);
    ff[0] = (1.0/(dr*dr))*(fr[0] + fl[0] -2.0*f[0]); //derivata precisa
//printf("\t\t%g : fr %g fl %g -2f %g /drdr= %g\n",x ,fr[0] ,fl[0] , f[0], ff[0] );

    (walker+i)->position_y[j] = y -dr;
    fl[1] =psi_test( params, walker, i  );
    (walker+i)->position_y[j] = y +dr;
    fr[1] =psi_test( params, walker, i  );
    (walker+i)->position_y[j] =y;
    f[1] = psi_test(params, walker , i);
    ff[1] = (1.0/(dr*dr))*(fr[1] + fl[1] -2.0*f[1]);
//printf("\t\t%g: fr %g fl %g -2f %g /drdr= %g\n",y ,fr[1] ,fl[1] , f[1], ff[1] );
    (walker+i)->position_z[j] = z -dr;
    fl[2] =psi_test( params, walker, i  );
    (walker+i)->position_z[j] = z +dr;
    fr[2] =psi_test( params, walker, i  );
    (walker+i)->position_z[j] =z;
    f[2] = psi_test(params, walker , i);
    ff[2] = (1.0/(dr*dr))*(fr[2] + fl[2] -2.0*f[2]);

    //printf("\t\t%g: fr %g fl %g -2f %g /drdr= %g\n\n",z, fr[2] ,fl[2] , f[2], ff[2] );
    return (ff[0] + ff[1] + ff[2]);
    }


/********************************************************************************/
