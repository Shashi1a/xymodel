# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include "mt19937ar.h"
# include "const.h"


int main(){

    int i,j ;
    long unsigned int seed[1000];
    unsigned long init[4] = {0x123, 0x234, 0x345, 0x456}, length = 4;
    char simtype[1024];
    int * right, * left, *up, *down;
    FILE * OUT ;

    double *enr_arr ;
    double *mag_arr;
    double *theta;
    double *mag_arr_x ;
    double *mag_arr_y;

    n_sites = L * L;

    // pi = 4*arctan(1);
    pi_val = 4.0 * atan(1.0);

    // declaring arrays for the neighbour table
    right = (int *)malloc(n_sites * sizeof(double));
    left = (int *)malloc(n_sites * sizeof(double));
    up = (int *)malloc(n_sites * sizeof(double));
    down = (int *)malloc(n_sites * sizeof(double));

    // declaring  array for the spin
    theta = (double *)malloc(n_sites * sizeof(double));


    // initializing the random number generator
    init_by_array(init, length);

    // calling the subroutine to set up the neighbour table
    neighbour(right, left, up, down);

    // this is for normal quench. 
    // 1. Equilibrate system at some high temperature
    // 2. Lower the temperature and make measurement
    if(quench_bool==1){
        sprintf(simtype,"quench");
        // random initial configuration

        latticeint(theta, right, left, up, down);

        // setting size of arrays based on the quench type
        enr_arr = (double *)malloc(meas_quench * sizeof(double));
        mag_arr = (double *)malloc(meas_quench * sizeof(double));
        mag_arr_x = (double *)malloc(meas_quench * sizeof(double));
        mag_arr_y = (double *)malloc(meas_quench * sizeof(double));

        // setting up initial temperature
        T_0 = 2.0 * T_kt ;
        beta = 1./T_0;
        for(i=0;i<equil_quench;i++){
            sweep(theta, right, left, up, down);
            enrmag_measurement(theta, right, left, up, down);
            //print_conf(theta, i,simtype);
        }
        
        // setting up the final temperature
        T_fin  = T_kt ;
        beta = 1. / T_fin;
        for(i=0;i<meas_quench;i++){
            sweep(theta, right, left, up, down);
            enrmag_measurement(theta, right, left, up, down);
            enr_arr[i] = enr;
            mag_arr_x[i] = mag_x;
            mag_arr_y[i] = mag_y;
            mag_arr[i] = sqrt((mag_x * mag_x) + (mag_y * mag_y));
            print_conf(theta, i,simtype);}
    }

    // we use an random configuration as  iniital state
    // set the temperature and start making measurements
    if(sdn_quench==1){
        sprintf(simtype,"sudden_quench");
        // initializing infinite temperature state
        latticeint(theta, right, left, up, down);

        // number of measurements equil to meas_sdn_quench
        enr_arr = (double *)malloc(meas_sdn_quench * sizeof(double));
        mag_arr = (double *)malloc(meas_sdn_quench * sizeof(double));
        mag_arr_x = (double *)malloc(meas_sdn_quench * sizeof(double));
        mag_arr_y = (double *)malloc(meas_sdn_quench * sizeof(double));

        // setting temperature to T_c and making measurements
        T_fin = 0.5 ;
        beta = 1. / T_fin;
        for (i = 0; i < meas_sdn_quench; i++){
            sweep(theta, right, left, up, down);
            //enrmag_measurement(theta, right, left, up, down);
            //enr_arr[i] = enr;
            //mag_arr_x[i] = mag_x;
            //mag_arr_y[i] = mag_y;
            //mag_arr[i] = sqrt((mag_x * mag_x) + (mag_y * mag_y));
            print_conf(theta, i,simtype);}
    }
    
    // use random configuration as an initial configuration
    // reduce temperature in a linear fashion wrt a rate
    if(lnr_quench == 1){
        sprintf(simtype,"linearquench");
        // initialize random initial state        
        latticeint(theta, right, left, up, down);

        // number of measurement equal to meas_lnr_quench
        enr_arr = (double *)malloc(meas_lnr_quench * sizeof(double));
        mag_arr = (double *)malloc(meas_lnr_quench * sizeof(double));
        mag_arr_x = (double *)malloc(meas_lnr_quench * sizeof(double));
        mag_arr_y = (double *)malloc(meas_lnr_quench * sizeof(double));

        // slowly initial temperature 
        T_0 = 2*T_kt;
        beta = 1./T_0 ;

        // relaxing the system at temperature T_0 before cooling it slowly
        for (i=0;i<equil_lnr_quench;i++){
            sweep(theta,right,left,up,down);
        }
        // slow cooling starting at temperature T_0
        for (i = 0; i < (2*rate_t)+1; i++){
            // t_i--> T_c (i = -rate_t)
            // t_i--> 0 (i=rate_t)
            t_i = T_kt*(1.0-(double)(i-rate_t)/rate_t);      
            beta = 1./t_i ;  
            sweep(theta, right, left, up, down);
            enrmag_measurement(theta, right, left, up, down);
            enr_arr[i] = enr;
            mag_arr_x[i] = mag_x;
            mag_arr_y[i] = mag_y;
            mag_arr[i] = sqrt((mag_x * mag_x) + (mag_y * mag_y));
            print_conf(theta,t_i,simtype);
        }
    }

    if (no_quench==1){
        sprintf(simtype,"noquench");
        // initializing the lattice, the theta are in range [-pi,pi]
        latticeint(theta,right,left,up,down);

        // number of measurements taken during equilibration cycle
        enr_arr = (double *)malloc(T_equil * sizeof(double));
        mag_arr = (double *)malloc(T_equil * sizeof(double));
        mag_arr_x = (double *)malloc(T_equil * sizeof(double));
        mag_arr_y = (double *)malloc(T_equil * sizeof(double));

        // call the monte-carlo sweep normal no quench
        // comment the below loop to only perform quenching
        beta = 1./Temp;  
        for(i = 0; i <T_equil;i++){
            sweep(theta,right,left,up,down);
            enrmag_measurement(theta,right,left,up,down);
            enr_arr[i] = enr ;
            mag_arr_x[i] = mag_x ;
            mag_arr_y[i] = mag_y ;
            mag_arr[i] = sqrt((mag_x*mag_x) + (mag_y*mag_y));
            //print_conf(theta, i,simtype);
            }

        // calling the subroutine to print the equilibration data
        print_eq(theta, enr_arr, mag_arr_x, mag_arr_y, mag_arr,simtype);

        // perform measurement
        for(i=0;i<n_meas;i++){
            sweep(theta,right,left,up,down);
            enrmag_measurement(theta,right,left,up,down);
            print_conf(theta,i,simtype);}
        }


    return 0;
}

/////////////////////////////////////////////////////////////////
/////////////// neighbour table for the system///////////////////
/////////////////////////////////////////////////////////////////

void neighbour(int *right, int *left, int *up, int *down){
        int site;
		int x,y;
		int rn,ln,un,dn;
		for(site = 0; site<n_sites;site++){
		    y = (int)((site) / L);
		    x = (site) % L;
		    rn = (x + 1 )% L;	
            ln = (x - 1 + L)% L;
		    un = (y + 1 )% L;	
            dn = (y - 1 + L)% L;
		    right[site]  = y * L + rn; 
            left[site]   = y * L + ln;
		    up[site]     = un * L + x;
            down[site]   = dn * L + x;}
        

}

/////////////////////////////////////////////////////////////////
/////////////// initializing the lattice randomly////////////////
/////////////////////////////////////////////////////////////////

void latticeint(double *theta,int*right,int*left,int*up,int*down){
    int i ;
    double rand;
    double randpi ; 
    printf("pi_val: %f \n",pi_val);
    
    //FILE * OUT ;
    //OUT = fopen("pi.dat","w");
    for (i = 0; i < n_sites;i++){
        rand = genrand_real1();
        randpi = (2.0*pi_val*rand)-pi_val;
        theta[i] = randpi;

        //fprintf(OUT,"%f \n",spin[i]);
        
    }
    //fclose (OUT);
}

///////////////////////////////////////////////////////////////////////
////////////////////////monte carlo sweep ////////////////////////////
///////////////////////////////////////////////////////////////////////

void sweep(double *theta, int *right, int *left, int *up, int *down){
    int i ;
    int si,ri,li,ui,di ;
    double f_ri ;

    double termir , termil , termiu , termid ;

    double e_u , e_v ;
    double delE ;
    double rand_n ;

    // the array to store the random sites that are to be flipped;
    // the array that store the value of random field;
    int * rand_site;
    double * rand_dir ; 
    FILE * OUT ;

    rand_site = (int*)malloc(n_sites*sizeof(int));
    rand_dir = (double*)malloc(n_sites*sizeof(double));
    //OUT = fopen("PI.dat","a+");
    // setting the arrays to store the random sites and random fields values.
    for (i=0;i<n_sites;i++){
        rand_site[i] = (int)(genrand_real2() * n_sites) ;
        rand_dir[i] = (2.0*pi_val*genrand_real1())-pi_val ;
    //    fprintf(OUT,"%lf\n",rand_dir[i]);
    }
    //fclose(OUT);
    // monte carlo to flip the variable at each site in the range [-pi,pi]
    for(i=0;i<n_sites;i++){
        si = rand_site[i];
        f_ri = rand_dir[i];
        
        // getting values of neighbouring sites
        ri = right[si];
        li = left[si];
        ui = up[si];
        di = down[si];

        // energy contribution of each neighbour
        termir = cos(theta[si] - theta[ri]);
        termil = cos(theta[si] - theta[li]);
        termiu = cos(theta[si] - theta[ui]);
        termid = cos(theta[si] - theta[di]);

        // energy of initial conf
        e_u = -J * (termir + termil + termiu + termid);
        
        // energy contribtion due to new direction
        termir = cos(f_ri - theta[ri]);
        termil = cos(f_ri - theta[li]);
        termiu = cos(f_ri - theta[ui]);
        termid = cos(f_ri - theta[di]);
        
        // energy of final conf
        e_v = - J * (termir + termil + termiu + termid);

        // energy difference 
        delE = e_v - e_u ;

        // change the value if it reduces the energy
        if (delE < 0.0){
            theta[si] = f_ri;
        }
        else{
            rand_n = genrand_real2();
            if (rand_n<exp(-beta*delE)){
                theta[si] = f_ri;
            }
            

        }

    }
}

////////////////////////////////////////////////////////
//////////////////energy calculation////////////////////
////////////////////////////////////////////////////////

void enrmag_measurement(double *theta,int*right,int*left,int*up,int*down){
    int i ;
    double si,li,ri,ui,di;
    double tempir, tempil, tempiu, tempid;
    enr = 0.0 ;
    mag_x = 0.0 ;
    mag_y = 0.0 ;
    for(i = 0;i<n_sites;i++){
        si = theta[i];
        ri = theta[right[i]];
        li = theta[left[i]];
        ui = theta[up[i]];
        di = theta[down[i]];

        tempir = cos(si-ri);
        tempil = cos(si-li);
        tempiu = cos(si-ui);
        tempid = cos(si-di);
        
        enr = enr -J *0.5* ( tempir + tempil + tempiu + tempid );
        mag_y = mag_y + sin(theta[i]);
        mag_x = mag_x + cos(theta[i]);
    }
    enr = enr/n_sites; 
    mag_x = mag_x/n_sites;
    mag_y = mag_y/n_sites ;
}

//////////////////////////////////////////////////////
////////////// modulus of a number ////////////////////
//////////////////////////////////////////////////////

double calc_modulus(double val_t){
    double ans_v ; 
    if  (val_t<0.0){
        ans_v = (-1) * fmod(fabs(val_t),pi_val);
        return ans_v;
    }
    else{
        ans_v = fmod(val_t,pi_val);
        return ans_v;
    }
}
///////////////////////////////////////////////////////
///////////////// to print equilibration data//////////
///////////////////////////////////////////////////////


void print_eq(double * theta,double *enr_arr,double*mag_arr_x,
        double*mag_arr_y,double*mag_arr,char simtype[])
    { 
        int i ;
        FILE * OUT ;
        char streq[1024] ;

        if (no_quench==1){
        sprintf(streq,"%s_Equildata_L%d_Temp%lf.dat",simtype,L,Temp);}
        else if(no_quench==0){
            sprintf(streq, "%s_Equildata_L%d_Temp%lf.dat",simtype,L, Temp);}
        OUT = fopen(streq, "w");

        for (i = 0; i < T_equil; i++)
        {
            fprintf(OUT, "%d \t %.15f \t %.15f \t %.15f \t %.15f \n", i, 
            enr_arr[i], mag_arr_x[i], mag_arr_y[i], mag_arr[i]);
        }
        fclose(OUT);
    }

////////////////////////////////////////////////////////////////
//////////////////print configurations//////////////////////////
////////////////////////////////////////////////////////////////


void print_conf(double *theta,double t_i,char simtype[]){

    int i;
    double j ; 
    FILE * OUT ;
    char strconf[1024];

    j = t_i ;
    if (no_quench==1){
        sprintf(strconf,"%s_sweepspinconf_L%d_tfinal%lf.dat",simtype,L,Temp);}
    else if(no_quench==0){
        sprintf(strconf, "%s_sweepspinconf_L%d.dat", simtype, L);
    }

    OUT = fopen(strconf,"a");
    fprintf(OUT,"%lf\t",j);
    for (i = 0 ;i <n_sites;i++){
        fprintf(OUT,"%lf\t",theta[i]);
    }
    fprintf(OUT, "\n");
    fclose(OUT);
}