/******************************************************************************************/
/******************************************************************************************/
int  L = 24 ; // linear dimension
int  n_sites; // number of sites

// parameters for no quench simulations
int no_quench = 1;     // set if you don't want any quench
int  T_equil = 100000; // time steps for the equilibration of the system
int n_meas = 10000 ; // no of measurement for no quench

// parameters for the normal quench simulations
int quench_bool = 0; // set this variable to perform normal quench T0 (2T_c)-->T_c
int equil_quench = 10; // equilibration cycle for the normal quench 
int meas_quench = 10 ;  // measurement cycle for the normal quench

// parameters for the sudden quench simulations
int sdn_quench = 0;            // set this variable to perform quench from T-->infinity to T_c
int meas_sdn_quench = 3000 ; // measurement cycle for the sudden quench

// parameters for the linear quench simulations
int lnr_quench = 0;         // set for linear quench
int equil_lnr_quench = 10 ; // equilibration cycle for the linear quench
int meas_lnr_quench = 10 ; // measurement cycle for the linear quench
int rate_t = 1024; // rate for the linear quench simulations

double J = 1.0 ;
double Temp = 0.90 ;
double T_kt = 0.89;
double T_fin ;
double T_0  ; 
double t_i ;
double beta;

double enr;
double mag_x;
double mag_y;
double pi_val  ; 


// subroutines for the simulations
void neighbour(int*right,int*left,int*up,int*down);
void latticeint(double *theta,int *right, int *left,int *up,int *down);
void sweep(double *theta,int *right,int *left,int *up,int *down);
void enrmag_measurement(double *theta, int *right, int *left, int *up, int *down);
double calc_modulus(double val_t);
void print_eq(double *theta, double *enr_arr, double *mag_arr_x,
              double *mag_arr_y, double *mag_arr,char simtype[]);
void print_conf(double *theta,double t_i, char simtype[]);
/****************************************************************************************/
/****************************************************************************************/
