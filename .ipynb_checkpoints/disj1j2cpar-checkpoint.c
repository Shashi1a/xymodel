#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "const.h"
#include "mt19937ar.h"

int
main(int argc , char *argv[])
	{

	int i,r,l,u,d,p,j,k,m,site;
	int size;
	int dis , box ;
	int dis_min ;
	int dis_max ;
	int 	*spin;
	int 	*sitepr;
	int 	*right, *left ,*up , *down;
	int 	*right_up, *right_down, *left_up, *left_down;
	double Tmax , Tmin , delT;
	double 	*Jleft,*Jright,*Jup,*Jdown;
	double *******P_add;
	long unsigned int seed[1000];
	unsigned long init[4]={0x123, 0x234, 0x345, 0x456}, length=4;


	double *mag_fe, *magfe_box, *magfe_box2, *magfe_box4;
	double *mag_st, *magst_box, *magst_box2, *magst_box4;
	double *magstn_box, *magstn_box2, *magstn_box4;
	double *energy, *energy_box, *energy_box2, *energy_box4;

	char string2[1024];
	FILE *result2;
	FILE * INP ;
	N_sites = L * L ;
	opti = 1;

/*************************generator stuff******************************************/


   	init_by_array(init, length);
   	for(i=0;i<n_dis;i++){
  	seed[i] = genrand_int32();/*printf("seed:%lu\n",seed[dis]);*/}



/***********************for neighbour table*****************************************/
	spin   = (int*)malloc(N_sites*sizeof(int));
	sitepr = (int*)malloc(N_sites*sizeof(int));

	right 	   = (int*)malloc(N_sites*sizeof(int));	left = (int*)malloc(N_sites*sizeof(int));
	up   	   = (int*)malloc(N_sites*sizeof(int));	down = (int*)malloc(N_sites*sizeof(int));

	right_up   = (int*)malloc(N_sites*sizeof(int));	right_down = (int*)malloc(N_sites*sizeof(int));
	left_up    = (int*)malloc(N_sites*sizeof(int));	left_down  = (int*)malloc(N_sites*sizeof(int));

	Jright  =   (double*)malloc(N_sites*sizeof(double));	Jleft = (double*)malloc(N_sites*sizeof(double));
	Jup     = (double*)malloc(N_sites*sizeof(double));  	Jdown  = (double*)malloc(N_sites*sizeof(double));
	P_add = (double*******)malloc(N_sites*sizeof(double******));

/*********************for probability table*****************************************/
 	for(i=0;i<N_sites;i++){
	P_add[i] = (double******)malloc(3*sizeof(double*****));
	for(j=0;j<3;j++){
	P_add[i][j] = (double*****)malloc(3*sizeof(double****));
	for(k =0;k<3;k++){
	P_add[i][j][k] = (double****)malloc(3*sizeof(double***));
	for(p=0;p<3;p++){
	P_add[i][j][k][p] = (double***)malloc(3*sizeof(double**));
	for(m=0;m<3;m++){
    P_add[i][j][k][p][m]=(double**)malloc(3*sizeof(double*));
   	for(l =0;l<3;l++){
   	P_add[i][j][k][p][m][l] = (double*)malloc(10*sizeof(double));}}}}}}

/**********************to store obsevables generated during the monte carlo step********/

	mag_fe = (double*)malloc(box_size*sizeof(double));
	mag_st = (double*)malloc(box_size*sizeof(double));
	energy = (double*)malloc(box_size*sizeof(double));

	magfe_box = (double*)malloc(N_box*sizeof(double));magfe_box2 = (double*)malloc(N_box*sizeof(double));magfe_box4 = (double*)malloc(N_box*sizeof(double));
	magst_box = (double*)malloc(N_box*sizeof(double));magst_box2 = (double*)malloc(N_box*sizeof(double));magst_box4 = (double*)malloc(N_box*sizeof(double));
	magstn_box = (double*)malloc(N_box*sizeof(double));magstn_box2 = (double*)malloc(N_box*sizeof(double));magstn_box4 = (double*)malloc(N_box*sizeof(double));
	energy_box = (double*)malloc(N_box*sizeof(double));energy_box2 = (double*)malloc(N_box*sizeof(double));energy_box4 = (double*)malloc(N_box*sizeof(double));

/****************************neighbour table******************************/
	neighbour(right,left,up,down,right_up,right_down,left_up,left_down);

/*	double *c_v, *c_v2; 		c_v = (double*)malloc(N_box*sizeof(double)); c_v2 = (double*)malloc(N_box*sizeof(double));
	double *chi_st, *chi_st2;	chi_st = (double*)malloc(N_box*sizeof(double)); chi_st2 = (double*)malloc(N_box*sizeof(double));
	double *chi_fe, *chi_fe2; 	chi_fe = (double*)malloc(N_box*sizeof(double)); chi_fe2 = (double*)malloc(N_box*sizeof(double));
	double *st_cum, *st_cum2; 	st_cum = (double*)malloc(N_box*sizeof(double));  st_cum2 = (double*)malloc(N_box*sizeof(double));
	double *fe_cum, *fe_cum2; 	fe_cum = (double*)malloc(N_box*sizeof(double)); fe_cum2 = (double*)malloc(N_box*sizeof(double));
	double *ener_cum, *ener_cum2; 	ener_cum = (double*)malloc(N_box*sizeof(double)); ener_cum2 = (double*)malloc(N_box*sizeof(double));


	double *resulte; 	resulte=(double*)malloc(2*sizeof(double));
	double *resultfe; 	resultfe=(double*)malloc(2*sizeof(double));
	double *resultst;   	resultst=(double*)malloc(2*sizeof(double));
	double *resultcv; 	    resultcv=(double*)malloc(2*sizeof(double));
	double *resultcumen; 	resultcumen=(double*)malloc(2*sizeof(double));
	double *resultchist; 	resultchist=(double*)malloc(2*sizeof(double));
	double *resultchife; 	resultchife=(double*)malloc(2*sizeof(double));
	double *resultcumfe; 	resultcumfe=(double*)malloc(2*sizeof(double));
	double *resultcumst; 	resultcumst=(double*)malloc(2*sizeof(double));*/

//	INP=fopen(argv[1],"r");
//	INP=fopen("inp.dat","r");
//	fscanf(INP,"%d\t %d\n \t %lf \t %lf \t %lf\n",&dis_min,&dis_max,&Tmax,&Tmin,&delT);
//	printf("dmn:%d\t dmx:%d\t tmx:%lf \t tmn:%lf\t dT:%lf\n",dis_min,dis_max,Tmax,Tmin,delT);
//	fclose(INP);
/*************************disorder configuration*************************/
	Tmax = 1.0 ;
	dis = 0;
		init_genrand(seed[dis]);avsite=0;
 		latticeint(sitepr,spin,right,left,up,down,right_up,right_down,left_up,left_down,Jright,Jleft,Jup,Jdown);
		measurement(spin,right,left,up,down,right_up,right_down,left_up,left_down,Jright,Jleft,Jup,Jdown);

		Temp=Tmax ;
			printf("%f\n",Temp);

			beta = 1./Temp;
			for(site=0;site<N_sites;site++){
			for(i=0;i<3;i++){
		for(r=0;r<3;r++){
		for(l=0;l<3;l++){
		for(u=0;u<3;u++){
		for(d=0;d<3;d++){
		for(p=0;p<10;p++){
		P_add[site][i][r][l][u][d][p]=exp(-2.*beta*(double)(i-1)*(Jright[site]*(double)(r-1) + Jleft[site]*(double)(l-1)+Jup[site]*(double)(u-1)+Jdown[site]*(double)(d-1)) - 2.*beta*(double)(i-1)*J2*(p-4));}}}}}}}
		printf("1");
/********************************loop for equilibration*******************************/
       	for(j = 0;j<T_equil; j++){
        	sweep(spin,right,left,up,down,right_up,right_down,left_down,left_up,Jright,Jleft,Jup,Jdown,P_add);}
				measurement(spin,right,left,up,down,right_up,right_down,left_up,left_down,Jright,Jleft,Jup,Jdown);
				printf("en:%lf \t mgn:%lf \t mgs:%lf  \t mgsdi:%lf \n",Ener,Magn_ferro,magx*magx+magy*magy,(magx*magx+magy*magy)/satmag);

/*********************measurement loop for montecarlo****************************/
	   			for(size = 0; size<box_size ; size++){

							sweep(spin,right,left,up,down,right_up,right_down,left_down,left_up,Jright,Jleft,Jup,Jdown,P_add);
							sprintf(string2,"config_L%d_Temp%f_J2%f.dat",L,Temp,fabs(J2));
							result2 = fopen(string2,"a");
						 	for(i = 0; i <N_sites; i++){
								fprintf(result2, "%d\t ",spin[i]);}
								fprintf(result2,"\n");
						  fclose(result2);}


    free(spin);
   	free(sitepr);
	free(right) ; free(left) ; free(up) ; free(down);
	free(right_up); free(right_down); free(left_up); free(left_down);
	free(Jleft); free(Jright); free(Jup); free(Jdown);
	free(P_add);

	free(mag_fe); free(magfe_box);  free(magfe_box2); free(magfe_box4);
	free(mag_st); free(magst_box); free(magst_box2); free(magst_box4);
	free(magstn_box); free(magstn_box2); free(magstn_box4);
	free(energy); free(energy_box); free(energy_box2); free(energy_box4);

	return 0 ;

	}

/***************************************************************************************************************************************/
void
neighbour(int*right,int*left,int*up,int*down,int*right_up,int*right_down,int*left_up,int*left_down)
	{
		int site;
		int in,jn;
		int rn,ln,un,dn;
		printf("%d\t %d\n",L,N_sites);
		for(site = 0; site<N_sites;site++){
		in = (site) / L;
		jn = (site) % L;
		rn = (jn + 1 + L)% L;	ln = (jn - 1 + L)% L;
		un = (in + 1 + L)% L;	dn = (in - 1 + L)% L;
		right[site]  = in * L + rn; left[site]   = in * L + ln;
		up[site]     = un * L + jn;	down[site]   = dn * L + jn;
		right_up[site] 	 =  L*un + rn;  right_down[site] =  L*dn + rn;
		left_up[site]    =  L*un + ln;  left_down[site]  =  L*dn + ln;}

	}
/***************************************************************************************************************************************/


/********************************************************************************************************************************************/
void
measurement(int *spin,int *right,
	int *left,int *up,int *down,int *right_up,int *right_down,int *left_up,int *left_down,double *Jright,double *Jleft,double *Jup,double *Jdown)
	{
		int i,x,y,mx,my;
		int ri,li,ui,di,rui,rdi,lui,ldi,si;
		double sum1;
		double sum2;
		double delta,delta1,delta2;
		magx=0.0;magy=0.0;Magn_ferro=0.0;Ener=0.0;

		for(i=0;i<N_sites;i++){
		ri=right[i];li=left[i];ui=up[i];di=down[i];
		rui=right_up[i];rdi=right_down[i];lui=left_up[i];ldi=left_down[i];
       	si = spin[i];
    		sum1   = Jright[i]*(double)spin[ri]+ Jleft[i]*(double)spin[li] + Jup[i]*(double)spin[ui] + Jdown[i]*(double)spin[di];
	 		sum2   = (double)spin[rui] +  (double)spin[rdi] +  (double)spin[lui] +  (double)spin[ldi];
	     	delta1 = (double)si * sum1;
	    	delta2 = (double)si * sum2;
	    	delta  = - 0.5 * delta1 -   0.5 * J2  * delta2;
	    	Ener +=   delta/(int)avsite;
	    	x  = i % L;y  = i / L;
	    	if(x % 2 != 1){mx = 1;}else{mx = -1;}
	    	if(y % 2 != 1){my = 1;}else{my = -1;}
	    	magx += (double)(si * mx)/(int)N_sites;
	    	magy += (double)(si * my)/(int)N_sites;
	    	Magn_ferro +=(double)si/(int)avsite;}

	}

/*****************************************************************************************************************************************/
void
sweep(int *spin,int *right,int *left,int *up,int *down,int *right_up,int *right_down,int *left_down,int *left_up,
	double *Jright,double *Jleft,double *Jup,double *Jdown,double *******P_add)
	{
		int ik;
		int site,sum2;
		int ri,si,li,ui,di,rui,rdi,lui,ldi;
		double sum1,delta,delta1,delta2,ku;
	    for(ik = 0; ik<N_sites; ik++){
	   	site = (int)N_sites*genrand_real2(); /*printf("%d\n",site);*/
	   	si = spin[site];
	   	if(si !=0){
	   	ri = right[site]; li = left[site]; ui = up[site]; di = down[site];
	  	rui = right_up[site]; rdi = right_down[site]; lui = left_up[site]; ldi = left_down[site];

	        sum1 = (double)spin[ri]*Jright[site] + (double)spin[li]*Jleft[site] + (double)spin[ui]*Jup[site] + (double)spin[di]*Jdown[site];
	        sum2 = spin[rui] + spin[rdi] + spin[lui] +	spin[ldi];   /*printf("%d\t",sum2);*/
	        delta1 = (double)si* sum1; /* printf("%lf\n",delta1);*/
	        delta2 = (double)si* (double)sum2;	/*printf("%lf\n",delta2);*/
	        delta  = 2.0* delta1 + 2.0*J2*delta2; /*printf("%lf\n",delta);*/
	    	if(delta < 0.)
        	{spin[site] = -spin[site];}
        	else{
			ku=genrand_real2();/*P_add[site][si+1][spin[ri]+1][spin[li]+1][spin[ui]+1][spin[di]+1][sum2+4]); */
	        if(ku<P_add[site][si+1][spin[ri]+1][spin[li]+1][spin[ui]+1][spin[di]+1][sum2+4])
	        {spin[site] = -spin[site];}
			else
			{spin[site] = spin[site];}}}}
	}

/***************************************************************************************************************************/
void
average(int box,double *mag_fe,double *mag_st,double *energy,double *magfe_box,double *magst_box,double *magstn_box,
	double *energy_box,double *magfe_box2,double *magst_box2,double *magstn_box2,double *energy_box2,double *magfe_box4,double *magst_box4,double *magstn_box4,double *energy_box4)
	{
		int i;
		double magf, mags, magsn, ener;
		double mag2f, mag2s, mags2n, ener2;
		double mag4f, mag4s, mags4n, ener4;
		magf = 0.0; mags = 0.0; magsn = 0.0; ener = 0.0;
		mag2f = 0.0; mag2s = 0.0; mags2n = 0.0; ener2 = 0.0;
		mag4f = 0.0; mag4s = 0.0; mags4n = 0.0; ener4 = 0.0;
		for(i=0;i<box_size;i++){
		magf += mag_fe[i]; mag2f += mag_fe[i]*mag_fe[i]; mag4f += mag_fe[i]*mag_fe[i]*mag_fe[i]*mag_fe[i];
		mags += sqrt(mag_st[i]); mag2s += mag_st[i] ; mag4s += mag_st[i] * mag_st[i];
		magsn += sqrt(mag_st[i]/satmag); mags2n += mag_st[i]/satmag; mags4n += (mag_st[i]/satmag)*(mag_st[i]/satmag);
		ener += energy[i]; ener2 += energy[i] * energy[i] ; ener4 += energy[i] * energy[i]*energy[i] * energy[i];}
		magfe_box[box] = magf/(double)box_size; magfe_box2[box] = mag2f/(double)box_size; magfe_box4[box] = mag4f/(double)box_size;
		magst_box[box] = mags/(double)box_size; magst_box2[box] = mag2s/(double)box_size; magst_box4[box] = mag4s/(double)box_size;
		magstn_box[box] = magsn/(double)box_size; magstn_box2[box] = mags2n/(double)box_size; magstn_box4[box] = mags4n/(double)box_size;
		energy_box[box] = ener/(double)box_size;energy_box2[box] = ener2/(double)box_size;energy_box4[box] = ener4/(double)box_size;
	}

/*************************************************************************************************************************
void
specific(double *energy_box,double *energy_box2,double *c_v,double *c_v2)
	{
		int i;
		for(i=0;i<N_box;i++){
		c_v[i] = avsite*beta*beta*(energy_box2[i] - energy_box[i]*energy_box[i]);
		c_v2[i] = c_v[i] * c_v[i];}
	}

**************************************************************************************************************************
void
susceptibility_st(double *magst_box,double *magst_box2,double *chi_st,double *chi_st2)
	{
		int i;
		for(i=0;i<N_box;i++){
		chi_st[i] = satmag*beta*(magst_box2[i]-magst_box[i]*magst_box[i]);
		chi_st2[i] = chi_st[i]*chi_st[i];}
	}

**************************************************************************************************************************
void
susceptibility_fe(double *magfe_box,double *magfe_box2,double *chi_fe,double *chi_fe2)
	{
		int i;
		for(i=0;i<N_box;i++){
		chi_fe[i] = avsite*beta*(magfe_box2[i] - magfe_box[i]*magfe_box[i]);
		chi_fe2[i] = chi_fe[i]*chi_fe[i];}
	}

************************************************************************************************************************
void
cumulant_fe(double *magfe_box2,double *magfe_box4,double *fe_cum,double *fe_cum2)
	{
		int i;
		for(i=0;i<N_box;i++){
		fe_cum[i] = 1. - 0.3333*(magfe_box4[i]/(magfe_box2[i]*magfe_box2[i]));
		fe_cum2[i] = fe_cum[i] * fe_cum[i];}
	}


**********************************************************************************************************************
void
cumulant_st(double *magst_box2,double *magst_box4,double *st_cum,double *st_cum2)
	{
		int i;
		for(i=0;i<N_box;i++){
		st_cum[i] = 2.0*(1.- 0.5*(magst_box4[i]/(magst_box2[i]*magst_box2[i])));
		st_cum2[i] = st_cum[i] * st_cum[i];}
	}
**********************************************************************************************************************
void
cumulant_en(double *energy_box2, double *energy_box4,double *ener_cum,double *ener_cum2)
	{
		int i;
		for(i=0;i<N_box;i++){
		ener_cum[i] = 1. - (1./3.)*(energy_box4[i]/(energy_box2[i]*energy_box2[i]));
		ener_cum2[i] = ener_cum[i]*ener_cum[i];}


	}
*********************************************************************************************************************
void
average1(double *resulte ,double *energy_box)
	{
		int i;
		double avg1=0.0; double avg2=0.0;
		double err;
		for(i=0;i<N_box;i++){
		avg1 += energy_box[i];
		avg2 += energy_box[i]*energy_box[i];}
		avg1 = avg1/(double)N_box;
		avg2 = avg2/(double)N_box;
		err = sqrt((1./(double)(N_box - 1))*(avg2 - avg1*avg1));
		resulte[0] = avg1;
		resulte[1] = err;

	}
*/


/*************************initializing the lattice**************************************************************************/
void
latticeint(int *sitepr,int *spin,int *right,
	int *left,int *up,int *down,int *right_up,int *right_down,int *left_up,int *left_down,double *Jright,double *Jleft,double *Jup,double *Jdown)
	{
		int i,xi,yi,indx,idx;
		double E1,Mf1,Ms1;
		double E2,Mf2,Ms2;
		double k;
		avsite=0;

		for(i = 0;i<N_sites;i++){
	    if(genrand_real2()<0.50){Jright[i] =  J1 - imp;}
        else{Jright[i] = J1 + imp;}
        if(genrand_real2()<0.50){Jup[i] =  J1 - imp;}
        else{Jup[i] = J1 + imp;}}

       	for(i =0;i<N_sites;i++){Jleft[i] = Jright[left[i]];
        Jdown[i] = Jup[down[i]]; /*printf("Jr[%d]:%lf \t Jl[%d]:%lf \t  Ju[%d]:%lf \t Jd[%d]:%lf \n",i,Jright[i],i,Jleft[i],i,Jup[i],i,Jdown[i]);*/}

		for(i=0;i<N_sites;i++){sitepr[i] = 1;}

		indx=0;
	    while(indx!=0){
	    	i=(int)N_sites*genrand_real2();
	    	if((sitepr[i] == 1)&&(sitepr[right[i]] == 1)&&(sitepr[left[i]] == 1)&&(sitepr[up[i]] == 1)&&(sitepr[down[i]] == 1)){
	    		sitepr[i] = 0;
	    		indx-=1;}}

	    for(i=0;i<N_sites;i++){ avsite+=sitepr[i];}

		for(i=0;i<N_sites;i++){
		xi=i%L;
		if(xi%2==0){spin[i] = -1 * sitepr[i];}
		else{spin[i] = 1*sitepr[i];}}
		measurement(spin,right,left,up,down,right_up,right_down,left_up,left_down,Jright,Jleft,Jup,Jdown);
		E1=Ener;
		Mf1=Magn_ferro;
		Ms1=magx*magx+magy*magy;



		for(i=0;i<N_sites;i++){
		yi=i/L;
		if(yi%2==0){spin[i] = -1*sitepr[i];}
		else{spin[i] = 1*sitepr[i];}}
		measurement(spin,right,left,up,down,right_up,right_down,left_up,left_down,Jright,Jleft,Jup,Jdown);
		E2=Ener;
		Mf2=Magn_ferro;
		Ms2=magx*magx+magy*magy;


		if(E1<E2){
		satmag=Ms1;
		for(i=0;i<N_sites;i++){
		xi=i%L;
		if(xi%2==0){spin[i] = -1 * sitepr[i];}
		else{spin[i] = 1*sitepr[i];}}}

		else{
		satmag=Ms2;
		for(i=0;i<N_sites;i++){
		yi=i/L;
		if(yi%2==0){spin[i] = -1*sitepr[i];}
		else{spin[i] = 1*sitepr[i];}}}

		if(opti==1){ /*option 1 for random configuration;*/
		for(i=0;i<N_sites;i++){
		k = genrand_real2();
		if(k<0.5){spin[i] = -1*sitepr[i];}
		else{spin[i] = 1*sitepr[i];}/*printf("%d\n",spin[i]);*/} }

		if(opti==2){/*option 2 is for ordered configuration*/
		for(i=0;i<N_sites;i++){
		spin[i] = 1*sitepr[i];/*printf("%d\n",spin[i]);*/}}

}
