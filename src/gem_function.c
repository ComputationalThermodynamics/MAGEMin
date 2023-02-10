/**
  Function to calculate chemical potential of endmembers/pure phases  
*/

#include <math.h>
#include <stdio.h>
#include <time.h>
#include <string.h>
#include <complex.h> 

#include "MAGEMin.h"
#include "gem_function.h"
#include "toolkit.h"

#define eps 1e-8

/**
  compute the Gibbs Free energy from the thermodynamic database
*/
PP_ref G_EM_function(		int 		 EM_database, 
							int 		 len_ox,
							double 		*bulk_rock, 
							double 		*apo, 
							double 		 P, 
							double 		 T, 
							char 		*name, 
							char		*state			
){
	/* Get thermodynamic data */
	struct EM_db EM_return;
	int i, p_id = find_EM_id(name);
	EM_return   = Access_EM_DB(p_id, EM_database);
	
	/* Get composition (in molar amount) */
	double composition[len_ox];
	for (i = 0; i < len_ox; i ++){
		composition[i] = EM_return.Comp[i];
	}
	
	/**
		NOTE: The function below is specific for tc_ds633 and might be different 
		for other databases. Ideally, we would therefore here call a seperate 
		routine depending on the EM_database.
    */

	/** if (EM_database == _tc_ds633_) { **/
	double t0, 	p0, 	R;
	double pth, theta, 	vv;
	double enthalpy, 	entropy, volume;
	double cpa, cpb, cpc, cpd;
	double alpha0, kappa0, kappa0p, kappa0pp, dkappa0dT;
	double cpterms, n;
	double vterm = 0.0;
	double ta = 0.0;
	double tb = 0.0;
	double tc = 0.0;
	double kbar2bar = 1e3;
	t0 = 298.15;
	p0 = 0.001;
	R  = 0.0083144; 
	
	enthalpy = EM_return.input_1[0];
	entropy  = EM_return.input_1[1];
	volume   = EM_return.input_1[2];
	
	cpa      = EM_return.input_2[0];
	cpb      = EM_return.input_2[1];
	cpc      = EM_return.input_2[2];
	cpd      = EM_return.input_2[3];
	
	alpha0   = EM_return.input_3[0];
	kappa0   = EM_return.input_3[1];
	kappa0p  = EM_return.input_3[2];
	kappa0pp = EM_return.input_3[3];	

	cpterms  = cpa* (T - t0) +          cpb* (pow(T,2.0) - pow(t0,2.0))/2.0 - 
                                        cpc* (1.0/T - 1.0/t0) + 
                                   2.0* cpd* (pow(T,0.5) - pow(t0,0.5))     - 
							   T* (2.0* cpa* (log(pow(T,0.5)) - log(pow(t0,0.5))) 
                               + cpb* (T - t0) - 
							   cpc/2.0* (pow(T,-2.) - pow(t0,-2.0)) - 2.0* cpd* (pow(T,-0.5) - pow(t0,-0.5)));
							   
	n        = EM_return.Comp[len_ox];
	
	char liq_tail[] = "L";
	if ( EndsWithTail(name, liq_tail) == 1 ) {
		dkappa0dT        = EM_return.input_3[4];
		pth              = 0.0;
		vv               = volume * exp(alpha0 * (T-t0));
		kappa0           = kappa0 + (dkappa0dT * (T-t0));
	}
	else {
		theta = (double)(round(10636/(entropy*1e3/n + 6.44)));
		pth   = theta* alpha0* kappa0 / (exp(theta/t0) * pow(theta/t0,2.0) / pow(exp(theta/t0) - 1.,2.0)) * (1./(exp(theta/(T)) - 1.) - 1./(exp(theta/t0) - 1.));
		vv    = volume;
	}
	/* EOS After Pitzer and Sterner, 1994 - API, The Journal of Chemical Physics */
	if (strcmp( name, "H2O") == 0 || strcmp( name, "CO2") == 0 ){

		double p_bar = 1000.*P; //in bar
		double c1, c2, c3, c4, c5, c6, c7, c8, c9, c10; 

		if (strcmp( name, "H2O") == 0){
			c1  =  0.24657688e6 / T + 0.51359951e2;
			c2  =  0.58638965e0 / T - 0.28646939e-2 + 0.31375577e-4 * T;
			c3  = -0.62783840e1 / T + 0.14791599e-1 + 0.35779579e-3 * T +  0.15432925e-7 * pow(T,2.0);
			c4  = -0.42719875e0 - 0.16325155e-4 * T;
			c5  =  0.56654978e4 / T - 0.16580167e2 + 0.76560762e-1 * T;
			c6  =  0.10917883e0;
			c7  =  0.38878656e13 / pow(T,4.0) - 0.13494878e9 / pow(T,2.0) + 0.30916564e6 / T + 0.75591105e1;
			c8  = -0.65537898e5 / T + 0.18810675e3;
			c9  = -0.14182435e14 / pow(T,4.0) + 0.18165390e9 / pow(T,2.0) - 0.19769068e6 / T - 0.23530318e2;
			c10 =  0.92093375e5 / T + 0.12246777e3;
		}
		else {	// can only be CO2
			c1  =  0.18261340e7 / T + 	0.79224365e2;
			c2  =  						0.66560660e-4 	+ 0.57152798e-5 * T + 0.30222363e-9 * pow(T,2.0);
			c3 	= 						0.59957845e-2 	+ 0.71669631e-4 * T + 0.62416103e-8 * pow(T,2.0);
			c4  = -0.13270279e1 / T +  -0.15210731e0  	+ 0.53654244e-3 * T - 0.71115142e-7 * pow(T,2.0);
			c5  =  0.12456776e0 / T +   0.49045367e1    + 0.98220560e-2 * T + 0.55962121e-5 * pow(T,2.0);
			c6  = 				     	0.75522299e0;
			c7  = -0.39344644e+12 / pow(T,4.0) + 0.90918237e8 / pow(T,2.0) + 0.42776716e6 / T - 0.22347856e2;
			c8 	=  0.40282608e3 / T +   0.11971627e3;
			c9  =  0.22995650e8 / pow(T,2.0) - 0.78971817e5 / T - 0.63376456e2;
			c10 =  0.95029765e5 / T + 0.18038071e2;
		}

		
		/* solve for volume at P, T */
		int    err,  k;
		double vsub, yr;
		double R1     = 83.144;
		double data[] = {R1,T,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,p_bar};
		
		double x1     = 3.0;
		double x2     = R1*T/P;
		
		double e      = 1e-14;
		int maxiter   = 500;
		int mode      = 0;               														/** Mode is used to send the right *data (see root_finding.c) */
		
		vsub          =  BrentRoots(x1,x2,data,e,mode,maxiter, &yr, &k, &err);
		
		double r      =   1.0/vsub;
		double Ares   =   R1*T*( c1*r + (1.0/(c2 + c3*r + c4*pow(r, 2.0) + c5*pow(r, 3.0) + c6*pow(r, 4.0)) - 1.0/c2) - c7/c8*(exp(-c8*r) - 1.0) - c9/c10*(exp(-c10*r) - 1.0) );
		vterm         =   (Ares + p_bar*vsub + R1*T*(log( R1*T / vsub ) - 1.0)) * 1e-4;	
	}
	/* here we use the CORK EOS to calculate G_O2 (see Holland & Powell, 1991) */
	else if(strcmp( name, "O2") == 0){
		
		double a0 	  =  5.45963e-5;
		double a1     = -8.63920e-6;
		double b0     =  9.18301e-4;
		double c0     = -3.30558e-5;
		double c1 	  =  2.30524e-6;
		double d0     =  6.93054e-7;
		double d1     = -8.38293e-8;

		double Tc     =  154.75;												//critical T constant for O2
		double Pc     =  0.05;													//critical P constant for O2
		
		double a	  = a0*pow(Tc,5.0/2.0)/Pc + a1*pow(Tc,3.0/2.0)/Pc*T;
		double b      = b0*Tc/Pc;
		double c      = c0*Tc/pow(Pc,3.0/2.0) + c1/pow(Pc,3.0/2.0)*T;
		double d      = d0*Tc/pow(Pc,2.0) + d1/pow(Pc,2.0)*T;
		
		vterm 		  = R*T/P + b - (a*R*sqrt(T))/((R*T+b*P)*(R*T + 2*b*P)) + c*sqrt(P) + d*P; 
	}
	else {
		ta     = (1. + kappa0p)/(1. + kappa0p + kappa0 * kappa0pp);
		tb     = (kappa0p + pow(kappa0p,2.0) - (kappa0 * kappa0pp))/(kappa0 * (1. + kappa0p));
		tc     = (1. + kappa0p + kappa0 * kappa0pp)/(kappa0p + pow(kappa0p,2.0) - kappa0 * kappa0pp);
		vterm  = vv*((P-p0)*(1.-ta)+ta*(-pow(1.+tb*(P-pth),(1.0-tc))+pow(1.0 + tb * (p0 - pth),(1.0 - tc)))/(tb* (tc - 1.)))/((1. - ta) + ta* pow(1. + tb * p0,(-tc)));

	}
	double gbase = (enthalpy - T*entropy + cpterms + vterm);	
	
	double landaut, smax, vmax, sfdh, sfdhv, sfw, sfwv, sfn, sffac; 
	double god, sod, q, v;
	double tc0, q20, q2;
	
	sfn 	= 0.0;
	vmax 	= 0.0;
	smax    = 0.0;
	landaut = 0.0;
	god     = 0.0;
	
	if ( EndsWithTail(name, liq_tail) != 1 ) {
		if (EM_return.input_3[4] == 0.0){
			landaut = 0.; smax = -0.001; vmax = 0.; sfdh = 0.; sfdhv = 0.;
			sfw     = 0.; sfwv = 0.;     sfn  = 0.; sffac = 0.;
		}
		else if (EM_return.input_3[4] == 1.0){
			landaut = EM_return.input_3[5]; smax = EM_return.input_3[6]; vmax = EM_return.input_3[7];
			sfdh = 0.; sfdhv = 0.; sfw = 0.; sfwv = 0.; sfn = 0.; sffac = 0.;
		}
		else if (EM_return.input_3[4] == 2.0){		
			landaut = 0.; smax = -0.001; vmax = 0.;
			sfdh = EM_return.input_3[5]; sfdhv = EM_return.input_3[6]; sfw = EM_return.input_3[7];
			sfwv = EM_return.input_3[8]; sfn   = EM_return.input_3[9]; sffac = EM_return.input_3[10];
		}
		//printf(" %4s %+10f\n",name,sfn);
		if (sfn > 0.){
			if (strcmp( state, "ordered") == 0 ){
				god = 0.;
			}
			else if (strcmp( state, "disordered") == 0 ){
				if (sffac < 0.){
					sod = sffac * R * (log(1./(sfn+1.)) + sfn*log(sfn/(sfn+1.)))*(1./sffac-sfn)/(sfn+1.);
				}
				else {
					sod = sffac * R * (log(1./(sfn+1.)) + sfn*log(sfn/(sfn+1.)));
				}
				god = sfdh + P*sfdhv + T*sod;
			}
			else if (strcmp( state, "equilibrium") == 0 ){
				if (sffac < 0.){
					/* solve for volume at P, T */
					int    err,  k;
					double vsub, yr;
					double data[] = {sfdh,P,sfdhv,sfw,T,sfwv,sfn,R,sffac};
					double x1     = eps;
					double x2     = 1.0-eps;
					double e      = 1e-12;
					int maxiter   = 500;
					int mode      = 1;               /** Mode is used to send the right *data (see root_finding.c) */
					
					q             =  BrentRoots(x1,x2,data,e,mode,maxiter,&yr,&k,&err);
					sod           = (((1. + sfn*q)*log((1. + sfn*q)/(sfn+1.)) + sfn*(1.-q)*log(sfn*(1.-q)/(sfn+1.)) - sffac*(sfn*(1.-q)*log((1.-q)/(sfn+1.)) + sfn*(sfn+q)*log((sfn+q)/(sfn+1.)) ))/(sfn+1.));
				}
				else {
					/* Test function to define min/max */
					v = eps;
					double v1 = ( sfdh + P*sfdhv + (sfw + P*sfwv)*(2.*v - 1.) + sffac*sfn/(sfn + 1.)*R*T * log(sfn*pow(1. - v,2.0)/((1. + sfn*v)*(sfn + v))) );
					v = 1.0-eps;
					double v2 = ( sfdh + P*sfdhv + (sfw + P*sfwv)*(2.*v - 1.) + sffac*sfn/(sfn + 1.)*R*T * log(sfn*pow(1. - v,2.0)/((1. + sfn*v)*(sfn + v))) );
					
					/* solve for volume at P, T */
					double x1, x2;
					int    err,  k;
					double vsub, yr;
					double data[] = {sfdh,P,sfdhv,sfw,sfwv,sffac,sfn,R,T};
					
					if (check_sign(v1, v2) == 1) {	x1     = eps;	x2     = 1.0 - eps;	}
					else {							x1     = 0.;	x2     = 1.0 - eps;	}
					
					double e     = 1e-12;
					int maxiter  = 500;
					int mode     = 2;               /** Mode is used to send the right *data (see root_finding.c) */
					q            =  BrentRoots(x1,x2,data,e,mode,maxiter,&yr,&k,&err);
					sod          =  (sffac*((1.+sfn*q)*log((1. + sfn*q)/(sfn + 1.)) + sfn*(1. - q)*log((1. - q)/(sfn + 1.)) + sfn*(1. - q)*log(sfn*(1. - q)/(sfn + 1.)) + sfn*(sfn + q)*log((sfn + q)/(sfn + 1.))) / (sfn + 1.));
				}
				god              = sfdh + P*sfdhv + q*(sfw - sfdh + P*(sfwv - sfdhv)) - pow(q,2.0)*(sfw + P*sfwv) + R*T*sod;
			}
			else {
				printf("wrong state (HAS TO BE: ordered, disordered or equilibrium)");
			}
		}
		else if (smax > 0.0){
			tc0 = landaut;
			q20 = sqrt(1.0 - t0 / tc0);
			if (strcmp( state, "ordered") == 0 ){
				god = smax*tc0*(-(2./3.) + q20*(1.0 - pow(q20,2.)/3.)) - T*smax*(q20 - 1.0) + P*vmax*(q20 - 1.0);
			}
			else if (strcmp( state, "disordered") == 0 ){
				god = smax*tc0*q20*(1.0 - pow(q20,2.)/3.) - T*smax*q20 + P*vmax*q20;
			}
			else if (strcmp( state, "equilibrium") == 0 ){
				if (vmax == 0){
					tc  = tc0;
				}
				else{
					tc  = tc0 + P * vmax / smax;
				}
				if(T >  tc){
					q2  = 0.0;
				}
				else{
					q2  = pow((tc - T) / tc0, 0.5);
				}
				god = smax*(tc0*(q20*(1.0 - (1./3.)*pow(q20, 2.0)) + (1./3.)*pow(q2, 3.0)) - q2*tc) - T*smax*(q20 - q2) + P*vmax*q20;
			}
			else{
				printf("wrong state (HAS TO BE: ordered, disordered or equilibrium)");
			}
		}
		else{
			god = 0.0;
		}
		gbase = gbase + god;
		
 	}
 	
	/* fill structure to send back to main */
	PP_ref PP_ref_db;
	
	/* Calculate normalizing factor using bulk-rock composition */
	double factor  = 0.0;
	
	/* Calculate the number of atoms in the bulk-rock composition */
	double fbc     = 0.0;
	for (i = 0; i < len_ox; i++){
		fbc += bulk_rock[i]*apo[i];
	}
	
	/* Calculate the number of atom in the solution */
	double ape = 0.0;
	for (i = 0; i < len_ox; i++){
		ape += composition[i]*apo[i];
	}
	
	/* Calculate normalizing factor */
	factor = fbc/ape;

	strcpy(PP_ref_db.Name, name);
	for (i = 0; i < len_ox; i++){
		PP_ref_db.Comp[i] = composition[i];
	}
	PP_ref_db.gbase   =  gbase;
	PP_ref_db.factor  =  factor;
	PP_ref_db.phase_shearModulus  =  (EM_return.input_4[0]*kbar2bar + (P - p0)*(EM_return.input_4[1])*kbar2bar + (T - t0)*(EM_return.input_4[2]))/kbar2bar;


	// printf(" %4s %+10f\n",name,PP_ref_db.gbase);
	// for (i = 0; i < len_ox; i++){
	// 	printf("%+10f",PP_ref_db.Comp[i]*PP_ref_db.factor); 
	// }
	// printf("\n");

	return (PP_ref_db);
}


