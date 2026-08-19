/*@ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 **
 **   Project      : MAGEMin
 **   License      : GNU GENERAL PUBLIC LICENSE Version 3, 29 June 2007
 **   Developers   : Nicolas Riel, Boris Kaus
 **   Contributors : Moccetti, N. B., Dominguez, H., Assunção J., Green E., Dolejš, D., Berlie N., and Rummel L.
 **   Organization : Institute of Geosciences, Johannes-Gutenberg University, Mainz
 **   Contact      : nriel[at]uni-mainz.de, kaus[at]uni-mainz.de
 **
 ** ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ @*/
/**
  Function to calculate chemical potential of endmembers/pure phases for thermocalc database
*/

#include <math.h>
#include <stdio.h>
#include <time.h>
#include <string.h>
#include <complex.h> 

#include "../MAGEMin.h"
#include "../toolkit.h"
#include "TC_gem_function.h"


#define eps 1e-8
#define max_ox 15


PP_ref TC_G_EM_function(	int 		 EM_dataset, 
							int 		 len_ox,
							int         *id,
							double 		*bulk_rock, 
							double 		*apo, 
							double 		 P, 
							double 		 T, 
							char 		*name, 
							char		*state			
){
	/* Get thermodynamic data */
	EM_db EM_return;
	int i, p_id = find_EM_id(name);
	EM_return   = Access_EM_DB(p_id, EM_dataset);
	
	/* Get composition (in molar amount) */
	double composition[len_ox];
	for (i = 0; i < len_ox; i ++){
		composition[i] = EM_return.Comp[id[i]];
	}

	double t0, 	p0, 	R;
	double pth, theta, 	vv;
	double enthalpy, 	entropy, volume;
	double cpa, cpb, cpc, cpd;
	double alpha0, kappa0, kappa0p, kappa0pp, dkappa0dT;
	double cpterms, n;
	double vterm 	= 0.0;
	double ta 		= 0.0;
	double tb 		= 0.0;
	double tc 		= 0.0;
	double kbar2bar = 1e3;
	double RTlnf 	= 0.0;
	double water_density_gcm3 = 0.0;		/** g/cm3, from the Pitzer & Sterner (1994) molar volume below; only set when name=="H2O" - used by G_DEW_function's solvent EOS (DEW2019) */
	t0 				= 298.15;
	p0 				= 0.001;
	R  				= 0.0083144; 
	
	enthalpy 		= EM_return.input_1[0];
	entropy  		= EM_return.input_1[1];
	volume   		= EM_return.input_1[2];
	
	cpa      		= EM_return.input_2[0];
	cpb      		= EM_return.input_2[1];
	cpc      		= EM_return.input_2[2];
	cpd      		= EM_return.input_2[3];
	
	alpha0   		= EM_return.input_3[0];
	kappa0   		= EM_return.input_3[1];
	kappa0p  		= EM_return.input_3[2];
	kappa0pp 		= EM_return.input_3[3];	

	cpterms  		= cpa* (T - t0) +   cpb* (pow(T,2.0) - pow(t0,2.0))/2.0 - 
                                        cpc* (1.0/T - 1.0/t0) + 
                                   2.0* cpd* (pow(T,0.5) - pow(t0,0.5))     - 
							   T* (2.0* cpa* (log(pow(T,0.5)) - log(pow(t0,0.5))) 
                               + cpb* (T - t0) - 
							   cpc/2.0* (pow(T,-2.) - pow(t0,-2.0)) - 2.0* cpd* (pow(T,-0.5) - pow(t0,-0.5)));
							   
	n        		= EM_return.Comp[max_ox];
	
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

        // printf("vterm: %f\n", vterm);
        // printf("Ares: %f\n", Ares);

		if (strcmp( name, "H2O") == 0){
			water_density_gcm3 = 18.015/vsub;		/** vsub in cm3/mol -> g/cm3 */
		}
	}
	/**
	 	here we use the CORK EOS to calculate G_O2 (see Holland & Powell, 1991) 
		Critical Temperature and pressures are taken from Holland & Powell 1991 for H2
	*/
	else if(strcmp( name, "H2") == 0){
		double Tc, Pc;
		// if(strcmp( name, "O2") == 0){
		// 	Tc     =  154.75;
		// 	Pc     =  0.05;	
		// }
		// else{ //(strcmp( name, "H2") == 0){	
			Tc     =  41.2;
			Pc     =  0.0211;	
		// }

		double a0 	  =  5.45963e-5;
		double a1     = -8.63920e-6;
		double b0     =  9.18301e-4;
		double c0     = -3.30558e-5;
		double c1 	  =  2.30524e-6;
		double d0     =  6.93054e-7;
		double d1     = -8.38293e-8;

		double a	  = a0*pow(Tc,5.0/2.0)/Pc + a1*pow(Tc,3.0/2.0)/Pc*T;
		double b      = b0*Tc/Pc;
		double c      = c0*Tc/pow(Pc,3.0/2.0) + c1/pow(Pc,3.0/2.0)*T;
		double d      = d0*Tc/pow(Pc,2.0) + d1/pow(Pc,2.0)*T;

		vterm 		  =(R*T/P + b - (a*R*sqrt(T))/((R*T+b*P)*(R*T + 2.0*b*P)) + c*sqrt(P) + d*P)/100.0; 
        // RTlnf        = R*T*log(1000.0 * P) + b*P + a/(b*sqrt(T)) * (log(R*T + b*P) - log(R*T + 2.0*P*b)) + 2.0/3.0*c*P*sqrt(P) + d/2.0*P*P;
	}
    else if(strcmp( name, "O2") == 0){
        vterm = 0.0;
    }
	else {
		ta     = (1. + kappa0p)/(1. + kappa0p + kappa0 * kappa0pp);
		tb     = (kappa0p + pow(kappa0p,2.0) - (kappa0 * kappa0pp))/(kappa0 * (1. + kappa0p));
		tc     = (1. + kappa0p + kappa0 * kappa0pp)/(kappa0p + pow(kappa0p,2.0) - kappa0 * kappa0pp);
		vterm  = vv*((P-p0)*(1.-ta)+ta*(-pow(1.+tb*(P-pth),(1.0-tc))+pow(1.0 + tb * (p0 - pth),(1.0 - tc)))/(tb* (tc - 1.)))/((1. - ta) + ta* pow(1. + tb * p0,(-tc)));
	}

    
	double gbase = (enthalpy - T*entropy + cpterms + vterm + RTlnf);	
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
	PP_ref_db.phase_density =  water_density_gcm3;	/** g/cm3, only meaningful for name=="H2O" (0.0 otherwise) */
	PP_ref_db.phase_shearModulus  =  (EM_return.input_4[0]*kbar2bar + (P - p0)*(EM_return.input_4[1])*kbar2bar + (T - t0)*(EM_return.input_4[2]))/kbar2bar;

	return (PP_ref_db);
}


/**************************************************************************************/
/**************************************************************************************/
/* DEW2019 (Deep Earth Water) aqueous species standard-state EOS                      */
/*--------------------------------------------------------------------------          */
/* Port of MAGEMin.jl's DEW_19_gbase.jl / DEW_19_a-x.jl (Sverjensky and co-workers,   */
/* Deep Earth Water model, https://www.dewcommunity.org). See                          */
/* tools/DEW_implementation_plan.md for the porting strategy.                         */
/*                                                                                     */
/* Water density (rho_w, g/cm3) is obtained from TC_G_EM_function's own H2O standard  */
/* state (Pitzer & Sterner, 1994 EOS) - this keeps DEW's solvent consistent with    */
/* the rest of the TC-group database's own H2O endmember, and its molar-volume-derived */
/* density is already in g/cm3.                                                        */
/**************************************************************************************/

/* Water saturation pressure (liquid-vapor equilibrium), Shock et al. (1992).
   T_celsius in Celsius, returns bar. */
double DEW_psat(double T_celsius){
    double result  =  1.44021565e+00;
    result += -2.75944904e-02 * T_celsius;
    result +=  3.50602876e-04 * T_celsius*T_celsius;
    result += -2.44834016e-06 * T_celsius*T_celsius*T_celsius;
    result +=  1.57085668e-08 * T_celsius*T_celsius*T_celsius*T_celsius;
    return result;
}

/* Solvent g-function, Shock et al. (1992), used inside the Born term. T in K, P in bar, rho_w in g/cm3 */
double DEW_gSolvent(double T, double P, double rho_w){
    double min_rho_w = 1.0;
    double TC = T - 273.15;

    if (rho_w >= 1.0){
        return 0.0;
    }
    else if ((P >= 500.0) && (rho_w <= min_rho_w)){
        return 0.0;
    }
    else if ((P < 500.0) && (P >= 220.46) && (TC >= 373.917)){
        return 0.0;
    }
    else if ((P < 220.46) && (P >= 1.00) && (P > DEW_psat(TC))){
        return 0.0;
    }
    else{
        double agP   = -2.037662,    agPP  =  5.747000e-3, agPPP = -6.557892e-6;
        double bgP   =  6.107361,    bgPP  = -1.074337e-2, bgPPP =  1.268348e-5;
        double ag1   =  3.66666e-16, ag2   = -1.504956e-10, ag3  =  5.01799e-14;

        double aG = agP + agPP*TC + agPPP*TC*TC;
        double bG = bgP + bgPP*TC + bgPPP*TC*TC;

        double f = 0.0;
        if ((P <= 1000.0) && (TC >= 155.0) && (TC <= 355.0)){
            f += ( pow((TC-155.0)/300.0, 4.8) + ag1*pow((TC-155.0)/300.0, 16.0) ) *
                 ( ag2*pow(1000.0-P, 3.0) + ag3*pow(1000.0-P, 4.0) );
        }

        return aG*pow(1.0 - rho_w, bG) - f;
    }
}

/* Johnson & Norton (1991) dielectric constant of water. T in K, rho_w in g/cm3 */
double DEW_JN_epsilon(double T, double rho_w){
    double T_hat = T / 298.15;
    double r     = rho_w;      /* already g/cm3, i.e. normalized to the 1000 kg/m3 reference */

    double k0 = 1.0;
    double k1 = 14.70333593 / T_hat;
    double k2 = 212.8462733 / T_hat - 115.4445173 + 19.55210915 * T_hat;
    double k3 = -83.3034798 / T_hat + 32.13240048 * T_hat - 6.69409865 * T_hat * T_hat;
    double k4 = -37.86202045 / (T_hat*T_hat) + 68.87359646 / T_hat - 27.29401652;

    return k0 + k1*r + k2*r*r + k3*r*r*r + k4*r*r*r*r;
}

/* Sverjensky et al. (2014) dielectric constant of water ("DM"). T in K, rho_w in g/cm3 */
double DEW_DM_epsilon(double T, double rho_w){
    double a1 = -1.57637700752506e-3, a2 =  6.81028783422197e-2, a3 =  7.54875480393944e-1;
    double b1 = -8.01665106535394e-5, b2 = -6.87161761831994e-2, b3 =  4.74797272182151e0;

    double aTK1 = -5.8274486041453000E-02, aTK2 =  2.2389805995733700E+00, aTK3 = -2.0249736922093000E+01;
    double bTK1 =  5.7128535105795700E-02, bTK2 = -2.2591436511452200E+00, bTK3 =  2.6398103834434400E+01;

    double TC = T - 273.15;

    double rho_wExponent, expExponent;
    if (TC > 0.0){
        rho_wExponent = a1*TC + a2*sqrt(TC) + a3;
        expExponent   = b1*TC + b2*sqrt(TC) + b3;
    }
    else{
        rho_wExponent = aTK1*T + aTK2*sqrt(T) + aTK3;
        expExponent   = bTK1*T + bTK2*sqrt(T) + bTK3;
    }

    return exp(expExponent) * pow(rho_w, rho_wExponent);
}

/* Smoothed (tanh-weighted around 5000 bar) blend of the DM/JN dielectric constant models. T in K, P in bar, rho_w in g/cm3 */
double DEW_epsilon(double T, double P, double rho_w){
    double DM_value = DEW_DM_epsilon(T, rho_w);
    double JN_value = DEW_JN_epsilon(T, rho_w);

    double width  = 1000.0;
    double centre = 5000.0;
    double weight, result = 0.0;

    if (P < centre){
        weight = 0.5 + tanh((P - centre)/width) / 2.0;
        if (weight < 0.001){ weight = 0.0; }
        result = weight*DM_value + (1.0 - weight)*JN_value;
    }
    else{
        weight = 0.5 - tanh((P - centre)/width) / 2.0;
        if (weight < 0.001){ weight = 0.0; }
        result = weight*JN_value + (1.0 - weight)*DM_value;
    }
    return result;
}

/* Born function B. T in K, P in bar, rho_w in g/cm3 */
double DEW_born_B(double T, double P, double rho_w){
    double Epsilon = DEW_epsilon(T, P, rho_w);
    return -1.0 / Epsilon;
}

/* Debye-Hückel A parameter, as used by the DEW aqueous speciation solver (aq_min_iterative port). */
double DEW_Agamma(double T, double P, double rho_w){
    double Epsilon = DEW_epsilon(T, P, rho_w);
    return 1.824829238e6 * sqrt(rho_w) / pow(Epsilon*T, 1.5);
}

/* Debye-Hückel B parameter, as used by the DEW aqueous speciation solver (aq_min_iterative port). */
double DEW_Bgamma(double T, double P, double rho_w){
    double Epsilon = DEW_epsilon(T, P, rho_w);
    return 50.29158649 * sqrt(rho_w) / sqrt(Epsilon*T);
}

/**
    Standard-state Gibbs free energy of a DEW2019 aqueous species (revised-HKF with
    variable theta/c1/c2/omega0 and the DEW Born-solvation term). Direct port of
    DEW_19_gbase.jl. T_r, P_r, Psi, eta and theta are universal solvent-EOS constants
    shared by every species (not stored per-species, see TC_endmembers.h).
*/
AQ_ref G_DEW_function(     int              len_ox,
                            int             *id,
                            double          *ElEntropy,
                            double           P,
                            double           T,
                            double           rho_w,
                            char            *name          ){

    /* Get thermodynamic data */
    DEW_db DEW_return;
    int i, p_id = find_DEW_id(name);
    DEW_return  = Access_DEW_DB(p_id);

    /* universal solvent-EOS constants, shared by every species */
    double T_r   = 298.15;
    double P_r   = 1.0;
    double Psi   = 2600.0;
    double eta   = 694657.0;
    double theta = 228.0;

    double G_ref  = DEW_return.input_1[0];
    double S_ref  = DEW_return.input_1[1];
    double a1     = DEW_return.input_1[2];
    double a2     = DEW_return.input_1[3];
    double a3     = DEW_return.input_2[0];
    double a4     = DEW_return.input_2[1];
    double c1     = DEW_return.input_2[2];
    double c2     = DEW_return.input_2[3];
    double rH     = DEW_return.input_3[0];
    double omega0 = DEW_return.input_3[1];
    double z      = DEW_return.input_3[2];
    double mu_Hp  = DEW_return.input_4[0];

    double Pbar = P * 1000.0;      /* P given in kbar, like every other TC_gem_function.c EOS */

    double theta2 = theta*theta;
    double theta3 = theta2*theta;
    double theta_inv_T = 1.0/(T - theta);

    double c1_theta3 = c1*theta3;
    double c2_theta  = c2*theta;
    double c_denom   = c1*theta2 + 2.0*c2;

    /* Term 1: G_ref */
    double gbase = G_ref;

    /* Term 2: pressure terms */
    double pressure_term = T*a1*theta_inv_T - a1*theta*theta_inv_T + a3*theta_inv_T;
    gbase += Pbar*pressure_term - P_r*pressure_term;

    /* Term 3: entropy term */
    gbase -= S_ref*(T - T_r);

    /* Term 4: log terms with c1, c2 */
    double log_arg_T  = T   + (-c1_theta3 - 2.0*c2_theta)/c_denom;
    double log_arg_Tr = T_r + (-c1_theta3 - 2.0*c2_theta)/c_denom;
    gbase += T*c2*log(log_arg_T)/theta2;

    /* complex reference-state term */
    double t_ = ( -T_r*c1*theta2*log(T_r) - T_r*c1*theta2 - T_r*c2*log(T_r)
                  + T_r*c2*log(T_r - c1_theta3/c_denom - 2.0*c2_theta/c_denom)
                  + c1*theta3*log(T_r) + c1*theta3 + c2*theta*log(T_r)
                  - c2*theta*log(T_r - c1_theta3/c_denom - 2.0*c2_theta/c_denom) + c2_theta );

    gbase -= T*t_/(T_r*theta2 - theta3);
    gbase -= T_r*c2*log(log_arg_Tr)/theta2;
    gbase += T_r*t_/(T_r*theta2 - theta3);

    /* Term 5: omega0 terms */
    gbase -= 5.7986499999999998e-5*omega0*(T - T_r) - 0.9872562762839302*omega0;

    /* Term 6: Born contribution */
    double gS = DEW_gSolvent(T, Pbar, rho_w);
    double born_contrib;
    if (z == 0.0){
        born_contrib = omega0;
    }
    else{
        born_contrib = -eta*z/(rH + gS) +
                        (eta*z/rH + omega0) / (1.0 + (z/rH + omega0/eta)*fabs(z)*gS/(z*z));
    }
    gbase -= (DEW_born_B(T, Pbar, rho_w) + 1.0)*born_contrib;

    /* Term 7: a2, a4 pressure-temperature term */
    double log_psi_P  = log(Pbar + Psi);
    double log_psi_Pr = log(P_r  + Psi);
    gbase += (T*a2 - a2*theta + a4)*(log_psi_P - log_psi_Pr)*theta_inv_T;

    /* Term 8: final c1, c2 terms */
    double log_arg_T2  = T   + (-c1_theta3 - c2_theta + theta*(c1*theta2 + c2))/c_denom;
    double log_arg_Tr2 = T_r + (-c1_theta3 - c2_theta + theta*(c1*theta2 + c2))/c_denom;

    gbase -= (T*c1*theta2   + T*c2)  *log(log_arg_T2) /theta2;
    gbase += (T_r*c1*theta2 + T_r*c2)*log(log_arg_Tr2)/theta2;

    /* Get composition + formation-reaction stoichiometry, sliced to the active oxide set */
    double comp[len_ox];
    double mu_comp[len_ox+1];
    for (i = 0; i < len_ox; i++){
        comp[i]    = DEW_return.Comp[id[i]];
        mu_comp[i] = DEW_return.MuComp[id[i]];
    }
    mu_comp[len_ox] = mu_Hp;

    /* HSC/SUPCRT conversion (mirrors DEW_19_gbase.jl's Sref block: MAGEMin's precomputed
       per-oxide reference-entropy table (z_b.ElEntropy) is already Tr-scaled AND already
       in kJ - unlike Julia's own raw per-element*Tr sum, which stays in J until the final
       /1000 - so cor must be added post-/1000, not folded into the J-scale gbase) */
    double cor = SUPCRT_to_HSC(ElEntropy, comp, len_ox);   /* kJ */

    double gbase_supcrt = gbase / 1000.0;
    double gbase_hsc    = gbase_supcrt + cor;

    /* special case for H+ / H2 -> gbase = 0 (matches DEW_19_gbase.jl) */
    if (strcmp(DEW_return.Name, "H+") == 0 || strcmp(DEW_return.Name, "H2") == 0){
        gbase_supcrt = 0.0;
        gbase_hsc    = 0.0;
    }

    /* fill structure to send back to main */
    AQ_ref AQ_ref_db;

    strcpy(AQ_ref_db.Name, name);
    for (i = 0; i < len_ox; i++){
        AQ_ref_db.Comp[i]   = comp[i];
        AQ_ref_db.MuComp[i] = mu_comp[i];
    }
    AQ_ref_db.MuComp[len_ox] = mu_comp[len_ox];

    AQ_ref_db.z             = z;
    AQ_ref_db.gbase_supcrt  = gbase_supcrt;
    AQ_ref_db.gbase         = gbase_hsc;
    AQ_ref_db.factor        = 1.0;

    return (AQ_ref_db);
}
