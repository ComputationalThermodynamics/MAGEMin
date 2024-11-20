/*@ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 **
 **   Project      : MAGEMin
 **   License      : GNU GENERAL PUBLIC LICENSE Version 3, 29 June 2007
 **   Developers   : Nicolas Riel, Boris Kaus
 **   Contributors : Dominguez, H., Assunção J., Green E., Berlie N., and Rummel L.
 **   Organization : Institute of Geosciences, Johannes-Gutenberg University, Mainz
 **   Contact      : nriel[at]uni-mainz.de, kaus[at]uni-mainz.de
 **
 ** ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ @*/
/**
  Function to calculate chemical potential of endmembers/pure phases for Stixrude database
*/

#include <math.h>
#include <stdio.h>
#include <time.h>
#include <string.h>
#include <complex.h> 

#include "../MAGEMin.h"
#include "../toolkit.h"
#include "SB_gem_function.h"

#define eps 1e-10
#define R 8.314462618  							// J/(mol·K)
#define LOG_EPS -36.04365338911715 				// Example value for log_eps
#define SQRT_EPS 1.4901161193847656e-08 		// Example value for sqrt_eps
#define val_infinity 19.4818182068004875 	// Example value for val_infinity
#define kbar2bar 1e3
#define bar2pa 1e5
#define T0 298.15
#define P0 0.001

double plg(double t) {
	double p4, dinc;

    double p0 	= exp(-t);
    double p1 	= 1.0;
    double p2 	= t * t;
    double p3 	= 2.0 * t;

    double plg 	= -2.1646464674222763831;

    int i 	= 1;
    while (i < 100000) {

        p4 		= (double)i;
        p1 		= p1 * p0;
        dinc 	= (p2 + (p3 + 2.0 / p4) / p4) * p1 / p4 / p4;
        plg 	= plg + dinc;

        if (fabs(dinc / (1.0 + fabs(plg))) < eps) {
            return plg;
        }

        i += 1;
    }

    return plg;
}

/**
	Calculate the thermal part of the Birch-Murnaghan equation of state.
	Returns the thermal part of the equation of state in [Pa].
	** 'Borrowed' from Burnman **
*/
double debye_temperature(double x, double grueneisen_0, double q_0, double Debye_0) {
    double f = 0.5 * (pow(x, 2.0 / 3.0) - 1.0);
    double a1_ii = 6.0 * grueneisen_0;  // EQ 47
    double a2_iikk = -12.0 * grueneisen_0 + 36.0 * pow(grueneisen_0, 2.0) - 18.0 * q_0 * grueneisen_0;  // EQ 47
    double nu_o_nu0_sq = 1.0 + a1_ii * f + 0.5 * a2_iikk * f * f;

    if (nu_o_nu0_sq > 0.0) {
        return Debye_0 * sqrt(nu_o_nu0_sq);
    } else {
        // printf("This volume (V = %.2f*V_0) exceeds the valid range of the thermal part of the slb equation of state.\n", 1.0 / x);
        return 0.0;
    }
}

/**
	Calculate the thermal part of the Birch-Murnaghan equation of state.
	Returns the thermal part of the equation of state in [Pa].
	** 'Borrowed' from Burnman **
*/
double grueneisen_parameter_slb(double V_0, double volume, double gruen_0, double q_0) {
    double x = V_0 / volume;
    double f = 0.5 * (pow(x, 2.0 / 3.0) - 1.0);
    double a1_ii = 6.0 * gruen_0;  // EQ 47
    double a2_iikk = -12.0 * gruen_0 + 36.0 * gruen_0 * gruen_0 - 18.0 * q_0 * gruen_0;  // EQ 47
    double nu_o_nu0_sq = 1.0 + a1_ii * f + 0.5 * a2_iikk * f * f;  // EQ 41
    return 1.0 / 6.0 / nu_o_nu0_sq * (2.0 * f + 1.0) * (a1_ii + a2_iikk * f);
}

/**
	Finite strain approximation for :math:`q`, the isotropic volume strain
	derivative of the grueneisen parameter.
	** 'Borrowed' from Burnman **
*/
double volume_dependent_q(double x, double grueneisen_0, double q_0) {
    double f = 0.5 * (pow(x, 2.0 / 3.0) - 1.0);
    double a1_ii = 6.0 * grueneisen_0;  // EQ 47
    double a2_iikk = -12.0 * grueneisen_0 + 36.0 * pow(grueneisen_0, 2.0) - 18.0 * q_0 * grueneisen_0;  // EQ 47
    double nu_o_nu0_sq = 1.0 + a1_ii * f + 0.5 * a2_iikk * f * f;  // EQ 41
    double gr = 1.0 / 6.0 / nu_o_nu0_sq * (2.0 * f + 1.0) * (a1_ii + a2_iikk * f);

    double q;
    if (fabs(grueneisen_0) < 1.0e-10) {
        q = 1.0 / 9.0 * (18.0 * gr - 6.0);
    } else {
        q = 1.0 / 9.0 * (18.0 * gr - 6.0 - 0.5 / nu_o_nu0_sq * (2.0 * f + 1.0) * (2.0 * f + 1.0) * a2_iikk / gr);
    }

    return q;
}


/**
	Finite strain approximation for eta_{s0}, the isotropic shear
	strain derivative of the grueneisen parameter.
	** 'Borrowed' from Burnman **
*/
double isotropic_eta_s(double x, double grueneisen_0, double eta_s_0, double q_0) {

    double f 			= 0.5 * (pow(x, 2.0 / 3.0) - 1.0);
    double a2_s 		= -2.0 * grueneisen_0 - 2.0 * eta_s_0;
    double a1_ii 		= 6.0 * grueneisen_0;
    double a2_iikk 		= -12.0 * grueneisen_0 + 36.0 * pow(grueneisen_0, 2.0) - 18.0 * q_0 * grueneisen_0;
    double nu_o_nu0_sq 	= 1.0 + a1_ii * f + 0.5 * a2_iikk * pow(f, 2.0);
    double gr 			= 1.0 / 6.0 / nu_o_nu0_sq * (2.0 * f + 1.0) * (a1_ii + a2_iikk * f);
    double eta_s 		= -gr - (0.5 * pow(nu_o_nu0_sq, -1.0) * pow((2.0 * f) + 1.0, 2.0) * a2_s);

    return eta_s;
}

/*
    Get the birch murnaghan shear modulus at a reference temperature, for a
    given volume. Returns shear modulus in [Pa] (the same units as in G_0).
    This uses a second order finite strain expansion
	** 'Borrowed' from Burnman **
*/
double shear_modulus_second_order(double volume, double V_0, double G_0, double Gprime_0, double K_0) {
    double x = V_0 / volume;
    double G = G_0 * pow(x, 5.0 / 3.0) * (1.0 - 0.5 * (pow(x, 2.0 / 3.0) - 1.0) * (5.0 - 3.0 * Gprime_0 * K_0 * 1e5 / G_0));

    return G;
}
/*
    compute the bulk modulus as per the third order
    birch-murnaghan equation of state.  Returns bulk
    modulus in the same units as the reference bulk
    modulus.  Pressure must be in [Pa].
	** 'Borrowed' from Burnman **
*/
double bulk_modulus(double volume, double V_0, double K_0, double Kprime_0) {
    double x = V_0 / volume;
    double f = 0.5 * (pow(x, 2.0 / 3.0) - 1.0);

    double K = pow(1.0 + 2.0 * f, 5.0 / 2.0) * (
        K_0
        + (3.0 * K_0 * Kprime_0 - 5 * K_0) * f
        + 27.0 / 2.0 * (K_0 * Kprime_0 - 4.0 * K_0) * f * f
    );

    return K;
}

/*
    Evaluate a Chebyshev series at points x.
    ** 'Borrowed' from Burnman **
	n = 17
*/
double chebyshev_representation[] = {
    2.707737068327440945 / 2.0,
    0.340068135211091751,
    -0.12945150184440869e-01,
    0.7963755380173816e-03,
    -0.546360009590824e-04,
    0.39243019598805e-05,
    -0.2894032823539e-06,
    0.217317613962e-07,
    -0.16542099950e-08,
    0.1272796189e-09,
    -0.987963460e-11,
    0.7725074e-12,
    -0.607797e-13,
    0.48076e-14,
    -0.3820e-15,
    0.305e-16,
    -0.24e-17
};

double _chebval(double x, const double* c, int len) {

    double c0, c1;

        double x2 = 2 * x;
        c0 = c[len - 2];
        c1 = c[len - 1];
        for (int i = 3; i <= len; i++) {
            double tmp = c0;
            c0 = c[len - i] - c1;
            c1 = tmp + c1 * x2;
        }
    
    return c0 + c1 * x;
}

/*
    Evaluate the Debye function using a Chebyshev series expansion coupled with
    asymptotic solutions of the function.  Shamelessly adapted from the GSL implementation
    of the same function (Itself adapted from Collected Algorithms from ACM).
    Should give the same result as debye_fn(x) to near machine-precision.
	** 'Borrowed' from Burnman **
*/
double debye_fn_cheb(double x) {
    
    double xcut = LOG_EPS;

    if (x < 2.0 * sqrt(2.0) * SQRT_EPS) {

        return 1.0 - 3.0 * x / 8.0 + x * x / 20.0;
    } else if (x <= 4.0) {
        double t = x * x / 8.0 - 1.0;
        double c = _chebval(t, chebyshev_representation, 17);

        return c - 0.375 * x;
    } else if (x < -(log(2.0) + LOG_EPS)) {
        int nexp = (int)floor(xcut / x);
        double ex = exp(-x);
        double xk = nexp * x;
        double rk = nexp;
        double sum = 0.0;
        for (int i = nexp; i > 0; i--) {
            double xk_inv = 1.0 / xk;
            sum *= ex;
            sum += (((6.0 * xk_inv + 6.0) * xk_inv + 3.0) * xk_inv + 1.0) / rk;
            rk -= 1.0;
            xk -= x;
        }

        return val_infinity / (x * x * x) - 3.0 * sum * ex;
    } else if (x < xcut) {
        double x3 = x * x * x;
        double sum = 6.0 + 6.0 * x + 3.0 * x * x + x3;

        return (val_infinity - 3.0 * sum * exp(-x)) / x3;
    } else {

        return ((val_infinity / x) / x) / x;
    }
}

/*
	Calculate the molar heat capacity at constant volume for a Debye solid.
	** 'Borrowed' from Burnman **
*/
double molar_heat_capacity_v(double T, double debye_T, double n) {
    if (T <= eps) {
        return 0.0;
    }
    double x = debye_T / T;
    double C_v = 3.0 * n * R * (4.0 * debye_fn_cheb(x) - 3.0 * x / (exp(x) - 1.0));
    return C_v;
}


/*
	Calculate the thermal energy for a Debye solid.
	** 'Borrowed' from Burnman **
*/
double thermal_energy(double T, double debye_T, double n) {
    if (T <= eps) {
        return 0.0;
    }
    double x = debye_T / T;
    double E_th = 3.0 * n * R * T * debye_fn_cheb(x);
    return E_th;
}

/* 
	Returns isothermal bulk modulus :math:`[Pa]`
	** 'Borrowed' from Burnman **
*/
double isothermal_bulk_modulus_reuss( 	double temperature, 
										double volume, 
										double T_0, 
										double V_0, 
										double grueneisen_0,
										double q_0, 
										double Debye_0, 
										double n, 
										double K_0, 
										double Kprime_0		){
											
    double debye_T = debye_temperature(V_0 / volume, grueneisen_0, q_0, Debye_0);
    double gr = grueneisen_parameter_slb(V_0, volume, grueneisen_0, q_0);

    // thermal energy at temperature T
    double E_th = thermal_energy(temperature, debye_T, n);

    // thermal energy at reference temperature
    double E_th_ref = thermal_energy(T_0, debye_T, n);

    // heat capacity at temperature T
    double C_v = molar_heat_capacity_v(temperature, debye_T, n);

    // heat capacity at reference temperature
    double C_v_ref = molar_heat_capacity_v(T_0, debye_T, n);

    double q = volume_dependent_q(V_0 / volume, grueneisen_0, q_0);

    double K = bulk_modulus(volume, V_0, K_0, Kprime_0)
        + (gr + 1.0 - q) * (gr / (volume/1e5)) * (E_th - E_th_ref)
        - (pow(gr, 2.0) / (volume/1e5)) * (C_v * temperature - C_v_ref * T_0);

    return K;
}

/*
	Returns shear modulus :math:`[Pa]`
	** 'Borrowed' from Burnman **
*/
double shear_modulus(	double temperature, 
						double volume, 
						double T_0, double V_0, 
						double grueneisen_0, double q_0, double Debye_0, double n, 
						double G_0, double Gprime_0, double K_0, double eta_s_0) {

    double debye_T 	= debye_temperature(V_0 / volume, grueneisen_0, q_0, Debye_0);
    double eta_s 	= isotropic_eta_s(V_0 / volume, grueneisen_0, eta_s_0, q_0);

    double E_th 	= thermal_energy(temperature, debye_T, n);
    double E_th_ref = thermal_energy(T_0, debye_T, n);

	return shear_modulus_second_order(volume, V_0, G_0, Gprime_0, K_0)
		- eta_s * (E_th - E_th_ref) / (volume/1e5);

}

PP_ref SB_G_EM_function(	int 		 EM_dataset, 
							int 		 len_ox,
							int         *id,
							double 		*bulk_rock, 
							double 		*apo, 
							double 		 Pkbar, 
							double 		 T, 
							char 		*name, 
							char		*state			
){

	/* Get thermodynamic data */
	EM_db_sb EM_return;
	int i, p_id = find_EM_id(name);
	EM_return   = Access_SB_EM_DB(p_id, EM_dataset);
	
	/* Get composition (in molar amount) */
	double composition[len_ox];
	for (i = 0; i < len_ox; i ++){
		composition[i] = EM_return.Comp[id[i]];
	}


	// double R  		= 8.31446261815324;
 	double P       = Pkbar * kbar2bar;
	/* declare the variables */
	double nr9, nr9T0, c1, c2, c3, aii, aiikk, as, aiikk2, aii2;
	double r23, r59, t1, t2, nr9t, tht, thT0;
	double b21, b22;
	double dfth, dfth0, root, V, V23;
	double f,df,d2f,dfc,d2fc,z,a2f,da,dtht,d2tht,dthT0,d2thT0,fpoly,fpoly0,etht,letht,d2fth,ethv,ethT0,lethT0,d2fth0,f1,df1,dv;
	double gamma, etas;
	double a,gbase;

	int    max_ite, itic, ibad, bad;

	max_ite 		= 128;

	double F0		=  EM_return.input_1[0];
	double n		=  EM_return.input_1[1];
	double V0 		= -EM_return.input_1[2];
	double K0 		=  EM_return.input_1[3];
	double Kp 		=  EM_return.input_1[4];
	double z00		=  EM_return.input_1[5];
	double gamma0 	=  EM_return.input_1[6];
	double q0 		=  EM_return.input_1[7];
	double etaS0	=  EM_return.input_1[8];
	double cme		=  EM_return.input_1[9];

	double g0 		=  EM_return.input_2[0]*1e9; // GPa to Pa
	double g0p 		=  EM_return.input_2[1];

    nr9 	= -9.0 * n * R;
    nr9T0 	= nr9 * T0;
    c1 		= -9.0 * (-V0) * K0;
    c2 		= Kp / 2.0 - 2.0;
    c3 		= 3.0 * c1 * c2;
    aii 	= 6.0 * gamma0;
	aiikk 	= -12.0*gamma0 + 36.0*pow(gamma0,2.0) - 18.0*q0*gamma0;
	as      = -(gamma0 + etaS0);
	aiikk2 	= 0.5 * aiikk;
    aii2 	= 0.5 * aii;
    // aiikk2 	= 0.5 * aii * (-2.0 + 6.0 * gamma0 - 3.0 * q0);
    // aii2 	= 3.0 * gamma0;
	b21		= (3.0*K0*g0p-5.0*g0);
	b22 	= ((6.0*g0p-24.0+4.5*Kp)*K0-14.0*g0);

    r23 = 2.0 / 3.0;
    r59 = 5.0 / 9.0;

    t1 		= z00 / T;
    t2 		= T / T0;
    nr9t 	= nr9 * T;

	tht 	= t1;
	thT0 	= tht * t2;

	dfth 	= nr9t * gamma0 / 	V0 * (3.0 * plg(tht) 	/ pow(tht,3.0) 	- log(1.0 - exp(-tht)));
	dfth0 	= nr9T0 * gamma0 / 	V0 * (3.0 * plg(thT0) 	/ pow(thT0,3.0) - log(1.0 - exp(-thT0)));

	root 	= K0 * ((2.0 + 2.0 * Kp) * (P + dfth - dfth0) + K0);

	V 		= 0.0;
	if (root > 0.0){
		V = ((2.0 + Kp) - sqrt(root) / K0) * V0 / (1.0 + Kp);
		if (V < V0 / 10.0 || V > V0 * 10.0){
			V = V0;
		}
	}
	else{
		V = V0;
	}

	ibad 	= 4;
	bad 	= 1;

	itic 	= 0;
	while (itic < max_ite){

		itic   += 1;

		V23 	= pow((V0 / V),r23);
		f 		= 0.5 * V23 - 0.5;
		df 		= -V23 / V / 3.0;
		d2f 	= r59 * V23 / pow(V,2.0);
		dfc 	= (c3 * f + c1) * f * df;
		d2fc 	= (2.0 * c3 * f + c1) * pow(df,2.0) + (c3 * f + c1) * f * d2f;

		z 		= 1.0 + (aii + aiikk2 * f) * f;

		if (z < 0.0 || V / V0 > 100.0 || V / V0 < 1e-2){
			break;
		// 	if (gv.verbose == 1){printf(" ERROR z or V/V0\n");}
		}

		root 	= sqrt(z);
		tht 	= t1 * root;
		thT0 	= tht * T / T0;
		a2f 	= aii2 + aiikk2 * f;
		da 		= a2f / root;
		dtht 	= t1 * da * df;
		d2tht 	= t1 * ((aiikk2 / root - pow(a2f,2.0) / pow(z,1.5)) * pow(df,2.0) + da * d2f);

		dthT0 	= dtht * t2;
		d2thT0 	= d2tht * t2;

		fpoly 	= 3.0 * plg(tht) / pow(tht,3.0);
		fpoly0 	= 3.0 * plg(thT0) / pow(thT0,3.0);
		etht 	= exp(-tht);

		if (1.0 - etht < 0.0){
			break;
		// 	printf("ERROR 1-etht\n");
		}

		letht 	= log(1.0 - etht);

		dfth 	= (letht - fpoly) * nr9t * dtht / tht;
		d2fth 	= ((4.0 * pow(dtht,2.0) / tht - d2tht) * (fpoly - letht) + pow(dtht,2.0) * etht / (1.0 - etht)) * nr9t / tht;
		ethT0 	= exp(-thT0);

		if (1.0 - ethT0 < 0.0){
			break;
		// 	printf("ERROR 1-thT0\n");
		}

		lethT0 	= log(1.0 - ethT0);

		dfth0 	= (lethT0 - fpoly0) * nr9T0 * dthT0 / thT0;
		d2fth0 	= ((4.0 * pow(dthT0,2.0) / thT0 - d2thT0) * (fpoly0 - lethT0) + pow(dthT0,2.0) * ethT0 / (1.0 - ethT0)) * nr9T0 / thT0;

		f1 		= -dfc - dfth + dfth0 - P;
		df1 	= -d2fc - d2fth + d2fth0;
		dv 		= f1 / df1;

		if (V - dv < 0.0){
			dv = V / 2.0;
		}

		V -= dv;
		
		if (itic >= max_ite || fabs(f1) > 1e40){
			if (fabs(f1 / P) < 0.0){
				ibad = 5;
				// if (gv.verbose == 1){printf("ERROR abs(f1/p)\n");}
			}
			break;
		}
		else if( fabs(dv / (1.0 + V)) < eps ){
				bad = 0;
				break;
		}
		
	}

	// if (bad == 1){
	// 	if (gv.verbose == 1){printf("ERROR bad\n");}
	// }

	/* get helmoltz energy:*/
	f 		= 0.5 * pow((V0 / V),r23) - 0.5;
	z 		= 1.0 + (aii + aiikk2 * f) * f;
	root 	= sqrt(z);

	tht 	= t1 * root;
	thT0 	= tht * t2;

	/* helmholtz energy */
	a 		= F0 + c1 * pow(f,2.0) * (0.5 + c2 * f) + nr9 * (T / pow(tht,3.0) * plg(tht) - T0 / pow(thT0,3.0) * plg(thT0));
	gbase 	= a + P * V - T * cme;


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
	PP_ref_db.gbase   =  gbase/kbar2bar;
	PP_ref_db.factor  =  factor;

	PP_ref_db.phase_shearModulus  = shear_modulus( 	T, 	V, 
													T0, V0, 
													gamma0, q0, z00,  n, 
													g0, g0p, K0, etaS0)/1e8;


	PP_ref_db.phase_bulkModulus  =  isothermal_bulk_modulus_reuss(	T, 	V, 
																	T0, V0, 
																	gamma0,
																	q0, 
																	z00, 
																	n, 
																	K0*bar2pa, 
																	Kp		)/1e8;

	// printf("gbase %4s %+10f\n",name,PP_ref_db.gbase);
	// printf("phase_shearModulus %4s %+10f\n",name,PP_ref_db.phase_shearModulus);
	// for (i = 0; i < len_ox; i++){
	// 	printf("%+10f",PP_ref_db.Comp[i]*PP_ref_db.factor); 
	// }
	// printf("\n");

	return (PP_ref_db);
}
