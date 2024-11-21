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
#include <stdbool.h>
#include <time.h>
#include <string.h>
#include <complex.h> 

#include "../MAGEMin.h"
#include "../toolkit.h"
#include "SB_gem_function.h"

#define eps 			1e-10
#define R 				8.314462618  				// J/(mol·K)
#define LOG_EPS 	   -36.04365338911715 			// Example value for log_eps
#define SQRT_EPS 		1.4901161193847656e-08 		// Example value for sqrt_eps
#define val_infinity 	19.4818182068004875 		// Example value for val_infinity
#define kbar2bar 		1e3
#define bar2pa 			1e5
#define T0 				298.15
#define P0 				0.001

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

/*
 	Returns molar volume. :math:`[m^3]`
	** 'Borrowed' from Burnman **
*/
// double volume(	double pressure, 
// 				double temperature, 
// 				double V_0, 
// 				double F_0, 
// 				double K_0, 
// 				double Kprime_0, 
// 				double n, 
// 				double grueneisen_0, 
// 				double q_0, 
// 				double Debye_0		){
					

//     double dV = 1.0e-2 * V_0;

//     double a1_ii = 6.0 * grueneisen_0;  // EQ 47
//     double a2_iikk = -12.0 * grueneisen_0 + 36.0 * pow(grueneisen_0, 2.0) - 18.0 * params[5] * grueneisen_0;  // EQ 47

//     double b_iikk = 9.0 * K_0;  // EQ 28
//     double b_iikkmm = 27.0 * K_0 * (Kprime_0 - 4.0);  // EQ 29


//     // Finding the volume at a given pressure requires a root-finding scheme.
//     // Here we use brentq to find the root.

//     // Root-finding using brentq requires bounds to be specified.
//     // We do this using a bracketing function.
//     double args[] = {
//         pressure,
//         temperature,
//         V_0,
//         Debye_0,
//         n,
//         a1_ii,
//         a2_iikk,
//         b_iikk,
//         b_iikkmm,
//     };

//     double sol;
//     // try {
//         // The first attempt to find a bracket for root finding uses V_0 as a starting point
//         sol = bracket(delta_pressure, V_0, dV, args);
//     // } catch (Exception) {
//     //     // At high temperature, the naive bracketing above may try a volume guess that exceeds the point at which the bulk modulus goes negative at that temperature.
//     //     // In this case, we try a more nuanced approach by first finding the volume at which the bulk modulus goes negative,
//     //     // and then either (a) raising an exception if the desired pressure is less than the pressure at that volume,
//     //     // or (b) using that pressure to create a better bracket for brentq.
//     //     double _K_T(double V, double T, double *params) {
//     //         return isothermal_bulk_modulus_reuss(0.0, T, V, params[0], params[1], params[2], params[3], grueneisen_0, params[5], K_0, Kprime_0, params[8], params[9], params[10]);
//     //     }

//     //     double sol_K_T = bracket(_K_T, V_0, dV, args);
//     //     double V_crit = brentq(_K_T, sol_K_T[0], sol_K_T[1], args);
//     //     double P_min = pressure(temperature, V_crit, params[0], params[1], params[2], params[3], grueneisen_0, params[5], K_0, Kprime_0, params[8], params[9], params[10]);
//     //     if (P_min > pressure) {
//     //         fprintf(stderr, "The desired pressure is not achievable at this temperature. The minimum pressure achievable is %.2e Pa.\n", P_min);
//     //         exit(EXIT_FAILURE);
//     //     } else {
//     //         try {
//     //             sol = bracket(delta_pressure, V_crit - dV, dV, args);
//     //         } catch (Exception) {
//     //             fprintf(stderr, "Cannot find a volume, perhaps you are outside of the range of validity for the equation of state?\n");
//     //             exit(EXIT_FAILURE);
//     //         }
//     //     }
//     // }

//     return brentq(delta_pressure, sol[0], sol[1], args);
// }

/*
 	Helmholtz free energy of lattice vibrations in the Debye model [J].
    It is important to note that this does NOT include the zero
    point energy for the lattice.  As long as you are
    calculating relative differences in F, this should cancel anyways.
	** 'Borrowed' from Burnman **
*/
double _helmholtz_free_energy(double T, double debye_T, double n) {
    if (T <= eps) {
        return 0.0;
    }
    double x = debye_T / T;
    double F = n * R * T * (3.0 * log(1.0 - exp(-x)) - debye_fn_cheb(x));
    return F;
}

/*
    Returns the Helmholtz free energy at the pressure and temperature
    of the mineral [J/mol]
	** 'Borrowed' from Burnman **
*/
double helmholtz_free_energy(double pressure, double temperature, double volume, double V_0, double T_0, double F_0, double K_0, double Kprime_0, double n, double grueneisen_0, double q_0, double Debye_0) {
    double x = V_0 / volume;
    double f = 0.5 * (pow(x, 2.0 / 3.0) - 1.0);
    double Debye_T = debye_temperature(V_0 / volume, grueneisen_0, q_0, Debye_0);

    double F_quasiharmonic = _helmholtz_free_energy(temperature, Debye_T, n)
        - _helmholtz_free_energy(T_0, Debye_T, n);

    double b_iikk = 9.0 * K_0;  // EQ 28
    double b_iikkmm = 27.0 * K_0 * (Kprime_0 - 4.0);  // EQ 29

    double F = F_0
        + 0.5 * b_iikk * f * f * V_0
        + (1.0 / 6.0) * V_0 * b_iikkmm * f * f * f
        + F_quasiharmonic;

    return F;
}

/*
	Returns Gibbs energy
	** 'Borrowed' from Perple_X **
*/
double compute_G0(	double t, 
					double p,
					double *V,

					double f0,
					double n,
					double v0,
					double k00,
					double k0p,
					double z00,
					double gamma0,
					double q0,
					double cme,
					double g0,
					double g0p ) {

    // Declare all variables at the beginning of the function
    double nr9, c1, c2, c3, aii, aiikk, aiikk2, aii2;
    double nr9t0, beta, gammel, t1, t2, nr9t, tht, tht0;
    double delt2, dfth, dfth0, root;
    double v, v23, f, df, d2f, dfc, d2fc, fel, dfel, d2fel;
    double z, a2f, da, dtht, d2tht, dtht0, d2tht0, fpoly;
    double fpoly0, etht, letht, d2fth, etht0, letht0, d2fth0;
    double f1, df1, dv, a;
	double G0;
    int itic, ibad;
    bool bad;

    // v0 = -EM_return.input_1[2];
    nr9 =  -9.0 * n * R;
    c1 = -9.0 * (-v0) * k00;
    c2 = k0p / 2.0 - 2.0;
    c3 = 3.0 * c1 * c2;
    aii = 6.0 * gamma0;
    aiikk = -12.0*gamma0 + 36.0*pow(gamma0,2.0) - 18.0*q0*gamma0;
    aiikk2 = 0.5 * aiikk;
    aii2 = 0.5 * aii;
    nr9t0 = nr9 * T0;
    // beta = (3.0*k00*g0p-5.0*g0);
    // gammel = ((6.0*g0p-24.0+4.5*k0p)*k00-14.0*g0);
      
    t1 = z00 / t;
    t2 = t / T0;
    nr9t = nr9 * t;

    // Initial guess for volume
    tht = t1;
    tht0 = tht * t2;
    delt2 = t * t - T0 * T0;

    dfth = nr9t * gamma0 / v0 * (3.0 * plg(tht) / (tht * tht * tht) - log(1.0 - exp(-tht)));
    dfth0 = nr9t0 * gamma0 / v0 * (3.0 * plg(tht0) / (tht0 * tht0 * tht0) - log(1.0 - exp(-tht0)));

    root = k00 * ((2.0 + 2.0 * k0p) * (p + dfth - dfth0) + k00);

    if (root > 0.0) {
        v = ((2.0 + k0p) - sqrt(root) / k00) * v0 / (1.0 + k0p);
        if (v < v0 / 10.0 || v > v0 * 10.0) v = v0;
    } else {
        v = v0;
    }

    itic = 0;
    ibad = 4;
    bad  = true;

    while (true) {
        itic++;

        v23 = pow(v0 / v, 2.0 / 3.0);
        f = 0.5 * v23 - 0.5;
        df = -v23 / (3.0 * v);
        d2f = (5.0 / 9.0) * v23 / (v * v);

        dfc = (c3 * f + c1) * f * df;
        d2fc = (2.0 * c3 * f + c1) * df * df + (c3 * f + c1) * f * d2f;

        // fel = -beta / 2.0 * pow(v / v0, gammel) * delt2;
        // dfel = fel * gammel / v;
        // d2fel = dfel * (gammel - 1.0) / v;

        z = 1.0 + (aii + aiikk2 * f) * f;

        if (z < 0.0 || v / v0 > 100.0 || v / v0 < 0.01) break;

        root = sqrt(z);

        tht = t1 * root;
        tht0 = tht * t2;

        a2f = aii2 + aiikk2 * f;
        da = a2f / root;
        dtht = t1 * da * df;
        d2tht = t1 * ((aiikk2 / root - a2f * a2f / pow(z, 1.5)) * df * df + da * d2f);

        dtht0 = dtht * t2;
        d2tht0 = d2tht * t2;

        fpoly = 3.0 * plg(tht) / (tht * tht * tht);
        fpoly0 = 3.0 * plg(tht0) / (tht0 * tht0 * tht0);

        etht = exp(-tht);

        if (1.0 - etht < 0.0) break;

        letht = log(1.0 - etht);

        dfth = (letht - fpoly) * nr9t * dtht / tht;
        d2fth = ((4.0 * dtht * dtht / tht - d2tht) * (fpoly - letht) + dtht * dtht * etht / (1.0 - etht)) * nr9t / tht;

        etht0 = exp(-tht0);

        if (1.0 - etht0 < 0.0) break;

        letht0 = log(1.0 - etht0);

        dfth0 = (letht0 - fpoly0) * nr9t0 * dtht0 / tht0;
        d2fth0 = ((4.0 * dtht0 * dtht0 / tht0 - d2tht0) * (fpoly0 - letht0) + dtht0 * dtht0 * etht0 / (1.0 - etht0)) * nr9t0 / tht0;

        f1 = -dfc - dfth + dfth0 -p;

        df1 = -d2fc - d2fth + d2fth0;

        dv = f1 / df1;

        if (v - dv < 0.0) dv = v / 2.0;

        v = v - dv;

        if (itic > 100 || fabs(f1) > 1e40) {
            if (fabs(f1 / p) < 1e-10) ibad = 5;
            break;
        } else if (fabs(dv / (1.0 + v)) < 1e-10) {
            bad = false;
            break;
        }
    }

	*V = v;

    f = 0.5 * pow(v0 / v, 2.0 / 3.0) - 0.5;
    z = 1.0 + (aii + aiikk2 * f) * f;
    root = sqrt(z);

    tht = t1 * root;
    tht0 = tht * t2;

    a = f0 + c1 * f * f * (0.5 + c2 * f) + nr9 * (t / (tht * tht * tht) * plg(tht) - T0 / (tht0 * tht0 * tht0) * plg(tht0));

    G0 = a + p * v - t * cme;

	return G0;
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

	gbase = compute_G0(	T, P, &V,
						F0,
						n,
						V0,
						K0,
						Kp,
						z00,
						gamma0,
						q0,
						cme,
						g0,
						g0p );

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
