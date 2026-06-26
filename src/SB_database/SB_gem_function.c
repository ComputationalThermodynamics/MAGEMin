/*@ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 **
 **   Project      : MAGEMin
 **   License      : GNU GENERAL PUBLIC LICENSE Version 3, 29 June 2007
 **   Developers   : Nicolas Riel, Boris Kaus
 **   Contributors : Nickolas B. Moccetti, Dominguez, H., Assunção J., Green E., Berlie N., and Rummel L.
 **   Organization : Institute of Geosciences, Johannes-Gutenberg University, Mainz
 **   Contact      : nriel[at]uni-mainz.de, kaus[at]uni-mainz.de
 **
 ** ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ @*/
/**
  Function to calculate chemical potential of endmembers/pure phases for Stixrude database
  Most of the functions have been borrowed from Burnman (R. Myhill) and Perple_X (J. Connolly)
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
#define T0 				300.0
#define P0 				0.001

/* 0 = legacy (Perple_X-style) solver, 1 = burnman-style (Brent volume solve + 3rd order shear),
   2 = same as 1 but with HeFESTo's analytic vibrational/spinodal volume bounds */
static int SB_eos_formulation = 0;

void SB_set_eos_formulation(int mode){
	SB_eos_formulation = mode;
}

/* 0 = compute_G0() reproduces its original (pre-correction) behavior exactly: a
   non-converged Newton iteration silently returns whatever volume it landed
   on. 1 = apply the fix found by comparing against Perple_X's own GSTX
   routine (rlib.f): honor the 'bad' convergence flag - which the original
   Fortran already uses to destabilize the phase on failure, but which was
   computed-and-discarded in this C port - by returning NAN instead (caught
   by SB_G_EM_function's centralized NaN/Inf guard), and tighten the loop's
   v/v0 sanity bound from the port-introduced [0.01,100] to [0.1,10], matching
   both Perple_X's own initial-guess bound and HeFESTo's analytic fallback
   range (sb_volume_bounds()). Off by default to keep default behavior
   unchanged; set via gv.SB_eos_cor / --SB_eos_cor. */
static int SB_eos_correction = 0;

void SB_set_eos_correction(int mode){
	SB_eos_correction = mode;
}

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
    Get the birch murnaghan shear modulus at a reference temperature, for a
    given volume. Returns shear modulus in [Pa] (the same units as in G_0).
    This uses a third order finite strain expansion. Every mineral in
    burnman's SLB_2011/SLB_2024 databases uses this order ("should be
    preferred, as it is more thermodynamically consistent" per burnman).
	** 'Borrowed' from Burnman **
*/
double shear_modulus_third_order(double volume, double V_0, double G_0, double Gprime_0, double K_0, double Kprime_0) {
    double x      = V_0 / volume;
    double f      = 0.5 * (pow(x, 2.0 / 3.0) - 1.0);
    double K_0_pa = K_0 * 1e5;

    double G = pow(1.0 + 2.0 * f, 5.0 / 2.0) * (
        G_0
        + (3.0 * K_0_pa * Gprime_0 - 5.0 * G_0) * f
        + (6.0 * K_0_pa * Gprime_0 - 24.0 * K_0_pa - 14.0 * G_0 + 4.5 * K_0_pa * Kprime_0) * f * f
    );

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

    /* MAGEMin's sb endmember database stores n (atoms per formula unit) with a
       sign baked in for compute_G0's own Newton-iteration algebra; thermal_energy()/
       molar_heat_capacity_v() are direct burnman ports and require the physical
       (positive) atom count - see the same fix in compute_G0_burnman(). */
    double n_abs = fabs(n);

    // thermal energy at temperature T
    double E_th = thermal_energy(temperature, debye_T, n_abs);

    // thermal energy at reference temperature
    double E_th_ref = thermal_energy(T_0, debye_T, n_abs);

    // heat capacity at temperature T
    double C_v = molar_heat_capacity_v(temperature, debye_T, n_abs);

    // heat capacity at reference temperature
    double C_v_ref = molar_heat_capacity_v(T_0, debye_T, n_abs);

    double q = volume_dependent_q(V_0 / volume, grueneisen_0, q_0);

    /* volume is in MAGEMin's native J/bar units; 1 J/bar = 1e-5 m^3 (1 bar = 1e5 Pa),
       so volume/1e5 converts to m^3 to match E_th/C_v's SI (J) units - NOT volume/1e4. */
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
						double G_0, double Gprime_0, double K_0, double Kprime_0, double eta_s_0) {

    double debye_T 	= debye_temperature(V_0 / volume, grueneisen_0, q_0, Debye_0);
    double eta_s 	= isotropic_eta_s(V_0 / volume, grueneisen_0, eta_s_0, q_0);

    /* see isothermal_bulk_modulus_reuss()/compute_G0_burnman() for why fabs(n) is needed */
    double n_abs = fabs(n);
    double E_th 	= thermal_energy(temperature, debye_T, n_abs);
    double E_th_ref = thermal_energy(T_0, debye_T, n_abs);

    double G_BM;
    if (SB_eos_formulation == 1 || SB_eos_formulation == 2){
        G_BM = shear_modulus_third_order(volume, V_0, G_0, Gprime_0, K_0, Kprime_0);
    }
    else{
        G_BM = shear_modulus_second_order(volume, V_0, G_0, Gprime_0, K_0);
    }

	/* volume is in MAGEMin's native J/bar units; volume/1e5 converts to m^3 (see
	   isothermal_bulk_modulus_reuss()) - NOT volume/1e4. */
	return G_BM - eta_s * (E_th - E_th_ref) / (volume/1e5);

}

/*
    Returns thermal expansivity. :math:`[1/K]`
    ** 'Borrowed' from Burnman **
*/
double thermal_expansivity(     double V_0, 
                                double volume, 
                                double grueneisen_0, 
                                double q_0, 
                                double Debye_0, 
                                double n, 
                                double K, 
                                double temperature		){

    double debye_T, C_v, gr_slb, alpha;
                                    
    debye_T  = debye_temperature(V_0 / volume, grueneisen_0, q_0, Debye_0);
    C_v      = molar_heat_capacity_v(temperature, debye_T, n);
    gr_slb   = grueneisen_parameter_slb(V_0, volume, grueneisen_0, q_0);
    alpha    = gr_slb * C_v / volume / K;
    
    return alpha;
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
	Pressure residual P(V,T) - P_target for the SLB EOS, matching burnman's
	_delta_pressure (slb.py). Used to solve for V via bracket + Brent's method
	instead of compute_G0's Newton iteration. Re-uses thermal_energy(), already
	defined above in this file (depends on the Chebyshev-Debye machinery), so
	this stays local to SB_gem_function.c rather than going through toolkit.c's
	AFunction()/BrentRoots() dispatcher.

	data = { P_target, T, V_0, T_0, Debye_0, n, a1_ii, a2_iikk, b_iikk, b_iikkmm, bel_0, gel }
	** 'Borrowed' from Burnman **
*/
double sb_pressure_residual(double v, double *data) {
	double P_target  = data[0];
	double T         = data[1];
	double V_0       = data[2];
	double T_0       = data[3];
	double Debye_0   = data[4];
	double n         = data[5];
	double a1_ii     = data[6];
	double a2_iikk   = data[7];
	double b_iikk    = data[8];
	double b_iikkmm  = data[9];
	double bel_0     = data[10];
	double gel       = data[11];

	double f           = 0.5 * (pow(V_0 / v, 2.0 / 3.0) - 1.0);
	double nu_o_nu0_sq = 1.0 + a1_ii * f + 0.5 * a2_iikk * f * f;

	/* Outside this range the finite-strain Debye temperature expansion goes
	   imaginary (same condition burnman's _debye_temperature() raises on) -
	   signal invalidity with NaN instead of silently substituting debye_T=0,
	   which would fabricate a wrong-but-finite root. sb_solve_volume/sb_brent_root
	   treat NaN as failure and the caller falls back to compute_G0(). */
	if (nu_o_nu0_sq <= 0.0) {
		return NAN;
	}
	double debye_T  = Debye_0 * sqrt(nu_o_nu0_sq);

	double E_th     = thermal_energy(T,   debye_T, n);
	double E_th_ref = thermal_energy(T_0, debye_T, n);
	double gr       = 1.0 / 6.0 / nu_o_nu0_sq * (2.0 * f + 1.0) * (a1_ii + a2_iikk * f);

	double Pel = 0.0;
	if (bel_0 != 0.0) {
		Pel = 0.5 * gel * bel_0 * pow(v / V_0, gel) * (T * T - T_0 * T_0) / v;
	}

	return (1.0 / 3.0) * pow(1.0 + 2.0 * f, 5.0 / 2.0) * (b_iikk * f + 0.5 * b_iikkmm * f * f)
		+ gr * (E_th - E_th_ref) / v + Pel - P_target;
}

/*
	Brent's method root-finder, adapted from BrentRoots (toolkit.c) but calling
	sb_pressure_residual() directly instead of dispatching through AFunction(mode,...).
	RootBracketed()/Minimum() are shared with toolkit.c via toolkit.h.
*/
double sb_brent_root(double x1, double x2, double *data, double Tolerance, int maxIterations, int *error) {

	double FPP = 1e-11, nearzero = 1e-40;
	double result, AA, BB, CC, DD, EE, FA, FB, FC, Tol1, PP, QQ, RR, SS, xm;
	int i;
	bool done;

	result = 0.0; EE = 0.0; CC = 0.0;
	i = 0; done = false; *error = 0;
	AA = x1; BB = x2; FA = sb_pressure_residual(AA,data); FB = sb_pressure_residual(BB,data);

	/* RootBracketed() treats NaN as "bracketed" (NaN fails both > 0 and < 0
	   tests), so a NaN residual (volume outside the EOS's valid range, see
	   sb_pressure_residual) must be caught explicitly here or it would be
	   chased by the iteration below and silently propagate a NaN/bogus root. */
	if (isnan(FA) || isnan(FB)) {
		*error = 1;
		return x1;
	}

	if (!(RootBracketed(FA,FB))) {
		*error = 1;
		return x1;
	}

	FC = FB;
	do {
		if (!(RootBracketed(FC,FB))) {
			CC = AA; FC = FA; DD = BB - AA; EE = DD;
		}
		if (fabs(FC) < fabs(FB)) {
			AA = BB; BB = CC; CC = AA;
			FA = FB; FB = FC; FC = FA;
		}
		Tol1 = 2.0 * FPP * fabs(BB) + 0.5 * Tolerance;
		xm = 0.5 * (CC-BB);
		if ((fabs(xm) <= Tol1) || (fabs(FA) < nearzero)) {
			result = BB;
			done = true;
		}
		else {
			if ((fabs(EE) >= Tol1) && (fabs(FA) > fabs(FB))) {
				SS = FB/FA;
				if (fabs(AA - CC) < nearzero) {
					PP = 2.0 * xm * SS;
					QQ = 1.0 - SS;
				}
				else {
					QQ = FA/FC;
					RR = FB/FC;
					PP = SS * (2.0 * xm * QQ * (QQ - RR) - (BB-AA) * (RR - 1.0));
					QQ = (QQ - 1.0) * (RR - 1.0) * (SS - 1.0);
				}
				if (PP > nearzero) QQ = -QQ;
				PP = fabs(PP);
				if ((2.0 * PP) < Minimum(3.0*xm*QQ-fabs(Tol1 * QQ), fabs(EE * QQ))) {
					EE = DD; DD = PP/QQ;
				}
				else {
					DD = xm; EE = DD;
				}
			}
			else {
				DD = xm; EE = DD;
			}
			AA = BB;
			FA = FB;
			if (fabs(DD) > Tol1) BB = BB + DD;
			else {
				if (xm > 0.0) BB = BB + fabs(Tol1);
				else BB = BB - fabs(Tol1);
			}
			FB = sb_pressure_residual(BB,data);
			if (isnan(FB)) {
				*error = 1;
				return x1;
			}
			i++;
		}
	} while ((!done) && (i < maxIterations));

	if (i >= maxIterations) *error = 2;
	return result;
}

/*
	Solves for the equilibrium volume via a golden-ratio bracket expansion
	followed by Brent's method, matching burnman's volume() (slb.py), instead
	of compute_G0's Newton iteration. Returns 0 on success, 1 if no bracket
	could be found within maxiter steps.
	** 'Borrowed' from Burnman **
*/
int sb_solve_volume(double *V, double P_target, double *data, double V_0) {

	double dx = fabs(1.0e-2 * V_0);
	double ratio = 1.618033988749895;
	int maxiter = 100;

	double x0 = V_0;
	double f0 = sb_pressure_residual(x0,data);
	double x_left = x0 - dx, x_right = x0 + dx;
	double f_left = sb_pressure_residual(x_left,data);
	double f_right = sb_pressure_residual(x_right,data);

	/* overshot zero, try making dx smaller */
	if ((f0 - f_left) * (f_right - f0) < 0.0) {
		while ((f0 - f_left) * (f_right - f0) < 0.0 && dx > 1.0e-12 * V_0) {
			dx /= ratio;
			x_left = x0 - dx; x_right = x0 + dx;
			f_left = sb_pressure_residual(x_left,data);
			f_right = sb_pressure_residual(x_right,data);
		}
	}

	for (int i = 0; i < maxiter; i++) {
		if (RootBracketed(f_left,f0)) {
			int error;
			*V = sb_brent_root(x_left, x0, data, 1.0e-12*V_0, 100, &error);
			return (error == 0) ? 0 : 1;
		}
		if (RootBracketed(f0,f_right)) {
			int error;
			*V = sb_brent_root(x0, x_right, data, 1.0e-12*V_0, 100, &error);
			return (error == 0) ? 0 : 1;
		}
		dx *= ratio;
		x_left = x0 - dx; x_right = x0 + dx;
		f_left = sb_pressure_residual(x_left,data);
		f_right = sb_pressure_residual(x_right,data);
		if (x_left <= 0.0) { return 1; }
	}

	return 1;
}

/*
	Volume bounds within which a real equilibrium volume can exist, ported
	from HeFESTo (L. Stixrude & C. Lithgow-Bertelloni's own reference
	implementation, parset.f) - used by compute_G0_burnman_bounded()/
	sb_solve_volume_bounded() (SB_eos=2) to keep the bracket search from
	ever evaluating the residual outside the EOS's valid range, rather
	than discovering invalidity after the fact (compute_G0_burnman()/
	SB_eos=1's approach).

	Two independent bounds are computed and combined:
	- "vibrational" limit: roots of nu_o_nu0_sq(f) = vsquaredminimum
	  (HeFESTo requires some safety margin above the literal zero-crossing),
	  i.e. (1-vsquaredminimum) + a1_ii*f + 0.5*a2_iikk*f^2 = 0, converted
	  to volume via V = V_0*(2f+1)^(-3/2) (inverse of f = 0.5*((V_0/V)^(2/3)-1)).
	- "spinodal" limit: roots of the cold (BM3) bulk modulus going
	  non-positive, 1 + (3Kp-5)*f + 13.5*(Kp-4)*f^2 = 0 - a mechanical
	  stability bound, independent of the thermal one.
	Both are clamped to [V_0/10, V_0*10] (HeFESTo's own fallback floor/ceiling),
	then combined by taking the tighter of the two (vl=max, vu=min).
*/
static void sb_volume_bounds(double a1_ii, double a2_iikk, double Kp, double V_0,
								double *vl, double *vu) {

	const double vsquaredminimum = 0.1;
	double vlow = 1.0e-15, vupp = 1.0e15;

	/* --- vibrational limit --- */
	double c = 1.0 - vsquaredminimum;
	double b = a1_ii;
	double a = 0.5 * a2_iikk;
	double det = b*b - 4.0*a*c;
	double fextremum = 0.0, vextremum = vupp;
	if (a != 0.0) {
		fextremum = -b / (2.0*a);
		if (2.0*fextremum + 1.0 > 0.0) {
			vextremum = V_0 * pow(2.0*fextremum + 1.0, -1.5);
		}
	}

	if (det >= 0.0) {
		if (a == 0.0) {
			double f1 = (b != 0.0) ? -c/b : 0.0;
			if (f1 > 0.0) vlow = V_0 * pow(2.0*f1 + 1.0, -1.5);
			if (f1 < 0.0) vupp = V_0 * pow(2.0*f1 + 1.0, -1.5);
		} else {
			double sq = sqrt(det);
			double f1 = (-b - sq) / (2.0*a);
			double f2 = (-b + sq) / (2.0*a);
			double fmax = (f1 > f2) ? f1 : f2;
			double fmin = (f1 < f2) ? f1 : f2;
			if (fmax < 0.0) {
				vupp = V_0 * pow(2.0*fmax + 1.0, -1.5);
			} else if (fmin > 0.0) {
				vlow = V_0 * pow(2.0*fmin + 1.0, -1.5);
			} else {
				vupp = V_0 * pow(2.0*fmin + 1.0, -1.5);
				vlow = V_0 * pow(2.0*fmax + 1.0, -1.5);
			}
		}
	}
	/* require vibrational frequency to increase with increasing f (decreasing V) */
	if (fextremum > 0.0) vlow = (vlow > vextremum) ? vlow : vextremum;
	if (fextremum < 0.0) vupp = (vupp < vextremum) ? vupp : vextremum;

	double vlim_lo = ((vlow > V_0/10.0) ? vlow : V_0/10.0);
	double vlim_hi = ((vupp < V_0*10.0) ? vupp : V_0*10.0);

	/* --- spinodal limit (cold bulk modulus stays positive) --- */
	double csp = 1.0;
	double bsp = 3.0*Kp - 5.0;
	double asp = 13.5*(Kp - 4.0);
	double detsp = bsp*bsp - 4.0*asp*csp;
	double vsplow = 1.0e-15, vspupp = 1.0e15;
	if (detsp >= 0.0) {
		if (asp == 0.0) {
			double f1 = (bsp != 0.0) ? -csp/bsp : 0.0;
			vspupp = V_0 * pow(2.0*f1 + 1.0, -1.5);
		} else {
			double sq = sqrt(detsp);
			double f1 = (-bsp - sq) / (2.0*asp);
			double f2 = (-bsp + sq) / (2.0*asp);
			double fmax = (f1 > f2) ? f1 : f2;
			double fmin = (f1 < f2) ? f1 : f2;
			if (fmax < 0.0) {
				vspupp = V_0 * pow(2.0*fmax + 1.0, -1.5);
			} else if (fmin > 0.0) {
				vsplow = V_0 * pow(2.0*fmin + 1.0, -1.5);
			} else {
				vspupp = V_0 * pow(2.0*fmin + 1.0, -1.5);
				vsplow = V_0 * pow(2.0*fmax + 1.0, -1.5);
			}
		}
	}
	double vsp_lo = ((vsplow > V_0/10.0) ? vsplow : V_0/10.0);
	double vsp_hi = ((vspupp < V_0*10.0) ? vspupp : V_0*10.0);

	*vl = (vlim_lo > vsp_lo) ? vlim_lo : vsp_lo;
	*vu = (vlim_hi < vsp_hi) ? vlim_hi : vsp_hi;
}

/*
	Same golden-ratio bracket expansion as sb_solve_volume(), but hard-clamped
	to [vl,vu] (sb_volume_bounds()) instead of expanding without limit -
	mirrors HeFESTo's cage.f. Returns 0 on success, 1 if no root exists
	within [vl,vu] (genuine vibrational/spinodal instability at this P,T).
*/
int sb_solve_volume_bounded(double *V, double P_target, double *data, double V_0, double vl, double vu) {

	double ratio = 1.618033988749895;
	int maxiter = 100;

	if (!(V_0 > vl && V_0 < vu)) { return 1; }

	double x0 = V_0;
	double f0 = sb_pressure_residual(x0,data);
	double dx = fabs(1.0e-2 * V_0);
	double x_left  = (x0 - dx > vl) ? x0 - dx : vl;
	double x_right = (x0 + dx < vu) ? x0 + dx : vu;
	double f_left  = sb_pressure_residual(x_left,data);
	double f_right = sb_pressure_residual(x_right,data);

	for (int i = 0; i < maxiter; i++) {
		if (!isnan(f_left) && RootBracketed(f_left,f0)) {
			int error;
			*V = sb_brent_root(x_left, x0, data, 1.0e-12*V_0, 100, &error);
			return (error == 0) ? 0 : 1;
		}
		if (!isnan(f_right) && RootBracketed(f0,f_right)) {
			int error;
			*V = sb_brent_root(x0, x_right, data, 1.0e-12*V_0, 100, &error);
			return (error == 0) ? 0 : 1;
		}
		if (x_left <= vl && x_right >= vu) { return 1; }
		dx *= ratio;
		x_left  = (x0 - dx > vl) ? x0 - dx : vl;
		x_right = (x0 + dx < vu) ? x0 + dx : vu;
		f_left  = sb_pressure_residual(x_left,data);
		f_right = sb_pressure_residual(x_right,data);
	}

	return 1;
}

/* forward declaration: compute_G0_burnman()/compute_G0_burnman_bounded() fall back to this if no volume bracket is found */
double compute_G0(	double t, double p, double *V,
					double f0, double n, double v0, double k00, double k0p,
					double z00, double gamma0, double q0, double cme,
					double g0, double g0p, int nativeFe);

/*
	Burnman-style alternative to compute_G0(): solves for volume via
	bracket+Brent (sb_solve_volume) instead of a Newton iteration, then
	assembles Gibbs energy the way burnman does it
	(gibbs_free_energy = helmholtz_free_energy + P*V), reusing the
	already-present helmholtz_free_energy(). Electronic (fel) and magnetic
	(Gmag) corrections for native Fe phases are identical to compute_G0's
	(same physics, solver-independent) and applied the same way.
*/
double compute_G0_burnman(	double t,
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
							double g0p,
							int nativeFe) {

	double b1 = 0.0, b2 = 0.0, Tc = 0.0, tau, delt2;
	double fel, Gmag, D, Klro, Ksro, v;

	if (nativeFe == 1) { // fea
		b1 = 0.00388; b2 = 1.47960; Tc = 1043.01; tau = t/Tc;
	} else if (nativeFe == 2) { // fee
		b1 = 0.00411; b2 = 1.69270;
	} else if (nativeFe == 3) { // feg
		b1 = 0.00375; b2 = 1.4796;
	}

	double a1_ii    = 6.0 * gamma0;
	double a2_iikk  = -12.0*gamma0 + 36.0*pow(gamma0,2.0) - 18.0*q0*gamma0;
	double b_iikk   = 9.0 * k00;
	double b_iikkmm = 27.0 * k00 * (k0p - 4.0);

	/* MAGEMin's sb endmember database stores n (atoms per formula unit) with a
	   sign baked in for compute_G0's own Newton-iteration algebra; thermal_energy()/
	   helmholtz_free_energy() are direct burnman ports and require the physical
	   (positive) atom count, so use fabs(n) here and below. */
	double n_abs = fabs(n);

	double data[12] = { p, t, v0, T0, z00, n_abs, a1_ii, a2_iikk, b_iikk, b_iikkmm, b1, b2 };

	if (sb_solve_volume(&v, p, data, v0) != 0) {
		/* fall back to the legacy Newton solver if no bracket was found
		   (e.g. far outside the validated P-T range) */
		return compute_G0(t, p, V, f0, n, v0, k00, k0p, z00, gamma0, q0, cme, g0, g0p, nativeFe);
	}

	*V = v;

	delt2 = t*t - T0*T0;
	fel   = 0.0;
	if (nativeFe > 0) {
		fel = -b1 / 2.0 * pow(v / v0, b2) * delt2;
	}

	double F  = helmholtz_free_energy(p, t, v, v0, T0, f0, k00, k0p, n_abs, gamma0, q0, z00) + fel;
	double G0 = F + p * v - t * cme;

	/* Magnetic contribution to Gibbs free energy (Chin-Hertzman-Sundman model) */
	if (nativeFe == 1) { // only fea is assumed magnetic
		D    = (518.0/1125.0)+(11692.0/15975.0)*(1.0/0.4 - 1.0);
		Klro = 9.46 / D;
		Ksro = (474.0/497.0)*(1.0/0.4 - 1.0)*Klro;
		if (tau <= 1.0) {
			Gmag = t/D*9.46*( 1.0 - (79.0*pow(tau, -1.0)/140/0.4 + 474.0/497.0*(1.0/0.4 - 1.0)) * (pow(tau, 3.0)/6 + pow(tau, 9.0)/135 + pow(tau, 15.0)/600) );
		} else {
			Gmag = -t/D*9.46*(pow(tau, -5.0)/10.0 + pow(tau, -15.0)/315.0 + pow(tau, -25.0)/1500.0);
		}
		G0 += Gmag;
	}

	return G0;
}

/*
	Same as compute_G0_burnman(), but the volume solve is hard-bounded by
	HeFESTo's analytic vibrational+spinodal limits (sb_volume_bounds(),
	sb_solve_volume_bounded()) instead of an unbounded bracket expansion -
	this is the only difference from compute_G0_burnman(); the EOS formula
	(BM3 + single-Debye MGD thermal model) is identical and already matches
	HeFESTo's own reduced form for standard (single-Debye-temperature)
	mantle minerals, so it is reused verbatim here.
*/
double compute_G0_burnman_bounded(	double t,
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
									double g0p,
									int nativeFe) {

	double b1 = 0.0, b2 = 0.0, Tc = 0.0, tau, delt2;
	double fel, Gmag, D, Klro, Ksro, v;

	if (nativeFe == 1) { // fea
		b1 = 0.00388; b2 = 1.47960; Tc = 1043.01; tau = t/Tc;
	} else if (nativeFe == 2) { // fee
		b1 = 0.00411; b2 = 1.69270;
	} else if (nativeFe == 3) { // feg
		b1 = 0.00375; b2 = 1.4796;
	}

	double n_abs = fabs(n);
	double a1_ii    = 6.0 * gamma0;
	double a2_iikk  = -12.0*gamma0 + 36.0*pow(gamma0,2.0) - 18.0*q0*gamma0;
	double b_iikk   = 9.0 * k00;
	double b_iikkmm = 27.0 * k00 * (k0p - 4.0);

	double vl, vu;
	sb_volume_bounds(a1_ii, a2_iikk, k0p, v0, &vl, &vu);

	double data[12] = { p, t, v0, T0, z00, n_abs, a1_ii, a2_iikk, b_iikk, b_iikkmm, b1, b2 };

	if (sb_solve_volume_bounded(&v, p, data, v0, vl, vu) != 0) {
		/* fall back to the legacy Newton solver if no bracket was found
		   within HeFESTo's own vibrational/spinodal bounds - a genuine
		   instability at this P,T, not just a search-range artifact.
		   (SB_G_EM_function guards against this fallback itself being
		   NaN/Inf, regardless of which SB_eos formulation produced it.) */
		return compute_G0(t, p, V, f0, n, v0, k00, k0p, z00, gamma0, q0, cme, g0, g0p, nativeFe);
	}

	*V = v;

	delt2 = t*t - T0*T0;
	fel   = 0.0;
	if (nativeFe > 0) {
		fel = -b1 / 2.0 * pow(v / v0, b2) * delt2;
	}

	double F  = helmholtz_free_energy(p, t, v, v0, T0, f0, k00, k0p, n_abs, gamma0, q0, z00) + fel;
	double G0 = F + p * v - t * cme;

	/* Magnetic contribution to Gibbs free energy (Chin-Hertzman-Sundman model) */
	if (nativeFe == 1) { // only fea is assumed magnetic
		D    = (518.0/1125.0)+(11692.0/15975.0)*(1.0/0.4 - 1.0);
		Klro = 9.46 / D;
		Ksro = (474.0/497.0)*(1.0/0.4 - 1.0)*Klro;
		if (tau <= 1.0) {
			Gmag = t/D*9.46*( 1.0 - (79.0*pow(tau, -1.0)/140/0.4 + 474.0/497.0*(1.0/0.4 - 1.0)) * (pow(tau, 3.0)/6 + pow(tau, 9.0)/135 + pow(tau, 15.0)/600) );
		} else {
			Gmag = -t/D*9.46*(pow(tau, -5.0)/10.0 + pow(tau, -15.0)/315.0 + pow(tau, -25.0)/1500.0);
		}
		G0 += Gmag;
	}

	return G0;
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
					double g0p,
                    int nativeFe) {

    // Declare all variables at the beginning of the function
    double nr9, c1, c2, c3, aii, aiikk, aiikk2, aii2;
    double nr9t0, beta, gammel, t1, t2, nr9t, tht, tht0;
    double delt2, dfth, dfth0, root;
    double v, v23, f, df, d2f, dfc, d2fc;
    double z, a2f, da, dtht, d2tht, dtht0, d2tht0, fpoly;
    double fpoly0, etht, letht, d2fth, etht0, letht0, d2fth0;
    double b1, b2, fel, dfel, d2fel, tau, Tc, Klro, Ksro, D;
    double f1, df1, dv, a;
	double G0, Gmag;
    int itic, ibad;
    bool bad;

    /* Electronic contributions to alpha-epsilon-gamma Fe */
    fel=0.0; dfel=0.0; d2fel=0.0;
    if (nativeFe == 1) { // fea
        b1 = 0.00388; b2 = 1.47960; Tc=1043.01; tau=t/Tc;
    } else if (nativeFe == 2) { // fee
        b1 = 0.00411; b2 = 1.69270;
    } else if (nativeFe == 3) { // feg
        b1 = 0.00375; b2 = 1.4796;
    }

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

        // Electronic contribution to the Helmholtz free energy
        if (nativeFe > 0) {
            fel = -b1 / 2.0 * pow(v / v0, b2) * delt2;
            dfel = fel * b2 / v;
            d2fel = dfel * (b2 - 1.0) / v;
        }

        z = 1.0 + (aii + aiikk2 * f) * f;

        /* Perple_X's own GSTX routine (the original this was ported from)
           has no v/v0 sanity range in this loop at all - only v<=0; the
           [0.01,100] range below was introduced during the C port and is
           looser than both Perple_X's own initial-guess bound (v0/10..v0*10,
           see above) and HeFESTo's analytic fallback range. Under the
           correction, tighten to match those. */
        double v_v0_hi = SB_eos_correction ? 10.0  : 100.0;
        double v_v0_lo = SB_eos_correction ? 0.1   : 0.01;
        if (z < 0.0 || v / v0 > v_v0_hi || v / v0 < v_v0_lo) break;

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

        f1 = -dfc - dfth + dfth0 - dfel -p;

        df1 = -d2fc - d2fth + d2fth0 - d2fel;

        dv = f1 / df1;

        if (v - dv < 0.0) dv = v / 2.0;

        v = v - dv;

        if (itic > 2048 || fabs(f1) > 1e40) {
            if (fabs(f1 / p) < 1e-10) { ibad = 5; bad = false; }
            break;
        } else if (fabs(dv / (1.0 + v)) < 1e-10) {
            bad = false;
            break;
        }
    }

	*V = v;

    /* compute_G0's own convergence flag ('bad') was tracked but never acted
       upon - a non-converged Newton iteration would silently return whatever
       volume/G0 it last landed on. Perple_X's original GSTX routine instead
       destabilizes the phase on failure; under the correction, return NAN
       (caught by SB_G_EM_function's centralized NaN/Inf guard, which applies
       the same large-penalty-G treatment) instead of trusting a divergent v.
       *V is also reset to v0 so any downstream property (shear/bulk modulus,
       reported volume) computed from it isn't built on a garbage value. */
    if (SB_eos_correction && bad) {
        *V = v0;
        return NAN;
    }

    // Recompute (converged) electronic contribution
    if (nativeFe > 0) {
        fel = -b1 / 2.0 * pow(v / v0, b2) * delt2;
    } else {
        fel = 0.0;
    }

    f = 0.5 * pow(v0 / v, 2.0 / 3.0) - 0.5;
    z = 1.0 + (aii + aiikk2 * f) * f;
    root = sqrt(z);

    tht = t1 * root;
    tht0 = tht * t2;

    a = f0 + c1 * f * f * (0.5 + c2 * f) + nr9 * (t / (tht * tht * tht) * plg(tht) - T0 / (tht0 * tht0 * tht0) * plg(tht0)) + fel;

    G0 = a + p * v - t * cme;

    // Magnetic contribution to Gibbs free energy
    if (nativeFe == 1) { // only fea is assumed magnetic
        // Tc=1043.01 K | SD=9.46 J/mol/K | tau=t/Tc | p=0.4 (nondimentional standard for bcc (fea); see Roslyakova et al. 2016)
        D = (518.0/1125.0)+(11692.0/15975.0)*(1.0/0.4 - 1.0);
        Klro = 9.46 / D;
        Ksro = (474.0/497.0)*(1.0/0.4 - 1.0)*Klro;
        if (tau <= 1.0) {
            Gmag = t/D*9.46*( 1.0 - (79.0*pow(tau, -1.0)/140/0.4 + 474.0/497.0*(1.0/0.4 - 1.0)) * (pow(tau, 3.0)/6 + pow(tau, 9.0)/135 + pow(tau, 15.0)/600) );
        } else {
            Gmag = -t/D*9.46*(pow(tau, -5.0)/10.0 + pow(tau, -15.0)/315.0 + pow(tau, -25.0)/1500.0);
        }
        G0 += Gmag;
    }

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
        composition[i] = EM_return.Comp[i];
    }

    double P       = Pkbar * kbar2bar;

    /* declare the variables */
    double nr9, nr9T0, c1, c2, c3, cpterms, t0, aii, aiikk, as, aiikk2, aii2;
    double r23, r59, t1, t2, nr9t, tht, thT0;
    double b21, b22;
    double dfth, dfth0, root, V, V23;
    double f,df,d2f,dfc,d2fc,z,a2f,da,dtht,d2tht,dthT0,d2thT0,fpoly,fpoly0,etht,letht,d2fth,ethv,ethT0,lethT0,d2fth0,f1,df1,dv;
    double gamma, etas;
    double a, gbase;

    int    max_ite, itic, ibad, bad;

    max_ite 		= 128;

    double F0		=  EM_return.input_1[0]; // or S0 for O2
    double n		=  EM_return.input_1[1]; // or c1 for O2
    double V0 		= -EM_return.input_1[2]; // or c2 for O2
    double K0 		=  EM_return.input_1[3]; // or c3 for O2
    double Kp 		=  EM_return.input_1[4]; // or c5 for O2
    double z00		=  EM_return.input_1[5];
    double gamma0 	=  EM_return.input_1[6];
    double q0 		=  EM_return.input_1[7];
    double etaS0	=  EM_return.input_1[8];
    double cme		=  EM_return.input_1[9];
    double g0 		=  EM_return.input_2[0]*1e9; // GPa to Pa
    double g0p 		=  EM_return.input_2[1];
    int   nativeFe   =  0;

    if (strcmp(name, "fea") == 0) {
        nativeFe = 1;
    } else if (strcmp(name, "fee") == 0) {
        nativeFe = 2;
    } else if (strcmp(name, "feg") == 0) {
        nativeFe = 3;
    }

    if (strcmp(name, "O2") == 0) {
        t0=298.15;
        /* Cp integration */
        cpterms  		= n* (T - t0) +   V0* (pow(T,2.0) - pow(t0,2.0))/2.0 - 
                                        K0* (1.0/T - 1.0/t0) + 
                                2.0* Kp* (pow(T,0.5) - pow(t0,0.5))     - 
                            T* (2.0* n* (log(pow(T,0.5)) - log(pow(t0,0.5))) 
                            + V0* (T - t0) - 
                            K0/2.0* (pow(T,-2.) - pow(t0,-2.0)) - 2.0* Kp* (pow(T,-0.5) - pow(t0,-0.5)));

        /* gbase = (enthalpy - T*entropy + cpterms + vterm + RTlnf) */
        gbase = (0.0 - (T-t0)*F0 + cpterms + 0.0 + 0.0);
    } else if (SB_eos_formulation == 1) {
        gbase = compute_G0_burnman(	T, P, &V,
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
                            g0p,
                            nativeFe);
    } else if (SB_eos_formulation == 2) {
        gbase = compute_G0_burnman_bounded(	T, P, &V,
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
                            g0p,
                            nativeFe);
    } else {
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
                            g0p,
                            nativeFe);
    }

    if (isnan(gbase) || isinf(gbase)) {
        /* Genuine volume/Gibbs solve failure (vibrational or spinodal
           instability at this P,T) rather than a search-range artifact -
           can occur with any of the three SB_eos formulations. Nothing
           downstream in the LP/PGE pipeline guards against NaN/Inf
           endmember G's (they compare as neither favorable nor
           unfavorable and can corrupt the whole-system solve), so -
           matching HeFESTo's own gspec.f philosophy of forcing G
           arbitrarily large for an unstable species rather than letting
           the failure propagate - substitute a large but finite penalty.
           1e9 J (-> 1e6 kJ after the /kbar2bar conversion below) matches
           the existing penalty-G convention in simplex_levelling.c. */
        gbase = 1.0e9;
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
	PP_ref_db.gbase   =  gbase/kbar2bar;
	PP_ref_db.factor  =  factor;

    if (strcmp(name, "O2") == 0) { // Gasseous phase
        PP_ref_db.phase_shearModulus = 0.0;
        PP_ref_db.phase_bulkModulus  = 0.0;
        PP_ref_db.phase_expansivity  = 0.0;
        PP_ref_db.phase_cp           = 0.0;
    } else {
	PP_ref_db.phase_shearModulus  = shear_modulus( 	T, 	V,
													T0, V0,
													gamma0, q0, z00,  n,
													g0, g0p, K0, Kp, etaS0)/1e8;


	PP_ref_db.phase_bulkModulus  =  isothermal_bulk_modulus_reuss(	T, 	V, 
																	T0, V0, 
																	gamma0,
																	q0, 
																	z00, 
																	n, 
																	K0*bar2pa, 
																	Kp		)/1e8;


    PP_ref_db.phase_expansivity =  thermal_expansivity(     V0/1e6, 
                                                            V/1e6, 
                                                            gamma0, 
                                                            q0, 
                                                            z00, 
                                                            -n, 
                                                            PP_ref_db.phase_bulkModulus*1e9, 
                                                            T		);  


    PP_ref_db.phase_cp = molar_heat_capacity_v(T, z00, -n)/1e3;
    }

    return (PP_ref_db);
    
}
