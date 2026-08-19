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
    Pitzer & Sterner (1994) real-gas EOS for pure H2O and pure CO2, ported
    from xMELTS-master/sources/fluid.c's fluidPhase(). That original
    function handles a general H2O-CO2 MIXTURE (composition x[2]); this
    port keeps only the two pure-component special cases (x=[1,0] and
    x=[0,1]) that "gh" actually needs for the liquid's H2O/CO2 standard
    states, dropping the cross/mixing terms (aMix, the general branch of
    idealGas(), and the higher-order (2nd/3rd) T,rho-derivatives that
    fluidPhase() only needed for its Cp/V outputs, not G).

    References (as cited in the original xMELTS source):
      H2O: Pitzer KS and Sterner SM (1994) Equations of state valid
           continuously from zero to extreme pressures for H2O and CO2.
           J Chem Phys 101: 3111-6
      CO2: Sterner SM and Pitzer KS (1994) An equation of state for
           carbon dioxide valid from zero to extreme pressures.
           Contr Mineral Petrol 117: 362-74

    The physics/equations themselves (virial-type Pitzer-Sterner
    coefficients, density solved by Newton iteration from a Redlich-Kwong
    initial guess, then G = A + P/rho with a final shift onto the
    Berman (1988) elemental reference-state scale) are a faithful,
    unabridged port - not a re-derivation.
*/
#include <math.h>
#include <float.h>

#include "GH_fluid_eos.h"

#define nCOEFF 11

/* Pitzer & Sterner (1994) coefficients. Columns per row: [unused, T^-4, T^-2, T^-1, const, T, T^2] */
static const double h2o_c[nCOEFF][7] = {
    { 0.0,            0.0,             0.0,            0.0,            0.0,            0.0,           0.0 },
    { 0.0,            0.0,             0.0,   0.24657688e+6,  0.51359951e+2,           0.0,           0.0 },
    { 0.0,            0.0,             0.0,   0.58638965e+0, -0.28646939e-2, 0.31375577e-4,           0.0 },
    { 0.0,            0.0,             0.0,  -0.62783840e+1,  0.14791599e-1, 0.35779579e-3, 0.15432925e-7 },
    { 0.0,            0.0,             0.0,            0.0, -0.42719875e+0, -0.16325155e-4,           0.0 },
    { 0.0,            0.0,             0.0,   0.56654978e+4, -0.16580167e+2, 0.76560762e-1,           0.0 },
    { 0.0,            0.0,             0.0,            0.0,  0.10917883e+0,           0.0,           0.0 },
    { 0.0,   0.38878656e+13,  -0.13494878e+9,  0.30916564e+6,  0.75591105e+1,           0.0,           0.0 },
    { 0.0,            0.0,             0.0,  -0.65537898e+5,  0.18810675e+3,           0.0,           0.0 },
    { 0.0,  -0.14182435e+14,   0.18165390e+9, -0.19769068e+6, -0.23530318e+2,           0.0,           0.0 },
    { 0.0,            0.0,             0.0,   0.92093375e+5,  0.12246777e+3,           0.0,           0.0 },
};

static const double co2_c[nCOEFF][7] = {
    { 0.0,            0.0,             0.0,            0.0,            0.0,            0.0,           0.0 },
    { 0.0,            0.0,             0.0,   0.18261340e+7,  0.79224365e+2,           0.0,           0.0 },
    { 0.0,            0.0,             0.0,            0.0,  0.66560660e-4, 0.57152798e-5, 0.30222363e-9 },
    { 0.0,            0.0,             0.0,            0.0,  0.59957845e-2, 0.71669631e-4, 0.62416103e-8 },
    { 0.0,            0.0,             0.0,  -0.13270279e+1, -0.15210731e+0, 0.53654244e-3, -0.71115142e-7 },
    { 0.0,            0.0,             0.0,   0.12456776e+0,  0.49045367e+1, 0.98220560e-2, 0.55962121e-5 },
    { 0.0,            0.0,             0.0,            0.0,  0.75522299e+0,           0.0,           0.0 },
    { 0.0,  -0.39344644e+12,   0.90918237e+8,  0.42776716e+6, -0.22347856e+2,           0.0,           0.0 },
    { 0.0,            0.0,             0.0,   0.40282608e+3,  0.11971627e+3,           0.0,           0.0 },
    { 0.0,            0.0,    0.22995650e+8, -0.78971817e+5, -0.63376456e+2,           0.0,           0.0 },
    { 0.0,            0.0,             0.0,   0.95029765e+5,  0.18038071e+2,           0.0,           0.0 },
};

/* Ideal-gas contribution coefficients (Cooper 1982, for H2O; Berman 1988 /
   Stull & Prophet, for CO2), same values fluidPhase()'s idealGas() uses. */
static const double H2O_b[9]    = { 0.0, 0.134865, -5.005143, 4.006320, 0.012436, 0.973150, 1.279500, 0.969560, 0.248730 };
static const double H2O_beta[9] = { 0.0, 0.0,       0.0,      0.0,      833.0,    2289.0,   5009.0,   5982.0,   17800.0 };

/**
    Redlich-Kwong initial density guess (Edminster 1968 cubic solution),
    ported verbatim from fluid.c's redlichKwong(). Returns the gas-branch
    root (index 0), which is what fluidPhase() always uses as its initial
    guess before Newton-refining against the full Pitzer-Sterner EOS.
*/
static double gh_redlichKwong_Z(double t, double p, double b, double a2b){
    double n, m, arg;

    if (a2b <= 0.0) a2b = 0.001;

    n = b*p*(a2b-b*p-1.0)/3.0 - a2b*b*p*b*p - 2.0/27.0;
    m = b*p*(a2b-b*p-1.0) - 1.0/3.0;

    arg = n*n/4.0 + m*m*m/27.0;
    if (arg > 0.0){
        double term1 = - n/2.0 + sqrt(arg);
        double term2 = - n/2.0 - sqrt(arg);
        double z = (term1 >= 0.0) ?  pow( term1, 1.0/3.0) : -pow(-term1, 1.0/3.0);
        z        += (term2 >= 0.0) ?  pow( term2, 1.0/3.0) : -pow(-term2, 1.0/3.0);
        z        += 1.0/3.0;
        return z;
    }
    else if (arg == 0.0){
        return 1.0;
    }
    else{
        double PI = acos(-1.0);
        double cosPhi = (n > 0) ? -sqrt(- (n*n/4.0) / (m*m*m/27.0)) : sqrt(- (n*n/4.0) / (m*m*m/27.0));
        double phi = acos(cosPhi);
        double r1  = 2.0*sqrt(-m/3.0)*cos(phi/3.0) + 1.0/3.0;
        double r2  = 2.0*sqrt(-m/3.0)*cos(phi/3.0 + 2.0*PI/3.0) + 1.0/3.0;
        double r3  = 2.0*sqrt(-m/3.0)*cos(phi/3.0 + 4.0*PI/3.0) + 1.0/3.0;
        double zgas = r1;
        if (r2 > zgas) zgas = r2;
        if (r3 > zgas) zgas = r3;
        return zgas;
    }
}

/**
    Pitzer & Sterner (1994) pure-component Gibbs energy, ported from
    fluidPhase()'s density-solve + residual-term + ideal-gas-term +
    "G = A + P/rho" logic, specialized to a single pure component
    (dropping the H2O-CO2 mixing terms entirely).
    T in K, Pbar in bar. Returns G in J/mol (Berman 1988 elemental
    reference scale, matching the rest of MAGEMin's endmembers).
*/
double GH_pitzer_sterner_G(int is_H2O, double T, double Pbar){

    const double R  = 83.14241; /* cm^3-bar/mol-K */
    const double P0 = 1.01325;  /* reference pressure (bar) */

    double a, b;
    if (is_H2O){
        a = (111.3057 + 50.70033*exp(-0.982646e-2*(T-273.15)))*10.0e5;
        b = 14.6;
    }
    else{
        a = (73.03 - 0.0714*(T-273.15) + 2.157e-5*(T-273.15)*(T-273.15))*10.0e5;
        b = 29.7;
    }

    double c[nCOEFF];
    for (int i = 1; i < nCOEFF; i++){
        const double *row = is_H2O ? h2o_c[i] : co2_c[i];
        c[i] = row[1]/(T*T*T*T) + row[2]/(T*T) + row[3]/T + row[4] + row[5]*T + row[6]*T*T;
    }

    double Zgas = gh_redlichKwong_Z(T, Pbar, b/(82.05*T), a/(b*82.05*pow(T,1.5)));
    double rh   = Pbar/(Zgas*R*T);

    double dp = DBL_MAX;
    for (int count = 1; count <= 200 && fabs(dp) > 10.0*DBL_EPSILON; count++){
        double temp1  = c[3] + 2.0*c[4]*rh + 3.0*c[5]*rh*rh + 4.0*c[6]*rh*rh*rh;
        double dtemp1 = 2.0*c[4] + 6.0*c[5]*rh + 12.0*c[6]*rh*rh;
        double temp2  = c[2] + c[3]*rh + c[4]*rh*rh + c[5]*rh*rh*rh + c[6]*rh*rh*rh*rh;
        double dtemp2 = c[3] + 2.0*c[4]*rh + 3.0*c[5]*rh*rh + 4.0*c[6]*rh*rh*rh;

        double pr  = rh + c[1]*rh*rh - rh*rh*temp1/(temp2*temp2)
                   + c[7]*rh*rh*exp(-c[8]*rh) + c[9]*rh*rh*exp(-c[10]*rh);
        pr *= R*T;
        pr -= Pbar;

        double dpr = 1.0 + 2.0*c[1]*rh - 2.0*rh*temp1/(temp2*temp2)
                   - rh*rh*(temp2*temp2*dtemp1 - temp1*2.0*temp2*dtemp2)/(temp2*temp2*temp2*temp2)
                   + 2.0*c[7]*rh*exp(-c[8]*rh) - c[7]*c[8]*rh*rh*exp(-c[8]*rh)
                   + 2.0*c[9]*rh*exp(-c[10]*rh) - c[9]*c[10]*rh*rh*exp(-c[10]*rh);
        dpr *= R*T;

        dp  = -pr/dpr;
        rh += dp;
        if (rh < 0.0) rh = 10.0*DBL_EPSILON;
    }

    double term1      = c[2] + c[3]*rh + c[4]*rh*rh + c[5]*rh*rh*rh + c[6]*rh*rh*rh*rh;
    double dterm1drh  = c[3] + 2.0*c[4]*rh + 3.0*c[5]*rh*rh + 4.0*c[6]*rh*rh*rh;

    double term2a     = c[7]/c[8];
    double dterm2bdrh = -c[8]*exp(-c[8]*rh);
    double term2      = term2a*(exp(-c[8]*rh) - 1.0);
    double dterm2drh  = term2a*dterm2bdrh;

    double term3a     = c[9]/c[10];
    double dterm3bdrh = -c[10]*exp(-c[10]*rh);
    double term3      = term3a*(exp(-c[10]*rh) - 1.0);
    double dterm3drh  = term3a*dterm3bdrh;

    double Ar     = c[1]*rh + 1.0/term1 - 1.0/c[2] - term2 - term3;
    double dArdrh = c[1] - dterm1drh/(term1*term1) - dterm2drh - dterm3drh;

    double Ai, dAidrh;
    if (is_H2O){
        double f = -1.0 + log(rh*R*T/P0) + H2O_b[1] + H2O_b[2]/T + H2O_b[3]*(1.0-log(T));
        for (int i = 4; i <= 8; i++){
            f += H2O_b[i]*log(1.0 - exp(-H2O_beta[i]/T));
        }
        Ai     = R*T*f;
        dAidrh = R*T/rh;
    }
    else{
        const double k0 = 93.0*10.0, k1 = -13.409e2*10.0, k2 = 1.238e5*10.0, k4 = -0.002876*10.0, k5 = 6336.2*10.0;
        const double Gref = -94265.0*4.184*10.0, Sref = 51.072*4.184*10.0, Tref = 298.15;

        Ai = -k2/(2.0*T) - 2.0*k1*T/sqrt(Tref) - k2*T/(2.0*Tref*Tref) - k5*T/Tref
           - k0*T*log(T) + k5*log(T) + R*T*log(R*T*rh) + 4.0*k1*sqrt(T) - k4*T*T/2.0
           + k0*T*log(Tref) - R*T*log(P0) - R*T - Sref*T + k0*T + k4*Tref*T + k2/Tref
           - k5*log(Tref) + Tref*Sref + Gref - Tref*k0 - 2.0*k1*sqrt(Tref) - k4*Tref*Tref/2.0 + k5;
        dAidrh = R*T/rh;
    }

    double A     = R*T*Ar + Ai;
    double dAdrh = R*T*dArdrh + dAidrh;
    double G     = (A + rh*dAdrh)/10.0; /* bar*cm3/mol -> J/mol */

    /* shift the EOS's own 298.15K/1bar value onto the Berman (1988)
       elemental reference scale, exactly as fluidPhase() does          */
    if (is_H2O){
        const double refH2O_g = -228538.00;
        G += refH2O_g - (-46493.8016496949);
    }
    else{
        const double refCO2_g = -394341.00;
        G += refCO2_g - (-394450.0);
    }

    return G;
}

/**
    Pure-H2O Pitzer-Sterner "h - T*s" recomputation, needed for the two
    real gibbs.c contexts that use this instead of fluidPhase()'s own raw
    *g output (pMELTS' standalone "water" phase, gibbs.c ~line 2336's
    "gH2O = hH2O - t*sH2O", commented there as "used gH2O for calibration";
    and, less obviously, pMELTS' own LIQUID "H2O" basis species too - gh's own H2O
    liquid-table row has all-zero Kress dV/dP terms, matching real
    xMELTS' own table, which means the *downstream*, unconditional
    Birch-Murnaghan block (gibbs.c ~line 1493) always takes its "no
    pressure correction" else-branch for H2O specifically, and THAT
    branch unconditionally recomputes gl = hl - t*sl from whatever hl/sl
    the H2O-specific branch left behind - so the liquid's own gl ALSO
    ends up as h-T*s, not the raw g the H2O branch itself computed. Both
    findings made 2026-07-16 by direct comparison against a real gibbs()
    dispatch (not assumed from the source alone) - see
    [[gh-spn-liq-gbase-verification]].

    Derivation (avoids porting fluidPhase()'s full dA/dT machinery):
    internally, *g=A+p/rh, *h=A+p/rh-t*dAdt, *s=-dAdt satisfy h-t*s=g
    EXACTLY before fluidPhase()'s own reference-state shift - so the
    only reason h-T*s differs from g in the FINAL (shifted) output is
    that fluid.c applies three INDEPENDENT additive shifts to g, h, s
    (its own calibrated constants, not a thermodynamically consistent
    triple: shift_g != shift_h - T*shift_s). Since GH_pitzer_sterner_G
    already returns g_raw+shift_g, this recovers h_shifted-T*s_shifted
    as g_raw + shift_h - T*shift_s = GH_pitzer_sterner_G(...) +
    (shift_h-shift_g) - T*shift_s, using fluid.c's own refH2O
    g/h/s constants (-228538.00/-241816.00/188.72) and its own
    298.15K/1bar EOS-only baseline values (-46493.8016496949/
    9430.96281231262/187.572582951064). Verified exact (<0.001 kJ) at
    two independent (T,P) against real fluidPhase() output directly.
*/
double GH_pitzer_sterner_H2O_hTs_G(double T, double Pbar){
    const double refH2O_g = -228538.00,   baseline_g = -46493.8016496949;
    const double refH2O_h = -241816.00,   baseline_h =   9430.96281231262;
    const double refH2O_s =     188.72,   baseline_s =    187.572582951064;

    double shift_g = refH2O_g - baseline_g;
    double shift_h = refH2O_h - baseline_h;
    double shift_s = refH2O_s - baseline_s;

    double g_current = GH_pitzer_sterner_G(1, T, Pbar);
    return g_current + (shift_h - shift_g) - T*shift_s;
}

/* Cooper (1982) H2O ideal-gas term, factored out of GH_pitzer_sterner_G's
   is_H2O branch above (same formula, unchanged) so GH_pitzer_sterner_mix_G
   below can evaluate it at the mixture's shared density. */
static void GH_ideal_gas_H2O(double T, double rh, double *Ai, double *dAidrh){
    const double R  = 83.14241;
    const double P0 = 1.01325;
    double f = -1.0 + log(rh*R*T/P0) + H2O_b[1] + H2O_b[2]/T + H2O_b[3]*(1.0-log(T));
    for (int i = 4; i <= 8; i++){
        f += H2O_b[i]*log(1.0 - exp(-H2O_beta[i]/T));
    }
    *Ai     = R*T*f;
    *dAidrh = R*T/rh;
}

/* Stull & Prophet CO2 ideal-gas term, factored out of GH_pitzer_sterner_G's
   !is_H2O branch above (same formula, unchanged). */
static void GH_ideal_gas_CO2(double T, double rh, double *Ai, double *dAidrh){
    const double R  = 83.14241;
    const double P0 = 1.01325;
    const double k0 = 93.0*10.0, k1 = -13.409e2*10.0, k2 = 1.238e5*10.0, k4 = -0.002876*10.0, k5 = 6336.2*10.0;
    const double Gref = -94265.0*4.184*10.0, Sref = 51.072*4.184*10.0, Tref = 298.15;

    *Ai = -k2/(2.0*T) - 2.0*k1*T/sqrt(Tref) - k2*T/(2.0*Tref*Tref) - k5*T/Tref
        - k0*T*log(T) + k5*log(T) + R*T*log(R*T*rh) + 4.0*k1*sqrt(T) - k4*T*T/2.0
        + k0*T*log(Tref) - R*T*log(P0) - R*T - Sref*T + k0*T + k4*Tref*T + k2/Tref
        - k5*log(Tref) + Tref*Sref + Gref - Tref*k0 - 2.0*k1*sqrt(Tref) - k4*Tref*Tref/2.0 + k5;
    *dAidrh = R*T/rh;
}

double GH_pitzer_sterner_mix_G(double x_h2o, double T, double Pbar, double *dGdx_h2o){

    const double R  = 83.14241; /* cm^3-bar/mol-K */
    const double x0 = x_h2o, x1 = 1.0 - x_h2o;

    /* van der Waals one-fluid Redlich-Kwong mixing, used only to seed the
       initial density guess - the only place a genuine H2O-CO2 cross
       term enters (ported verbatim from fluidPhase()).                  */
    double aCO2  = (73.03 - 0.0714*(T-273.15) + 2.157e-5*(T-273.15)*(T-273.15))*10.0e5;
    double bCO2  = 29.7;
    double aH2O  = (111.3057 + 50.70033*exp(-0.982646e-2*(T-273.15)))*10.0e5;
    double bH2O  = 14.6;
    double a0CO2 = 46.0e6;
    double a0H2O = exp(4.881243 + 0.1823047e-2*(T-273.15) - 0.1712269e-4*(T-273.15)*(T-273.15)
                       + 6.479419e-8*(T-273.15)*(T-273.15)*(T-273.15))*10.0e5;
    double k     = exp(-11.071 + 5953.0/T - 2.746e6/(T*T) + 4.646e8/(T*T*T));
    double aMix  = sqrt(a0H2O*a0CO2) + k*82.05*82.05*pow(T, 2.5);

    double aRK = x0*x0*aH2O + x1*x1*aCO2 + 2.0*x0*x1*aMix;
    double bRK = x0*bH2O + x1*bCO2;

    double Zgas = gh_redlichKwong_Z(T, Pbar, bRK/(82.05*T), aRK/(bRK*82.05*pow(T,1.5)));
    double rh   = Pbar/(Zgas*R*T);

    /* mixture virial coefficients: mole-fraction-linear combination of the
       pure H2O/CO2 coefficients ("an initial guess at the functional
       form" per fluidPhase()'s own comment - ported as-is). c[] mixes;
       dc[] = c_H2O-c_CO2 is the x0-derivative of that linear mix, needed
       below for dG/dx_h2o.                                              */
    double c[nCOEFF], dc[nCOEFF];
    for (int i = 1; i < nCOEFF; i++){
        const double *rowH = h2o_c[i];
        const double *rowC = co2_c[i];
        double cH2O = rowH[1]/(T*T*T*T) + rowH[2]/(T*T) + rowH[3]/T + rowH[4] + rowH[5]*T + rowH[6]*T*T;
        double cCO2 = rowC[1]/(T*T*T*T) + rowC[2]/(T*T) + rowC[3]/T + rowC[4] + rowC[5]*T + rowC[6]*T*T;
        c[i]  = x0*cH2O + x1*cCO2;
        dc[i] = cH2O - cCO2;
    }

    double dp = DBL_MAX;
    for (int count = 1; count <= 200 && fabs(dp) > 10.0*DBL_EPSILON; count++){
        double temp1  = c[3] + 2.0*c[4]*rh + 3.0*c[5]*rh*rh + 4.0*c[6]*rh*rh*rh;
        double dtemp1 = 2.0*c[4] + 6.0*c[5]*rh + 12.0*c[6]*rh*rh;
        double temp2  = c[2] + c[3]*rh + c[4]*rh*rh + c[5]*rh*rh*rh + c[6]*rh*rh*rh*rh;
        double dtemp2 = c[3] + 2.0*c[4]*rh + 3.0*c[5]*rh*rh + 4.0*c[6]*rh*rh*rh;

        double pr  = rh + c[1]*rh*rh - rh*rh*temp1/(temp2*temp2)
                   + c[7]*rh*rh*exp(-c[8]*rh) + c[9]*rh*rh*exp(-c[10]*rh);
        pr *= R*T;
        pr -= Pbar;

        double dpr = 1.0 + 2.0*c[1]*rh - 2.0*rh*temp1/(temp2*temp2)
                   - rh*rh*(temp2*temp2*dtemp1 - temp1*2.0*temp2*dtemp2)/(temp2*temp2*temp2*temp2)
                   + 2.0*c[7]*rh*exp(-c[8]*rh) - c[7]*c[8]*rh*rh*exp(-c[8]*rh)
                   + 2.0*c[9]*rh*exp(-c[10]*rh) - c[9]*c[10]*rh*rh*exp(-c[10]*rh);
        dpr *= R*T;

        dp  = -pr/dpr;
        rh += dp;
        if (rh < 0.0) rh = 10.0*DBL_EPSILON;
    }

    double term1      = c[2] + c[3]*rh + c[4]*rh*rh + c[5]*rh*rh*rh + c[6]*rh*rh*rh*rh;
    double dterm1drh  = c[3] + 2.0*c[4]*rh + 3.0*c[5]*rh*rh + 4.0*c[6]*rh*rh*rh;
    double term1_x    = dc[2] + dc[3]*rh + dc[4]*rh*rh + dc[5]*rh*rh*rh + dc[6]*rh*rh*rh*rh;

    double term2a     = c[7]/c[8];
    double term2a_x   = dc[7]/c[8] - c[7]*dc[8]/(c[8]*c[8]);
    double term2b     = exp(-c[8]*rh) - 1.0;
    double dterm2bdrh = -c[8]*exp(-c[8]*rh);
    double term2b_x   = -dc[8]*rh*exp(-c[8]*rh);
    double term2      = term2a*term2b;
    double dterm2drh  = term2a*dterm2bdrh;
    double term2_x    = term2a_x*term2b + term2a*term2b_x;

    double term3a     = c[9]/c[10];
    double term3a_x   = dc[9]/c[10] - c[9]*dc[10]/(c[10]*c[10]);
    double term3b     = exp(-c[10]*rh) - 1.0;
    double dterm3bdrh = -c[10]*exp(-c[10]*rh);
    double term3b_x   = -dc[10]*rh*exp(-c[10]*rh);
    double term3      = term3a*term3b;
    double dterm3drh  = term3a*dterm3bdrh;
    double term3_x    = term3a_x*term3b + term3a*term3b_x;

    double Ar     = c[1]*rh + 1.0/term1 - 1.0/c[2] - term2 - term3;
    double dArdrh = c[1] - dterm1drh/(term1*term1) - dterm2drh - dterm3drh;
    double Ar_x   = dc[1]*rh - term1_x/(term1*term1) + dc[2]/(c[2]*c[2]) - term2_x - term3_x;

    double Ai_h2o, dAidrh_h2o, Ai_co2, dAidrh_co2;
    GH_ideal_gas_H2O(T, rh, &Ai_h2o, &dAidrh_h2o);
    GH_ideal_gas_CO2(T, rh, &Ai_co2, &dAidrh_co2);
    double Ai      = x0*Ai_h2o + x1*Ai_co2;
    double dAidrh  = x0*dAidrh_h2o + x1*dAidrh_co2;
    double Ai_x    = Ai_h2o - Ai_co2;

    double A      = R*T*Ar   + Ai;
    double dAdrh  = R*T*dArdrh + dAidrh;
    double A_x    = R*T*Ar_x + Ai_x;   /* dA/dx_h2o at fixed rho (envelope theorem) */

    double G      = (A + rh*dAdrh)/10.0; /* bar*cm3/mol -> J/mol */
    double dGdx   = A_x/10.0;

    /* reference-state shift, weighted by composition (reduces exactly to
       GH_pitzer_sterner_G's own shift at x_h2o=1 or x_h2o=0)             */
    const double refH2O_g = -228538.00;
    const double refCO2_g = -394341.00;
    G    += x0*(refH2O_g - (-46493.8016496949)) + x1*(refCO2_g - (-394450.0));
    dGdx += (refH2O_g - (-46493.8016496949)) - (refCO2_g - (-394450.0));

    *dGdx_h2o = dGdx;
    return G;
}

/* ============================================================================
   Haar (1984) H2O EOS at fixed p=1 bar, ported from water.c's whaar().
   Value + 1st-rho-derivative only (drops T-derivative and 2nd-rho-
   derivative bookkeeping only needed for real whaar()'s unused Cp/V
   outputs). See GH_fluid_eos.h for verification/provenance notes.
   ============================================================================ */
static void GH_kubik(double b, double c, double d, double *x1, double *x2, double *x2i, double *x3){
    double q, p, r, phi3, ff;
    const double pi = 3.14159263538979;
    *x2 = 0.0; *x2i = 0.0; *x3 = 0.0;
    if (c == 0.0 && d == 0.0) { *x1 = -b; return; }
    q = (2.0*b*b*b/27.0 - b*c/3.0 + d)/2.0;
    p = (3.0*c - b*b)/9.0;
    ff = fabs(p);
    r = sqrt(ff);
    ff = r*q;
    if (ff < 0.0) r = -r;
    ff = q/(r*r*r);
    if (p > 0.0) {
        phi3 = log(ff + sqrt(ff*ff+1.0))/3.0;
        *x1 = -r*(exp(phi3) - exp(-phi3)) - b/3.0;
        *x2i = 1.0;
    } else if (q*q + p*p*p > 0.0) {
        phi3 = log(ff + sqrt(ff*ff-1.0))/3.0;
        *x1 = -r*(exp(phi3) + exp(-phi3)) - b/3.0;
        *x2i = 1.0;
    } else {
        phi3 = atan(sqrt(1.0-ff*ff)/ff)/3.0;
        *x1 = -2.0*r*cos(phi3) - b/3.0;
        *x2 = 2.0*r*cos(pi/3.0-phi3) - b/3.0;
        *x2i = 0.0;
        *x3 = 2.0*r*cos(pi/3.0+phi3) - b/3.0;
    }
}

static double GH_psat2(double t){
    double w, wsq, v, ff;
    int i;
    static const double a[9] = { 0.0,
        -7.8889166, 2.5514255, -6.716169, 33.239495,
        -105.38479, 174.35319, -148.39348, 48.631602
    };
    if (t <= 314.0) return exp(6.3573118-8858.843/t + 607.56335/pow(t, 0.6));
    v = t/647.25;
    w = fabs(1.0-v);
    wsq = sqrt(w);
    ff = 0.0;
    for(i=1;i<=8;i++) { ff = ff + a[i]*w; w = w*wsq; }
    return 220.93*exp(ff/v);
}

double GH_haar_H2O_G(double t, double p){
    static const double gi[41] = { 0.0,
        -.53062968529023e4,  .22744901424408e5,  .78779333020687e4,
        -.69830527374994e3,  .17863832875422e6, -.39514731563338e6,
         .33803884280753e6, -.13855050202703e6, -.25637436613260e7,
         .48212575981415e7, -.34183016969660e7,  .12223156417448e7,
         .11797433655832e8, -.21734810110373e8,  .10829952168620e8,
        -.25441998064049e7, -.31377774947767e8,  .52911910757704e8,
        -.13802577177877e8, -.25109914369001e7,  .46561826115608e8,
        -.72752773275387e8,  .41774246148294e7,  .14016358244614e8,
        -.31555231392127e8,  .47929666384584e8,  .40912664781209e7,
        -.13626369388386e8,  .69625220862664e7, -.10834900096447e8,
        -.22722827401688e7,  .38365486000660e7,  .68833257944332e5,
         .21757245522644e6, -.26627944829770e5, -.70730418082074e6,
        -.225e1,            -1.68e1,             .055e1,
        -93.0e1
    };
    static const int ki[41] = { 0,
        1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5,
        6, 6, 6, 6, 7, 7, 7, 7, 9, 9, 9, 9, 3, 3, 1, 5, 2, 2, 2, 4
    };
    static const int li[41] = { 0,
        1, 2, 4, 6, 1, 2, 4, 6, 1, 2, 4, 6, 1, 2, 4, 6, 1, 2, 4, 6,
        1, 2, 4, 6, 1, 2, 4, 6, 1, 2, 4, 6, 0, 3, 3, 3, 0, 2, 0, 0
    };
    static const double rhoi[41] = { 0.0,
        0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0, 0.319, 0.310, 0.310, 1.55
    };
    static const double ttti[41] = { 0.0,
        0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0, 640.0, 640.0, 641.6, 270.0
    };
    static const double alpi[41] = { 0.0,
        0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0, 34.0, 40.0, 30.0, 1050.0
    };
    static const double beti[41] = { 0.0,
        0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0, 2.0e4, 2.0e4, 4.0e4, 25.0
    };
    static const double bi[6]  = { 0.7478629, -0.3540782, 0.0, 0.007159876, 0.0, -0.003528426 };
    static const double bbi[6] = { 1.1278334, -0.5944001, -5.010996, 0.0, 0.63684256, 0.0 };
    static const double ci[19] = { 0.0,
        .19730271018e2,      .209662681977e2,   -.483429455355,
        .605743189245e1,   22.56023885,        -9.87532442,
       -.43135538513e1,      .458155781,        -.47754901883e-1,
        .41238460633e-2,    -.27929052852e-3,    .14481695261e-4,
       -.56473658748e-6,     .16200446e-7,      -.3303822796e-9,
        .451916067368e-11,  -.370734122708e-13,  .137546068238e-15
    };

    const double r     = 4.6152;
    const double gref  = -54955.23970014;
    const double t0    = 647.073;
    const double rr    = 8.31441;
    const double alpha = 11.0;
    const double beta  = 133.0/3.0;
    const double gammaC= 7.0/2.0;
    const double P0    = 1.01325;
    /* p is now a caller-supplied parameter (was hardcoded to 1.0) - real
       gibbs.c calls whaar() at TWO different pressures depending on which
       standard state is being built: the liquid's own H2O basis species
       always uses the fixed 1-bar reference (matches the original port),
       while the standalone "water" pure phase (GH_gem_function.c's own
       "water" branch) calls whaar() at the ACTUAL pressure (capped at
       10000 bar) - see GH_wdh78_G below for the correction above that cap. */

    double taui[7];
    int i, count;

    taui[0] = 1.0; taui[1] = t/t0; for(i=2;i<=6;i++) taui[i] = taui[i-1]*taui[1];

    double b  = bi[1]*log(taui[1]) + bi[0] + bi[3]/taui[3] + bi[5]/taui[5];
    double bb = bbi[0] + bbi[1]/taui[1] + bbi[2]/taui[2] + bbi[4]/taui[4];

    double ps = 220.55;
    if (t <= 647.25) ps = GH_psat2(t);

    double ark = 1.279186e8 - 2.241415e4 * t;
    double brk = 1.428062e1 + 6.092237e-4 * t;
    double oft = ark/(p*sqrt(t));
    double buk = -10.0*rr*t/p;
    double cuk = oft - brk*brk + brk*buk;
    double duk = - brk*oft;
    double x1, x2, x2i, x3;
    GH_kubik(buk, cuk, duk, &x1, &x2, &x2i, &x3);

    double vol;
    if (x2i != 0.0) vol = x1;
    else            vol = (p < ps) ? fmax(x1, fmax(x2, x3)) : fmin(x1, fmin(x2, x3));

    double rhn = (vol <= 0.0) ? 1.9 : (1.0/vol)*18.0152;

    double dp = DBL_MAX, dr = DBL_MAX, rh = rhn;
    for (count=1; count<=100 && (dp > 10.0*DBL_EPSILON || dr > 10.0*DBL_EPSILON); count++){
        double y, ermi[10], prv, dpr;
        rh = rhn;
        if (rh <= 0.0) rh = 1.e-8;
        if (rh >  1.9) rh = 1.9;
        y = rh*b/4.0;
        ermi[0] = 1.0; ermi[1] = 1.0-exp(-rh);
        for (i=2;i<=9;i++) ermi[i] = ermi[i-1]*ermi[1];

        prv = 0.0; dpr = 0.0;
        for (i=1; i<=36; i++) {
            prv = prv + gi[i]/taui[li[i]]*ermi[ki[i]-1];
            dpr = dpr + (2.0+rh*(ki[i]*exp(-rh)-1.0)/ermi[1])*gi[i]/taui[li[i]]*ermi[ki[i]-1];
        }
        for(i=37; i<=40; i++) {
            double del, tau, abc, q10, qm;
            del = rh/rhoi[i] - 1.0;
            tau = t/ttti[i]  - 1.0;
            abc = -alpi[i] * pow(del, (double) ki[i]) - beti[i] * tau*tau;
            q10 = (abc > -100.00) ? gi[i] * pow(del, (double) li[i]) * exp(abc) : 0.0;
            qm = li[i]/del - ki[i]*alpi[i]*pow(del, (double) (ki[i]-1));
            prv = prv + q10*qm*rh*rh/rhoi[i];
            dpr = dpr + (q10*qm*rh*rh/rhoi[i]) * (2.0/rh+qm/rhoi[i])
                      - rh*rh/(rhoi[i]*rhoi[i])*q10*
                        (li[i]/del/del + ki[i]*(ki[i]-1)*alpi[i]*pow(del, (double) (ki[i]-2)));
        }

        prv = rh*(rh*exp(-rh)*prv + r*t*((1.0 + alpha*y + beta*y*y)/pow(1.0-y,3.0)
                                + 4.0*y*(bb/b - gammaC)));
        dpr = rh*exp(-rh)*dpr
            + r*t*( (1.0 + 2.0*alpha*y + 3.0*beta*y*y)/pow(1.0-y,3.0)
                   + 3.0*y*(1.0 + alpha*y + beta*y*y)/pow(1.0-y,4.0)
                   + 2.0*4.0*y*(bb/b - gammaC));

        if (dpr <= 0.0) rhn *= (p <= ps) ? 0.95 : 1.05;
        else {
            double x;
            if (dpr < 0.01) dpr = 0.01;
            x = (p - prv)/dpr;
            if (fabs(x) > 0.1) x = 0.1*x/fabs(x);
            rhn = rh + x;
        }
        dp = fabs(1.0 - prv/p);
        dr = fabs(1.0 - rhn/rh);
    }
    rh = rhn;

    /* base function Z, value + rho-derivative only */
    double y = rh*b/4.0;
    double dydrh = b/4.0;

    double Z = - log(1.0-y) - (beta-1.0)/(1.0-y)
             + (alpha+beta+1.0)/(2.0*(1.0-y)*(1.0-y)) + 4*y*(bb/b - gammaC)
             - (alpha-beta+3.0)/2.0 + log(rh*r*t/P0);
    double dZdrh = 1.0/rh + 4.0*(bb/b-gammaC)*dydrh + dydrh/(1.0-y)
             + (alpha+beta+1.0)*dydrh/pow(1.0-y,3.0)
             - (beta-1.0)*dydrh/pow(1.0-y,2.0);

    double Ab = r*t*Z;
    double dAbdrh = r*t*dZdrh;

    /* residual function Ar, value + rho-derivative only */
    double ermi[10];
    ermi[0] = 1.0; ermi[1] = 1.0-exp(-rh);
    for (i=2; i<=9; i++) ermi[i] = ermi[i-1]*ermi[1];

    double Ar = 0.0, dArdrh = 0.0;
    for (i=1; i<=36; i++){
        Ar     += gi[i]/ki[i]/taui[li[i]]*ermi[ki[i]];
        dArdrh += gi[i]/taui[li[i]]*ermi[ki[i]-1]*exp(-rh);
    }
    for(i=37; i<=40; i++){
        double del = rh/rhoi[i] - 1.0;
        double tau = t/ttti[i] - 1.0;
        double Q = -alpi[i]*pow(del, (double) ki[i]) - beti[i]*tau*tau;
        double dQdrh = (ki[i] == 0) ? 0.0 : -alpi[i]*ki[i]*pow(del, (double) (ki[i]-1))/rhoi[i];
        double expQ = (Q > -100.0) ? exp(Q) : 0.0;
        Ar += gi[i]*pow(del, (double) li[i])*expQ;
        dArdrh += (li[i] == 0) ? gi[i]*expQ*dQdrh
                : gi[i]*li[i]*pow(del, (double) (li[i]-1))*expQ/rhoi[i]
                  + gi[i]*pow(del, (double) li[i])*expQ*dQdrh;
    }

    /* ideal gas Ai, T-only value (no rho-dependence) */
    double tr = t/1.0e2;
    double Zi = 1.0 + (ci[1]/tr + ci[2])*log(tr);
    for (i=3; i<=18; i++) Zi += ci[i]*pow(tr, (double) (i-6));
    double Ai = -r*t*Zi;

    double A     = Ab + Ar + Ai;
    double dAdrh = dAbdrh + dArdrh;

    double pcalc = rh*rh*dAdrh;
    double gH2O  = A + pcalc/rh;

    gH2O *= 1.80152; /* -> J/mol */

    /* Berman shift, MODE_xMELTS/MODE__MELTS/MODE__MELTSandCO2/MODE__MELTSandCO2_H2O branch (the only one gh needs) */
    gH2O += -285829.96 - (298.15*69.9146) - gref;

    return gH2O;
}

/* ============================================================================
   wdh78 high-pressure correction for water.c's whaar() EOS, ported from
   water.c's own wdh78() ("Returns difference in thermodynamic properties
   of water between p,t and 10kb,t"). Real gibbs.c's standalone "water"
   pure-phase branch calls whaar() at min(P,10000) then, if the ACTUAL
   pressure exceeds 10000 bar, adds this delta to extend the EOS beyond
   whaar()'s own reliable range - found missing entirely from gh (which
   only had the liquid's fixed-1-bar H2O, never the standalone "water"
   phase) during the 2026-07-15 grid sweep. Value-only port (drops the
   h/s/cp/v/dvdt/etc: outputs real wdh78() also computes, none needed for
   a standard-state G).
   ============================================================================ */
static double GH_wdh78_poly(double t_celsius, double p_bar){
    static const double a[5][5] = {
        { -5.6130073e+04,  3.8101798e-01, -2.1167697e-06,  2.0266445e-11, -8.3225572e-17 },
        { -1.5285559e+01,  1.3752390e-04, -1.5586868e-09,  6.6329577e-15,  0.0           },
        { -2.6092451e-02,  3.5988857e-08, -2.7916588e-14,  0.0,            0.0           },
        {  1.7140501e-05, -1.6860893e-11,  0.0,            0.0,            0.0           },
        { -6.0126987e-09,  0.0,            0.0,            0.0,            0.0           }
    };
    double g = 0.0;
    for (int j = 0; j < 5; j++){
        for (int l = 0; l < 5-j; l++){
            g += a[j][l]*pow(t_celsius, (double) j)*pow(p_bar, (double) l);
        }
    }
    return g;
}

double GH_wdh78_G(double t, double p){
    double t_celsius = t - 273.15;
    double g10kb = GH_wdh78_poly(t_celsius, 10000.0);
    double g_p   = GH_wdh78_poly(t_celsius, p);
    return (g_p - g10kb) * 4.184; /* cal -> J, matching real wdh78()'s own *4.184 */
}

/* ============================================================================
   Duan (1992) pure-CO2 EOS at fixed p=1 bar, ported from fluidPhase.c's
   duanCO2Driver()+idealGasCO2(), specialized to the pure endpoint and
   value-only. See GH_fluid_eos.h for verification/provenance notes.
   ============================================================================ */
static const double GH_CO2Tc = 304.1282;
static const double GH_CO2Pc = 73.773;
#define GH_CO2Vc (8.314467*GH_CO2Tc/GH_CO2Pc)

static const double GH_CO2La1  =  1.14400435E-01;
static const double GH_CO2La2  = -9.38526684E-01;
static const double GH_CO2La3  =  7.21857006E-01;
static const double GH_CO2La4  =  8.81072902E-03;
static const double GH_CO2La5  =  6.36473911E-02;
static const double GH_CO2La6  = -7.70822213E-02;
static const double GH_CO2La7  =  9.01506064E-04;
static const double GH_CO2La8  = -6.81834166E-03;
static const double GH_CO2La9  =  7.32364258E-03;
static const double GH_CO2La10 = -1.10288237E-04;
static const double GH_CO2La11 =  1.26524193E-03;
static const double GH_CO2La12 = -1.49730823E-03;
static const double GH_CO2La   =  7.81940730E-03;
static const double GH_CO2Lb   = -4.22918013E+00;
static const double GH_CO2Lc   =  1.58500000E-01;

/* Berman(1988)-style ideal-gas Cp/H/S polynomial coefficients for CO2,
   ported from fluidPhase.c's idealCoeff[][CO2] column (index 1)          */
static const double GH_idealCoeff_CO2[13] = {
    -1.8188731,       12.903022,       -9.6634864,       4.2251879,
    -1.0421640,       0.12683515,      -0.0049939675,    2.4950242,
    -0.82723750,      0.15372481,      -0.015861243,     0.00086017150,
    -0.000019222165
};

double GH_duan_CO2_G(double t){
    const double R_duan = 8.314467;
    const double p = 1.0; /* fixed reference pressure - matches propertiesOfPureCO2(t, 1.0, ...) call site in real gibbs.c */
    double CO2Tr = t/GH_CO2Tc;

    double bEnd = GH_CO2La1 + GH_CO2La2/CO2Tr/CO2Tr + GH_CO2La3/CO2Tr/CO2Tr/CO2Tr;
    double cEnd = GH_CO2La4 + GH_CO2La5/CO2Tr/CO2Tr + GH_CO2La6/CO2Tr/CO2Tr/CO2Tr;
    double dEnd = GH_CO2La7 + GH_CO2La8/CO2Tr/CO2Tr + GH_CO2La9/CO2Tr/CO2Tr/CO2Tr;
    double eEnd = GH_CO2La10 + GH_CO2La11/CO2Tr/CO2Tr + GH_CO2La12/CO2Tr/CO2Tr/CO2Tr;
    double fEnd = GH_CO2La/CO2Tr/CO2Tr/CO2Tr;

    double bv = bEnd*GH_CO2Vc;
    double cv = cEnd*GH_CO2Vc*GH_CO2Vc;
    double dv = dEnd*GH_CO2Vc*GH_CO2Vc*GH_CO2Vc*GH_CO2Vc;
    double ev = eEnd*GH_CO2Vc*GH_CO2Vc*GH_CO2Vc*GH_CO2Vc*GH_CO2Vc;
    double fv = fEnd*GH_CO2Vc*GH_CO2Vc;
    double beta = GH_CO2Lb;
    double gammav = GH_CO2Lc*GH_CO2Vc*GH_CO2Vc;

    double bvPrime = 2.0*bv;
    double cvPrime = 3.0*cv;
    double dvPrime = 5.0*dv;
    double evPrime = 6.0*ev;
    double fvPrime = 2.0*fv;
    double betaPrime = beta;
    double gammavPrime = 3.0*gammav;

    double v, z;
    int iter = 0;
    double delv = 1.0, vPrevious = 1.0, delvPrevious = 1.0;
    v = R_duan*t/p;
    while (iter < 200) {
        z = 1.0 + bv/v + cv/v/v + dv/v/v/v/v + ev/v/v/v/v/v + (fv/v/v) * (beta + gammav/v/v) * exp(-gammav/v/v);
        delv = z*R_duan*t/p - v;
        if ( ((iter > 1) && (delv*delvPrevious < 0.0)) || (fabs(delv) < v*100.0*DBL_EPSILON) ) break;
        vPrevious = v;
        delvPrevious = delv;
        v = (z*R_duan*t/p + v)/2.0;
        iter++;
    }
    if (fabs(delv) > v*100.0*DBL_EPSILON && iter < 200) {
        double dx;
        double rtb = (delv < 0.0) ? (dx = vPrevious-v,v) : (dx = v-vPrevious,vPrevious);
        iter = 0;
        while (iter < 200) {
            v = rtb + (dx *= 0.5);
            z = 1.0 + bv/v + cv/v/v + dv/v/v/v/v + ev/v/v/v/v/v + (fv/v/v) * (beta + gammav/v/v) * exp(-gammav/v/v);
            delv = z*R_duan*t/p - v;
            if (delv <= 0.0) rtb = v;
            if ( (fabs(dx) < 100.0*DBL_EPSILON) || (delv == 0.0) ) break;
            iter++;
        }
    }

    double lnPhiCO2 = 0.0;
    lnPhiCO2 += -log(z);
    lnPhiCO2 += bvPrime/v;
    lnPhiCO2 += cvPrime/2.0/v/v;
    lnPhiCO2 += dvPrime/4.0/v/v/v/v;
    lnPhiCO2 += evPrime/5.0/v/v/v/v/v;
    lnPhiCO2 += ((fvPrime*beta + betaPrime*fv)/2.0/gammav)*(1.0-exp(-gammav/v/v));
    lnPhiCO2 += ((fvPrime*gammav+gammavPrime*fv-fv*beta*(gammavPrime-gammav))/2.0/gammav/gammav)
               *(1.0 - (gammav/v/v + 1.0)*exp(-gammav/v/v));
    lnPhiCO2 += ((gammavPrime-gammav)*fv/2.0/gammav/gammav)*(-2.0 + (gammav*gammav/v/v/v/v + 2.0*gammav/v/v + 2.0)*exp(-gammav/v/v));

    double phi = exp(lnPhiCO2);

    /* idealGasCO2: h0(t), s0(t) integrated Cp polynomial */
    double s0=0.0, h0=0.0;
    int i;
    for (i=0; i<7; i++) h0 += GH_idealCoeff_CO2[i]*pow(t/1000.0, (double) (i+1))/((double) (i+1));
    h0 += GH_idealCoeff_CO2[7]*log(t/1000.0);
    for (i=8; i<13; i++) h0 += GH_idealCoeff_CO2[i]/pow(t/1000.0, (double) (i-7))/((double) (7-i));
    s0 = GH_idealCoeff_CO2[0]*log(t/1000.0);
    for (i=1; i<7; i++) s0 += GH_idealCoeff_CO2[i]*pow(t/1000.0, (double) i)/((double) i);
    for (i=7; i<13; i++) s0 += GH_idealCoeff_CO2[i]/pow(t/1000.0, (double) (i-6))/((double) (6-i));
    h0 *= 8.31451*1000.0;
    s0 *= 8.31451;
    h0 += -385358.2260;
    s0 +=  210.0304;

    const double R_outer = 8.3143;
    double g = h0 - t*s0 + R_outer*t*log(phi*p);
    return g;
}

/* ============================================================================
   Duan & Zhang (2006) H2O-CO2 EOS - the ACTUAL rhyolite-MELTS "fluid"
   solution-phase model. Ported 2026-07-17 from real
   xMELTS-master/sources/fluidPhase.c (byte-identical in
   xMELTS-FinalH2O-CO2Model): the real "fluid" phase's function pointers
   (sol_struct_data.h: testFlu/conFlu/actFlu/gmixFlu/...) all resolve into
   fluidPhase.c, whose engine is duan()/duanH2O()/duanCO2() (DZ2006 virial
   EOS with composition-dependent mixing rules + fugacity coefficients),
   with components "h2oduan"/"co2duan" whose standard states are
   propertiesOfPureH2O/CO2() (ideal-gas h,s + RT*ln(phi*p)) - NOT
   fluid.c's fluidPhase() (Pitzer-Sterner + the crude linear c[i] mixing),
   which real MELTS only uses for the LIQUID's H2O/CO2 standard states via
   gibbs.c. gh's earlier "fl" used the Pitzer-Sterner model by mistake
   (whose huge, unphysical Gex produced a spurious H2O-CO2 two-fluid
   split); this port replaces it.

   Value + composition-derivative (Prime) slices only, transcribed
   VERBATIM from fluidPhase.c's BVc/CVc/DVc/EVc/FVc/Beta/GammaVc
   AndDerivative general-composition branches (the pure-endpoint branches
   there are algebraically identical special cases - verified against
   GH_duan_CO2_G's existing pure-CO2 shortcuts above), the duanDriver
   volume iteration + lnPhi expressions, and idealGasH2O/CO2. All
   T-derivative machinery (only needed for real MELTS' Cp/V outputs)
   dropped. Coefficient switch at p<=2000 bar between the low-P (La) and
   high-P (Ha) sets, exactly as real duan((p<=2000.0)?1:0, ...).

   Assembly (matching real MELTS exactly, R = 8.3143 J/K = silmin.h's R):
     g_pure_i  = h0_i(T) - T*s0_i(T) + R*T*ln(phi_pure_i * p)     [propertiesOfPure*]
     G_mix     = R*T*sum_i x_i*ln(x_i * phi_i^mix / phi_i^pure)   [gmixFlu]
     G_total   = sum_i x_i*g_pure_i + G_mix                        (J/mol)
     dG/dx_h2o = mu_H2O - mu_CO2 (Gibbs-Duhem-exact, actFlu's mu form)
   ============================================================================ */

/* H2O critical constants + DZ2006 coefficient tables (CO2's low-P set and
   Tc/Pc already exist above as GH_CO2La1..12 / GH_CO2Tc / GH_CO2Pc) */
static const double GH_H2OTc = 647.25;
static const double GH_H2OPc = 221.19;
#define GH_H2OVc (8.314467*GH_H2OTc/GH_H2OPc)

/* pure EOS terms, 0 to 0.2 GPa (fluidPhase.c "La" set) */
static const double GH_H2OLa1  =  4.38269941E-02;
static const double GH_H2OLa2  = -1.68244362E-01;
static const double GH_H2OLa3  = -2.36923373E-01;
static const double GH_H2OLa4  =  1.13027462E-02;
static const double GH_H2OLa5  = -7.67764181E-02;
static const double GH_H2OLa6  =  9.71820593E-02;
static const double GH_H2OLa7  =  6.62674916E-05;
static const double GH_H2OLa8  =  1.06637349E-03;
static const double GH_H2OLa9  = -1.23265258E-03;
static const double GH_H2OLa10 = -8.93953948E-06;
static const double GH_H2OLa11 = -3.88124606E-05;
static const double GH_H2OLa12 =  5.61510206E-05;
static const double GH_H2OLa   =  7.51274488E-03;
static const double GH_H2OLb   =  2.51598931E+00;
static const double GH_H2OLc   =  3.94000000E-02;

/* pure EOS terms, 0.2 to 10 GPa (fluidPhase.c "Ha" set) */
static const double GH_H2OHa1  =  4.68071541E-02;
static const double GH_H2OHa2  = -2.81275941E-01;
static const double GH_H2OHa3  = -2.43926365E-01;
static const double GH_H2OHa4  =  1.10016958E-02;
static const double GH_H2OHa5  = -3.86603525E-02;
static const double GH_H2OHa6  =  9.30095461E-02;
static const double GH_H2OHa7  = -1.15747171E-05;
static const double GH_H2OHa8  =  4.19873848E-04;
static const double GH_H2OHa9  = -5.82739501E-04;
static const double GH_H2OHa10 =  1.00936000E-06;
static const double GH_H2OHa11 = -1.01713593E-05;
static const double GH_H2OHa12 =  1.63934213E-05;
static const double GH_H2OHa   = -4.49505919E-02;
static const double GH_H2OHb   = -3.15028174E-01;
static const double GH_H2OHc   =  1.25000000E-02;

static const double GH_CO2Ha1  =  5.72573440E-03;
static const double GH_CO2Ha2  =  7.94836769E+00;
static const double GH_CO2Ha3  = -3.84236281E+01;
static const double GH_CO2Ha4  =  3.71600369E-02;
static const double GH_CO2Ha5  = -1.92888994E+00;
static const double GH_CO2Ha6  =  6.64254770E+00;
static const double GH_CO2Ha7  = -7.02203950E-06;
static const double GH_CO2Ha8  =  1.77093234E-02;
static const double GH_CO2Ha9  = -4.81892026E-02;
static const double GH_CO2Ha10 =  3.88344869E-06;
static const double GH_CO2Ha11 = -5.54833167E-04;
static const double GH_CO2Ha12 =  1.70489748E-03;
static const double GH_CO2Ha   = -4.13039220E-01;
static const double GH_CO2Hb   = -8.47988634E+00;
static const double GH_CO2Hc   =  2.80000000E-02;

/* idealCoeff[][H2O] column (index 0) - the CO2 column exists above */
static const double GH_idealCoeff_H2O[13] = {
     3.10409601236035e+01, -3.91422080460869e+01,  3.79695277233575e+01,
    -2.18374910952284e+01,  7.42251494566339e+00, -1.38178929609470e+00,
     1.08807067571454e-01, -1.20771176848589e+01,  3.39105078851732e+00,
    -5.84520979955060e-01,  5.89930846488082e-02, -3.12970001415882e-03,
     6.57460740981757e-05
};

#define iW 0    /* H2O index, real fluidPhase.c's "#define H2O 0" */
#define iC 1    /* CO2 index, real fluidPhase.c's "#define CO2 1" */

/* verbatim fluidPhase.c powSum() - signed-cube-root weighted mean */
static double GH_powSum(double a, double fa, double b, double fb) {
    double sum = 0.0;
    sum += (a >= 0.0) ? fa*pow(a, 1.0/3.0) : -fa*pow(-a, 1.0/3.0);
    sum += (b >= 0.0) ? fb*pow(b, 1.0/3.0) : -fb*pow(-b, 1.0/3.0);
    return (sum >= 0.0) ? pow(sum/(fa+fb), 3.0) : -pow(-sum/(fa+fb), 3.0);
}

/* value + Prime slices of BVc/CVc/DVc/EVc/FVc/Beta/GammaVc AndDerivative,
   general-composition branches, transcribed verbatim */
static void GH_dz_coeffs(int useLowPcoeff, double t, const double x[2],
        double *bv, double bvP[2], double *cv, double cvP[2],
        double *dv, double dvP[2], double *ev, double evP[2],
        double *fv, double fvP[2], double *beta, double betaP[2],
        double *gammav, double gammavP[2]){

    double H2OTr = t/GH_H2OTc;
    double CO2Tr = t/GH_CO2Tc;
    double bEnd[2], cEnd[2], dEnd[2], eEnd[2], fEnd[2], betaEnd[2], gammaEnd[2];
    double k1, k2, k3;

    if (useLowPcoeff) {
        bEnd[iW] = GH_H2OLa1 + GH_H2OLa2/H2OTr/H2OTr + GH_H2OLa3/H2OTr/H2OTr/H2OTr;
        bEnd[iC] = GH_CO2La1 + GH_CO2La2/CO2Tr/CO2Tr + GH_CO2La3/CO2Tr/CO2Tr/CO2Tr;
        k1 = 3.131 - 5.0624e-3*t + 1.8641e-6*t*t - 31.409/t;
        cEnd[iW] = GH_H2OLa4 + GH_H2OLa5/H2OTr/H2OTr + GH_H2OLa6/H2OTr/H2OTr/H2OTr;
        cEnd[iC] = GH_CO2La4 + GH_CO2La5/CO2Tr/CO2Tr + GH_CO2La6/CO2Tr/CO2Tr/CO2Tr;
        k2 = -46.646 + 4.2877e-2*t - 1.0892e-5*t*t + 1.5782e4/t;
        dEnd[iW] = GH_H2OLa7 + GH_H2OLa8/H2OTr/H2OTr + GH_H2OLa9/H2OTr/H2OTr/H2OTr;
        dEnd[iC] = GH_CO2La7 + GH_CO2La8/CO2Tr/CO2Tr + GH_CO2La9/CO2Tr/CO2Tr/CO2Tr;
        eEnd[iW] = GH_H2OLa10 + GH_H2OLa11/H2OTr/H2OTr + GH_H2OLa12/H2OTr/H2OTr/H2OTr;
        eEnd[iC] = GH_CO2La10 + GH_CO2La11/CO2Tr/CO2Tr + GH_CO2La12/CO2Tr/CO2Tr/CO2Tr;
        fEnd[iW] = GH_H2OLa/H2OTr/H2OTr/H2OTr;
        fEnd[iC] = GH_CO2La/CO2Tr/CO2Tr/CO2Tr;
        betaEnd[iW] = GH_H2OLb;
        betaEnd[iC] = GH_CO2Lb;
        gammaEnd[iW] = GH_H2OLc;
        gammaEnd[iC] = GH_CO2Lc;
        k3 = 0.9;
    } else {
        bEnd[iW] = GH_H2OHa1 + GH_H2OHa2/H2OTr/H2OTr + GH_H2OHa3/H2OTr/H2OTr/H2OTr;
        bEnd[iC] = GH_CO2Ha1 + GH_CO2Ha2/CO2Tr/CO2Tr + GH_CO2Ha3/CO2Tr/CO2Tr/CO2Tr;
        k1 = 9.034 - 7.9212e-3*t + 2.3285e-6*t*t - 2.4221e3/t;
        cEnd[iW] = GH_H2OHa4 + GH_H2OHa5/H2OTr/H2OTr + GH_H2OHa6/H2OTr/H2OTr/H2OTr;
        cEnd[iC] = GH_CO2Ha4 + GH_CO2Ha5/CO2Tr/CO2Tr + GH_CO2Ha6/CO2Tr/CO2Tr/CO2Tr;
        k2 = -1.068 + 1.8756e-3*t - 4.9371e-7*t*t + 6.6180e2/t;
        dEnd[iW] = GH_H2OHa7 + GH_H2OHa8/H2OTr/H2OTr + GH_H2OHa9/H2OTr/H2OTr/H2OTr;
        dEnd[iC] = GH_CO2Ha7 + GH_CO2Ha8/CO2Tr/CO2Tr + GH_CO2Ha9/CO2Tr/CO2Tr/CO2Tr;
        eEnd[iW] = GH_H2OHa10 + GH_H2OHa11/H2OTr/H2OTr + GH_H2OHa12/H2OTr/H2OTr/H2OTr;
        eEnd[iC] = GH_CO2Ha10 + GH_CO2Ha11/CO2Tr/CO2Tr + GH_CO2Ha12/CO2Tr/CO2Tr/CO2Tr;
        fEnd[iW] = GH_H2OHa/H2OTr/H2OTr/H2OTr;
        fEnd[iC] = GH_CO2Ha/CO2Tr/CO2Tr/CO2Tr;
        betaEnd[iW] = GH_H2OHb;
        betaEnd[iC] = GH_CO2Hb;
        gammaEnd[iW] = GH_H2OHc;
        gammaEnd[iC] = GH_CO2Hc;
        k3 = 1.0;
    }

    /* BVc general-x branch */
    *bv  = 0.0;
    *bv += bEnd[iW]*GH_H2OVc*x[iW]*x[iW];
    *bv += 2.0*GH_powSum(bEnd[iW], 1.0, bEnd[iC], 1.0)*k1*GH_powSum(GH_H2OVc, 1.0, GH_CO2Vc, 1.0)*x[iW]*x[iC];
    *bv += bEnd[iC]*GH_CO2Vc*x[iC]*x[iC];

    bvP[iW]  = 0.0;
    bvP[iW] += 2.0*bEnd[iW]*GH_H2OVc*x[iW];
    bvP[iW] += 2.0*GH_powSum(bEnd[iW], 1.0, bEnd[iC], 1.0)*k1*GH_powSum(GH_H2OVc, 1.0, GH_CO2Vc, 1.0)*x[iC];

    bvP[iC]  = 0.0;
    bvP[iC] += 2.0*GH_powSum(bEnd[iW], 1.0, bEnd[iC], 1.0)*k1*GH_powSum(GH_H2OVc, 1.0, GH_CO2Vc, 1.0)*x[iW];
    bvP[iC] += 2.0*bEnd[iC]*GH_CO2Vc*x[iC];

    /* CVc general-x branch */
    *cv  = 0.0;
    *cv += cEnd[iW]*GH_H2OVc*GH_H2OVc*x[iW]*x[iW]*x[iW];
    *cv += 3.0*GH_powSum(cEnd[iW], 2.0, cEnd[iC], 1.0)*k2*pow( GH_powSum(GH_H2OVc, 2.0, GH_CO2Vc, 1.0), 2.0)*x[iW]*x[iW]*x[iC];
    *cv += 3.0*GH_powSum(cEnd[iW], 1.0, cEnd[iC], 2.0)*k2*pow( GH_powSum(GH_H2OVc, 1.0, GH_CO2Vc, 2.0), 2.0)*x[iW]*x[iC]*x[iC];
    *cv += cEnd[iC]*GH_CO2Vc*GH_CO2Vc*x[iC]*x[iC]*x[iC];

    cvP[iW]  = 0.0;
    cvP[iW] += 3.0*cEnd[iW]*GH_H2OVc*GH_H2OVc*x[iW]*x[iW];
    cvP[iW] += 6.0*GH_powSum(cEnd[iW], 2.0, cEnd[iC], 1.0)*k2*pow(GH_powSum(GH_H2OVc, 2.0, GH_CO2Vc, 1.0), 2.0)*x[iW]*x[iC];
    cvP[iW] += 3.0*GH_powSum(cEnd[iW], 1.0, cEnd[iC], 2.0)*k2*pow(GH_powSum(GH_H2OVc, 1.0, GH_CO2Vc, 2.0), 2.0)*x[iC]*x[iC];

    cvP[iC]  = 0.0;
    cvP[iC] += 3.0*GH_powSum(cEnd[iW], 2.0, cEnd[iC], 1.0)*k2*pow(GH_powSum(GH_H2OVc, 2.0, GH_CO2Vc, 1.0), 2.0)*x[iW]*x[iW];
    cvP[iC] += 6.0*GH_powSum(cEnd[iW], 1.0, cEnd[iC], 2.0)*k2*pow(GH_powSum(GH_H2OVc, 1.0, GH_CO2Vc, 2.0), 2.0)*x[iW]*x[iC];
    cvP[iC] += 3.0*cEnd[iC]*GH_CO2Vc*GH_CO2Vc*x[iC]*x[iC];

    /* DVc general-x branch (no k factor) */
    *dv  = 0.0;
    *dv += dEnd[iW]*GH_H2OVc*GH_H2OVc*GH_H2OVc*GH_H2OVc*x[iW]*x[iW]*x[iW]*x[iW]*x[iW];
    *dv +=  5.0*GH_powSum(dEnd[iW], 4.0, dEnd[iC], 1.0)*pow(GH_powSum(GH_H2OVc, 4.0, GH_CO2Vc, 1.0), 4.0)*x[iW]*x[iW]*x[iW]*x[iW]*x[iC];
    *dv += 10.0*GH_powSum(dEnd[iW], 3.0, dEnd[iC], 2.0)*pow(GH_powSum(GH_H2OVc, 3.0, GH_CO2Vc, 2.0), 4.0)*x[iW]*x[iW]*x[iW]*x[iC]*x[iC];
    *dv += 10.0*GH_powSum(dEnd[iW], 2.0, dEnd[iC], 3.0)*pow(GH_powSum(GH_H2OVc, 2.0, GH_CO2Vc, 3.0), 4.0)*x[iW]*x[iW]*x[iC]*x[iC]*x[iC];
    *dv +=  5.0*GH_powSum(dEnd[iW], 1.0, dEnd[iC], 4.0)*pow(GH_powSum(GH_H2OVc, 1.0, GH_CO2Vc, 4.0), 4.0)*x[iW]*x[iC]*x[iC]*x[iC]*x[iC];
    *dv += dEnd[iC]*GH_CO2Vc*GH_CO2Vc*GH_CO2Vc*GH_CO2Vc*x[iC]*x[iC]*x[iC]*x[iC]*x[iC];

    dvP[iW]  =  0.0;
    dvP[iW] +=  5.0*dEnd[iW]*GH_H2OVc*GH_H2OVc*GH_H2OVc*GH_H2OVc*x[iW]*x[iW]*x[iW]*x[iW];
    dvP[iW] += 20.0*GH_powSum(dEnd[iW], 4.0, dEnd[iC], 1.0)*pow(GH_powSum(GH_H2OVc, 4.0, GH_CO2Vc, 1.0), 4.0)*x[iW]*x[iW]*x[iW]*x[iC];
    dvP[iW] += 30.0*GH_powSum(dEnd[iW], 3.0, dEnd[iC], 2.0)*pow(GH_powSum(GH_H2OVc, 3.0, GH_CO2Vc, 2.0), 4.0)*x[iW]*x[iW]*x[iC]*x[iC];
    dvP[iW] += 20.0*GH_powSum(dEnd[iW], 2.0, dEnd[iC], 3.0)*pow(GH_powSum(GH_H2OVc, 2.0, GH_CO2Vc, 3.0), 4.0)*x[iW]*x[iC]*x[iC]*x[iC];
    dvP[iW] +=  5.0*GH_powSum(dEnd[iW], 1.0, dEnd[iC], 4.0)*pow(GH_powSum(GH_H2OVc, 1.0, GH_CO2Vc, 4.0), 4.0)*x[iC]*x[iC]*x[iC]*x[iC];

    dvP[iC]  =  0.0;
    dvP[iC] +=  5.0*GH_powSum(dEnd[iW], 4.0, dEnd[iC], 1.0)*pow(GH_powSum(GH_H2OVc, 4.0, GH_CO2Vc, 1.0), 4.0)*x[iW]*x[iW]*x[iW]*x[iW];
    dvP[iC] += 20.0*GH_powSum(dEnd[iW], 3.0, dEnd[iC], 2.0)*pow(GH_powSum(GH_H2OVc, 3.0, GH_CO2Vc, 2.0), 4.0)*x[iW]*x[iW]*x[iW]*x[iC];
    dvP[iC] += 30.0*GH_powSum(dEnd[iW], 2.0, dEnd[iC], 3.0)*pow(GH_powSum(GH_H2OVc, 2.0, GH_CO2Vc, 3.0), 4.0)*x[iW]*x[iW]*x[iC]*x[iC];
    dvP[iC] += 20.0*GH_powSum(dEnd[iW], 1.0, dEnd[iC], 4.0)*pow(GH_powSum(GH_H2OVc, 1.0, GH_CO2Vc, 4.0), 4.0)*x[iW]*x[iC]*x[iC]*x[iC];
    dvP[iC] +=  5.0*dEnd[iC]*GH_CO2Vc*GH_CO2Vc*GH_CO2Vc*GH_CO2Vc*x[iC]*x[iC]*x[iC]*x[iC];

    /* EVc general-x branch (no k factor) */
    *ev  = 0.0;
    *ev += eEnd[iW]*GH_H2OVc*GH_H2OVc*GH_H2OVc*GH_H2OVc*GH_H2OVc*x[iW]*x[iW]*x[iW]*x[iW]*x[iW]*x[iW];
    *ev +=  6.0*GH_powSum(eEnd[iW], 5.0, eEnd[iC], 1.0)*pow(GH_powSum(GH_H2OVc, 5.0, GH_CO2Vc, 1.0), 5.0)*x[iW]*x[iW]*x[iW]*x[iW]*x[iW]*x[iC];
    *ev += 15.0*GH_powSum(eEnd[iW], 4.0, eEnd[iC], 2.0)*pow(GH_powSum(GH_H2OVc, 4.0, GH_CO2Vc, 2.0), 5.0)*x[iW]*x[iW]*x[iW]*x[iW]*x[iC]*x[iC];
    *ev += 20.0*GH_powSum(eEnd[iW], 3.0, eEnd[iC], 3.0)*pow(GH_powSum(GH_H2OVc, 3.0, GH_CO2Vc, 3.0), 5.0)*x[iW]*x[iW]*x[iW]*x[iC]*x[iC]*x[iC];
    *ev += 15.0*GH_powSum(eEnd[iW], 2.0, eEnd[iC], 4.0)*pow(GH_powSum(GH_H2OVc, 2.0, GH_CO2Vc, 4.0), 5.0)*x[iW]*x[iW]*x[iC]*x[iC]*x[iC]*x[iC];
    *ev +=  6.0*GH_powSum(eEnd[iW], 1.0, eEnd[iC], 5.0)*pow(GH_powSum(GH_H2OVc, 1.0, GH_CO2Vc, 5.0), 5.0)*x[iW]*x[iC]*x[iC]*x[iC]*x[iC]*x[iC];
    *ev += eEnd[iC]*GH_CO2Vc*GH_CO2Vc*GH_CO2Vc*GH_CO2Vc*GH_CO2Vc*x[iC]*x[iC]*x[iC]*x[iC]*x[iC]*x[iC];

    evP[iW]  =  0.0;
    evP[iW] +=  6.0*eEnd[iW]*GH_H2OVc*GH_H2OVc*GH_H2OVc*GH_H2OVc*GH_H2OVc*x[iW]*x[iW]*x[iW]*x[iW]*x[iW];
    evP[iW] += 30.0*GH_powSum(eEnd[iW], 5.0, eEnd[iC], 1.0)*pow(GH_powSum(GH_H2OVc, 5.0, GH_CO2Vc, 1.0), 5.0)*x[iW]*x[iW]*x[iW]*x[iW]*x[iC];
    evP[iW] += 60.0*GH_powSum(eEnd[iW], 4.0, eEnd[iC], 2.0)*pow(GH_powSum(GH_H2OVc, 4.0, GH_CO2Vc, 2.0), 5.0)*x[iW]*x[iW]*x[iW]*x[iC]*x[iC];
    evP[iW] += 60.0*GH_powSum(eEnd[iW], 3.0, eEnd[iC], 3.0)*pow(GH_powSum(GH_H2OVc, 3.0, GH_CO2Vc, 3.0), 5.0)*x[iW]*x[iW]*x[iC]*x[iC]*x[iC];
    evP[iW] += 30.0*GH_powSum(eEnd[iW], 2.0, eEnd[iC], 4.0)*pow(GH_powSum(GH_H2OVc, 2.0, GH_CO2Vc, 4.0), 5.0)*x[iW]*x[iC]*x[iC]*x[iC]*x[iC];
    evP[iW] +=  6.0*GH_powSum(eEnd[iW], 1.0, eEnd[iC], 5.0)*pow(GH_powSum(GH_H2OVc, 1.0, GH_CO2Vc, 5.0), 5.0)*x[iC]*x[iC]*x[iC]*x[iC]*x[iC];

    evP[iC]  =  0.0;
    evP[iC] +=  6.0*GH_powSum(eEnd[iW], 5.0, eEnd[iC], 1.0)*pow(GH_powSum(GH_H2OVc, 5.0, GH_CO2Vc, 1.0), 5.0)*x[iW]*x[iW]*x[iW]*x[iW]*x[iW];
    evP[iC] += 30.0*GH_powSum(eEnd[iW], 4.0, eEnd[iC], 2.0)*pow(GH_powSum(GH_H2OVc, 4.0, GH_CO2Vc, 2.0), 5.0)*x[iW]*x[iW]*x[iW]*x[iW]*x[iC];
    evP[iC] += 60.0*GH_powSum(eEnd[iW], 3.0, eEnd[iC], 3.0)*pow(GH_powSum(GH_H2OVc, 3.0, GH_CO2Vc, 3.0), 5.0)*x[iW]*x[iW]*x[iW]*x[iC]*x[iC];
    evP[iC] += 60.0*GH_powSum(eEnd[iW], 2.0, eEnd[iC], 4.0)*pow(GH_powSum(GH_H2OVc, 2.0, GH_CO2Vc, 4.0), 5.0)*x[iW]*x[iW]*x[iC]*x[iC]*x[iC];
    evP[iC] += 30.0*GH_powSum(eEnd[iW], 1.0, eEnd[iC], 5.0)*pow(GH_powSum(GH_H2OVc, 1.0, GH_CO2Vc, 5.0), 5.0)*x[iW]*x[iC]*x[iC]*x[iC]*x[iC];
    evP[iC] +=  6.0*eEnd[iC]*GH_CO2Vc*GH_CO2Vc*GH_CO2Vc*GH_CO2Vc*GH_CO2Vc*x[iC]*x[iC]*x[iC]*x[iC]*x[iC];

    /* FVc general-x branch (no k factor) */
    *fv  = 0.0;
    *fv += fEnd[iW]*GH_H2OVc*GH_H2OVc*x[iW]*x[iW];
    *fv += 2.0*GH_powSum(fEnd[iW], 1.0, fEnd[iC], 1.0)*pow(GH_powSum(GH_H2OVc, 1.0, GH_CO2Vc, 1.0), 2.0)*x[iW]*x[iC];
    *fv += fEnd[iC]*GH_CO2Vc*GH_CO2Vc*x[iC]*x[iC];

    fvP[iW]  = 0.0;
    fvP[iW] += 2.0*fEnd[iW]*GH_H2OVc*GH_H2OVc*x[iW];
    fvP[iW] += 2.0*GH_powSum(fEnd[iW], 1.0, fEnd[iC], 1.0)*pow(GH_powSum(GH_H2OVc, 1.0, GH_CO2Vc, 1.0), 2.0)*x[iC];

    fvP[iC]  = 0.0;
    fvP[iC] += 2.0*GH_powSum(fEnd[iW], 1.0, fEnd[iC], 1.0)*pow(GH_powSum(GH_H2OVc, 1.0, GH_CO2Vc, 1.0), 2.0)*x[iW];
    fvP[iC] += 2.0*fEnd[iC]*GH_CO2Vc*GH_CO2Vc*x[iC];

    /* Beta (linear mixing, no endpoint branching beyond the table switch) */
    *beta  = betaEnd[iW]*x[iW] + betaEnd[iC]*x[iC];
    betaP[iW] = betaEnd[iW];
    betaP[iC] = betaEnd[iC];

    /* GammaVc general-x branch */
    *gammav  = 0.0;
    *gammav += gammaEnd[iW]*GH_H2OVc*GH_H2OVc*x[iW]*x[iW]*x[iW];
    *gammav += 3.0*GH_powSum(gammaEnd[iW], 2.0, gammaEnd[iC], 1.0)*k3*pow(GH_powSum(GH_H2OVc, 2.0, GH_CO2Vc, 1.0), 2.0)*x[iW]*x[iW]*x[iC];
    *gammav += 3.0*GH_powSum(gammaEnd[iW], 1.0, gammaEnd[iC], 2.0)*k3*pow(GH_powSum(GH_H2OVc, 1.0, GH_CO2Vc, 2.0), 2.0)*x[iW]*x[iC]*x[iC];
    *gammav += gammaEnd[iC]*GH_CO2Vc*GH_CO2Vc*x[iC]*x[iC]*x[iC];

    gammavP[iW]  = 0.0;
    gammavP[iW] += 3.0*gammaEnd[iW]*GH_H2OVc*GH_H2OVc*x[iW]*x[iW];
    gammavP[iW] += 6.0*GH_powSum(gammaEnd[iW], 2.0, gammaEnd[iC], 1.0)*k3*pow(GH_powSum(GH_H2OVc, 2.0, GH_CO2Vc, 1.0), 2.0)*x[iW]*x[iC];
    gammavP[iW] += 3.0*GH_powSum(gammaEnd[iW], 1.0, gammaEnd[iC], 2.0)*k3*pow(GH_powSum(GH_H2OVc, 1.0, GH_CO2Vc, 2.0), 2.0)*x[iC]*x[iC];

    gammavP[iC]  = 0.0;
    gammavP[iC] += 3.0*GH_powSum(gammaEnd[iW], 2.0, gammaEnd[iC], 1.0)*k3*pow(GH_powSum(GH_H2OVc, 2.0, GH_CO2Vc, 1.0), 2.0)*x[iW]*x[iW];
    gammavP[iC] += 6.0*GH_powSum(gammaEnd[iW], 1.0, gammaEnd[iC], 2.0)*k3*pow(GH_powSum(GH_H2OVc, 1.0, GH_CO2Vc, 2.0), 2.0)*x[iW]*x[iC];
    gammavP[iC] += 3.0*gammaEnd[iC]*GH_CO2Vc*GH_CO2Vc*x[iC]*x[iC];
}

/* duanDriver's volume iteration + lnPhi expressions, verbatim (value only) */
static void GH_dz_phi_driver(int useLowPcoeff, double t, double p, const double x[2],
                      double *lnPhiH2O, double *lnPhiCO2){
    double bv, cv, dv, ev, fv, beta, gammav;
    double bvP[2], cvP[2], dvP[2], evP[2], fvP[2], betaP[2], gammavP[2];
    double v, z;

    GH_dz_coeffs(useLowPcoeff, t, x, &bv, bvP, &cv, cvP, &dv, dvP, &ev, evP, &fv, fvP, &beta, betaP, &gammav, gammavP);

    {
        int iter = 0;
        double delv = 1.0, vPrevious = 1.0, delvPrevious = 1.0;

        v = 8.314467*t/p;
        while (iter < 200) {
            z = 1.0 + bv/v + cv/v/v + dv/v/v/v/v + ev/v/v/v/v/v + (fv/v/v) * (beta + gammav/v/v) * exp(-gammav/v/v);
            delv = z*8.314467*t/p - v;
            if ( ((iter > 1) && (delv*delvPrevious < 0.0)) || (fabs(delv) < v*100.0*DBL_EPSILON) ) break;
            vPrevious = v;
            delvPrevious = delv;
            v = (z*8.314467*t/p + v)/2.0;
            iter++;
        }
        if (iter < 200 && fabs(delv) > v*100.0*DBL_EPSILON) {
            double dx;
            double rtb = (delv < 0.0) ? (dx = vPrevious-v,v) : (dx = v-vPrevious,vPrevious);
            iter = 0;
            while (iter < 200) {
                v = rtb + (dx *= 0.5);
                z = 1.0 + bv/v + cv/v/v + dv/v/v/v/v + ev/v/v/v/v/v + (fv/v/v) * (beta + gammav/v/v) * exp(-gammav/v/v);
                delv = z*8.314467*t/p - v;
                if (delv <= 0.0) rtb = v;
                if ( (fabs(dx) < 100.0*DBL_EPSILON) || (delv == 0.0) ) break;
                iter++;
            }
        }
    }

    *lnPhiH2O  = 0.0;
    *lnPhiH2O += -log(z);
    *lnPhiH2O += bvP[iW]/v;
    *lnPhiH2O += cvP[iW]/2.0/v/v;
    *lnPhiH2O += dvP[iW]/4.0/v/v/v/v;
    *lnPhiH2O += evP[iW]/5.0/v/v/v/v/v;
    *lnPhiH2O += ((fvP[iW]*beta + betaP[iW]*fv)/2.0/gammav)*(1.0-exp(-gammav/v/v));
    *lnPhiH2O += ((fvP[iW]*gammav+gammavP[iW]*fv-fv*beta*(gammavP[iW]-gammav))/2.0/gammav/gammav)
                *(1.0 - (gammav/v/v + 1.0)*exp(-gammav/v/v));
    *lnPhiH2O += ((gammavP[iW]-gammav)*fv/2.0/gammav/gammav)*(-2.0 + (gammav*gammav/v/v/v/v + 2.0*gammav/v/v + 2.0)*exp(-gammav/v/v));

    *lnPhiCO2  = 0.0;
    *lnPhiCO2 += -log(z);
    *lnPhiCO2 += bvP[iC]/v;
    *lnPhiCO2 += cvP[iC]/2.0/v/v;
    *lnPhiCO2 += dvP[iC]/4.0/v/v/v/v;
    *lnPhiCO2 += evP[iC]/5.0/v/v/v/v/v;
    *lnPhiCO2 += ((fvP[iC]*beta + betaP[iC]*fv)/2.0/gammav)*(1.0-exp(-gammav/v/v));
    *lnPhiCO2 += ((fvP[iC]*gammav+gammavP[iC]*fv-fv*beta*(gammavP[iC]-gammav))/2.0/gammav/gammav)
                *(1.0 - (gammav/v/v + 1.0)*exp(-gammav/v/v));
    *lnPhiCO2 += ((gammavP[iC]-gammav)*fv/2.0/gammav/gammav)*(-2.0 + (gammav*gammav/v/v/v/v + 2.0*gammav/v/v + 2.0)*exp(-gammav/v/v));
}

/* real duan()/duanH2O()/duanCO2() wrapper behavior: above 2000 bar the
   high-P (Ha) coefficient set is SPLICED onto the low-P (La) one for
   continuity - lnPhi(p>2000) = lnPhi_Ha(p) + lnPhi_La(2000) - lnPhi_Ha(2000),
   evaluated at the same composition (fluidPhase.c lines 1788-1810/1968-1982).
   Without this the two coefficient sets disagree at 2000 bar by a
   T-dependent constant (measured up to ~9.3 kJ at 873 K) - found via the
   first benchmark pass, which was exact below 2000 bar and offset by
   exactly that constant above. */
static void GH_dz_phi(double t, double p, const double x[2],
                      double *lnPhiH2O, double *lnPhiCO2){
    int useLow = (p <= 2000.0) ? 1 : 0;
    GH_dz_phi_driver(useLow, t, p, x, lnPhiH2O, lnPhiCO2);
    if (!useLow){
        double loW, loC, hiW, hiC;
        GH_dz_phi_driver(1, t, 2000.0, x, &loW, &loC);
        GH_dz_phi_driver(0, t, 2000.0, x, &hiW, &hiC);
        *lnPhiH2O += loW - hiW;
        *lnPhiCO2 += loC - hiC;
    }
}

/* idealGasH2O/idealGasCO2 h0(t), s0(t) (values only), verbatim */
static void GH_dz_ideal_hs(int is_h2o, double t, double *h0, double *s0){
    const double *ic = is_h2o ? GH_idealCoeff_H2O : GH_idealCoeff_CO2;
    int i;
    *h0 = 0.0;
    for (i=0; i<7; i++) *h0 += ic[i]*pow(t/1000.0, (double) (i+1))/((double) (i+1));
    *h0 += ic[7]*log(t/1000.0);
    for (i=8; i<13; i++) *h0 += ic[i]/pow(t/1000.0, (double) (i-7))/((double) (7-i));
    *s0 = ic[0]*log(t/1000.0);
    for (i=1; i<7; i++)  *s0 += ic[i]*pow(t/1000.0, (double) i)/((double) i);
    for (i=7; i<13; i++) *s0 += ic[i]/pow(t/1000.0, (double) (i-6))/((double) (6-i));
    *h0 *= 8.31451*1000.0;
    *s0 *= 8.31451;
    if (is_h2o){ *h0 += -355665.4136; *s0 +=  359.6505; }
    else       { *h0 += -385358.2260; *s0 +=  210.0304; }
}

/* == real propertiesOfPureH2O/CO2's g output (J/mol): ideal-gas h-T*s +
   R*T*ln(phi_pure*p), phi_pure from the DZ2006 EOS at the pure endpoint */
double GH_duan_pure_G(int is_h2o, double t, double p){
    const double R_ = 8.3143;
    double x[2] = { is_h2o ? 1.0 : 0.0, is_h2o ? 0.0 : 1.0 };
    double lnPhiW, lnPhiC, h0, s0;
    GH_dz_phi(t, p, x, &lnPhiW, &lnPhiC);
    GH_dz_ideal_hs(is_h2o, t, &h0, &s0);
    return h0 - t*s0 + R_*t*log(exp(is_h2o ? lnPhiW : lnPhiC)*p);
}

/* Full DZ2006 fluid G(x) in J/mol (x_h2o = mole fraction H2O), plus
   dG/dx_h2o = mu_H2O - mu_CO2 (Gibbs-Duhem-exact). Endpoint guards match
   real gmixFlu's 100*DBL_EPSILON zero-mixing cutoffs. */
double GH_duan_mix_G(double x_h2o, double t, double p, double *dGdx_h2o){
    const double R_ = 8.3143;
    double gW = GH_duan_pure_G(1, t, p);
    double gC = GH_duan_pure_G(0, t, p);
    double x[2] = { x_h2o, 1.0 - x_h2o };

    if ((fabs(x[iC]) < 100.0*DBL_EPSILON) || (fabs(x[iW]) < 100.0*DBL_EPSILON)){
        if (dGdx_h2o != (double *) 0) *dGdx_h2o = gW - gC;
        return x[iW]*gW + x[iC]*gC;
    }

    double lnPhiW, lnPhiC, lnPhiWpure, lnPhiCpure, h0, s0, dum;
    GH_dz_phi(t, p, x, &lnPhiW, &lnPhiC);
    {
        double xp[2] = {1.0, 0.0};
        GH_dz_phi(t, p, xp, &lnPhiWpure, &dum);
    }
    {
        double xp[2] = {0.0, 1.0};
        GH_dz_phi(t, p, xp, &dum, &lnPhiCpure);
    }

    double muW = gW + R_*t*(log(x[iW]) + lnPhiW - lnPhiWpure);
    double muC = gC + R_*t*(log(x[iC]) + lnPhiC - lnPhiCpure);

    if (dGdx_h2o != (double *) 0) *dGdx_h2o = muW - muC;
    return x[iW]*muW + x[iC]*muC;
}

/* Mixing-only chemical potentials of the DZ2006 fluid (J/mol):
   muGex_i = R*T*ln(x_i * phi_i^mix / phi_i^pure)  - actFlu's mu[] exactly
   (fluidPhase.c line 2486-2487). The pure-component standard states
   GH_duan_pure_G belong in gbase[] (so the LP/NLopt hyperplane rotation
   applies to them); these mixing parts are what obj_gh_fluid adds on top.
   Endpoint guards match real gmixFlu/actFlu's 100*DBL_EPSILON cutoffs. */
void GH_duan_mix_muGex(double x_h2o, double t, double p, double *muW, double *muC){
    const double R_ = 8.3143;
    double x[2] = { x_h2o, 1.0 - x_h2o };

    if ((fabs(x[iC]) < 100.0*DBL_EPSILON) || (fabs(x[iW]) < 100.0*DBL_EPSILON)){
        *muW = 0.0; *muC = 0.0;
        return;
    }

    double lnPhiW, lnPhiC, lnPhiWpure, lnPhiCpure, dum;
    GH_dz_phi(t, p, x, &lnPhiW, &lnPhiC);
    {
        double xp[2] = {1.0, 0.0};
        GH_dz_phi(t, p, xp, &lnPhiWpure, &dum);
    }
    {
        double xp[2] = {0.0, 1.0};
        GH_dz_phi(t, p, xp, &dum, &lnPhiCpure);
    }
    *muW = R_*t*(log(x[iW]) + lnPhiW - lnPhiWpure);
    *muC = R_*t*(log(x[iC]) + lnPhiC - lnPhiCpure);
}
