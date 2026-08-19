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
    Closed-form standard-state Gibbs energy for Ghiorso/MELTS liquid basis
    species (research group "gh"). Unlike the THERMOCALC ("tc") and
    Stixrude ("sb") endmember EOS, which need an iterative volume/EOS
    solve, the MELTS liquid volume model used here (a linear-in-(T,P)
    "Kress" polynomial, see GH_endmembers.h) integrates in closed form, so
    no Newton iteration is needed.

    G(T,P) is built the way MELTS itself builds it: from a *solid*
    reference state at Tr=298.15K, Pr=1bar (H, S, 4-term Berman Cp),
    integrated up to the species' fusion temperature Tfus, plus a fusion
    correction (deltaH_fus = Tfus*Sfus), plus the liquid's own constant Cp
    integrated from Tfus to T, plus a pressure correction from the Kress
    liquid-volume polynomial integrated from Pr to P.

    All internal arithmetic is done in MELTS' native units (J, bar, K,
    matching xMELTS' R=8.3143 J/K and Tr=298.15K/Pr=1bar reference state);
    the result is converted to MAGEMin's convention (kJ, kbar) only in the
    final assignment to PP_ref_db.gbase, mirroring how SB_G_EM_function
    converts its own bar/J-based EOS result via "/kbar2bar" at the very end.
*/

#include <math.h>
#include <string.h>
#include <float.h>
#include <stdio.h>

#include "../MAGEMin.h"
#include "GH_endmembers.h"
#include "GH_PP_endmembers.h"
#include "GH_fluid_eos.h"
#include "GH_gem_function.h"

#define GH_Tr    298.15   /* reference temperature (K), matches xMELTS TR      */
#define GH_Pr    1.0      /* reference pressure (bar), matches xMELTS PR      */
#define GH_kbar2bar 1000.0

/** GH_G_EM_function's own "EM_database" parameter is NOT the calibration-
    family id (0=xMELTS/1=rMELTS/2=pMELTS) despite the name - it's fed
    gv.EM_dataset (always 1 for gh, a leftover from tc/sb's own, unrelated
    "endmember dataset version" concept) all the way through the shared
    get_em_data->G_EM_function->GH_G_EM_function call chain, which has
    ~990 call sites across tc/sb/gh and was not worth widening just for
    this. Real gv.EM_database is threaded in here instead the same way
    GH_cpx_SS2/GH_spn_multistart_flag are: a plain global (not static,
    since GH_SS_objective_init_function that sets it lives in a different
    file, gh_objective_functions.c) set once at solution-phase init time.
    Found while wiring pMELTS liquid support, 2026-07-14 - see
    [[gh-multicalibration-xmelts-rmelts-pmelts]]. */
int GH_actual_EM_database = 0;

/** See GH_gem_function.h for the full explanation. Default 0 = pure
    standalone "water" phase semantics; set to 1 only around the liquid's
    own "H2O" basis-species lookups (gh_gss_function.c). */
int GH_H2O_liquid_context = 0;

/**
    Integrate the Berman (1988) solid Cp polynomial (k0,k1,k2,k3) from Tr to
    T, adding the order-disorder lambda-transition contribution (Tt, deltaH,
    l1, l2) when the phase has one (Tt != 0) - ported from xMELTS'
    sources/gibbs.c (the generic CP_BERMAN branch of intCpsolid, shared by
    every Berman-calibrated solid that isn't given its own dedicated
    pressure-shifted transition curve, e.g. nepheline/kalsilite there and
    quartz/tridymite/cristobalite here, see GH_PP_endmembers.h).
*/
static void GH_berman_HS(  double T, double H0, double S0, const double cp[8],
                            double *H_out, double *S_out               ){
    double k0 = cp[0], k1 = cp[1], k2 = cp[2], k3 = cp[3];
    double Tt = cp[4], deltaH = cp[5], l1 = cp[6], l2 = cp[7];

    double H = H0 + k0*(T-GH_Tr) + 2.0*k1*(sqrt(T)-sqrt(GH_Tr))
             - k2*(1.0/T - 1.0/GH_Tr) - 0.5*k3*(1.0/(T*T) - 1.0/(GH_Tr*GH_Tr));
    double S = S0 + k0*log(T/GH_Tr)
             - 2.0*k1*(1.0/sqrt(T) - 1.0/sqrt(GH_Tr))
             - 0.5*k2*(1.0/(T*T) - 1.0/(GH_Tr*GH_Tr))
             - (1.0/3.0)*k3*(1.0/(T*T*T) - 1.0/(GH_Tr*GH_Tr*GH_Tr));

    if (Tt != 0.0){
        if (T > Tt){
            H += deltaH + 0.5*l1*l1*(Tt*Tt-GH_Tr*GH_Tr)
               + (2.0/3.0)*l1*l2*(Tt*Tt*Tt-GH_Tr*GH_Tr*GH_Tr)
               + 0.25*l2*l2*(Tt*Tt*Tt*Tt-GH_Tr*GH_Tr*GH_Tr*GH_Tr);
            S += deltaH/Tt + l1*l1*(Tt-GH_Tr) + l1*l2*(Tt*Tt-GH_Tr*GH_Tr)
               + (1.0/3.0)*l2*l2*(Tt*Tt*Tt-GH_Tr*GH_Tr*GH_Tr);
        }
        else{
            H += 0.5*l1*l1*(T*T-GH_Tr*GH_Tr)
               + (2.0/3.0)*l1*l2*(T*T*T-GH_Tr*GH_Tr*GH_Tr)
               + 0.25*l2*l2*(T*T*T*T-GH_Tr*GH_Tr*GH_Tr*GH_Tr);
            S += l1*l1*(T-GH_Tr) + l1*l2*(T*T-GH_Tr*GH_Tr)
               + (1.0/3.0)*l2*l2*(T*T*T-GH_Tr*GH_Tr*GH_Tr);
        }
    }
    *H_out = H;
    *S_out = S;
}

/**
    Berman (1988) solid volume-EOS pressure integral (dG = V dP from Pr to
    P), ported from xMELTS' sources/gibbs.c (intEOSsolid, EOS_BERMAN case):
    V(T,P) = V0*(1 + v1*(P-Pr) + v2*(P-Pr)^2 + v3*(T-Tr) + v4*(T-Tr)^2).
    Closed form, no iterative solve needed (unlike sb's Birch-Murnaghan-
    Debye or tc's Modified Tait volume models).
*/
static double GH_berman_EOS_dG(double T, double P, double V0, const double eos[4]){
    double v1 = eos[0], v2 = eos[1], v3 = eos[2], v4 = eos[3];
    return V0*( (v1/2.0-v2)*(P*P-GH_Pr*GH_Pr) + v2*(P*P*P-GH_Pr*GH_Pr*GH_Pr)/3.0
              + (1.0-v1+v2+v3*(T-GH_Tr)+v4*(T-GH_Tr)*(T-GH_Tr))*(P-GH_Pr) );
}

/**
    Ideal-gas standard state for O2, ported from xMELTS' sources/gibbs.c
    (the unconditional "ideal gas treatment" fallback in gibbs(), used
    whenever none of that file's SPECIAL_O2/CONSTANT_P_O2/ZERO_O2 build
    macros are set - the vanilla case). H0=0, S0=205.15 J/K are O2's own
    elemental reference values (O2 gas *is* the reference state for
    oxygen), not the Tr,Pr=298.15K/1bar solid convention used elsewhere.
*/
static double GH_O2_G(double T, double P){
    double R = 8.3143;
    double H = 23.10248*(T-GH_Tr) + 2.0*804.8876*(sqrt(T)-sqrt(GH_Tr))
             - 1762835.0*(1.0/T - 1.0/GH_Tr) - 18172.9196*log(T/GH_Tr)
             + 0.5*0.002676*(T*T-GH_Tr*GH_Tr);
    double S = 205.15 + 23.10248*log(T/GH_Tr)
             - 2.0*804.8876*(1.0/sqrt(T) - 1.0/sqrt(GH_Tr))
             - 0.5*1762835.0*(1.0/(T*T) - 1.0/(GH_Tr*GH_Tr))
             + 18172.9196*(1.0/T - 1.0/GH_Tr) + 0.002676*(T-GH_Tr)
             - R*log(P/GH_Pr);
    return H - T*S;
}

/**
    Real MELTS treatment for the SiO2 polymorphs quartz, tridymite and
    cristobalite (ported from xMELTS' sources/gibbs.c special cases for
    each name, VALUE only - no lambdaV volume-derivative terms, since
    those never feed back into gs/hs in the source either, and MAGEMin
    gets volume/moduli from finite differences on G, not analytics).
    Below the pressure-shifted transition temperature cp_t = Tt + dTdP*
    (P-1), integrates the alpha-phase Cp + lambda transition (delt =
    Tt-cp_t shifted, reducing exactly to the plain lambda formula used
    elsewhere in gh when delt=0, e.g. tridymite at any P since its
    dTdP=0); above cp_t, switches to the beta phase's own H0/S0 (Berman
    1988) reintegrated with the SAME Cp coefficients from Tr, and its own
    V0/eos_berman for the pressure term.

    Two quartz-only quirks confirmed against source (2026-07-14, found by
    diffing gh's gbase directly against a real xMELTS Test_gibbs run,
    which showed quartz off by a near-constant ~1.29 kJ across a wide T/P
    sweep while chromite/hercynite/magnetite/spinel/ulvospinel matched to
    <0.001 kJ - ruling out the shared Berman/EOS engine and pointing at
    something quartz-specific):
    1. QUARTZ_ADJUSTMENT = -1291.0 J: a real, hardcoded empirical
       calibration correction in gibbs.c's RHYOLITE_ADJUSTMENTS branch
       (the one gh replicates - WETTING_ANGLE_CORR, the other additive
       term in the same source line, is always 0.0 so doesn't need
       porting), applied unconditionally to quartz's G regardless of
       alpha/beta branch. Not present for cristobalite/tridymite.
    2. The alpha-branch lambda-correction reference temperature: quartz's
       own source hardcodes the literal constant 373.0 (K) instead of Tr
       (298.15 K) - unlike tridymite/cristobalite, which both genuinely
       use Tr there. Kept as a quartz-only literal below, not routed
       through GH_Tr, to stay faithful to that source-level asymmetry.
*/
static double GH_SiO2_polymorph_G(int is_quartz, double T, double P,
                                   const PP_db_gh *alpha, const PP_db_gh_beta *beta){
    double k0 = alpha->cp_berman[0], k1 = alpha->cp_berman[1];
    double k2 = alpha->cp_berman[2], k3 = alpha->cp_berman[3];
    double Tt = alpha->cp_berman[4], l1 = alpha->cp_berman[6], l2 = alpha->cp_berman[7];
    double cp_t = Tt + beta->dTdP*(P-1.0);
    const double QUARTZ_ADJUSTMENT = -1291.0; /* J, RHYOLITE_ADJUSTMENTS branch of gibbs.c */

    double H, S, V0; const double *eos;
    if (T > cp_t){
        H = beta->H + k0*(T-GH_Tr) + 2.0*k1*(sqrt(T)-sqrt(GH_Tr))
          - k2*(1.0/T-1.0/GH_Tr) - 0.5*k3*(1.0/(T*T)-1.0/(GH_Tr*GH_Tr));
        S = beta->S + k0*log(T/GH_Tr)
          - 2.0*k1*(1.0/sqrt(T)-1.0/sqrt(GH_Tr))
          - 0.5*k2*(1.0/(T*T)-1.0/(GH_Tr*GH_Tr))
          - (1.0/3.0)*k3*(1.0/(T*T*T)-1.0/(GH_Tr*GH_Tr*GH_Tr));
        V0 = beta->V; eos = beta->eos_berman;
    } else {
        H = alpha->H + k0*(T-GH_Tr) + 2.0*k1*(sqrt(T)-sqrt(GH_Tr))
          - k2*(1.0/T-1.0/GH_Tr) - 0.5*k3*(1.0/(T*T)-1.0/(GH_Tr*GH_Tr));
        S = alpha->S + k0*log(T/GH_Tr)
          - 2.0*k1*(1.0/sqrt(T)-1.0/sqrt(GH_Tr))
          - 0.5*k2*(1.0/(T*T)-1.0/(GH_Tr*GH_Tr))
          - (1.0/3.0)*k3*(1.0/(T*T*T)-1.0/(GH_Tr*GH_Tr*GH_Tr));

        double delt = Tt - cp_t;
        double x1 = l1*l1*delt + 2.0*l1*l2*delt*delt + l2*l2*delt*delt*delt;
        double x2 = l1*l1 + 4.0*l1*l2*delt + 3.0*l2*l2*delt*delt;
        double x3 = 2.0*l1*l2 + 3.0*l2*l2*delt;
        double x4 = l2*l2;
        double Tr_shift = (is_quartz ? 373.0 : GH_Tr) - delt;
        H += x1*(T-Tr_shift) + x2*(T*T-Tr_shift*Tr_shift)/2.0
           + x3*(T*T*T-Tr_shift*Tr_shift*Tr_shift)/3.0
           + x4*(T*T*T*T-Tr_shift*Tr_shift*Tr_shift*Tr_shift)/4.0;
        S += x1*(log(T)-log(Tr_shift)) + x2*(T-Tr_shift)
           + x3*(T*T-Tr_shift*Tr_shift)/2.0 + x4*(T*T*T-Tr_shift*Tr_shift*Tr_shift)/3.0;

        V0 = alpha->V; eos = alpha->eos_berman;
    }
    double G = (H - T*S) + GH_berman_EOS_dG(T, P, V0, eos);
    if (is_quartz){ G += QUARTZ_ADJUSTMENT; }
    return G;
}

/**
    Sanidine (K-feldspar) Al-Si disordering correction, ported from
    xMELTS' sources/gibbs.c "sanidine" special case (value only - the
    extra Cp/pressure-cross terms in the source only refine hs/ss/cps for
    T<1436K and never feed back into gs beyond what's captured here via
    the td=min(1436,T) cap). Zero below Tr (T<=Tr).
*/
static double GH_sanidine_ordering_G(double T, double P){
    double td = (T < 1436.0) ? T : 1436.0;
    double d0=282.98, d1=-4.83e3, d2=36.21e5, d3=-15.733e-2, d4=34.770e-6, d5=41.063e4;

    if (T <= GH_Tr) return 0.0;

    double dhdis = d0*(td-GH_Tr) + 2.0*d1*(sqrt(td)-sqrt(GH_Tr)) - d2*(1.0/td-1.0/GH_Tr)
                 + d3*(td*td-GH_Tr*GH_Tr)/2.0 + d4*(td*td*td-GH_Tr*GH_Tr*GH_Tr)/3.0;
    double dsdis = d0*(log(td)-log(GH_Tr)) - 2.0*d1*(1.0/sqrt(td)-1.0/sqrt(GH_Tr))
                 - d2*(1.0/(td*td)-1.0/(GH_Tr*GH_Tr))/2.0 + d3*(td-GH_Tr)
                 + d4*(td*td-GH_Tr*GH_Tr)/2.0;
    double dvdis = dhdis/d5;

    return dhdis - T*dsdis + dvdis*(P-1.0);
}

/**
    Albite Al-Si tricritical ordering correction (Salje 1985), ported
    from xMELTS' sources/albite.c: a 2-order-parameter (q, qod) Landau
    model, Newton-solved here directly (2x2 analytic Hessian, same
    iteration albite.c's own order() performs) since this is a single
    pure-phase standard-state evaluation, not a joint solution-phase
    optimization - no NLopt/order-parameter machinery needed. Returns 0
    when the disordered solution (q=qod=0) is the equilibrium state
    (matches albite.c's "skip" test), i.e. above ~tc=1251K.
*/
static double GH_albite_ordering_G(double T, double P){
    const double a0=5.479, b=6854.0, aod0=41.620, bod=-9301.0, cod=43600.0;
    const double d0=-2.171, d1=-3.043, d2=-0.001569, d3=0.000002109;
    const double tc=1251.0, tod=824.1;
    const double vscale = 335282.925;

    double s0 = 0.6, s1 = 0.9;
    for (int iter = 0; iter < 200; iter++){
        double pf = 1.0 + (P-1.0)/vscale;
        double dgds0 = (-a0*tc*s0 + b*s0*s0*s0 + (d0-d2*T*T-2.0*d3*T*T*T)*s1)*pf
                     + T*(a0*s0 + (d1+2.0*d2*T+3.0*d3*T*T)*s1);
        double dgds1 = (-aod0*tod*s1 + bod*s1*s1*s1 + cod*s1*s1*s1*s1*s1
                        + (d0-d2*T*T-2.0*d3*T*T*T)*s0)*pf
                     + T*(aod0*s1 + (d1+2.0*d2*T+3.0*d3*T*T)*s0);

        double j00 = (-a0*tc + 3.0*b*s0*s0)*pf + T*a0;
        double j01 = (d0-d2*T*T-2.0*d3*T*T*T)*pf + T*(d1+2.0*d2*T+3.0*d3*T*T);
        double j11 = (-aod0*tod + 3.0*bod*s1*s1 + 5.0*cod*s1*s1*s1*s1)*pf + T*aod0;

        double det = j00*j11 - j01*j01;
        if (fabs(det) < 1e-30) break;
        double inv00 =  j11/det, inv01 = -j01/det, inv11 =  j00/det;

        double ds0 = inv00*dgds0 + inv01*dgds1;
        double ds1 = inv01*dgds0 + inv11*dgds1;
        s0 -= ds0; s1 -= ds1;

        if (fabs(ds0) < 1e-10 && fabs(ds1) < 1e-10) break;
    }

    if (fabs(s0) < 1.4901161193847656e-8 && fabs(s1) < 1.4901161193847656e-8){
        return 0.0;   /* disordered (high albite/monalbite): no correction */
    }

    double pf = 1.0 + (P-1.0)/vscale;
    double S = -(0.5*a0*s0*s0 + 0.5*aod0*s1*s1 + (d1+2.0*d2*T+3.0*d3*T*T)*s0*s1);
    double H = -0.5*a0*tc*s0*s0 + 0.25*b*s0*s0*s0*s0
             - 0.5*aod0*tod*s1*s1 + 0.25*bod*s1*s1*s1*s1
             + (cod/6.0)*s1*s1*s1*s1*s1*s1
             + (d0-d2*T*T-2.0*d3*T*T*T)*s0*s1;
    double V = H/vscale;
    return H - T*S + (P-1.0)*V;
}

/**
    Rhombohedral-oxide (geikielite-hematite-ilmenite-pyrophanite-corundum,
    Ghiorso & Evans 2008) intracrystalline-ordering correction for the
    ilmenite/geikielite/pyrophanite standard states, ported from xMELTS'
    sources/rhomsghiorso.c (pureOrder + IL_G/GK_G/PY_G macros): a single
    Landau-type order parameter s in [0,1) on that endmember's own A/B
    cation sites, bounded-Newton-solved (same iteration pureOrder itself
    performs, start s=0.98). dh/wh/dv/wv are the (Hordering, Wordering,
    dVordering, WVordering) constants for ilmenite/geikielite/pyrophanite -
    numerically identical across all three in xMELTS' own calibration, but
    kept as distinct call-site arguments to stay faithful to the 3
    independently-named macro sets (IL_, GK_, PY_) in the source.
*/
static double GH_rhm_pure_order_G(double T, double P, double dh, double wh, double dv, double wv){
    double R = 8.3143;
    double s = 0.98;
    for (int iter = 0; iter < 1000; iter++){
        double dgds   = R*T*(log(1.0+s) - log(1.0-s))
                      - 2.0*(dh+(P-1.0)*dv)*s + (wh+(P-1.0)*wv)*(2.0*s - 4.0*s*s*s);
        double d2gds2 = R*T*(1.0/(1.0+s) + 1.0/(1.0-s))
                      - 2.0*(dh+(P-1.0)*dv) + (wh+(P-1.0)*wv)*(2.0 - 12.0*s*s);
        double sNew   = s - dgds/d2gds2;
        sNew = (sNew > 1.0 - 10.0*DBL_EPSILON) ? 1.0 - 10.0*DBL_EPSILON : sNew;
        sNew = (sNew < 0.0) ? 0.0 : sNew;
        int converged = (fabs(sNew - s) <= 10.0*DBL_EPSILON);
        s = sNew;
        if (converged) break;
    }
    double S = -R*((1.0+s)*log(1.0+s) + (1.0-s)*log(1.0-s) - 2.0*log(2.0));
    double H = (dh+(P-1.0)*dv)*(1.0 - s*s) + (wh+(P-1.0)*wv)*s*s*(1.0 - s*s);
    return H - T*S;
}

/**
    Rhombohedral-oxide short-range-order correction for the hematite/
    corundum standard states used inside that solution (xMELTS' HM_G/CR_G
    macros): fSRO(T,order) is a natural cubic spline through 8 calibrated
    points (SRO600..SRO1300), but every one of those calibration values is
    the identical constant 0.0730205 - the spline is exactly flat, so
    fSRO(T,0)==0.0730205 and all its T-derivatives are 0 for any T,
    collapsing the general spline machinery to this closed form (P-
    independent, matching DHM_GDP=DCR_GDP=0 in the source).
*/
static double GH_rhm_SRO_G(double T){
    const double SROconst = 0.0730205;
    double R = 8.3143;
    return -T*(SROconst*2.0*R*log(2.0));
}

/**
    Self-consistent thermal Vinet finite-strain volume EOS, ported from
    xMELTS' sources/gibbs.c (the EOS_VINET branch of gibbs()) - a genuinely
    different volume-integration engine from the Berman polynomial EOS used
    everywhere else in gh, needed only for the na-nepheline term inside the
    vc-nepheline/ca-nepheline "phantom reaction" standard states (see
    GH_vcneph_G/GH_caneph_G below). Only the Gibbs energy contribution is
    ported (matching gh's general policy of getting volume/moduli from
    finite differences on G, not analytic derivatives) - the h/s/cp/v/dvdt/
    dvdp/etc terms the source also computes are not needed here.
*/
static double GH_vinet_EOS_dG(double T, double P, double V0, double alpha, double K, double Kp){
    double Tr = GH_Tr, Pr = GH_Pr;
    double eta = 3.0*(Kp-1.0)/2.0;
    double x = 1.0, x0 = 1.0;
    double fn, dfn;
    int iter;

    iter = 0;
    do {
        fn  = x*x*(P/10000.0) - 3.0*K*(1.0-x)*exp(eta*(1.0-x)) - x*x*alpha*K*(T-Tr);
        dfn = 2.0*x*(P/10000.0) + 3.0*K*(1.0+eta*(1.0-x))*exp(eta*(1.0-x)) - 2.0*alpha*K*(T-Tr);
        x = x - fn/dfn;
        iter++;
    } while ((iter < 500) && (fn*fn > DBL_EPSILON));

    iter = 0;
    do {
        fn  = x0*x0*(Pr/10000.0) - 3.0*K*(1.0-x0)*exp(eta*(1.0-x0)) - x0*x0*alpha*K*(T-Tr);
        dfn = 2.0*x0*(Pr/10000.0) + 3.0*K*(1.0+eta*(1.0-x0))*exp(eta*(1.0-x0)) - 2.0*alpha*K*(T-Tr);
        x0 = x0 - fn/dfn;
        iter++;
    } while ((iter < 500) && (fn*fn > DBL_EPSILON));

    double a  =  (9.0*V0*K/(eta*eta))*(1.0 - eta*(1.0-x))*exp(eta*(1.0-x));
           a +=  V0*(T-Tr)*K*alpha*(x*x*x - 1.0) - 9.0*V0*K/(eta*eta);
           a -=  (9.0*V0*K/(eta*eta))*(1.0 - eta*(1.0-x0))*exp(eta*(1.0-x0));
           a -=  V0*(T-Tr)*K*alpha*(x0*x0*x0 - 1.0) - 9.0*V0*K/(eta*eta);

    return -a*10000.0 + P*V0*x*x*x - Pr*V0*x0*x0*x0;
}

/**
    Na-nepheline evaluated with a phase-specific override volume V0 and the
    Vinet EOS above instead of its own standalone Berman EOS (ported from
    xMELTS' gibbs.c "vc-nepheline"/"ca-nepheline" branches' inner
    gibbs(t,p,"na-nepheline",&tempRef,...) call - same H0/S0/Cp as
    na-nepheline's own standalone entry ("nane" in GH_PP_endmembers.c),
    only V0 and the EOS differ).
*/
static double GH_naneph_vinet_G(double T, double P, double V0){
    double H_T, S_T;
    const double cp[8] = { 205.24*4.0, -7.599E2*4.0, -108.383E5*4.0, 208.182E7*4.0,
                            467.0, 241.0*4.0, -50.249E-2*2.0, 165.95E-5*2.0 };
    GH_berman_HS(T, -2093004.0*4.0, 124.641*4.0, cp, &H_T, &S_T);
    return (H_T - T*S_T) + GH_vinet_EOS_dG(T, P, V0, 31.802e-6, 48.7805, 1.4747);
}

/**
    Vc-nepheline (Na3Al3Si5O16) real standard state, ported from xMELTS'
    gibbs.c: a "phantom" 2:1 reaction [2*high-albite + na-nepheline]/2,
    where "high-albite" is the PLAIN Berman H/S/Cp evaluation of "vcne"'s
    own table row (GH_PP_endmembers.c - real xMELTS bakes high-albite's
    calibration directly into that row, with V=0/EOS=0 so the ordinary
    generic Berman dispatch already gives exactly the "high-albite, no
    pressure term" piece) plus 0.5x na-nepheline's own Vinet-EOS-corrected
    term computed here. Added on top of the generic gbase_J in
    GH_G_EM_function, mirroring albite/sanidine's ordering-correction
    dispatch pattern (not a replacement, an addition).
*/
static double GH_vcneph_G(double T, double P){
    return 0.5*GH_naneph_vinet_G(T, P, 5.434181*8.0);
}

/**
    Ca-nepheline (CaNa2Al4Si4O16) real standard state, ported from xMELTS'
    gibbs.c: same "phantom" 2:1 reaction construction as vc-nepheline
    (with "high-anorthite" = "cane"'s own plain-Berman table row, V=0/
    EOS=0), plus a real +23096 J enthalpy correction and a real +15.8765
    J/K zero-point entropy correction (both applied to the high-anorthite
    side before averaging, per the source) plus 0.5x na-nepheline's own
    Vinet-EOS-corrected term.
*/
static double GH_caneph_G(double T, double P){
    return 23096.0 - T*15.8765 + 0.5*GH_naneph_vinet_G(T, P, 5.433181*8.0);
}

/**
    pMELTS liquid pressure correction, ported from xMELTS' sources/gibbs.c
    (the "if (calculationMode == MODE_pMELTS)" block, ~line 1493-1595):
    a genuinely different volume-integration engine from the Kress
    polynomial used for xMELTS/rMELTS - NOT an addition on top of the
    Kress dG_P term, a full REPLACEMENT of it (confirmed by tracing
    gibbs.c's own control flow: the Kress branch a few lines above is
    gated on MODE__MELTS/MODE__MELTSandCO2/MODE__MELTSandCO2_H2O or
    MODE_xMELTS - MODE_pMELTS matches none of those, so hl/sl reaching
    this block are still the plain Pr-reference (Tr-to-Tfus-to-T
    integrated) values, exactly matching gh's own G_Pr before this call).

    A 3rd-order Birch-Murnaghan EOS with Kp hardcoded to 5.0 for pMELTS
    (real source: "double Kp = (calculationMode == MODE_pMELTS) ? 5.0 :
    liquid->eos.Kress.d2vdp2;"), volume solved by Newton's method from the
    reference volume v0 (the Kress dV/dT-corrected but NOT dV/dP-corrected
    volume at Pr), then G = G_Pr + p*v - pr*v0 + minusIntPdV (a BM3
    pressure-integral closed form). Value-only port (drops the extensive
    T-derivative bookkeeping real gibbs.c also carries for Cp/dCp/dT,
    matching gh's general policy elsewhere).

    Trigger condition (real source: "p > pr && (dvdp+d2vdtp*(t-trl)) !=
    0.0") - when false (essentially only if a pMELTS endmember's Kress
    dV/dP terms are literally zero, or in the P<=Pr edge case never hit by
    gh's own P>=1kbar grid), real gibbs.c applies NO pressure correction
    at all (not even the trivial Vl*(P-Pr) Kress linear term) - a real,
    distinct-from-Kress fallback behavior, ported here as returning 0.0.

    Found missing entirely 2026-07-15/16 during the 10x10 grid sweep -
    caught because the earlier Phase 3 pMELTS liquid port verification was
    only ever checked at one low-pressure point (~2kbar) where this
    correction is still negligible. Grows to multi-kJ magnitude by
    25kbar. See [[gh-multicalibration-xmelts-rmelts-pmelts]] and
    [[gh-spn-liq-gbase-verification]].
*/
static double GH_pmelts_liq_BM3_dG(double T, double P, double Vl,
                                    double dvdt, double dvdp, double d2vdtp){
    const double Trl = 1673.0;
    double v0     = Vl + dvdt*(T - Trl);
    double dv0dp  = dvdp + d2vdtp*(T - Trl);

    if (!(P > GH_Pr) || dv0dp == 0.0){
        return 0.0;
    }

    double K  = -v0/dv0dp;
    const double Kp = 5.0;

    double v = v0*0.99, vLast = v0;
    int iter = 0;
    while (fabs(v - vLast) > 10.0*DBL_EPSILON && iter < 1000){
        double r23 = pow(v0/v, 2.0/3.0), r43 = pow(v0/v, 4.0/3.0);
        double r53 = pow(v0/v, 5.0/3.0), r73 = pow(v0/v, 7.0/3.0);
        double rm13 = pow(v0/v, -1.0/3.0);
        double bracket = 1.0 - 0.75*(4.0-Kp)*(r23 - 1.0);

        double fn  = 1.5*K*(r73 - r53)*bracket - P;
        double dfn = 1.5*K*((7.0/3.0)*r43 - (5.0/3.0)*r23)*(-v0/(v*v))*bracket
                   + 1.5*K*(r73 - r53)*(-0.75*(4.0-Kp)*(2.0/3.0)*rm13)*(-v0/(v*v));

        vLast = v;
        v    += -fn/dfn;
        if (v > v0)      v = v0;
        if (v < 0.01*v0) v = 0.01*v0;
        iter++;
    }

    double f = (pow(v0/v, 2.0/3.0) - 1.0)/2.0;
    double minusIntPdV = 4.5*K*(Kp-4.0)*v0*f*f*f + 4.5*K*v0*f*f;

    return P*v - GH_Pr*v0 + minusIntPdV;
}

static PP_ref GH_pack_PP_ref(char *name, int len_ox, const double *Comp,
                              double *bulk_rock, double *apo, double gbase_J){
    PP_ref PP_ref_db;
    strcpy(PP_ref_db.Name, name);
    for (int i = 0; i < len_ox; i++){
        PP_ref_db.Comp[i] = Comp[i];
    }
    double fbc = 0.0;
    for (int i = 0; i < len_ox; i++){
        fbc += bulk_rock[i]*apo[i];
    }
    double ape = 0.0;
    for (int i = 0; i < len_ox; i++){
        ape += Comp[i]*apo[i];
    }
    PP_ref_db.gbase             = gbase_J/GH_kbar2bar;   /* J -> kJ */
    PP_ref_db.factor             = fbc/ape;
    PP_ref_db.phase_shearModulus = 0.0;
    PP_ref_db.phase_bulkModulus  = 0.0;
    PP_ref_db.phase_expansivity  = 0.0;
    PP_ref_db.phase_cp           = 0.0;
    return PP_ref_db;
}

PP_ref GH_G_EM_function(   int          EM_database,
                            int          len_ox,
                            int         *id,
                            double      *bulk_rock,
                            double      *apo,
                            double       Pkbar,
                            double       T,
                            char        *name,
                            char        *state          ){

    double P  = Pkbar * GH_kbar2bar;   /* kbar -> bar */

    /* Common rock-forming pure phases (research group "gh"'s own thermodynamic
       basis, ported from xMELTS' solid-phase database - see
       GH_PP_endmembers.h), not part of the 13-species liquid endmember table
       looked up via find_EM_id() below: must be special-cased before that
       lookup, since they aren't in gh's endmember hashtable.                  */
    /* SiO2 polymorphs: real MELTS dual alpha/beta reference-state
       treatment (pressure-shifted transition), not the generic Berman-
       solid path below - see GH_SiO2_polymorph_G(). */
    int beta_id = GH_find_SiO2_beta_id(name);
    if (beta_id >= 0 && (strcmp(name,"q")==0 || strcmp(name,"crst")==0 || strcmp(name,"trd")==0)){
        int pp_id_sio2 = GH_find_PP_id(name);
        PP_db_gh alpha = Access_GH_PP_DB(pp_id_sio2);
        PP_db_gh_beta beta = Access_GH_SiO2_beta_DB(beta_id);
        double gbase_J = GH_SiO2_polymorph_G(strcmp(name,"q")==0, T, P, &alpha, &beta);
        return GH_pack_PP_ref(name, len_ox, alpha.Comp, bulk_rock, apo, gbase_J);
    }

    int pp_id = GH_find_PP_id(name);
    if (pp_id >= 0){
        PP_db_gh PP_return = Access_GH_PP_DB(pp_id);
        double H_T, S_T;
        GH_berman_HS(T, PP_return.H, PP_return.S, PP_return.cp_berman, &H_T, &S_T);
        double gbase_J = (H_T - T*S_T) + GH_berman_EOS_dG(T, P, PP_return.V, PP_return.eos_berman);
        /* real MELTS Al-Si order-disorder corrections for these two
           feldspar endmember standard states (sources/gibbs.c/albite.c) */
        if (strcmp(name, "ab") == 0){
            gbase_J += GH_albite_ordering_G(T, P);
        }
        else if (strcmp(name, "san") == 0){
            gbase_J += GH_sanidine_ordering_G(T, P);
        }
        /* rhm-oxide (geikielite-hematite-ilmenite-pyrophanite-corundum,
           Ghiorso & Evans 2008) intracrystalline-ordering / short-range-
           order corrections - see GH_rhm_pure_order_G/GH_rhm_SRO_G above.
           dh/wh/dv/wv are numerically identical across ilm/gei/pyr in
           xMELTS' own calibration (dhilm=dhgei=dhpyr=17477.0 J,
           whilm=whgei=whpyr=3189.0 J, dvilm=dvgei=dvpyr=0.010758 J/bar,
           wvilm=wvgei=wvpyr=0.035089 J/bar). */
        else if (strcmp(name, "ilm") == 0 || strcmp(name, "gei") == 0 || strcmp(name, "pyr") == 0){
            gbase_J += GH_rhm_pure_order_G(T, P, 17477.0, 3189.0, 0.010758, 0.035089);
        }
        else if (strcmp(name, "hem") == 0 || strcmp(name, "crn") == 0){
            gbase_J += GH_rhm_SRO_G(T);
        }
        /* nepheline/kalsilite "phantom reaction" standard states - see
           GH_vcneph_G/GH_caneph_G above */
        else if (strcmp(name, "vcne") == 0){
            gbase_J += GH_vcneph_G(T, P);
        }
        else if (strcmp(name, "cane") == 0){
            gbase_J += GH_caneph_G(T, P);
        }
        return GH_pack_PP_ref(name, len_ox, PP_return.Comp, bulk_rock, apo, gbase_J);
    }
    if (strcmp(name, "O2") == 0){
        double O2_Comp[16] = {0,0,0,0,0,0,0,0,2.0,0,0,0,0,0,0,0};   /* O2 = 2 "O" (excess-oxygen axis) */
        double gbase_J = GH_O2_G(T, P);
        return GH_pack_PP_ref(name, len_ox, O2_Comp, bulk_rock, apo, gbase_J);
    }

    /* endmember lookup (implicit declaration of find_EM_id, resolved at
       link time against src/hash_init.h's definition - same pattern used
       by SB_G_EM_function/TC_G_EM_function, neither of which includes
       hash_init.h directly either)                                        */
    int p_id            = find_EM_id(name);
    EM_db_gh EM_return   = Access_GH_EM_DB(GH_actual_EM_database, p_id);

    /* H2O and CO2 liquid standard states. NOT GH_pitzer_sterner_G/
       fluidPhase() - that model is a separate real-gas EOS real gibbs.c
       only uses for pMELTS' H2O (fixed-9550bar branch, not yet ported)
       and for gh's own "fl" mixed-fluid solution phase (GH_fluid_eos.c).
       For xMELTS/rMELTS, real gibbs.c's actual liquid H2O/CO2 branches
       (MODE_xMELTS/MODE__MELTSandCO2_H2O, ~line 1018/1063) build G from a
       Robie-type ideal-gas polynomial + a real EOS at a FIXED 1-bar
       reference (Haar 1984 for H2O via whaar(), Duan 1992 for CO2 via
       propertiesOfPureCO2()) plus a simple linear-in-(T,P) volume
       correction (Oaks-Lange for H2O, a Kress-like linear term for CO2) -
       found 2026-07-15 when live-verifying gh's liq gbase against real
       MELTS output showed GH_pitzer_sterner_G was never actually right
       here (off by ~2.7-10.3 kJ for CO2, ~47-50 kJ for H2O, both real,
       non-constant-offset discrepancies). GH_haar_H2O_G/GH_duan_CO2_G
       (GH_fluid_eos.c) port the real whaar()/propertiesOfPureCO2()
       formulas faithfully (verified bit-exact against a standalone
       harness linked against real xMELTS). phase->h/phase->s in the real
       formula (the liquid table's own H2O/CO2 "ref" H0/S0) are exactly
       0 in real xMELTS' own table too (confirmed via the same harness),
       matching gh's existing all-zero EM_return.H0/S0 for these two -
       no separate correction term needed for that part.
       See [[gh-multicalibration-xmelts-rmelts-pmelts]]. */
    /* Standalone "water" pure phase (gv.PP_list's "H2O" entry, plus
       toolkit.c's system_aH2O activity calc, both reached via this same
       name-dispatch since gh's PP list keeps the "H2O" string rather than
       renaming to real xMELTS' own distinct "water" - see
       GH_H2O_liquid_context's header comment). Real gibbs.c's own
       "water" branch (line ~2326) is whaar()-at-actual-P(capped 10000)+
       wdh78()-above-that-cap for xMELTS/rMELTS - no Robie/Oaks-Lange
       wrapper (that wrapper is specific to the liquid's own "H2O" branch
       below). Found missing 2026-07-15 during the 10x10 grid sweep (gh
       previously used GH_pitzer_sterner_G here unconditionally, off by
       several kJ) - see [[gh-spn-liq-gbase-verification]].

       pMELTS' own "water" branch: real gibbs.c calls fluidPhase(t,p,x,...)
       at the ACTUAL pressure (unlike the liquid's own H2O, which uses a
       fixed 9550 bar reference - see the liquid-context branch below) and
       explicitly OVERRIDES its own *g output with a *h-t*s* recomputation
       ("used gH2O for calibration", gibbs.c's own comment) - these two
       disagree by ~70 kJ (confirmed via a standalone fluidPhase() harness
       at T=1273.15K/P=5000bar: g=-352.33 kJ vs h-Ts=-422.99 kJ), because
       the Berman reference-state shifts fluid.c applies to g vs (h,s) are
       independently calibrated constants, not a thermodynamically
       consistent triple - shifting g alone does not preserve g=h-Ts once
       shifted. Fixed 2026-07-16 via GH_pitzer_sterner_H2O_hTs_G
       (GH_fluid_eos.c), which recovers h-T*s from GH_pitzer_sterner_G's
       raw g plus fluid.c's own known shift constants, without needing to
       port fluidPhase()'s full dA/dT machinery - see its header comment
       for the derivation, verified exact against a real fluidPhase()
       harness. */
    if (strcmp(name, "H2O") == 0 && !GH_H2O_liquid_context){
        double gbase_J;
        if (GH_actual_EM_database == 2){
            /* No 10000 bar cap here (unlike the whaar()-based xMELTS/
               rMELTS branch below) - real gibbs.c's pMELTS "water" branch
               calls fluidPhase() at the actual p directly, uncapped
               (confirmed by reading the source: only the whaar() branch
               has "pHaar = (p<=10000.0)?p:10000.0"). */
            gbase_J = GH_pitzer_sterner_H2O_hTs_G(T, P);
        }
        else {
            double Pcap = (P <= 10000.0) ? P : 10000.0;
            gbase_J = GH_haar_H2O_G(T, Pcap);
            if (P > 10000.0){
                gbase_J += GH_wdh78_G(T, P);
            }
        }
        return GH_pack_PP_ref(name, len_ox, EM_return.Comp, bulk_rock, apo, gbase_J);
    }

    if (strcmp(name, "H2O") == 0 || strcmp(name, "CO2") == 0){
        double gbase_J;
        if (GH_actual_EM_database == 2 && strcmp(name, "CO2") == 0){
            /* pMELTS' own liquid table has no real-gas EOS treatment for
               CO2 - real gibbs.c leaves it as a literal placeholder zero
               under MODE_pMELTS (verified directly via a standalone
               gibbs() harness returning exactly G=0.0 for pMELTS "CO2") -
               a real, confirmed-correct source behavior, not a porting
               gap. */
            gbase_J = 0.0;
        }
        else if (GH_actual_EM_database == 2){
            /* pMELTS' own liquid "H2O" branch (real gibbs.c line ~939-956,
               MODE_pMELTS): calls fluidPhase(t, 9550.0, x=[1,0], ...) at a
               FIXED 9550 bar reference (NOT the actual P), computing
               gl=gH2O+r*t*(a/t+b) with a=phase->h/r=0, b=-phase->s/r=0 for
               gh's all-zero H2O ref data - so THIS branch's own gl is just
               gH2O (fluidPhase's raw g). BUT gh's H2O liquid-table row also
               has all-zero Kress dV/dP terms (matching real xMELTS' own
               table), which means the *downstream*, unconditional
               Birch-Murnaghan block (gibbs.c ~line 1493) always takes its
               "dv0dp==0" else-branch for H2O, and THAT branch
               unconditionally overwrites gl = hl - t*sl using the hl/sl
               this H2O branch left behind (hl=hH2O, sl=sH2O since a=b=0) -
               so the FINAL liquid gbase is actually hH2O-T*sH2O, not the
               raw gH2O. Confirmed 2026-07-16 by direct comparison against
               a real gibbs() dispatch call (not assumed from reading the
               H2O branch in isolation, which was my first, wrong attempt
               at this fix - g and h-Ts differ by ~70kJ here, same root
               cause as the standalone "water" phase below).
               GH_pitzer_sterner_H2O_hTs_G (GH_fluid_eos.c) computes this
               h-T*s value without needing to port fluidPhase()'s dA/dT
               machinery - see its header comment for the derivation.
               Fixed 2026-07-16, closing the last pMELTS gap from the
               10x10 grid sweep - see [[gh-spn-liq-gbase-verification]]. */
            gbase_J = GH_pitzer_sterner_H2O_hTs_G(T, 9550.0);
        }
        else if (strcmp(name, "H2O") == 0){
            const double r_ = 8.3143;
            double a = -33676.0 + 2783.6851512128/r_;   /* phase->h = 0 */
            double b = 18.3527  - 2.3838467967178/r_;   /* phase->s = 0 */
            double gRobie = r_*T*(2.9147*log(T) - 9.6863e-4*T + 6.8593e-8*T*T
                                + 77.8899/sqrt(T) - 28954.8/T - 2263.27/(T*T) - 15.8997);
            const double vOaksLange = 2.775, dvdtOaksLange = 1.086e-3, dvdpOaksLange = -0.382e-4;
            gbase_J = r_*T*(a/T + b) + GH_haar_H2O_G(T, 1.0) - gRobie
                    + vOaksLange*(P-1.0) + dvdtOaksLange*(T-1673.15)*(P-1.0)
                    + 0.5*dvdpOaksLange*(P-1.0)*(P-1.0);
        }
        else {
            const double hCO2 = -630.93193811701, sCO2 = -109.39331414050;
            const double vCO2 = 4.0157994267547, dvCO2dt = 1.213189e-3, dvCO2dp = -0.4267387e-4;
            const double pr_ref = 1.0, trl = 1673.0;
            gbase_J = GH_duan_CO2_G(T) + hCO2 - T*sCO2
                    + vCO2*(P-pr_ref) + dvCO2dt*(T-trl)*(P-pr_ref)
                    + dvCO2dp*(P*P/2.0 - pr_ref*pr_ref/2.0 - pr_ref*(P-pr_ref));
        }
        return GH_pack_PP_ref(name, len_ox, EM_return.Comp, bulk_rock, apo, gbase_J);
    }

    /* Liquid SiO2: real gibbs.c gives it its own dedicated, directly-
       calibrated empirical H0/S0/Cp(T) polynomial (h0_sio2=-901554.0,
       s0_sio2=48.475, Cp=al+bl*T+cl/T^2+dl/sqrt(T)) integrated straight
       from Tr, with a glass-transition-like threshold tg_sio2=1480K above
       which a simpler constant-Cp=81.373 form takes over - NOT the
       generic solid(beta-cristobalite)+fusion-correction+constant-Cpl
       construction used below for every other oxide liquid component.
       Found missing 2026-07-14 while checking gh's liq gbase against
       real MELTS at a genuinely spinel-stable T/P/bulk: after fixing the
       separate GH_liq_Trl bug (see below), 10 of 11 oxide endmembers
       matched to machine precision but SiO2 still showed a real ~0.2 kJ
       gap - traced to this entirely separate model, confirmed by
       reproducing gibbs.c's own SiO2 branch by hand and matching real
       MELTS' computed G to 6 decimal places. Uses the SAME Kress EOS
       pressure/volume terms (same Vl/kress[] data) as the generic path -
       only the T-only H0/S0/Cp part differs. See
       [[gh-spn-liq-gbase-verification]]. */
    if (strcmp(name, "SiO2") == 0){
        const double h0_sio2 = -901554.0, s0_sio2 = 48.475;
        const double al_sio2 = 127.200, bl_sio2 = -10.777e-3, cl_sio2 = 4.3127e5, dl_sio2 = -1463.8;
        const double tg_sio2 = 1480.0, cp_sio2 = 81.373;

        double hl, sl;
        if (T >= tg_sio2){
            hl = h0_sio2 + al_sio2*(tg_sio2-GH_Tr) + bl_sio2*(tg_sio2*tg_sio2-GH_Tr*GH_Tr)/2.0
               - cl_sio2*(1.0/tg_sio2-1.0/GH_Tr) + 2.0*dl_sio2*(sqrt(tg_sio2)-sqrt(GH_Tr));
            sl = s0_sio2 + al_sio2*log(tg_sio2/GH_Tr) + bl_sio2*(tg_sio2-GH_Tr)
               - (cl_sio2/2.0)*(1.0/(tg_sio2*tg_sio2)-1.0/(GH_Tr*GH_Tr)) - 2.0*dl_sio2*(1.0/sqrt(tg_sio2)-1.0/sqrt(GH_Tr));
            hl += (T-tg_sio2)*cp_sio2;
            sl += cp_sio2*log(T/tg_sio2);
        } else {
            hl = h0_sio2 + al_sio2*(T-GH_Tr) + bl_sio2*(T*T-GH_Tr*GH_Tr)/2.0
               - cl_sio2*(1.0/T-1.0/GH_Tr) + 2.0*dl_sio2*(sqrt(T)-sqrt(GH_Tr));
            sl = s0_sio2 + al_sio2*log(T/GH_Tr) + bl_sio2*(T-GH_Tr)
               - (cl_sio2/2.0)*(1.0/(T*T)-1.0/(GH_Tr*GH_Tr)) - 2.0*dl_sio2*(1.0/sqrt(T)-1.0/sqrt(GH_Tr));
        }

        const double GH_liq_Trl_sio2 = 1673.0;
        double dT = T - GH_liq_Trl_sio2, dP = P - GH_Pr;
        double Vl = EM_return.Vl, dvdt = EM_return.kress[0], dvdp = EM_return.kress[1];
        double d2vdtp = EM_return.kress[2], d2vdp2 = EM_return.kress[3];
        double gbase_J = (hl - T*sl) + Vl*dP + dvdt*dT*dP + 0.5*dvdp*dP*dP
                       + 0.5*d2vdtp*dT*dP*dP + (d2vdp2/6.0)*dP*dP*dP;

        return GH_pack_PP_ref(name, len_ox, EM_return.Comp, bulk_rock, apo, gbase_J);
    }

    double Tfus = EM_return.Tfus;
    double Sfus = EM_return.Sfus;
    double Cpl  = EM_return.Cpl;
    double Vl    = EM_return.Vl;
    double dvdt   = EM_return.kress[0];
    double dvdp   = EM_return.kress[1];
    double d2vdtp = EM_return.kress[2];
    double d2vdp2 = EM_return.kress[3];

    /* integrate the SOLID's Berman Cp (+ lambda transition, if any) from Tr to Tfus */
    double H_sol_Tfus, S_sol_Tfus;
    GH_berman_HS(Tfus, EM_return.H, EM_return.S, EM_return.cp_berman, &H_sol_Tfus, &S_sol_Tfus);

    /* fusion correction: deltaG_fus(Tfus) = 0 => deltaH_fus = Tfus*deltaS_fus */
    double H_liq_Tfus = H_sol_Tfus + Tfus*Sfus;
    double S_liq_Tfus = S_sol_Tfus + Sfus;

    /* integrate the LIQUID's constant Cp from Tfus to T */
    double H_liq_T = H_liq_Tfus + Cpl*(T - Tfus);
    double S_liq_T = S_liq_Tfus + Cpl*log(T/Tfus);

    double G_Pr = H_liq_T - T*S_liq_T;

    /* pressure correction: integral of the Kress liquid-volume polynomial
       V(T,P) = Vl + dvdt*(T-Trl) + dvdp*(P-Pr) + d2vdtp*(T-Trl)*(P-Pr) + 0.5*d2vdp2*(P-Pr)^2
       from Pr to P, at fixed T. Trl=1673.0K is real gibbs.c's own liquid-
       specific reference temperature for this EOS ("static const double
       trl = 1673.0;" in gibbs.c) - NOT the general Tr=298.15K used for
       the Berman solid-state integrals above. Confirmed missing 2026-07-
       14: found by calling xMELTS' own gibbs() directly (not re-deriving
       the formula) on each of the 13 liquid endmembers at T=1129.10C/
       P=2kbar (a real, MELTS-confirmed spinel-stable point) and comparing
       - every one of the 11 oxide endmembers (all but H2O/CO2, which use
       a separate Pitzer-Sterner path) was off by 0.04-3.3 kJ, and the
       predicted magnitude from just this one wrong constant (recomputed
       standalone with each endmember's own dvdt/d2vdtp) matched every
       observed diff to <0.001 kJ - see [[gh-spn-liq-gbase-verification]].
       An earlier, purely algebraic/symbolic check of this same formula
       (done in an earlier session) verified the *structure* of the 4-term
       Kress expansion reduces correctly, but never cross-checked that its
       "Tr" numerically matched gibbs.c's own liquid-specific trl constant
       - a purely structural check can miss a wrong constant like this. */
    const double GH_liq_Trl = 1673.0;
    double dT = T - GH_liq_Trl;
    double dP = P - GH_Pr;
    double dG_P;
    if (GH_actual_EM_database == 2){
        /* pMELTS: Birch-Murnaghan replaces the Kress pressure integral
           entirely - see GH_pmelts_liq_BM3_dG's header comment. */
        dG_P = GH_pmelts_liq_BM3_dG(T, P, Vl, dvdt, dvdp, d2vdtp);
    }
    else {
        dG_P =  Vl*dP
              + dvdt*dT*dP
              + 0.5*dvdp*dP*dP
              + 0.5*d2vdtp*dT*dP*dP
              + (d2vdp2/6.0)*dP*dP*dP;
    }

    double gbase_J = G_Pr + dG_P;

    PP_ref PP_ref_db = GH_pack_PP_ref(name, len_ox, EM_return.Comp, bulk_rock, apo, gbase_J);
    PP_ref_db.phase_cp = Cpl/GH_kbar2bar;

    return PP_ref_db;
}
