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
    Objective function for the Stage-A Ghiorso/MELTS liquid ("liq"): a
    13-endmember symmetric regular solution (ideal mole-fraction mixing +
    Margules excess energy), following the same formulation as
    SB_database's solution models (e.g. obj_sb11_ol): n_xeos == n_em,
    p[i] = x[i] directly (no reduced/rotated basis), with the Sigma(p)=1
    closure enforced by an NLopt equality constrain.

    rMELTS only (GH_actual_EM_database==1): also includes real rMELTS'
    embedded liquid speciation reaction CaSiO3 + CO2 <-> SiO2 + CaCO3
    (hidden species s), FIXED 2026-07-14 - see [[gh-gexcess-verification]].
    A first attempt at this same fix (earlier the same day, see
    [[gh-spn-liq-gbase-verification]] point 5) used a piecewise
    boundary-check-then-Newton-loop structure that broke SLSQP's line
    search during pseudocompound refinement (gradient non-smoothness at
    the s=0/s>0 transition). This version instead solves s via a Newton
    iteration that is always evaluated (no discrete branch other than a
    single r9,r13>TINY guard for the genuinely-zero-reactant case, where
    s->0 continuously anyway), matching real rhyolite-MELTS' own order()
    Newton solve in structure (not bisection). Derivation verified exactly
    (to floating-point precision) against real gmixLiq_CO2_H2O across
    multiple T/P/compositions via a from-scratch harness calling the real
    compiled library directly.
*/
#include <math.h>
#include <complex.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#include "../MAGEMin.h"
#include "gh_objective_functions.h"
#include "GH_fluid_eos.h"
#include "GH_gem_function.h"

/** Mirrors gv.gh_multistart_order (obj_gh_spn's own signature is fixed by
    NLopt's callback interface and doesn't receive gv directly) - copied
    once in GH_SS_objective_init_function, which does receive gv. Default 0
    matches real xMELTS' own order() function exactly (single physically-
    motivated starting guess, no multi-start) - see MAGEMin.h. */
static int GH_spn_multistart_flag = 0;

/** cpx's SS2 Taylor coefficient is xMELTS-only in real clinopyroxene.c
    ("#define SS2 ((calculationMode == MODE_xMELTS) ? 0.25*(S027) : 0.0)")
    - set once by GH_SS_objective_init_function per gv.EM_database, read
    by the huge cpx constant table further down. Default matches xMELTS
    (EM_database==0) in case this is ever read before init runs. */
static double GH_cpx_SS2 = -1.08018328;
/** Same role as GH_cpx_SS2 but for opx: real orthopyroxene.c has the
    IDENTICAL "#define SS2 ((calculationMode==MODE_xMELTS)?0.25*(S027):0.0)"
    macro text - only the VALUE of S027 differs (oS027, since opx's own
    file-local "clino" always resolves FALSE), giving 0.25*oS027 =
    -0.57977688 instead of cpx's 0.25*cS027 = -1.08018328. See obj_gh_opx's
    header comment. */
static double GH_opx_SS2 = -0.57977688;

/** rMELTS-only liquid speciation reaction CaSiO3 + CO2 <-> SiO2 + CaCO3
    (hidden species, extent s). W(CaCO3,X) for gh's own 13 rMELTS liq
    endmembers (SiO2,TiO2,Al2O3,Fe2O3,MgCr2O4,Fe2SiO4,MnSi0.5O2,Mg2SiO4,
    CaSiO3,Na2SiO3,KAlSiO4,CO2,H2O in that order - EM_tmp in
    gh_gss_function.c), extracted directly from real
    param_struct_data_CO2_H2O.h's self-labeled "W(CaCO3,X)" entries. */
static const double GH_rmelts_Wcc[13] = {
     63281.2,  /* SiO2     */
    -79202.7,  /* TiO2     */
     46716.1,  /* Al2O3    */
     65508.7,  /* Fe2O3    */
         0.0,  /* MgCr2O4  */
    -72996.8,  /* Fe2SiO4  */
         0.0,  /* MnSi0.5O2*/
    -24872.6,  /* Mg2SiO4  */
     37534.1,  /* CaSiO3   */
   -311011.3,  /* Na2SiO3  */
    -27865.0,  /* KAlSiO4  */
         0.0,  /* CO2      */
      7873.0   /* H2O      */
};
/** CaCO3's own standard-state Gibbs energy is a real xMELTS "phantom
    reaction" (sources/gibbs.c's getCaCO3properties): CaSiO3 + CO2 - SiO2,
    plus a small empirical correction (gibbs.c's hCorr/vCorr constants,
    applied to G only, not S). Verified exact (floating-point precision)
    against real gibbs(t,p,"CaCO3",...) at 4 independent T/P points. P
    here must be in BAR (pr=1 bar reference), not kbar - gh's own d->P is
    in kbar, hence the *1000.0. */
#define GH_RMELTS_CACO3_HCORR (-17574.497522747)
#define GH_RMELTS_CACO3_VCORR (-1.9034060173857)

/** Flat-array index of W(i,k), i<k, for gh's fixed 13-endmember rMELTS
    liq enumeration ("for i: for k=i+1..12: it++") - same convention the
    main Margules loop below uses, extracted as a formula so the CaCO3
    correction can look up arbitrary pairs without re-walking the loop. */
static inline int GH_rmelts_widx(int i, int k){
    return i*13 - i*(i+1)/2 + (k-i-1);
}

/** Hex(s) = sum_{i<k,13 endmembers} W(i,k)*x_i(s)*x_k(s)
           + sum_k W(CaCO3,k)*x_k(s)*s
    where x_SiO2=p[0]+s, x_CaSiO3=p[8]-s, x_CO2=p[11]-s, all other x_i=p[i]
    unchanged, and s is CaCO3's own (hidden-species) mole fraction. Exactly
    quadratic in s since every x_i(s) is affine in s - evaluated directly
    (not via a hand-simplified closed form) to avoid the repeated algebra
    mistakes documented in [[gh-gexcess-verification]]'s derivation notes. */
static double GH_rmelts_Hex(const double *p, const double *W, double s){
    double xa[13];
    for (int i = 0; i < 13; i++) xa[i] = p[i];
    xa[0]  += s;
    xa[8]  -= s;
    xa[11] -= s;
    double Hex = 0.0;
    for (int i = 0; i < 13; i++){
        for (int k = i+1; k < 13; k++){
            Hex += W[GH_rmelts_widx(i,k)]*xa[i]*xa[k];
        }
    }
    for (int k = 0; k < 13; k++){
        Hex += GH_rmelts_Wcc[k]*xa[k]*s;
    }
    return Hex;
}

/** dHex(s)/dp_j at fixed s (envelope-theorem gradient piece - s* is a
    stationary point of Gex over s, so the total derivative of the
    correction w.r.t. p_j equals this explicit partial derivative; no
    implicit differentiation through s*(p) is needed). */
static double GH_rmelts_dHexdp(const double *p, const double *W, double s, int j){
    double xa[13];
    for (int i = 0; i < 13; i++) xa[i] = p[i];
    xa[0]  += s;
    xa[8]  -= s;
    xa[11] -= s;
    double dj = 0.0;
    for (int k = 0; k < 13; k++){
        if (k == j) continue;
        double Wjk = (j < k) ? W[GH_rmelts_widx(j,k)] : W[GH_rmelts_widx(k,j)];
        dj += Wjk*xa[k];
    }
    dj += GH_rmelts_Wcc[j]*s;
    return dj;
}

/** Ideal configurational entropy (sum x_i*log(x_i)) over gh's 13 basis
    endmembers PLUS the hidden CaCO3 species (x_CaCO3=s) - real rMELTS'
    own gmixLiq_CO2_H2O sums entropy over all NE species, not just the
    NA basis ones, when s>0. Does NOT include the separate H2O-specific
    extra term (obj_gh_liq's own Sconfig_h2o, unaffected by s). */
static double GH_rmelts_entropy(const double *p, const double *d_em, double s){
    double xa[13];
    for (int i = 0; i < 13; i++) xa[i] = p[i] + d_em[i];
    xa[0]  += s;
    xa[8]  -= s;
    xa[11] -= s;
    double config = 0.0;
    for (int i = 0; i < 13; i++) if (xa[i] > 0.0) config += xa[i]*creal(clog(xa[i]));
    if (s > 0.0) config += s*creal(clog(s));
    return config;
}

double obj_gh_liq(unsigned n, const double *x, double *grad, void *SS_ref_db){
    SS_ref *d = (SS_ref *) SS_ref_db;

    int n_em    = d->n_em;
    double R    = d->R;
    double T    = d->T;
    double *p   = d->p;
    double *gb  = d->gb_lvl;
    double *mu_Gex = d->mu_Gex;

    for (int i = 0; i < n_em; i++){
        p[i] = x[i];
    }

    /* symmetric regular-solution excess energy (identical pattern to obj_sb11_ol,
       generalized from 2 to n_em endmembers)                                     */
    for (int i = 0; i < n_em; i++){
        double Gex = 0.0;
        int it = 0;
        for (int j = 0; j < n_em; j++){
            double tmp = d->eye[i][j] - p[j];
            for (int k = j+1; k < n_em; k++){
                Gex -= tmp*(d->eye[i][k]-p[k])*(d->W[it]);
                it += 1;
            }
        }
        mu_Gex[i] = Gex/1000.0;   /* J -> kJ */
        d->sf[i]  = p[i];         /* mixing fraction, used only for the generic sf>=0 check */
    }

    d->sum_apep = 0.0;
    for (int i = 0; i < n_em; i++){
        d->sum_apep += d->ape[i]*p[i];
    }
    d->factor = d->fbc/d->sum_apep;

    /* d_em[i]=1 for a boiled-out endmember makes log(p[i]+d_em[i])=log(1)=0
       regardless of p[i]; a non-boiled-out endmember never reaches p[i]=0
       (see this file's header comment) so the log below is unconditional. */
    double Sconfig = 0.0;
    double dSi[13];
    for (int i = 0; i < n_em; i++){
        Sconfig += p[i]*creal(clog(p[i] + d->d_em[i]));
        dSi[i]   = R*T*(creal(clog(p[i] + d->d_em[i]))+1.0);
    }
    Sconfig *= R*T;

    /* Extra H2O-specific ideal-entropy term, found and fixed 2026-07-14
       via [[gh-gexcess-verification]]'s follow-up liq investigation: real
       xMELTS' own gmixLiq_v34 (sources/liquid_v34.c) adds
       "R*t*(x[NA-1]*log(x[NA-1]) + (1-x[NA-1])*log(1-x[NA-1]))" AFTER its
       plain per-component ideal-mixing sum, where x[NA-1] is always H2O
       (the last of the NA liquid components in xMELTS' own ordering) -
       a real, deliberate extra term (not present for any other
       component), missing entirely from this port before now. H2O is
       likewise always the LAST endmember in gh's own EM_list for every
       liq calibration this function serves (xMELTS/rMELTS index 12 of 13,
       pMELTS index 11 of 12 - both confirmed directly against
       G_SS_gh_liq_function's EM_tmp/EM_tmp_pMELTS arrays), so
       p[n_em-1]/d_em[n_em-1] is always H2O here too - safe to use
       generically rather than needing a per-calibration branch. Confirmed
       numerically: this term's magnitude matches the previously-observed
       H2O-sensitive residual almost exactly (e.g. -3.982 kJ predicted vs
       -3.981 kJ observed at x_H2O=0.1, T=1200C) at a pure SiO2-H2O binary
       composition where every other possible source of diff is excluded. */
    double x_h2o = p[n_em-1] + d->d_em[n_em-1];
    /* y*log(y) -> 0 as y -> 0+ (standard convention); guarded because
       H2O being boiled out of the bulk (d_em[n_em-1]=1) pins x_h2o to
       exactly 1.0, which would otherwise hit 0*log(0) = NaN in the
       (1.0-x_h2o) term below - unlike the generic Sconfig loop above,
       this H2O-specific binary-entropy form isn't protected by the
       d_em offset trick at BOTH of its own boundaries. */
    double log_x_h2o   = (x_h2o > 0.0) ? creal(clog(x_h2o))       : 0.0;
    double log_1mx_h2o = (x_h2o < 1.0) ? creal(clog(1.0-x_h2o))   : 0.0;
    double Sconfig_h2o = R*T*(x_h2o*log_x_h2o + (1.0-x_h2o)*log_1mx_h2o);
    dSi[n_em-1] += R*T*(log_x_h2o - log_1mx_h2o);

    /* rMELTS-only: CaSiO3 + CO2 <-> SiO2 + CaCO3 speciation reaction, see
       this function's header comment. A first attempt (same day, see
       [[gh-spn-liq-gbase-verification]] point 5) used a discrete
       "p[CaSiO3]>TINY && p[CO2]>TINY" value-zeroing branch here and broke
       SLSQP's line search during pseudocompound refinement (confirmed by
       direct pre/post-fix comparison at T=1300C/P=8kbar/test=0: pre-fix
       gives a clean 2-instance liq split, that branch gave 6 fragmented
       liq instances) - the branch discarded a real, non-negligible
       correction over a wide range of small-but-nonzero compositions,
       creating a genuine value/gradient jump at the threshold. Fixed here
       by reparametrizing s = min(p[CaSiO3],p[CO2]) * sigma, sigma in
       (0,1): this keeps x_CaSiO3=p[CaSiO3]-s and x_CO2=p[CO2]-s
       nonnegative BY CONSTRUCTION for any sigma in (0,1) (no separate
       epsilon-floor needed), and s->0 EXACTLY (not approximately) and
       CONTINUOUSLY as either reactant's mole fraction ->0, since s is
       their product's multiple - so the only remaining branch below
       (fmin(r9,r13) > 0.0) fires solely at the literal, measure-zero
       exact-zero boundary where the correction's true value already IS
       zero in the limit, not an approximation of it. */
    double dGex_rmelts = 0.0;
    double dSi_rmelts[13] = {0.0};
    if (GH_actual_EM_database == 1 && n_em == 13){
        /* d->W[]/GH_rmelts_Wcc[] are raw J (same convention mu_Gex's own
           /1000.0 conversion already relies on); d->R is kJ-scaled
           (0.0083144) for direct use against df_raw. This whole block
           works in J internally (Rgas = d->R*1000.0, matching e.g.
           obj_gh_fsp's own convention for the same reason) and converts
           only the final correction back to kJ at the end. */
        const int iSi = 0, iCa = 8, iCO2 = 11;
        double Rgas = R*1000.0;
        double Pbar = d->P*1000.0;
        double C = GH_RMELTS_CACO3_HCORR + GH_RMELTS_CACO3_VCORR*(Pbar-1.0);

        double r9 = p[iCa], r13 = p[iCO2], x0 = p[iSi];
        double m = fmin(r9, r13);
        double s = 0.0;

        if (m > 0.0){
            double Hex0  = GH_rmelts_Hex(p, d->W,  0.0);
            double Hex_p = GH_rmelts_Hex(p, d->W,  1.0);
            double Hex_m = GH_rmelts_Hex(p, d->W, -1.0);
            double Bcoef = 0.5*(Hex_p - Hex_m);
            double Dcoef = 0.5*(Hex_p + Hex_m) - Hex0;

            double sigma = 0.5;
            double lo = 1.0e-14, hi = 1.0 - 1.0e-14;
            for (int iter = 0; iter < 100; iter++){
                double sv = m*sigma;
                double xSi = x0+sv, xCa = r9-sv, xCO2v = r13-sv;
                double f   = Bcoef + 2.0*Dcoef*sv + Rgas*T*(creal(clog(xSi))-creal(clog(xCa))-creal(clog(xCO2v))+creal(clog(sv))) + C;
                double dfds = 2.0*Dcoef + Rgas*T*(1.0/xSi + 1.0/xCa + 1.0/xCO2v + 1.0/sv);
                double fp  = dfds*m;   /* df/dsigma = df/ds * ds/dsigma, ds/dsigma = m */
                double sigma_new = sigma - f/fp;
                if (!(sigma_new > lo && sigma_new < hi)) sigma_new = 0.5*(sigma + (f > 0.0 ? lo : hi));
                if (fabs(sigma_new - sigma) < 1.0e-15){ sigma = sigma_new; break; }
                sigma = sigma_new;
            }
            s = m*sigma;

            double Hex_s = Hex0 + Bcoef*s + Dcoef*s*s;
            double ent_s = GH_rmelts_entropy(p, d->d_em, s);
            double ent_0 = GH_rmelts_entropy(p, d->d_em, 0.0);
            dGex_rmelts = ((Hex_s - Hex0) + Rgas*T*(ent_s - ent_0) + s*C) / 1000.0;

            if (grad){
                double xSi_s = x0+s, xCa_s = r9-s, xCO2_s = r13-s;
                for (int j = 0; j < 13; j++){
                    dSi_rmelts[j] = (GH_rmelts_dHexdp(p, d->W, s, j) - GH_rmelts_dHexdp(p, d->W, 0.0, j)) / 1000.0;
                }
                dSi_rmelts[iSi]  += R*T*(creal(clog(xSi_s))  - creal(clog(x0)));
                dSi_rmelts[iCa]  += R*T*(creal(clog(xCa_s))  - creal(clog(r9)));
                dSi_rmelts[iCO2] += R*T*(creal(clog(xCO2_s)) - creal(clog(r13)));
            }
        }
    }

    d->df_raw = 0.0;
    for (int i = 0; i < n_em; i++){
        d->df_raw += (mu_Gex[i] + gb[i])*p[i];
    }
    d->df_raw += Sconfig + Sconfig_h2o + dGex_rmelts;
    d->df = d->df_raw * d->factor;

    if (grad){
        for (int i = 0; i < n_em; i++){
            grad[i] = (mu_Gex[i] + gb[i] + dSi[i] + dSi_rmelts[i])*d->factor
                      - (d->df_raw*d->factor*(d->ape[i]/d->sum_apep));
        }
    }

    return d->df;
}

/**
    Olivine (Fo-Fa): a genuine reduction of xMELTS' full 6-component,
    order-parameter-bearing olivine model (sources/olivine.c: NR=5, NS=4,
    Taylor-expansion in r[]/s[] "Darken equation" coefficients), not an
    ad hoc simplification. Substituting r0=r2=r3=-1 (no Mn/Co), r4=0 (no
    Ca) and evaluating the reduced 2-variable (r1=2*xFa-1, s1=M1/M2 Fe-Mg
    ordering) Taylor polynomial gives:
        H(r1,s1) = C0 + C2*r1^2 + C5*s1^2   (C1=C3=C4=0 identically)
    with C0 = -C2 (verified numerically: H=0 at both pure endmembers) and
    C5 >= 0 (0 at Pr=1bar, growing slowly with P). Since both H and the
    ideal site-mixing entropy S(s1) are even functions of s1 minimized/
    maximized at s1=0 for every composition, the equilibrium order
    parameter is s1*=0 identically across the whole Fo-Fa join (real,
    not assumed: physically consistent with olivine showing negligible
    M1/M2 Fe-Mg ordering, unlike e.g. orthopyroxene) - so no NLopt order-
    parameter dimension is needed here, only the disordered 2-site
    entropy already used, PLUS the corrected effective regular-solution
    energy W(P) = 4*C0(P) = 20300 + 0.015*(Pbar-1) J (NOT the raw
    WH1MGFE/WH2MGFE=5075 J ordering-reaction parameter previously used
    by mistake - that constant is the M1/M2 EXCHANGE energy entering the
    Taylor coefficients, not the bulk disordered-limit regular-solution
    W; verified by direct numerical evaluation of the full xMELTS
    coefficient set reduced to the Fo-Fa limit).
*/
double obj_gh_ol(unsigned n, const double *x, double *grad, void *SS_ref_db){
    SS_ref *d = (SS_ref *) SS_ref_db;

    int n_em    = d->n_em;
    double R    = d->R;
    double T    = d->T;
    double *p   = d->p;
    double *gb  = d->gb_lvl;
    double *mu_Gex = d->mu_Gex;

    /* W(P) already computed in G_SS_gh_ol_function (gh_gss_function.c) */
    double W = d->W[0];

    for (int i = 0; i < n_em; i++){
        p[i] = x[i];
    }

    /* mu_i = W*p_j^2 (j != i), the same eye[][]-trick chemical potential
       used everywhere else, now with the corrected W */
    mu_Gex[0] = W*p[1]*p[1]/1000.0;
    mu_Gex[1] = W*p[0]*p[0]/1000.0;
    d->sf[0] = p[0]; d->sf[1] = p[1];

    d->sum_apep = 0.0;
    for (int i = 0; i < n_em; i++){
        d->sum_apep += d->ape[i]*p[i];
    }
    d->factor = d->fbc/d->sum_apep;

    /* unconditional log() - see obj_gh_liq's header comment */
    double Sconfig = 2.0*R*T*(p[0]*creal(clog(p[0] + d->d_em[0])) + p[1]*creal(clog(p[1] + d->d_em[1])));

    d->df_raw = 0.0;
    for (int i = 0; i < n_em; i++){
        d->df_raw += (mu_Gex[i] + gb[i])*p[i];
    }
    d->df_raw += Sconfig;
    d->df = d->df_raw * d->factor;

    if (grad){
        double dS0 = 2.0*R*T*(creal(clog(p[0] + d->d_em[0]))+1.0);
        double dS1 = 2.0*R*T*(creal(clog(p[1] + d->d_em[1]))+1.0);
        grad[0] = (dS0 + mu_Gex[0] + gb[0])*d->factor - (d->df_raw*d->factor*(d->ape[0]/d->sum_apep));
        grad[1] = (dS1 + mu_Gex[1] + gb[1])*d->factor - (d->df_raw*d->factor*(d->ape[1]/d->sum_apep));
    }
    return d->df;
}

/**
    Real mixed H2O-CO2 fluid (Pitzer & Sterner 1994), from xMELTS'
    sources/fluid.c fluidPhase() - see GH_pitzer_sterner_mix_G's header
    comment (GH_fluid_eos.h) for the actual mixing physics. Unlike every
    other gh phase, there is no stored-W[] Margules Gex here: the mixture
    EOS returns the TOTAL molar Gibbs energy G(x,T,P) directly (a real
    fluid, not a lattice solution). Gex is defined the same way it behaves
    everywhere else in gh - the mixture's G minus the ideal mechanical
    mixture of the two pure-fluid gb[] values (gb[0]=H2O, gb[1]=CO2,
    already at the current T,P via get_em_data in
    G_SS_gh_fluid_function) - then folded into mu_Gex[]/df via the exact
    same NA=2,NR=1 "delta_ik-p_k" construction obj_gh_ol above uses for
    its own (much simpler, W*p0*p1) Gex. dGex/dx0 is analytic via the
    envelope theorem baked into GH_pitzer_sterner_mix_G itself.
*/
double obj_gh_fluid(unsigned n, const double *x, double *grad, void *SS_ref_db){
    SS_ref *d = (SS_ref *) SS_ref_db;

    double T    = d->T;
    double Pbar = d->P*1000.0;
    double *p   = d->p;
    double *gb  = d->gb_lvl;
    double *mu_Gex = d->mu_Gex;

    for (int i = 0; i < 2; i++){ p[i] = x[i]; }
    double x0 = p[0], x1 = p[1];

    /* Duan & Zhang (2006) - the ACTUAL rhyolite-MELTS "fluid" phase model
       (fluidPhase.c's gmixFlu/actFlu engine + "h2oduan"/"co2duan"
       propertiesOfPureH2O/CO2 standard states), ported+verified exact
       (0.000000 J worst diff vs the real library for G AND dG/dx over
       873-1673K x 500-9000bar x full composition range, incl. the 2000-bar
       low-P/high-P coefficient splice). Replaced GH_pitzer_sterner_mix_G
       (+ an h-Ts reference correction) 2026-07-17: that model is only real
       MELTS' LIQUID H2O/CO2 standard-state helper (fluid.c), never its
       fluid phase - its crude linear-c[i] mixing carried a huge unphysical
       Gex (~+25 kJ at x=0.5) that produced a spurious pure-H2O + pure-CO2
       two-fluid split; DZ2006 mixes near-ideally and gives the single
       mixed fluid real rhyolite-MELTS predicts.

       Wiring fixed at the same time: the old code back-computed
       Gex = Gtot - x*gb and then re-added gb, which cancels ALGEBRAICALLY -
       df was the absolute Gtot(x), so the NLopt minimization was blind to
       the rotated hyperplane (Gamma) and always slid to the lowest-absolute-G
       endpoint (pure CO2), leaving raw levelling PCs to patch the LP (the
       "fl pseudocompound that did not minimize" symptom). The old model's
       huge Gex barrier masked this by trapping NLopt on each side of the
       gap. Now the DZ2006 PURE standard states live in gbase[] (set in
       G_SS_gh_fluid_function, so rotate_hyperplane applies Gamma to them
       exactly like every other gh phase) and this objective adds only the
       MIXING chemical potentials muGex_i = R*T*ln(x_i*phi_i/phi_i_pure)
       (actFlu's mu[] exactly). Gradient form is the standard gh pattern -
       exact by Gibbs-Duhem (sum x*dmu = 0 at constant T,P). */
    double muW, muC;
    GH_duan_mix_muGex(x0, T, Pbar, &muW, &muC);
    mu_Gex[0] = muW/1000.0;   /* J -> kJ */
    mu_Gex[1] = muC/1000.0;
    d->sf[0] = p[0]; d->sf[1] = p[1];

    d->sum_apep = 0.0;
    for (int i = 0; i < 2; i++){ d->sum_apep += d->ape[i]*p[i]; }
    d->factor = d->fbc/d->sum_apep;

    d->df_raw = (mu_Gex[0]+gb[0])*p[0] + (mu_Gex[1]+gb[1])*p[1];
    d->df = d->df_raw * d->factor;

    if (grad){
        grad[0] = (mu_Gex[0]+gb[0])*d->factor - (d->df_raw*d->factor*(d->ape[0]/d->sum_apep));
        grad[1] = (mu_Gex[1]+gb[1])*d->factor - (d->df_raw*d->factor*(d->ape[1]/d->sum_apep));
    }
    return d->df;
}

/**
    Biotite (Ann-Phl): identical pattern, single symmetric W term, no
    site-multiplicity factor on the ideal entropy (matches xMELTS'
    sources/biotite.c calibration as-is - see G_SS_gh_bi_init_function).
*/
double obj_gh_bi(unsigned n, const double *x, double *grad, void *SS_ref_db){
    SS_ref *d = (SS_ref *) SS_ref_db;

    int n_em    = d->n_em;
    double R    = d->R;
    double T    = d->T;
    double *p   = d->p;
    double *gb  = d->gb_lvl;
    double *mu_Gex = d->mu_Gex;

    for (int i = 0; i < n_em; i++){
        p[i] = x[i];
    }

    for (int i = 0; i < n_em; i++){
        double Gex = 0.0;
        int it = 0;
        for (int j = 0; j < n_em; j++){
            double tmp = d->eye[i][j] - p[j];
            for (int k = j+1; k < n_em; k++){
                Gex -= tmp*(d->eye[i][k]-p[k])*(d->W[it]);
                it += 1;
            }
        }
        mu_Gex[i] = Gex/1000.0;
        d->sf[i]  = p[i];
    }

    d->sum_apep = 0.0;
    for (int i = 0; i < n_em; i++){
        d->sum_apep += d->ape[i]*p[i];
    }
    d->factor = d->fbc/d->sum_apep;

    /* unconditional log() - see obj_gh_liq's header comment */
    double Sconfig = R*T*(p[0]*creal(clog(p[0] + d->d_em[0])) + p[1]*creal(clog(p[1] + d->d_em[1])));

    d->df_raw = 0.0;
    for (int i = 0; i < n_em; i++){
        d->df_raw += (mu_Gex[i] + gb[i])*p[i];
    }
    d->df_raw += Sconfig;
    d->df = d->df_raw * d->factor;

    if (grad){
        double dS0 = R*T*(creal(clog(p[0] + d->d_em[0]))+1.0);
        double dS1 = R*T*(creal(clog(p[1] + d->d_em[1]))+1.0);
        grad[0] = (dS0 + mu_Gex[0] + gb[0])*d->factor - (d->df_raw*d->factor*(d->ape[0]/d->sum_apep));
        grad[1] = (dS1 + mu_Gex[1] + gb[1])*d->factor - (d->df_raw*d->factor*(d->ape[1]/d->sum_apep));
    }
    return d->df;
}

/**
    Feldspar (Ab-An-San): Elkins & Grove (1990) asymmetric ternary
    Margules, ported from xMELTS' sources/feldspar.c H/S/V macros
    (real published W parameters). Unlike the symmetric phases above,
    Wij != Wji here, so the simple eye[][] trick (which only supports one
    W per pair) cannot represent this - instead the raw partial
    derivatives of H, S_excess and V (each expanded as an explicit
    polynomial in the 3 independent p[i], no variable eliminated) are
    combined into the standard partial-molar transform
    mu_i = Gex + sum_k (delta_ik - p_k) * dGex/dp_k, which is the general
    n-component form of the same identity the eye[][] trick specializes
    to for simple symmetric regular solutions (verified: for Gex=W*x0*x1
    it reduces to exactly W*x1^2, matching obj_gh_ol/obj_sb11_ol). This
    keeps the same "direct p=x + Sigma(p)=1 equality constraint" NLopt
    formulation as every other gh/sb phase - only the excess-energy
    formula itself is richer. xMELTS' own Al-Si/albite-monalbite
    order-disorder correction (sources/gibbs.c) is not reproduced, same
    simplification as these endmembers' own standard states.
*/
double obj_gh_fsp(unsigned n, const double *x, double *grad, void *SS_ref_db){
    SS_ref *d = (SS_ref *) SS_ref_db;

    double R = d->R;
    double T = d->T;
    double Pbar = d->P*1000.0;   /* kbar -> bar, matching W calibration units */
    double *p   = d->p;
    double *gb  = d->gb_lvl;
    double *mu_Gex = d->mu_Gex;

    for (int i = 0; i < 3; i++){ p[i] = x[i]; }
    double xab = p[0], xan = p[1], xor_ = p[2];

    const double whabor=18810.0, wsabor=10.3, wvabor=0.4602;
    const double whorab=27320.0, wsorab=10.3, wvorab=0.3264;
    const double whaban=7924.0,  whanab=0.0;
    const double whanor=38974.0, wvanor=-0.1037;
    const double whoran=40317.0;
    const double whabanor=12545.0, wvabanor=-1.095;

    double KH = 0.5*(whaban+whanab+whabor+whorab+whanor+whoran) + whabanor;
    double KS = 0.5*(wsabor+wsorab);
    double KV = 0.5*(wvabor+wvorab+wvanor) + wvabanor;

    double H  = whaban*xab*xan*xan + whanab*xab*xab*xan + whabor*xab*xor_*xor_
              + whorab*xab*xab*xor_ + whanor*xan*xor_*xor_ + whoran*xan*xan*xor_
              + xab*xan*xor_*KH;
    double Sx = wsabor*xab*xor_*xor_ + wsorab*xab*xab*xor_ + xab*xan*xor_*KS;
    double V  = wvabor*xab*xor_*xor_ + wvorab*xab*xab*xor_ + wvanor*xan*xor_*xor_
              + xab*xan*xor_*KV;

    double dH[3], dS[3], dV[3];
    dH[0] = whaban*xan*xan + 2.0*whanab*xab*xan + whabor*xor_*xor_ + 2.0*whorab*xab*xor_ + xan*xor_*KH;
    dH[1] = 2.0*whaban*xab*xan + whanab*xab*xab + whanor*xor_*xor_ + 2.0*whoran*xan*xor_ + xab*xor_*KH;
    dH[2] = 2.0*whabor*xab*xor_ + whorab*xab*xab + 2.0*whanor*xan*xor_ + whoran*xan*xan + xab*xan*KH;

    dS[0] = wsabor*xor_*xor_ + 2.0*wsorab*xab*xor_ + xan*xor_*KS;
    dS[1] = xab*xor_*KS;
    dS[2] = 2.0*wsabor*xab*xor_ + wsorab*xab*xab + xab*xan*KS;

    dV[0] = 2.0*wvorab*xab*xor_ + wvabor*xor_*xor_ + xan*xor_*KV;
    dV[1] = wvanor*xor_*xor_ + xab*xor_*KV;
    dV[2] = 2.0*wvabor*xab*xor_ + wvorab*xab*xab + 2.0*wvanor*xan*xor_ + xab*xan*KV;

    double Gex = H - T*Sx + (Pbar-1.0)*V;
    double dGex[3];
    for (int k = 0; k < 3; k++){ dGex[k] = dH[k] - T*dS[k] + (Pbar-1.0)*dV[k]; }

    for (int i = 0; i < 3; i++){
        double mu = Gex;
        for (int k = 0; k < 3; k++){
            double delta_ik = (i==k) ? 1.0 : 0.0;
            mu += (delta_ik - p[k])*dGex[k];
        }
        mu_Gex[i] = mu/1000.0;   /* J -> kJ */
        d->sf[i]  = p[i];
    }

    d->sum_apep = 0.0;
    for (int i = 0; i < 3; i++){ d->sum_apep += d->ape[i]*p[i]; }
    d->factor = d->fbc/d->sum_apep;

    /* unconditional log() - see obj_gh_liq's header comment */
    double Sconfig = 0.0;
    double dSi[3];
    for (int i = 0; i < 3; i++){
        Sconfig += p[i]*creal(clog(p[i] + d->d_em[i]));
        dSi[i]   = R*T*(creal(clog(p[i] + d->d_em[i]))+1.0);
    }
    Sconfig *= R*T;

    d->df_raw = 0.0;
    for (int i = 0; i < 3; i++){ d->df_raw += (mu_Gex[i] + gb[i])*p[i]; }
    d->df_raw += Sconfig;
    d->df = d->df_raw * d->factor;

    if (grad){
        for (int i = 0; i < 3; i++){
            grad[i] = (mu_Gex[i] + gb[i] + dSi[i])*d->factor - (d->df_raw*d->factor*(d->ape[i]/d->sum_apep));
        }
    }
    return d->df;
}

/**
    Garnet (Gr-Py-Alm): Berman & Koziol (1991) asymmetric ternary
    Margules, ported from xMELTS' sources/garnet.c H/S/V macros (real
    published W parameters, "1==Grossular, 2==Pyrope, 3==Almandine").
    Same mu-transform approach as obj_gh_fsp above, since Wij != Wji here
    too. Unlike feldspar, garnet's own convention uses a plain P (not
    P-1) in its volume term (see sources/garnet.c's G macro) - kept as-is.
*/
double obj_gh_g(unsigned n, const double *x, double *grad, void *SS_ref_db){
    SS_ref *d = (SS_ref *) SS_ref_db;

    double R = d->R;
    double T = d->T;
    double Pbar = d->P*1000.0;   /* kbar -> bar, matching W calibration units */
    double *p   = d->p;
    double *gb  = d->gb_lvl;
    double *mu_Gex = d->mu_Gex;

    for (int i = 0; i < 3; i++){ p[i] = x[i]; }
    double gr = p[0], py = p[1], al = p[2];

    const double WH112=21560.0, WS112=18.79, WV112=0.10;
    const double WH122=69200.0, WS122=18.79, WV122=0.10;
    const double WH113=20320.0, WS113=5.08,  WV113=0.17;
    const double WH133=2620.0,  WS133=5.08,  WV133=0.09;
    const double WH223=230.0,   WS223=0.0,   WV223=0.01;
    const double WH233=3720.0,  WS233=0.0,   WV233=0.06;
    const double WH123=0.0,     WS123=0.0,   WV123=0.0;

    double KH = 0.5*(WH112+WH122+WH113+WH133+WH223+WH233) + WH123;
    double KS = 0.5*(WS112+WS122+WS113+WS133+WS223+WS233) + WS123;
    double KV = 0.5*(WV112+WV122+WV113+WV133+WV223+WV233) + WV123;

    double H  = WH112*gr*gr*py + WH122*gr*py*py + WH113*gr*gr*al + WH133*gr*al*al
              + WH223*py*py*al + WH233*py*al*al + gr*py*al*KH;
    double Sx = WS112*gr*gr*py + WS122*gr*py*py + WS113*gr*gr*al + WS133*gr*al*al
              + WS223*py*py*al + WS233*py*al*al + gr*py*al*KS;
    double V  = WV112*gr*gr*py + WV122*gr*py*py + WV113*gr*gr*al + WV133*gr*al*al
              + WV223*py*py*al + WV233*py*al*al + gr*py*al*KV;

    double dH[3], dS[3], dV[3];
    dH[0] = 2.0*WH112*gr*py + WH122*py*py + 2.0*WH113*gr*al + WH133*al*al + py*al*KH;
    dH[1] = WH112*gr*gr + 2.0*WH122*gr*py + 2.0*WH223*py*al + WH233*al*al + gr*al*KH;
    dH[2] = WH113*gr*gr + 2.0*WH133*gr*al + WH223*py*py + 2.0*WH233*py*al + gr*py*KH;

    dS[0] = 2.0*WS112*gr*py + WS122*py*py + 2.0*WS113*gr*al + WS133*al*al + py*al*KS;
    dS[1] = WS112*gr*gr + 2.0*WS122*gr*py + 2.0*WS223*py*al + WS233*al*al + gr*al*KS;
    dS[2] = WS113*gr*gr + 2.0*WS133*gr*al + WS223*py*py + 2.0*WS233*py*al + gr*py*KS;

    dV[0] = 2.0*WV112*gr*py + WV122*py*py + 2.0*WV113*gr*al + WV133*al*al + py*al*KV;
    dV[1] = WV112*gr*gr + 2.0*WV122*gr*py + 2.0*WV223*py*al + WV233*al*al + gr*al*KV;
    dV[2] = WV113*gr*gr + 2.0*WV133*gr*al + WV223*py*py + 2.0*WV233*py*al + gr*py*KV;

    double Gex = H - T*Sx + Pbar*V;
    double dGex[3];
    for (int k = 0; k < 3; k++){ dGex[k] = dH[k] - T*dS[k] + Pbar*dV[k]; }

    for (int i = 0; i < 3; i++){
        double mu = Gex;
        for (int k = 0; k < 3; k++){
            double delta_ik = (i==k) ? 1.0 : 0.0;
            mu += (delta_ik - p[k])*dGex[k];
        }
        mu_Gex[i] = mu/1000.0;   /* J -> kJ */
        d->sf[i]  = p[i];
    }

    d->sum_apep = 0.0;
    for (int i = 0; i < 3; i++){ d->sum_apep += d->ape[i]*p[i]; }
    d->factor = d->fbc/d->sum_apep;

    /* unconditional log() - see obj_gh_liq's header comment */
    double Sconfig = 0.0;
    double dSi[3];
    for (int i = 0; i < 3; i++){
        Sconfig += p[i]*creal(clog(p[i] + d->d_em[i]));
        dSi[i]   = 3.0*R*T*(creal(clog(p[i] + d->d_em[i]))+1.0);
    }
    Sconfig *= 3.0*R*T;

    d->df_raw = 0.0;
    for (int i = 0; i < 3; i++){ d->df_raw += (mu_Gex[i] + gb[i])*p[i]; }
    d->df_raw += Sconfig;
    d->df = d->df_raw * d->factor;

    if (grad){
        for (int i = 0; i < 3; i++){
            grad[i] = (mu_Gex[i] + gb[i] + dSi[i])*d->factor - (d->df_raw*d->factor*(d->ape[i]/d->sum_apep));
        }
    }
    return d->df;
}

/**
    Hornblende (Parg-Fparg-Mhst): ported from xMELTS' sources/hornblende.c.
    Unlike fsp/g above, this phase has NO separate ideal Sigma(p)*log(p)
    term added on top of Gex - the site-mixing entropy S below (M12 site,
    multiplicity 4, Mg-Fe2+; M3 site, multiplicity 1, Al-Fe3+) already IS
    the complete configurational entropy of the model (real crystallographic
    sites, not a "treat 3 endmembers as ideally-mixing molecules" fiction),
    so folding it into Gex=H-T*S and using the same
    mu_i=Gex+sum_k(delta_ik-p_k)*dGex/dp_k transform already gives the
    complete df_raw with no extra Sconfig addition (verified: the standard
    partial-molar identity sum_i mu_Gex[i]*p[i]=Gex holds regardless of
    what Gex represents, so no double-counting). H/W are raw J (matching
    fsp/g's own convention); Rgas=d->R*1000 (the true gas constant, not the
    kJ-scaled d->R used for the simpler phases' R*T*p*log(p) term) is used
    so the entropy term combines correctly with the raw-J H/W before the
    final /1000 J->kJ conversion. Only the two non-reference endmembers'
    own site fractions (xFe2M12=p[1], xFe3M3=p[2]) can hit log(0) if
    boiled out (bulk missing FeO or, for mhst, FeO+O) - guarded the same
    way as the simpler phases, via +d_em at the log() call site; their
    complements (xMgM12, xAlM3) are bounds-protected already (never reach
    exactly 1) so need no such guard.
*/
double obj_gh_hb(unsigned n, const double *x, double *grad, void *SS_ref_db){
    SS_ref *d = (SS_ref *) SS_ref_db;

    double T    = d->T;
    double Rgas = d->R*1000.0;
    double *p   = d->p;
    double *gb  = d->gb_lvl;
    double *mu_Gex = d->mu_Gex;

    for (int i = 0; i < 3; i++){ p[i] = x[i]; }

    const double WFEMG = 28116.48, WFEAL = 16780.0;

    double xMgM12  = 1.0-p[1];
    double xFe2M12 = p[1] + d->d_em[1];
    double xAlM3   = 1.0-p[2];
    double xFe3M3  = p[2] + d->d_em[2];

    double H  = WFEMG*p[1] - WFEMG*p[1]*p[1] + WFEAL*p[2] - WFEAL*p[2]*p[2];
    double S  = -Rgas*(4.0*xMgM12*creal(clog(xMgM12)) + 4.0*xFe2M12*creal(clog(xFe2M12)) + xFe3M3*creal(clog(xFe3M3)) + xAlM3*creal(clog(xAlM3)));
    double Gex = H - T*S;   /* V=0 identically for this phase */

    double dGex[3];
    dGex[0] = 0.0;   /* pargasite is the reference endmember (r0,r1 = x_fparg,x_mhst in xMELTS' own basis) */
    dGex[1] = WFEMG*(1.0-2.0*p[1]) + 4.0*Rgas*T*creal(clog(xFe2M12/xMgM12));
    dGex[2] = WFEAL*(1.0-2.0*p[2]) +     Rgas*T*creal(clog(xFe3M3/xAlM3));

    for (int i = 0; i < 3; i++){
        double mu = Gex;
        for (int k = 0; k < 3; k++){
            double delta_ik = (i==k) ? 1.0 : 0.0;
            mu += (delta_ik - p[k])*dGex[k];
        }
        mu_Gex[i] = mu/1000.0;
        d->sf[i]  = p[i];
    }

    d->sum_apep = 0.0;
    for (int i = 0; i < 3; i++){ d->sum_apep += d->ape[i]*p[i]; }
    d->factor = d->fbc/d->sum_apep;

    d->df_raw = 0.0;
    for (int i = 0; i < 3; i++){ d->df_raw += (mu_Gex[i] + gb[i])*p[i]; }
    d->df = d->df_raw * d->factor;

    if (grad){
        for (int i = 0; i < 3; i++){
            grad[i] = (mu_Gex[i] + gb[i])*d->factor - (d->df_raw*d->factor*(d->ape[i]/d->sum_apep));
        }
    }
    return d->df;
}

/**
    Leucite (Lc-Anl-Nlc): ported from xMELTS' sources/leucite.c. Same
    "no separate Sconfig" architecture as obj_gh_hb above - the two-site
    (L: K-Na-H2O; S: K-Na-vacancy, multiplicity 3/2) entropy S below is
    the complete configurational entropy already. p[0]=lc, p[1]=anl are
    xMELTS' own independent r0,r1 (per leucite.c's FR0/FR1 macros testing
    i==0/i==1); p[2]=nlc is the dependent/reference endmember here (note:
    the opposite convention from hornblende, where index 0 was the
    reference - this is a per-source-file choice, not a fixed pattern).
    Known limitation (not yet handled, deferred): unlike obj_gh_hb, the
    site fractions here are PRODUCTS of p[0],p[1] (e.g. xKL=p0*(1-p1)), so
    a simple +d_em at the log() site isn't sufficient by itself to safely
    boil out K2O or H2O absence (p0->0 or p1->0) - not guarded yet, since
    none of gh's current test bulks zero out K2O or H2O; revisit if that
    changes.
*/
double obj_gh_lc(unsigned n, const double *x, double *grad, void *SS_ref_db){
    SS_ref *d = (SS_ref *) SS_ref_db;

    double T    = d->T;
    double Rgas = d->R*1000.0;
    double *p   = d->p;
    double *gb  = d->gb_lvl;
    double *mu_Gex = d->mu_Gex;

    for (int i = 0; i < 3; i++){ p[i] = x[i]; }
    double p0 = p[0], p1 = p[1];

    const double WNAK = 7000.0, WNAH = 7000.0, DGR = 53000.0;

    double a  = p0*(1.0-p1);
    double b  = (1.0-p0)*(1.0-p1);
    double c  = p1;
    double dd = (2.0/3.0)*p0*p1;
    double e  = (2.0/3.0)*(1.0-p0)*p1;
    double f  = 1.0-(2.0/3.0)*p1;

    double K = 1.5*creal(clog(2.0/3.0)) + 0.5*creal(clog(0.5));

    double H = WNAK*p0 + WNAH*p1 - WNAK*p0*p0 - WNAH*p1*p1 + DGR*p0*p1;
    double S = -Rgas*( a*creal(clog(a)) + b*creal(clog(b)) + c*creal(clog(c))
                      + 1.5*(dd*creal(clog(dd)) + e*creal(clog(e)) + f*creal(clog(f)))
                      - p1*K );

    double Gex = H - T*S;   /* V=0 identically for this phase */

    double dH0 = WNAK - 2.0*WNAK*p0 + DGR*p1;
    double dH1 = WNAH - 2.0*WNAH*p1 + DGR*p0;

    double dS0 = -Rgas*( (1.0-p1)*creal(clog(a/b)) + p1*creal(clog(dd/e)) );
    double dS1 = -Rgas*( p0*creal(clog(dd/a)) + (1.0-p0)*creal(clog(e/b)) + creal(clog(c/f)) - K );

    double dGex[3];
    dGex[0] = dH0 - T*dS0;
    dGex[1] = dH1 - T*dS1;
    dGex[2] = 0.0;   /* na-leucite is the dependent/reference endmember here */

    for (int i = 0; i < 3; i++){
        double mu = Gex;
        for (int k = 0; k < 3; k++){
            double delta_ik = (i==k) ? 1.0 : 0.0;
            mu += (delta_ik - p[k])*dGex[k];
        }
        mu_Gex[i] = mu/1000.0;
        d->sf[i]  = p[i];
    }

    d->sum_apep = 0.0;
    for (int i = 0; i < 3; i++){ d->sum_apep += d->ape[i]*p[i]; }
    d->factor = d->fbc/d->sum_apep;

    d->df_raw = 0.0;
    for (int i = 0; i < 3; i++){ d->df_raw += (mu_Gex[i] + gb[i])*p[i]; }
    d->df = d->df_raw * d->factor;

    if (grad){
        for (int i = 0; i < 3; i++){
            grad[i] = (mu_Gex[i] + gb[i])*d->factor - (d->df_raw*d->factor*(d->ape[i]/d->sum_apep));
        }
    }
    return d->df;
}

/**
    Melilite (Ak-Geh-Fak-Na): ported from xMELTS' sources/melilite.c. This
    is gh's first phase with a genuine internal order parameter (NS=1) -
    a NEW pattern relative to every phase above. p[0]=ak (dependent
    endmember, xMELTS' own convention - same as obj_gh_hb, opposite of
    obj_gh_lc), p[1]=geh=r0, p[2]=fak=r1, p[3]=na=r2. s represents Al/Si
    ordering between the T1 and T2 tetrahedral sites within the gehlenite
    component (xAl3T1=r0-s, xAl3T2=(r0+s)/2 - conserves total Al while
    redistributing it between sites).

    s is found by an embedded 1D Newton iteration on dG/ds=0 (using the
    analytic d2G/ds2 as the Newton denominator), run to convergence BEFORE
    computing G/grad for the current trial p[] - i.e. s is treated as
    "solved instantly" at every NLopt evaluation, not as its own NLopt
    dimension. By the envelope theorem this means the composition
    gradient dG/dp_k, evaluated at the converged s, needs NO ds/dp_k
    chain-rule correction (verified directly in xMELTS' own source: the
    analogous SECOND-derivative mask in melilite.c's gmixMel does exactly
    this) - so DGDR0/DGDR1/DGDR2 below are the same partial derivatives
    used elsewhere, just evaluated at s=s*.

    Like obj_gh_hb/obj_gh_lc, there is no separate ideal Sconfig term -
    the 3-site (Ca/Na octahedral, Mg/Fe/Al/Si T1, Al/Si T2) entropy S
    below is already the complete configurational entropy.

    ak/fak/na-melilite (p[0]/r1/r2) all get the same +d_em boiled-out
    guard as obj_gh_hb's xFe2M12/xFe3M3 (see xMg2T1/xNaOct/xFe2T1 below).
    geh (r0) is the one exception, left as a documented limitation: its
    boil-out would force r0=0, which degenerates the T1/T2 order-
    parameter's own valid s-range - a harder case than a simple additive
    guard fixes (see the comment at xCaOct/xNaOct/xFe2T1 below).

    Total molar G = sum(p[i]*gbase[i]) + Gex, and Gex here MUST be
    Gmix(r) = G(r,s*) - ENDMEMBERS(r), matching real MELTS' gmixMel's own
    "*gmix = G; ... *gmix -= ENDMEMBERS;" - a proper Gibbs-energy-of-
    mixing function must vanish at every pure-endmember composition, and
    the raw joint G(r,s*) alone does not (gh's gbase[] does NOT carry any
    per-endmember ordering correction for ak/geh/fak/na-melilite, unlike
    rhm, so unlike obj_gh_rhm this ENDMEMBERS(r) subtraction is the ONLY
    place gehlenite's own T1/T2 Landau ordering correction (real
    melilite.c's GH_G macro - AK_G/FE_AK_G/NA_G are all exactly 0, only
    geh's own correction is nonzero) enters gh's model at all. Missing
    entirely before 2026-07-14 - a real bug, found via
    [[gh-gexcess-verification]] and fixed via the embedded 1D Newton solve
    in GH_mel_mix_geh_order_G below (same DGH_GDS0/D2GH_GDS0S0 formulas as
    real melilite.c's own pureOrder, T,P-only - no r/s dependence, so its
    contribution to dGex is just -geh_G on index 1, geh's own p[] slot).
*/
static double GH_mel_mix_geh_order_G(double T){
    const double DG22p=12000.0, W22p=38354.0;
    double R = 8.3143;
    double s = 0.1;
    for (int iter = 0; iter < 1000; iter++){
        double dgds   = R*T*(-2.0*creal(clog(1.0-s)) + creal(clog(s)) + creal(clog(1.0+s))) + DG22p + W22p*(1.0-2.0*s);
        double d2gds2 = R*T*(2.0/(1.0-s) + 1.0/s + 1.0/(1.0+s)) - 2.0*W22p;
        double sNew   = s - dgds/d2gds2;
        sNew = (sNew > 1.0 - 2.220446049250313e-16) ? 1.0 - 2.220446049250313e-16 : sNew;
        sNew = (sNew < 2.220446049250313e-16) ? 2.220446049250313e-16 : sNew;
        int converged = (fabs(sNew - s) <= 10.0*2.220446049250313e-16);
        s = sNew;
        if (converged) break;
    }
    double S = -R*( 2.0*(1.0-s)*creal(clog(1.0-s)) + s*creal(clog(s)) + (1.0+s)*creal(clog(1.0+s)) - 2.0*creal(clog(2.0)) );
    double H = DG22p*s + W22p*s*(1.0-s);
    return H - T*S;
}

double obj_gh_mel(unsigned n, const double *x, double *grad, void *SS_ref_db){
    SS_ref *d = (SS_ref *) SS_ref_db;

    double T    = d->T;
    double Rgas = d->R*1000.0;
    double *p   = d->p;
    double *gb  = d->gb_lvl;
    double *mu_Gex = d->mu_Gex;

    for (int i = 0; i < 4; i++){ p[i] = x[i]; }
    double r0 = p[1], r1 = p[2], r2 = p[3];
    double xMg2T1 = p[0] + d->d_em[0];   /* = x_ak, gh's own independent variable - NOT recomputed as 1-r0-r1-r2, matching the "p[i]=x[i] direct" convention every other gh phase already uses. +d_em guards log(xMg2T1) if ak is boiled out (bounds_ref pins p[0] to exactly 0) - mirrors obj_gh_hb's xFe2M12/xFe3M3 pattern. */

    const double DG22p=12000.0, W12=9000.0, W12p=13000.0, W13=16000.0, W14=0.0,
                 W22p=38354.0,  W23=9000.0, W2p3=13000.0, W24=9000.0,  W2p4=13000.0, W34=0.0;

    /* --- embedded 1D solve for the order parameter s ---
       Melilite's W22p term produces genuine order-disorder instability
       (negative-curvature regions in s), where plain damped Newton can
       cycle forever between an interior candidate and a clamp boundary
       instead of converging (verified empirically via a standalone FD
       gradient check, scratchpad verify_mel_grad.c). Bisection on the
       sign of dG/ds is unconditionally convergent for a bounded 1D root,
       and naturally falls back to a boundary solution when dG/ds doesn't
       change sign anywhere across the box (G monotonic in s throughout,
       so the constrained optimum is a boundary, not an interior root). */
    double eps_s = 1.0e-10;
    double s_lo = fmax(fmax(r0-(1.0-eps_s), eps_s-r2), 2.0*eps_s-r0);
    double s_hi = fmin(fmin(r0-eps_s, (1.0-eps_s)-r2), 2.0*(1.0-eps_s)-r0);

    double dgds_lo, dgds_hi;
    {
        double xAl3T1 = r0-s_lo, xSi4T1 = r2+s_lo, xAl3T2 = 0.5*(r0+s_lo), xSi4T2 = 1.0-0.5*(r0+s_lo);
        dgds_lo = Rgas*T*(-creal(clog(xAl3T1))+creal(clog(xSi4T1))+creal(clog(xAl3T2))-creal(clog(xSi4T2)))
                + (DG22p+W12p-W12) - 2.0*W22p*s_lo + (W22p-W12p+W12)*r0
                + (W2p3-W23-W12p+W12)*r1 + (W2p4-W24-W12p+W12)*r2;
    }
    {
        double xAl3T1 = r0-s_hi, xSi4T1 = r2+s_hi, xAl3T2 = 0.5*(r0+s_hi), xSi4T2 = 1.0-0.5*(r0+s_hi);
        dgds_hi = Rgas*T*(-creal(clog(xAl3T1))+creal(clog(xSi4T1))+creal(clog(xAl3T2))-creal(clog(xSi4T2)))
                + (DG22p+W12p-W12) - 2.0*W22p*s_hi + (W22p-W12p+W12)*r0
                + (W2p3-W23-W12p+W12)*r1 + (W2p4-W24-W12p+W12)*r2;
    }

    double s;
    if (dgds_lo > 0.0 && dgds_hi > 0.0){
        s = s_lo;
    }
    else if (dgds_lo < 0.0 && dgds_hi < 0.0){
        s = s_hi;
    }
    else {
        double a = s_lo, b = s_hi, fa = dgds_lo;
        for (int it = 0; it < 80; it++){
            double mid = 0.5*(a+b);
            double xAl3T1 = r0-mid, xSi4T1 = r2+mid, xAl3T2 = 0.5*(r0+mid), xSi4T2 = 1.0-0.5*(r0+mid);
            double fm = Rgas*T*(-creal(clog(xAl3T1))+creal(clog(xSi4T1))+creal(clog(xAl3T2))-creal(clog(xSi4T2)))
                      + (DG22p+W12p-W12) - 2.0*W22p*mid + (W22p-W12p+W12)*r0
                      + (W2p3-W23-W12p+W12)*r1 + (W2p4-W24-W12p+W12)*r2;
            if ((fa > 0.0) == (fm > 0.0)){ a = mid; fa = fm; }
            else { b = mid; }
            if (b-a < 1.0e-14){ break; }
        }
        s = 0.5*(a+b);
    }

    /* xNaOct/xFe2T1 get the same +d_em guard as obj_gh_hb's xFe2M12/xFe3M3:
       if na-melilite (index 3) or fak (index 2) is boiled out, bounds_ref
       pins p[3]/p[2] to exactly 0, and without the guard log(xNaOct)/
       log(xFe2T1) would evaluate log(0)=NaN. xCaOct=1-r2 needs no guard
       (it only approaches 0 as r2->1, a legitimate composition, not a
       boiled-out one). geh (index 1, r0) has no such guard here: its
       boil-out would force r0=0, which (see the s_lo/s_hi formulas above)
       makes the T1/T2-site order-parameter's own valid range degenerate/
       empty - a genuinely harder case than a simple additive guard can
       fix, left as a documented limitation like leucite's boil-out gap.
       In practice all of ak/geh/fak/na-melilite's defining oxides (CaO,
       MgO, FeO, Al2O3, SiO2, Na2O) are on gh's "always floored to >=1e-4"
       list (toolkit.c), so none of these boil-out branches currently
       trigger for any reachable bulk composition - this guard is
       defensive/for-correctness, matching the established pattern,
       not a fix for an observed failure. */
    double xCaOct = 1.0-r2, xNaOct = r2 + d->d_em[3], xFe2T1 = r1 + d->d_em[2];
    double xAl3T1 = r0-s, xSi4T1 = r2+s, xAl3T2 = 0.5*(r0+s), xSi4T2 = 1.0-0.5*(r0+s);

    double H = (DG22p+W12p-W12)*s + W12*r0*(1.0-r0) - W22p*s*s + W13*r1*(1.0-r1)
             + (W22p-W12p+W12)*r0*s + (W23-W12-W13)*r0*r1 + (W2p3-W23-W12p+W12)*r1*s
             + W14*r2*(1.0-r2) + (W24-W14-W12)*r0*r2 + (W34-W14-W13)*r1*r2
             + (W2p4-W24-W12p+W12)*r2*s;

    double S = -Rgas*( 2.0*xCaOct*creal(clog(xCaOct)) + 2.0*xNaOct*creal(clog(xNaOct))
                      + xMg2T1*creal(clog(xMg2T1)) + xFe2T1*creal(clog(xFe2T1))
                      + xAl3T1*creal(clog(xAl3T1)) + xSi4T1*creal(clog(xSi4T1))
                      + 2.0*xAl3T2*creal(clog(xAl3T2)) + 2.0*xSi4T2*creal(clog(xSi4T2)) );

    double Graw = H - T*S;   /* V=0 identically for this phase */

    /* ENDMEMBERS(r) = p[1]*geh_G(T) - only gehlenite has a nonzero
       individual ordering correction among ak/geh/fak/na (see this
       function's header comment). */
    double geh_G = GH_mel_mix_geh_order_G(T);
    double Gex = Graw - p[1]*geh_G;

    /* NOTE: xMELTS' own DGDR0/DGDR1/DGDR2 macros (reduced r-basis) carry a
       "-log(xMg2T1)" term, valid there because xMg2T1=1-r0-r1-r2 is
       algebraically DEPENDENT on r0,r1,r2 in that basis. gh instead treats
       xMg2T1=p[0] as its own independent NLopt variable (matching every
       other gh phase's direct p=x convention, Sigma(p)=1 enforced
       externally by NLopt) - so that log() term does not belong in these
       partials, and dropping it removes a hidden "-1" that used to cancel
       a "+1" elsewhere in the reduced-basis derivation; each partial needs
       an explicit "+1" (i.e. +Rgas*T) to compensate. Correspondingly
       dGex[0] is no longer 0: Gex genuinely depends on p[0] through the
       entropy term. Verified against central finite differences (see
       scratchpad verify_mel_grad.c). */
    double dGdr0 = Rgas*T*(creal(clog(xAl3T1))+creal(clog(xAl3T2))-creal(clog(xSi4T2))+1.0)
                 + W12*(1.0-2.0*r0) + (W22p-W12p+W12)*s
                 + (W23-W12-W13)*r1 + (W24-W14-W12)*r2;
    double dGdr1 = Rgas*T*(creal(clog(xFe2T1))+1.0)
                 + W13*(1.0-2.0*r1) + (W23-W12-W13)*r0
                 + (W2p3-W23-W12p+W12)*s + (W34-W14-W13)*r2;
    double dGdr2 = Rgas*T*(-2.0*creal(clog(xCaOct))+2.0*creal(clog(xNaOct))+creal(clog(xSi4T1))+1.0)
                 + W14*(1.0-2.0*r2) + (W24-W14-W12)*r0
                 + (W34-W14-W13)*r1 + (W2p4-W24-W12p+W12)*s;

    double dGex[4];
    dGex[0] = Rgas*T*(creal(clog(xMg2T1))+1.0);
    dGex[1] = dGdr0 - geh_G;   /* -d(ENDMEMBERS)/dp[1], geh's own correction */
    dGex[2] = dGdr1;
    dGex[3] = dGdr2;

    for (int i = 0; i < 4; i++){
        double mu = Gex;
        for (int k = 0; k < 4; k++){
            double delta_ik = (i==k) ? 1.0 : 0.0;
            mu += (delta_ik - p[k])*dGex[k];
        }
        mu_Gex[i] = mu/1000.0;
        d->sf[i]  = p[i];
    }

    d->sum_apep = 0.0;
    for (int i = 0; i < 4; i++){ d->sum_apep += d->ape[i]*p[i]; }
    d->factor = d->fbc/d->sum_apep;

    d->df_raw = 0.0;
    for (int i = 0; i < 4; i++){ d->df_raw += (mu_Gex[i] + gb[i])*p[i]; }
    d->df = d->df_raw * d->factor;

    if (grad){
        for (int i = 0; i < 4; i++){
            grad[i] = (mu_Gex[i] + gb[i])*d->factor - (d->df_raw*d->factor*(d->ape[i]/d->sum_apep));
        }
    }
    return d->df;
}

/**
    Cummingtonite (cumm-grun, Fe-Mg amphibole): ported from xMELTS'
    sources/cummingtonite.c (Ghiorso et al. 1994/1995). NR=1, NS=2 - gh's
    first phase needing a genuine 2D embedded order-parameter solve (M4,
    M1+M3, M2 site Fe-Mg ordering), generalizing melilite's 1D pattern.

    xMELTS' own r-basis (r=2*p[1]-1, i.e. r=-1 at pure cumm, r=+1 at pure
    grun) was converted to gh's direct p[0]=x_cumm, p[1]=x_grun convention
    by substituting r=2*p[1]-1 throughout H(r,s0,s1) and DGDR0 - verified
    against xMELTS' own native r-basis formulas by direct algebraic
    substitution AND numerically (both H and dGex/dp1 match to ~1e-6 via
    central finite differences across p1 in [0.05,0.8], T in [600,1100]K -
    see scratchpad check_cum_basis.c / check_cum_full.c). p[0] (cumm) is
    the dependent endmember (dGex[0]=0), matching every gh phase's
    convention of picking one endmember as the partial-molar reference.

    The embedded solve for s0 (M1+M3 ordering) and s1 (M2 ordering) uses
    a damped 2x2 Newton step: at each iteration, solve the 2x2 linear
    system d2G/ds2 * delta_s = -dG/ds via Cramer's rule, then shrink the
    step length lambda so every one of the 3 independent site fractions
    (xfe2m4, xfe2m13, xfe2m2 - their Mg-complements follow automatically)
    stays within (eps, 1-eps). This step-limiting logic is ported directly
    from xMELTS' own order() function in cummingtonite.c (not invented) -
    unlike melilite's simpler 1D case (where plain Newton could oscillate
    forever near a curvature sign change), this 2D Landau-expansion model
    is well-behaved once step-limited, confirmed convergent and FD-matched
    at every test composition/temperature tried, including near-endmember
    compositions (p1=0.05, p1=0.8).

    By the envelope theorem (verified numerically here as everywhere else
    in gh), the composition gradient dGex/dp1 needs no ds/dp1 correction -
    it is evaluated at the converged s0*,s1* using the same partial
    derivative formula xMELTS itself uses (DGDR0, chain-ruled by the
    constant factor dr/dp1=2 from the r=2*p1-1 substitution).
*/
double obj_gh_cum(unsigned n, const double *x, double *grad, void *SS_ref_db){
    SS_ref *d = (SS_ref *) SS_ref_db;

    double T    = d->T;
    double Rgas = d->R*1000.0;
    double *p   = d->p;
    double *gb  = d->gb_lvl;
    double *mu_Gex = d->mu_Gex;

    for (int i = 0; i < 2; i++){ p[i] = x[i]; }
    double p0 = p[0], p1 = p[1];

    const double HEX=-25856.0, HDEX=5905.0, WH123=5115.64, WH13=4592.11, WH2=4324.28, WH4=3559.52,
                 HX123_4=12102.8, HX134_2=15113.3, HX13_24=3296.10;

    const double H0c    = 5.0*WH123/4.0+WH4/2.0+HX123_4/4.0;
    const double HS13   = 3.0*HEX/5.0+HDEX;
    const double HS2    = 2.0*HEX/5.0-HDEX;
    const double HRS13  = -15.0*WH123/7.0+3.0*WH13-6.0*WH4/7.0-3.0*HX123_4/7.0+HX13_24/2.0;
    const double HRS2   = -10.0*WH123/7.0+2.0*WH2-4.0*WH4/7.0-2.0*HX123_4/7.0+HX134_2/2.0;
    const double HS13S13 = -45.0*WH123/49.0-3.0*WH13/7.0-18.0*WH4/49.0-9.0*HX123_4/49.0+3.0*HX13_24/7.0;
    const double HS13S2  = -60.0*WH123/49.0+12.0*WH13/7.0+12.0*WH2/7.0-24.0*WH4/49.0+25.0*HX123_4/98.0-3.0*HX13_24/14.0-HX134_2/14.0;
    const double HS2S2   = -20.0*WH123/49.0-6.0*WH2/7.0-8.0*WH4/49.0-4.0*HX123_4/49.0+2.0*HX134_2/7.0;

    double r = 2.0*p1-1.0;

    /* --- embedded 2D solve for the order parameters s0 (M1+M3), s1 (M2) ---
       Damped Newton with step-length limiting (ported from xMELTS' own
       cummingtonite.c order() function) rather than plain Newton, keeping
       every site fraction inside (eps,1-eps) at each step. */
    double eps_s = 1.0e-8;
    double s0 = 0.0, s1 = 0.0;
    for (int iter = 0; iter < 200; iter++){
        double xfe2m4  = p1+3.0*s0/7.0+2.0*s1/7.0;
        double xfe2m13 = p1-4.0*s0/7.0+2.0*s1/7.0;
        double xfe2m2  = p1+3.0*s0/7.0-5.0*s1/7.0;
        double xmg2m4 = 1.0-xfe2m4, xmg2m13 = 1.0-xfe2m13, xmg2m2 = 1.0-xfe2m2;

        double dgds0 = Rgas*T*( (6.0/7.0)*creal(clog(xfe2m4/xmg2m4)) - (12.0/7.0)*creal(clog(xfe2m13/xmg2m13)) + (6.0/7.0)*creal(clog(xfe2m2/xmg2m2)) )
                     + HS13 + HRS13*r + 2.0*HS13S13*s0 + HS13S2*s1;
        double dgds1 = Rgas*T*( (4.0/7.0)*creal(clog(xfe2m4/xmg2m4)) + (6.0/7.0)*creal(clog(xfe2m13/xmg2m13)) - (10.0/7.0)*creal(clog(xfe2m2/xmg2m2)) )
                     + HS2 + HRS2*r + HS13S2*s0 + 2.0*HS2S2*s1;

        double H00 = Rgas*T*( (18.0/49.0)/xfe2m4 + (18.0/49.0)/xmg2m4 + (48.0/49.0)/xfe2m13 + (48.0/49.0)/xmg2m13
                            + (18.0/49.0)/xfe2m2 + (18.0/49.0)/xmg2m2 ) + 2.0*HS13S13;
        double H01 = Rgas*T*( (12.0/49.0)/xfe2m4 + (12.0/49.0)/xmg2m4 - (24.0/49.0)/xfe2m13 - (24.0/49.0)/xmg2m13
                            - (30.0/49.0)/xfe2m2 - (30.0/49.0)/xmg2m2 ) + HS13S2;
        double H11 = Rgas*T*( (8.0/49.0)/xfe2m4 + (8.0/49.0)/xmg2m4 + (12.0/49.0)/xfe2m13 + (12.0/49.0)/xmg2m13
                            + (50.0/49.0)/xfe2m2 + (50.0/49.0)/xmg2m2 ) + 2.0*HS2S2;

        double det = H00*H11-H01*H01;
        double ds0 = -( H11*dgds0-H01*dgds1)/det;
        double ds1 = -(-H01*dgds0+H00*dgds1)/det;

        double dxfe2m4  = 3.0*ds0/7.0+2.0*ds1/7.0;
        double dxfe2m13 = -4.0*ds0/7.0+2.0*ds1/7.0;
        double dxfe2m2  = 3.0*ds0/7.0-5.0*ds1/7.0;
        double lambda = 1.0;
        if (xfe2m4+lambda*dxfe2m4 < eps_s){ lambda = fmin(lambda,(eps_s-xfe2m4)/dxfe2m4); }
        if (xfe2m4+lambda*dxfe2m4 > 1.0-eps_s){ lambda = fmin(lambda,(1.0-eps_s-xfe2m4)/dxfe2m4); }
        if (xfe2m13+lambda*dxfe2m13 < eps_s){ lambda = fmin(lambda,(eps_s-xfe2m13)/dxfe2m13); }
        if (xfe2m13+lambda*dxfe2m13 > 1.0-eps_s){ lambda = fmin(lambda,(1.0-eps_s-xfe2m13)/dxfe2m13); }
        if (xfe2m2+lambda*dxfe2m2 < eps_s){ lambda = fmin(lambda,(eps_s-xfe2m2)/dxfe2m2); }
        if (xfe2m2+lambda*dxfe2m2 > 1.0-eps_s){ lambda = fmin(lambda,(1.0-eps_s-xfe2m2)/dxfe2m2); }
        if (lambda < 0.0){ lambda = 0.0; }
        if (lambda > 1.0){ lambda = 1.0; }

        double s0n = s0 + lambda*ds0, s1n = s1 + lambda*ds1;
        if (fabs(s0n-s0) < 1.0e-14 && fabs(s1n-s1) < 1.0e-14){ s0 = s0n; s1 = s1n; break; }
        s0 = s0n; s1 = s1n;
    }

    double xfe2m4  = p1+3.0*s0/7.0+2.0*s1/7.0;
    double xfe2m13 = p1-4.0*s0/7.0+2.0*s1/7.0;
    double xfe2m2  = p1+3.0*s0/7.0-5.0*s1/7.0;
    double xmg2m4 = 1.0-xfe2m4, xmg2m13 = 1.0-xfe2m13, xmg2m2 = 1.0-xfe2m2;

    double H = 4.0*H0c*p0*p1 + (HS13-HRS13)*s0 + (HS2-HRS2)*s1
             + 2.0*HRS13*p1*s0 + 2.0*HRS2*p1*s1
             + HS13S13*s0*s0 + HS13S2*s0*s1 + HS2S2*s1*s1;

    double S = -Rgas*( 2.0*xfe2m4*creal(clog(xfe2m4)) + 2.0*xmg2m4*creal(clog(xmg2m4))
                      + 3.0*xfe2m13*creal(clog(xfe2m13)) + 3.0*xmg2m13*creal(clog(xmg2m13))
                      + 2.0*xfe2m2*creal(clog(xfe2m2)) + 2.0*xmg2m2*creal(clog(xmg2m2)) );

    double Gex = H - T*S;   /* V=0 identically for this phase */

    double dGdp1 = 4.0*H0c*(p0-p1) + 2.0*HRS13*s0 + 2.0*HRS2*s1
                 + Rgas*T*( 2.0*creal(clog(xfe2m4/xmg2m4)) + 3.0*creal(clog(xfe2m13/xmg2m13)) + 2.0*creal(clog(xfe2m2/xmg2m2)) );

    double dGex[2];
    dGex[0] = 0.0;   /* cummingtonite is the dependent endmember */
    dGex[1] = dGdp1;

    for (int i = 0; i < 2; i++){
        double mu = Gex;
        for (int k = 0; k < 2; k++){
            double delta_ik = (i==k) ? 1.0 : 0.0;
            mu += (delta_ik - p[k])*dGex[k];
        }
        mu_Gex[i] = mu/1000.0;
        d->sf[i]  = p[i];
    }

    d->sum_apep = 0.0;
    for (int i = 0; i < 2; i++){ d->sum_apep += d->ape[i]*p[i]; }
    d->factor = d->fbc/d->sum_apep;

    d->df_raw = 0.0;
    for (int i = 0; i < 2; i++){ d->df_raw += (mu_Gex[i] + gb[i])*p[i]; }
    d->df = d->df_raw * d->factor;

    if (grad){
        for (int i = 0; i < 2; i++){
            grad[i] = (mu_Gex[i] + gb[i])*d->factor - (d->df_raw*d->factor*(d->ape[i]/d->sum_apep));
        }
    }
    return d->df;
}

/**
    Spinel (Chr-Herc-Mt-Spl-Usp): from xMELTS' sources/spinel.c (Sack &
    Ghiorso 1991, Contributions to Mineralogy and Petrology 106:474-505;
    also American Mineralogist). NR=4, NS=3 - gh's first 3D embedded
    order-parameter solve (tetrahedral/octahedral Mg-Fe2+-Al-Fe3+
    distribution across two sites; Cr3+ and Ti4+ are always fully
    octahedral by construction, no order parameter needed for them).

    p[0]=chromite(=r1), p[2]=magnetite(=r3), p[3]=spinel(=r0),
    p[4]=ulvospinel(=r2) are DIRECTLY xMELTS' own r0..r3 (a pure relabeling,
    not a basis transform - unlike cummingtonite, xMELTS' own r-basis here
    already IS 4 of the 5 endmembers' mole fractions). p[1]=hercynite is
    the dependent endmember (dGex[1]=0), matching xMELTS' own x[1]=1-r0-r1-r2-r3
    convention. Verified against xMELTS' own H, S, DGDR and D2GDS-S formulas by
    direct substitution AND numerically (dGex/dp matches central finite
    differences to ~1e-6 relative, 60+ compositions incl. near-endmember
    corners - see scratchpad check_spn_full.c/check_spn_gs.c).

    Chromite (p[0]=r1) plays a dual role in xMELTS' own model: it is both
    an independent composition variable AND (via "s3=x3=r1", not a genuine
    free order parameter) contributes extra Taylor terms (gs3*r1 etc.)
    that would otherwise belong to a 4th order parameter. Since these are
    just additional r1-dependent terms in H (not a new degree of freedom),
    they fall out automatically from differentiating the H formula below
    w.r.t. p[0] in the ordinary way - no special handling needed.

    The order-parameter solve uses a joint damped 3D Newton step (all
    three s's updated together each iteration via the real d2G/ds2
    Hessian, one shared step-length damping factor computed from all 8
    site fractions' feasibility bounds) - ported directly from xMELTS'
    own order() function (GH_spn_solve_s_from), not the per-variable
    Gauss-Seidel/bisection scheme an earlier version of this port used.
    That Gauss-Seidel substitute was numerically robust in isolation
    (FD-matched across 60+ compositions) but introduced small solver-
    path-dependent inconsistencies between the returned G and its
    analytic gradient for very close-by trial compositions, enough to
    defeat NLopt's gradient-based SLSQP: diagnosed via
    nlopt_get_numevals/status logging showing ~90% of spn's local
    optimizations silently hitting the 1024-eval cap without converging
    (a ~30x slowdown, invisible without instrumentation - the run still
    reported Status 0 with a plausible-looking G). The real joint Newton
    fixes this (same diagnostic: most calls now reach FTOL_REACHED) and
    is more faithful, including inheriting MELTS' own known limitation at
    extreme near-single-endmember compositions (can still fail to
    converge/land on a different local minimum there) - per this
    project's fidelity principle, reproducing that is correct, not a
    regression to engineer around with different-physics robustness code.
*/
static void GH_narrow_bound(double val0, double slope, double eps, double *lo, double *hi){
    if (slope > 0.0){
        double a = (eps-val0)/slope, b = (1.0-eps-val0)/slope;
        if (a > *lo){ *lo = a; }
        if (b < *hi){ *hi = b; }
    } else if (slope < 0.0){
        double a = (1.0-eps-val0)/slope, b = (eps-val0)/slope;
        if (a > *lo){ *lo = a; }
        if (b < *hi){ *hi = b; }
    }
}
/* Site fractions from the order parameters, clamped into (eps,1-eps) -
   matching xMELTS' own order() guard (spinel.c lines 2221-2241) exactly:
   clamp AFTER computing the raw affine value, BEFORE using it in any
   log()/1/x term, every time it's evaluated (not once at the end). */
static void GH_spn_site_fracs(double p0,double p2,double p3,double p4,double s0,double s1,double s2,double d_em2,
    double *xmg2tet,double *xfe2tet,double *xal3tet,double *xfe3tet,
    double *xmg2oct,double *xfe2oct,double *xal3oct,double *xfe3oct){
    double eps = 2.220446049250313e-16; /* DBL_EPSILON, matching xMELTS' own guard exactly */
    double v;
    v = (p3+s0)/2.0;                                        *xmg2tet = (v<=0.0)?eps:((v>=1.0)?1.0-eps:v);
    v = p4-0.5*p3-0.5*s0+s1+p0+s2;                           *xfe2tet = (v<=0.0)?eps:((v>=1.0)?1.0-eps:v);
    v = 1.0-p0-p4-p2-s1;                                     *xal3tet = (v<=0.0)?eps:((v>=1.0)?1.0-eps:v);
    v = p2-s2+d_em2;                                         *xfe3tet = (v<=0.0)?eps:((v>=1.0)?1.0-eps:v);
    v = (p3-s0)/4.0;                                         *xmg2oct = (v<=0.0)?eps:((v>=1.0)?1.0-eps:v);
    v = (2.0-p3+s0-2.0*s1-2.0*p0-2.0*s2)/4.0;                *xfe2oct = (v<=0.0)?eps:((v>=1.0)?1.0-eps:v);
    v = (1.0-p0-p4-p2+s1)/2.0;                               *xal3oct = (v<=0.0)?eps:((v>=1.0)?1.0-eps:v);
    v = (p2+s2)/2.0+d_em2;                                   *xfe3oct = (v<=0.0)?eps:((v>=1.0)?1.0-eps:v);
}
static void GH_spn_dgds(double p0,double p2,double p3,double p4,double s0,double s1,double s2,double T,double Rgas,
                         const double *gc, double d_em2, double *dgds0,double *dgds1,double *dgds2){
    /* gc[] packs: gs1,gs2,gs3,hs4,ss4,gx2s1,gx2s2,gx2s3,gx2s4,gx3s1,gx3s2,gx3s3,gx3s4,
       gx4s1,gx4s2,gx4s3,gx4s4,gx5s1,gx5s2,gx5s3,gx5s4,gs1s1,gs1s2,gs1s3,gs1s4,
       gs2s2,gs2s3,gs2s4,gs3s3,gs3s4,gs4s4  (31 entries, indices below) */
    double gs1=gc[0],gs2=gc[1],hs4=gc[3],ss4=gc[4];
    double gx2s1=gc[5],gx2s2=gc[6],gx2s4=gc[8];
    double gx3s1=gc[9],gx3s2=gc[10],gx3s4=gc[12];
    double gx4s1=gc[13],gx4s2=gc[14],gx4s4=gc[16];
    double gx5s1=gc[17],gx5s2=gc[18],gx5s4=gc[20];
    double gs1s1=gc[21],gs1s2=gc[22],gs1s3=gc[23],gs1s4=gc[24];
    double gs2s2=gc[25],gs2s3=gc[26],gs2s4=gc[27],gs3s4=gc[29],gs4s4=gc[30];

    double xmg2tet,xfe2tet,xal3tet,xfe3tet,xmg2oct,xfe2oct,xal3oct,xfe3oct;
    GH_spn_site_fracs(p0,p2,p3,p4,s0,s1,s2,d_em2,&xmg2tet,&xfe2tet,&xal3tet,&xfe3tet,&xmg2oct,&xfe2oct,&xal3oct,&xfe3oct);

    *dgds0 = gs1 + gx2s1*p3 + gx3s1*p0 + gx4s1*p4 + gx5s1*p2
           + 2.0*gs1s1*s0 + gs1s2*s1 + gs1s3*p0 + gs1s4*s2
           + 0.5*Rgas*T*(creal(clog(xmg2tet/xfe2tet)) + creal(clog(xfe2oct/xmg2oct)));
    *dgds1 = gs2 + gx2s2*p3 + gx3s2*p0 + gx4s2*p4 + gx5s2*p2
           + gs1s2*s0 + 2.0*gs2s2*s1 + gs2s3*p0 + gs2s4*s2
           + Rgas*T*(creal(clog(xfe2tet/xal3tet)) + creal(clog(xal3oct/xfe2oct)));
    *dgds2 = hs4 - T*ss4 + gx2s4*p3 + gx3s4*p0 + gx4s4*p4 + gx5s4*p2
           + gs1s4*s0 + gs2s4*s1 + gs3s4*p0 + 2.0*gs4s4*s2
           + Rgas*T*(creal(clog(xfe2tet/xfe3tet)) + creal(clog(xfe3oct/xfe2oct)));
}
/* Hessian d2G/ds_i ds_j, from xMELTS' own D2GDS0S0/D2GDS0S1/D2GDS0S2/
   D2GDS1S1/D2GDS1S2/D2GDS2S2 macros (spinel.c lines 1660-1674) - not
   previously ported (gh's spn used a bisection scheme needing only first
   derivatives); needed now for the real joint Newton step below. */
static void GH_spn_d2gds(double p0,double p2,double p3,double p4,double s0,double s1,double s2,double T,double Rgas,
                          const double *gc, double d_em2,
                          double *d00,double *d01,double *d02,double *d11,double *d12,double *d22){
    double gs1s1=gc[21],gs1s2=gc[22],gs1s4=gc[24],gs2s2=gc[25],gs2s4=gc[27],gs4s4=gc[30];
    double xmg2tet,xfe2tet,xal3tet,xfe3tet,xmg2oct,xfe2oct,xal3oct,xfe3oct;
    GH_spn_site_fracs(p0,p2,p3,p4,s0,s1,s2,d_em2,&xmg2tet,&xfe2tet,&xal3tet,&xfe3tet,&xmg2oct,&xfe2oct,&xal3oct,&xfe3oct);

    *d00 = 2.0*gs1s1 + 0.25*Rgas*T*(1.0/xmg2tet + 1.0/xfe2tet + 0.5/xfe2oct + 0.5/xmg2oct);
    *d01 = gs1s2 + 0.5*Rgas*T*(-1.0/xfe2tet - 0.5/xfe2oct);
    *d02 = gs1s4 + 0.5*Rgas*T*(-1.0/xfe2tet - 0.5/xfe2oct);
    *d11 = 2.0*gs2s2 + Rgas*T*(1.0/xfe2tet + 1.0/xal3tet + 0.5/xal3oct + 0.5/xfe2oct);
    *d12 = gs2s4 + Rgas*T*(1.0/xfe2tet + 0.5/xfe2oct);
    *d22 = 2.0*gs4s4 + Rgas*T*(1.0/xfe2tet + 1.0/xfe3tet + 0.5/xfe3oct + 0.5/xfe2oct);
}
static double GH_spn_G_at(double p0,double p2,double p3,double p4,double s0,double s1,double s2,double T,double Rgas,
                           const double *gc,
                           double gx2x2,double gx2x3,double gx2x4,double gx2x5,
                           double gx3x3,double gx3x4,double gx3x5,double gx4x4,double gx4x5,double gx5x5,
                           double WV1,double WV2,double Pval,double d_em0,double d_em2,double d_em4){
    double xmg2tet=(p3+s0)/2.0, xfe2tet=p4-0.5*p3-0.5*s0+s1+p0+s2, xal3tet=1.0-p0-p4-p2-s1, xfe3tet=p2-s2+d_em2;
    double xmg2oct=(p3-s0)/4.0, xfe2oct=(2.0-p3+s0-2.0*s1-2.0*p0-2.0*s2)/4.0, xal3oct=(1.0-p0-p4-p2+s1)/2.0, xfe3oct=(p2+s2)/2.0+d_em2;
    double xcr3oct=p0+d_em0, xti4oct=(p4+d_em4)/2.0;

    double H = gx2x2*p3*p3 + gx2x3*p3*p0 + gx2x4*p3*p4 + gx2x5*p3*p2
             + gc[5]*p3*s0 + gc[6]*p3*s1 + gc[7]*p3*p0 + gc[8]*p3*s2
             + gx3x3*p0*p0 + gx3x4*p0*p4 + gx3x5*p0*p2 + gc[9]*p0*s0 + gc[10]*p0*s1
             + gc[11]*p0*p0 + gc[12]*p0*s2
             + gx4x4*p4*p4 + gx4x5*p4*p2 + gc[13]*p4*s0 + gc[14]*p4*s1 + gc[15]*p4*p0 + gc[16]*p4*s2
             + gx5x5*p2*p2 + gc[17]*p2*s0 + gc[18]*p2*s1 + gc[19]*p2*p0 + gc[20]*p2*s2
             + gc[0]*s0 + gc[1]*s1 + gc[2]*p0 + gc[3]*s2
             + gc[21]*s0*s0 + gc[22]*s0*s1 + gc[23]*s0*p0 + gc[24]*s0*s2
             + gc[25]*s1*s1 + gc[26]*s1*p0 + gc[27]*s1*s2
             + gc[28]*p0*p0 + gc[29]*p0*s2 + gc[30]*s2*s2;

    double S = -Rgas*( xmg2tet*creal(clog(xmg2tet)) + xfe2tet*creal(clog(xfe2tet))
                      + xal3tet*creal(clog(xal3tet)) + xfe3tet*creal(clog(xfe3tet))
                      + 2.0*xmg2oct*creal(clog(xmg2oct)) + 2.0*xfe2oct*creal(clog(xfe2oct))
                      + 2.0*xal3oct*creal(clog(xal3oct)) + 2.0*xfe3oct*creal(clog(xfe3oct))
                      + 2.0*xcr3oct*creal(clog(xcr3oct)) + 2.0*xti4oct*creal(clog(xti4oct)) ) + gc[4]*s2;

    double Gvol = p4*p2*(WV1*p2+WV2*p4)*(Pval-1.0);
    return H - T*S + Gvol;
}
/* Joint damped Newton solve for s0,s1,s2 - ported directly from xMELTS'
   own order() function (spinel.c lines 2203-2325), not a gh-invented
   scheme: at each iteration, solve the 3x3 linear system
   d2gds2 * delta = -dgds (Cramer's rule instead of xMELTS' gaussj, same
   math) for a joint step in all three order parameters together, then
   compute ONE shared step-length damping factor lambda (starting at 1.0,
   tightened in turn by each of the 8 site fractions' own linear
   feasibility bound - the exact sequence and arithmetic xMELTS uses, not
   a "min over all constraints" rewrite) and apply lambda*delta to all
   three simultaneously. This replaces an earlier per-variable Gauss-
   Seidel/bisection scheme that, while numerically robust in isolation,
   produced small solver-path-dependent inconsistencies between the
   returned G and its analytic gradient for very close-by trial
   compositions - enough to defeat NLopt's gradient-based SLSQP, which
   hit its 1024-eval cap on ~90% of calls without converging (diagnosed
   via nlopt_get_numevals/status logging). The real joint Newton matches
   what MELTS itself actually computes, including its own known limits
   (documented: can fail to converge or land on a different local minimum
   at extreme near-single-endmember compositions) - per this project's
   own fidelity principle, that is the correct behavior to reproduce, not
   something to engineer around with different-physics robustness code. */
static void GH_spn_solve_s_from(double p0,double p2,double p3,double p4,double T,double Rgas,const double *gc,double d_em2,
                                 double s0init,double s1init,double s2init,double *s0o,double *s1o,double *s2o){
    const int maxIter = 200; /* xMELTS' own MAX_ITER (spinel.c) */
    double s0=s0init, s1=s1init, s2=s2init;
    double s0Old=2.0, s1Old=2.0, s2Old=2.0; /* force >=1 iteration, matching xMELTS' sOld[i]=2.0 init */
    int iter = 0;
    while ( (fabs(s0-s0Old) > 10.0*2.220446049250313e-16 ||
             fabs(s1-s1Old) > 10.0*2.220446049250313e-16 ||
             fabs(s2-s2Old) > 10.0*2.220446049250313e-16) && iter < maxIter ){
        double xmg2tet,xfe2tet,xal3tet,xfe3tet,xmg2oct,xfe2oct,xal3oct,xfe3oct;
        GH_spn_site_fracs(p0,p2,p3,p4,s0,s1,s2,d_em2,&xmg2tet,&xfe2tet,&xal3tet,&xfe3tet,&xmg2oct,&xfe2oct,&xal3oct,&xfe3oct);

        double dgds0,dgds1,dgds2;
        GH_spn_dgds(p0,p2,p3,p4,s0,s1,s2,T,Rgas,gc,d_em2,&dgds0,&dgds1,&dgds2);
        double d00,d01,d02,d11,d12,d22;
        GH_spn_d2gds(p0,p2,p3,p4,s0,s1,s2,T,Rgas,gc,d_em2,&d00,&d01,&d02,&d11,&d12,&d22);

        s0Old=s0; s1Old=s1; s2Old=s2;

        /* solve [[d00,d01,d02],[d01,d11,d12],[d02,d12,d22]]*delta = -[dgds0,dgds1,dgds2] */
        double r0=-dgds0, r1=-dgds1, r2=-dgds2;
        double det  = d00*(d11*d22-d12*d12) - d01*(d01*d22-d12*d02) + d02*(d01*d12-d11*d02);
        double det0 =  r0*(d11*d22-d12*d12) - d01*( r1*d22-d12* r2) + d02*( r1*d12-d11* r2);
        double det1 = d00*( r1*d22- r2*d12) -  r0*(d01*d22-d12*d02) + d02*(d01* r2- r1*d02);
        double det2 = d00*(d11* r2-d12* r1) - d01*(d01* r2-d12* r0) +  r0*(d01*d12-d11*d02);
        double deltaS0 = det0/det, deltaS1 = det1/det, deltaS2 = det2/det;

        double lambda = 1.0;
        double val, slope;
        /* xmg2tet = (p3+s0)/2 -> d/dstep = deltaS0/2 */
        val=xmg2tet; slope=deltaS0/2.0;
        if      (slope != 0.0 && val+lambda*slope < 0.0) lambda = -val/slope;
        else if (slope != 0.0 && val+lambda*slope > 1.0) lambda = (1.0-val)/slope;
        /* xfe2tet: -deltaS0/2+deltaS1+deltaS2 */
        val=xfe2tet; slope=-deltaS0/2.0+deltaS1+deltaS2;
        if      (slope != 0.0 && val+lambda*slope < 0.0) lambda = -val/slope;
        else if (slope != 0.0 && val+lambda*slope > 1.0) lambda = (1.0-val)/slope;
        /* xal3tet: -deltaS1 */
        val=xal3tet; slope=-deltaS1;
        if      (slope != 0.0 && val+lambda*slope < 0.0) lambda = -val/slope;
        else if (slope != 0.0 && val+lambda*slope > 1.0) lambda = (1.0-val)/slope;
        /* xfe3tet: -deltaS2 */
        val=xfe3tet; slope=-deltaS2;
        if      (slope != 0.0 && val+lambda*slope < 0.0) lambda = -val/slope;
        else if (slope != 0.0 && val+lambda*slope > 1.0) lambda = (1.0-val)/slope;
        /* xmg2oct: -deltaS0/4 */
        val=xmg2oct; slope=-deltaS0/4.0;
        if      (slope != 0.0 && val+lambda*slope < 0.0) lambda = -val/slope;
        else if (slope != 0.0 && val+lambda*slope > 1.0) lambda = (1.0-val)/slope;
        /* xfe2oct: deltaS0/4-deltaS1/2-deltaS2/2 */
        val=xfe2oct; slope=deltaS0/4.0-deltaS1/2.0-deltaS2/2.0;
        if      (slope != 0.0 && val+lambda*slope < 0.0) lambda = -val/slope;
        else if (slope != 0.0 && val+lambda*slope > 1.0) lambda = (1.0-val)/slope;
        /* xal3oct: deltaS1/2 */
        val=xal3oct; slope=deltaS1/2.0;
        if      (slope != 0.0 && val+lambda*slope < 0.0) lambda = -val/slope;
        else if (slope != 0.0 && val+lambda*slope > 1.0) lambda = (1.0-val)/slope;
        /* xfe3oct: deltaS2/2 */
        val=xfe3oct; slope=deltaS2/2.0;
        if      (slope != 0.0 && val+lambda*slope < 0.0) lambda = -val/slope;
        else if (slope != 0.0 && val+lambda*slope > 1.0) lambda = (1.0-val)/slope;

        if (lambda < 1.0){
            s0 = s0Old + lambda*deltaS0; s1 = s1Old + lambda*deltaS1; s2 = s2Old + lambda*deltaS2;
        } else {
            s0 = s0Old + deltaS0; s1 = s1Old + deltaS1; s2 = s2Old + deltaS2;
        }
        iter++;
    }
    *s0o=s0; *s1o=s1; *s2o=s2;
}
/** Spinel's own 5 individual pure-endmember ordering corrections
    (CR_G/HC_G/MT_G/SP_G/UV_G in real spinel.c - chromite/hercynite/
    magnetite/spinel/ulvospinel, matching gh's own p[]/EM_list order
    directly), ported for the same reason and via the same file-local-
    duplication convention as GH_rhm_mix_endmember_order_G/
    GH_mel_mix_geh_order_G above: obj_gh_spn's gbase[] carries NO
    per-endmember ordering correction at all (confirmed via
    GH_gem_function.c's own gbase construction - spn's 5 endmembers get
    only the plain Berman EOS, no correction branch), so ENDMEMBERS(r) =
    sum(p_i * that endmember's own correction) is the ONLY place these 5
    corrections enter gh's model, exactly like melilite's gehlenite. All
    5 corrections are T-only (P-independent - DCR_GDP=DHC_GDP=DMT_GDP=
    DSP_GDP=DUV_GDP=0 in the source). CR_G/UV_G are closed-form (no order
    parameter); SP_G/HC_G/MT_G each have their OWN independent 1D Newton-
    solved order parameter (s0/s1/s2 respectively - 3 SEPARATE
    single-variable problems, not a joint 3x3 solve like
    GH_spn_solve_s_from below, which is for the mixed-composition Gex
    itself, a genuinely different problem) - ported from real spinel.c's
    own separate "pureOrder" function. Missing entirely before
    2026-07-14 - a real bug (the largest of the 3 found in
    [[gh-gexcess-verification]], ~30-90 kJ, since these 5 corrections are
    all tens of kJ, unlike rhm's sub-kJ ones), fixed by subtracting
    ENDMEMBERS(r) in obj_gh_spn below, same pattern as rhm/mel. */
static void GH_spn_mix_endmember_G(double T, double *CR_G, double *HC_G, double *MT_G, double *SP_G, double *UV_G){
    const double H11=-8.7e3*4.184, W11=4.5e3*4.184, W14=20.8e3*4.184,
                 W15=10.0e3*4.184, W15P=11.7e3*4.184,
                 W22=3.6e3*4.184, H24=6.55e3*4.184, W24U=12.6e3*4.184, W2P4U=10.9e3*4.184,
                 H55=6.25e3*4.184, S55=0.0, W55=0.0,
                 HEX=-3.6e3*4.184, HX=2.4e3*4.184, WOCT=2.0e3*4.184, WTET=2.0e3*4.184,
                 H33=-20.0e3*4.184, W33=7.0e3*4.184, W13=10.0e3*4.184, W13P=11.3e3*4.184,
                 W1P4=12.4e3*4.184;
    const double R = 8.3143;
    const double eps = 2.220446049250313e-16;

    const double gs1 = 0.5*(-WOCT + W24U - W14 + 0.5*HX + 0.5*HEX + H24);
    const double gs2 = W11 + H11;
    const double gs3 = W13P - W13 + H33;
    const double hs4 = W15P - W15 + H55;
    const double gx2x2 = -0.25*(WTET+WOCT+HX);
    const double gx2s1 = 0.5*(WOCT-WTET);
    const double gx2s2 = 0.5*(WTET-W22+W11-WOCT+2.0*W2P4U-2.0*W1P4-HEX+2.0*H24);
    const double gx3x3 = -W13;
    const double gx3s3 = W33-W13P+W13;
    const double gx4x4 = -W14;
    const double gx5x5 = -W15;
    const double gx5s4 = W55-W15P+W15;
    const double gs1s1 = 0.25*(-WTET-WOCT+HX);
    const double gs1s2 = 0.5*(WTET-W22+W11+WOCT-HX);
    const double gs2s2 = -W11;
    const double gs3s3 = -W33;
    const double gs4s4 = -W55;

    double s0=0.5, s1=0.9, s2=0.1;
    for (int iter = 0; iter < 1000; iter++){
        double dgds0 = gs1 + gs2/2.0 + gx2s1 + gx2s2/2.0 + 2.0*gs1s1*s0
                     + gs1s2*(0.5+s0) + gs2s2*(1.0+s0)/2.0
                     + R*T*(0.5*creal(clog(1.0+s0)) + 0.5*creal(clog(3.0+s0)) - creal(clog(1.0-s0)));
        double dgds1 = gs2 + 2.0*gs2s2*s1
                     + R*T*(creal(clog(s1)) + creal(clog(1.0+s1)) - 2.0*creal(clog(1.0-s1)));
        double dgds2 = hs4 - T*S55 + gx5s4 + 2.0*gs4s4*s2
                     + R*T*(creal(clog(s2)) - 2.0*creal(clog(1.0-s2)) + creal(clog(1.0+s2)));
        double d00 = 2.0*gs1s1 + gs1s2 + gs2s2/2.0 + R*T*(0.5/(1.0+s0) + 0.5/(3.0+s0) + 1.0/(1.0-s0));
        double d11 = 2.0*gs2s2 + R*T*(1.0/s1 + 1.0/(1.0+s1) + 2.0/(1.0-s1));
        double d22 = 2.0*gs4s4 + R*T*(1.0/s2 + 2.0/(1.0-s2) + 1.0/(1.0+s2));
        double s0N = s0 - dgds0/d00, s1N = s1 - dgds1/d11, s2N = s2 - dgds2/d22;
        s0N = (s0N < 1.0-eps) ? s0N : 1.0-eps;   s1N = (s1N < 1.0-eps) ? s1N : 1.0-eps;   s2N = (s2N < 1.0-eps) ? s2N : 1.0-eps;
        s0N = (s0N > -1.0+eps) ? s0N : -1.0+eps; s1N = (s1N > eps) ? s1N : eps;           s2N = (s2N > eps) ? s2N : eps;
        int converged = (fabs(s0N-s0)<=10.0*eps) && (fabs(s1N-s1)<=10.0*eps) && (fabs(s2N-s2)<=10.0*eps);
        s0=s0N; s1=s1N; s2=s2N;
        if (converged) break;
    }

    double SP_S = -R*(0.5*(1.0+s0)*creal(clog(1.0+s0))+(1.0-s0)*creal(clog(1.0-s0))+0.5*(3.0+s0)*creal(clog(3.0+s0))-5.0*creal(clog(2.0)));
    double SP_H = gs1*s0 + gs2*(1.0+s0)/2.0 + gx2x2 + gx2s1*s0 + gx2s2*(1.0+s0)/2.0
                + gs1s1*s0*s0 + gs1s2*s0*(1.0+s0)/2.0 + gs2s2*(1.0+s0)*(1.0+s0)/4.0;
    *SP_G = SP_H - T*SP_S;

    double HC_S = -R*(s1*creal(clog(s1))+2.0*(1.0-s1)*creal(clog(1.0-s1))+(1.0+s1)*creal(clog(1.0+s1))-2.0*creal(clog(2.0)));
    double HC_H = gs2*s1 + gs2s2*s1*s1;
    *HC_G = HC_H - T*HC_S;

    *CR_G = gs3 + gx3x3 + gx3s3 + gs3s3;   /* CR_S=0.0 identically */

    *UV_G = gx4x4 - T*(R*2.0*creal(clog(2.0)));

    double MT_S = -R*(s2*creal(clog(s2))+2.0*(1.0-s2)*creal(clog(1.0-s2))+(1.0+s2)*creal(clog(1.0+s2))-2.0*creal(clog(2.0))) + S55*s2;
    double MT_H = hs4*s2 + gx5x5 + gx5s4*s2 + gs4s4*s2*s2;
    *MT_G = MT_H - T*MT_S;
}

double obj_gh_spn(unsigned n, const double *x, double *grad, void *SS_ref_db){
    SS_ref *d = (SS_ref *) SS_ref_db;

    double T    = d->T;
    double Rgas = d->R*1000.0;
    double *p   = d->p;
    double *gb  = d->gb_lvl;
    double *mu_Gex = d->mu_Gex;

    for (int i = 0; i < 5; i++){ p[i] = x[i]; }
    double p0 = p[0], p2 = p[2], p3 = p[3], p4 = p[4];   /* p[1]=hercynite, dependent */

    const double H11=-8.7e3*4.184, W11=4.5e3*4.184, W14=20.8e3*4.184, W1P4=12.4e3*4.184,
                 W15=10.0e3*4.184, W1P5=14.4e3*4.184, W15P=11.7e3*4.184, W1P5P=7.0e3*4.184,
                 W22=3.6e3*4.184, H24=6.55e3*4.184, W24U=12.6e3*4.184, W2P4U=10.9e3*4.184,
                 H25=8.05e3*4.184, W25PU=15.3e3*4.184, W2P5U=14.4e3*4.184, W45=6.0e3*4.184,
                 W45P=2.0e3*4.184, W4U5PU=1.8e3*4.184, H55=6.25e3*4.184, S55=0.0, W55=0.0,
                 HEX=-3.6e3*4.184, HX=2.4e3*4.184, WOCT=2.0e3*4.184, WTET=2.0e3*4.184,
                 H33=-20.0e3*4.184, H23=0.0, W33=7.0e3*4.184, W13=10.0e3*4.184,
                 W1P3=15.2e3*4.184, W13P=11.3e3*4.184, W1P3P=5.9e3*4.184, W2P3U=10.0e3*4.184,
                 W23PU=9.7e3*4.184, W34=12.5e3*4.184, W3P4=10.0e3*4.184, W3PU4U=10.4e3*4.184,
                 W35=0.0, W3P5=10.0e3*4.184, W35P=8.0e3*4.184, W3P5P=0.0;
    const double WV1=-0.1250, WV2=0.1018;

    const double gx2x2 = -0.25*(WTET+WOCT+HX);
    const double gx2x3 = 0.5*(W2P3U-W22-W1P3+W11+2.0*H23);
    const double gx2x4 = 0.5*(WTET-W24U+W14+0.5*HX-0.5*HEX+H24);
    const double gx2x5 = 0.5*(-W22+W11+W2P5U-W1P5+2.0*H25);
    const double gx3x3 = -W13;
    const double gx3x4 = W34-W14-W13;
    const double gx3x5 = W35-W15-W13;
    const double gx4x4 = -W14;
    const double gx4x5 = W45-W15-W14;
    const double gx5x5 = -W15;

    double gc[31];
    gc[0]  = 0.5*(-WOCT+W24U-W14+0.5*HX+0.5*HEX+H24);              /* gs1 */
    gc[1]  = W11+H11;                                              /* gs2 */
    gc[2]  = W13P-W13+H33;                                         /* gs3 */
    gc[3]  = W15P-W15+H55;                                         /* hs4 */
    gc[4]  = S55;                                                  /* ss4 */
    gc[5]  = 0.5*(WOCT-WTET);                                      /* gx2s1 */
    gc[6]  = 0.5*(WTET-W22+W11-WOCT+2.0*W2P4U-2.0*W1P4-HEX+2.0*H24); /* gx2s2 */
    gc[7]  = 0.5*(WTET-WOCT+2.0*W3PU4U-2.0*W3P4-W2P3U-W23PU+W22+W1P3+W13P-W11-HEX+2.0*H24-2.0*H23); /* gx2s3 */
    gc[8]  = 0.5*(WTET-WOCT-W11+W22+2.0*W4U5PU-2.0*W45P-W2P5U-W25PU+W1P5+W15P-HEX-2.0*H25+2.0*H24);  /* gx2s4 */
    gc[9]  = 0.5*(W2P3U-W22-W1P3+W11);                             /* gx3s1 */
    gc[10] = W1P3-W13-W11;                                         /* gx3s2 */
    gc[11] = W33-W13P+W13;                                         /* gx3s3 */
    gc[12] = W35P-W35-W15P+W15;                                    /* gx3s4 */
    gc[13] = 0.5*(WTET-W24U+W14-0.5*HX+0.5*HEX-H24);               /* gx4s1 */
    gc[14] = -W11+W1P4-W14;                                        /* gx4s2 */
    gc[15] = W3P4-W34-W13P+W13;                                    /* gx4s3 */
    gc[16] = W45P-W45-W15P+W15;                                    /* gx4s4 */
    gc[17] = 0.5*(-W22+W11+W2P5U-W1P5);                            /* gx5s1 */
    gc[18] = -W11+W1P5-W15;                                        /* gx5s2 */
    gc[19] = W3P5-W35-W13P+W13;                                    /* gx5s3 */
    gc[20] = W55-W15P+W15;                                         /* gx5s4 */
    gc[21] = 0.25*(-WTET-WOCT+HX);                                 /* gs1s1 */
    gc[22] = 0.5*(WTET-W22+W11+WOCT-HX);                           /* gs1s2 */
    gc[23] = 0.5*(WTET+WOCT-W2P3U-W23PU+W22+W1P3+W13P-W11-HX);     /* gs1s3 */
    gc[24] = 0.5*(WTET+WOCT+W22-W11-W2P5U-W25PU+W1P5+W15P-HX);     /* gs1s4 */
    gc[25] = -W11;                                                 /* gs2s2 */
    gc[26] = W1P3P-W1P3-W13P+W13;                                  /* gs2s3 */
    gc[27] = W1P5P-W1P5-W15P+W15;                                  /* gs2s4 */
    gc[28] = -W33;                                                 /* gs3s3 */
    gc[29] = W3P5P-W3P5-W35P+W35;                                  /* gs3s4 */
    gc[30] = -W55;                                                 /* gs4s4 */

    /* --- embedded 3D solve for order parameters s0(M1+M3-like tet/oct split
       for Mg-Fe2+), s1 (Al distribution), s2 (Fe3+ distribution) ---
       Default (gv.gh_multistart_order=0): ONE physically-motivated starting
       guess - the exact same "distribute cations proportional to available
       site ratio" formula xMELTS' own order() function uses - fed into the
       real joint damped 3D Newton solve (GH_spn_solve_s_from), matching
       xMELTS' own order() exactly (see this function's broader header
       comment for why an earlier Gauss-Seidel/bisection substitute was
       replaced).

       Known limitation (confirmed via finite differences, not a
       transcription error - every one of the 31 Taylor coefficients and
       ~35 raw W-parameters was cross-checked term-by-term against xMELTS'
       own spinel.c): this 3D order-parameter problem has genuine multiple
       local minima at extreme, near-single-endmember compositions (e.g.
       ~99.7% pure MgAl2O4 "spl") - a real feature of the non-convex energy
       landscape, not unique to this port. Real MELTS' own single Newton
       path can land on any of these depending on where it starts, same as
       here; this is confined to a composition regime (>99% single
       endmember) far from where any real spinel-bearing bulk would land.

       gv.gh_multistart_order=1 (legacy/opt-in, NOT the default): try
       several starting points spanning the feasible s-cube and keep
       whichever converges to the lowest G - finds a lower, more globally-
       optimal G at those same extreme corners (with a more physically
       "normal"-structure order parameter), but can then disagree with what
       real MELTS itself would compute there, since MELTS never multi-
       starts either. See scratchpad check_spn_multistart2.c. */
    double s0, s1, s2;
    {
        double totAl=2.0*(1.0-p0-p4-p2), totCr=2.0*p0, totFe3=2.0*p2, totMg=p3, totTi=p4;
        double ratio = 2.0-totCr-totTi;
        double xmg2oct0=totMg*ratio/(1.0+ratio);
        double xal3oct0=totAl*ratio/(1.0+ratio), xfe3oct0=totFe3*ratio/(1.0+ratio);
        double xmg2tet0=totMg-xmg2oct0, xal3tet0=totAl-xal3oct0, xfe3tet0=totFe3-xfe3oct0;
        xmg2oct0/=2.0; xal3oct0/=2.0; xfe3oct0/=2.0;
        s0 = xmg2tet0-2.0*xmg2oct0;
        s1 = xal3oct0-xal3tet0/2.0;
        s2 = xfe3oct0-xfe3tet0/2.0;
    }
    if (!GH_spn_multistart_flag){
        double cs0, cs1, cs2;
        GH_spn_solve_s_from(p0,p2,p3,p4,T,Rgas,gc,d->d_em[2],s0,s1,s2,&cs0,&cs1,&cs2);
        s0 = cs0; s1 = cs1; s2 = cs2;
    }
    else {
        const double starts[7][3] = {
            {0.0,0.0,0.0}, {0.8,0.0,0.0}, {-0.8,0.0,0.0}, {0.0,0.8,0.0}, {0.0,-0.8,0.0},
            {0.5,0.5,0.0}, {-0.5,-0.5,0.0}
        };
        double bestG = 1.0e300, bs0=s0, bs1=s1, bs2=s2;
        double cand_starts[8][3];
        cand_starts[0][0]=s0; cand_starts[0][1]=s1; cand_starts[0][2]=s2;
        for (int i = 0; i < 7; i++){ cand_starts[i+1][0]=starts[i][0]; cand_starts[i+1][1]=starts[i][1]; cand_starts[i+1][2]=starts[i][2]; }
        for (int i = 0; i < 8; i++){
            double cs0, cs1, cs2;
            GH_spn_solve_s_from(p0,p2,p3,p4,T,Rgas,gc,d->d_em[2],cand_starts[i][0],cand_starts[i][1],cand_starts[i][2],&cs0,&cs1,&cs2);
            double G = GH_spn_G_at(p0,p2,p3,p4,cs0,cs1,cs2,T,Rgas,gc,
                                    gx2x2,gx2x3,gx2x4,gx2x5,gx3x3,gx3x4,gx3x5,gx4x4,gx4x5,gx5x5,
                                    WV1,WV2,d->P,d->d_em[0],d->d_em[2],d->d_em[4]);
            if (G < bestG){ bestG = G; bs0 = cs0; bs1 = cs1; bs2 = cs2; }
        }
        s0 = bs0; s1 = bs1; s2 = bs2;
    }

    /* xcr3oct/xti4oct/xfe3tet/xfe3oct get the same +d_em guard used
       throughout gh (hb, melilite): if chromite, ulvospinel or magnetite is
       boiled out (bounds_ref pins its p to exactly 0), log() of its site
       fraction would otherwise be log(0)=NaN. The +d_em2 guard here keeps
       the S/H/dGdp block finite regardless of whatever s2 the Newton solve
       converges to in that degenerate case, consistent with z_em[2]=0
       already excluding magnetite from the real mass balance. */
    double xmg2tet=(p3+s0)/2.0, xfe2tet=p4-0.5*p3-0.5*s0+s1+p0+s2, xal3tet=1.0-p0-p4-p2-s1, xfe3tet=p2-s2+d->d_em[2];
    double xmg2oct=(p3-s0)/4.0, xfe2oct=(2.0-p3+s0-2.0*s1-2.0*p0-2.0*s2)/4.0, xal3oct=(1.0-p0-p4-p2+s1)/2.0, xfe3oct=(p2+s2)/2.0+d->d_em[2];
    double xcr3oct=p0+d->d_em[0], xti4oct=(p4+d->d_em[4])/2.0;

    double H = gx2x2*p3*p3 + gx2x3*p3*p0 + gx2x4*p3*p4 + gx2x5*p3*p2
             + gc[5]*p3*s0 + gc[6]*p3*s1 + gc[7]*p3*p0 + gc[8]*p3*s2
             + gx3x3*p0*p0 + gx3x4*p0*p4 + gx3x5*p0*p2 + gc[9]*p0*s0 + gc[10]*p0*s1
             + gc[11]*p0*p0 + gc[12]*p0*s2
             + gx4x4*p4*p4 + gx4x5*p4*p2 + gc[13]*p4*s0 + gc[14]*p4*s1 + gc[15]*p4*p0 + gc[16]*p4*s2
             + gx5x5*p2*p2 + gc[17]*p2*s0 + gc[18]*p2*s1 + gc[19]*p2*p0 + gc[20]*p2*s2
             + gc[0]*s0 + gc[1]*s1 + gc[2]*p0 + gc[3]*s2
             + gc[21]*s0*s0 + gc[22]*s0*s1 + gc[23]*s0*p0 + gc[24]*s0*s2
             + gc[25]*s1*s1 + gc[26]*s1*p0 + gc[27]*s1*s2
             + gc[28]*p0*p0 + gc[29]*p0*s2 + gc[30]*s2*s2;

    double S = -Rgas*( xmg2tet*creal(clog(xmg2tet)) + xfe2tet*creal(clog(xfe2tet))
                      + xal3tet*creal(clog(xal3tet)) + xfe3tet*creal(clog(xfe3tet))
                      + 2.0*xmg2oct*creal(clog(xmg2oct)) + 2.0*xfe2oct*creal(clog(xfe2oct))
                      + 2.0*xal3oct*creal(clog(xal3oct)) + 2.0*xfe3oct*creal(clog(xfe3oct))
                      + 2.0*xcr3oct*creal(clog(xcr3oct)) + 2.0*xti4oct*creal(clog(xti4oct)) ) + gc[4]*s2;

    double Gvol = p4*p2*(WV1*p2+WV2*p4)*((d->P)-1.0);
    double Graw = H - T*S + Gvol;

    /* ENDMEMBERS(r) = sum(p_i * that endmember's own individual ordering
       correction) - see GH_spn_mix_endmember_G's header comment. p[1]=herc
       is genuinely present here (unlike in Graw, which only depends on
       p[1] implicitly as the dependent endmember), so its own dGex
       contribution is nonzero below, same pattern as rhm's hem fix. */
    double CR_G, HC_G, MT_G, SP_G, UV_G;
    GH_spn_mix_endmember_G(T, &CR_G, &HC_G, &MT_G, &SP_G, &UV_G);
    double ENDMEMBERS = p0*CR_G + p[1]*HC_G + p2*MT_G + p3*SP_G + p4*UV_G;
    double Gex = Graw - ENDMEMBERS;

    /* dGdp0: chromite (p0=r1) plays a dual role in xMELTS' own model (see
       this function's header comment) - its "gx3-family" terms (p0 as an
       ordinary composition variable, e.g. gx3x3*p0^2) and its "gs3-family"
       terms (p0 standing in for what would otherwise be order parameter
       s3, e.g. gs3s3*p0^2 via gc[28] - a DIFFERENT W-parameter, W33 not
       W13) are both genuinely present as distinct additive terms in H -
       but each term differentiates ONCE, not once per "role". An earlier
       version of this formula double-counted gc[2] and gc[7]*p3 while
       omitting gx2x3*p3 entirely (the term symmetric to dGdp3's own
       gx2x3*p0, from the same gx2x3*p3*p0 cross term in H) - a real bug,
       not a deliberate double-block sum: confirmed by finite-difference
       (grad[0] was off by ~19 kJ at every tested composition, the other
       four gradients were already exact) and by symmetry with dGdp2/
       dGdp3/dGdp4, which already collect each H cross-term exactly once.
       This is very likely why spinel's NLopt optimizations needed ~20x
       more evaluations than cpx's to (mostly not) converge - SLSQP was
       being fed a wrong gradient, not fighting genuine landscape
       non-smoothness. */
    double dGdp0 = gc[2] + gc[7]*p3 + gx2x3*p3 + 2.0*gx3x3*p0 + gx3x4*p4 + gx3x5*p2 + gc[9]*s0 + gc[10]*s1 + 2.0*gc[11]*p0 + gc[12]*s2
                 + gc[15]*p4 + gc[19]*p2 + gc[23]*s0 + gc[26]*s1 + 2.0*gc[28]*p0 + gc[29]*s2
                 + Rgas*T*(creal(clog(xfe2tet/xal3tet)) + 2.0*creal(clog(xcr3oct)) - creal(clog(xfe2oct)) - creal(clog(xal3oct)));

    /* dGdp2/dGdp4 pressure (Gvol) cross terms: xMELTS' own DGDR3/DGDR2
       macros (spinel.c ~1595-1600) differentiate Gvol =
       p2*p4*(WV1*p2+WV2*p4)*(P-1) asymmetrically - the "extra" cross term
       beyond the naive product-rule piece carries WV1 for dGdp2 (d/d[mt])
       and WV2 for dGdp4 (d/d[usp]), confirmed independently by direct
       calculus on this file's own Gvol formula in GH_spn_G_at. An earlier
       version had these swapped (WV2 in dGdp2, WV1 in dGdp4) - a real,
       wrong-sign gradient bug (WV1=-0.1250 vs WV2=+0.1018, opposite sign)
       wherever mt and usp are both non-negligible, exactly the composition
       regime a MgCr2O4-poor, Fe-Ti-rich spinel lands in. Found 2026-07-13
       while investigating gh predicting spinel stable ~280C above real
       xMELTS' liquidus for the same bulk/T/P. */
    double dGdp2 = gx5x5*2.0*p2 + gx2x5*p3 + gx3x5*p0 + gx4x5*p4 + gc[17]*s0 + gc[18]*s1 + gc[19]*p0 + gc[20]*s2
                 + Rgas*T*(creal(clog(xfe3tet/xal3tet)) + creal(clog(xfe3oct/xal3oct)))
                 + p4*(WV1*p2+WV2*p4)*((d->P)-1.0) + p2*p4*WV1*((d->P)-1.0);
    double dGdp3 = gx2x2*2.0*p3 + gx2x3*p0 + gx2x4*p4 + gx2x5*p2 + gc[5]*s0 + gc[6]*s1 + gc[7]*p0 + gc[8]*s2
                 + 0.5*Rgas*T*(creal(clog(xmg2tet/xfe2tet)) + creal(clog(xmg2oct/xfe2oct)));
    double dGdp4 = gx4x4*2.0*p4 + gx2x4*p3 + gx3x4*p0 + gx4x5*p2 + gc[13]*s0 + gc[14]*s1 + gc[15]*p0 + gc[16]*s2
                 + Rgas*T*(creal(clog(xfe2tet/xal3tet)) + creal(clog(xti4oct/xal3oct)))
                 + p2*(WV1*p2+WV2*p4)*((d->P)-1.0) + p2*p4*WV2*((d->P)-1.0);

    double dGex[5];
    dGex[0] = dGdp0 - CR_G;
    dGex[1] = -HC_G;   /* hercynite is dependent in the r-basis but present in ENDMEMBERS */
    dGex[2] = dGdp2 - MT_G;
    dGex[3] = dGdp3 - SP_G;
    dGex[4] = dGdp4 - UV_G;

    for (int i = 0; i < 5; i++){
        double mu = Gex;
        for (int k = 0; k < 5; k++){
            double delta_ik = (i==k) ? 1.0 : 0.0;
            mu += (delta_ik - p[k])*dGex[k];
        }
        mu_Gex[i] = mu/1000.0;
        d->sf[i]  = p[i];
    }

    d->sum_apep = 0.0;
    for (int i = 0; i < 5; i++){ d->sum_apep += d->ape[i]*p[i]; }
    d->factor = d->fbc/d->sum_apep;

    d->df_raw = 0.0;
    for (int i = 0; i < 5; i++){ d->df_raw += (mu_Gex[i] + gb[i])*p[i]; }
    d->df = d->df_raw * d->factor;

    if (grad){
        for (int i = 0; i < 5; i++){
            grad[i] = (mu_Gex[i] + gb[i])*d->factor - (d->df_raw*d->factor*(d->ape[i]/d->sum_apep));
        }
    }
    return d->df;
}

/**
    Rhombohedral oxide (Geikielite-Hematite-Ilmenite-Pyrophanite-Corundum,
    Ghiorso & Evans 2008): from xMELTS' sources/rhomsghiorso.c. NR=4
    (r0=Xilm, r1=Xgei, r2=Xpyr, r3=Xcrn, Xhem=1-r0-r1-r2-r3 dependent),
    NS=3 joint order parameters (s0=ilmenite's, s1=geikielite's,
    s2=pyrophanite's own A/B-site cation-ordering degree - a genuinely
    different, JOINT solve from the 3 independent per-pure-endmember
    order-parameter solves already folded into gei/ilm/pyr's own gbase,
    see GH_rhm_pure_order_G in GH_gem_function.c). p[]/EM_list order
    (gei,hem,ilm,pyr,crn) is a pure relabeling of xMELTS' own r-basis:
    p[0]=r1, p[1]=Xhem (dependent), p[2]=r0, p[3]=r2, p[4]=r3 - same
    "direct r-basis" pattern spn's chr/herc/mt/spl/usp already uses, no
    basis rotation needed.

    W-parameter table below is ported verbatim (names and numeric values)
    from rhomsghiorso.c's file-scope statics. The joint order-parameter
    solve (GH_rhm_solve_s_from) ports xMELTS' own order() function
    EXACTLY, per this project's solver-fidelity principle (see spn's own
    header comment above for why a plausible-looking substitute silently
    broke NLopt convergence there): a full (UNDAMPED) joint 3-variable
    Newton step each iteration (3x3 Hessian solved via Cramer's rule - no
    shared linear-algebra utility exists in this codebase, same "same
    math, no gaussj" approach spn's own joint solve uses), a simple
    box-clamp of s[i] into [-r[i]+DBL_EPSILON, r[i]-DBL_EPSILON] after
    each step (MIN then MAX, in that order - matching the source exactly,
    including its degenerate-but-well-defined behavior of pinning s to
    +DBL_EPSILON when r[i]=0, e.g. a boiled-out endmember), starting
    guess s[i]=0.9*r[i], convergence tolerance 10*DBL_EPSILON per
    component, MAX_ITER=100 - deliberately NOT spn's own per-site-fraction
    damped-lambda line search, since that scheme is specific to
    spinel.c's order() and has no counterpart in rhomsghiorso.c's.

    The short-range-order spline fSRO(T,order) is a natural cubic spline
    through 8 calibrated points that are all the identical constant
    0.0730205 in xMELTS' own calibration - exactly flat, so
    fSRO(T,0)==0.0730205 and every T-derivative is 0 for any T (see
    GH_rhm_SRO_G in GH_gem_function.c for the same simplification applied
    to hematite/corundum's own standard states).

    Total molar G = sum(p[i]*gbase[i]) + Gex(r,s*(r)), where gbase[i]
    already carries each endmember's own ordering/SRO correction (matching
    xMELTS' own G_total = ENDMEMBERS(r) + Gmix(r) construction, see
    gmixMsg's *gmix = G - ENDMEMBERS in the source). Gex here MUST
    therefore be Gmix(r) = G(r,s*) - ENDMEMBERS(r), not the raw joint
    G(r,s*) - since gbase[] already carries ENDMEMBERS(r)'s per-endmember
    pieces once, adding the un-subtracted G(r,s*) on top double-counts
    them. (An earlier version of this function skipped the subtraction,
    reasoning gbase[] already had it covered - true for the total's
    ENDMEMBERS(r) term, but wrong for Gex, which must independently
    vanish at every pure endmember composition like any proper mixing
    energy; this is a real bug, found via [[gh-gexcess-verification]] and
    fixed 2026-07-14. dGex now includes -ENDMEMBERS(r)'s own p-gradient,
    see GH_rhm_mix_endmember_order_G/GH_rhm_mix_endmember_SRO_G above.)
    dGdr from GH_rhm_dGdr uses the envelope theorem (DGDR0..3 evaluated at
    the converged s*, no ds/dr needed, since dGex/ds=0 there by
    construction of the joint solve) plus the same delta_ik-p_k tangent-
    plane trick used throughout gh (liq/ol/spn) to turn a single scalar
    Gex + its analytic p-gradient into per-endmember mu_Gex[] with
    sum(p[i]*mu_Gex[i])==Gex exactly.
*/
static const double GH_rhm_dvilm=0.010758, GH_rhm_dvgei=0.010758, GH_rhm_dvpyr=0.010758;
static const double GH_rhm_wvilm=0.035089, GH_rhm_wvgei=0.035089, GH_rhm_wvpyr=0.035089;
static const double GH_rhm_dwvhmilm=0.013701, GH_rhm_dwvhmgei=0.013701, GH_rhm_dwvhmpyr=0.013701;
static const double GH_rhm_wvhmilm2=0.0, GH_rhm_wvhmgei2=0.0, GH_rhm_wvhmpyr2=0.0;
static const double GH_rhm_dwvcrnilm=0.013701, GH_rhm_dwvcrngei=0.013701, GH_rhm_dwvcrnpyr=0.013701;
static const double GH_rhm_wvhmilm=-0.11764, GH_rhm_wvhmgei=-0.11764, GH_rhm_wvhmpyr=-0.11764;
static const double GH_rhm_wvilmgei=0.0, GH_rhm_wvilmpyr=0.0, GH_rhm_wvgeipyr=0.0;
static const double GH_rhm_dhilm=17477.0, GH_rhm_dhgei=17477.0, GH_rhm_dhpyr=17477.0;
static const double GH_rhm_whilm=3189.0, GH_rhm_whgei=3189.0, GH_rhm_whpyr=3189.0;
static const double GH_rhm_dwhhmilm=-5626.63, GH_rhm_dwhhmgei=-5626.63, GH_rhm_dwhhmpyr=-5626.63;
static const double GH_rhm_whhmilm2=-833.14, GH_rhm_whhmgei2=-833.14, GH_rhm_whhmpyr2=-833.14;
static const double GH_rhm_dwhcrnilm=0.0, GH_rhm_dwhcrngei=0.0, GH_rhm_dwhcrnpyr=0.0;
static const double GH_rhm_whmcrn=69000.0;
static const double GH_rhm_whhmilm=22535.6, GH_rhm_whhmgei=22535.6, GH_rhm_whhmpyr=22535.6;
static const double GH_rhm_wcrnilm=22535.6, GH_rhm_wcrngei=22535.6, GH_rhm_wcrnpyr=22535.6;
static const double GH_rhm_whilmgei=2600.0, GH_rhm_whilmpyr=2200.0, GH_rhm_whgeipyr=2600.0;
static const double GH_rhm_whilmgeiT=88099.9, GH_rhm_whilmpyrT=30244.0, GH_rhm_whgeipyrT=2600.0;
/* whilmilmgei=whilmgeigei=whilmilmpyr=whilmpyrpyr=whgeigeipyr=whgeipyrpyr=0.0
   in xMELTS' own calibration - the H/DGDR/D2GDS2 terms they multiply are
   simply omitted below rather than carried as always-0 additions. */
static const double GH_rhm_SROconst = 0.0730205;
static const double GH_rhm_R = 8.3143;

static void GH_rhm_site_fracs(const double r[4], const double s[3],
    double *xfe2a,double *xmg2a,double *xmn2a,double *xti4a,double *xal3a,double *xfe3a,
    double *xfe2b,double *xmg2b,double *xmn2b,double *xti4b,double *xal3b,double *xfe3b){
    double eps = 2.220446049250313e-16;
    double v;
    v=(r[0]+s[0])/2.0;                              *xfe2a=(v<=0.0)?eps:v;
    v=(r[1]+s[1])/2.0;                              *xmg2a=(v<=0.0)?eps:v;
    v=(r[2]+s[2])/2.0;                              *xmn2a=(v<=0.0)?eps:v;
    v=(r[0]-s[0]+r[1]-s[1]+r[2]-s[2])/2.0;          *xti4a=(v<=0.0)?eps:v;
    v=r[3];                                          *xal3a=(v<=0.0)?eps:v;
    v=1.0-r[0]-r[1]-r[2]-r[3];                      *xfe3a=(v<=0.0)?eps:v;

    v=(r[0]-s[0])/2.0;                              *xfe2b=(v<=0.0)?eps:v;
    v=(r[1]-s[1])/2.0;                              *xmg2b=(v<=0.0)?eps:v;
    v=(r[2]-s[2])/2.0;                              *xmn2b=(v<=0.0)?eps:v;
    v=(r[0]+s[0]+r[1]+s[1]+r[2]+s[2])/2.0;          *xti4b=(v<=0.0)?eps:v;
    v=r[3];                                          *xal3b=(v<=0.0)?eps:v;
    v=1.0-r[0]-r[1]-r[2]-r[3];                      *xfe3b=(v<=0.0)?eps:v;
}
static void GH_rhm_ID_fracs(const double r[4],
    double *xfe2ID,double *xmg2ID,double *xmn2ID,double *xti4ID,double *xal3ID,double *xfe3ID){
    double eps = 2.220446049250313e-16;
    *xfe2ID = (r[0]/2.0 > 0.0) ? r[0]/2.0 : eps;
    *xmg2ID = (r[1]/2.0 > 0.0) ? r[1]/2.0 : eps;
    *xmn2ID = (r[2]/2.0 > 0.0) ? r[2]/2.0 : eps;
    *xti4ID = ((r[0]+r[1]+r[2])/2.0 > 0.0) ? (r[0]+r[1]+r[2])/2.0 : eps;
    *xal3ID = (r[3] > 0.0) ? r[3] : eps;
    *xfe3ID = (1.0-r[0]-r[1]-r[2]-r[3] > 0.0) ? 1.0-r[0]-r[1]-r[2]-r[3] : eps;
}
static void GH_rhm_dgds(const double r[4], const double s[3], double t, double p,
                         double *dgds0, double *dgds1, double *dgds2){
    double xfe2a,xmg2a,xmn2a,xti4a,xal3a,xfe3a,xfe2b,xmg2b,xmn2b,xti4b,xal3b,xfe3b;
    GH_rhm_site_fracs(r,s,&xfe2a,&xmg2a,&xmn2a,&xti4a,&xal3a,&xfe3a,&xfe2b,&xmg2b,&xmn2b,&xti4b,&xal3b,&xfe3b);
    double R = GH_rhm_R;
    double r0=r[0],r1=r[1],r2=r[2],r3=r[3], s0=s[0],s1=s[1],s2=s[2];
    double sumr = r0+r1+r2+r3;

    *dgds0 = 0.5*R*t*( creal(clog(xfe2a)) - creal(clog(xti4a)) - creal(clog(xfe2b)) + creal(clog(xti4b)) )
           + (GH_rhm_whilmgei-GH_rhm_whilmgeiT)*s1/2.0 + (GH_rhm_whilmpyr-GH_rhm_whilmpyrT)*s2/2.0
           - 2.0*(GH_rhm_dhilm+(p-1.0)*GH_rhm_dvilm+((GH_rhm_dwhhmilm+(p-1.0)*GH_rhm_dwvhmilm)/2.0+(GH_rhm_whhmilm2+(p-1.0)*GH_rhm_wvhmilm2)/4.0)*(1.0-sumr)
                    -(GH_rhm_dwhcrnilm+(p-1.0)*GH_rhm_dwvcrnilm)*r3/2.0)*s0
           + 2.0*(GH_rhm_whilm+(p-1.0)*GH_rhm_wvilm)*s0*(r0*r0-2.0*s0*s0);

    *dgds1 = 0.5*R*t*( creal(clog(xmg2a)) - creal(clog(xti4a)) - creal(clog(xmg2b)) + creal(clog(xti4b)) )
           + (GH_rhm_whilmgei-GH_rhm_whilmgeiT)*s0/2.0 + (GH_rhm_whgeipyr-GH_rhm_whgeipyrT)*s2/2.0
           - 2.0*(GH_rhm_dhgei+(p-1.0)*GH_rhm_dvgei+((GH_rhm_dwhhmgei+(p-1.0)*GH_rhm_dwvhmgei)/2.0+(GH_rhm_whhmgei2+(p-1.0)*GH_rhm_wvhmgei2)/4.0)*(1.0-sumr)
                    -(GH_rhm_dwhcrngei+(p-1.0)*GH_rhm_dwvcrngei)*r3/2.0)*s1
           + 2.0*(GH_rhm_whgei+(p-1.0)*GH_rhm_wvgei)*s1*(r1*r1-2.0*s1*s1);

    *dgds2 = 0.5*R*t*( creal(clog(xmn2a)) - creal(clog(xti4a)) - creal(clog(xmn2b)) + creal(clog(xti4b)) )
           + (GH_rhm_whilmpyr-GH_rhm_whilmpyrT)*s0/2.0 + (GH_rhm_whgeipyr-GH_rhm_whgeipyrT)*s1/2.0
           - 2.0*(GH_rhm_dhpyr+(p-1.0)*GH_rhm_dvpyr+((GH_rhm_dwhhmpyr+(p-1.0)*GH_rhm_dwvhmpyr)/2.0+(GH_rhm_whhmpyr2+(p-1.0)*GH_rhm_wvhmpyr2)/4.0)*(1.0-sumr)
                    -(GH_rhm_dwhcrnpyr+(p-1.0)*GH_rhm_dwvcrnpyr)*r3/2.0)*s2
           + 2.0*(GH_rhm_whpyr+(p-1.0)*GH_rhm_wvpyr)*s2*(r2*r2-2.0*s2*s2);
}
static void GH_rhm_d2gds2(const double r[4], const double s[3], double t, double p,
                           double *d00,double *d01,double *d02,double *d11,double *d12,double *d22){
    double xfe2a,xmg2a,xmn2a,xti4a,xal3a,xfe3a,xfe2b,xmg2b,xmn2b,xti4b,xal3b,xfe3b;
    GH_rhm_site_fracs(r,s,&xfe2a,&xmg2a,&xmn2a,&xti4a,&xal3a,&xfe3a,&xfe2b,&xmg2b,&xmn2b,&xti4b,&xal3b,&xfe3b);
    double R = GH_rhm_R;
    double r0=r[0],r1=r[1],r2=r[2],r3=r[3], s0=s[0],s1=s[1],s2=s[2];
    double sumr = r0+r1+r2+r3;

    *d00 = 0.25*R*t*( 1.0/xfe2a + 1.0/xti4a + 1.0/xfe2b + 1.0/xti4b )
         - 2.0*(GH_rhm_dhilm+(p-1.0)*GH_rhm_dvilm+((GH_rhm_dwhhmilm+(p-1.0)*GH_rhm_dwvhmilm)/2.0+(GH_rhm_whhmilm2+(p-1.0)*GH_rhm_wvhmilm2)/4.0)*(1.0-sumr)
                  -(GH_rhm_dwhcrnilm+(p-1.0)*GH_rhm_dwvcrnilm)*r3/2.0)
         + 2.0*(GH_rhm_whilm+(p-1.0)*GH_rhm_wvilm)*(r0*r0-6.0*s0*s0);
    *d01 = 0.25*R*t*( 1.0/xti4a + 1.0/xti4b ) + (GH_rhm_whilmgei-GH_rhm_whilmgeiT)/2.0;
    *d02 = 0.25*R*t*( 1.0/xti4a + 1.0/xti4b ) + (GH_rhm_whilmpyr-GH_rhm_whilmpyrT)/2.0;
    *d11 = 0.25*R*t*( 1.0/xmg2a + 1.0/xti4a + 1.0/xmg2b + 1.0/xti4b )
         - 2.0*(GH_rhm_dhgei+(p-1.0)*GH_rhm_dvgei+((GH_rhm_dwhhmgei+(p-1.0)*GH_rhm_dwvhmgei)/2.0+(GH_rhm_whhmgei2+(p-1.0)*GH_rhm_wvhmgei2)/4.0)*(1.0-sumr)
                  -(GH_rhm_dwhcrngei+(p-1.0)*GH_rhm_dwvcrngei)*r3/2.0)
         + 2.0*(GH_rhm_whgei+(p-1.0)*GH_rhm_wvgei)*(r1*r1-6.0*s1*s1);
    *d12 = 0.25*R*t*( 1.0/xti4a + 1.0/xti4b ) + (GH_rhm_whgeipyr-GH_rhm_whgeipyrT)/2.0;
    *d22 = 0.25*R*t*( 1.0/xmn2a + 1.0/xti4a + 1.0/xmn2b + 1.0/xti4b )
         - 2.0*(GH_rhm_dhpyr+(p-1.0)*GH_rhm_dvpyr+((GH_rhm_dwhhmpyr+(p-1.0)*GH_rhm_dwvhmpyr)/2.0+(GH_rhm_whhmpyr2+(p-1.0)*GH_rhm_wvhmpyr2)/4.0)*(1.0-sumr)
                  -(GH_rhm_dwhcrnpyr+(p-1.0)*GH_rhm_dwvcrnpyr)*r3/2.0)
         + 2.0*(GH_rhm_whpyr+(p-1.0)*GH_rhm_wvpyr)*(r2*r2-6.0*s2*s2);
}
/* Joint order-parameter solve, ported directly from xMELTS' own order()
   function (rhomsghiorso.c) - see this block's header comment above for
   why this must match the source's exact algorithm (undamped full Newton
   step + simple box clamp), not a different-but-plausible substitute. */
static void GH_rhm_solve_s_from(const double r[4], double t, double p,
                                 double *s0o, double *s1o, double *s2o){
    const int maxIter = 100; /* xMELTS' own MAX_ITER (rhomsghiorso.c) */
    double eps = 2.220446049250313e-16;
    double s[3] = { 0.9*r[0], 0.9*r[1], 0.9*r[2] };
    double sOld[3] = { 2.0, 2.0, 2.0 };
    int iter = 0;
    while ( (fabs(s[0]-sOld[0]) > 10.0*eps ||
             fabs(s[1]-sOld[1]) > 10.0*eps ||
             fabs(s[2]-sOld[2]) > 10.0*eps) && iter < maxIter ){
        double dgds0,dgds1,dgds2;
        GH_rhm_dgds(r, s, t, p, &dgds0, &dgds1, &dgds2);
        double d00,d01,d02,d11,d12,d22;
        GH_rhm_d2gds2(r, s, t, p, &d00,&d01,&d02,&d11,&d12,&d22);

        sOld[0]=s[0]; sOld[1]=s[1]; sOld[2]=s[2];

        /* solve [[d00,d01,d02],[d01,d11,d12],[d02,d12,d22]]*delta = -[dgds0,dgds1,dgds2]
           via Cramer's rule (no shared linear-algebra utility exists in
           this codebase - same "same math, no gaussj" approach spn's own
           joint solve uses). */
        double g0=-dgds0, g1=-dgds1, g2=-dgds2;
        double det  = d00*(d11*d22-d12*d12) - d01*(d01*d22-d12*d02) + d02*(d01*d12-d11*d02);
        double det0 =  g0*(d11*d22-d12*d12) - d01*( g1*d22-d12* g2) + d02*( g1*d12-d11* g2);
        double det1 = d00*( g1*d22- g2*d12) -  g0*(d01*d22-d12*d02) + d02*(d01* g2- g1*d02);
        double det2 = d00*(d11* g2-d12* g1) - d01*(d01* g2-d12* g0) +  g0*(d01*d12-d11*d02);
        double delta0 = det0/det, delta1 = det1/det, delta2 = det2/det;

        s[0] = sOld[0] + delta0;
        s[1] = sOld[1] + delta1;
        s[2] = sOld[2] + delta2;

        /* box-clamp into (-r[i],r[i]) - MIN then MAX, in that order,
           matching xMELTS' own order() exactly (including its degenerate-
           but-well-defined behavior of pinning s to +eps when r[i]=0). */
        for (int i = 0; i < 3; i++){
            s[i] = (s[i] < r[i]-eps) ? s[i] : r[i]-eps;
            s[i] = (s[i] > -r[i]+eps) ? s[i] : -r[i]+eps;
        }
        iter++;
    }
    *s0o = s[0]; *s1o = s[1]; *s2o = s[2];
}
static double GH_rhm_S(const double r[4], const double s[3], double t){
    double xfe2a,xmg2a,xmn2a,xti4a,xal3a,xfe3a,xfe2b,xmg2b,xmn2b,xti4b,xal3b,xfe3b;
    GH_rhm_site_fracs(r,s,&xfe2a,&xmg2a,&xmn2a,&xti4a,&xal3a,&xfe3a,&xfe2b,&xmg2b,&xmn2b,&xti4b,&xal3b,&xfe3b);
    double xfe2ID,xmg2ID,xmn2ID,xti4ID,xal3ID,xfe3ID;
    GH_rhm_ID_fracs(r,&xfe2ID,&xmg2ID,&xmn2ID,&xti4ID,&xal3ID,&xfe3ID);
    double R = GH_rhm_R;

    return -R*( xfe2a*creal(clog(xfe2a)) + xmg2a*creal(clog(xmg2a)) + xmn2a*creal(clog(xmn2a)) + xti4a*creal(clog(xti4a))
              + xfe2b*creal(clog(xfe2b)) + xmg2b*creal(clog(xmg2b)) + xmn2b*creal(clog(xmn2b)) + xti4b*creal(clog(xti4b))
              - 2.0*xfe2ID*creal(clog(xfe2ID)) - 2.0*xmg2ID*creal(clog(xmg2ID)) - 2.0*xmn2ID*creal(clog(xmn2ID)) - 2.0*xti4ID*creal(clog(xti4ID)) )
         - (1.0-GH_rhm_SROconst)*2.0*R*( xfe3ID*creal(clog(xfe3ID)) + xal3ID*creal(clog(xal3ID)) + xfe2ID*creal(clog(xfe2ID))
                                        + xmg2ID*creal(clog(xmg2ID)) + xmn2ID*creal(clog(xmn2ID)) + xti4ID*creal(clog(xti4ID)) )
         + GH_rhm_SROconst*2.0*R*creal(clog(2.0));
}
static double GH_rhm_H(const double r[4], const double s[3], double p){
    double r0=r[0],r1=r[1],r2=r[2],r3=r[3], s0=s[0],s1=s[1],s2=s[2];
    double sumr = r0+r1+r2+r3;
    return
          ((GH_rhm_whhmilm+(p-1.0)*GH_rhm_wvhmilm)+(GH_rhm_dwhhmilm+(p-1.0)*GH_rhm_dwvhmilm)*(1.0-2.0*r0-r1-r2-r3))*r0*(1.0-sumr)
        + ((GH_rhm_whhmgei+(p-1.0)*GH_rhm_wvhmgei)+(GH_rhm_dwhhmgei+(p-1.0)*GH_rhm_dwvhmgei)*(1.0-r0-2.0*r1-r2-r3))*r1*(1.0-sumr)
        + ((GH_rhm_whhmpyr+(p-1.0)*GH_rhm_wvhmpyr)+(GH_rhm_dwhhmpyr+(p-1.0)*GH_rhm_dwvhmpyr)*(1.0-r0-r1-2.0*r2-r3))*r2*(1.0-sumr)
        +  GH_rhm_whmcrn*r3*(1.0-sumr)
        + (GH_rhm_whilmgei+(p-1.0)*GH_rhm_wvilmgei+GH_rhm_whilmgeiT)*r0*r1/2.0 + (GH_rhm_whilmgei-GH_rhm_whilmgeiT)*s0*s1/2.0
        + (GH_rhm_whilmpyr+(p-1.0)*GH_rhm_wvilmpyr+GH_rhm_whilmpyrT)*r0*r2/2.0 + (GH_rhm_whilmpyr-GH_rhm_whilmpyrT)*s0*s2/2.0
        + (GH_rhm_whgeipyr+(p-1.0)*GH_rhm_wvgeipyr+GH_rhm_whgeipyrT)*r1*r2/2.0 + (GH_rhm_whgeipyr-GH_rhm_whgeipyrT)*s1*s2/2.0
        + (GH_rhm_wcrnilm+(GH_rhm_dwhcrnilm+(p-1.0)*GH_rhm_dwvcrnilm)*(r0-r3))*r0*r3
        + (GH_rhm_wcrngei+(GH_rhm_dwhcrngei+(p-1.0)*GH_rhm_dwvcrngei)*(r1-r3))*r1*r3
        + (GH_rhm_wcrnpyr+(GH_rhm_dwhcrnpyr+(p-1.0)*GH_rhm_dwvcrnpyr)*(r2-r3))*r2*r3
        + (GH_rhm_dhilm+(p-1.0)*GH_rhm_dvilm+((GH_rhm_dwhhmilm+(p-1.0)*GH_rhm_dwvhmilm)/2.0+(GH_rhm_whhmilm2+(p-1.0)*GH_rhm_wvhmilm2)/4.0)*(1.0-sumr)
                    -(GH_rhm_dwhcrnilm+(p-1.0)*GH_rhm_dwvcrnilm)*r3/2.0)*(r0*r0-s0*s0)
        + (GH_rhm_dhgei+(p-1.0)*GH_rhm_dvgei+((GH_rhm_dwhhmgei+(p-1.0)*GH_rhm_dwvhmgei)/2.0+(GH_rhm_whhmgei2+(p-1.0)*GH_rhm_wvhmgei2)/4.0)*(1.0-sumr)
                    -(GH_rhm_dwhcrngei+(p-1.0)*GH_rhm_dwvcrngei)*r3/2.0)*(r1*r1-s1*s1)
        + (GH_rhm_dhpyr+(p-1.0)*GH_rhm_dvpyr+((GH_rhm_dwhhmpyr+(p-1.0)*GH_rhm_dwvhmpyr)/2.0+(GH_rhm_whhmpyr2+(p-1.0)*GH_rhm_wvhmpyr2)/4.0)*(1.0-sumr)
                    -(GH_rhm_dwhcrnpyr+(p-1.0)*GH_rhm_dwvcrnpyr)*r3/2.0)*(r2*r2-s2*s2)
        + (GH_rhm_whilm+(p-1.0)*GH_rhm_wvilm)*s0*s0*(r0*r0-s0*s0)
        + (GH_rhm_whgei+(p-1.0)*GH_rhm_wvgei)*s1*s1*(r1*r1-s1*s1)
        + (GH_rhm_whpyr+(p-1.0)*GH_rhm_wvpyr)*s2*s2*(r2*r2-s2*s2);
}
static void GH_rhm_dGdr(const double r[4], const double s[3], double t, double p,
                         double *dgdr0,double *dgdr1,double *dgdr2,double *dgdr3){
    double xfe2a,xmg2a,xmn2a,xti4a,xal3a,xfe3a,xfe2b,xmg2b,xmn2b,xti4b,xal3b,xfe3b;
    GH_rhm_site_fracs(r,s,&xfe2a,&xmg2a,&xmn2a,&xti4a,&xal3a,&xfe3a,&xfe2b,&xmg2b,&xmn2b,&xti4b,&xal3b,&xfe3b);
    double xfe2ID,xmg2ID,xmn2ID,xti4ID,xal3ID,xfe3ID;
    GH_rhm_ID_fracs(r,&xfe2ID,&xmg2ID,&xmn2ID,&xti4ID,&xal3ID,&xfe3ID);
    double R = GH_rhm_R;
    double r0=r[0],r1=r[1],r2=r[2],r3=r[3], s0=s[0],s1=s[1],s2=s[2];
    double sumr = r0+r1+r2+r3;

    *dgdr0 = R*t*( 0.5*creal(clog(xfe2a)) + 0.5*creal(clog(xti4a)) + 0.5*creal(clog(xfe2b)) + 0.5*creal(clog(xti4b)) - creal(clog(xfe2ID)) - creal(clog(xti4ID)) )
           + (1.0-GH_rhm_SROconst)*2.0*R*t*(- creal(clog(xfe3ID)) + 0.5*creal(clog(xfe2ID)) + 0.5*creal(clog(xti4ID)) )
           + ((GH_rhm_whhmilm+(p-1.0)*GH_rhm_wvhmilm)+(GH_rhm_dwhhmilm+(p-1.0)*GH_rhm_dwvhmilm)*(1.0-2.0*r0-r1-r2-r3))*(1.0-2.0*r0-r1-r2-r3)
           - ((GH_rhm_whhmgei+(p-1.0)*GH_rhm_wvhmgei)+(GH_rhm_dwhhmgei+(p-1.0)*GH_rhm_dwvhmgei)*(1.0-r0-2.0*r1-r2-r3))*r1
           - ((GH_rhm_whhmpyr+(p-1.0)*GH_rhm_wvhmpyr)+(GH_rhm_dwhhmpyr+(p-1.0)*GH_rhm_dwvhmpyr)*(1.0-r0-r1-2.0*r2-r3))*r2
           - 2.0*(GH_rhm_dwhhmilm+(p-1.0)*GH_rhm_dwvhmilm)*r0*(1.0-sumr)
           -     (GH_rhm_dwhhmgei+(p-1.0)*GH_rhm_dwvhmgei)*r1*(1.0-sumr)
           -     (GH_rhm_dwhhmpyr+(p-1.0)*GH_rhm_dwvhmpyr)*r2*(1.0-sumr)
           -  GH_rhm_whmcrn*r3
           + (GH_rhm_whilmgei+(p-1.0)*GH_rhm_wvilmgei+GH_rhm_whilmgeiT)*r1/2.0
           + (GH_rhm_whilmpyr+(p-1.0)*GH_rhm_wvilmpyr+GH_rhm_whilmpyrT)*r2/2.0
           + (GH_rhm_wcrnilm+(GH_rhm_dwhcrnilm+(p-1.0)*GH_rhm_dwvcrnilm)*(r0-r3))*r3 + (GH_rhm_dwhcrnilm+(p-1.0)*GH_rhm_dwvcrnilm)*r0*r3
           + 2.0*(GH_rhm_dhilm+(p-1.0)*GH_rhm_dvilm+((GH_rhm_dwhhmilm+(p-1.0)*GH_rhm_dwvhmilm)/2.0+(GH_rhm_whhmilm2+(p-1.0)*GH_rhm_wvhmilm2)/4.0)*(1.0-sumr)
                      -(GH_rhm_dwhcrnilm+(p-1.0)*GH_rhm_dwvcrnilm)*r3/2.0)*r0
           - ((GH_rhm_dwhhmilm+(p-1.0)*GH_rhm_dwvhmilm)/2.0+(GH_rhm_whhmilm2+(p-1.0)*GH_rhm_wvhmilm2)/4.0)*(r0*r0-s0*s0)
           - ((GH_rhm_dwhhmgei+(p-1.0)*GH_rhm_dwvhmgei)/2.0+(GH_rhm_whhmgei2+(p-1.0)*GH_rhm_wvhmgei2)/4.0)*(r1*r1-s1*s1)
           - ((GH_rhm_dwhhmpyr+(p-1.0)*GH_rhm_dwvhmpyr)/2.0+(GH_rhm_whhmpyr2+(p-1.0)*GH_rhm_wvhmpyr2)/4.0)*(r2*r2-s2*s2)
           + 2.0*(GH_rhm_whilm+(p-1.0)*GH_rhm_wvilm)*s0*s0*r0;

    *dgdr1 = R*t*( 0.5*creal(clog(xmg2a)) + 0.5*creal(clog(xti4a)) + 0.5*creal(clog(xmg2b)) + 0.5*creal(clog(xti4b)) - creal(clog(xmg2ID)) - creal(clog(xti4ID)) )
           + (1.0-GH_rhm_SROconst)*2.0*R*t*(- creal(clog(xfe3ID)) + 0.5*creal(clog(xmg2ID)) + 0.5*creal(clog(xti4ID)) )
           - ((GH_rhm_whhmilm+(p-1.0)*GH_rhm_wvhmilm)+(GH_rhm_dwhhmilm+(p-1.0)*GH_rhm_dwvhmilm)*(1.0-2.0*r0-r1-r2-r3))*r0
           + ((GH_rhm_whhmgei+(p-1.0)*GH_rhm_wvhmgei)+(GH_rhm_dwhhmgei+(p-1.0)*GH_rhm_dwvhmgei)*(1.0-r0-2.0*r1-r2-r3))*(1.0-r0-2.0*r1-r2-r3)
           - ((GH_rhm_whhmpyr+(p-1.0)*GH_rhm_wvhmpyr)+(GH_rhm_dwhhmpyr+(p-1.0)*GH_rhm_dwvhmpyr)*(1.0-r0-r1-2.0*r2-r3))*r2
           -     (GH_rhm_dwhhmilm+(p-1.0)*GH_rhm_dwvhmilm)*r0*(1.0-sumr)
           - 2.0*(GH_rhm_dwhhmgei+(p-1.0)*GH_rhm_dwvhmgei)*r1*(1.0-sumr)
           -     (GH_rhm_dwhhmpyr+(p-1.0)*GH_rhm_dwvhmpyr)*r2*(1.0-sumr)
           -  GH_rhm_whmcrn*r3
           + (GH_rhm_whilmgei+(p-1.0)*GH_rhm_wvilmgei+GH_rhm_whilmgeiT)*r0/2.0
           + (GH_rhm_whgeipyr+(p-1.0)*GH_rhm_wvgeipyr+GH_rhm_whgeipyrT)*r2/2.0
           + (GH_rhm_wcrngei+(GH_rhm_dwhcrngei+(p-1.0)*GH_rhm_dwvcrngei)*(r1-r3))*r3 + (GH_rhm_dwhcrngei+(p-1.0)*GH_rhm_dwvcrngei)*r1*r3
           - ((GH_rhm_dwhhmilm+(p-1.0)*GH_rhm_dwvhmilm)/2.0+(GH_rhm_whhmilm2+(p-1.0)*GH_rhm_wvhmilm2)/4.0)*(r0*r0-s0*s0)
           + 2.0*(GH_rhm_dhgei+(p-1.0)*GH_rhm_dvgei+((GH_rhm_dwhhmgei+(p-1.0)*GH_rhm_dwvhmgei)/2.0+(GH_rhm_whhmgei2+(p-1.0)*GH_rhm_wvhmgei2)/4.0)*(1.0-sumr)
                      -(GH_rhm_dwhcrngei+(p-1.0)*GH_rhm_dwvcrngei)*r3/2.0)*r1
           - ((GH_rhm_dwhhmgei+(p-1.0)*GH_rhm_dwvhmgei)/2.0+(GH_rhm_whhmgei2+(p-1.0)*GH_rhm_wvhmgei2)/4.0)*(r1*r1-s1*s1)
           - ((GH_rhm_dwhhmpyr+(p-1.0)*GH_rhm_dwvhmpyr)/2.0+(GH_rhm_whhmpyr2+(p-1.0)*GH_rhm_wvhmpyr2)/4.0)*(r2*r2-s2*s2)
           + 2.0*(GH_rhm_whgei+(p-1.0)*GH_rhm_wvgei)*s1*s1*r1;

    *dgdr2 = R*t*( 0.5*creal(clog(xmn2a)) + 0.5*creal(clog(xti4a)) + 0.5*creal(clog(xmn2b)) + 0.5*creal(clog(xti4b)) - creal(clog(xmn2ID)) - creal(clog(xti4ID)) )
           + (1.0-GH_rhm_SROconst)*2.0*R*t*(- creal(clog(xfe3ID)) + 0.5*creal(clog(xmn2ID)) + 0.5*creal(clog(xti4ID)) )
           - ((GH_rhm_whhmilm+(p-1.0)*GH_rhm_wvhmilm)+(GH_rhm_dwhhmilm+(p-1.0)*GH_rhm_dwvhmilm)*(1.0-2.0*r0-r1-r2-r3))*r0
           - ((GH_rhm_whhmgei+(p-1.0)*GH_rhm_wvhmgei)+(GH_rhm_dwhhmgei+(p-1.0)*GH_rhm_dwvhmgei)*(1.0-r0-2.0*r1-r2-r3))*r1
           + ((GH_rhm_whhmpyr+(p-1.0)*GH_rhm_wvhmpyr)+(GH_rhm_dwhhmpyr+(p-1.0)*GH_rhm_dwvhmpyr)*(1.0-r0-r1-2.0*r2-r3))*(1.0-r0-r1-2.0*r2-r3)
           -     (GH_rhm_dwhhmilm+(p-1.0)*GH_rhm_dwvhmilm)*r0*(1.0-sumr)
           -     (GH_rhm_dwhhmgei+(p-1.0)*GH_rhm_dwvhmgei)*r1*(1.0-sumr)
           - 2.0*(GH_rhm_dwhhmpyr+(p-1.0)*GH_rhm_dwvhmpyr)*r2*(1.0-sumr)
           -  GH_rhm_whmcrn*r3
           + (GH_rhm_whilmpyr+(p-1.0)*GH_rhm_wvilmpyr+GH_rhm_whilmpyrT)*r0/2.0
           + (GH_rhm_whgeipyr+(p-1.0)*GH_rhm_wvgeipyr+GH_rhm_whgeipyrT)*r1/2.0
           + (GH_rhm_wcrnpyr+(GH_rhm_dwhcrnpyr+(p-1.0)*GH_rhm_dwvcrnpyr)*(r2-r3))*r3 + (GH_rhm_dwhcrnpyr+(p-1.0)*GH_rhm_dwvcrnpyr)*r2*r3
           - ((GH_rhm_dwhhmilm+(p-1.0)*GH_rhm_dwvhmilm)/2.0+(GH_rhm_whhmilm2+(p-1.0)*GH_rhm_wvhmilm2)/4.0)*(r0*r0-s0*s0)
           - ((GH_rhm_dwhhmgei+(p-1.0)*GH_rhm_dwvhmgei)/2.0+(GH_rhm_whhmgei2+(p-1.0)*GH_rhm_wvhmgei2)/4.0)*(r1*r1-s1*s1)
           + 2.0*(GH_rhm_dhpyr+(p-1.0)*GH_rhm_dvpyr+((GH_rhm_dwhhmpyr+(p-1.0)*GH_rhm_dwvhmpyr)/2.0+(GH_rhm_whhmpyr2+(p-1.0)*GH_rhm_wvhmpyr2)/4.0)*(1.0-sumr)
                      -(GH_rhm_dwhcrnpyr+(p-1.0)*GH_rhm_dwvcrnpyr)*r3/2.0)*r2
           - ((GH_rhm_dwhhmpyr+(p-1.0)*GH_rhm_dwvhmpyr)/2.0+(GH_rhm_whhmpyr2+(p-1.0)*GH_rhm_wvhmpyr2)/4.0)*(r2*r2-s2*s2)
           + 2.0*(GH_rhm_whpyr+(p-1.0)*GH_rhm_wvpyr)*s2*s2*r2;

    *dgdr3 = (1.0-GH_rhm_SROconst)*2.0*R*t*(- creal(clog(xfe3ID)) + creal(clog(xal3ID)) )
           - ((GH_rhm_whhmilm+(p-1.0)*GH_rhm_wvhmilm)+(GH_rhm_dwhhmilm+(p-1.0)*GH_rhm_dwvhmilm)*(1.0-2.0*r0-r1-r2-r3))*r0
           - ((GH_rhm_whhmgei+(p-1.0)*GH_rhm_wvhmgei)+(GH_rhm_dwhhmgei+(p-1.0)*GH_rhm_dwvhmgei)*(1.0-r0-2.0*r1-r2-r3))*r1
           - ((GH_rhm_whhmpyr+(p-1.0)*GH_rhm_wvhmpyr)+(GH_rhm_dwhhmpyr+(p-1.0)*GH_rhm_dwvhmpyr)*(1.0-r0-r1-2.0*r2-r3))*r2
           - (GH_rhm_dwhhmilm+(p-1.0)*GH_rhm_dwvhmilm)*r0*(1.0-sumr)
           - (GH_rhm_dwhhmgei+(p-1.0)*GH_rhm_dwvhmgei)*r1*(1.0-sumr)
           - (GH_rhm_dwhhmpyr+(p-1.0)*GH_rhm_dwvhmpyr)*r2*(1.0-sumr)
           +  GH_rhm_whmcrn*(1.0-r0-r1-r2-2.0*r3)
           + (GH_rhm_wcrnilm+(GH_rhm_dwhcrnilm+(p-1.0)*GH_rhm_dwvcrnilm)*(r0-r3))*r0 - (GH_rhm_dwhcrnilm+(p-1.0)*GH_rhm_dwvcrnilm)*r0*r3
           + (GH_rhm_wcrngei+(GH_rhm_dwhcrngei+(p-1.0)*GH_rhm_dwvcrngei)*(r1-r3))*r1 - (GH_rhm_dwhcrngei+(p-1.0)*GH_rhm_dwvcrngei)*r1*r3
           + (GH_rhm_wcrnpyr+(GH_rhm_dwhcrnpyr+(p-1.0)*GH_rhm_dwvcrnpyr)*(r2-r3))*r2 - (GH_rhm_dwhcrnpyr+(p-1.0)*GH_rhm_dwvcrnpyr)*r2*r3
           - ((GH_rhm_dwhhmilm+(p-1.0)*GH_rhm_dwvhmilm)/2.0 + (GH_rhm_whhmilm2+(p-1.0)*GH_rhm_wvhmilm2)/4.0 + (GH_rhm_dwhcrnilm+(p-1.0)*GH_rhm_dwvcrnilm)/2.0)*(r0*r0-s0*s0)
           - ((GH_rhm_dwhhmgei+(p-1.0)*GH_rhm_dwvhmgei)/2.0 + (GH_rhm_whhmgei2+(p-1.0)*GH_rhm_wvhmgei2)/4.0 + (GH_rhm_dwhcrngei+(p-1.0)*GH_rhm_dwvcrngei)/2.0)*(r1*r1-s1*s1)
           - ((GH_rhm_dwhhmpyr+(p-1.0)*GH_rhm_dwvhmpyr)/2.0 + (GH_rhm_whhmpyr2+(p-1.0)*GH_rhm_wvhmpyr2)/4.0 + (GH_rhm_dwhcrnpyr+(p-1.0)*GH_rhm_dwvcrnpyr)/2.0)*(r2*r2-s2*s2);
}
/** Duplicated, file-local copies of GH_gem_function.c's own
    GH_rhm_pure_order_G/GH_rhm_SRO_G (same formulas, same constants -
    matching this file's own established convention of keeping its rhm
    helpers, e.g. GH_rhm_H/GH_rhm_S/GH_rhm_dGdr above, self-contained
    rather than sharing statics across files). Needed to subtract
    ENDMEMBERS(r) = sum(p_i * that endmember's own individual ordering
    correction) from Gex below - see obj_gh_rhm's header comment for why
    this subtraction was missing (a real bug, found and fixed 2026-07-14):
    gh's gbase[] already carries each endmember's own correction via these
    same two functions (GH_gem_function.c), so obj_gh_rhm reporting the
    full G(r,s*) - ENDMEMBERS(r) mixing-only quantity (matching real
    MELTS' gmixMsg's own "*gmix -= ENDMEMBERS") is required to avoid
    double-counting it, exactly the pattern cpx's Emix/es_end already
    uses. Verified against real xMELTS' gmixMsg (gexcess_verification/). */
static double GH_rhm_mix_endmember_order_G(double T, double P, double dh, double wh, double dv, double wv){
    double R = 8.3143;
    double eps = 2.220446049250313e-16; /* DBL_EPSILON */
    double s = 0.98;
    for (int iter = 0; iter < 1000; iter++){
        double dgds   = R*T*(creal(clog(1.0+s)) - creal(clog(1.0-s)))
                      - 2.0*(dh+(P-1.0)*dv)*s + (wh+(P-1.0)*wv)*(2.0*s - 4.0*s*s*s);
        double d2gds2 = R*T*(1.0/(1.0+s) + 1.0/(1.0-s))
                      - 2.0*(dh+(P-1.0)*dv) + (wh+(P-1.0)*wv)*(2.0 - 12.0*s*s);
        double sNew   = s - dgds/d2gds2;
        sNew = (sNew > 1.0 - 10.0*eps) ? 1.0 - 10.0*eps : sNew;
        sNew = (sNew < 0.0) ? 0.0 : sNew;
        int converged = (fabs(sNew - s) <= 10.0*eps);
        s = sNew;
        if (converged) break;
    }
    double S = -R*((1.0+s)*creal(clog(1.0+s)) + (1.0-s)*creal(clog(1.0-s)) - 2.0*creal(clog(2.0)));
    double H = (dh+(P-1.0)*dv)*(1.0 - s*s) + (wh+(P-1.0)*wv)*s*s*(1.0 - s*s);
    return H - T*S;
}
static double GH_rhm_mix_endmember_SRO_G(double T){
    const double SROconst = 0.0730205;
    double R = 8.3143;
    return -T*(SROconst*2.0*R*creal(clog(2.0)));
}

double obj_gh_rhm(unsigned n, const double *x, double *grad, void *SS_ref_db){
    SS_ref *d = (SS_ref *) SS_ref_db;

    double T    = d->T;
    double Pbar = d->P*1000.0;   /* kbar -> bar, matching xMELTS' own P=Pr=1bar reference */
    double *p   = d->p;
    double *gb  = d->gb_lvl;
    double *mu_Gex = d->mu_Gex;

    for (int i = 0; i < 5; i++){ p[i] = x[i]; }
    /* r0=Xilm=p[2], r1=Xgei=p[0], r2=Xpyr=p[3], r3=Xcrn=p[4]; p[1]=Xhem dependent */
    double r[4] = { p[2], p[0], p[3], p[4] };

    double s0, s1, s2;
    GH_rhm_solve_s_from(r, T, Pbar, &s0, &s1, &s2);
    double s[3] = { s0, s1, s2 };

    double Graw = GH_rhm_H(r, s, Pbar) - T*GH_rhm_S(r, s, T);

    double dgdr0, dgdr1, dgdr2, dgdr3;
    GH_rhm_dGdr(r, s, T, Pbar, &dgdr0, &dgdr1, &dgdr2, &dgdr3);

    /* ENDMEMBERS(r): gei/ilm/pyr share the same ordering correction
       (dh/wh/dv/wv identical across the three in xMELTS' calibration,
       matching GH_gem_function.c's own gbase construction); hem/crn share
       the flat-SRO-spline correction. p[1]=hem is genuinely present in
       ENDMEMBERS (unlike in Graw, which only depends on r0..r3), so its
       own dGex contribution is nonzero here even though it's the
       dependent endmember in the r-basis. */
    double corr_ord = GH_rhm_mix_endmember_order_G(T, Pbar, 17477.0, 3189.0, 0.010758, 0.035089);
    double corr_sro = GH_rhm_mix_endmember_SRO_G(T);
    double ENDMEMBERS = (p[0]+p[2]+p[3])*corr_ord + (p[1]+p[4])*corr_sro;
    double Gex = Graw - ENDMEMBERS;

    double dGex[5];
    dGex[0] = dgdr1 - corr_ord;   /* gei */
    dGex[1] = -corr_sro;         /* hem, dependent in r-basis but present in ENDMEMBERS */
    dGex[2] = dgdr0 - corr_ord;   /* ilm */
    dGex[3] = dgdr2 - corr_ord;   /* pyr */
    dGex[4] = dgdr3 - corr_sro;   /* crn */

    for (int i = 0; i < 5; i++){
        double mu = Gex;
        for (int k = 0; k < 5; k++){
            double delta_ik = (i==k) ? 1.0 : 0.0;
            mu += (delta_ik - p[k])*dGex[k];
        }
        mu_Gex[i] = mu/1000.0;   /* J -> kJ */
        d->sf[i]  = p[i];
    }

    d->sum_apep = 0.0;
    for (int i = 0; i < 5; i++){ d->sum_apep += d->ape[i]*p[i]; }
    d->factor = d->fbc/d->sum_apep;

    d->df_raw = 0.0;
    for (int i = 0; i < 5; i++){ d->df_raw += (mu_Gex[i] + gb[i])*p[i]; }
    d->df = d->df_raw * d->factor;

    if (grad){
        for (int i = 0; i < 5; i++){
            grad[i] = (mu_Gex[i] + gb[i])*d->factor - (d->df_raw*d->factor*(d->ape[i]/d->sum_apep));
        }
    }
    return d->df;
}

/**
    Nepheline (Sack & Ghiorso 1995): from xMELTS' sources/nepheline.c.
    NR=3 (r0=X(k-nepheline), r1=X(vc-nepheline), r2=X(ca-nepheline),
    Xna-nepheline=1-r0-r1-r2 dependent), NS=1 (single K/Na large-vs-small
    site-ordering parameter s0). p[]/EM_list order (nane,kne,vcne,cane) is
    a direct r-basis relabeling: p[0]=Xna-nepheline dependent, p[1]=r0,
    p[2]=r1, p[3]=r2 - same pattern rhm/spn already use.

    W-parameter table ported verbatim (names/values) from nepheline.c's
    file-scope #defines. Of the 15 originally-named W/A parameters, 9
    (DWKNALS, DWKNASS, DWVNASS, DWVNALS, DWVKLS, AWLS, AWSS, AWVNA, AWVK)
    and 3 more (WCANA, WCAK, WPLAG) are all 0.0 in xMELTS' own calibration -
    the (large) polynomial terms they multiply in the source's H/DGDR/
    D2GDR2 macros are simply omitted below rather than carried as always-0
    additions (same simplification approach used for rhm's 6 zero cross-
    terms).

    The order-parameter solve (GH_nph_solve_s0) ports xMELTS' own order()
    function (nepheline.c) as faithfully as NS=1 allows: a Newton step each
    iteration (trivial 1x1 "matrix solve", no Cramer's rule needed), then a
    damped line-search (lambda halved) against ALL SIX site fractions'
    [0,1] (or [0,1/3] for xcass) feasibility bounds simultaneously - unlike
    rhm's simple box-clamp, this matches nepheline.c's own algorithm
    exactly, including its distinctive "re-project s from the (possibly
    clamped) site fractions after each step" behavior (s = xkls-xkss,
    recomputed from clamped site fractions rather than used directly from
    the Newton step - see nepheline.c order(), lines ~874-909).

    A PENALTY/r1 barrier term (PENALTY=sqrt(DBL_EPSILON)) is added to G to
    keep the vc-nepheline mole fraction away from exactly 0, ported
    unconditionally from nepheline.c's own #define G (no DBL_MAX special
    case needed here, unlike kalsilite.c's copy of the same guard - gh's
    NLopt bounds already keep every xeos, including r1=p[2], at or above
    eps > 0, so r1 is never exactly 0 in gh's own evaluation).

    Total molar G = sum(p[i]*gbase[i]) + Gex(r,s0*(r)), where gbase[i]
    already carries each endmember's own standard-state correction (na-
    nepheline/k-nepheline plain Berman; vc-nepheline/ca-nepheline their
    real Vinet-EOS "phantom reaction" correction - see GH_vcneph_G/
    GH_caneph_G in GH_gem_function.c) - matching xMELTS' own construction
    exactly, since nepheline.c's H/S/V macros contain NO endmember-linear
    terms (DH1..DH4 are defined in nepheline.c but never actually used
    there - confirmed by direct inspection - they're purely informational,
    matching the values baked into the endmember table's own H ref
    entries instead). dGex/dp uses the envelope theorem (DGDR0..2 at the
    converged s0*, no ds/dr needed) plus the same delta_ik-p_k tangent-
    plane trick used throughout gh.
*/
static const double GH_nph_HEX=-31744.63, GH_nph_SEX=-20.92, GH_nph_VEX=-1.05;
static const double GH_nph_HX=-13893.18, GH_nph_SX=12.55, GH_nph_VX=0.0;
static const double GH_nph_H23=73520.0, GH_nph_H24=42560.0, GH_nph_S23=0.0, GH_nph_S24=0.0;
static const double GH_nph_WHNAKLS=6861.76, GH_nph_WVNAKLS=0.33, GH_nph_WHNAKSS=51002.96, GH_nph_WVNAKSS=0.54;
static const double GH_nph_WVN=14644.0, GH_nph_WVK=8368.0;
static const double GH_nph_S4=-15.8765;
static const double GH_nph_R = 8.3143;
static const double GH_nph_PENALTY = 1.4901161193847656e-08; /* sqrt(DBL_EPSILON) */

static void GH_nph_site_fracs(double r0, double r1, double r2, double s0,
    double *xkls, double *xvcls, double *xnals, double *xkss, double *xcass, double *xnass){
    double eps = 2.220446049250313e-16;
    double v;
    v = r0 + 3.0*s0/4.0;                    *xkls  = (v<eps)?eps:((v>1.0-eps)?1.0-eps:v);
    v = r1 + r2;                            *xvcls = (v<eps)?eps:((v>1.0-eps)?1.0-eps:v);
    v = 1.0 - r0 - r1 - r2 - 3.0*s0/4.0;    *xnals = (v<eps)?eps:((v>1.0-eps)?1.0-eps:v);
    v = r0 - s0/4.0;                        *xkss  = (v<eps)?eps:((v>1.0-eps)?1.0-eps:v);
    v = r2/3.0;                             *xcass = (v<eps)?eps:((v>1.0/3.0-eps)?1.0/3.0-eps:v);
    v = 1.0 - r0 - r2/3.0 + s0/4.0;         *xnass = (v<eps)?eps:((v>1.0-eps)?1.0-eps:v);
}
static double GH_nph_dgds0(double r0, double r1, double r2, double s0, double t, double p){
    double xkls,xvcls,xnals,xkss,xcass,xnass;
    GH_nph_site_fracs(r0,r1,r2,s0,&xkls,&xvcls,&xnals,&xkss,&xcass,&xnass);
    double R = GH_nph_R;
    double GEX = GH_nph_HEX - t*GH_nph_SEX + (p-1.0)*GH_nph_VEX;
    double GX  = GH_nph_HX  - t*GH_nph_SX  + (p-1.0)*GH_nph_VX;
    double WNAKLS = GH_nph_WHNAKLS + (p-1.0)*GH_nph_WVNAKLS;
    double WNAKSS = GH_nph_WHNAKSS + (p-1.0)*GH_nph_WVNAKSS;
    double G23 = GH_nph_H23 - t*GH_nph_S23;

    return R*t*(0.75*creal(clog(xkls)) - 0.75*creal(clog(xnals)) - 0.75*creal(clog(xkss)) + 0.75*creal(clog(xnass)))
         + 0.25*(2.0*GEX+GX+3.0*WNAKLS-WNAKSS)
         + (1.0/8.0)*(3.0*GX-9.0*WNAKLS-WNAKSS)*s0
         - 0.5*(GX+3.0*WNAKLS-WNAKSS)*r0
         + 0.125*(GX-GEX-2.0*G23-6.0*WNAKLS-6.0*GH_nph_WVN+6.0*GH_nph_WVK)*r1
         + (18.0*GH_nph_WVN-18.0*GH_nph_WVK-18.0*WNAKLS+4.0*WNAKSS-9.0*GEX-3.0*GX+6.0*G23-36.0*(GH_nph_H24-t*GH_nph_S24))*r2/24.0;
}
static double GH_nph_d2gds0s0(double r0, double r1, double r2, double s0, double t, double p){
    double xkls,xvcls,xnals,xkss,xcass,xnass;
    GH_nph_site_fracs(r0,r1,r2,s0,&xkls,&xvcls,&xnals,&xkss,&xcass,&xnass);
    double R = GH_nph_R;
    double GX  = GH_nph_HX  - t*GH_nph_SX  + (p-1.0)*GH_nph_VX;
    double WNAKLS = GH_nph_WHNAKLS + (p-1.0)*GH_nph_WVNAKLS;
    double WNAKSS = GH_nph_WHNAKSS + (p-1.0)*GH_nph_WVNAKSS;

    return R*t*(0.75*0.75/xkls + 0.75*0.75/xnals + 3.0*0.25*0.25/xkss + 3.0*0.25*0.25/xnass)
         + (1.0/8.0)*(3.0*GX-9.0*WNAKLS-WNAKSS);
}
/* Joint (here: scalar, NS=1) order-parameter solve, ported directly from
   xMELTS' own order() function (nepheline.c) - see this block's header
   comment for why the damped line-search-against-site-fractions and the
   "re-project s from the clamped site fractions" step must both be
   replicated exactly, not simplified to a plain box-clamp. */
static double GH_nph_solve_s0(double r0, double r1, double r2, double t, double p){
    const int MAX_ITER = 200; /* xMELTS' own MAX_ITER (nepheline.c) */
    double eps = 2.220446049250313e-16;

    double sNew = -4.0*r0*(r1+2.0*r2/3.0)/(4.0-r1-2.0*r2);
    double xkls,xvcls,xnals,xkss,xcass,xnass;
    GH_nph_site_fracs(r0,r1,r2,sNew,&xkls,&xvcls,&xnals,&xkss,&xcass,&xnass);
    sNew = xkls - xkss;

    double sOld = 2.0;
    int iter = 0;
    while ((fabs(sNew-sOld) > 10.0*eps) && (iter < MAX_ITER)){
        double s0 = sNew;
        double dgds  = GH_nph_dgds0(r0,r1,r2,s0,t,p);
        double d2gds2 = GH_nph_d2gds0s0(r0,r1,r2,s0,t,p);
        sOld = s0;

        double sCorr = -dgds/d2gds2;
        double lambda = 1.0;
        double sTrial = sOld + lambda*sCorr;
        double xkls_,xvcls_,xnals_,xkss_,xcass_,xnass_;
        GH_nph_site_fracs(r0,r1,r2,sTrial,&xkls_,&xvcls_,&xnals_,&xkss_,&xcass_,&xnass_);
        /* re-derive the UNCLAMPED site fractions (site_fracs already
           clamps) to test feasibility, matching the source's own
           unclamped xkls/xvcls/... comparison in its while() guard */
        double rxkls  = r0 + 3.0*sTrial/4.0;
        double rxvcls = r1 + r2;
        double rxnals = 1.0 - r0 - r1 - r2 - 3.0*sTrial/4.0;
        double rxkss  = r0 - sTrial/4.0;
        double rxcass = r2/3.0;
        double rxnass = 1.0 - r0 - r2/3.0 + sTrial/4.0;
        while ((rxkls<0.0 || rxkls>1.0 || rxvcls<0.0 || rxvcls>1.0 ||
                rxnals<0.0 || rxnals>1.0 || rxkss<0.0  || rxkss>1.0 ||
                rxcass<0.0 || rxcass>1.0/3.0 || rxnass<0.0 || rxnass>1.0)
               && lambda > eps){
            lambda /= 2.0;
            sTrial = sOld + lambda*sCorr;
            rxkls  = r0 + 3.0*sTrial/4.0;
            rxvcls = r1 + r2;
            rxnals = 1.0 - r0 - r1 - r2 - 3.0*sTrial/4.0;
            rxkss  = r0 - sTrial/4.0;
            rxcass = r2/3.0;
            rxnass = 1.0 - r0 - r2/3.0 + sTrial/4.0;
        }

        GH_nph_site_fracs(r0,r1,r2,sTrial,&xkls,&xvcls,&xnals,&xkss,&xcass,&xnass);
        sNew = xkls - xkss;
        iter++;
    }
    return sNew;
}
static double GH_nph_S(double r0, double r1, double r2, double s0, double t){
    double xkls,xvcls,xnals,xkss,xcass,xnass;
    GH_nph_site_fracs(r0,r1,r2,s0,&xkls,&xvcls,&xnals,&xkss,&xcass,&xnass);
    double R = GH_nph_R;

    return GH_nph_S4*r2 - R*(xkls*creal(clog(xkls)) + xvcls*creal(clog(xvcls)) +
             xnals*creal(clog(xnals)) - (1.0-3.0*xcass)*creal(clog(1.0-3.0*xcass)) +
             3.0*xkss*creal(clog(xkss)) + 3.0*xcass*creal(clog(xcass)) + 3.0*xnass*creal(clog(xnass)))
         + 0.25*(2.0*GH_nph_SEX+GH_nph_SX)*s0 + GH_nph_SX*r0*(1.0-r0)
         + (3.0/16.0)*GH_nph_SX*s0*s0 - 0.5*GH_nph_SX*r0*s0
         + 0.5*(2.0*GH_nph_S23+GH_nph_SEX-GH_nph_SX)*r0*r1 + 0.125*(GH_nph_SX-GH_nph_SEX-2.0*GH_nph_S23)*r1*s0
         + (6.0*GH_nph_S23+3.0*GH_nph_SEX-15.0*GH_nph_SX)*r0*r2/18.0
         + (-9.0*GH_nph_SEX-3.0*GH_nph_SX+6.0*GH_nph_S23-36.0*GH_nph_S24)*r2*s0/24.0;
}
static double GH_nph_H(double r0, double r1, double r2, double s0, double p){
    return 0.25*(2.0*GH_nph_HEX+GH_nph_HX+3.0*GH_nph_WHNAKLS-GH_nph_WHNAKSS)*s0
         + (GH_nph_HX+GH_nph_WHNAKLS+GH_nph_WHNAKSS)*r0*(1.0-r0)
         + (1.0/16.0)*(3.0*GH_nph_HX-9.0*GH_nph_WHNAKLS-GH_nph_WHNAKSS)*s0*s0
         - 0.5*(GH_nph_HX+3.0*GH_nph_WHNAKLS-GH_nph_WHNAKSS)*r0*s0 + GH_nph_WVN*r1*(1.0-r1)
         + 0.5*(2.0*GH_nph_H23+GH_nph_HEX-GH_nph_HX-2.0*GH_nph_WHNAKLS-2.0*GH_nph_WVN+2.0*GH_nph_WVK)*r0*r1
         + 0.125*(GH_nph_HX-GH_nph_HEX-2.0*GH_nph_H23-6.0*GH_nph_WHNAKLS-6.0*GH_nph_WVN+6.0*GH_nph_WVK)*r1*s0
         + (6.0*GH_nph_H23+3.0*GH_nph_HEX-15.0*GH_nph_HX-18.0*GH_nph_WHNAKLS-4.0*GH_nph_WHNAKSS+18.0*GH_nph_WVN-18.0*GH_nph_WVK)*r0*r2/18.0
         + (-GH_nph_WVN)*r1*r2
         + (18.0*GH_nph_WVN-18.0*GH_nph_WVK-18.0*GH_nph_WHNAKLS+4.0*GH_nph_WHNAKSS-9.0*GH_nph_HEX-3.0*GH_nph_HX+6.0*GH_nph_H23-36.0*GH_nph_H24)*r2*s0/24.0;
}
static double GH_nph_V(double r0, double r1, double r2, double s0, double p){
    return 0.25*(2.0*GH_nph_VEX+GH_nph_VX+3.0*GH_nph_WVNAKLS-GH_nph_WVNAKSS)*s0
         + (GH_nph_VX+GH_nph_WVNAKLS+GH_nph_WVNAKSS)*r0*(1.0-r0)
         + (1.0/16.0)*(3.0*GH_nph_VX-9.0*GH_nph_WVNAKLS-GH_nph_WVNAKSS)*s0*s0
         - 0.5*(GH_nph_VX+3.0*GH_nph_WVNAKLS-GH_nph_WVNAKSS)*r0*s0
         + 0.5*(GH_nph_VEX-GH_nph_VX-2.0*GH_nph_WVNAKLS)*r0*r1
         + 0.125*(GH_nph_VX-GH_nph_VEX-6.0*GH_nph_WVNAKLS)*r1*s0
         + (3.0*GH_nph_VEX-15.0*GH_nph_VX-18.0*GH_nph_WVNAKLS-4.0*GH_nph_WVNAKSS)*r0*r2/18.0
         + (-18.0*GH_nph_WVNAKLS+4.0*GH_nph_WVNAKSS-9.0*GH_nph_VEX-3.0*GH_nph_VX)*r2*s0/24.0;
}
static void GH_nph_dGdr(double r0, double r1, double r2, double s0, double t, double p,
                         double *dgdr0, double *dgdr1, double *dgdr2){
    double xkls,xvcls,xnals,xkss,xcass,xnass;
    GH_nph_site_fracs(r0,r1,r2,s0,&xkls,&xvcls,&xnals,&xkss,&xcass,&xnass);
    double R = GH_nph_R;
    double GEX = GH_nph_HEX - t*GH_nph_SEX + (p-1.0)*GH_nph_VEX;
    double GX  = GH_nph_HX  - t*GH_nph_SX  + (p-1.0)*GH_nph_VX;
    double WNAKLS = GH_nph_WHNAKLS + (p-1.0)*GH_nph_WVNAKLS;
    double WNAKSS = GH_nph_WHNAKSS + (p-1.0)*GH_nph_WVNAKSS;
    double G23 = GH_nph_H23 - t*GH_nph_S23;
    double G24 = GH_nph_H24 - t*GH_nph_S24;
    double G4  = -t*GH_nph_S4;

    *dgdr0 = R*t*(creal(clog(xkls)) - creal(clog(xnals)) + 3.0*creal(clog(xkss)) - 3.0*creal(clog(xnass)))
           + (GX+WNAKLS+WNAKSS)*(1.0-2.0*r0)
           - 0.5*(GX+3.0*WNAKLS-WNAKSS)*s0
           + 0.5*(2.0*G23+GEX-GX-2.0*WNAKLS-2.0*GH_nph_WVN+2.0*GH_nph_WVK)*r1
           + (6.0*G23+3.0*GEX-15.0*GX-18.0*WNAKLS-4.0*WNAKSS+18.0*GH_nph_WVN-18.0*GH_nph_WVK)*r2/18.0;

    *dgdr1 = R*t*(creal(clog(xvcls)) - creal(clog(xnals))) + GH_nph_WVN*(1.0-2.0*r1)
           + 0.5*(2.0*G23+GEX-GX-2.0*WNAKLS-2.0*GH_nph_WVN+2.0*GH_nph_WVK)*r0
           + 0.125*(GX-GEX-2.0*G23-6.0*WNAKLS-6.0*GH_nph_WVN+6.0*GH_nph_WVK)*s0
           + (-GH_nph_WVN)*r2
           - GH_nph_PENALTY/(r1*r1);

    *dgdr2 = G4 + R*t*(creal(clog(xvcls)) - creal(clog(xnals)) + creal(clog(1.0-3.0*xcass)) + creal(clog(xcass)) - creal(clog(xnass)) + 1.0)
           + (6.0*G23+3.0*GEX-15.0*GX-18.0*WNAKLS-4.0*WNAKSS+18.0*GH_nph_WVN-18.0*GH_nph_WVK)*r0/18.0
           + (-GH_nph_WVN)*r1
           + (18.0*GH_nph_WVN-18.0*GH_nph_WVK-18.0*WNAKLS+4.0*WNAKSS-9.0*GEX-3.0*GX+6.0*G23-36.0*G24)*s0/24.0;
}
double obj_gh_nph(unsigned n, const double *x, double *grad, void *SS_ref_db){
    SS_ref *d = (SS_ref *) SS_ref_db;

    double T    = d->T;
    double Pbar = d->P*1000.0;
    double *p   = d->p;
    double *gb  = d->gb_lvl;
    double *mu_Gex = d->mu_Gex;

    for (int i = 0; i < 4; i++){ p[i] = x[i]; }
    double r0 = p[1], r1 = p[2], r2 = p[3];

    double s0 = GH_nph_solve_s0(r0, r1, r2, T, Pbar);

    double Gex = GH_nph_H(r0,r1,r2,s0,Pbar) - T*GH_nph_S(r0,r1,r2,s0,T) + (Pbar-1.0)*GH_nph_V(r0,r1,r2,s0,Pbar)
               + GH_nph_PENALTY/r1;

    double dgdr0, dgdr1, dgdr2;
    GH_nph_dGdr(r0, r1, r2, s0, T, Pbar, &dgdr0, &dgdr1, &dgdr2);

    double dGex[4];
    dGex[0] = 0.0;    /* na-nepheline, dependent endmember */
    dGex[1] = dgdr0;
    dGex[2] = dgdr1;
    dGex[3] = dgdr2;

    for (int i = 0; i < 4; i++){
        double mu = Gex;
        for (int k = 0; k < 4; k++){
            double delta_ik = (i==k) ? 1.0 : 0.0;
            mu += (delta_ik - p[k])*dGex[k];
        }
        mu_Gex[i] = mu/1000.0;   /* J -> kJ */
        d->sf[i]  = p[i];
    }

    d->sum_apep = 0.0;
    for (int i = 0; i < 4; i++){ d->sum_apep += d->ape[i]*p[i]; }
    d->factor = d->fbc/d->sum_apep;

    d->df_raw = 0.0;
    for (int i = 0; i < 4; i++){ d->df_raw += (mu_Gex[i] + gb[i])*p[i]; }
    d->df = d->df_raw * d->factor;

    if (grad){
        for (int i = 0; i < 4; i++){
            grad[i] = (mu_Gex[i] + gb[i])*d->factor - (d->df_raw*d->factor*(d->ape[i]/d->sum_apep));
        }
    }
    return d->df;
}

/**
    Kalsilite (Sack & Ghiorso 1995): from xMELTS' sources/kalsilite.c.
    Shares nepheline's NR=3 r-basis (r0=X(k-nepheline "knk"), r1=X(vc-
    nepheline), r2=X(ca-nepheline)) and its PENALTY/r1 barrier, but has NO
    order parameter at all (kalsilite.c defines no NS/order() function) -
    a direct algebraic ideal+regular solution over 4 sites (xk,xvc,xna,
    xca), all plain functions of r with no embedded solve needed.
*/
static void GH_kls_site_fracs(double r0, double r1, double r2,
    double *xk, double *xvc, double *xna, double *xca){
    double eps = 2.220446049250313e-16;
    double v;
    v = r0;                            *xk  = (v<eps)?eps:((v>1.0-eps)?1.0-eps:v);
    v = (r1+r2)/4.0;                   *xvc = (v<eps)?eps:((v>1.0-eps)?1.0-eps:v);
    v = 1.0 - r0 - r1/4.0 - r2/2.0;    *xna = (v<eps)?eps:((v>1.0-eps)?1.0-eps:v);
    v = r2/4.0;                        *xca = (v<eps)?eps:((v>1.0-eps)?1.0-eps:v);
}
double obj_gh_kls(unsigned n, const double *x, double *grad, void *SS_ref_db){
    SS_ref *d = (SS_ref *) SS_ref_db;

    const double DH1=32145.67, DS1=10.46, DV1=0.0;
    const double DH2=-5648.4,  DS2=0.0,   DV2=-0.04;
    const double DH3=14644.0,  DS3=-18.7, DV3=0.0;
    const double DH4=35184.0,  DS4=-7.17, DV4=0.0;
    const double WKALS=29288.0, WKVKALS=30334.0, WNVKALS=14646.0, WKCKALS=14644.0;
    const double S4=-15.8765;
    const double R = 8.3143;
    const double PENALTY = 1.4901161193847656e-08; /* sqrt(DBL_EPSILON) */

    double T    = d->T;
    double Pbar = d->P*1000.0;
    double *p   = d->p;
    double *gb  = d->gb_lvl;
    double *mu_Gex = d->mu_Gex;

    for (int i = 0; i < 4; i++){ p[i] = x[i]; }
    double r0 = p[1], r1 = p[2], r2 = p[3];
    double sumr = r0+r1+r2;

    double xk,xvc,xna,xca;
    GH_kls_site_fracs(r0,r1,r2,&xk,&xvc,&xna,&xca);

    double S = DS1*(1.0-sumr) + DS2*r0 + DS3*r1 + (DS4+S4)*r2
             - 4.0*R*(xk*creal(clog(xk)) + (xvc-xca)*creal(clog(xvc-xca)) + xna*creal(clog(xna)) + xca*creal(clog(xca)));
    double H = DH1*(1.0-sumr) + DH2*r0 + DH3*r1 + DH4*r2
             + WKALS*r0*(1.0-sumr) + WNVKALS*r1*(1.0-sumr) + WKVKALS*r0*r1 + WKCKALS*r0*r2;
    double V = DV1*(1.0-sumr) + DV2*r0 + DV3*r1 + DV4*r2;
    double G4 = -T*S4;

    double Gex = H - T*S + (Pbar-1.0)*V + PENALTY/r1;

    double dgdr0 = - (DH1-T*DS1+(Pbar-1.0)*DV1) + (DH2-T*DS2+(Pbar-1.0)*DV2) + R*T*4.0*(creal(clog(xk))-creal(clog(xna)))
                 + WKALS*(1.0-2.0*r0-r1-r2) - WNVKALS*r1 + WKVKALS*r1 + WKCKALS*r2;
    double dgdr1 = - (DH1-T*DS1+(Pbar-1.0)*DV1) + (DH3-T*DS3+(Pbar-1.0)*DV3) + R*T*(creal(clog(xvc-xca))-creal(clog(xna)))
                 - WKALS*r0 + WNVKALS*(1.0-r0-2.0*r1-r2) + WKVKALS*r0
                 - PENALTY/(r1*r1);
    double dgdr2 = - (DH1-T*DS1+(Pbar-1.0)*DV1) + (DH4-T*DS4+(Pbar-1.0)*DV4) + G4
                 + R*T*(creal(clog(xca))-2.0*creal(clog(xna))-1.0)
                 - WKALS*r0 - WNVKALS*r1 + WKCKALS*r0;

    double dGex[4];
    dGex[0] = 0.0;    /* na-nepheline, dependent endmember */
    dGex[1] = dgdr0;
    dGex[2] = dgdr1;
    dGex[3] = dgdr2;

    for (int i = 0; i < 4; i++){
        double mu = Gex;
        for (int k = 0; k < 4; k++){
            double delta_ik = (i==k) ? 1.0 : 0.0;
            mu += (delta_ik - p[k])*dGex[k];
        }
        mu_Gex[i] = mu/1000.0;   /* J -> kJ */
        d->sf[i]  = p[i];
    }

    d->sum_apep = 0.0;
    for (int i = 0; i < 4; i++){ d->sum_apep += d->ape[i]*p[i]; }
    d->factor = d->fbc/d->sum_apep;

    d->df_raw = 0.0;
    for (int i = 0; i < 4; i++){ d->df_raw += (mu_Gex[i] + gb[i])*p[i]; }
    d->df = d->df_raw * d->factor;

    if (grad){
        for (int i = 0; i < 4; i++){
            grad[i] = (mu_Gex[i] + gb[i])*d->factor - (d->df_raw*d->factor*(d->ape[i]/d->sum_apep));
        }
    }
    return d->df;
}

/**
    Clinopyroxene / Orthopyroxene (Di-Cen-Hed-CaTs(Al,"alumino-buffonite")-
    CaTs(Fe3+,"buffonite")-Ess-Jd): from xMELTS' sources/clinopyroxene.c
    and orthopyroxene.c (Sack & Ghiorso 1993). NR=6, NS=2 - gh's largest
    embedded order-parameter phase (Step 5, the last planned addition).

    KEY FINDING (confirmed by diffing the two ~5400-line source files,
    not assumed): clinopyroxene.c and orthopyroxene.c are the SAME model.
    `ISCLINO` (the compile flag that would let xMELTS' own isClino()
    compare a clino-vs-ortho G and pick the lower) is `#undef`'d
    unconditionally in both files, and BOTH files separately hardcode
    `static const int clino = TRUE;` in their own order()/pureOrder()
    functions - so both registered MELTS phases ("clinopyroxene" and
    "orthopyroxene") always evaluate the identical clino-branch Taylor
    coefficients. The only differences between the two files are
    testCpx/testOpx's composition-validity bounds and conCpx/conOpx's
    default Fe2+/Fe3+ probe-analysis partitioning - neither relevant to
    gh's direct-composition p=x approach. So gh registers "cpx" and "opx"
    as two separate solution phases sharing the exact same objective
    function (obj_gh_opx is a thin wrapper around obj_gh_cpx below),
    letting MAGEMin's own multi-phase competition pick whichever is
    stable - reproducing real MELTS' behavior (which itself never
    actually calls isClino() in production either, since ISCLINO is
    always undefined) without needing any new "compare and choose"
    architecture.

    All 245 H/S/V Taylor coefficients below were extracted by compiling
    xMELTS' own #define macros (clinopyroxene.c lines 211-893) with
    clino=1 and calculationMode=MODE_xMELTS hardcoded, then printing every
    coefficient's numeric value directly - not manually re-derived, which
    eliminates transcription/arithmetic risk for a ~100-term formula (see
    scratchpad gen_cpx_constants.c). This also settled what the "p"-
    prefixed constants seen during initial inspection actually are:
    clino=TRUE collapses H1..H7/S1..S7/V1..V7 (the ortho-vertex offsets)
    to exactly zero and activates pH1..pH7 (Ghiorso's "sliding standard
    state for pigeonite", added in xMELTS v3.1-1) - fixed constants folded
    into the clino-branch coefficients, NOT a third structural state
    requiring its own isClino-style comparison.

    G(r,s) (the "G" expression below) is xMELTS' own raw Taylor+ideal-
    entropy sum - unlike cum/spn, it does NOT already vanish at every
    non-dependent pure endmember (only at diopside, r=0,s=0, where
    H0=S0=V0=0 identically when clino - verified analytically and
    numerically, not assumed). xMELTS' own gmixCpx explicitly subtracts a
    mechanical-mixture "ENDMEMBERS" reference (each endmember's own
    G(r=e_i,s*) via purePyx) before returning gmix; ported the same way:
    `ends[i]` solves the SAME 2D order-parameter problem at each of the 7
    pure-endmember r-vectors, and Gex = G(r,s*) - sum_i p_i*ends[i].
    By the envelope theorem, dGex[i] = dG/dp_i - ends[i] (ends[i] doesn't
    depend on the trial composition, only on T,P).

    p<->r basis: p[0]=di,p[1]=cen,p[2]=hed,p[3]=cats(Al-Tschermak),
    p[4]=buff(Fe3+-Tschermak),p[5]=ess,p[6]=jd - gh's direct mole-fraction
    convention. r0=p2, r1=p3+p6/2, r2=p4-p6/2, r3=p5+p6/2, r4=p6, r5=p1 -
    a genuine LINEAR SUBSTITUTION (not a relabeling like spn's r-basis
    was), extracted directly from xMELTS' own FR2..FR7/ENDMEMBERS macros
    (clinopyroxene.c lines 1924-1948: e.g. jadeite's own r-vector is
    r=(0,0.5,-0.5,0.5,1,0), not a unit vector). di (p[0]) is the
    dependent endmember.

    Order-parameter separability: s0 (Al3+/Fe3+ ordering on M1+tet) and
    s1 (Fe2+/Mg ordering on M1+M2) affect two COMPLETELY DISJOINT sets of
    site fractions (confirmed directly from xMELTS' own order() variable
    assignments) - so (unlike spn's 3-way-coupled s0/s1/s2) each 1D
    sub-problem's bounds only need the OTHER order parameter as a fixed
    linear offset. Solved via the same "bounds + bisection-on-sign,
    Gauss-Seidel cycle to self-consistency" pattern as spn
    (GH_cpx_bounds_s0/s1, GH_cpx_bisect1d, GH_cpx_solve_s_from), default
    single-start (gv.gh_multistart_order=0) matching real MELTS' own
    single Newton path; the starting guess is xMELTS' own ab-initio
    formula from order() (site-ratio distribution), not invented.

    Independently eps-clamping the 12 elementary site fractions does not
    guarantee every SUM/DIFFERENCE of them used as a log() argument in
    SIC/DGDR1/2/4/5 stays within (0,1) - GH_cpx_composites clamps those 7
    composite quantities too (mirrors xMELTS' own order() guard "if
    xmg2m1+xfe2m1(+xna1m2) >= 1.0, decrement xmg2m1/xfe2m1"). Verified
    (standalone finite-difference script, scratchpad test_cpx2.c): Gex is
    ~0 (to float noise) at all 7 pure endmembers and dGex matches central
    finite differences to ~1e-6 relative, across a T/P/composition sweep
    including near-endmember and near-Ti/Na-free compositions.
*/
    static const double H0 = 0, H027 = -16466.0265632, HS1 = -2677.12;
    static const double HS1S1 = -8917.5, HS1S2 = 7804.93, HS2 = -9811.9766408;
    static const double HS2S2 = 7460.595, HX2 = 7029.12, HX2S1 = -3819.645;
    static const double HX2S2 = -3117.08, HX2S2S2 = 169.975, HX2X2 = -7029.12;
    static const double HX2X2S2 = -339.95, HX2X2X7 = 339.95, HX2X3 = -8948.56;
    static const double HX2X3S2 = 5177.7, HX2X3X7 = -679.9, HX2X4 = -18247.52;
    static const double HX2X4S2 = 5177.7, HX2X4X7 = -679.9, HX2X5 = -20147.715;
    static const double HX2X5S2 = 5177.7, HX2X5X7 = -679.9, HX2X6 = -5424.3875;
    static const double HX2X6S2 = 2588.85, HX2X6X7 = -339.95, HX2X7 = 657.9867184;
    static const double HX2X7S2 = 10695.35, HX2X7X7 = -1012.005, HX3 = 16318;
    static const double HX3S1 = -19744.44, HX3S2 = 7537.4816796, HX3S2S2 = -5628.7875;
    static const double HX3X3 = -16318, HX3X3S2 = 1464.4, HX3X3X7 = -1464.4;
    static const double HX3X4 = -18357, HX3X4S2 = 2366.575, HX3X4X7 = -3491.025;
    static const double HX3X5 = -5465.44, HX3X5S2 = 3491.025, HX3X5X7 = -2366.575;
    static const double HX3X6 = 8160.4, HX3X6S2 = 1464.4, HX3X6X7 = -1464.4;
    static const double HX3X7 = 22496.3516796, HX3X7S2 = 1234.28, HX3X7X7 = -1923.3325;
    static const double HX4 = 18828, HX4S1 = 3347.06, HX4S2 = 7968.9216796;
    static const double HX4S2S2 = -5628.7875, HX4X4 = -18828, HX4X4S2 = 1464.4;
    static const double HX4X4X7 = -1464.4, HX4X5 = -15951.94, HX4X5S2 = 3491.025;
    static const double HX4X5X7 = -2366.575, HX4X6 = 23252.15, HX4X6S2 = 2026.625;
    static const double HX4X6X7 = -902.175, HX4X7 = 27147.9216796, HX4X7S2 = 1234.28;
    static const double HX4X7X7 = -1923.3325, HX5 = 23597.12, HX5S1 = -2065.5;
    static const double HX5S2 = -854.3466408, HX5S2S2 = -5347.675, HX5X5 = -9937;
    static const double HX5X5S2 = 2588.85, HX5X5X7 = -339.95, HX5X6 = -16226.25;
    static const double HX5X6S2 = 2588.85, HX5X6X7 = -339.95, HX5X7 = 17668.5133592;
    static const double HX5X7S2 = 3697.61, HX5X7X7 = 258.885, HX6 = -13220.68;
    static const double HX6S1 = -7322, HX6S2 = 7437.9366796, HX6S2S2 = -2673.8375;
    static const double HX6X6 = 8735.875, HX6X6S2 = 787.76875, HX6X6X7 = 55.56875;
    static const double HX6X7 = -651.4583204, HX6X7S2 = 1848.805, HX6X7X7 = 129.4425;
    static const double HX7 = 11610.1033592, HX7S1 = 4240.94, HX7S2 = 5460.6166408;
    static const double HX7S2S2 = -10695.35, HX7X7 = -17564.4583592, HX7X7S2 = -4246.76;
    static const double HX7X7X7 = 590.99, S0 = 0, S027 = -4.32073312;
    static const double SS1 = -0, SS1S1 = 0, SS1S2 = 0;
    /* SS2 is NOT calibration-agnostic like every other constant in this
       table: real clinopyroxene.c defines it as
       "#define SS2 ((calculationMode == MODE_xMELTS) ? 0.25*(S027) : 0.0)"
       - i.e. this one term (0.25*S027 = 0.25*(-4.32073312) = -1.08018328,
       matching the value below) is xMELTS-only; rMELTS/pMELTS both use
       0.0 instead. Every OTHER constant in this table is a genuine
       compile-time constant (fine as "static const" at this file scope);
       this one alone needs a runtime value, set by
       GH_SS_objective_init_function into the file-scope GH_cpx_SS2
       variable declared near GH_spn_multistart_flag above - so every
       formula below that used the extracted "SS2" constant now reads
       GH_cpx_SS2 directly instead (a file-scope "static const double SS2
       = GH_cpx_SS2;" here would itself violate C's compile-time-constant-
       initializer rule for statics, hence no local alias). Found while
       wiring multi-calibration support, 2026-07-14. */
    static const double SS2S2 = 0, SX2 = 0;
    static const double SX2S1 = 0, SX2S2 = 0, SX2X2 = 0;
    static const double SX2X3 = 0, SX2X4 = 0, SX2X5 = 0;
    static const double SX2X6 = 0, SX2X7 = -5.8088564, SX2X7X7 = 3.64848984;
    static const double SX3 = 0, SX3S1 = 0, SX3S2 = -0.54009164;
    static const double SX3X3 = 0, SX3X4 = 0, SX3X5 = 0;
    static const double SX3X6 = 0, SX3X7 = 0.35523636, SX3X7S2 = 0.91212246;
    static const double SX3X7X7 = 0.91212246, SX4 = 0, SX4S1 = 0;
    static const double SX4S2 = -0.54009164, SX4X4 = 0, SX4X5 = 0;
    static const double SX4X6 = 0, SX4X7 = 0.35523636, SX4X7S2 = 0.91212246;
    static const double SX4X7X7 = 0.91212246, SX5 = 0, SX5S1 = 0;
    static const double SX5S2 = -1.08018328, SX5X5 = 0, SX5X6 = 0;
    static const double SX5X7 = -0.18485528, SX5X7S2 = 1.82424492, SX5X7X7 = 1.82424492;
    static const double SX6 = 0, SX6S1 = 0, SX6S2 = -0.54009164;
    static const double SX6X6 = 0, SX6X7 = -0.09242764, SX6X7S2 = 0.91212246;
    static const double SX6X7X7 = 0.91212246, SX7 = -3.97551128, SX7S1 = 0;
    static const double SX7S2 = 2.9044282, SX7X7 = 5.7997562, SX7X7S2 = -1.82424492;
    static const double SX7X7X7 = -1.82424492, V0 = 0, V027 = -0.0385316;
    static const double VS1 = 0, VS1S1 = 0, VS1S2 = 0;
    static const double VS2 = -0.0064949, VS2S2 = -0.089956, VX2 = 0;
    static const double VX2S1 = 0, VX2S2 = 0.087864, VX2S2S2 = -0.005753;
    static const double VX2X2 = -0, VX2X2S2 = 0.011506, VX2X2X7 = -0.011506;
    static const double VX2X3 = -0, VX2X3S2 = -0.043932, VX2X3X7 = 0.023012;
    static const double VX2X4 = -0, VX2X4S2 = -0.043932, VX2X4X7 = 0.023012;
    static const double VX2X5 = -0, VX2X5S2 = -0.043932, VX2X5X7 = 0.023012;
    static const double VX2X6 = -0, VX2X6S2 = -0.021966, VX2X6X7 = 0.011506;
    static const double VX2X7 = -0.0862098, VX2X7S2 = -0.09937, VX2X7X7 = 0.08420296532;
    static const double VX3 = 0, VX3S1 = 0, VX3S2 = 0.00878155;
    static const double VX3S2S2 = 0.0509925, VX3X3 = 0, VX3X3S2 = -0.016736;
    static const double VX3X3X7 = 0.016736, VX3X4 = 0, VX3X4S2 = -0.030857;
    static const double VX3X4X7 = 0.036087, VX3X5 = 0, VX3X5S2 = -0.036087;
    static const double VX3X5X7 = 0.030857, VX3X6 = 0, VX3X6S2 = -0.016736;
    static const double VX3X6X7 = 0.016736, VX3X7 = -0.07179645, VX3X7S2 = -0.02928800867;
    static const double VX3X7X7 = 0.02013549133, VX4 = 0, VX4S1 = 0;
    static const double VX4S2 = 0.00878155, VX4S2S2 = 0.0509925, VX4X4 = 0;
    static const double VX4X4S2 = -0.016736, VX4X4X7 = 0.016736, VX4X5 = 0;
    static const double VX4X5S2 = -0.036087, VX4X5X7 = 0.030857, VX4X6 = 0;
    static const double VX4X6S2 = -0.019351, VX4X6X7 = 0.014121, VX4X7 = -0.07179645;
    static const double VX4X7S2 = -0.02928800867, VX4X7X7 = 0.02013549133, VX5 = 0;
    static const double VX5S1 = 0, VX5S2 = 0.0091951, VX5S2S2 = 0.049685;
    static const double VX5X5 = 0, VX5X5S2 = -0.021966, VX5X5X7 = 0.011506;
    static const double VX5X6 = 0, VX5X6S2 = -0.021966, VX5X6X7 = 0.011506;
    static const double VX5X7 = -0.0713829, VX5X7S2 = -0.03451801734, VX5X7X7 = 0.01621298266;
    static const double VX6 = 0, VX6S1 = 0, VX6S2 = 0.00459755;
    static const double VX6S2S2 = 0.0248425, VX6X6 = 0, VX6X6S2 = -0.00614525;
    static const double VX6X6X7 = 0.00222275, VX6X7 = -0.03569145, VX6X7S2 = -0.01725900867;
    static const double VX6X7X7 = 0.00810649133, VX7 = -0.0106789, VX7S1 = 0;
    static const double VX7S2 = 0.01904688266, VX7S2S2 = 0.09937, VX7X7 = 0.04624288266;
    static const double VX7X7S2 = -0.01255198266, VX7X7X7 = -0.04497798266;

static void GH_cpx_site_fracs(double r0,double r1,double r2,double r3,double r4,double r5,double s0,double s1,
    double *xal3m1,double *xfe2m1,double *xfe3m1,double *xmg2m1,double *xti4m1,
    double *xca2m2,double *xfe2m2,double *xmg2m2,double *xna1m2,
    double *xal3tet,double *xfe3tet,double *xsi4tet){
    double eps = 1.0e-9, v;
    v = (r1+r2)/2.0;                       *xti4m1  = (v<eps)?eps:((v>1.0-eps)?1.0-eps:v);
    v = 1.0-r4-r5;                         *xca2m2  = (v<eps)?eps:((v>1.0-eps)?1.0-eps:v);
    v = r4;                                *xna1m2  = (v<eps)?eps:((v>1.0-eps)?1.0-eps:v);
    v = (4.0-2.0*r1-2.0*r2-2.0*r3+r4)/4.0; *xsi4tet = (v<eps)?eps:((v>1.0-eps)?1.0-eps:v);
    v = (2.0*r3+r4-2.0*s0)/4.0;            *xal3m1  = (v<eps)?eps:((v>1.0-eps)?1.0-eps:v);
    v = (2.0*r0-r5-s1)/2.0;                *xfe2m1  = (v<eps)?eps:((v>1.0-eps)?1.0-eps:v);
    v = (2.0*r3+r4+2.0*s0)/4.0;            *xfe3m1  = (v<eps)?eps:((v>1.0-eps)?1.0-eps:v);
    v = 1.0-r0-r3-0.5*(r1+r2+r4-r5-s1);    *xmg2m1  = (v<eps)?eps:((v>1.0-eps)?1.0-eps:v);
    v = (r5+s1)/2.0;                       *xfe2m2  = (v<eps)?eps:((v>1.0-eps)?1.0-eps:v);
    v = (r5-s1)/2.0;                       *xmg2m2  = (v<eps)?eps:((v>1.0-eps)?1.0-eps:v);
    v = (4.0*r1+2.0*r3-r4+2.0*s0)/8.0;     *xal3tet = (v<eps)?eps:((v>1.0-eps)?1.0-eps:v);
    v = (4.0*r2+2.0*r3-r4-2.0*s0)/8.0;     *xfe3tet = (v<eps)?eps:((v>1.0-eps)?1.0-eps:v);
}

static void GH_cpx_composites(double xmg2m1,double xfe2m1,double xti4m1,double xna1m2,double xsi4tet,
    double *xm1_MgFeNoTi,double *xm1_notTi,double *xm1_MgFe,double *xtet_notSi,
    double *xrest_notNa,double *xm1_notMgFe,double *xm2_notNa){
    double eps = 1.0e-9, v;
    v = xmg2m1+xfe2m1-xti4m1;     *xm1_MgFeNoTi = (v<eps)?eps:((v>1.0-eps)?1.0-eps:v);
    v = 1.0-xti4m1;               *xm1_notTi    = (v<eps)?eps:((v>1.0-eps)?1.0-eps:v);
    v = xmg2m1+xfe2m1;            *xm1_MgFe     = (v<eps)?eps:((v>1.0-eps)?1.0-eps:v);
    v = 1.0-xsi4tet;              *xtet_notSi   = (v<eps)?eps:((v>1.0-eps)?1.0-eps:v);
    v = 1.0-xmg2m1-xfe2m1-xna1m2; *xrest_notNa  = (v<eps)?eps:((v>1.0-eps)?1.0-eps:v);
    v = 1.0-xmg2m1-xfe2m1;        *xm1_notMgFe  = (v<eps)?eps:((v>1.0-eps)?1.0-eps:v);
    v = 1.0-xna1m2;               *xm2_notNa    = (v<eps)?eps:((v>1.0-eps)?1.0-eps:v);
}

static void GH_cpx_dgds(double r0,double r1,double r2,double r3,double r4,double r5,double s0,double s1,
                         double T,double Rgas,double Pv,double *dgds0,double *dgds1){
    double xal3m1,xfe2m1,xfe3m1,xmg2m1,xti4m1,xca2m2,xfe2m2,xmg2m2,xna1m2,xal3tet,xfe3tet,xsi4tet;
    GH_cpx_site_fracs(r0,r1,r2,r3,r4,r5,s0,s1,&xal3m1,&xfe2m1,&xfe3m1,&xmg2m1,&xti4m1,&xca2m2,&xfe2m2,&xmg2m2,&xna1m2,&xal3tet,&xfe3tet,&xsi4tet);
    *dgds0 = 0.5*Rgas*T*(creal(clog(xfe3m1))-creal(clog(xal3m1))+creal(clog(xal3tet))-creal(clog(xfe3tet))) +                ((HS1)-T*(SS1)+Pv*(VS1)) +                ((HX2S1)-T*(SX2S1)+Pv*(VX2S1))*r0 +                ((HX3S1)-T*(SX3S1)+Pv*(VX3S1))*r1 +                ((HX4S1)-T*(SX4S1)+Pv*(VX4S1))*r2 +                ((HX5S1)-T*(SX5S1)+Pv*(VX5S1))*r3 +                ((HX6S1)-T*(SX6S1)+Pv*(VX6S1))*r4 +                ((HX7S1)-T*(SX7S1)+Pv*(VX7S1))*r5 +                ((HS1S1)-T*(SS1S1)+Pv*(VS1S1))*s0*2.0 +                ((HS1S2)-T*(SS1S2)+Pv*(VS1S2))*s1;
    *dgds1 = 0.5*Rgas*T*(creal(clog(xmg2m1))-creal(clog(xfe2m1))+creal(clog(xfe2m2))-creal(clog(xmg2m2))) +                ((HS2)-T*(GH_cpx_SS2)+Pv*(VS2)) +                ((HX2S2)-T*(SX2S2)+Pv*(VX2S2))*r0 +                ((HX3S2)-T*(SX3S2)+Pv*(VX3S2))*r1 +                ((HX4S2)-T*(SX4S2)+Pv*(VX4S2))*r2 +                ((HX5S2)-T*(SX5S2)+Pv*(VX5S2))*r3 +                ((HX6S2)-T*(SX6S2)+Pv*(VX6S2))*r4 +                ((HX7S2)-T*(SX7S2)+Pv*(VX7S2))*r5 +                ((HS1S2)-T*(SS1S2)+Pv*(VS1S2))*s0 +                ((HS2S2)-T*(SS2S2)+Pv*(VS2S2))*s1*2.0 +                ((HX2X2S2)+Pv*(VX2X2S2))*r0*r0 +                ((HX2X3S2)+Pv*(VX2X3S2))*r0*r1 +                ((HX2X4S2)+Pv*(VX2X4S2))*r0*r2 +                ((HX2X5S2)+Pv*(VX2X5S2))*r0*r3 +                ((HX2X6S2)+Pv*(VX2X6S2))*r0*r4 +                ((HX2X7S2)+Pv*(VX2X7S2))*r0*r5 +                ((HX2S2S2)+Pv*(VX2S2S2))*r0*s1*2.0 +                ((HX3X3S2)+Pv*(VX3X3S2))*r1*r1 +                ((HX3X4S2)+Pv*(VX3X4S2))*r1*r2 +                ((HX3X5S2)+Pv*(VX3X5S2))*r1*r3 +                ((HX3X6S2)+Pv*(VX3X6S2))*r1*r4 +                ((HX3X7S2)-T*(SX3X7S2)+Pv*(VX3X7S2))*r1*r5 +                ((HX3S2S2)+Pv*(VX3S2S2))*r1*s1*2.0 +                ((HX4X4S2)+Pv*(VX4X4S2))*r2*r2 +                ((HX4X5S2)+Pv*(VX4X5S2))*r2*r3 +                ((HX4X6S2)+Pv*(VX4X6S2))*r2*r4 +                ((HX4X7S2)-T*(SX4X7S2)+Pv*(VX4X7S2))*r2*r5 +                ((HX4S2S2)+Pv*(VX4S2S2))*r2*s1*2.0 +                ((HX5X5S2)+Pv*(VX5X5S2))*r3*r3 +                ((HX5X6S2)+Pv*(VX5X6S2))*r3*r4 +                ((HX5X7S2)-T*(SX5X7S2)+Pv*(VX5X7S2))*r3*r5 +                ((HX5S2S2)+Pv*(VX5S2S2))*r3*s1*2.0 +                ((HX6X6S2)+Pv*(VX6X6S2))*r4*r4 +                ((HX6X7S2)-T*(SX6X7S2)+Pv*(VX6X7S2))*r4*r5 +                ((HX6S2S2)+Pv*(VX6S2S2))*r4*s1*2.0 +                ((HX7X7S2)-T*(SX7X7S2)+Pv*(VX7X7S2))*r5*r5 +                ((HX7S2S2)+Pv*(VX7S2S2))*r5*s1*2.0;
}

static void GH_cpx_bounds_s0(double r1,double r2,double r3,double r4,double eps,double *lo,double *hi){
    *lo = -1.0e6; *hi = 1.0e6;
    double val, slope;
    val = (2.0*r3+r4)/4.0;      slope = -0.5;  GH_narrow_bound(val,slope,eps,lo,hi); /* xal3m1  */
    val = (2.0*r3+r4)/4.0;      slope =  0.5;  GH_narrow_bound(val,slope,eps,lo,hi); /* xfe3m1  */
    val = (4.0*r1+2.0*r3-r4)/8.0; slope = 0.25; GH_narrow_bound(val,slope,eps,lo,hi); /* xal3tet */
    val = (4.0*r2+2.0*r3-r4)/8.0; slope = -0.25; GH_narrow_bound(val,slope,eps,lo,hi); /* xfe3tet */
}
static void GH_cpx_bounds_s1(double r0,double r1,double r2,double r3,double r4,double r5,double eps,double *lo,double *hi){
    *lo = -1.0e6; *hi = 1.0e6;
    double val, slope;
    val = (2.0*r0-r5)/2.0;                        slope = -0.5; GH_narrow_bound(val,slope,eps,lo,hi); /* xfe2m1 */
    val = 1.0-r0-r3-0.5*(r1+r2+r4-r5);            slope =  0.5; GH_narrow_bound(val,slope,eps,lo,hi); /* xmg2m1 */
    val = r5/2.0;                                  slope =  0.5; GH_narrow_bound(val,slope,eps,lo,hi); /* xfe2m2 */
    val = r5/2.0;                                  slope = -0.5; GH_narrow_bound(val,slope,eps,lo,hi); /* xmg2m2 */
}
/* Hessian d2G/ds_i ds_j, from xMELTS' own D2GDS0S0/D2GDS0S1/D2GDS1S1
   macros (clinopyroxene.c lines 2702-2720) - needed for the real joint
   Newton solve below (not previously ported; an earlier version of this
   port used a per-variable Gauss-Seidel/bisection cycle needing only
   first derivatives). */
static void GH_cpx_d2gds(double r0,double r1,double r2,double r3,double r4,double r5,double s0,double s1,
                          double T,double Rgas,double Pv,double *d00,double *d01,double *d11){
    double xal3m1,xfe2m1,xfe3m1,xmg2m1,xti4m1,xca2m2,xfe2m2,xmg2m2,xna1m2,xal3tet,xfe3tet,xsi4tet;
    GH_cpx_site_fracs(r0,r1,r2,r3,r4,r5,s0,s1,&xal3m1,&xfe2m1,&xfe3m1,&xmg2m1,&xti4m1,&xca2m2,&xfe2m2,&xmg2m2,&xna1m2,&xal3tet,&xfe3tet,&xsi4tet);
    *d00 = Rgas*T*(0.25/xfe3m1+0.25/xal3m1+0.125/xal3tet+0.125/xfe3tet) + 2.0*((HS1S1)-T*(SS1S1)+Pv*(VS1S1));
    *d01 = ((HS1S2)-T*(SS1S2)+Pv*(VS1S2));
    *d11 = Rgas*T*(0.25/xmg2m1+0.25/xfe2m1+0.25/xfe2m2+0.25/xmg2m2) + 2.0*((HS2S2)-T*(SS2S2)+Pv*(VS2S2))
         + 2.0*r0*((HX2S2S2)+Pv*(VX2S2S2)) + 2.0*r1*((HX3S2S2)+Pv*(VX3S2S2)) + 2.0*r2*((HX4S2S2)+Pv*(VX4S2S2))
         + 2.0*r3*((HX5S2S2)+Pv*(VX5S2S2)) + 2.0*r4*((HX6S2S2)+Pv*(VX6S2S2)) + 2.0*r5*((HX7S2S2)+Pv*(VX7S2S2));
}
/* Joint damped Newton solve for s0,s1 - ported directly from xMELTS' own
   order() (clinopyroxene.c lines 3601-3675), which uses a genuine joint
   2x2 Newton step (via gaussj, i.e. the real Hessian above) followed by
   hard-clamping s0/s1 into their r-dependent feasible bounds - not the
   per-variable Gauss-Seidel/bisection cycle an earlier version of this
   port used (same story, and same fix, as spinel's solver - see
   obj_gh_spn's header comment). Reuses GH_cpx_bounds_s0/s1 (already
   verified equivalent to xMELTS' explicit MAX/MIN bound formulas) for
   the clamp step instead of re-deriving xMELTS' explicit expressions. */
static void GH_cpx_solve_s_from(double r0,double r1,double r2,double r3,double r4,double r5,double T,double Rgas,double Pv,
                                 double s0init,double s1init,double *s0o,double *s1o){
    const int maxIter = 1000; /* xMELTS' own maxIter (clinopyroxene.c) */
    double eps_s = 1.0e-8;
    double s0=s0init, s1=s1init;
    double s0Old=2.0, s1Old=2.0; /* force >=1 iteration, matching xMELTS' sOld[i]=2.0 init */
    int iter = 0;
    while ( (fabs(s0-s0Old) > 10.0*2.220446049250313e-16 ||
             fabs(s1-s1Old) > 10.0*2.220446049250313e-16) && iter < maxIter ){
        double dgds0, dgds1;
        GH_cpx_dgds(r0,r1,r2,r3,r4,r5,s0,s1,T,Rgas,Pv,&dgds0,&dgds1);
        double d00, d01, d11;
        GH_cpx_d2gds(r0,r1,r2,r3,r4,r5,s0,s1,T,Rgas,Pv,&d00,&d01,&d11);

        s0Old = s0; s1Old = s1;

        double det = d00*d11 - d01*d01;
        double deltaS0 = -( d11*dgds0 - d01*dgds1)/det;
        double deltaS1 = -(-d01*dgds0 + d00*dgds1)/det;

        double sN0 = s0Old + deltaS0, sN1 = s1Old + deltaS1;
        double lo, hi;
        GH_cpx_bounds_s0(r1,r2,r3,r4,eps_s,&lo,&hi);
        if (sN0 < lo){ sN0 = lo; } if (sN0 > hi){ sN0 = hi; }
        GH_cpx_bounds_s1(r0,r1,r2,r3,r4,r5,eps_s,&lo,&hi);
        if (sN1 < lo){ sN1 = lo; } if (sN1 > hi){ sN1 = hi; }
        s0 = sN0; s1 = sN1;
        iter++;
    }
    *s0o = s0; *s1o = s1;
}
/* xMELTS' own ab-initio starting guess (order()'s "distribute cations
   proportional to available site" formula), used both for the trial
   composition and every fixed endmember corner below. */
static void GH_cpx_guess_s(double r0,double r1,double r2,double r3,double r4,double r5,double *s0g,double *s1g){
    double totAl=r1+r3, totCa=1.0-r4-r5, totFe2=r0, totFe3=r2+r3;
    double totMg=1.0-r0-0.5*(r1+r2)-r3-0.5*r4+r5, totNa=r4, totTi=0.5*(r1+r2);
    double totM2=totCa+totNa, totM1=totTi+(totMg+totFe2-(1.0-totM2));
    *s0g = (totFe3+totAl != 0.0) ? (1.0-totM1)*(totFe3-totAl)/(totFe3+totAl) : 0.0;
    *s1g = (totFe2+totMg != 0.0) ? (1.0-totM2)*(totFe2-totMg)/(totFe2+totMg) : 0.0;
}
static double GH_cpx_G_at(double r0,double r1,double r2,double r3,double r4,double r5,double s0,double s1,
                           double T,double Rgas,double Pv){
    double xal3m1,xfe2m1,xfe3m1,xmg2m1,xti4m1,xca2m2,xfe2m2,xmg2m2,xna1m2,xal3tet,xfe3tet,xsi4tet;
    GH_cpx_site_fracs(r0,r1,r2,r3,r4,r5,s0,s1,&xal3m1,&xfe2m1,&xfe3m1,&xmg2m1,&xti4m1,&xca2m2,&xfe2m2,&xmg2m2,&xna1m2,&xal3tet,&xfe3tet,&xsi4tet);
    double xm1_MgFeNoTi,xm1_notTi,xm1_MgFe,xtet_notSi,xrest_notNa,xm1_notMgFe,xm2_notNa;
    GH_cpx_composites(xmg2m1,xfe2m1,xti4m1,xna1m2,xsi4tet,&xm1_MgFeNoTi,&xm1_notTi,&xm1_MgFe,&xtet_notSi,&xrest_notNa,&xm1_notMgFe,&xm2_notNa);
    double SIC = -Rgas*(xmg2m1*creal(clog(xmg2m1))                + xfe2m1*creal(clog(xfe2m1))                + xal3m1*creal(clog(xal3m1))                + xfe3m1*creal(clog(xfe3m1))                + xti4m1*creal(clog(xti4m1))                + (xm1_MgFeNoTi)*creal(clog(xm1_MgFeNoTi))                - (xm1_notTi)*creal(clog(xm1_notTi))                - (xm1_MgFe)*creal(clog(xm1_MgFe))                - 2.0*(xtet_notSi)*creal(clog(xtet_notSi))                + 2.0*xal3tet*creal(clog(xal3tet))                + 2.0*xfe3tet*creal(clog(xfe3tet))                + xca2m2*creal(clog(xca2m2))                + xna1m2*creal(clog(xna1m2))                + xmg2m2*creal(clog(xmg2m2))                + xfe2m2*creal(clog(xfe2m2))                + (xrest_notNa)*creal(clog(xrest_notNa))                - (xm1_notMgFe)*creal(clog(xm1_notMgFe))                - (xm2_notNa)*creal(clog(xm2_notNa)) );
    return -T*(SIC) + (H0)-T*(S0)+Pv*(V0) +              ((HX2)-T*(SX2)+Pv*(VX2))*r0 +              ((HX3)-T*(SX3)+Pv*(VX3))*r1 +              ((HX4)-T*(SX4)+Pv*(VX4))*r2 +              ((HX5)-T*(SX5)+Pv*(VX5))*r3 +              ((HX6)-T*(SX6)+Pv*(VX6))*r4 +              ((HX7)-T*(SX7)+Pv*(VX7))*r5 +              ((HS1)-T*(SS1)+Pv*(VS1))*s0 +              ((HS2)-T*(GH_cpx_SS2)+Pv*(VS2))*s1 +              ((HX2X2)-T*(SX2X2)+Pv*(VX2X2))*r0*r0 +              ((HX2X3)-T*(SX2X3)+Pv*(VX2X3))*r0*r1 +              ((HX2X4)-T*(SX2X4)+Pv*(VX2X4))*r0*r2 +              ((HX2X5)-T*(SX2X5)+Pv*(VX2X5))*r0*r3 +              ((HX2X6)-T*(SX2X6)+Pv*(VX2X6))*r0*r4 +              ((HX2X7)-T*(SX2X7)+Pv*(VX2X7))*r0*r5 +              ((HX2S1)-T*(SX2S1)+Pv*(VX2S1))*r0*s0 +              ((HX2S2)-T*(SX2S2)+Pv*(VX2S2))*r0*s1 +              ((HX3X3)-T*(SX3X3)+Pv*(VX3X3))*r1*r1 +              ((HX3X4)-T*(SX3X4)+Pv*(VX3X4))*r1*r2 +              ((HX3X5)-T*(SX3X5)+Pv*(VX3X5))*r1*r3 +              ((HX3X6)-T*(SX3X6)+Pv*(VX3X6))*r1*r4 +              ((HX3X7)-T*(SX3X7)+Pv*(VX3X7))*r1*r5 +              ((HX3S1)-T*(SX3S1)+Pv*(VX3S1))*r1*s0 +              ((HX3S2)-T*(SX3S2)+Pv*(VX3S2))*r1*s1 +              ((HX4X4)-T*(SX4X4)+Pv*(VX4X4))*r2*r2 +              ((HX4X5)-T*(SX4X5)+Pv*(VX4X5))*r2*r3 +              ((HX4X6)-T*(SX4X6)+Pv*(VX4X6))*r2*r4 +              ((HX4X7)-T*(SX4X7)+Pv*(VX4X7))*r2*r5 +              ((HX4S1)-T*(SX4S1)+Pv*(VX4S1))*r2*s0 +              ((HX4S2)-T*(SX4S2)+Pv*(VX4S2))*r2*s1 +              ((HX5X5)-T*(SX5X5)+Pv*(VX5X5))*r3*r3 +              ((HX5X6)-T*(SX5X6)+Pv*(VX5X6))*r3*r4 +              ((HX5X7)-T*(SX5X7)+Pv*(VX5X7))*r3*r5 +              ((HX5S1)-T*(SX5S1)+Pv*(VX5S1))*r3*s0 +              ((HX5S2)-T*(SX5S2)+Pv*(VX5S2))*r3*s1 +              ((HX6X6)-T*(SX6X6)+Pv*(VX6X6))*r4*r4 +              ((HX6X7)-T*(SX6X7)+Pv*(VX6X7))*r4*r5 +              ((HX6S1)-T*(SX6S1)+Pv*(VX6S1))*r4*s0 +              ((HX6S2)-T*(SX6S2)+Pv*(VX6S2))*r4*s1 +              ((HX7X7)-T*(SX7X7)+Pv*(VX7X7))*r5*r5 +              ((HX7S1)-T*(SX7S1)+Pv*(VX7S1))*r5*s0 +              ((HX7S2)-T*(SX7S2)+Pv*(VX7S2))*r5*s1 +              ((HS1S1)-T*(SS1S1)+Pv*(VS1S1))*s0*s0 +              ((HS1S2)-T*(SS1S2)+Pv*(VS1S2))*s0*s1 +              ((HS2S2)-T*(SS2S2)+Pv*(VS2S2))*s1*s1 +              ((HX2X2X7)+Pv*(VX2X2X7))*r0*r0*r5 +              ((HX2X2S2)+Pv*(VX2X2S2))*r0*r0*s1 +              ((HX2X3X7)+Pv*(VX2X3X7))*r0*r1*r5 +              ((HX2X3S2)+Pv*(VX2X3S2))*r0*r1*s1 +              ((HX2X4X7)+Pv*(VX2X4X7))*r0*r2*r5 +              ((HX2X4S2)+Pv*(VX2X4S2))*r0*r2*s1 +              ((HX2X5X7)+Pv*(VX2X5X7))*r0*r3*r5 +              ((HX2X5S2)+Pv*(VX2X5S2))*r0*r3*s1 +              ((HX2X6X7)+Pv*(VX2X6X7))*r0*r4*r5 +              ((HX2X6S2)+Pv*(VX2X6S2))*r0*r4*s1 +              ((HX2X7X7)-T*(SX2X7X7)+Pv*(VX2X7X7))*r0*r5*r5 +              ((HX2X7S2)+Pv*(VX2X7S2))*r0*r5*s1 +              ((HX2S2S2)+Pv*(VX2S2S2))*r0*s1*s1 +              ((HX3X3X7)+Pv*(VX3X3X7))*r1*r1*r5 +              ((HX3X3S2)+Pv*(VX3X3S2))*r1*r1*s1 +              ((HX3X4X7)+Pv*(VX3X4X7))*r1*r2*r5 +              ((HX3X4S2)+Pv*(VX3X4S2))*r1*r2*s1 +              ((HX3X5X7)+Pv*(VX3X5X7))*r1*r3*r5 +              ((HX3X5S2)+Pv*(VX3X5S2))*r1*r3*s1 +              ((HX3X6X7)+Pv*(VX3X6X7))*r1*r4*r5 +              ((HX3X6S2)+Pv*(VX3X6S2))*r1*r4*s1 +              ((HX3X7X7)-T*(SX3X7X7)+Pv*(VX3X7X7))*r1*r5*r5 +              ((HX3X7S2)-T*(SX3X7S2)+Pv*(VX3X7S2))*r1*r5*s1 +              ((HX3S2S2)+Pv*(VX3S2S2))*r1*s1*s1 +              ((HX4X4X7)+Pv*(VX4X4X7))*r2*r2*r5 +              ((HX4X4S2)+Pv*(VX4X4S2))*r2*r2*s1 +              ((HX4X5X7)+Pv*(VX4X5X7))*r2*r3*r5 +              ((HX4X5S2)+Pv*(VX4X5S2))*r2*r3*s1 +              ((HX4X6X7)+Pv*(VX4X6X7))*r2*r4*r5 +              ((HX4X6S2)+Pv*(VX4X6S2))*r2*r4*s1 +              ((HX4X7X7)-T*(SX4X7X7)+Pv*(VX4X7X7))*r2*r5*r5 +              ((HX4X7S2)-T*(SX4X7S2)+Pv*(VX4X7S2))*r2*r5*s1 +              ((HX4S2S2)+Pv*(VX4S2S2))*r2*s1*s1 +              ((HX5X5X7)+Pv*(VX5X5X7))*r3*r3*r5 +              ((HX5X5S2)+Pv*(VX5X5S2))*r3*r3*s1 +              ((HX5X6X7)+Pv*(VX5X6X7))*r3*r4*r5 +              ((HX5X6S2)+Pv*(VX5X6S2))*r3*r4*s1 +              ((HX5X7X7)-T*(SX5X7X7)+Pv*(VX5X7X7))*r3*r5*r5 +              ((HX5X7S2)-T*(SX5X7S2)+Pv*(VX5X7S2))*r3*r5*s1 +              ((HX5S2S2)+Pv*(VX5S2S2))*r3*s1*s1 +              ((HX6X6X7)+Pv*(VX6X6X7))*r4*r4*r5 +              ((HX6X6S2)+Pv*(VX6X6S2))*r4*r4*s1 +              ((HX6X7X7)-T*(SX6X7X7)+Pv*(VX6X7X7))*r4*r5*r5 +              ((HX6X7S2)-T*(SX6X7S2)+Pv*(VX6X7S2))*r4*r5*s1 +              ((HX6S2S2)+Pv*(VX6S2S2))*r4*s1*s1 +              ((HX7X7X7)-T*(SX7X7X7)+Pv*(VX7X7X7))*r5*r5*r5 +              ((HX7X7S2)-T*(SX7X7S2)+Pv*(VX7X7S2))*r5*r5*s1 +              ((HX7S2S2)+Pv*(VX7S2S2))*r5*s1*s1;
}
/* Solve the embedded order-parameter problem at a FIXED r-vector (used
   only for the trial composition itself - see GH_cpx_pure_ES_G below for
   why the 7 pure-endmember reference values do NOT need this) and return
   G(r,s*) there. */
static double GH_cpx_solve_and_G(double r0,double r1,double r2,double r3,double r4,double r5,
                                  double T,double Rgas,double Pv,double *s0o,double *s1o){
    double s0g,s1g; GH_cpx_guess_s(r0,r1,r2,r3,r4,r5,&s0g,&s1g);
    double cs0,cs1; GH_cpx_solve_s_from(r0,r1,r2,r3,r4,r5,T,Rgas,Pv,s0g,s1g,&cs0,&cs1);
    *s0o = cs0; *s1o = cs1;
    return GH_cpx_G_at(r0,r1,r2,r3,r4,r5,cs0,cs1,T,Rgas,Pv);
}

/**
    Pure-endmember reference values ("ends[]"), matching xMELTS' own
    purePyx (NOT a re-run of the solid-solution's order() problem at each
    corner - a brute-force earlier version of this port did that, and got
    numerically correct but needlessly expensive results). Direct
    inspection of purePyx shows DI_G/EN_G/HD_G/CA_G/CF_G/JD_G are CLOSED
    FORM (s fixed at implicit exact values baked into the macro, no
    solve), and only ES_G needs an iterative solve - of a *different*,
    single-variable Landau order parameter specific to pure essenite
    (xMELTS' own separate `pureOrder` function), not the solid-solution's
    s0/s1.

    Better still: DI_G=EN_G=HD_G=CA_G=CF_G=JD_G are all EXACTLY 0, for
    every T,P - not approximately, verified by direct numeric summation
    of the extracted Taylor coefficients (H/S/V parts each individually
    cancel to 0.0, JD to -3.6e-12 = float noise). This isn't a
    coincidence: clino=TRUE zeroes every vertex H_i/S_i/V_i (i=1..7), and
    these six closed-form macros are literally that vertex's own
    H_i-T*S_i+(P-1)*V_i. Only essenite has a genuine nonzero reference
    value, from its own intrinsic Al3+/Fe3+ ordering entropy (the log(1+-s)
    term in ES_G/ES_S, which survives even though essenite's own H5/S5/V5
    are also 0). So Emix = p[5]*ends[5] - only essenite's own mole
    fraction matters.
*/
static double GH_cpx_pure_ES_dgds(double s,double T,double Rgas,double Pv){
    return Rgas*T*(creal(clog(1.0+s))-creal(clog(1.0-s))) + (HS1) + (HX5S1) + (HS1S1)*s*2.0
         - T*((SS1)+(SX5S1)+(SS1S1)*s*2.0) + Pv*((VS1)+(VX5S1)+(VS1S1)*s*2.0);
}
static double GH_cpx_pure_ES_G_at_s(double s,double T,double Rgas,double Pv){
    return Rgas*T*((1.0-s)*creal(clog(1.0-s)) + (1.0+s)*creal(clog(1.0+s)) - 2.0*creal(clog(2.0)))
         + (H0)+(HX5)+(HS1)*s+(HX5X5)+(HX5S1)*s+(HS1S1)*s*s
         - T*((S0)+(SX5)+(SS1)*s+(SX5X5)+(SX5S1)*s+(SS1S1)*s*s)
         + Pv*((V0)+(VX5)+(VS1)*s+(VX5X5)+(VS1S1)*s*s+(VX5S1)*s);
}
static double GH_cpx_pure_ES_G(double T,double Rgas,double Pv){
    double eps = 1.0e-8, lo = -1.0+eps, hi = 1.0-eps;
    double fa = GH_cpx_pure_ES_dgds(lo,T,Rgas,Pv), fb = GH_cpx_pure_ES_dgds(hi,T,Rgas,Pv);
    double a = lo, b = hi;
    if (fa > 0.0 && fb > 0.0){ a = b = lo; }
    else if (fa < 0.0 && fb < 0.0){ a = b = hi; }
    else {
        for (int it = 0; it < 80; it++){
            double mid = 0.5*(a+b), fm = GH_cpx_pure_ES_dgds(mid,T,Rgas,Pv);
            if ((fa > 0.0) == (fm > 0.0)){ a = mid; fa = fm; } else { b = mid; }
            if (b-a < 1.0e-14){ break; }
        }
    }
    return GH_cpx_pure_ES_G_at_s(0.5*(a+b),T,Rgas,Pv);
}

double obj_gh_cpx(unsigned n, const double *x, double *grad, void *SS_ref_db){
    SS_ref *d = (SS_ref *) SS_ref_db;

    double T    = d->T;
    double Rgas = d->R*1000.0;
    double Pv   = d->P - 1.0;
    double *p   = d->p;
    double *gb  = d->gb_lvl;
    double *mu_Gex = d->mu_Gex;

    for (int i = 0; i < 7; i++){ p[i] = x[i]; }
    double p1=p[1], p2=p[2], p3=p[3], p4=p[4], p5=p[5], p6=p[6];

    /* r-basis (xMELTS' own convention - see this function's header comment) */
    double r0=p2, r1=p3+0.5*p6, r2=p4-0.5*p6, r3=p5+0.5*p6, r4=p6, r5=p1;

    double s0, s1;
    double Graw = GH_cpx_solve_and_G(r0,r1,r2,r3,r4,r5,T,Rgas,Pv,&s0,&s1);

    double xal3m1,xfe2m1,xfe3m1,xmg2m1,xti4m1,xca2m2,xfe2m2,xmg2m2,xna1m2,xal3tet,xfe3tet,xsi4tet;
    GH_cpx_site_fracs(r0,r1,r2,r3,r4,r5,s0,s1,&xal3m1,&xfe2m1,&xfe3m1,&xmg2m1,&xti4m1,&xca2m2,&xfe2m2,&xmg2m2,&xna1m2,&xal3tet,&xfe3tet,&xsi4tet);
    double xm1_MgFeNoTi,xm1_notTi,xm1_MgFe,xtet_notSi,xrest_notNa,xm1_notMgFe,xm2_notNa;
    GH_cpx_composites(xmg2m1,xfe2m1,xti4m1,xna1m2,xsi4tet,&xm1_MgFeNoTi,&xm1_notTi,&xm1_MgFe,&xtet_notSi,&xrest_notNa,&xm1_notMgFe,&xm2_notNa);

    double dgdr0 = Rgas*T*(creal(clog(xfe2m1)) - creal(clog(xmg2m1))) +                ((HX2)-T*(SX2)+Pv*(VX2)) +                ((HX2X2)-T*(SX2X2)+Pv*(VX2X2))*r0*2.0 +                ((HX2X3)-T*(SX2X3)+Pv*(VX2X3))*r1 +                ((HX2X4)-T*(SX2X4)+Pv*(VX2X4))*r2 +                ((HX2X5)-T*(SX2X5)+Pv*(VX2X5))*r3 +                ((HX2X6)-T*(SX2X6)+Pv*(VX2X6))*r4 +                ((HX2X7)-T*(SX2X7)+Pv*(VX2X7))*r5 +                ((HX2S1)-T*(SX2S1)+Pv*(VX2S1))*s0 +                ((HX2S2)-T*(SX2S2)+Pv*(VX2S2))*s1 +                ((HX2X2X7)+Pv*(VX2X2X7))*r0*r5*2.0 +                ((HX2X2S2)+Pv*(VX2X2S2))*r0*s1*2.0 +                ((HX2X3X7)+Pv*(VX2X3X7))*r1*r5 +                ((HX2X3S2)+Pv*(VX2X3S2))*r1*s1 +                ((HX2X4X7)+Pv*(VX2X4X7))*r2*r5 +                ((HX2X4S2)+Pv*(VX2X4S2))*r2*s1 +                ((HX2X5X7)+Pv*(VX2X5X7))*r3*r5 +                ((HX2X5S2)+Pv*(VX2X5S2))*r3*s1 +                ((HX2X6X7)+Pv*(VX2X6X7))*r4*r5 +                ((HX2X6S2)+Pv*(VX2X6S2))*r4*s1 +                ((HX2X7X7)-T*(SX2X7X7)+Pv*(VX2X7X7))*r5*r5 +                ((HX2X7S2)+Pv*(VX2X7S2))*r5*s1 +                ((HX2S2S2)+Pv*(VX2S2S2))*s1*s1;
    double dgdr1 = Rgas*T*(0.5*creal(clog(xti4m1)) - 0.5*creal(clog(xmg2m1)) + 0.5*creal(clog(xm1_notTi))                  - creal(clog(xm1_MgFeNoTi)) + 0.5*creal(clog(xm1_MgFe))                  + 0.5*creal(clog(xrest_notNa)) - creal(clog(xtet_notSi))                  - 0.5*creal(clog(xm1_notMgFe)) + creal(clog(xal3tet))) +                ((HX3)-T*(SX3)+Pv*(VX3)) +                ((HX2X3)-T*(SX2X3)+Pv*(VX2X3))*r0 +                ((HX3X3)-T*(SX3X3)+Pv*(VX3X3))*r1*2.0 +                ((HX3X4)-T*(SX3X4)+Pv*(VX3X4))*r2 +                ((HX3X5)-T*(SX3X5)+Pv*(VX3X5))*r3 +                ((HX3X6)-T*(SX3X6)+Pv*(VX3X6))*r4 +                ((HX3X7)-T*(SX3X7)+Pv*(VX3X7))*r5 +                ((HX3S1)-T*(SX3S1)+Pv*(VX3S1))*s0 +                ((HX3S2)-T*(SX3S2)+Pv*(VX3S2))*s1 +                ((HX2X3X7)+Pv*(VX2X3X7))*r0*r5 +                ((HX2X3S2)+Pv*(VX2X3S2))*r0*s1 +                ((HX3X3X7)+Pv*(VX3X3X7))*r1*r5*2.0 +                ((HX3X3S2)+Pv*(VX3X3S2))*r1*s1*2.0 +                ((HX3X4X7)+Pv*(VX3X4X7))*r2*r5 +                ((HX3X4S2)+Pv*(VX3X4S2))*r2*s1 +                ((HX3X5X7)+Pv*(VX3X5X7))*r3*r5 +                ((HX3X5S2)+Pv*(VX3X5S2))*r3*s1 +                ((HX3X6X7)+Pv*(VX3X6X7))*r4*r5 +                ((HX3X6S2)+Pv*(VX3X6S2))*r4*s1 +                ((HX3X7X7)-T*(SX3X7X7)+Pv*(VX3X7X7))*r5*r5 +                ((HX3X7S2)-T*(SX3X7S2)+Pv*(VX3X7S2))*r5*s1 +                ((HX3S2S2)+Pv*(VX3S2S2))*s1*s1;
    double dgdr2 = Rgas*T*(0.5*creal(clog(xti4m1)) - 0.5*creal(clog(xmg2m1)) + 0.5*creal(clog(xm1_notTi))                  - creal(clog(xm1_MgFeNoTi)) + 0.5*creal(clog(xm1_MgFe))                  + 0.5*creal(clog(xrest_notNa)) - creal(clog(xtet_notSi))                  - 0.5*creal(clog(xm1_notMgFe)) + creal(clog(xfe3tet))) +                ((HX4)-T*(SX4)+Pv*(VX4)) +                ((HX2X4)-T*(SX2X4)+Pv*(VX2X4))*r0 +                ((HX3X4)-T*(SX3X4)+Pv*(VX3X4))*r1 +                ((HX4X4)-T*(SX4X4)+Pv*(VX4X4))*r2*2.0 +                ((HX4X5)-T*(SX4X5)+Pv*(VX4X5))*r3 +                ((HX4X6)-T*(SX4X6)+Pv*(VX4X6))*r4 +                ((HX4X7)-T*(SX4X7)+Pv*(VX4X7))*r5 +                ((HX4S1)-T*(SX4S1)+Pv*(VX4S1))*s0 +                ((HX4S2)-T*(SX4S2)+Pv*(VX4S2))*s1 +                ((HX2X4X7)+Pv*(VX2X4X7))*r0*r5 +                ((HX2X4S2)+Pv*(VX2X4S2))*r0*s1 +                ((HX3X4X7)+Pv*(VX3X4X7))*r1*r5 +                ((HX3X4S2)+Pv*(VX3X4S2))*r1*s1 +                ((HX4X4X7)+Pv*(VX4X4X7))*r2*r5*2.0 +                ((HX4X4S2)+Pv*(VX4X4S2))*r2*s1*2.0 +                ((HX4X5X7)+Pv*(VX4X5X7))*r3*r5 +                ((HX4X5S2)+Pv*(VX4X5S2))*r3*s1 +                ((HX4X6X7)+Pv*(VX4X6X7))*r4*r5 +                ((HX4X6S2)+Pv*(VX4X6S2))*r4*s1 +                ((HX4X7X7)-T*(SX4X7X7)+Pv*(VX4X7X7))*r5*r5 +                ((HX4X7S2)-T*(SX4X7S2)+Pv*(VX4X7S2))*r5*s1 +                ((HX4S2S2)+Pv*(VX4S2S2))*s1*s1;
    double dgdr3 = Rgas*T*(0.5*creal(clog(xal3m1)) + 0.5*creal(clog(xfe3m1)) - creal(clog(xmg2m1))                  - creal(clog(xm1_MgFeNoTi)) + creal(clog(xm1_MgFe))                  - creal(clog(xtet_notSi)) + 0.5*creal(clog(xal3tet)) + 0.5*creal(clog(xfe3tet))                  + creal(clog(xrest_notNa)) - creal(clog(xm1_notMgFe))) +                ((HX5)-T*(SX5)+Pv*(VX5)) +                ((HX2X5)-T*(SX2X5)+Pv*(VX2X5))*r0 +                ((HX3X5)-T*(SX3X5)+Pv*(VX3X5))*r1 +                ((HX4X5)-T*(SX4X5)+Pv*(VX4X5))*r2 +                ((HX5X5)-T*(SX5X5)+Pv*(VX5X5))*r3*2.0 +                ((HX5X6)-T*(SX5X6)+Pv*(VX5X6))*r4 +                ((HX5X7)-T*(SX5X7)+Pv*(VX5X7))*r5 +                ((HX5S1)-T*(SX5S1)+Pv*(VX5S1))*s0 +                ((HX5S2)-T*(SX5S2)+Pv*(VX5S2))*s1 +                ((HX2X5X7)+Pv*(VX2X5X7))*r0*r5 +                ((HX2X5S2)+Pv*(VX2X5S2))*r0*s1 +                ((HX3X5X7)+Pv*(VX3X5X7))*r1*r5 +                ((HX3X5S2)+Pv*(VX3X5S2))*r1*s1 +                ((HX4X5X7)+Pv*(VX4X5X7))*r2*r5 +                ((HX4X5S2)+Pv*(VX4X5S2))*r2*s1 +                ((HX5X5X7)+Pv*(VX5X5X7))*r3*r5*2.0 +                ((HX5X5S2)+Pv*(VX5X5S2))*r3*s1*2.0 +                ((HX5X6X7)+Pv*(VX5X6X7))*r4*r5 +                ((HX5X6S2)+Pv*(VX5X6S2))*r4*s1 +                ((HX5X7X7)-T*(SX5X7X7)+Pv*(VX5X7X7))*r5*r5 +                ((HX5X7S2)-T*(SX5X7S2)+Pv*(VX5X7S2))*r5*s1 +                ((HX5S2S2)+Pv*(VX5S2S2))*s1*s1;
    double dgdr4 = Rgas*T*(0.25*creal(clog(xal3m1)) + 0.25*creal(clog(xfe3m1)) - 0.5*creal(clog(xmg2m1))                  - 0.5*creal(clog(xm1_MgFeNoTi)) + 0.5*creal(clog(xm1_MgFe))                  + 0.5*creal(clog(xtet_notSi)) - 0.25*creal(clog(xal3tet))                  - 0.25*creal(clog(xfe3tet)) - creal(clog(xca2m2)) + creal(clog(xna1m2))                  - 0.5*creal(clog(xrest_notNa))                  - 0.5*creal(clog(xm1_notMgFe)) + creal(clog(xm2_notNa))) +                ((HX6)-T*(SX6)+Pv*(VX6)) +                ((HX2X6)-T*(SX2X6)+Pv*(VX2X6))*r0 +                ((HX3X6)-T*(SX3X6)+Pv*(VX3X6))*r1 +                ((HX4X6)-T*(SX4X6)+Pv*(VX4X6))*r2 +                ((HX5X6)-T*(SX5X6)+Pv*(VX5X6))*r3 +                ((HX6X6)-T*(SX6X6)+Pv*(VX6X6))*r4*2.0 +                ((HX6X7)-T*(SX6X7)+Pv*(VX6X7))*r5 +                ((HX6S1)-T*(SX6S1)+Pv*(VX6S1))*s0 +                ((HX6S2)-T*(SX6S2)+Pv*(VX6S2))*s1 +                ((HX2X6X7)+Pv*(VX2X6X7))*r0*r5 +                ((HX2X6S2)+Pv*(VX2X6S2))*r0*s1 +                ((HX3X6X7)+Pv*(VX3X6X7))*r1*r5 +                ((HX3X6S2)+Pv*(VX3X6S2))*r1*s1 +                ((HX4X6X7)+Pv*(VX4X6X7))*r2*r5 +                ((HX4X6S2)+Pv*(VX4X6S2))*r2*s1 +                ((HX5X6X7)+Pv*(VX5X6X7))*r3*r5 +                ((HX5X6S2)+Pv*(VX5X6S2))*r3*s1 +                ((HX6X6X7)+Pv*(VX6X6X7))*r4*r5*2.0 +                ((HX6X6S2)+Pv*(VX6X6S2))*r4*s1*2.0 +                ((HX6X7X7)-T*(SX6X7X7)+Pv*(VX6X7X7))*r5*r5 +                ((HX6X7S2)-T*(SX6X7S2)+Pv*(VX6X7S2))*r5*s1 +                ((HX6S2S2)+Pv*(VX6S2S2))*s1*s1;
    double dgdr5 = 0.5*Rgas*T*(creal(clog(xmg2m1)) - 2.0*creal(clog(xca2m2))                  + creal(clog(xmg2m2)) + creal(clog(xfe2m2/xfe2m1))) +                ((HX7)-T*(SX7)+Pv*(VX7)) +                ((HX2X7)-T*(SX2X7)+Pv*(VX2X7))*r0 +                ((HX3X7)-T*(SX3X7)+Pv*(VX3X7))*r1 +                ((HX4X7)-T*(SX4X7)+Pv*(VX4X7))*r2 +                ((HX5X7)-T*(SX5X7)+Pv*(VX5X7))*r3 +                ((HX6X7)-T*(SX6X7)+Pv*(VX6X7))*r4 +                ((HX7X7)-T*(SX7X7)+Pv*(VX7X7))*r5*2.0 +                ((HX7S1)-T*(SX7S1)+Pv*(VX7S1))*s0 +                ((HX7S2)-T*(SX7S2)+Pv*(VX7S2))*s1 +                ((HX2X2X7)+Pv*(VX2X2X7))*r0*r0 +                ((HX2X3X7)+Pv*(VX2X3X7))*r0*r1 +                ((HX2X4X7)+Pv*(VX2X4X7))*r0*r2 +                ((HX2X5X7)+Pv*(VX2X5X7))*r0*r3 +                ((HX2X6X7)+Pv*(VX2X6X7))*r0*r4 +                ((HX2X7X7)-T*(SX2X7X7)+Pv*(VX2X7X7))*r0*r5*2.0 +                ((HX2X7S2)+Pv*(VX2X7S2))*r0*s1 +                ((HX3X3X7)+Pv*(VX3X3X7))*r1*r1 +                ((HX3X4X7)+Pv*(VX3X4X7))*r1*r2 +                ((HX3X5X7)+Pv*(VX3X5X7))*r1*r3 +                ((HX3X6X7)+Pv*(VX3X6X7))*r1*r4 +                ((HX3X7X7)-T*(SX3X7X7)+Pv*(VX3X7X7))*r1*r5*2.0 +                ((HX3X7S2)-T*(SX3X7S2)+Pv*(VX3X7S2))*r1*s1 +                ((HX4X4X7)+Pv*(VX4X4X7))*r2*r2 +                ((HX4X5X7)+Pv*(VX4X5X7))*r2*r3 +                ((HX4X6X7)+Pv*(VX4X6X7))*r2*r4 +                ((HX4X7X7)-T*(SX4X7X7)+Pv*(VX4X7X7))*r2*r5*2.0 +                ((HX4X7S2)-T*(SX4X7S2)+Pv*(VX4X7S2))*r2*s1 +                ((HX5X5X7)+Pv*(VX5X5X7))*r3*r3 +                ((HX5X6X7)+Pv*(VX5X6X7))*r3*r4 +                ((HX5X7X7)-T*(SX5X7X7)+Pv*(VX5X7X7))*r3*r5*2.0 +                ((HX5X7S2)-T*(SX5X7S2)+Pv*(VX5X7S2))*r3*s1 +                ((HX6X6X7)+Pv*(VX6X6X7))*r4*r4 +                ((HX6X7X7)-T*(SX6X7X7)+Pv*(VX6X7X7))*r4*r5*2.0 +                ((HX6X7S2)-T*(SX6X7S2)+Pv*(VX6X7S2))*r4*s1 +                ((HX7X7X7)-T*(SX7X7X7)+Pv*(VX7X7X7))*r5*r5*3.0 +                ((HX7X7S2)-T*(SX7X7S2)+Pv*(VX7X7S2))*r5*s1*2.0 +                ((HX7S2S2)+Pv*(VX7S2S2))*s1*s1;

    /* ends[]: see GH_cpx_pure_ES_G's header comment - DI/EN/HD/CA/CF/JD are
       exactly 0 (proven, T,P-independent), only essenite is nonzero and
       needs a (cheap, 1D) solve. */
    double es_end = GH_cpx_pure_ES_G(T,Rgas,Pv);

    double Emix = p5*es_end;
    double Gex = Graw - Emix;

    double dGex[7];
    dGex[0] = 0.0;
    dGex[1] = dgdr5;
    dGex[2] = dgdr0;
    dGex[3] = dgdr1;
    dGex[4] = dgdr2;
    dGex[5] = dgdr3 - es_end;
    dGex[6] = 0.5*dgdr1 - 0.5*dgdr2 + 0.5*dgdr3 + dgdr4;

    for (int i = 0; i < 7; i++){
        double mu = Gex;
        for (int k = 0; k < 7; k++){
            double delta_ik = (i==k) ? 1.0 : 0.0;
            mu += (delta_ik - p[k])*dGex[k];
        }
        mu_Gex[i] = mu/1000.0;
        d->sf[i]  = p[i];
    }

    d->sum_apep = 0.0;
    for (int i = 0; i < 7; i++){ d->sum_apep += d->ape[i]*p[i]; }
    d->factor = d->fbc/d->sum_apep;

    d->df_raw = 0.0;
    for (int i = 0; i < 7; i++){ d->df_raw += (mu_Gex[i] + gb[i])*p[i]; }
    d->df = d->df_raw * d->factor;

    if (grad){
        for (int i = 0; i < 7; i++){
            grad[i] = (mu_Gex[i] + gb[i])*d->factor - (d->df_raw*d->factor*(d->ape[i]/d->sum_apep));
        }
    }
    return d->df;
}

/**
    Orthopyroxene: a full, independent port of real MELTS'
    orthopyroxene.c - NOT the same energetics as clinopyroxene (an earlier
    version of this function just called obj_gh_cpx, reasoning both source
    files hardcode clino=TRUE; that was wrong - see below - a real bug,
    found via [[gh-gexcess-verification]] and fixed 2026-07-14).

    Root cause of the earlier mistake: both clinopyroxene.c and
    orthopyroxene.c gate their own isClino() behind "#ifdef ISCLINO",
    immediately followed by "#undef ISCLINO" in the SAME file before
    isClino() is defined - so the conditional branch is dead code in BOTH
    files, and each file's own "#else" fallback fires unconditionally:
    clinopyroxene.c's isClino() always returns TRUE, orthopyroxene.c's
    always returns FALSE. Since real cpx and real opx are literally the
    SAME source-code machinery (same macro names, same formula structure)
    parameterized by "clino", real cpx and real opx genuinely have
    DIFFERENT calibrated energetics - not just the "clino ? cWHFEMG :
    oWHFEMG"-style W-parameters an initial inspection suggests, but a much
    larger effect: real clinopyroxene.c ALSO zeroes 7 "vertex" corrections
    H1..H7/S1..S7/V1..V7 when clino=TRUE and activates a SEPARATE set of
    "pigeonite" corrections pH1..pH7 instead (Ghiorso's "sliding standard
    state for pigeonite", v3.1-1) - for opx (clino=FALSE), this reverses:
    H1..H7 (built from DH1..DH4/DcTOoH3..H6) become the active corrections
    and pH1..pH7 vanish. This affects essentially all ~245 Taylor
    coefficients below, not just a handful of W-parameters.

    All 244 H/S/V Taylor coefficients below were generated the same way as
    obj_gh_cpx's own (see that function's header comment) - mechanically
    transcribing every relevant #define macro from clinopyroxene.c (text-
    identical to orthopyroxene.c) into a small Python script and evaluating
    numerically with clino=False, NOT hand-re-derived. That transcription
    was cross-validated first: run with clino=True, all 244 outputs
    matched obj_gh_cpx's own already-verified constants exactly (0
    mismatches) - see scratchpad cpx_opx_taylor.py, gexcess_verification/.

    Unlike cpx (where DI_G=EN_G=HD_G=CA_G=CF_G=JD_G=0 exactly, since clino
    zeroes every H_i/S_i/V_i vertex, leaving only essenite's own genuine
    Al3+/Fe3+ ordering entropy nonzero), for opx ALL 7 pure-endmember
    reference values ("ends[]", matching real orthopyroxene.c's own
    purePyx) are nonzero, since H1..H7 are now active - see GH_opx_ends
    below. Emix is therefore sum_i(p_i*ends[i]) over all 7 endmembers, not
    just p[5]*es_end.

    p<->r basis and the joint s0/s1 order-parameter machinery
    (site-fractions, composite-quantity clamping, bounds, starting guess)
    are pure composition/site-fraction bookkeeping, NOT calibration-
    specific - reused unchanged from GH_cpx_site_fracs/GH_cpx_composites/
    GH_cpx_bounds_s0/GH_cpx_bounds_s1/GH_cpx_guess_s.
*/
/** opx's own Taylor-coefficient block, needed by GH_opx_dgds/GH_opx_d2gds/
    GH_opx_G_at/GH_opx_ends below. Declared as function-local statics
    (duplicated per function, not a shared file-scope block) because C has
    no namespacing - a file-scope "static const double HX2 = ...;" here
    would collide with cpx's own already-declared file-scope HX2. Same
    244 values as GH_opx_ends' own copy (see obj_gh_opx's header comment
    for how they were generated/verified) - only the subset each function
    actually needs is declared, to keep -Wunused-variable quiet. */
static void GH_opx_dgds(double r0,double r1,double r2,double r3,double r4,double r5,double s0,double s1,
                         double T,double Rgas,double Pv,double *dgds0,double *dgds1){
    static const double HS1 = -2677.12, HS2 = -7635.799999999999, HS2S2 = 15637.7;
    static const double HX2S1 = -3819.6450000000004, HX2S2 = -10773.800000000001, HX2S2S2 = 209.19999999999982;
    static const double HX3S1 = -19744.440000000002, HX3S2 = 5934.734999999998;
    static const double HX4S1 = 3347.0599999999986, HX4S2 = 6366.1749999999965;
    static const double HX5S1 = -2065.5, HX5S2 = -4133.060000000004;
    static const double HX6S1 = -7322, HX6S2 = 5798.579999999997;
    static const double HX7S1 = 4240.940000000002, HX7S2 = 2782.3600000000006, HX7S2S2 = -20501.6;
    static const double HS1S1 = -8917.5, HS1S2 = 7804.93;
    static const double HX2X2S2 = -418.39999999999964, HX2X3S2 = 10041.6, HX2X4S2 = 10041.6;
    static const double HX2X5S2 = 10041.6, HX2X6S2 = 5020.8, HX2X7S2 = 20501.6;
    static const double HX3X3S2 = 2719.6, HX3X4S2 = 4288.599999999999, HX3X5S2 = 6589.8;
    static const double HX3X6S2 = 2719.6, HX3X7S2 = 1255.199999999999;
    static const double HX4X4S2 = 2719.6, HX4X5S2 = 6589.8, HX4X6S2 = 3870.2, HX4X7S2 = 1255.199999999999;
    static const double HX5X5S2 = 5020.8, HX5X6S2 = 5020.8, HX5X7S2 = 5857.5999999999985;
    static const double HX6X6S2 = 1542.85, HX6X7S2 = 2928.7999999999993;
    static const double VS1 = 0, VS2 = -0.09523640000000068, VS2S2 = -0.14957800000000002;
    static const double VX2S1 = 0, VX2S2 = 0.142256, VX2S2S2 = -0.007322;
    static const double VX3S1 = 0, VX3S2 = 0.07267179999999968;
    static const double VX4S1 = 0, VX4S2 = 0.07267179999999968;
    static const double VX5S1 = 0, VX5S2 = 0.10141159999999935;
    static const double VX6S1 = 0, VX6S2 = 0.050705799999999676;
    static const double VX7S1 = 0, VX7S2 = 0.022016400000000658, VX7S2S2 = 0.15690000000000004;
    static const double VS1S1 = 0, VS1S2 = 0;
    static const double VX2X2S2 = 0.014644, VX2X3S2 = -0.071128, VX2X4S2 = -0.071128;
    static const double VX2X5S2 = -0.071128, VX2X6S2 = -0.035564, VX2X7S2 = -0.15690000000000004;
    static const double VX3X3S2 = -0.025104, VX3X4S2 = -0.044978000000000004, VX3X5S2 = -0.055438;
    static const double VX3X6S2 = -0.025104, VX3X7S2 = -0.043932;
    static const double VX4X4S2 = -0.025104, VX4X5S2 = -0.055438, VX4X6S2 = -0.030334, VX4X7S2 = -0.043932;
    static const double VX5X5S2 = -0.035564, VX5X6S2 = -0.035564, VX5X7S2 = -0.064852;
    static const double VX6X6S2 = -0.010198500000000001, VX6X7S2 = -0.032426;
    /* all exactly 0 in opx's Taylor set (SS1/SX*S1/SS2S2/SX*X7S2, matching
       cpx's own identical zero pattern - see cpx_opx_taylor.py output) */
    static const double SS1_OPX = 0, SX2S1_OPX = 0, SX3S1_OPX = 0, SX4S1_OPX = 0;
    static const double SX5S1_OPX = 0, SX6S1_OPX = 0, SX7S1_OPX = 0, SS2S2_OPX = 0;
    static const double SX3X7S2_OPX = 0, SX4X7S2_OPX = 0, SX5X7S2_OPX = 0;
    static const double SX6X7S2_OPX = 0, SX7X7S2_OPX = 0;
    double xal3m1,xfe2m1,xfe3m1,xmg2m1,xti4m1,xca2m2,xfe2m2,xmg2m2,xna1m2,xal3tet,xfe3tet,xsi4tet;
    GH_cpx_site_fracs(r0,r1,r2,r3,r4,r5,s0,s1,&xal3m1,&xfe2m1,&xfe3m1,&xmg2m1,&xti4m1,&xca2m2,&xfe2m2,&xmg2m2,&xna1m2,&xal3tet,&xfe3tet,&xsi4tet);
    *dgds0 = 0.5*Rgas*T*(creal(clog(xfe3m1))-creal(clog(xal3m1))+creal(clog(xal3tet))-creal(clog(xfe3tet))) +                ((HS1)-T*(SS1_OPX)+Pv*(VS1)) +                ((HX2S1)-T*(SX2S1_OPX)+Pv*(VX2S1))*r0 +                ((HX3S1)-T*(SX3S1_OPX)+Pv*(VX3S1))*r1 +                ((HX4S1)-T*(SX4S1_OPX)+Pv*(VX4S1))*r2 +                ((HX5S1)-T*(SX5S1_OPX)+Pv*(VX5S1))*r3 +                ((HX6S1)-T*(SX6S1_OPX)+Pv*(VX6S1))*r4 +                ((HX7S1)-T*(SX7S1_OPX)+Pv*(VX7S1))*r5 +                ((HS1S1)-T*(SS1S1)+Pv*(VS1S1))*s0*2.0 +                ((HS1S2)-T*(SS1S2)+Pv*(VS1S2))*s1;
    *dgds1 = 0.5*Rgas*T*(creal(clog(xmg2m1))-creal(clog(xfe2m1))+creal(clog(xfe2m2))-creal(clog(xmg2m2))) +                ((HS2)-T*(GH_opx_SS2)+Pv*(VS2)) +                ((HX2S2)-T*(SX2S2)+Pv*(VX2S2))*r0 +                ((HX3S2)-T*(SX3S2)+Pv*(VX3S2))*r1 +                ((HX4S2)-T*(SX4S2)+Pv*(VX4S2))*r2 +                ((HX5S2)-T*(SX5S2)+Pv*(VX5S2))*r3 +                ((HX6S2)-T*(SX6S2)+Pv*(VX6S2))*r4 +                ((HX7S2)-T*(SX7S2)+Pv*(VX7S2))*r5 +                ((HS1S2)-T*(SS1S2)+Pv*(VS1S2))*s0 +                ((HS2S2)-T*(SS2S2_OPX)+Pv*(VS2S2))*s1*2.0 +                ((HX2X2S2)+Pv*(VX2X2S2))*r0*r0 +                ((HX2X3S2)+Pv*(VX2X3S2))*r0*r1 +                ((HX2X4S2)+Pv*(VX2X4S2))*r0*r2 +                ((HX2X5S2)+Pv*(VX2X5S2))*r0*r3 +                ((HX2X6S2)+Pv*(VX2X6S2))*r0*r4 +                ((HX2X7S2)+Pv*(VX2X7S2))*r0*r5 +                ((HX2S2S2)+Pv*(VX2S2S2))*r0*s1*2.0 +                ((HX3X3S2)+Pv*(VX3X3S2))*r1*r1 +                ((HX3X4S2)+Pv*(VX3X4S2))*r1*r2 +                ((HX3X5S2)+Pv*(VX3X5S2))*r1*r3 +                ((HX3X6S2)+Pv*(VX3X6S2))*r1*r4 +                ((HX3X7S2)-T*(SX3X7S2_OPX)+Pv*(VX3X7S2))*r1*r5 +                ((HX3S2S2)+Pv*(VX3S2S2))*r1*s1*2.0 +                ((HX4X4S2)+Pv*(VX4X4S2))*r2*r2 +                ((HX4X5S2)+Pv*(VX4X5S2))*r2*r3 +                ((HX4X6S2)+Pv*(VX4X6S2))*r2*r4 +                ((HX4X7S2)-T*(SX4X7S2_OPX)+Pv*(VX4X7S2))*r2*r5 +                ((HX4S2S2)+Pv*(VX4S2S2))*r2*s1*2.0 +                ((HX5X5S2)+Pv*(VX5X5S2))*r3*r3 +                ((HX5X6S2)+Pv*(VX5X6S2))*r3*r4 +                ((HX5X7S2)-T*(SX5X7S2_OPX)+Pv*(VX5X7S2))*r3*r5 +                ((HX5S2S2)+Pv*(VX5S2S2))*r3*s1*2.0 +                ((HX6X6S2)+Pv*(VX6X6S2))*r4*r4 +                ((HX6X7S2)-T*(SX6X7S2_OPX)+Pv*(VX6X7S2))*r4*r5 +                ((HX6S2S2)+Pv*(VX6S2S2))*r4*s1*2.0 +                ((HX7X7S2)-T*(SX7X7S2_OPX)+Pv*(VX7X7S2))*r5*r5 +                ((HX7S2S2)+Pv*(VX7S2S2))*r5*s1*2.0;
}

static void GH_opx_d2gds(double r0,double r1,double r2,double r3,double r4,double r5,double s0,double s1,
                          double T,double Rgas,double Pv,double *d00,double *d01,double *d11){
    static const double HS1S1 = -8917.5, HS1S2 = 7804.93, HS2S2 = 15637.7;
    static const double VS1S1 = 0, VS1S2 = 0, VS2S2 = -0.14957800000000002;
    static const double HX2S2S2 = 209.19999999999982, HX3S2S2 = -10826.1, HX4S2S2 = -10826.1;
    static const double HX5S2S2 = -10250.8, HX6S2S2 = -5125.4, HX7S2S2 = -20501.6;
    static const double VX2S2S2 = -0.007322, VX3S2S2 = 0.08106500000000001, VX4S2S2 = 0.08106500000000001;
    static const double VX5S2S2 = 0.07845000000000002, VX6S2S2 = 0.03922500000000001, VX7S2S2 = 0.15690000000000004;
    static const double SS2S2_OPX = 0;
    double xal3m1,xfe2m1,xfe3m1,xmg2m1,xti4m1,xca2m2,xfe2m2,xmg2m2,xna1m2,xal3tet,xfe3tet,xsi4tet;
    GH_cpx_site_fracs(r0,r1,r2,r3,r4,r5,s0,s1,&xal3m1,&xfe2m1,&xfe3m1,&xmg2m1,&xti4m1,&xca2m2,&xfe2m2,&xmg2m2,&xna1m2,&xal3tet,&xfe3tet,&xsi4tet);
    *d00 = Rgas*T*(0.25/xfe3m1+0.25/xal3m1+0.125/xal3tet+0.125/xfe3tet) + 2.0*((HS1S1)-T*(0.0)+Pv*(VS1S1));
    *d01 = ((HS1S2)-T*(0.0)+Pv*(VS1S2));
    *d11 = Rgas*T*(0.25/xmg2m1+0.25/xfe2m1+0.25/xfe2m2+0.25/xmg2m2) + 2.0*((HS2S2)-T*(SS2S2_OPX)+Pv*(VS2S2))
         + 2.0*r0*((HX2S2S2)+Pv*(VX2S2S2)) + 2.0*r1*((HX3S2S2)+Pv*(VX3S2S2)) + 2.0*r2*((HX4S2S2)+Pv*(VX4S2S2))
         + 2.0*r3*((HX5S2S2)+Pv*(VX5S2S2)) + 2.0*r4*((HX6S2S2)+Pv*(VX6S2S2)) + 2.0*r5*((HX7S2S2)+Pv*(VX7S2S2));
}
static void GH_opx_solve_s_from(double r0,double r1,double r2,double r3,double r4,double r5,double T,double Rgas,double Pv,
                                 double s0init,double s1init,double *s0o,double *s1o){
    const int maxIter = 1000;
    double eps_s = 1.0e-8;
    double s0=s0init, s1=s1init;
    double s0Old=2.0, s1Old=2.0;
    int iter = 0;
    while ( (fabs(s0-s0Old) > 10.0*2.220446049250313e-16 ||
             fabs(s1-s1Old) > 10.0*2.220446049250313e-16) && iter < maxIter ){
        double dgds0, dgds1;
        GH_opx_dgds(r0,r1,r2,r3,r4,r5,s0,s1,T,Rgas,Pv,&dgds0,&dgds1);
        double d00, d01, d11;
        GH_opx_d2gds(r0,r1,r2,r3,r4,r5,s0,s1,T,Rgas,Pv,&d00,&d01,&d11);

        s0Old = s0; s1Old = s1;

        double det = d00*d11 - d01*d01;
        double deltaS0 = -( d11*dgds0 - d01*dgds1)/det;
        double deltaS1 = -(-d01*dgds0 + d00*dgds1)/det;

        double sN0 = s0Old + deltaS0, sN1 = s1Old + deltaS1;
        double lo, hi;
        GH_cpx_bounds_s0(r1,r2,r3,r4,eps_s,&lo,&hi);
        if (sN0 < lo){ sN0 = lo; } if (sN0 > hi){ sN0 = hi; }
        GH_cpx_bounds_s1(r0,r1,r2,r3,r4,r5,eps_s,&lo,&hi);
        if (sN1 < lo){ sN1 = lo; } if (sN1 > hi){ sN1 = hi; }
        s0 = sN0; s1 = sN1;
        iter++;
    }
    *s0o = s0; *s1o = s1;
}
static double GH_opx_G_at(double r0,double r1,double r2,double r3,double r4,double r5,double s0,double s1,
                           double T,double Rgas,double Pv){
    static const double H0 = 2331.7432000000003, HS1 = -2677.12, HS2 = -7635.799999999999;
    static const double HS1S1 = -8917.5, HS1S2 = 7804.93, HS2S2 = 15637.7;
    static const double HX2 = 8086.416799999999, HX2S1 = -3819.6450000000004, HX2S2 = -10773.800000000001;
    static const double HX2X2 = -8368, HX2X3 = -9618, HX2X4 = -18916.96, HX2X5 = -21486.595, HX2X6 = -6093.8275, HX2X7 = 7050.039999999999;
    static const double HX3 = 39105.4968, HX3S1 = -19744.440000000002, HX3S2 = 5934.734999999998;
    static const double HX3X3 = -16318.000000000002, HX3X4 = -18357, HX3X5 = -5465.440000000004, HX3X6 = 8160.4000000000015, HX3X7 = -2484.9950000000044;
    static const double HX4 = 41615.4968, HX4S1 = 3347.0599999999986, HX4S2 = 6366.1749999999965;
    static const double HX4X4 = -18828, HX4X5 = -15951.939999999999, HX4X6 = 23252.149999999994, HX4X7 = 2166.574999999997;
    static const double HX5 = 46384.6168, HX5S1 = -2065.5, HX5S2 = -4133.060000000004;
    static const double HX5X5 = -9937, HX5X6 = -16226.249999999996, HX5X7 = -8988.800000000007;
    static const double HX6 = -1826.9316, HX6S1 = -7322, HX6S2 = 5798.579999999997;
    static const double HX6X6 = 8735.874999999996, HX6X7 = -13980.114999999998;
    static const double HX7 = 11203.4968, HX7S1 = 4240.940000000002, HX7S2 = 2782.3600000000006, HX7X7 = -28168.78;
    static const double HX2X2X7 = 418.39999999999964, HX2X2S2 = -418.39999999999964;
    static const double HX2X3X7 = -836.7999999999993, HX2X3S2 = 10041.6;
    static const double HX2X4X7 = -836.7999999999993, HX2X4S2 = 10041.6;
    static const double HX2X5X7 = -836.7999999999993, HX2X5S2 = 10041.6;
    static const double HX2X6X7 = -418.39999999999964, HX2X6S2 = 5020.8;
    static const double HX2X7X7 = -2301.199999999997, HX2X7S2 = 20501.6, HX2S2S2 = 209.19999999999982;
    static const double HX3X3X7 = -2719.6, HX3X3S2 = 2719.6;
    static const double HX3X4X7 = -6589.8, HX3X4S2 = 4288.599999999999;
    static const double HX3X5X7 = -4288.599999999999, HX3X5S2 = 6589.8;
    static const double HX3X6X7 = -2719.6, HX3X6S2 = 2719.6;
    static const double HX3X7X7 = -4236.3, HX3X7S2 = 1255.199999999999, HX3S2S2 = -10826.1;
    static const double HX4X4X7 = -2719.6, HX4X4S2 = 2719.6;
    static const double HX4X5X7 = -4288.599999999999, HX4X5S2 = 6589.8;
    static const double HX4X6X7 = -1568.9999999999998, HX4X6S2 = 3870.2;
    static const double HX4X7X7 = -4236.3, HX4X7S2 = 1255.199999999999, HX4S2S2 = -10826.1;
    static const double HX5X5X7 = -418.39999999999964, HX5X5S2 = 5020.8;
    static const double HX5X6X7 = -418.39999999999964, HX5X6S2 = 5020.8;
    static const double HX5X7X7 = -209.19999999999982, HX5X7S2 = 5857.5999999999985, HX5S2S2 = -10250.8;
    static const double HX6X6X7 = 183.05000000000007, HX6X6S2 = 1542.85;
    static const double HX6X7X7 = -104.59999999999991, HX6X7S2 = 2928.7999999999993, HX6S2S2 = -5125.4;
    static const double HX7X7X7 = 1255.199999999999, HX7X7S2 = -8368, HX7S2S2 = -20501.6;
    static const double S0 = 1.3807200000000002;
    static const double SX2 = -0.8284320000000002, SX3 = -3.38072, SX4 = -3.38072, SX5 = -3.38072, SX6 = -1.69036, SX7 = -4.78469688;
    static const double SS1 = 0, SS1S1 = 0, SS1S2 = 0;
    static const double SX2S1 = 0, SX2S2 = 0, SX2X2 = 0, SX2X3 = 0, SX2X4 = 0, SX2X5 = 0, SX2X6 = 0, SX2X7 = -1.15955376;
    static const double SX3S1 = 0, SX3S2 = -0.28988844, SX3X3 = 0, SX3X4 = 0, SX3X5 = 0, SX3X6 = 0, SX3X7 = -0.28988844;
    static const double SX4S1 = 0, SX4S2 = -0.28988844, SX4X4 = 0, SX4X5 = 0, SX4X6 = 0, SX4X7 = -0.28988844;
    static const double SX5S1 = 0, SX5S2 = -0.57977688, SX5X5 = 0, SX5X6 = 0, SX5X7 = -0.57977688;
    static const double SX6S1 = 0, SX6S2 = -0.28988844, SX6X6 = 0, SX6X7 = -0.28988844;
    static const double SX7X7 = 0.57977688, SX7S2 = 0.57977688;
    static const double SX2X2X7 = 0, SX2X2S2 = 0, SX2X3X7 = 0, SX2X3S2 = 0, SX2X4X7 = 0, SX2X4S2 = 0, SX2X5X7 = 0, SX2X5S2 = 0;
    static const double SX2X6X7 = 0, SX2X6S2 = 0, SX2X7X7 = 0, SX2X7S2 = 0, SX2S2S2 = 0;
    static const double SX3X3X7 = 0, SX3X3S2 = 0, SX3X4X7 = 0, SX3X4S2 = 0, SX3X5X7 = 0, SX3X5S2 = 0, SX3X6X7 = 0, SX3X6S2 = 0;
    static const double SX3X7X7 = 0, SX3X7S2 = 0, SX3S2S2 = 0;
    static const double SX4X4X7 = 0, SX4X4S2 = 0, SX4X5X7 = 0, SX4X5S2 = 0, SX4X6X7 = 0, SX4X6S2 = 0;
    static const double SX4X7X7 = 0, SX4X7S2 = 0, SX4S2S2 = 0;
    static const double SX5X5X7 = 0, SX5X5S2 = 0, SX5X6X7 = 0, SX5X6S2 = 0, SX5X7X7 = 0, SX5X7S2 = 0, SX5S2S2 = 0;
    static const double SX6X6X7 = 0, SX6X6S2 = 0, SX6X7X7 = 0, SX6X7S2 = 0, SX6S2S2 = 0;
    static const double SX7X7X7 = 0, SX7X7S2 = 0, SX7S2S2 = 0;
    static const double V0 = 0.075312;
    static const double VS1 = 0, VS1S1 = 0, VS1S2 = 0, VS2 = -0.09523640000000068, VS2S2 = -0.14957800000000002;
    static const double VX2 = 0.03710980000000141, VX2S1 = 0, VX2S2 = 0.142256;
    static const double VX2X2 = -0.014121, VX2X3 = -0.0070605, VX2X4 = -0.0070605, VX2X5 = -0.014121, VX2X6 = -0.0070605, VX2X7 = -0.01317580000000132;
    static const double VX3 = -0.126602, VX3S1 = 0, VX3S2 = 0.07267179999999968;
    static const double VX3X3 = 0, VX3X4 = 0, VX3X5 = 0, VX3X6 = 0, VX3X7 = -0.06121620000000034;
    static const double VX4 = -0.126602, VX4S1 = 0, VX4S2 = 0.07267179999999968;
    static const double VX4X4 = 0, VX4X5 = 0, VX4X6 = 0, VX4X7 = -0.06121620000000034;
    static const double VX5 = -0.126602, VX5S1 = 0, VX5S2 = 0.10141159999999935;
    static const double VX5X5 = 0, VX5X6 = 0, VX5X7 = -0.03247640000000067;
    static const double VX6 = -0.063301, VX6S1 = 0, VX6S2 = 0.050705799999999676;
    static const double VX6X6 = 0, VX6X7 = -0.016238200000000334;
    static const double VX7 = -0.14879160000000066, VX7S1 = 0, VX7S2 = 0.022016400000000658, VX7X7 = -0.01250159999999934;
    static const double VX2X2S2 = 0.014644, VX2X2X7 = -0.014644;
    static const double VX2X3S2 = -0.071128, VX2X3X7 = 0.029288;
    static const double VX2X4S2 = -0.071128, VX2X4X7 = 0.029288;
    static const double VX2X5S2 = -0.071128, VX2X5X7 = 0.029288;
    static const double VX2X6S2 = -0.035564, VX2X6X7 = 0.014644;
    static const double VX2X7S2 = -0.15690000000000004, VX2X7X7 = 0.080542, VX2S2S2 = -0.007322;
    static const double VX3X3S2 = -0.025104, VX3X3X7 = 0.025104;
    static const double VX3X4S2 = -0.044978000000000004, VX3X4X7 = 0.055438;
    static const double VX3X5S2 = -0.055438, VX3X5X7 = 0.044978000000000004;
    static const double VX3X6S2 = -0.025104, VX3X6X7 = 0.025104;
    static const double VX3X7S2 = -0.043932, VX3X7X7 = 0.025627, VX3S2S2 = 0.08106500000000001;
    static const double VX4X4S2 = -0.025104, VX4X4X7 = 0.025104;
    static const double VX4X5S2 = -0.055438, VX4X5X7 = 0.044978000000000004;
    static const double VX4X6S2 = -0.030334, VX4X6X7 = 0.019874000000000003;
    static const double VX4X7S2 = -0.043932, VX4X7X7 = 0.025627, VX4S2S2 = 0.08106500000000001;
    static const double VX5X5S2 = -0.035564, VX5X5X7 = 0.014644;
    static const double VX5X6S2 = -0.035564, VX5X6X7 = 0.014644;
    static const double VX5X7S2 = -0.064852, VX5X7X7 = 0.007322, VX5S2S2 = 0.07845000000000002;
    static const double VX6X6S2 = -0.010198500000000001, VX6X6X7 = 0.0023534999999999997;
    static const double VX6X7S2 = -0.032426, VX6X7X7 = 0.003661, VX6S2S2 = 0.03922500000000001;
    static const double VX7X7S2 = 0.012552000000000008, VX7X7X7 = -0.043932, VX7S2S2 = 0.15690000000000004;

    double xal3m1,xfe2m1,xfe3m1,xmg2m1,xti4m1,xca2m2,xfe2m2,xmg2m2,xna1m2,xal3tet,xfe3tet,xsi4tet;
    GH_cpx_site_fracs(r0,r1,r2,r3,r4,r5,s0,s1,&xal3m1,&xfe2m1,&xfe3m1,&xmg2m1,&xti4m1,&xca2m2,&xfe2m2,&xmg2m2,&xna1m2,&xal3tet,&xfe3tet,&xsi4tet);
    double xm1_MgFeNoTi,xm1_notTi,xm1_MgFe,xtet_notSi,xrest_notNa,xm1_notMgFe,xm2_notNa;
    GH_cpx_composites(xmg2m1,xfe2m1,xti4m1,xna1m2,xsi4tet,&xm1_MgFeNoTi,&xm1_notTi,&xm1_MgFe,&xtet_notSi,&xrest_notNa,&xm1_notMgFe,&xm2_notNa);
    double SIC = -Rgas*(xmg2m1*creal(clog(xmg2m1))                + xfe2m1*creal(clog(xfe2m1))                + xal3m1*creal(clog(xal3m1))                + xfe3m1*creal(clog(xfe3m1))                + xti4m1*creal(clog(xti4m1))                + (xm1_MgFeNoTi)*creal(clog(xm1_MgFeNoTi))                - (xm1_notTi)*creal(clog(xm1_notTi))                - (xm1_MgFe)*creal(clog(xm1_MgFe))                - 2.0*(xtet_notSi)*creal(clog(xtet_notSi))                + 2.0*xal3tet*creal(clog(xal3tet))                + 2.0*xfe3tet*creal(clog(xfe3tet))                + xca2m2*creal(clog(xca2m2))                + xna1m2*creal(clog(xna1m2))                + xmg2m2*creal(clog(xmg2m2))                + xfe2m2*creal(clog(xfe2m2))                + (xrest_notNa)*creal(clog(xrest_notNa))                - (xm1_notMgFe)*creal(clog(xm1_notMgFe))                - (xm2_notNa)*creal(clog(xm2_notNa)) );
    return -T*(SIC) + (H0)-T*(S0)+Pv*(V0) +              ((HX2)-T*(SX2)+Pv*(VX2))*r0 +              ((HX3)-T*(SX3)+Pv*(VX3))*r1 +              ((HX4)-T*(SX4)+Pv*(VX4))*r2 +              ((HX5)-T*(SX5)+Pv*(VX5))*r3 +              ((HX6)-T*(SX6)+Pv*(VX6))*r4 +              ((HX7)-T*(SX7)+Pv*(VX7))*r5 +              ((HS1)-T*(SS1)+Pv*(VS1))*s0 +              ((HS2)-T*(GH_opx_SS2)+Pv*(VS2))*s1 +              ((HX2X2)-T*(SX2X2)+Pv*(VX2X2))*r0*r0 +              ((HX2X3)-T*(SX2X3)+Pv*(VX2X3))*r0*r1 +              ((HX2X4)-T*(SX2X4)+Pv*(VX2X4))*r0*r2 +              ((HX2X5)-T*(SX2X5)+Pv*(VX2X5))*r0*r3 +              ((HX2X6)-T*(SX2X6)+Pv*(VX2X6))*r0*r4 +              ((HX2X7)-T*(SX2X7)+Pv*(VX2X7))*r0*r5 +              ((HX2S1)-T*(SX2S1)+Pv*(VX2S1))*r0*s0 +              ((HX2S2)-T*(SX2S2)+Pv*(VX2S2))*r0*s1 +              ((HX3X3)-T*(SX3X3)+Pv*(VX3X3))*r1*r1 +              ((HX3X4)-T*(SX3X4)+Pv*(VX3X4))*r1*r2 +              ((HX3X5)-T*(SX3X5)+Pv*(VX3X5))*r1*r3 +              ((HX3X6)-T*(SX3X6)+Pv*(VX3X6))*r1*r4 +              ((HX3X7)-T*(SX3X7)+Pv*(VX3X7))*r1*r5 +              ((HX3S1)-T*(SX3S1)+Pv*(VX3S1))*r1*s0 +              ((HX3S2)-T*(SX3S2)+Pv*(VX3S2))*r1*s1 +              ((HX4X4)-T*(SX4X4)+Pv*(VX4X4))*r2*r2 +              ((HX4X5)-T*(SX4X5)+Pv*(VX4X5))*r2*r3 +              ((HX4X6)-T*(SX4X6)+Pv*(VX4X6))*r2*r4 +              ((HX4X7)-T*(SX4X7)+Pv*(VX4X7))*r2*r5 +              ((HX4S1)-T*(SX4S1)+Pv*(VX4S1))*r2*s0 +              ((HX4S2)-T*(SX4S2)+Pv*(VX4S2))*r2*s1 +              ((HX5X5)-T*(SX5X5)+Pv*(VX5X5))*r3*r3 +              ((HX5X6)-T*(SX5X6)+Pv*(VX5X6))*r3*r4 +              ((HX5X7)-T*(SX5X7)+Pv*(VX5X7))*r3*r5 +              ((HX5S1)-T*(SX5S1)+Pv*(VX5S1))*r3*s0 +              ((HX5S2)-T*(SX5S2)+Pv*(VX5S2))*r3*s1 +              ((HX6X6)-T*(SX6X6)+Pv*(VX6X6))*r4*r4 +              ((HX6X7)-T*(SX6X7)+Pv*(VX6X7))*r4*r5 +              ((HX6S1)-T*(SX6S1)+Pv*(VX6S1))*r4*s0 +              ((HX6S2)-T*(SX6S2)+Pv*(VX6S2))*r4*s1 +              ((HX7X7)-T*(SX7X7)+Pv*(VX7X7))*r5*r5 +              ((HX7S1)-T*(SX7S1)+Pv*(VX7S1))*r5*s0 +              ((HX7S2)-T*(SX7S2)+Pv*(VX7S2))*r5*s1 +              ((HS1S1)-T*(SS1S1)+Pv*(VS1S1))*s0*s0 +              ((HS1S2)-T*(SS1S2)+Pv*(VS1S2))*s0*s1 +              ((HS2S2)-T*(SS2S2)+Pv*(VS2S2))*s1*s1 +              ((HX2X2X7)+Pv*(VX2X2X7))*r0*r0*r5 +              ((HX2X2S2)+Pv*(VX2X2S2))*r0*r0*s1 +              ((HX2X3X7)+Pv*(VX2X3X7))*r0*r1*r5 +              ((HX2X3S2)+Pv*(VX2X3S2))*r0*r1*s1 +              ((HX2X4X7)+Pv*(VX2X4X7))*r0*r2*r5 +              ((HX2X4S2)+Pv*(VX2X4S2))*r0*r2*s1 +              ((HX2X5X7)+Pv*(VX2X5X7))*r0*r3*r5 +              ((HX2X5S2)+Pv*(VX2X5S2))*r0*r3*s1 +              ((HX2X6X7)+Pv*(VX2X6X7))*r0*r4*r5 +              ((HX2X6S2)+Pv*(VX2X6S2))*r0*r4*s1 +              ((HX2X7X7)-T*(SX2X7X7)+Pv*(VX2X7X7))*r0*r5*r5 +              ((HX2X7S2)+Pv*(VX2X7S2))*r0*r5*s1 +              ((HX2S2S2)+Pv*(VX2S2S2))*r0*s1*s1 +              ((HX3X3X7)+Pv*(VX3X3X7))*r1*r1*r5 +              ((HX3X3S2)+Pv*(VX3X3S2))*r1*r1*s1 +              ((HX3X4X7)+Pv*(VX3X4X7))*r1*r2*r5 +              ((HX3X4S2)+Pv*(VX3X4S2))*r1*r2*s1 +              ((HX3X5X7)+Pv*(VX3X5X7))*r1*r3*r5 +              ((HX3X5S2)+Pv*(VX3X5S2))*r1*r3*s1 +              ((HX3X6X7)+Pv*(VX3X6X7))*r1*r4*r5 +              ((HX3X6S2)+Pv*(VX3X6S2))*r1*r4*s1 +              ((HX3X7X7)-T*(SX3X7X7)+Pv*(VX3X7X7))*r1*r5*r5 +              ((HX3X7S2)-T*(SX3X7S2)+Pv*(VX3X7S2))*r1*r5*s1 +              ((HX3S2S2)+Pv*(VX3S2S2))*r1*s1*s1 +              ((HX4X4X7)+Pv*(VX4X4X7))*r2*r2*r5 +              ((HX4X4S2)+Pv*(VX4X4S2))*r2*r2*s1 +              ((HX4X5X7)+Pv*(VX4X5X7))*r2*r3*r5 +              ((HX4X5S2)+Pv*(VX4X5S2))*r2*r3*s1 +              ((HX4X6X7)+Pv*(VX4X6X7))*r2*r4*r5 +              ((HX4X6S2)+Pv*(VX4X6S2))*r2*r4*s1 +              ((HX4X7X7)-T*(SX4X7X7)+Pv*(VX4X7X7))*r2*r5*r5 +              ((HX4X7S2)-T*(SX4X7S2)+Pv*(VX4X7S2))*r2*r5*s1 +              ((HX4S2S2)+Pv*(VX4S2S2))*r2*s1*s1 +              ((HX5X5X7)+Pv*(VX5X5X7))*r3*r3*r5 +              ((HX5X5S2)+Pv*(VX5X5S2))*r3*r3*s1 +              ((HX5X6X7)+Pv*(VX5X6X7))*r3*r4*r5 +              ((HX5X6S2)+Pv*(VX5X6S2))*r3*r4*s1 +              ((HX5X7X7)-T*(SX5X7X7)+Pv*(VX5X7X7))*r3*r5*r5 +              ((HX5X7S2)-T*(SX5X7S2)+Pv*(VX5X7S2))*r3*r5*s1 +              ((HX5S2S2)+Pv*(VX5S2S2))*r3*s1*s1 +              ((HX6X6X7)+Pv*(VX6X6X7))*r4*r4*r5 +              ((HX6X6S2)+Pv*(VX6X6S2))*r4*r4*s1 +              ((HX6X7X7)-T*(SX6X7X7)+Pv*(VX6X7X7))*r4*r5*r5 +              ((HX6X7S2)-T*(SX6X7S2)+Pv*(VX6X7S2))*r4*r5*s1 +              ((HX6S2S2)+Pv*(VX6S2S2))*r4*s1*s1 +              ((HX7X7X7)-T*(SX7X7X7)+Pv*(VX7X7X7))*r5*r5*r5 +              ((HX7X7S2)-T*(SX7X7S2)+Pv*(VX7X7S2))*r5*r5*s1 +              ((HX7S2S2)+Pv*(VX7S2S2))*r5*s1*s1;
}
static double GH_opx_solve_and_G(double r0,double r1,double r2,double r3,double r4,double r5,
                                  double T,double Rgas,double Pv,double *s0o,double *s1o){
    double s0g,s1g; GH_cpx_guess_s(r0,r1,r2,r3,r4,r5,&s0g,&s1g);
    double cs0,cs1; GH_opx_solve_s_from(r0,r1,r2,r3,r4,r5,T,Rgas,Pv,s0g,s1g,&cs0,&cs1);
    *s0o = cs0; *s1o = cs1;
    return GH_opx_G_at(r0,r1,r2,r3,r4,r5,cs0,cs1,T,Rgas,Pv);
}

double obj_gh_opx(unsigned n, const double *x, double *grad, void *SS_ref_db){
    SS_ref *d = (SS_ref *) SS_ref_db;

    double T    = d->T;
    double Rgas = d->R*1000.0;
    double Pv   = d->P - 1.0;
    double *p   = d->p;
    double *gb  = d->gb_lvl;
    double *mu_Gex = d->mu_Gex;

    for (int i = 0; i < 7; i++){ p[i] = x[i]; }
    double p1=p[1], p2=p[2], p3=p[3], p4=p[4], p5=p[5], p6=p[6];

    /* SAME r-basis as cpx (pure composition bookkeeping, not calibration-
       specific) - see obj_gh_cpx's header comment. */
    double r0=p2, r1=p3+0.5*p6, r2=p4-0.5*p6, r3=p5+0.5*p6, r4=p6, r5=p1;

    double s0, s1;
    double Graw = GH_opx_solve_and_G(r0,r1,r2,r3,r4,r5,T,Rgas,Pv,&s0,&s1);

    double xal3m1,xfe2m1,xfe3m1,xmg2m1,xti4m1,xca2m2,xfe2m2,xmg2m2,xna1m2,xal3tet,xfe3tet,xsi4tet;
    GH_cpx_site_fracs(r0,r1,r2,r3,r4,r5,s0,s1,&xal3m1,&xfe2m1,&xfe3m1,&xmg2m1,&xti4m1,&xca2m2,&xfe2m2,&xmg2m2,&xna1m2,&xal3tet,&xfe3tet,&xsi4tet);
    double xm1_MgFeNoTi,xm1_notTi,xm1_MgFe,xtet_notSi,xrest_notNa,xm1_notMgFe,xm2_notNa;
    GH_cpx_composites(xmg2m1,xfe2m1,xti4m1,xna1m2,xsi4tet,&xm1_MgFeNoTi,&xm1_notTi,&xm1_MgFe,&xtet_notSi,&xrest_notNa,&xm1_notMgFe,&xm2_notNa);

    static const double HX2 = 8086.416799999999, HX2S1 = -3819.6450000000004, HX2S2 = -10773.800000000001;
    static const double HX2X2 = -8368, HX2X3 = -9618, HX2X4 = -18916.96, HX2X5 = -21486.595, HX2X6 = -6093.8275, HX2X7 = 7050.039999999999;
    static const double HX3 = 39105.4968, HX3S1 = -19744.440000000002, HX3S2 = 5934.734999999998;
    static const double HX3X3 = -16318.000000000002, HX3X4 = -18357, HX3X5 = -5465.440000000004, HX3X6 = 8160.4000000000015, HX3X7 = -2484.9950000000044;
    static const double HX4 = 41615.4968, HX4S1 = 3347.0599999999986, HX4S2 = 6366.1749999999965;
    static const double HX4X4 = -18828, HX4X5 = -15951.939999999999, HX4X6 = 23252.149999999994, HX4X7 = 2166.574999999997;
    static const double HX5 = 46384.6168, HX5S1 = -2065.5, HX5S2 = -4133.060000000004;
    static const double HX5X5 = -9937, HX5X6 = -16226.249999999996, HX5X7 = -8988.800000000007;
    static const double HX6 = -1826.9316, HX6S1 = -7322, HX6S2 = 5798.579999999997;
    static const double HX6X6 = 8735.874999999996, HX6X7 = -13980.114999999998;
    static const double HX7 = 11203.4968, HX7S1 = 4240.940000000002, HX7S2 = 2782.3600000000006, HX7X7 = -28168.78;
    static const double HX2X2X7 = 418.39999999999964, HX2X2S2 = -418.39999999999964;
    static const double HX2X3X7 = -836.7999999999993, HX2X3S2 = 10041.6;
    static const double HX2X4X7 = -836.7999999999993, HX2X4S2 = 10041.6;
    static const double HX2X5X7 = -836.7999999999993, HX2X5S2 = 10041.6;
    static const double HX2X6X7 = -418.39999999999964, HX2X6S2 = 5020.8;
    static const double HX2X7X7 = -2301.199999999997, HX2X7S2 = 20501.6, HX2S2S2 = 209.19999999999982;
    static const double HX3X3X7 = -2719.6, HX3X3S2 = 2719.6;
    static const double HX3X4X7 = -6589.8, HX3X4S2 = 4288.599999999999;
    static const double HX3X5X7 = -4288.599999999999, HX3X5S2 = 6589.8;
    static const double HX3X6X7 = -2719.6, HX3X6S2 = 2719.6;
    static const double HX3X7X7 = -4236.3, HX3X7S2 = 1255.199999999999, HX3S2S2 = -10826.1;
    static const double HX4X4X7 = -2719.6, HX4X4S2 = 2719.6;
    static const double HX4X5X7 = -4288.599999999999, HX4X5S2 = 6589.8;
    static const double HX4X6X7 = -1568.9999999999998, HX4X6S2 = 3870.2;
    static const double HX4X7X7 = -4236.3, HX4X7S2 = 1255.199999999999, HX4S2S2 = -10826.1;
    static const double HX5X5X7 = -418.39999999999964, HX5X5S2 = 5020.8;
    static const double HX5X6X7 = -418.39999999999964, HX5X6S2 = 5020.8;
    static const double HX5X7X7 = -209.19999999999982, HX5X7S2 = 5857.5999999999985, HX5S2S2 = -10250.8;
    static const double HX6X6X7 = 183.05000000000007, HX6X6S2 = 1542.85;
    static const double HX6X7X7 = -104.59999999999991, HX6X7S2 = 2928.7999999999993, HX6S2S2 = -5125.4;
    static const double HX7X7X7 = 1255.199999999999, HX7X7S2 = -8368, HX7S2S2 = -20501.6;
    static const double SX2 = -0.8284320000000002, SX3 = -3.38072, SX4 = -3.38072, SX5 = -3.38072, SX6 = -1.69036, SX7 = -4.78469688;
    static const double SX2S1 = 0, SX2S2 = 0, SX2X2 = 0, SX2X3 = 0, SX2X4 = 0, SX2X5 = 0, SX2X6 = 0, SX2X7 = -1.15955376;
    static const double SX3S1 = 0, SX3S2 = -0.28988844, SX3X3 = 0, SX3X4 = 0, SX3X5 = 0, SX3X6 = 0, SX3X7 = -0.28988844;
    static const double SX4S1 = 0, SX4S2 = -0.28988844, SX4X4 = 0, SX4X5 = 0, SX4X6 = 0, SX4X7 = -0.28988844;
    static const double SX5S1 = 0, SX5S2 = -0.57977688, SX5X5 = 0, SX5X6 = 0, SX5X7 = -0.57977688;
    static const double SX6S1 = 0, SX6S2 = -0.28988844, SX6X6 = 0, SX6X7 = -0.28988844;
    static const double SX7X7 = 0.57977688, SX7S2 = 0.57977688;
    static const double SX2X2X7 = 0, SX2X2S2 = 0, SX2X3X7 = 0, SX2X3S2 = 0, SX2X4X7 = 0, SX2X4S2 = 0, SX2X5X7 = 0, SX2X5S2 = 0;
    static const double SX2X6X7 = 0, SX2X6S2 = 0, SX2X7X7 = 0, SX2X7S2 = 0, SX2S2S2 = 0;
    static const double SX3X3X7 = 0, SX3X3S2 = 0, SX3X4X7 = 0, SX3X4S2 = 0, SX3X5X7 = 0, SX3X5S2 = 0, SX3X6X7 = 0, SX3X6S2 = 0;
    static const double SX3X7X7 = 0, SX3X7S2 = 0, SX3S2S2 = 0;
    static const double SX4X4X7 = 0, SX4X4S2 = 0, SX4X5X7 = 0, SX4X5S2 = 0, SX4X6X7 = 0, SX4X6S2 = 0;
    static const double SX4X7X7 = 0, SX4X7S2 = 0, SX4S2S2 = 0;
    static const double SX5X5X7 = 0, SX5X5S2 = 0, SX5X6X7 = 0, SX5X6S2 = 0, SX5X7X7 = 0, SX5X7S2 = 0, SX5S2S2 = 0;
    static const double SX6X6X7 = 0, SX6X6S2 = 0, SX6X7X7 = 0, SX6X7S2 = 0, SX6S2S2 = 0;
    static const double SX7X7X7 = 0, SX7X7S2 = 0, SX7S2S2 = 0;
    static const double VX2 = 0.03710980000000141, VX2S1 = 0, VX2S2 = 0.142256;
    static const double VX2X2 = -0.014121, VX2X3 = -0.0070605, VX2X4 = -0.0070605, VX2X5 = -0.014121, VX2X6 = -0.0070605, VX2X7 = -0.01317580000000132;
    static const double VX3 = -0.126602, VX3S1 = 0, VX3S2 = 0.07267179999999968;
    static const double VX3X3 = 0, VX3X4 = 0, VX3X5 = 0, VX3X6 = 0, VX3X7 = -0.06121620000000034;
    static const double VX4 = -0.126602, VX4S1 = 0, VX4S2 = 0.07267179999999968;
    static const double VX4X4 = 0, VX4X5 = 0, VX4X6 = 0, VX4X7 = -0.06121620000000034;
    static const double VX5 = -0.126602, VX5S1 = 0, VX5S2 = 0.10141159999999935;
    static const double VX5X5 = 0, VX5X6 = 0, VX5X7 = -0.03247640000000067;
    static const double VX6 = -0.063301, VX6S1 = 0, VX6S2 = 0.050705799999999676;
    static const double VX6X6 = 0, VX6X7 = -0.016238200000000334;
    static const double VX7 = -0.14879160000000066, VX7S1 = 0, VX7S2 = 0.022016400000000658, VX7X7 = -0.01250159999999934;
    static const double VX2X2S2 = 0.014644, VX2X2X7 = -0.014644;
    static const double VX2X3S2 = -0.071128, VX2X3X7 = 0.029288;
    static const double VX2X4S2 = -0.071128, VX2X4X7 = 0.029288;
    static const double VX2X5S2 = -0.071128, VX2X5X7 = 0.029288;
    static const double VX2X6S2 = -0.035564, VX2X6X7 = 0.014644;
    static const double VX2X7S2 = -0.15690000000000004, VX2X7X7 = 0.080542, VX2S2S2 = -0.007322;
    static const double VX3X3S2 = -0.025104, VX3X3X7 = 0.025104;
    static const double VX3X4S2 = -0.044978000000000004, VX3X4X7 = 0.055438;
    static const double VX3X5S2 = -0.055438, VX3X5X7 = 0.044978000000000004;
    static const double VX3X6S2 = -0.025104, VX3X6X7 = 0.025104;
    static const double VX3X7S2 = -0.043932, VX3X7X7 = 0.025627, VX3S2S2 = 0.08106500000000001;
    static const double VX4X4S2 = -0.025104, VX4X4X7 = 0.025104;
    static const double VX4X5S2 = -0.055438, VX4X5X7 = 0.044978000000000004;
    static const double VX4X6S2 = -0.030334, VX4X6X7 = 0.019874000000000003;
    static const double VX4X7S2 = -0.043932, VX4X7X7 = 0.025627, VX4S2S2 = 0.08106500000000001;
    static const double VX5X5S2 = -0.035564, VX5X5X7 = 0.014644;
    static const double VX5X6S2 = -0.035564, VX5X6X7 = 0.014644;
    static const double VX5X7S2 = -0.064852, VX5X7X7 = 0.007322, VX5S2S2 = 0.07845000000000002;
    static const double VX6X6S2 = -0.010198500000000001, VX6X6X7 = 0.0023534999999999997;
    static const double VX6X7S2 = -0.032426, VX6X7X7 = 0.003661, VX6S2S2 = 0.03922500000000001;
    static const double VX7X7S2 = 0.012552000000000008, VX7X7X7 = -0.043932, VX7S2S2 = 0.15690000000000004;

    double dgdr0 = Rgas*T*(creal(clog(xfe2m1)) - creal(clog(xmg2m1))) +                ((HX2)-T*(SX2)+Pv*(VX2)) +                ((HX2X2)-T*(SX2X2)+Pv*(VX2X2))*r0*2.0 +                ((HX2X3)-T*(SX2X3)+Pv*(VX2X3))*r1 +                ((HX2X4)-T*(SX2X4)+Pv*(VX2X4))*r2 +                ((HX2X5)-T*(SX2X5)+Pv*(VX2X5))*r3 +                ((HX2X6)-T*(SX2X6)+Pv*(VX2X6))*r4 +                ((HX2X7)-T*(SX2X7)+Pv*(VX2X7))*r5 +                ((HX2S1)-T*(SX2S1)+Pv*(VX2S1))*s0 +                ((HX2S2)-T*(SX2S2)+Pv*(VX2S2))*s1 +                ((HX2X2X7)+Pv*(VX2X2X7))*r0*r5*2.0 +                ((HX2X2S2)+Pv*(VX2X2S2))*r0*s1*2.0 +                ((HX2X3X7)+Pv*(VX2X3X7))*r1*r5 +                ((HX2X3S2)+Pv*(VX2X3S2))*r1*s1 +                ((HX2X4X7)+Pv*(VX2X4X7))*r2*r5 +                ((HX2X4S2)+Pv*(VX2X4S2))*r2*s1 +                ((HX2X5X7)+Pv*(VX2X5X7))*r3*r5 +                ((HX2X5S2)+Pv*(VX2X5S2))*r3*s1 +                ((HX2X6X7)+Pv*(VX2X6X7))*r4*r5 +                ((HX2X6S2)+Pv*(VX2X6S2))*r4*s1 +                ((HX2X7X7)-T*(SX2X7X7)+Pv*(VX2X7X7))*r5*r5 +                ((HX2X7S2)+Pv*(VX2X7S2))*r5*s1 +                ((HX2S2S2)+Pv*(VX2S2S2))*s1*s1;
    double dgdr1 = Rgas*T*(0.5*creal(clog(xti4m1)) - 0.5*creal(clog(xmg2m1)) + 0.5*creal(clog(xm1_notTi))                  - creal(clog(xm1_MgFeNoTi)) + 0.5*creal(clog(xm1_MgFe))                  + 0.5*creal(clog(xrest_notNa)) - creal(clog(xtet_notSi))                  - 0.5*creal(clog(xm1_notMgFe)) + creal(clog(xal3tet))) +                ((HX3)-T*(SX3)+Pv*(VX3)) +                ((HX2X3)-T*(SX2X3)+Pv*(VX2X3))*r0 +                ((HX3X3)-T*(SX3X3)+Pv*(VX3X3))*r1*2.0 +                ((HX3X4)-T*(SX3X4)+Pv*(VX3X4))*r2 +                ((HX3X5)-T*(SX3X5)+Pv*(VX3X5))*r3 +                ((HX3X6)-T*(SX3X6)+Pv*(VX3X6))*r4 +                ((HX3X7)-T*(SX3X7)+Pv*(VX3X7))*r5 +                ((HX3S1)-T*(SX3S1)+Pv*(VX3S1))*s0 +                ((HX3S2)-T*(SX3S2)+Pv*(VX3S2))*s1 +                ((HX2X3X7)+Pv*(VX2X3X7))*r0*r5 +                ((HX2X3S2)+Pv*(VX2X3S2))*r0*s1 +                ((HX3X3X7)+Pv*(VX3X3X7))*r1*r5*2.0 +                ((HX3X3S2)+Pv*(VX3X3S2))*r1*s1*2.0 +                ((HX3X4X7)+Pv*(VX3X4X7))*r2*r5 +                ((HX3X4S2)+Pv*(VX3X4S2))*r2*s1 +                ((HX3X5X7)+Pv*(VX3X5X7))*r3*r5 +                ((HX3X5S2)+Pv*(VX3X5S2))*r3*s1 +                ((HX3X6X7)+Pv*(VX3X6X7))*r4*r5 +                ((HX3X6S2)+Pv*(VX3X6S2))*r4*s1 +                ((HX3X7X7)-T*(SX3X7X7)+Pv*(VX3X7X7))*r5*r5 +                ((HX3X7S2)-T*(SX3X7S2)+Pv*(VX3X7S2))*r5*s1 +                ((HX3S2S2)+Pv*(VX3S2S2))*s1*s1;
    double dgdr2 = Rgas*T*(0.5*creal(clog(xti4m1)) - 0.5*creal(clog(xmg2m1)) + 0.5*creal(clog(xm1_notTi))                  - creal(clog(xm1_MgFeNoTi)) + 0.5*creal(clog(xm1_MgFe))                  + 0.5*creal(clog(xrest_notNa)) - creal(clog(xtet_notSi))                  - 0.5*creal(clog(xm1_notMgFe)) + creal(clog(xfe3tet))) +                ((HX4)-T*(SX4)+Pv*(VX4)) +                ((HX2X4)-T*(SX2X4)+Pv*(VX2X4))*r0 +                ((HX3X4)-T*(SX3X4)+Pv*(VX3X4))*r1 +                ((HX4X4)-T*(SX4X4)+Pv*(VX4X4))*r2*2.0 +                ((HX4X5)-T*(SX4X5)+Pv*(VX4X5))*r3 +                ((HX4X6)-T*(SX4X6)+Pv*(VX4X6))*r4 +                ((HX4X7)-T*(SX4X7)+Pv*(VX4X7))*r5 +                ((HX4S1)-T*(SX4S1)+Pv*(VX4S1))*s0 +                ((HX4S2)-T*(SX4S2)+Pv*(VX4S2))*s1 +                ((HX2X4X7)+Pv*(VX2X4X7))*r0*r5 +                ((HX2X4S2)+Pv*(VX2X4S2))*r0*s1 +                ((HX3X4X7)+Pv*(VX3X4X7))*r1*r5 +                ((HX3X4S2)+Pv*(VX3X4S2))*r1*s1 +                ((HX4X4X7)+Pv*(VX4X4X7))*r2*r5*2.0 +                ((HX4X4S2)+Pv*(VX4X4S2))*r2*s1*2.0 +                ((HX4X5X7)+Pv*(VX4X5X7))*r3*r5 +                ((HX4X5S2)+Pv*(VX4X5S2))*r3*s1 +                ((HX4X6X7)+Pv*(VX4X6X7))*r4*r5 +                ((HX4X6S2)+Pv*(VX4X6S2))*r4*s1 +                ((HX4X7X7)-T*(SX4X7X7)+Pv*(VX4X7X7))*r5*r5 +                ((HX4X7S2)-T*(SX4X7S2)+Pv*(VX4X7S2))*r5*s1 +                ((HX4S2S2)+Pv*(VX4S2S2))*s1*s1;
    double dgdr3 = Rgas*T*(0.5*creal(clog(xal3m1)) + 0.5*creal(clog(xfe3m1)) - creal(clog(xmg2m1))                  - creal(clog(xm1_MgFeNoTi)) + creal(clog(xm1_MgFe))                  - creal(clog(xtet_notSi)) + 0.5*creal(clog(xal3tet)) + 0.5*creal(clog(xfe3tet))                  + creal(clog(xrest_notNa)) - creal(clog(xm1_notMgFe))) +                ((HX5)-T*(SX5)+Pv*(VX5)) +                ((HX2X5)-T*(SX2X5)+Pv*(VX2X5))*r0 +                ((HX3X5)-T*(SX3X5)+Pv*(VX3X5))*r1 +                ((HX4X5)-T*(SX4X5)+Pv*(VX4X5))*r2 +                ((HX5X5)-T*(SX5X5)+Pv*(VX5X5))*r3*2.0 +                ((HX5X6)-T*(SX5X6)+Pv*(VX5X6))*r4 +                ((HX5X7)-T*(SX5X7)+Pv*(VX5X7))*r5 +                ((HX5S1)-T*(SX5S1)+Pv*(VX5S1))*s0 +                ((HX5S2)-T*(SX5S2)+Pv*(VX5S2))*s1 +                ((HX2X5X7)+Pv*(VX2X5X7))*r0*r5 +                ((HX2X5S2)+Pv*(VX2X5S2))*r0*s1 +                ((HX3X5X7)+Pv*(VX3X5X7))*r1*r5 +                ((HX3X5S2)+Pv*(VX3X5S2))*r1*s1 +                ((HX4X5X7)+Pv*(VX4X5X7))*r2*r5 +                ((HX4X5S2)+Pv*(VX4X5S2))*r2*s1 +                ((HX5X5X7)+Pv*(VX5X5X7))*r3*r5*2.0 +                ((HX5X5S2)+Pv*(VX5X5S2))*r3*s1*2.0 +                ((HX5X6X7)+Pv*(VX5X6X7))*r4*r5 +                ((HX5X6S2)+Pv*(VX5X6S2))*r4*s1 +                ((HX5X7X7)-T*(SX5X7X7)+Pv*(VX5X7X7))*r5*r5 +                ((HX5X7S2)-T*(SX5X7S2)+Pv*(VX5X7S2))*r5*s1 +                ((HX5S2S2)+Pv*(VX5S2S2))*s1*s1;
    double dgdr4 = Rgas*T*(0.25*creal(clog(xal3m1)) + 0.25*creal(clog(xfe3m1)) - 0.5*creal(clog(xmg2m1))                  - 0.5*creal(clog(xm1_MgFeNoTi)) + 0.5*creal(clog(xm1_MgFe))                  + 0.5*creal(clog(xtet_notSi)) - 0.25*creal(clog(xal3tet))                  - 0.25*creal(clog(xfe3tet)) - creal(clog(xca2m2)) + creal(clog(xna1m2))                  - 0.5*creal(clog(xrest_notNa))                  - 0.5*creal(clog(xm1_notMgFe)) + creal(clog(xm2_notNa))) +                ((HX6)-T*(SX6)+Pv*(VX6)) +                ((HX2X6)-T*(SX2X6)+Pv*(VX2X6))*r0 +                ((HX3X6)-T*(SX3X6)+Pv*(VX3X6))*r1 +                ((HX4X6)-T*(SX4X6)+Pv*(VX4X6))*r2 +                ((HX5X6)-T*(SX5X6)+Pv*(VX5X6))*r3 +                ((HX6X6)-T*(SX6X6)+Pv*(VX6X6))*r4*2.0 +                ((HX6X7)-T*(SX6X7)+Pv*(VX6X7))*r5 +                ((HX6S1)-T*(SX6S1)+Pv*(VX6S1))*s0 +                ((HX6S2)-T*(SX6S2)+Pv*(VX6S2))*s1 +                ((HX2X6X7)+Pv*(VX2X6X7))*r0*r5 +                ((HX2X6S2)+Pv*(VX2X6S2))*r0*s1 +                ((HX3X6X7)+Pv*(VX3X6X7))*r1*r5 +                ((HX3X6S2)+Pv*(VX3X6S2))*r1*s1 +                ((HX4X6X7)+Pv*(VX4X6X7))*r2*r5 +                ((HX4X6S2)+Pv*(VX4X6S2))*r2*s1 +                ((HX5X6X7)+Pv*(VX5X6X7))*r3*r5 +                ((HX5X6S2)+Pv*(VX5X6S2))*r3*s1 +                ((HX6X6X7)+Pv*(VX6X6X7))*r4*r5*2.0 +                ((HX6X6S2)+Pv*(VX6X6S2))*r4*s1*2.0 +                ((HX6X7X7)-T*(SX6X7X7)+Pv*(VX6X7X7))*r5*r5 +                ((HX6X7S2)-T*(SX6X7S2)+Pv*(VX6X7S2))*r5*s1 +                ((HX6S2S2)+Pv*(VX6S2S2))*s1*s1;
    double dgdr5 = 0.5*Rgas*T*(creal(clog(xmg2m1)) - 2.0*creal(clog(xca2m2))                  + creal(clog(xmg2m2)) + creal(clog(xfe2m2/xfe2m1))) +                ((HX7)-T*(SX7)+Pv*(VX7)) +                ((HX2X7)-T*(SX2X7)+Pv*(VX2X7))*r0 +                ((HX3X7)-T*(SX3X7)+Pv*(VX3X7))*r1 +                ((HX4X7)-T*(SX4X7)+Pv*(VX4X7))*r2 +                ((HX5X7)-T*(SX5X7)+Pv*(VX5X7))*r3 +                ((HX6X7)-T*(SX6X7)+Pv*(VX6X7))*r4 +                ((HX7X7)-T*(SX7X7)+Pv*(VX7X7))*r5*2.0 +                ((HX7S1)-T*(SX7S1)+Pv*(VX7S1))*s0 +                ((HX7S2)-T*(SX7S2)+Pv*(VX7S2))*s1 +                ((HX2X2X7)+Pv*(VX2X2X7))*r0*r0 +                ((HX2X3X7)+Pv*(VX2X3X7))*r0*r1 +                ((HX2X4X7)+Pv*(VX2X4X7))*r0*r2 +                ((HX2X5X7)+Pv*(VX2X5X7))*r0*r3 +                ((HX2X6X7)+Pv*(VX2X6X7))*r0*r4 +                ((HX2X7X7)-T*(SX2X7X7)+Pv*(VX2X7X7))*r0*r5*2.0 +                ((HX2X7S2)+Pv*(VX2X7S2))*r0*s1 +                ((HX3X3X7)+Pv*(VX3X3X7))*r1*r1 +                ((HX3X4X7)+Pv*(VX3X4X7))*r1*r2 +                ((HX3X5X7)+Pv*(VX3X5X7))*r1*r3 +                ((HX3X6X7)+Pv*(VX3X6X7))*r1*r4 +                ((HX3X7X7)-T*(SX3X7X7)+Pv*(VX3X7X7))*r1*r5*2.0 +                ((HX3X7S2)-T*(SX3X7S2)+Pv*(VX3X7S2))*r1*s1 +                ((HX4X4X7)+Pv*(VX4X4X7))*r2*r2 +                ((HX4X5X7)+Pv*(VX4X5X7))*r2*r3 +                ((HX4X6X7)+Pv*(VX4X6X7))*r2*r4 +                ((HX4X7X7)-T*(SX4X7X7)+Pv*(VX4X7X7))*r2*r5*2.0 +                ((HX4X7S2)-T*(SX4X7S2)+Pv*(VX4X7S2))*r2*s1 +                ((HX5X5X7)+Pv*(VX5X5X7))*r3*r3 +                ((HX5X6X7)+Pv*(VX5X6X7))*r3*r4 +                ((HX5X7X7)-T*(SX5X7X7)+Pv*(VX5X7X7))*r3*r5*2.0 +                ((HX5X7S2)-T*(SX5X7S2)+Pv*(VX5X7S2))*r3*s1 +                ((HX6X6X7)+Pv*(VX6X6X7))*r4*r4 +                ((HX6X7X7)-T*(SX6X7X7)+Pv*(VX6X7X7))*r4*r5*2.0 +                ((HX6X7S2)-T*(SX6X7S2)+Pv*(VX6X7S2))*r4*s1 +                ((HX7X7X7)-T*(SX7X7X7)+Pv*(VX7X7X7))*r5*r5*3.0 +                ((HX7X7S2)-T*(SX7X7S2)+Pv*(VX7X7S2))*r5*s1*2.0 +                ((HX7S2S2)+Pv*(VX7S2S2))*s1*s1;

    /* ends[]: real orthopyroxene.c's own "purePyx"/"pureOrder" functions
       (computing the ENDMEMBERS reference values) hardcode
       "static const int clino = TRUE;" INTERNALLY, regardless of the
       surrounding opx (clino=FALSE) mixing-model context - confirmed by
       direct source reading (orthopyroxene.c lines 1385, 1470). This
       means real opx's endmember reference corrections use CPX's own
       Taylor coefficients, not opx's - so DI_G=EN_G=HD_G=CA_G=CF_G=JD_G=0
       exactly for opx too (same as cpx), and only essenite is nonzero,
       using the SAME GH_cpx_pure_ES_G already used by obj_gh_cpx (no
       separate opx version needed). An earlier version of this function
       used a from-scratch 7-endmember "GH_opx_ends" (opx Taylor
       coefficients) reasoning the whole ENDMEMBERS calc should follow the
       overall clino=FALSE context - that was wrong: it's a real, easy-to-
       miss exception in the source (confirmed the hard way, via a
       ~10 kJ residual against real MELTS at typical/pure compositions
       that didn't resolve until this was found). GH_opx_ends is kept
       (unused) as a documented dead end, not deleted, in case a future
       session needs the derivation again. */
    double es_end = GH_cpx_pure_ES_G(T, Rgas, Pv);
    double Emix = p5*es_end;
    double Gex = Graw - Emix;

    double dGex[7];
    dGex[0] = 0.0;
    dGex[1] = dgdr5;
    dGex[2] = dgdr0;
    dGex[3] = dgdr1;
    dGex[4] = dgdr2;
    dGex[5] = dgdr3 - es_end;
    dGex[6] = 0.5*dgdr1 - 0.5*dgdr2 + 0.5*dgdr3 + dgdr4;

    for (int i = 0; i < 7; i++){
        double mu = Gex;
        for (int k = 0; k < 7; k++){
            double delta_ik = (i==k) ? 1.0 : 0.0;
            mu += (delta_ik - p[k])*dGex[k];
        }
        mu_Gex[i] = mu/1000.0;
        d->sf[i]  = p[i];
    }

    d->sum_apep = 0.0;
    for (int i = 0; i < 7; i++){ d->sum_apep += d->ape[i]*p[i]; }
    d->factor = d->fbc/d->sum_apep;

    d->df_raw = 0.0;
    for (int i = 0; i < 7; i++){ d->df_raw += (mu_Gex[i] + gb[i])*p[i]; }
    d->df = d->df_raw * d->factor;

    if (grad){
        for (int i = 0; i < 7; i++){
            grad[i] = (mu_Gex[i] + gb[i])*d->factor - (d->df_raw*d->factor*(d->ape[i]/d->sum_apep));
        }
    }
    return d->df;
}

void GH_SS_objective_init_function(    obj_type            *SS_objective,
                                        global_variable      gv                  ){
    GH_spn_multistart_flag = gv.gh_multistart_order;
    GH_cpx_SS2 = (gv.EM_database == 0) ? -1.08018328 : 0.0;
    GH_opx_SS2 = (gv.EM_database == 0) ? -0.57977688 : 0.0;
    GH_actual_EM_database = gv.EM_database;
    for (int iss = 0; iss < gv.len_ss; iss++){
        if (strcmp( gv.SS_list[iss], "liq") == 0 ){
            SS_objective[iss] = obj_gh_liq;
        }
        else if (strcmp( gv.SS_list[iss], "ol") == 0 ){
            SS_objective[iss] = obj_gh_ol;
        }
        else if (strcmp( gv.SS_list[iss], "fsp") == 0 ){
            SS_objective[iss] = obj_gh_fsp;
        }
        else if (strcmp( gv.SS_list[iss], "bi") == 0 ){
            SS_objective[iss] = obj_gh_bi;
        }
        else if (strcmp( gv.SS_list[iss], "g") == 0 ){
            SS_objective[iss] = obj_gh_g;
        }
        else if (strcmp( gv.SS_list[iss], "hb") == 0 ){
            SS_objective[iss] = obj_gh_hb;
        }
        else if (strcmp( gv.SS_list[iss], "lc") == 0 ){
            SS_objective[iss] = obj_gh_lc;
        }
        else if (strcmp( gv.SS_list[iss], "mel") == 0 ){
            SS_objective[iss] = obj_gh_mel;
        }
        else if (strcmp( gv.SS_list[iss], "cum") == 0 ){
            SS_objective[iss] = obj_gh_cum;
        }
        else if (strcmp( gv.SS_list[iss], "spn") == 0 ){
            SS_objective[iss] = obj_gh_spn;
        }
        else if (strcmp( gv.SS_list[iss], "cpx") == 0 ){
            SS_objective[iss] = obj_gh_cpx;
        }
        else if (strcmp( gv.SS_list[iss], "opx") == 0 ){
            SS_objective[iss] = obj_gh_opx;
        }
        else if (strcmp( gv.SS_list[iss], "fl") == 0 ){
            SS_objective[iss] = obj_gh_fluid;
        }
        else if (strcmp( gv.SS_list[iss], "rhm") == 0 ){
            SS_objective[iss] = obj_gh_rhm;
        }
        else if (strcmp( gv.SS_list[iss], "nph") == 0 ){
            SS_objective[iss] = obj_gh_nph;
        }
        else if (strcmp( gv.SS_list[iss], "kls") == 0 ){
            SS_objective[iss] = obj_gh_kls;
        }
    }
}

void GH_PC_init(                       PC_type             *PC_read,
                                        global_variable      gv                  ){
    for (int iss = 0; iss < gv.len_ss; iss++){
        if (strcmp( gv.SS_list[iss], "liq") == 0 ){
            PC_read[iss] = obj_gh_liq;
        }
        else if (strcmp( gv.SS_list[iss], "ol") == 0 ){
            PC_read[iss] = obj_gh_ol;
        }
        else if (strcmp( gv.SS_list[iss], "fsp") == 0 ){
            PC_read[iss] = obj_gh_fsp;
        }
        else if (strcmp( gv.SS_list[iss], "bi") == 0 ){
            PC_read[iss] = obj_gh_bi;
        }
        else if (strcmp( gv.SS_list[iss], "g") == 0 ){
            PC_read[iss] = obj_gh_g;
        }
        else if (strcmp( gv.SS_list[iss], "hb") == 0 ){
            PC_read[iss] = obj_gh_hb;
        }
        else if (strcmp( gv.SS_list[iss], "lc") == 0 ){
            PC_read[iss] = obj_gh_lc;
        }
        else if (strcmp( gv.SS_list[iss], "mel") == 0 ){
            PC_read[iss] = obj_gh_mel;
        }
        else if (strcmp( gv.SS_list[iss], "cum") == 0 ){
            PC_read[iss] = obj_gh_cum;
        }
        else if (strcmp( gv.SS_list[iss], "spn") == 0 ){
            PC_read[iss] = obj_gh_spn;
        }
        else if (strcmp( gv.SS_list[iss], "cpx") == 0 ){
            PC_read[iss] = obj_gh_cpx;
        }
        else if (strcmp( gv.SS_list[iss], "opx") == 0 ){
            PC_read[iss] = obj_gh_opx;
        }
        else if (strcmp( gv.SS_list[iss], "fl") == 0 ){
            PC_read[iss] = obj_gh_fluid;
        }
        else if (strcmp( gv.SS_list[iss], "rhm") == 0 ){
            PC_read[iss] = obj_gh_rhm;
        }
        else if (strcmp( gv.SS_list[iss], "nph") == 0 ){
            PC_read[iss] = obj_gh_nph;
        }
        else if (strcmp( gv.SS_list[iss], "kls") == 0 ){
            PC_read[iss] = obj_gh_kls;
        }
    }
}

/**
    p-to-xeos map for the "pure endmember as LP-swap candidate" mechanism
    (swap_pure_endmembers, simplex_levelling.c): unlike tc's per-phase
    algebra (e.g. p2x_ig_fsp drops a reference endmember to build a reduced
    basis), every gh phase is direct p=x with no basis reduction, so the
    map is just identity + clamp to bounds. One shared function suffices
    for all 5 phases (liq/ol/fsp/bi/g).
*/
void p2x_gh_generic(void *SS_ref_db, double eps){
    SS_ref *d = (SS_ref *) SS_ref_db;
    for (int i = 0; i < d->n_xeos; i++){
        d->iguess[i] = d->p[i];
        if (d->iguess[i] < d->bounds[i][0]){ d->iguess[i] = d->bounds[i][0]; }
        if (d->iguess[i] > d->bounds[i][1]){ d->iguess[i] = d->bounds[i][1]; }
    }
}

void GH_P2X_init(                      P2X_type            *P2X_read,
                                        global_variable      gv                  ){
    for (int iss = 0; iss < gv.len_ss; iss++){
        if (strcmp( gv.SS_list[iss], "liq") == 0 ||
            strcmp( gv.SS_list[iss], "ol")  == 0 ||
            strcmp( gv.SS_list[iss], "fsp") == 0 ||
            strcmp( gv.SS_list[iss], "bi")  == 0 ||
            strcmp( gv.SS_list[iss], "g")   == 0 ||
            strcmp( gv.SS_list[iss], "hb")  == 0 ||
            strcmp( gv.SS_list[iss], "lc")  == 0 ||
            strcmp( gv.SS_list[iss], "mel") == 0 ||
            strcmp( gv.SS_list[iss], "cum") == 0 ||
            strcmp( gv.SS_list[iss], "spn") == 0 ||
            strcmp( gv.SS_list[iss], "cpx") == 0 ||
            strcmp( gv.SS_list[iss], "opx") == 0 ||
            strcmp( gv.SS_list[iss], "fl") == 0 ||
            strcmp( gv.SS_list[iss], "rhm") == 0 ||
            strcmp( gv.SS_list[iss], "nph") == 0 ||
            strcmp( gv.SS_list[iss], "kls") == 0){
            P2X_read[iss] = p2x_gh_generic;
        }
    }
}

/** Evaluate a candidate pseudocompound's G and reconstruct its oxide composition
    (mirrors SB_PC_function's structure). */
SS_ref GH_PC_function(     global_variable      gv,
                            PC_type             *PC_read,
                            SS_ref               SS_ref_db,
                            bulk_info            z_b,
                            int                  ph_id                  ){

    double G0 = (*PC_read[ph_id])( SS_ref_db.n_xeos,
                                    SS_ref_db.iguess,
                                    SS_ref_db.dfx,
                                    &SS_ref_db                     );
    SS_ref_db.df = G0;

    for (int j = 0; j < gv.len_ox; j++){
        SS_ref_db.ss_comp[j] = 0.0;
    }
    for (int i = 0; i < SS_ref_db.n_em; i++){
        for (int j = 0; j < gv.len_ox; j++){
            SS_ref_db.ss_comp[j] += SS_ref_db.Comp[i][j]*SS_ref_db.p[i]*SS_ref_db.z_em[i];
        }
    }

    SS_ref_db.sf_ok = 1;

    return SS_ref_db;
}
