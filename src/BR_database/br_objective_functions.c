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
/* Source: Pourteau et al. (2014), Contrib Mineral Petrol -- Ctd/Car are
   both "(IDEAL) M(1):Fe,Mg", W=0, so this is a plain ideal 1-site binary. */
#include <complex.h>
#include <string.h>

#include "br_objective_functions.h"

static double obj_br_ideal_binary(unsigned n, const double *x, double *grad, void *SS_ref_db){
    (void) n;
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

double obj_br_ctd(unsigned n, const double *x, double *grad, void *SS_ref_db){
    return obj_br_ideal_binary(n, x, grad, SS_ref_db);
}

double obj_br_car(unsigned n, const double *x, double *grad, void *SS_ref_db){
    return obj_br_ideal_binary(n, x, grad, SS_ref_db);
}

/* Chl: (-SITE,MARGULES) M4(1):Al,Fe,Mg - M23(4):Fe,Mg,Al - M1(1):Fe,Mg,Al,v - T1(2):Si,Al.
   M4 is Al in all 5 endmembers (invariant, no entropy/energy contribution).
   M23 and T1 site fractions are linear functions of the 5 endmember
   fractions (Tschermak-coupled to M1, not independent x-eos). Only M1
   carries Margules non-ideality (SITEMARG table); M23/T1 mix ideally.
   Endmember order: dph, clin, feam, sud, ames. */
double obj_br_chl(unsigned n, const double *x, double *grad, void *SS_ref_db){
    (void) n;
    SS_ref *d = (SS_ref *) SS_ref_db;

    int n_em = d->n_em;
    double R = d->R;
    double T = d->T;
    double *p  = d->p;
    double *gb = d->gb_lvl;
    double *W  = d->W;
    double eps = 1e-10;

    for (int i = 0; i < n_em; i++){ p[i] = x[i]; }

    double x23_Fe = p[0] + p[2];
    double x23_Mg = p[1] + p[4] + 0.5*p[3];
    double x23_Al = 0.5*p[3];

    double x1_Fe = p[0];
    double x1_Mg = p[1];
    double x1_Al = p[2] + p[4];
    double x1_v  = p[3];

    double xT1_Al = 0.5*p[0] + 0.5*p[1] + p[2] + 0.5*p[3] + p[4];
    double xT1_Si = 1.0 - xT1_Al;

    d->sf[0] = x23_Fe; d->sf[1] = x23_Mg; d->sf[2] = x23_Al;
    d->sf[3] = x1_Fe;  d->sf[4] = x1_Mg;  d->sf[5] = x1_Al; d->sf[6] = x1_v;
    d->sf[7] = xT1_Al; d->sf[8] = xT1_Si;

    double S_terms =
        4.0*( x23_Fe*creal(clog(x23_Fe+eps)) + x23_Mg*creal(clog(x23_Mg+eps)) + x23_Al*creal(clog(x23_Al+eps)) ) +
        1.0*( x1_Fe*creal(clog(x1_Fe+eps)) + x1_Mg*creal(clog(x1_Mg+eps)) + x1_Al*creal(clog(x1_Al+eps)) + x1_v*creal(clog(x1_v+eps)) ) +
        2.0*( xT1_Al*creal(clog(xT1_Al+eps)) + xT1_Si*creal(clog(xT1_Si+eps)) );

    /* W order: FeMg, FeAl, Fev, MgAl, Mgv, Alv */
    double G_xs_J = W[0]*x1_Fe*x1_Mg + W[1]*x1_Fe*x1_Al + W[2]*x1_Fe*x1_v
                  + W[3]*x1_Mg*x1_Al + W[4]*x1_Mg*x1_v  + W[5]*x1_Al*x1_v;

    d->sum_apep = 0.0;
    for (int i = 0; i < n_em; i++){ d->sum_apep += d->ape[i]*p[i]; }
    d->factor = d->fbc/d->sum_apep;

    d->df_raw = 0.0;
    for (int i = 0; i < n_em; i++){ d->df_raw += gb[i]*p[i]; }
    d->df_raw += R*T*S_terms;
    d->df_raw += G_xs_J/1000.0;
    d->df = d->df_raw * d->factor;

    if (grad){
        /* dS/d(site fraction), each = ln(x+eps)+1 */
        double dS_23Fe = creal(clog(x23_Fe+eps))+1.0;
        double dS_23Mg = creal(clog(x23_Mg+eps))+1.0;
        double dS_23Al = creal(clog(x23_Al+eps))+1.0;
        double dS_1Fe  = creal(clog(x1_Fe+eps))+1.0;
        double dS_1Mg  = creal(clog(x1_Mg+eps))+1.0;
        double dS_1Al  = creal(clog(x1_Al+eps))+1.0;
        double dS_1v   = creal(clog(x1_v+eps))+1.0;
        double dS_T1Al = creal(clog(xT1_Al+eps))+1.0;
        double dS_T1Si = creal(clog(xT1_Si+eps))+1.0;

        /* dGxs/d(M1 site fraction) */
        double dG_1Fe = W[0]*x1_Mg + W[1]*x1_Al + W[2]*x1_v;
        double dG_1Mg = W[0]*x1_Fe + W[3]*x1_Al + W[4]*x1_v;
        double dG_1Al = W[1]*x1_Fe + W[3]*x1_Mg + W[5]*x1_v;
        double dG_1v  = W[2]*x1_Fe + W[4]*x1_Mg + W[5]*x1_Al;

        /* Jacobian d(site fraction)/dp[k], k = dph,clin,feam,sud,ames */
        static const double d23Fe_dp[5] = {1,0,1,0,0};
        static const double d23Mg_dp[5] = {0,1,0,0.5,1};
        static const double d23Al_dp[5] = {0,0,0,0.5,0};
        static const double d1Fe_dp[5]  = {1,0,0,0,0};
        static const double d1Mg_dp[5]  = {0,1,0,0,0};
        static const double d1Al_dp[5]  = {0,0,1,0,1};
        static const double d1v_dp[5]   = {0,0,0,1,0};
        static const double dT1Al_dp[5] = {0.5,0.5,1,0.5,1};

        for (int k = 0; k < n_em; k++){
            double dS =
                4.0*( dS_23Fe*d23Fe_dp[k] + dS_23Mg*d23Mg_dp[k] + dS_23Al*d23Al_dp[k] ) +
                1.0*( dS_1Fe*d1Fe_dp[k] + dS_1Mg*d1Mg_dp[k] + dS_1Al*d1Al_dp[k] + dS_1v*d1v_dp[k] ) +
                2.0*( dS_T1Al - dS_T1Si )*dT1Al_dp[k];

            double dG =
                dG_1Fe*d1Fe_dp[k] + dG_1Mg*d1Mg_dp[k] + dG_1Al*d1Al_dp[k] + dG_1v*d1v_dp[k];

            grad[k] = (gb[k] + R*T*dS + dG/1000.0)*d->factor
                    - (d->df_raw*d->factor*(d->ape[k]/d->sum_apep));
        }
    }

    return d->df;
}

/* MicaDubacq: (-SITE,MARGULES) I(1):Na,K,v,h - M1(1):v - M2(2):Fe,Mg,Al - T1(2):Si,Al.
   `h` is never occupied by any of the 5 real endmembers (confirmed from
   Pourteau's own endmember table) and M1 is `v` in all 5 - both invariant,
   no entropy/energy contribution. I-site uses Pourteau's doubled-species
   (subregular, Redlich-Kister-style) asymmetric Margules convention - NOT
   van Laar - confirmed by the exact term count (6 binary + 1 ternary = 7,
   matching a complete 3-component asymmetric expansion) and by TC's own
   feldspar model needing a *per-endmember* alpha array that this table
   simply doesn't provide. M2 uses the same simple symmetric form as Chl's
   M1 (single line per pair, WH/WS/WV). Endmember order: mu, pa, cel, fcel,
   prlph. */
double obj_br_mica(unsigned n, const double *x, double *grad, void *SS_ref_db){
    (void) n;
    SS_ref *d = (SS_ref *) SS_ref_db;

    int n_em = d->n_em;
    double R = d->R;
    double T = d->T;
    double *p  = d->p;
    double *gb = d->gb_lvl;
    double *W  = d->W;
    double eps = 1e-10;

    for (int i = 0; i < n_em; i++){ p[i] = x[i]; }

    double x_K  = p[0] + p[2] + p[3];
    double x_Na = p[1];
    double x_v  = p[4];

    double x_M2_Al = p[0] + p[1] + 0.5*p[2] + 0.5*p[3] + p[4];
    double x_M2_Mg = 0.5*p[2];
    double x_M2_Fe = 0.5*p[3];

    double x_T1_Al = 0.5*p[0] + 0.5*p[1];
    double x_T1_Si = 1.0 - x_T1_Al;

    d->sf[0] = x_K;  d->sf[1] = x_Na;    d->sf[2] = x_v;
    d->sf[3] = x_M2_Al; d->sf[4] = x_M2_Mg; d->sf[5] = x_M2_Fe;
    d->sf[6] = x_T1_Al; d->sf[7] = x_T1_Si;

    double S_terms =
        1.0*( x_K*creal(clog(x_K+eps)) + x_Na*creal(clog(x_Na+eps)) + x_v*creal(clog(x_v+eps)) ) +
        2.0*( x_M2_Al*creal(clog(x_M2_Al+eps)) + x_M2_Mg*creal(clog(x_M2_Mg+eps)) + x_M2_Fe*creal(clog(x_M2_Fe+eps)) ) +
        2.0*( x_T1_Al*creal(clog(x_T1_Al+eps)) + x_T1_Si*creal(clog(x_T1_Si+eps)) );

    /* I-site subregular asymmetric term: Gex_ij = xi*xj*(xi*Wiij + xj*Wijj) */
    double G_xs_I_J =
        x_K*x_Na*(x_K*W[0] + x_Na*W[1])
      + x_Na*x_v*(x_Na*W[2] + x_v*W[3])
      + x_K*x_v*(x_K*W[4] + x_v*W[5])
      + x_K*x_Na*x_v*W[6];

    /* M2-site symmetric term */
    double G_xs_M2_J = W[7]*x_M2_Mg*x_M2_Al + W[8]*x_M2_Fe*x_M2_Al;

    d->sum_apep = 0.0;
    for (int i = 0; i < n_em; i++){ d->sum_apep += d->ape[i]*p[i]; }
    d->factor = d->fbc/d->sum_apep;

    d->df_raw = 0.0;
    for (int i = 0; i < n_em; i++){ d->df_raw += gb[i]*p[i]; }
    d->df_raw += R*T*S_terms;
    d->df_raw += (G_xs_I_J + G_xs_M2_J)/1000.0;
    d->df = d->df_raw * d->factor;

    if (grad){
        double dS_K  = creal(clog(x_K+eps))+1.0;
        double dS_Na = creal(clog(x_Na+eps))+1.0;
        double dS_v  = creal(clog(x_v+eps))+1.0;
        double dS_M2Al = creal(clog(x_M2_Al+eps))+1.0;
        double dS_M2Mg = creal(clog(x_M2_Mg+eps))+1.0;
        double dS_M2Fe = creal(clog(x_M2_Fe+eps))+1.0;
        double dS_T1Al = creal(clog(x_T1_Al+eps))+1.0;
        double dS_T1Si = creal(clog(x_T1_Si+eps))+1.0;

        /* dGxs_I/d(I site fraction), from Gex_I = W0 K^2 Na + W1 K Na^2 + W2 Na^2 v
           + W3 Na v^2 + W4 K^2 v + W5 K v^2 + W6 K Na v */
        double dGI_K  = 2.0*W[0]*x_K*x_Na + W[1]*x_Na*x_Na + 2.0*W[4]*x_K*x_v + W[5]*x_v*x_v + W[6]*x_Na*x_v;
        double dGI_Na = W[0]*x_K*x_K + 2.0*W[1]*x_K*x_Na + 2.0*W[2]*x_Na*x_v + W[3]*x_v*x_v + W[6]*x_K*x_v;
        double dGI_v  = W[2]*x_Na*x_Na + 2.0*W[3]*x_Na*x_v + W[4]*x_K*x_K + 2.0*W[5]*x_K*x_v + W[6]*x_K*x_Na;

        /* dGxs_M2/d(M2 site fraction) */
        double dGM2_Al = W[7]*x_M2_Mg + W[8]*x_M2_Fe;
        double dGM2_Mg = W[7]*x_M2_Al;
        double dGM2_Fe = W[8]*x_M2_Al;

        /* Jacobian d(site fraction)/dp[k], k = mu,pa,cel,fcel,prlph */
        static const double dK_dp[5]    = {1,0,1,1,0};
        static const double dNa_dp[5]   = {0,1,0,0,0};
        static const double dv_dp[5]    = {0,0,0,0,1};
        static const double dM2Al_dp[5] = {1,1,0.5,0.5,1};
        static const double dM2Mg_dp[5] = {0,0,0.5,0,0};
        static const double dM2Fe_dp[5] = {0,0,0,0.5,0};
        static const double dT1Al_dp[5] = {0.5,0.5,0,0,0};

        for (int k = 0; k < n_em; k++){
            double dS =
                1.0*( dS_K*dK_dp[k] + dS_Na*dNa_dp[k] + dS_v*dv_dp[k] ) +
                2.0*( dS_M2Al*dM2Al_dp[k] + dS_M2Mg*dM2Mg_dp[k] + dS_M2Fe*dM2Fe_dp[k] ) +
                2.0*( dS_T1Al - dS_T1Si )*dT1Al_dp[k];

            double dG =
                dGI_K*dK_dp[k] + dGI_Na*dNa_dp[k] + dGI_v*dv_dp[k] +
                dGM2_Al*dM2Al_dp[k] + dGM2_Mg*dM2Mg_dp[k] + dGM2_Fe*dM2Fe_dp[k];

            grad[k] = (gb[k] + R*T*dS + dG/1000.0)*d->factor
                    - (d->df_raw*d->factor*(d->ape[k]/d->sum_apep));
        }
    }

    return d->df;
}

double obj_br_ilm(unsigned n, const double *x, double *grad, void *SS_ref_db){
    return obj_br_ideal_binary(n, x, grad, SS_ref_db);
}

double obj_br_ep(unsigned n, const double *x, double *grad, void *SS_ref_db){
    return obj_br_ideal_binary(n, x, grad, SS_ref_db);
}

double obj_br_spl(unsigned n, const double *x, double *grad, void *SS_ref_db){
    return obj_br_ideal_binary(n, x, grad, SS_ref_db);
}

/* Same 2-endmember, single-site regular-solution formalism as
   obj_br_ideal_binary, generalized to an arbitrary site multiplicity
   (TALC/BIOTITE M(3), OLIVINE/OPX M(2), AMPH M(4), STAU T(4)). */
static double obj_br_ideal_binary_mult(unsigned n, const double *x, double *grad, void *SS_ref_db, double mult){
    (void) n;
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

    double Sconfig = R*T*mult*(p[0]*creal(clog(p[0] + d->d_em[0])) + p[1]*creal(clog(p[1] + d->d_em[1])));

    d->df_raw = 0.0;
    for (int i = 0; i < n_em; i++){
        d->df_raw += (mu_Gex[i] + gb[i])*p[i];
    }
    d->df_raw += Sconfig;
    d->df = d->df_raw * d->factor;

    if (grad){
        double dS0 = R*T*mult*(creal(clog(p[0] + d->d_em[0]))+1.0);
        double dS1 = R*T*mult*(creal(clog(p[1] + d->d_em[1]))+1.0);
        grad[0] = (dS0 + mu_Gex[0] + gb[0])*d->factor - (d->df_raw*d->factor*(d->ape[0]/d->sum_apep));
        grad[1] = (dS1 + mu_Gex[1] + gb[1])*d->factor - (d->df_raw*d->factor*(d->ape[1]/d->sum_apep));
    }
    return d->df;
}

double obj_br_talc(unsigned n, const double *x, double *grad, void *SS_ref_db){
    return obj_br_ideal_binary_mult(n, x, grad, SS_ref_db, 3.0);
}

double obj_br_bt(unsigned n, const double *x, double *grad, void *SS_ref_db){
    return obj_br_ideal_binary_mult(n, x, grad, SS_ref_db, 3.0);
}

double obj_br_ol(unsigned n, const double *x, double *grad, void *SS_ref_db){
    return obj_br_ideal_binary_mult(n, x, grad, SS_ref_db, 2.0);
}

double obj_br_opx(unsigned n, const double *x, double *grad, void *SS_ref_db){
    return obj_br_ideal_binary_mult(n, x, grad, SS_ref_db, 2.0);
}

double obj_br_crd(unsigned n, const double *x, double *grad, void *SS_ref_db){
    return obj_br_ideal_binary_mult(n, x, grad, SS_ref_db, 2.0);
}

double obj_br_amph(unsigned n, const double *x, double *grad, void *SS_ref_db){
    return obj_br_ideal_binary_mult(n, x, grad, SS_ref_db, 4.0);
}

double obj_br_stau(unsigned n, const double *x, double *grad, void *SS_ref_db){
    return obj_br_ideal_binary_mult(n, x, grad, SS_ref_db, 4.0);
}

/* Garnet: (MARGULES,SITE) M(3):Mg,Fe, asymmetric subregular binary
   (py, alm), same Gex_ij = xi*xj*(xi*Wiij + xj*Wijj) form as one pair of
   MicaDubacq's I-site, standalone. */
double obj_br_grt(unsigned n, const double *x, double *grad, void *SS_ref_db){
    (void) n;
    SS_ref *d = (SS_ref *) SS_ref_db;

    double R = d->R;
    double T = d->T;
    double *p  = d->p;
    double *gb = d->gb_lvl;
    double *W  = d->W;
    double eps = 1e-10;
    double mult = 3.0;

    p[0] = x[0]; p[1] = x[1];
    double x1 = p[0], x2 = p[1];

    d->sf[0] = x1; d->sf[1] = x2;

    double S_terms = mult*( x1*creal(clog(x1+eps)) + x2*creal(clog(x2+eps)) );

    double G_xs = x1*x2*(x1*W[0] + x2*W[1]);

    d->sum_apep = 0.0;
    for (int i = 0; i < d->n_em; i++){ d->sum_apep += d->ape[i]*p[i]; }
    d->factor = d->fbc/d->sum_apep;

    d->df_raw = gb[0]*p[0] + gb[1]*p[1] + R*T*S_terms + G_xs/1000.0;
    d->df = d->df_raw * d->factor;

    if (grad){
        double dS1 = mult*(creal(clog(x1+eps))+1.0);
        double dS2 = mult*(creal(clog(x2+eps))+1.0);
        double dG1 = 2.0*W[0]*x1*x2 + W[1]*x2*x2;
        double dG2 = W[0]*x1*x1 + 2.0*W[1]*x1*x2;

        grad[0] = (gb[0] + R*T*dS1 + dG1/1000.0)*d->factor - (d->df_raw*d->factor*(d->ape[0]/d->sum_apep));
        grad[1] = (gb[1] + R*T*dS2 + dG2/1000.0)*d->factor - (d->df_raw*d->factor*(d->ape[1]/d->sum_apep));
    }

    return d->df;
}

/* Omphacite: (MARGULES,IDEAL), plain ideal ternary entropy (no site
   splitting) + 3-pair subregular Margules, no ternary correction term.
   Endmember order: di, jd, hed. */
double obj_br_omph(unsigned n, const double *x, double *grad, void *SS_ref_db){
    (void) n;
    SS_ref *d = (SS_ref *) SS_ref_db;

    double R = d->R;
    double T = d->T;
    double *p  = d->p;
    double *gb = d->gb_lvl;
    double *W  = d->W;
    double eps = 1e-10;

    for (int i = 0; i < 3; i++){ p[i] = x[i]; }
    double b = p[0]; /* di */
    double a = p[1]; /* jd */
    double c = p[2]; /* hed */

    d->sf[0] = b; d->sf[1] = a; d->sf[2] = c;

    double S_terms = ( a*creal(clog(a+eps)) + b*creal(clog(b+eps)) + c*creal(clog(c+eps)) );

    /* W[0..1]=Jd-Di (Wjj*d,Wj*dd), W[2..3]=Di-Hd, W[4..5]=Jd-Hd */
    double G_xs = W[0]*a*a*b + W[1]*a*b*b
                + W[2]*b*b*c + W[3]*b*c*c
                + W[4]*a*a*c + W[5]*a*c*c;

    d->sum_apep = 0.0;
    for (int i = 0; i < 3; i++){ d->sum_apep += d->ape[i]*p[i]; }
    d->factor = d->fbc/d->sum_apep;

    d->df_raw = gb[0]*p[0] + gb[1]*p[1] + gb[2]*p[2] + R*T*S_terms + G_xs/1000.0;
    d->df = d->df_raw * d->factor;

    if (grad){
        double dSa = creal(clog(a+eps))+1.0;
        double dSb = creal(clog(b+eps))+1.0;
        double dSc = creal(clog(c+eps))+1.0;

        double dGb = W[0]*a*a + 2.0*W[1]*a*b + 2.0*W[2]*b*c + W[3]*c*c;
        double dGa = 2.0*W[0]*a*b + W[1]*b*b + 2.0*W[4]*a*c + W[5]*c*c;
        double dGc = W[2]*b*b + 2.0*W[3]*b*c + W[4]*a*a + 2.0*W[5]*a*c;

        grad[0] = (gb[0] + R*T*dSb + dGb/1000.0)*d->factor - (d->df_raw*d->factor*(d->ape[0]/d->sum_apep));
        grad[1] = (gb[1] + R*T*dSa + dGa/1000.0)*d->factor - (d->df_raw*d->factor*(d->ape[1]/d->sum_apep));
        grad[2] = (gb[2] + R*T*dSc + dGc/1000.0)*d->factor - (d->df_raw*d->factor*(d->ape[2]/d->sum_apep));
    }

    return d->df;
}

/* AMPHx: (IDEAL,MARGULES) A(1):v,Na - M1(2):Mg,Al, reciprocal 2-site
   entropy + symmetric ternary Margules on raw endmember fractions
   (unlike Chl/Mica, whose Margules acts on site fractions).
   Endmember order: tr, tsch, parg. */
double obj_br_amphx(unsigned n, const double *x, double *grad, void *SS_ref_db){
    (void) n;
    SS_ref *d = (SS_ref *) SS_ref_db;

    double R = d->R;
    double T = d->T;
    double *p  = d->p;
    double *gb = d->gb_lvl;
    double *W  = d->W;
    double eps = 1e-10;

    for (int i = 0; i < 3; i++){ p[i] = x[i]; }
    double p_tr = p[0], p_tsch = p[1], p_parg = p[2];

    double x_Av  = p_tr + p_tsch;
    double x_ANa = p_parg;
    double x_M1Mg = p_tr + 0.5*p_parg;
    double x_M1Al = p_tsch + 0.5*p_parg;

    d->sf[0] = x_Av; d->sf[1] = x_ANa; d->sf[2] = x_M1Mg; d->sf[3] = x_M1Al;

    double S_terms =
        1.0*( x_Av*creal(clog(x_Av+eps)) + x_ANa*creal(clog(x_ANa+eps)) ) +
        2.0*( x_M1Mg*creal(clog(x_M1Mg+eps)) + x_M1Al*creal(clog(x_M1Al+eps)) );

    double G_xs = W[0]*p_tr*p_tsch + W[1]*p_tr*p_parg + W[2]*p_tsch*p_parg;

    d->sum_apep = 0.0;
    for (int i = 0; i < 3; i++){ d->sum_apep += d->ape[i]*p[i]; }
    d->factor = d->fbc/d->sum_apep;

    d->df_raw = gb[0]*p[0] + gb[1]*p[1] + gb[2]*p[2] + R*T*S_terms + G_xs/1000.0;
    d->df = d->df_raw * d->factor;

    if (grad){
        double dS_Av   = creal(clog(x_Av+eps))+1.0;
        double dS_ANa  = creal(clog(x_ANa+eps))+1.0;
        double dS_M1Mg = creal(clog(x_M1Mg+eps))+1.0;
        double dS_M1Al = creal(clog(x_M1Al+eps))+1.0;

        static const double dAv_dp[3]   = {1,1,0};
        static const double dANa_dp[3]  = {0,0,1};
        static const double dM1Mg_dp[3] = {1,0,0.5};
        static const double dM1Al_dp[3] = {0,1,0.5};

        double dG_tr   = W[0]*p_tsch + W[1]*p_parg;
        double dG_tsch = W[0]*p_tr   + W[2]*p_parg;
        double dG_parg = W[1]*p_tr   + W[2]*p_tsch;
        double dG[3] = { dG_tr, dG_tsch, dG_parg };

        for (int k = 0; k < 3; k++){
            double dS =
                1.0*( dS_Av*dAv_dp[k] + dS_ANa*dANa_dp[k] ) +
                2.0*( dS_M1Mg*dM1Mg_dp[k] + dS_M1Al*dM1Al_dp[k] );

            grad[k] = (gb[k] + R*T*dS + dG[k]/1000.0)*d->factor
                    - (d->df_raw*d->factor*(d->ape[k]/d->sum_apep));
        }
    }

    return d->df;
}

/* Feldspar: (EXT,MARGULES), plain ideal ternary entropy (no site given) +
   3-pair subregular + ternary Margules, Fuhrman & Lindsley (1988).
   Endmember order: ab, kfs, an. Al-Si order-disorder (D1/D2) dropped -
   fixed-G0 endmembers, per the decided simplification. */
double obj_br_fsp(unsigned n, const double *x, double *grad, void *SS_ref_db){
    (void) n;
    SS_ref *d = (SS_ref *) SS_ref_db;

    double R = d->R;
    double T = d->T;
    double *p  = d->p;
    double *gb = d->gb_lvl;
    double *W  = d->W;
    double eps = 1e-10;

    for (int i = 0; i < 3; i++){ p[i] = x[i]; }
    double px = p[0]; /* ab */
    double py = p[1]; /* kfs */
    double pz = p[2]; /* an */

    d->sf[0] = px; d->sf[1] = py; d->sf[2] = pz;

    double S_terms = ( px*creal(clog(px+eps)) + py*creal(clog(py+eps)) + pz*creal(clog(pz+eps)) );

    /* W[0..1]=Ab-Kfs (W_ab2kfs,W_abkfs2), W[2..3]=Ab-An, W[4..5]=An-Kfs, W[6]=ternary */
    double G_xs = W[0]*px*px*py + W[1]*px*py*py
                + W[2]*px*px*pz + W[3]*px*pz*pz
                + W[4]*py*pz*pz + W[5]*py*py*pz
                + W[6]*px*py*pz;

    d->sum_apep = 0.0;
    for (int i = 0; i < 3; i++){ d->sum_apep += d->ape[i]*p[i]; }
    d->factor = d->fbc/d->sum_apep;

    d->df_raw = gb[0]*p[0] + gb[1]*p[1] + gb[2]*p[2] + R*T*S_terms + G_xs/1000.0;
    d->df = d->df_raw * d->factor;

    if (grad){
        double dSx = creal(clog(px+eps))+1.0;
        double dSy = creal(clog(py+eps))+1.0;
        double dSz = creal(clog(pz+eps))+1.0;

        double dGx = 2.0*W[0]*px*py + W[1]*py*py + 2.0*W[2]*px*pz + W[3]*pz*pz + W[6]*py*pz;
        double dGy = W[0]*px*px + 2.0*W[1]*px*py + W[4]*pz*pz + 2.0*W[5]*py*pz + W[6]*px*pz;
        double dGz = W[2]*px*px + 2.0*W[3]*px*pz + 2.0*W[4]*py*pz + W[5]*py*py + W[6]*px*py;

        grad[0] = (gb[0] + R*T*dSx + dGx/1000.0)*d->factor - (d->df_raw*d->factor*(d->ape[0]/d->sum_apep));
        grad[1] = (gb[1] + R*T*dSy + dGy/1000.0)*d->factor - (d->df_raw*d->factor*(d->ape[1]/d->sum_apep));
        grad[2] = (gb[2] + R*T*dSz + dGz/1000.0)*d->factor - (d->df_raw*d->factor*(d->ape[2]/d->sum_apep));
    }

    return d->df;
}

void BR_SS_objective_init_function(obj_type *SS_objective, global_variable gv){
    for (int iss = 0; iss < gv.len_ss; iss++){
        if (strcmp( gv.SS_list[iss], "ctd") == 0 ){
            SS_objective[iss] = obj_br_ctd;
        }
        else if (strcmp( gv.SS_list[iss], "car") == 0 ){
            SS_objective[iss] = obj_br_car;
        }
        else if (strcmp( gv.SS_list[iss], "chl") == 0 ){
            SS_objective[iss] = obj_br_chl;
        }
        else if (strcmp( gv.SS_list[iss], "mica") == 0 ){
            SS_objective[iss] = obj_br_mica;
        }
        else if (strcmp( gv.SS_list[iss], "talc") == 0 ){
            SS_objective[iss] = obj_br_talc;
        }
        else if (strcmp( gv.SS_list[iss], "ilm") == 0 ){
            SS_objective[iss] = obj_br_ilm;
        }
        else if (strcmp( gv.SS_list[iss], "bt") == 0 ){
            SS_objective[iss] = obj_br_bt;
        }
        else if (strcmp( gv.SS_list[iss], "ol") == 0 ){
            SS_objective[iss] = obj_br_ol;
        }
        else if (strcmp( gv.SS_list[iss], "ep") == 0 ){
            SS_objective[iss] = obj_br_ep;
        }
        else if (strcmp( gv.SS_list[iss], "opx") == 0 ){
            SS_objective[iss] = obj_br_opx;
        }
        else if (strcmp( gv.SS_list[iss], "amph") == 0 ){
            SS_objective[iss] = obj_br_amph;
        }
        else if (strcmp( gv.SS_list[iss], "spl") == 0 ){
            SS_objective[iss] = obj_br_spl;
        }
        else if (strcmp( gv.SS_list[iss], "stau") == 0 ){
            SS_objective[iss] = obj_br_stau;
        }
        else if (strcmp( gv.SS_list[iss], "crd") == 0 ){
            SS_objective[iss] = obj_br_crd;
        }
        else if (strcmp( gv.SS_list[iss], "grt") == 0 ){
            SS_objective[iss] = obj_br_grt;
        }
        else if (strcmp( gv.SS_list[iss], "omph") == 0 ){
            SS_objective[iss] = obj_br_omph;
        }
        else if (strcmp( gv.SS_list[iss], "amphx") == 0 ){
            SS_objective[iss] = obj_br_amphx;
        }
        else if (strcmp( gv.SS_list[iss], "fsp") == 0 ){
            SS_objective[iss] = obj_br_fsp;
        }
    }
}

void BR_PC_init(PC_type *PC_read, global_variable gv){
    for (int iss = 0; iss < gv.len_ss; iss++){
        if (strcmp( gv.SS_list[iss], "ctd") == 0 ){
            PC_read[iss] = obj_br_ctd;
        }
        else if (strcmp( gv.SS_list[iss], "car") == 0 ){
            PC_read[iss] = obj_br_car;
        }
        else if (strcmp( gv.SS_list[iss], "chl") == 0 ){
            PC_read[iss] = obj_br_chl;
        }
        else if (strcmp( gv.SS_list[iss], "mica") == 0 ){
            PC_read[iss] = obj_br_mica;
        }
        else if (strcmp( gv.SS_list[iss], "talc") == 0 ){
            PC_read[iss] = obj_br_talc;
        }
        else if (strcmp( gv.SS_list[iss], "ilm") == 0 ){
            PC_read[iss] = obj_br_ilm;
        }
        else if (strcmp( gv.SS_list[iss], "bt") == 0 ){
            PC_read[iss] = obj_br_bt;
        }
        else if (strcmp( gv.SS_list[iss], "ol") == 0 ){
            PC_read[iss] = obj_br_ol;
        }
        else if (strcmp( gv.SS_list[iss], "ep") == 0 ){
            PC_read[iss] = obj_br_ep;
        }
        else if (strcmp( gv.SS_list[iss], "opx") == 0 ){
            PC_read[iss] = obj_br_opx;
        }
        else if (strcmp( gv.SS_list[iss], "amph") == 0 ){
            PC_read[iss] = obj_br_amph;
        }
        else if (strcmp( gv.SS_list[iss], "spl") == 0 ){
            PC_read[iss] = obj_br_spl;
        }
        else if (strcmp( gv.SS_list[iss], "stau") == 0 ){
            PC_read[iss] = obj_br_stau;
        }
        else if (strcmp( gv.SS_list[iss], "crd") == 0 ){
            PC_read[iss] = obj_br_crd;
        }
        else if (strcmp( gv.SS_list[iss], "grt") == 0 ){
            PC_read[iss] = obj_br_grt;
        }
        else if (strcmp( gv.SS_list[iss], "omph") == 0 ){
            PC_read[iss] = obj_br_omph;
        }
        else if (strcmp( gv.SS_list[iss], "amphx") == 0 ){
            PC_read[iss] = obj_br_amphx;
        }
        else if (strcmp( gv.SS_list[iss], "fsp") == 0 ){
            PC_read[iss] = obj_br_fsp;
        }
    }
}

void p2x_br_generic(void *SS_ref_db, double eps){
    (void) eps;
    SS_ref *d = (SS_ref *) SS_ref_db;
    for (int i = 0; i < d->n_xeos; i++){
        d->iguess[i] = d->p[i];
        if (d->iguess[i] < d->bounds[i][0]){ d->iguess[i] = d->bounds[i][0]; }
        if (d->iguess[i] > d->bounds[i][1]){ d->iguess[i] = d->bounds[i][1]; }
    }
}

void BR_P2X_init(P2X_type *P2X_read, global_variable gv){
    for (int iss = 0; iss < gv.len_ss; iss++){
        if (strcmp( gv.SS_list[iss], "ctd") == 0 ||
            strcmp( gv.SS_list[iss], "car") == 0 ||
            strcmp( gv.SS_list[iss], "chl") == 0 ||
            strcmp( gv.SS_list[iss], "mica") == 0 ||
            strcmp( gv.SS_list[iss], "talc") == 0 ||
            strcmp( gv.SS_list[iss], "ilm") == 0 ||
            strcmp( gv.SS_list[iss], "bt") == 0 ||
            strcmp( gv.SS_list[iss], "ol") == 0 ||
            strcmp( gv.SS_list[iss], "ep") == 0 ||
            strcmp( gv.SS_list[iss], "opx") == 0 ||
            strcmp( gv.SS_list[iss], "amph") == 0 ||
            strcmp( gv.SS_list[iss], "spl") == 0 ||
            strcmp( gv.SS_list[iss], "stau") == 0 ||
            strcmp( gv.SS_list[iss], "crd") == 0 ||
            strcmp( gv.SS_list[iss], "grt") == 0 ||
            strcmp( gv.SS_list[iss], "omph") == 0 ||
            strcmp( gv.SS_list[iss], "amphx") == 0 ||
            strcmp( gv.SS_list[iss], "fsp") == 0 ){
            P2X_read[iss] = p2x_br_generic;
        }
    }
}

SS_ref BR_PC_function(global_variable gv, PC_type *PC_read, SS_ref SS_ref_db, bulk_info z_b, int ph_id){
    (void) z_b;
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
