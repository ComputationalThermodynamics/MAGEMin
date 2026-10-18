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

#include <math.h>
#include <float.h>
#include <complex.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#include "../MAGEMin.h"
#include "gh_objective_functions.h"

#define GH_R_J  8.3143
#define GH_R_KJ 0.0083143
#include "GH_fluid_eos.h"
#include "GH_gem_function.h"

static inline void gex_sym_n2(SS_ref *d, const double *p, int jmax, double *Gex)
{
    int n = d->n_em, it = 0;
    double Q = 0.0;
    for (int i = 0; i < n; i++){ Gex[i] = 0.0; }
    for (int j = 0; j < jmax; j++){
        for (int k = j+1; k < n; k++){
            double w = d->W[it++];
            Gex[j] += w*p[k];
            Gex[k] += w*p[j];
            Q      += w*p[j]*p[k];
        }
    }
    for (int i = 0; i < n; i++){ Gex[i] -= Q; }
}

static double GH_df_assemble(SS_ref *d, const double *p, const double *gb, const double *mu_Gex, int n, double *grad){
    d->sum_apep = 0.0;
    for (int i = 0; i < n; i++){ d->sum_apep += d->ape[i]*p[i]; }
    d->factor = d->fbc/d->sum_apep;

    d->df_raw = 0.0;
    for (int i = 0; i < n; i++){ d->df_raw += (mu_Gex[i] + gb[i])*p[i]; }
    d->df = d->df_raw * d->factor;

    if (grad){
        for (int i = 0; i < n; i++){
            grad[i] = (mu_Gex[i] + gb[i])*d->factor - (d->df_raw*d->factor*(d->ape[i]/d->sum_apep));
        }
    }
    return d->df;
}

static int GH_nograd(const double *grad, const double *p, int n, double norm){
    if (grad){ return 0; }
    double S = 0.0;
    for (int i = 0; i < n; i++){ S += p[i]; }
    return fabs(S - norm) <= 1.0e-12*fabs(norm);
}

static double GH_df_fast(SS_ref *d, const double *p, const double *gb, int n, double Gex, double Sconfig){
    d->sum_apep = 0.0;
    double S = 0.0;
    d->df_raw = 0.0;
    for (int i = 0; i < n; i++){
        d->sum_apep += d->ape[i]*p[i];
        d->df_raw   += gb[i]*p[i];
        S           += p[i];
    }
    d->factor  = d->fbc/d->sum_apep;
    d->df_raw += S*Gex/1000.0 + Sconfig;
    d->df      = d->df_raw * d->factor;
    return d->df;
}

static void GH_mu_from_dGex(SS_ref *d, const double *p, double norm, int n, double Gex, const double *dGex){
    for (int i = 0; i < n; i++){
        double mu = Gex;
        for (int k = 0; k < n; k++){
            double delta_ik = (i==k) ? 1.0 : 0.0;
            mu += (delta_ik - p[k]/norm)*dGex[k];
        }
        d->mu_Gex[i] = mu/1000.0;
        d->sf[i]     = p[i];
    }
}

typedef struct { double h, s, v; } GH_hsv;

static inline double GH_hsv_G(GH_hsv c, double T, double dP){
    return c.h - T*c.s + dP*c.v;
}

static void GH_ternary_cubic(const GH_hsv W[6], GH_hsv W123, double T, double dP, const double *p, double *G, double *dG){
    GH_hsv K = { W123.h, W123.s, W123.v };
    double w[6];
    for (int i = 0; i < 6; i++){
        K.h += 0.5*W[i].h; K.s += 0.5*W[i].s; K.v += 0.5*W[i].v;
        w[i] = GH_hsv_G(W[i], T, dP);
    }
    double k = GH_hsv_G(K, T, dP);
    double a = p[0], b = p[1], c = p[2];
    *G = w[0]*a*a*b + w[1]*a*b*b + w[2]*a*a*c + w[3]*a*c*c + w[4]*b*b*c + w[5]*b*c*c + a*b*c*k;
    dG[0] = 2.0*w[0]*a*b + w[1]*b*b + 2.0*w[2]*a*c + w[3]*c*c + b*c*k;
    dG[1] = w[0]*a*a + 2.0*w[1]*a*b + 2.0*w[4]*b*c + w[5]*c*c + a*c*k;
    dG[2] = w[2]*a*a + 2.0*w[3]*a*c + w[4]*b*b + 2.0*w[5]*b*c + a*b*k;
}

static void GH_invert3(const double H[3][3], double inv[3][3]){
    double c00 = H[1][1]*H[2][2]-H[1][2]*H[2][1], c01 = H[1][2]*H[2][0]-H[1][0]*H[2][2], c02 = H[1][0]*H[2][1]-H[1][1]*H[2][0];
    double det = H[0][0]*c00 + H[0][1]*c01 + H[0][2]*c02;
    inv[0][0] = c00/det; inv[0][1] = (H[0][2]*H[2][1]-H[0][1]*H[2][2])/det; inv[0][2] = (H[0][1]*H[1][2]-H[0][2]*H[1][1])/det;
    inv[1][0] = c01/det; inv[1][1] = (H[0][0]*H[2][2]-H[0][2]*H[2][0])/det; inv[1][2] = (H[0][2]*H[1][0]-H[0][0]*H[1][2])/det;
    inv[2][0] = c02/det; inv[2][1] = (H[0][1]*H[2][0]-H[0][0]*H[2][1])/det; inv[2][2] = (H[0][0]*H[1][1]-H[0][1]*H[1][0])/det;
}

static __attribute__((noinline)) double rlog(double complex z)
{
#if defined(__APPLE__)
    return (cimag(z) == 0.0) ? log(fabs(creal(z))) : creal(clog(z));
#else
    return creal(clog(z));
#endif
}

static const double GH_rmelts_Wcc[13] = {
     63281.2,
    -79202.7,
     46716.1,
     65508.7,
         0.0,
    -72996.8,
         0.0,
    -24872.6,
     37534.1,
   -311011.3,
    -27865.0,
         0.0,
      7873.0
};

static const GH_hsv GH_RMELTS_CACO3_C = { -17574.497522747, 0.0, -1.9034060173857 };

static double GH_config_entropy(const SS_ref *d, const double *p, int n, double mult, double R, double T, double *dSi){
    double Sconfig = 0.0;
    for (int i = 0; i < n; i++){
        double lg = rlog((p[i] + d->d_em[i]));
        Sconfig += p[i]*lg;
        dSi[i]   = mult*R*T*(lg + ((d->d_em[i] == 0.0) ? 1.0 : p[i]/(p[i] + d->d_em[i])));
    }
    Sconfig *= mult*R*T;
    return Sconfig;
}

static double GH_df_assemble_cfg(SS_ref *d, const double *p, const double *gb, const double *mu_Gex, int n, double *grad,
                                 double Sconfig, const double *dSi){
    d->sum_apep = 0.0;
    for (int i = 0; i < n; i++){ d->sum_apep += d->ape[i]*p[i]; }
    d->factor = d->fbc/d->sum_apep;

    d->df_raw = 0.0;
    for (int i = 0; i < n; i++){ d->df_raw += (mu_Gex[i] + gb[i])*p[i]; }
    d->df_raw += Sconfig;
    d->df = d->df_raw * d->factor;

    if (grad){
        for (int i = 0; i < n; i++){
            grad[i] = (mu_Gex[i] + gb[i] + dSi[i])*d->factor - (d->df_raw*d->factor*(d->ape[i]/d->sum_apep));
        }
    }
    return d->df;
}

static inline int GH_rmelts_widx(int i, int k){
    return i*13 - i*(i+1)/2 + (k-i-1);
}

static inline double GH_rmelts_Wsym(const double *W, int j, int k){
    return (j < k) ? W[GH_rmelts_widx(j,k)] : W[GH_rmelts_widx(k,j)];
}

static void GH_rmelts_Wv(const double *W, double Wv[13]){
    for (int j = 0; j < 13; j++){
        Wv[j] = ((j != 0) ? GH_rmelts_Wsym(W, j, 0) : 0.0)
              - ((j != 8) ? GH_rmelts_Wsym(W, j, 8) : 0.0)
              - ((j != 11) ? GH_rmelts_Wsym(W, j, 11) : 0.0);
    }
}

static inline double GH_xlogx(double x){
    return (x > 0.0) ? x*rlog((x)) : 0.0;
}

int GH_sf_constrained(const char *name){
    return strcmp(name, "cpx") == 0 || strcmp(name, "opx") == 0 || strcmp(name, "spn") == 0 || strcmp(name, "ol") == 0
        || strcmp(name, "hb") == 0 || strcmp(name, "nph") == 0 || strcmp(name, "kls") == 0 || strcmp(name, "amp") == 0
        || strcmp(name, "lc") == 0;
}

int GH_n_sf(const char *name, int n_em){
    if (strcmp(name, "cpx") == 0) return 22;
    if (strcmp(name, "opx") == 0) return 25;
    if (strcmp(name, "spn") == 0) return 10;
    if (strcmp(name, "ol")  == 0) return 8;
    if (strcmp(name, "hb")  == 0) return 4;
    if (strcmp(name, "nph") == 0 || strcmp(name, "kls") == 0) return 7;
    if (strcmp(name, "amp") == 0) return 6;
    if (strcmp(name, "lc")  == 0) return 4;
    return n_em;
}

static void GH_sf_pyx(const double *p, double *sf, int opx){
    double S   = p[0]+p[1]+p[2]+p[3]+p[4]+p[5]+p[6];
    double Fe2 = p[2], Na = p[6], Al = p[3]+p[5]+p[6], Fe3 = p[4]+p[5], Ti = p[3]+p[4];
    double Ca  = p[0]+p[2]+p[3]+p[4]+p[5], Mg = p[0]+2.0*p[1]+0.5*(p[3]+p[4]);
    sf[0]  = 0.5*Fe2;        sf[1]  = 0.5*(2.0*S - Fe2);
    sf[2]  = Na;             sf[3]  = S - Na;
    sf[4]  = 0.5*Al;         sf[5]  = 0.5*(2.0*S - Al);
    sf[6]  = 0.5*Fe3;        sf[7]  = 0.5*(2.0*S - Fe3);
    sf[8]  = Ti;             sf[9]  = S - Ti;
    sf[10] = Ca;             sf[11] = S - Ca;
    sf[12] = 0.5*Mg;         sf[13] = 0.5*(2.0*S - Mg);
    sf[14] = 0.25*(p[6] + p[1]);
    sf[15] = 0.25*(S - p[6] - p[1]);
    if (opx){
        double Mgo = S + p[1] - p[2] - 0.5*(p[3] + p[4]) - p[5] - 1.5*p[6];
        sf[16] = 0.5*p[1];
        sf[17] = 0.25*(S - p[3] - p[4] - p[5] - p[6]);
        sf[18] = 0.25*(p[5] + 0.5*(p[3] + p[4]));
        sf[19] = 0.25*(p[5] + p[6] + 0.5*(p[3] + p[4]));
        sf[20] = 0.25*(S - p[0] - p[1] - p[2] - p[6]);
        sf[21] = 0.25*Mgo;
        sf[22] = 0.25*(2.0*S - Mgo);
        sf[23] = 0.25*(0.5*(p[3] + p[4] + p[5]) + 0.75*p[6]);
        sf[24] = 0.25*(p[5] + p[6]);
    }
    else{
        sf[16] = 0.25*p[1];
        sf[17] = 0.25*(p[5] + p[6]);
        sf[18] = 0.25*(0.5*(p[3] + p[4]) + p[5]);
        sf[19] = 0.25*(S - p[3] - p[4] - p[5] - p[6]);
        sf[20] = 0.25*(p[3] + p[4] + p[5]);
        sf[21] = 0.25*(S - p[3] - p[4] - p[5]);
    }
}

static void GH_sf_spn(const double *p, double *sf){
    double S = p[0]+p[1]+p[2]+p[3]+p[4];
    sf[0] = p[0];              sf[1] = S - p[0];
    sf[2] = p[2];              sf[3] = S - p[2];
    sf[4] = p[4];              sf[5] = S - p[4];
    sf[6] = 0.5*p[3];          sf[7] = 0.5*(S + p[4] - p[3]);
    sf[8] = p[1] + p[3];       sf[9] = 0.5*(S - p[1]);
}

static void GH_sf_ol(const double *p, double *sf){
    double S = p[0]+p[1]+p[2]+p[3];
    sf[0] = p[2];              sf[1] = S - p[2];
    sf[2] = p[1];              sf[3] = S - p[1];
    sf[4] = p[3];              sf[5] = S - p[3];
    sf[6] = p[0] + 0.5*p[3];   sf[7] = 0.5*(S - p[0]);
}

static void GH_sf_hb(const double *p, double *sf){
    sf[0] = p[1];
    sf[1] = p[2];
    sf[2] = p[0] + p[1];
    sf[3] = p[0] + p[2];
}

static void GH_sf_nph(const double *p, double *sf){
    double S = p[0]+p[1]+p[2]+p[3];
    sf[0] = 0.5*(p[0] + 0.75*p[2] + 0.5*p[3]);
    sf[1] = 0.25*(S - p[0]);
    sf[2] = p[1];
    sf[3] = S - p[1] - 0.25*p[2] - 0.5*p[3];
    sf[4] = p[2];
    sf[5] = p[3];
    sf[6] = S - p[2] - p[3];
}

static void GH_sf_amp(const double *p, double *sf){
    double S = p[0]+p[1]+p[2];
    sf[0] = p[1];
    sf[1] = S - 2.0*p[2]/7.0 - p[1];
    sf[2] = p[2];
    sf[3] = S - p[2];
    sf[4] = p[0] + 5.0*p[2]/7.0;
    sf[5] = S - p[0] - p[2];
}

void GH_sf_of_p(const char *name, const double *p, double *sf){
    if      (strcmp(name, "cpx") == 0){ GH_sf_pyx(p, sf, 0); }
    else if (strcmp(name, "opx") == 0){ GH_sf_pyx(p, sf, 1); }
    else if (strcmp(name, "spn") == 0){ GH_sf_spn(p, sf); }
    else if (strcmp(name, "ol")  == 0){ GH_sf_ol(p, sf); }
    else if (strcmp(name, "hb")  == 0){ GH_sf_hb(p, sf); }
    else if (strcmp(name, "nph") == 0 || strcmp(name, "kls") == 0){ GH_sf_nph(p, sf); }
    else if (strcmp(name, "amp") == 0){ GH_sf_amp(p, sf); }
    else if (strcmp(name, "lc")  == 0){ sf[0] = p[0]; sf[1] = p[1]; sf[2] = p[2] + p[0]; sf[3] = p[2] + p[1]; }
}

double obj_gh_liq(unsigned n, const double *x, double *grad, void *SS_ref_db){
    SS_ref *d = (SS_ref *) SS_ref_db;

    int n_em    = d->n_em;
    double R    = GH_R_KJ;
    double T    = d->T;
    double *p   = d->p;
    double *gb  = d->gb_lvl;
    double *mu_Gex = d->mu_Gex;

    for (int i = 0; i < n_em; i++){
        p[i] = x[i];
    }

    double Gex_n2[n_em];
    gex_sym_n2(d, p, n_em, Gex_n2);
    for (int i = 0; i < n_em; i++){
        double Gex = Gex_n2[i];
        mu_Gex[i] = Gex/1000.0;
        d->sf[i]  = p[i];
    }

    double dSi[13];
    double Sconfig = GH_config_entropy(d, p, n_em, 1.0, R, T, dSi);

    double x_h2o = p[n_em-1] + d->d_em[n_em-1];

    double log_x_h2o   = (x_h2o > 0.0) ? rlog((x_h2o))       : 0.0;
    double log_1mx_h2o = (x_h2o < 1.0) ? rlog((1.0-x_h2o))   : 0.0;
    double Sconfig_h2o = 0.0;
    if (d->d_em[n_em-1] == 0.0){
        Sconfig_h2o = R*T*(x_h2o*log_x_h2o + (1.0-x_h2o)*log_1mx_h2o);
        dSi[n_em-1] += R*T*(log_x_h2o - log_1mx_h2o);
    }

    double dGex_rmelts = 0.0;
    double dSi_rmelts[13] = {0.0};
    if (d->EM_database == 1 && n_em == 13){

        const int iSi = 0, iCa = 8, iCO2 = 11;
        double Rgas = R*1000.0;
        double Pbar = d->P*1000.0;
        double C = GH_hsv_G(GH_RMELTS_CACO3_C, 0.0, Pbar-1.0);

        double r9 = p[iCa], r13 = p[iCO2], x0 = p[iSi];
        double a0 = x0 + d->d_em[iSi], a9 = r9 + d->d_em[iCa], a13 = r13 + d->d_em[iCO2];
        double m = fmin(r9, r13);
        double s = 0.0;

        if (m > 0.0 && d->d_em[iCa] == 0.0 && d->d_em[iCO2] == 0.0){
            double Wv[13];
            GH_rmelts_Wv(d->W, Wv);
            double Bcoef = 0.0;
            for (int k = 0; k < 13; k++){ Bcoef += (Wv[k] + GH_rmelts_Wcc[k])*p[k]; }
            double Dcoef = GH_rmelts_Wsym(d->W, iCa, iCO2) - GH_rmelts_Wsym(d->W, iSi, iCa) - GH_rmelts_Wsym(d->W, iSi, iCO2)
                         + GH_rmelts_Wcc[iSi] - GH_rmelts_Wcc[iCa] - GH_rmelts_Wcc[iCO2];

            double sigma = 0.5;
            double lo = 1.0e-14, hi = 1.0 - 1.0e-14;
            if (2.0*Dcoef + Rgas*T*(1.0/m + 1.0/(a0+m) + 1.0/a9 + 1.0/a13) > 0.0){
                double z = -(Bcoef + C)/(Rgas*T);
                z = (z > 600.0) ? 600.0 : ((z < -600.0) ? -600.0 : z);
                double K  = exp(z);
                double qa = 1.0 - K, qb = a0 + K*(a9 + a13), qc = K*a9*a13;
                double disc = qb*qb + 4.0*qa*qc;
                if (disc >= 0.0){
                    double s0 = 2.0*qc/(qb + sqrt(disc));
                    double sg = s0/m;
                    if (sg > lo && sg < hi){ sigma = sg; }
                }
            }
            for (int iter = 0; iter < 100; iter++){
                double sv = m*sigma;
                double xSi = a0+sv, xCa = a9-sv, xCO2v = a13-sv;
                double f   = Bcoef + 2.0*Dcoef*sv + Rgas*T*rlog(((xSi*sv)/(xCa*xCO2v))) + C;
                double dfds = 2.0*Dcoef + Rgas*T*(1.0/xSi + 1.0/xCa + 1.0/xCO2v + 1.0/sv);
                double fp  = dfds*m;
                double sigma_new = sigma - f/fp;
                if (!(sigma_new > lo && sigma_new < hi)) sigma_new = 0.5*(sigma + (f > 0.0 ? lo : hi));
                if (fabs(sigma_new - sigma) < 1.0e-15){ sigma = sigma_new; break; }
                sigma = sigma_new;
            }
            s = m*sigma;

            double dent = GH_xlogx(a0+s) - GH_xlogx(a0) + GH_xlogx(a9-s) - GH_xlogx(a9)
                        + GH_xlogx(a13-s) - GH_xlogx(a13) + GH_xlogx(s);
            dGex_rmelts = (Bcoef*s + Dcoef*s*s + Rgas*T*dent + s*C) / 1000.0;

            if (grad){
                double xSi_s = a0+s, xCa_s = a9-s, xCO2_s = a13-s;
                for (int j = 0; j < 13; j++){
                    dSi_rmelts[j] = s*(Wv[j] + GH_rmelts_Wcc[j]) / 1000.0;
                }
                dSi_rmelts[iSi]  += R*T*(rlog((xSi_s))  - rlog((a0)));
                dSi_rmelts[iCa]  += R*T*(rlog((xCa_s))  - rlog((a9)));
                dSi_rmelts[iCO2] += R*T*(rlog((xCO2_s)) - rlog((a13)));
            }
        }
    }

    double dSi_tot[13];
    for (int i = 0; i < n_em; i++){ dSi_tot[i] = dSi[i] + dSi_rmelts[i]; }
    return GH_df_assemble_cfg(d, p, gb, mu_Gex, n_em, grad, Sconfig + Sconfig_h2o + dGex_rmelts, dSi_tot);
}

typedef struct { GH_hsv ex, x, w1, w2; } GH_ol_pair;
typedef struct { GH_hsv w2ca; double f; } GH_ol_cat;
enum { GH_OL_MGMN, GH_OL_MGFE, GH_OL_MGCO, GH_OL_MGNI, GH_OL_MNFE, GH_OL_MNCO, GH_OL_MNNI, GH_OL_FECO, GH_OL_FENI, GH_OL_CONI, GH_OL_NPAIR };
enum { GH_OL_MG, GH_OL_MN, GH_OL_FE, GH_OL_CO, GH_OL_NI, GH_OL_NCAT };
static const GH_ol_pair GH_OL_PAIR[GH_OL_NPAIR] = {
    [GH_OL_MGMN] = { .ex = {15800.0, 0.0, 0.0}, .x = {8750.0, 0.0, 0.0}, .w1 = {6625.0, 0.0, 0.0}, .w2 = {6625.0, 0.0, 0.0} },
    [GH_OL_MGFE] = { .ex = {0.0, 0.0, 0.0}, .x = {10150.0, 0.0, 0.015}, .w1 = {5075.0, 0.0, 0.0}, .w2 = {5075.0, 0.0, 0.0} },
    [GH_OL_MGCO] = { .ex = {-15000.0, 0.0, 0.0}, .x = {3000.0, 0.0, 0.0}, .w1 = {1500.0, 0.0, 0.0}, .w2 = {1500.0, 0.0, 0.0} },
    [GH_OL_MGNI] = { .ex = {-19750.0, 0.0, 0.0}, .x = {2200.0, 0.0, 0.0}, .w1 = {-600.0, 0.0, 0.0}, .w2 = {2800.0, 0.0, 0.0} },
    [GH_OL_MNFE] = { .ex = {-11800.0, 0.0, 0.0}, .x = {500.0, 0.0, 0.0}, .w1 = {1750.0, 0.0, 0.0}, .w2 = {1750.0, 0.0, 0.0} },
    [GH_OL_MNCO] = { .ex = {0.0, 0.0, 0.0}, .x = {0.0, 0.0, 0.0}, .w1 = {0.0, 0.0, 0.0}, .w2 = {0.0, 0.0, 0.0} },
    [GH_OL_MNNI] = { .ex = {0.0, 0.0, 0.0}, .x = {0.0, 0.0, 0.0}, .w1 = {0.0, 0.0, 0.0}, .w2 = {0.0, 0.0, 0.0} },
    [GH_OL_FECO] = { .ex = {0.0, 0.0, 0.0}, .x = {3000.0, 0.0, 0.0}, .w1 = {1500.0, 0.0, 0.0}, .w2 = {1500.0, 0.0, 0.0} },
    [GH_OL_FENI] = { .ex = {-20000.0, 0.0, 0.0}, .x = {10000.0, 0.0, 0.0}, .w1 = {5000.0, 0.0, 0.0}, .w2 = {5000.0, 0.0, 0.0} },
    [GH_OL_CONI] = { .ex = {0.0, 0.0, 0.0}, .x = {0.0, 0.0, 0.0}, .w1 = {0.0, 0.0, 0.0}, .w2 = {0.0, 0.0, 0.0} },
};
static const GH_ol_cat GH_OL_CAT[GH_OL_NCAT] = {
    [GH_OL_MG] = { .w2ca = {34500.0, 0.0, 0.35}, .f = 0.0 },
    [GH_OL_MN] = { .w2ca = {16000.0, 0.0, 0.0}, .f = 9500.0 },
    [GH_OL_FE] = { .w2ca = {21900.0, 0.0, 0.0}, .f = 9500.0 },
    [GH_OL_CO] = { .w2ca = {30000.0, 0.0, 0.0}, .f = 0.0 },
    [GH_OL_NI] = { .w2ca = {40000.0, 0.0, 0.0}, .f = 0.0 },
};

enum { GH_OL_R1, GH_OL_R2, GH_OL_R3, GH_OL_R4, GH_OL_R5, GH_OL_S1, GH_OL_S2, GH_OL_S3, GH_OL_S4, GH_OL_NZ };

typedef struct {
    double g0;
    double lin[GH_OL_NZ];
    double quad[GH_OL_NZ][GH_OL_NZ];
} GH_ol_coef;

typedef struct {
    double m1[4], m2[4];
    double ca2, mg1, mg2;
} GH_ol_sf;

static void GH_ol_coef_at(double p, GH_ol_coef *c){
    enum { MGMN = GH_OL_MGMN, MGFE = GH_OL_MGFE, MGCO = GH_OL_MGCO, MGNI = GH_OL_MGNI, MNFE = GH_OL_MNFE,
           MNCO = GH_OL_MNCO, MNNI = GH_OL_MNNI, FECO = GH_OL_FECO, FENI = GH_OL_FENI, CONI = GH_OL_CONI };
    enum { MG = GH_OL_MG, MN = GH_OL_MN, FE = GH_OL_FE, CO = GH_OL_CO, NI = GH_OL_NI };
    enum { R1 = GH_OL_R1, R2 = GH_OL_R2, R3 = GH_OL_R3, R4 = GH_OL_R4, R5 = GH_OL_R5,
           S1 = GH_OL_S1, S2 = GH_OL_S2, S3 = GH_OL_S3, S4 = GH_OL_S4 };

    double ex[GH_OL_NPAIR], x[GH_OL_NPAIR], w1[GH_OL_NPAIR], w2[GH_OL_NPAIR];
    double A[GH_OL_NPAIR], B[GH_OL_NPAIR], C[GH_OL_NPAIR], D[GH_OL_NPAIR];
    double w2ca[GH_OL_NCAT], f[GH_OL_NCAT];

    for (int k = 0; k < GH_OL_NPAIR; k++){
        ex[k] = GH_hsv_G(GH_OL_PAIR[k].ex, 0.0, p - 1.0);
        x[k]  = GH_hsv_G(GH_OL_PAIR[k].x, 0.0, p - 1.0);
        w1[k] = GH_hsv_G(GH_OL_PAIR[k].w1, 0.0, p - 1.0);
        w2[k] = GH_hsv_G(GH_OL_PAIR[k].w2, 0.0, p - 1.0);
        A[k]  = x[k] + w1[k] + w2[k];
        B[k]  = x[k] - w1[k] - w2[k];
        C[k]  = ex[k] + w1[k] - w2[k];
        D[k]  = ex[k] - w1[k] + w2[k];
    }
    for (int k = 0; k < GH_OL_NCAT; k++){
        w2ca[k] = GH_hsv_G(GH_OL_CAT[k].w2ca, 0.0, p - 1.0);
        f[k]    = GH_OL_CAT[k].f;
    }

    c->g0 = 0.25*(A[MNFE] + A[MNCO] + A[MNNI] + A[FECO] + A[FENI] + A[CONI] - 2.0*A[MGMN] - 2.0*A[MGFE]
            - 2.0*A[MGCO] - 2.0*A[MGNI]);
    c->lin[R1] = 0.25*(-3.0*A[MGMN] - A[MGFE] - A[MGCO] - A[MGNI] + A[MNFE] + A[MNCO] + A[MNNI]);
    c->lin[R2] = 0.25*(-A[MGMN] - 3.0*A[MGFE] - A[MGCO] - A[MGNI] + A[MNFE] + A[FECO] + A[FENI]);
    c->lin[R3] = 0.25*(-A[MGMN] - A[MGFE] - 3.0*A[MGCO] - A[MGNI] + A[MNCO] + A[FECO] + A[CONI]);
    c->lin[R4] = 0.25*(-A[MGMN] - A[MGFE] - A[MGCO] - 3.0*A[MGNI] + A[MNNI] + A[FENI] + A[CONI]);
    c->lin[R5] = -w2ca[MG] - 0.25*((f[MN] + ex[MGMN] + x[MGMN] - 2.0*w2ca[MN] + 2.0*w2[MGMN]) + (f[FE]
            + ex[MGFE] + x[MGFE] - 2.0*w2ca[FE] + 2.0*w2[MGFE]) + (f[CO] + ex[MGCO] + x[MGCO]
            - 2.0*w2ca[CO] + 2.0*w2[MGCO]) + (f[NI] + ex[MGNI] + x[MGNI] - 2.0*w2ca[NI] + 2.0*w2[MGNI]));
    c->lin[S1] = 0.25*((ex[MGMN] + 3.0*w1[MGMN] - 3.0*w2[MGMN]) - D[MGFE] - D[MGCO] - D[MGNI] + D[MNFE]
            + D[MNCO] + D[MNNI]);
    c->lin[S2] = 0.25*(-D[MGMN] + (ex[MGFE] + 3.0*w1[MGFE] - 3.0*w2[MGFE]) - D[MGCO] - D[MGNI] - C[MNFE]
            + D[FECO] + D[FENI]);
    c->lin[S3] = 0.25*(D[MGMN] + D[MGFE] - (ex[MGCO] + 3.0*w1[MGCO] - 3.0*w2[MGCO]) + D[MGNI] + C[MNCO]
            + C[FECO] - D[CONI]);
    c->lin[S4] = 0.25*(D[MGMN] + D[MGFE] + D[MGCO] - (ex[MGNI] + 3.0*w1[MGNI] - 3.0*w2[MGNI]) + C[MNNI]
            + C[FENI] + C[CONI]);
    c->quad[R1][R1] = -0.25*A[MGMN];
    c->quad[R1][R2] = 0.25*(A[MNFE] - A[MGMN] - A[MGFE]);
    c->quad[R1][R3] = 0.25*(A[MNCO] - A[MGMN] - A[MGCO]);
    c->quad[R1][R4] = 0.25*(A[MNNI] - A[MGMN] - A[MGNI]);
    c->quad[R1][R5] = -0.25*(f[MN] + ex[MGMN] + x[MGMN] + 2.0*w2ca[MG] - 2.0*w2ca[MN] + 2.0*w2[MGMN]);
    c->quad[R1][S1] = 0.5*(w1[MGMN] - w2[MGMN]);
    c->quad[R1][S2] = 0.25*(-C[MGMN] + C[MGFE] - C[MNFE]);
    c->quad[R1][S3] = 0.25*(D[MGMN] - C[MGCO] + C[MNCO]);
    c->quad[R1][S4] = 0.25*(D[MGMN] - C[MGNI] + C[MNNI]);
    c->quad[R2][R2] = -0.25*A[MGFE];
    c->quad[R2][R3] = 0.25*(A[FECO] - A[MGFE] - A[MGCO]);
    c->quad[R2][R4] = 0.25*(A[FENI] - A[MGFE] - A[MGNI]);
    c->quad[R2][R5] = -0.25*(f[FE] + ex[MGFE] + x[MGFE] + 2.0*w2ca[MG] - 2.0*w2ca[FE] + 2.0*w2[MGFE]);
    c->quad[R2][S1] = 0.25*(C[MGMN] - D[MGFE] + D[MNFE]);
    c->quad[R2][S2] = 0.5*(w1[MGFE] - w2[MGFE]);
    c->quad[R2][S3] = 0.25*(D[MGFE] - C[MGCO] + C[FECO]);
    c->quad[R2][S4] = 0.25*(D[MGFE] - C[MGNI] + C[FENI]);
    c->quad[R3][R3] = -0.25*A[MGCO];
    c->quad[R3][R4] = 0.25*(A[CONI] - A[MGCO] - A[MGNI]);
    c->quad[R3][R5] = -0.25*(f[CO] + ex[MGCO] + x[MGCO] + 2.0*w2ca[MG] - 2.0*w2ca[CO] + 2.0*w2[MGCO]);
    c->quad[R3][S1] = 0.25*(C[MGMN] - D[MGCO] + D[MNCO]);
    c->quad[R3][S2] = 0.25*(C[MGFE] - D[MGCO] + D[FECO]);
    c->quad[R3][S3] = 0.5*(-w1[MGCO] + w2[MGCO]);
    c->quad[R3][S4] = 0.25*(C[MGCO] - D[MGNI] + D[CONI]);
    c->quad[R4][R4] = -0.25*A[MGNI];
    c->quad[R4][R5] = -0.25*(f[NI] + ex[MGNI] + x[MGNI] + 2.0*w2ca[MG] - 2.0*w2ca[NI] + 2.0*w2[MGNI]);
    c->quad[R4][S1] = 0.25*(C[MGMN] - D[MGNI] + D[MNNI]);
    c->quad[R4][S2] = 0.25*(C[MGFE] - D[MGNI] + D[FENI]);
    c->quad[R4][S3] = 0.25*(-D[MGCO] + D[MGNI] - D[CONI]);
    c->quad[R4][S4] = 0.5*(-w1[MGNI] + w2[MGNI]);
    c->quad[R5][R5] = -1.0*w2ca[MG];
    c->quad[R5][S1] = 0.25*(f[MN] + ex[MGMN] + x[MGMN] - 2.0*w2ca[MG] + 2.0*w2ca[MN] - 2.0*w2[MGMN]);
    c->quad[R5][S2] = 0.25*(f[FE] + ex[MGFE] + x[MGFE] - 2.0*w2ca[MG] + 2.0*w2ca[FE] - 2.0*w2[MGFE]);
    c->quad[R5][S3] = -0.25*(f[CO] + ex[MGCO] + x[MGCO] - 2.0*w2ca[MG] + 2.0*w2ca[CO] - 2.0*w2[MGCO]);
    c->quad[R5][S4] = -0.25*(f[NI] + ex[MGNI] + x[MGNI] - 2.0*w2ca[MG] + 2.0*w2ca[NI] - 2.0*w2[MGNI]);
    c->quad[S1][S1] = 0.25*B[MGMN];
    c->quad[S1][S2] = 0.25*(B[MGMN] + B[MGFE] - B[MNFE]);
    c->quad[S1][S3] = 0.25*(-B[MGMN] - B[MGCO] + B[MNCO]);
    c->quad[S1][S4] = 0.25*(-B[MGMN] - B[MGNI] + B[MNNI]);
    c->quad[S2][S2] = 0.25*B[MGFE];
    c->quad[S2][S3] = 0.25*(-B[MGFE] - B[MGCO] + B[FECO]);
    c->quad[S2][S4] = 0.25*(-B[MGFE] - B[MGNI] + B[FENI]);
    c->quad[S3][S3] = 0.25*B[MGCO];
    c->quad[S3][S4] = 0.25*(B[MGCO] + B[MGNI] - B[CONI]);
    c->quad[S4][S4] = 0.25*B[MGNI];

    for (int i = 0; i < GH_OL_NZ; i++){
        for (int j = 0; j < i; j++){ c->quad[i][j] = c->quad[j][i]; }
    }
}

static void GH_ol_site_fracs(const double r[5], const double s[4], GH_ol_sf *f){
    static const double sgn[4] = {-1.0, -1.0, 1.0, 1.0};
    for (int i = 0; i < 4; i++){
        f->m1[i] = (r[i] + sgn[i]*s[i] + 1.0)/2.0;
        f->m2[i] = (r[i] - sgn[i]*s[i] + 1.0)/2.0;
    }
    f->ca2 = r[4];
    f->mg1 = 1.0 - f->m1[0] - f->m1[1] - f->m1[2] - f->m1[3];
    f->mg2 = 1.0 - f->m2[0] - f->m2[1] - f->m2[2] - f->m2[3] - f->ca2;
}

static void GH_ol_clamp_sf(GH_ol_sf *f, double eps){
    for (int i = 0; i < 4; i++){
        if (f->m1[i] <= eps) f->m1[i] = eps;
        if (f->m2[i] <= eps) f->m2[i] = eps;
    }
    if (f->ca2 <= eps) f->ca2 = eps;
    if (f->mg1 <= eps) f->mg1 = eps;
    if (f->mg2 <= eps) f->mg2 = eps;
}

static void GH_ol_pin_absent(GH_ol_sf *f, const int absent[4], double eps){
    for (int i = 0; i < 4; i++){
        if (absent[i]){ f->m1[i] = eps; f->m2[i] = eps; }
    }
}

static int GH_ol_sf_positive(const GH_ol_sf *f, const int absent[4], const double r[5]){
    for (int i = 0; i < 4; i++){
        if (!absent[i] && (!(f->m1[i] > 0.0) || !(f->m2[i] > 0.0))) return 0;
    }
    double mg_free = 1.0 - 0.5*(r[0]+1.0) - 0.5*(r[1]+1.0) - 0.5*(r[2]+1.0) - 0.5*(r[3]+1.0);
    if (r[4] > 0.0 && !(f->ca2 > 0.0)) return 0;
    if (mg_free > 0.0 && !(f->mg1 > 0.0)) return 0;
    if (2.0*mg_free - r[4] > 0.0 && !(f->mg2 > 0.0)) return 0;
    if (!(f->mg1 >= -1.0e-12) || !(f->mg2 >= -1.0e-12)) return 0;
    return 1;
}

static void GH_ol_pack_z(const double r[5], const double s[4], double z[GH_OL_NZ]){
    for (int i = 0; i < 5; i++){ z[i] = r[i]; }
    for (int i = 0; i < 4; i++){ z[5+i] = s[i]; }
}

static double GH_ol_entropy(const GH_ol_sf *f){
    return -GH_R_J*( f->m1[0]*log(f->m1[0]) + f->m2[0]*log(f->m2[0]) + f->m1[1]*log(f->m1[1]) + f->m2[1]*log(f->m2[1]) +
                     f->m1[2]*log(f->m1[2]) + f->m2[2]*log(f->m2[2]) + f->m1[3]*log(f->m1[3]) + f->m2[3]*log(f->m2[3]) +
                     f->ca2*log(f->ca2) + f->mg1*log(f->mg1) + f->mg2*log(f->mg2));
}

static double GH_ol_gibbs(const GH_ol_coef *c, const double z[GH_OL_NZ], const GH_ol_sf *f, double t){
    double h = c->g0;
    for (int i = 0; i < GH_OL_NZ; i++){ h += c->lin[i]*z[i]; }
    for (int i = 0; i < GH_OL_NZ; i++){
        for (int j = i; j < GH_OL_NZ; j++){ h += c->quad[i][j]*z[i]*z[j]; }
    }
    return h - t*GH_ol_entropy(f);
}

static double GH_ol_grad_poly(const GH_ol_coef *c, const double z[GH_OL_NZ], int i){
    double g = c->lin[i] + 2.0*c->quad[i][i]*z[i];
    for (int j = 0; j < GH_OL_NZ; j++){
        if (j != i){ g += c->quad[i][j]*z[j]; }
    }
    return g;
}

static void GH_ol_dg_dr(const GH_ol_coef *c, const double z[GH_OL_NZ], const GH_ol_sf *f, double t, double dr[5]){
    for (int i = 0; i < 4; i++){
        dr[i] = GH_ol_grad_poly(c, z, i) + 0.5*GH_R_J*t*(log(f->m1[i]*f->m2[i]/f->mg1/f->mg2));
    }
    dr[4] = GH_ol_grad_poly(c, z, 4) + GH_R_J*t*(log(f->ca2/f->mg2));
}

static double GH_ol_dg_ds(const GH_ol_coef *c, const double z[GH_OL_NZ], const GH_ol_sf *f, double t, int i){
    double ord = (i < 2) ? log(f->m2[i]*f->mg1/f->m1[i]/f->mg2) : log(f->m1[i]*f->mg2/f->m2[i]/f->mg1);
    return GH_ol_grad_poly(c, z, 5+i) + 0.5*GH_R_J*t*ord;
}

static double GH_ol_d2g_ds2(const GH_ol_coef *c, const GH_ol_sf *f, double t, int i, int j){
    if (i == j){
        return 2.0*c->quad[5+i][5+i] + 0.25*GH_R_J*t*(1.0/f->m1[i] + 1.0/f->m2[i] + 1.0/f->mg1 + 1.0/f->mg2);
    }
    double sg = ((i < 2) == (j < 2)) ? 1.0 : -1.0;
    return c->quad[5+i][5+j] + sg*0.25*GH_R_J*t*(1.0/f->mg1 + 1.0/f->mg2);
}

static void GH_ol_newton_terms(const GH_ol_coef *c, const double z[GH_OL_NZ], const GH_ol_sf *f, double t,
                               const int absent[4], double dgds[4], double H[4][4]){
    for (int i = 0; i < 4; i++){
        dgds[i] = absent[i] ? 0.0 : GH_ol_dg_ds(c, z, f, t, i);
        H[i][i] = absent[i] ? 1.0 : GH_ol_d2g_ds2(c, f, t, i, i);
        for (int j = i+1; j < 4; j++){
            H[i][j] = (!absent[i] && !absent[j]) ? GH_ol_d2g_ds2(c, f, t, i, j) : 0.0;
            H[j][i] = H[i][j];
        }
    }
}

static int GH_ol_invert4(double a[4][4]){
    int indxc[4], indxr[4], ipiv[4] = {0,0,0,0};
    int irow = 0, icol = 0;
    for (int i = 0; i < 4; i++){
        double big = 0.0;
        for (int j = 0; j < 4; j++){
            if (ipiv[j] != 1){
                for (int k = 0; k < 4; k++){
                    if (ipiv[k] == 0){
                        if (fabs(a[j][k]) >= big){ big = fabs(a[j][k]); irow = j; icol = k; }
                    }
                }
            }
        }
        ++(ipiv[icol]);
        if (irow != icol){
            for (int l = 0; l < 4; l++){ double tmp = a[irow][l]; a[irow][l] = a[icol][l]; a[icol][l] = tmp; }
        }
        indxr[i] = irow; indxc[i] = icol;
        if (a[icol][icol] == 0.0){ return 1; }
        double pivinv = 1.0/a[icol][icol];
        a[icol][icol] = 1.0;
        for (int l = 0; l < 4; l++){ a[icol][l] *= pivinv; }
        for (int ll = 0; ll < 4; ll++){
            if (ll != icol){
                double dum = a[ll][icol];
                a[ll][icol] = 0.0;
                for (int l = 0; l < 4; l++){ a[ll][l] -= a[icol][l]*dum; }
            }
        }
    }
    for (int l = 3; l >= 0; l--){
        if (indxr[l] != indxc[l]){
            for (int k = 0; k < 4; k++){ double tmp = a[k][indxr[l]]; a[k][indxr[l]] = a[k][indxc[l]]; a[k][indxc[l]] = tmp; }
        }
    }
    return 0;
}

static double GH_ol_G_pinned(const GH_ol_coef *c, double t, const double r[5], double u, double a){
    double s[4] = {u + a, u - a, 0.0, 0.0};
    GH_ol_sf f;
    double z[GH_OL_NZ];
    GH_ol_site_fracs(r, s, &f);
    GH_ol_clamp_sf(&f, DBL_EPSILON);
    GH_ol_pack_z(r, s, z);
    return GH_ol_gibbs(c, z, &f, t);
}
static void GH_ol_safe_order(const GH_ol_coef *c, double t, const double r[5], const int absent[4],
                             const double sMin[4], const double sMax[4], double sg[4]){
    const double eps = DBL_EPSILON;
    double s[4] = {0.0, 0.0, 0.0, 0.0};
    GH_ol_sf f0;
    GH_ol_site_fracs(r, s, &f0);
    if (!GH_ol_sf_positive(&f0, absent, r)){
        for (int i = 0; i < 2; i++){ if (!absent[i]){ s[i] = -0.5*r[4]; } }
    }
    double mg_free = 1.0 - 0.5*(r[0]+1.0) - 0.5*(r[1]+1.0) - 0.5*(r[2]+1.0) - 0.5*(r[3]+1.0);
    int pinned = (2.0*mg_free - r[4] <= 1.0e-12) && !absent[0] && !absent[1] && absent[2] && absent[3];
    if (pinned){
        double u  = s[0];
        double lo = fmax(fmax(-(r[0]+1.0) - u, u - (r[1]+1.0)), fmax(sMin[0] - u, u - sMax[1]));
        double hi = fmin(fmin( (r[0]+1.0) - u, u + (r[1]+1.0)), fmin(sMax[0] - u, u - sMin[1]));
        lo += 1.0e-12; hi -= 1.0e-12;
        const double gr = 0.6180339887498949;
        double a = lo, b = hi;
        double x1 = b - gr*(b - a), x2 = a + gr*(b - a);
        double f1 = GH_ol_G_pinned(c, t, r, u, x1), f2 = GH_ol_G_pinned(c, t, r, u, x2);
        for (int it = 0; it < 200 && (b - a) > 1.0e-13; it++){
            if (f1 < f2){ b = x2; x2 = x1; f2 = f1; x1 = b - gr*(b - a); f1 = GH_ol_G_pinned(c, t, r, u, x1); }
            else        { a = x1; x1 = x2; f1 = f2; x2 = a + gr*(b - a); f2 = GH_ol_G_pinned(c, t, r, u, x2); }
        }
        double am = 0.5*(a + b);
        sg[0] = u + am; sg[1] = u - am; sg[2] = 0.0; sg[3] = 0.0;
        return;
    }
    for (int iter = 0; iter < 500; iter++){
        GH_ol_sf f;
        double z[GH_OL_NZ], dgds[4], H[4][4];
        GH_ol_site_fracs(r, s, &f);
        GH_ol_clamp_sf(&f, eps);
        GH_ol_pack_z(r, s, z);
        GH_ol_newton_terms(c, z, &f, t, absent, dgds, H);
        double gmax = 0.0;
        for (int i = 0; i < 4; i++){ if (fabs(dgds[i]) > gmax){ gmax = fabs(dgds[i]); } }
        if (gmax < 1.0e-9){ break; }
        double ds[4] = {0.0, 0.0, 0.0, 0.0};
        GH_ol_invert4(H);
        for (int i = 0; i < 4; i++){
            for (int j = 0; j < 4; j++){ ds[i] -= H[i][j]*dgds[j]; }
            if (absent[i]){ ds[i] = 0.0; }
        }
        double lambda = 1.0;
        for (int ls = 0; ls < 60; ls++){
            double sn[4];
            GH_ol_sf fn;
            int okp = 1;
            for (int i = 0; i < 4; i++){
                sn[i] = s[i] + lambda*ds[i];
                if (!absent[i] && (sn[i] < sMin[i] || sn[i] > sMax[i])){ okp = 0; }
            }
            GH_ol_site_fracs(r, sn, &fn);
            if (okp && !GH_ol_sf_positive(&fn, absent, r)){ okp = 0; }
            if (okp){ break; }
            lambda *= 0.5;
        }
        double dmax = 0.0;
        for (int i = 0; i < 4; i++){
            s[i] += lambda*ds[i];
            if (fabs(lambda*ds[i]) > dmax){ dmax = fabs(lambda*ds[i]); }
        }
        if (dmax < 100.0*eps){ break; }
    }
    for (int i = 0; i < 4; i++){ sg[i] = s[i]; }
}

static void GH_ol_s_bounds(const double r[5], double sMin[4], double sMax[4]){
    sMax[0] = fmin(1.0-r[0], 1.0+r[0]);
    sMax[0] = fmin(sMax[0], 1.0-r[0]-2.0*r[4]);
    sMin[0] = fmax(r[0]-1.0, -r[0]-1.0);
    sMax[1] = fmin(1.0-r[1], 1.0+r[1]);
    sMax[1] = fmin(sMax[1], 1.0-r[1]-2.0*r[4]);
    sMin[1] = fmax(r[1]-1.0, -r[1]-1.0);
    sMax[2] = fmin(1.0-r[2], 1.0+r[2]);
    sMin[2] = fmax(r[2]-1.0, -r[2]-1.0);
    sMin[2] = fmax(sMin[2], r[2]-1.0-2.0*r[4]);
    sMax[3] = fmin(1.0-r[3], 1.0+r[3]);
    sMin[3] = fmax(r[3]-1.0, -r[3]-1.0);
    sMin[3] = fmax(sMin[3], r[3]-1.0-2.0*r[4]);
}

static void GH_ol_order(const GH_ol_coef *c, double t, const double r[5], double s_out[4], GH_ol_sf *sf){
    const double eps = DBL_EPSILON;
    GH_ol_sf f = {{0.0}, {0.0}, 0.0, 0.0, 0.0};
    double sOld[4], sNew[4], sMax[4], sMin[4];
    int iter = 0;

    for (int i = 0; i < 4; i++){ sOld[i] = 2.0; }

    GH_ol_s_bounds(r, sMin, sMax);

    sNew[0] = (r[0] > 0.0) ? sMax[0]/3.0  : 0.35*sMin[0]+0.65*sMax[0];
    sNew[1] = (r[1] > 0.0) ? sMax[1]/15.0 : 0.45*sMin[1]+0.55*sMax[1];
    sNew[2] = (r[2] > 0.0) ? sMax[2]/3.0  : 0.25*sMin[2]+0.75*sMax[2];
    sNew[3] = (r[3] > 0.0) ? sMax[3]/2.0  : 0.25*sMin[3]+0.75*sMax[3];

    int absent[4];
    for (int i = 0; i < 4; i++){ absent[i] = (fabs(r[i]+1.0) < 10.0*eps); }

    while ((fabs(sNew[0]-sOld[0]) > 100.0*eps || fabs(sNew[1]-sOld[1]) > 100.0*eps ||
            fabs(sNew[2]-sOld[2]) > 100.0*eps || fabs(sNew[3]-sOld[3]) > 100.0*eps) && iter < 1000){
        double s[4] = {sNew[0], sNew[1], sNew[2], sNew[3]};
        double z[GH_OL_NZ], dgds[4], d2gds2[4][4];

        for (int i = 0; i < 4; i++){ if (r[i] == -1.0) s[i] = 0.0; }

        GH_ol_site_fracs(r, s, &f);
        GH_ol_clamp_sf(&f, eps);
        GH_ol_pack_z(r, s, z);
        GH_ol_newton_terms(c, z, &f, t, absent, dgds, d2gds2);
        GH_ol_pin_absent(&f, absent, eps);

        for (int i = 0; i < 4; i++){ sOld[i] = s[i]; }

        GH_ol_invert4(d2gds2);

        for (int i = 0; i < 4; i++){ if (absent[i]) d2gds2[i][i] = 0.0; }

        for (int i = 0; i < 4; i++){
            for (int j = 0; j < 4; j++){ s[i] += -d2gds2[i][j]*dgds[j]; }
        }
        for (int i = 0; i < 4; i++){
            if (!absent[i]){ s[i] = fmax(s[i], sMin[i]+3.0*eps); s[i] = fmin(s[i], sMax[i]-3.0*eps); }
        }

        for (int i = 0; i < 4; i++){ sNew[i] = s[i]; }
        iter++;
    }

    GH_ol_sf raw;
    GH_ol_site_fracs(r, sOld, &raw);
    if (iter >= 1000 || !GH_ol_sf_positive(&raw, absent, r)){
        double sg[4] = {0.0, 0.0, 0.0, 0.0};
        GH_ol_safe_order(c, t, r, absent, sMin, sMax, sg);
        for (int i = 0; i < 4; i++){ sOld[i] = sg[i]; }
        GH_ol_site_fracs(r, sg, &f);
        GH_ol_clamp_sf(&f, eps);
        GH_ol_pin_absent(&f, absent, eps);
    }
    for (int i = 0; i < 4; i++){ s_out[i] = sOld[i]; }
    *sf = f;
}

static double GH_ol_gmix(double t, double p, const double r[5], double dx[5]){
    GH_ol_coef c;
    GH_ol_sf f;
    double s[4], z[GH_OL_NZ];
    GH_ol_coef_at(p, &c);
    GH_ol_order(&c, t, r, s, &f);
    GH_ol_pack_z(r, s, z);
    if (dx){ GH_ol_dg_dr(&c, z, &f, t, dx); }
    return GH_ol_gibbs(&c, z, &f, t);
}

double obj_gh_ol(unsigned n, const double *x, double *grad, void *SS_ref_db){
    SS_ref *d = (SS_ref *) SS_ref_db;

    double T    = d->T;
    double Pbar = d->P*1000.0;
    double *p   = d->p;
    double *gb  = d->gb_lvl;
    double *mu_Gex = d->mu_Gex;

    for (int i = 0; i < 4; i++){ p[i] = x[i]; }
    double sum = p[0] + p[1] + p[2] + p[3];
    double r[5] = { 2.0*p[2]/sum - 1.0, 2.0*p[1]/sum - 1.0, -1.0, -1.0, p[3]/sum };
    for (int i = 0; i < 2; i++){ if (fabs(r[i] + 1.0) < 1.0e-12){ r[i] = -1.0; } }

    double dx[5];
    double Gex = GH_ol_gmix(T, Pbar, r, dx);

    double dGex[4];
    dGex[0] = 0.0;
    dGex[1] = 2.0*dx[1];
    dGex[2] = 2.0*dx[0];
    dGex[3] = dx[4];

    GH_mu_from_dGex(d, p, sum, 4, Gex, dGex);

    GH_sf_ol(p, d->sf);
    return GH_df_assemble(d, p, gb, mu_Gex, 4, grad);
}

double obj_gh_fluid(unsigned n, const double *x, double *grad, void *SS_ref_db){
    SS_ref *d = (SS_ref *) SS_ref_db;

    double T    = d->T;
    double Pbar = d->P*1000.0;
    double *p   = d->p;
    double *gb  = d->gb_lvl;
    double *mu_Gex = d->mu_Gex;

    for (int i = 0; i < 2; i++){ p[i] = x[i]; }
    double x0 = p[0], x1 = p[1];

    double muW, muC;
    GH_duan_mix_muGex(x0, T, Pbar, &muW, &muC);
    mu_Gex[0] = muW/1000.0;
    mu_Gex[1] = muC/1000.0;
    d->sf[0] = p[0]; d->sf[1] = p[1];

    return GH_df_assemble(d, p, gb, mu_Gex, 2, grad);
}

double obj_gh_bi(unsigned n, const double *x, double *grad, void *SS_ref_db){
    SS_ref *d = (SS_ref *) SS_ref_db;

    int n_em    = d->n_em;
    double R    = GH_R_KJ;
    double T    = d->T;
    double *p   = d->p;
    double *gb  = d->gb_lvl;
    double *mu_Gex = d->mu_Gex;

    for (int i = 0; i < n_em; i++){
        p[i] = x[i];
    }

    double Gex_n2[n_em];
    gex_sym_n2(d, p, n_em, Gex_n2);
    for (int i = 0; i < n_em; i++){
        double Gex = Gex_n2[i];
        mu_Gex[i] = Gex/1000.0;
        d->sf[i]  = p[i];
    }

    double dSi[2];
    double Sconfig = GH_config_entropy(d, p, 2, 1.0, R, T, dSi);
    return GH_df_assemble_cfg(d, p, gb, mu_Gex, n_em, grad, Sconfig, dSi);
}

double obj_gh_fsp(unsigned n, const double *x, double *grad, void *SS_ref_db){
    SS_ref *d = (SS_ref *) SS_ref_db;

    double R = GH_R_KJ;
    double T = d->T;
    double Pbar = d->P*1000.0;
    double *p   = d->p;
    double *gb  = d->gb_lvl;
    double *mu_Gex = d->mu_Gex;

    for (int i = 0; i < 3; i++){ p[i] = x[i]; }

    const GH_hsv Wfsp[6] = {
        { 0.0,     0.0,  0.0     },
        { 7924.0,  0.0,  0.0     },
        { 27320.0, 10.3, 0.3264  },
        { 18810.0, 10.3, 0.4602  },
        { 40317.0, 0.0,  0.0     },
        { 38974.0, 0.0, -0.1037  }
    };
    const GH_hsv Wfsp123 = { 12545.0, 0.0, -1.095 };

    double Gex, dGex[3];
    GH_ternary_cubic(Wfsp, Wfsp123, T, Pbar-1.0, p, &Gex, dGex);

    GH_mu_from_dGex(d, p, 1.0, 3, Gex, dGex);

    double dSi[3];
    double Sconfig = GH_config_entropy(d, p, 3, 1.0, R, T, dSi);
    return GH_df_assemble_cfg(d, p, gb, mu_Gex, 3, grad, Sconfig, dSi);
}

double obj_gh_g(unsigned n, const double *x, double *grad, void *SS_ref_db){
    SS_ref *d = (SS_ref *) SS_ref_db;

    double R = GH_R_KJ;
    double T = d->T;
    double Pbar = d->P*1000.0;
    double *p   = d->p;
    double *gb  = d->gb_lvl;
    double *mu_Gex = d->mu_Gex;

    for (int i = 0; i < 3; i++){ p[i] = x[i]; }

    const GH_hsv Wg[6] = {
        { 21560.0, 18.79, 0.10 },
        { 69200.0, 18.79, 0.10 },
        { 20320.0, 5.08,  0.17 },
        { 2620.0,  5.08,  0.09 },
        { 230.0,   0.0,   0.01 },
        { 3720.0,  0.0,   0.06 }
    };
    const GH_hsv Wg123 = { 0.0, 0.0, 0.0 };

    double Gex, dGex[3];
    GH_ternary_cubic(Wg, Wg123, T, Pbar, p, &Gex, dGex);

    GH_mu_from_dGex(d, p, 1.0, 3, Gex, dGex);

    double dSi[3];
    double Sconfig = GH_config_entropy(d, p, 3, 3.0, R, T, dSi);
    return GH_df_assemble_cfg(d, p, gb, mu_Gex, 3, grad, Sconfig, dSi);
}

double obj_gh_hb(unsigned n, const double *x, double *grad, void *SS_ref_db){
    SS_ref *d = (SS_ref *) SS_ref_db;

    double T    = d->T;
    double Rgas = GH_R_J;
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
    double S  = -Rgas*(4.0*xMgM12*rlog((xMgM12)) + 4.0*xFe2M12*rlog((xFe2M12)) + xFe3M3*rlog((xFe3M3)) + xAlM3*rlog((xAlM3)));
    double Gex = H - T*S;

    double dGex[3];
    dGex[0] = 0.0;
    dGex[1] = WFEMG*(1.0-2.0*p[1]) + 4.0*Rgas*T*rlog((xFe2M12/xMgM12));
    dGex[2] = WFEAL*(1.0-2.0*p[2]) +     Rgas*T*rlog((xFe3M3/xAlM3));

    GH_mu_from_dGex(d, p, 1.0, 3, Gex, dGex);

    GH_sf_hb(p, d->sf);
    return GH_df_assemble(d, p, gb, mu_Gex, 3, grad);
}

double obj_gh_lc(unsigned n, const double *x, double *grad, void *SS_ref_db){
    SS_ref *d = (SS_ref *) SS_ref_db;

    double T    = d->T;
    double Rgas = GH_R_J;
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

    double K = 1.5*rlog((2.0/3.0)) + 0.5*rlog((0.5));

    double H = WNAK*p0 + WNAH*p1 - WNAK*p0*p0 - WNAH*p1*p1 + DGR*p0*p1;
    double S = -Rgas*( a*rlog((a)) + b*rlog((b)) + c*rlog((c))
                      + 1.5*(dd*rlog((dd)) + e*rlog((e)) + f*rlog((f)))
                      - p1*K );

    double Gex = H - T*S;

    double dH0 = WNAK - 2.0*WNAK*p0 + DGR*p1;
    double dH1 = WNAH - 2.0*WNAH*p1 + DGR*p0;

    double dS0 = -Rgas*( (1.0-p1)*rlog((a/b)) + p1*rlog((dd/e)) );
    double dS1 = -Rgas*( p0*rlog((dd/a)) + (1.0-p0)*rlog((e/b)) + rlog((c/f)) - K );

    double dGex[3];
    dGex[0] = dH0 - T*dS0;
    dGex[1] = dH1 - T*dS1;
    dGex[2] = 0.0;

    GH_mu_from_dGex(d, p, 1.0, 3, Gex, dGex);

    GH_sf_of_p("lc", p, d->sf);
    return GH_df_assemble(d, p, gb, mu_Gex, 3, grad);
}

static double GH_mel_mix_geh_order_G(double T){
    const double DG22p=12000.0, W22p=38354.0;
    double R = 8.3143;
    double s = 0.1;
    for (int iter = 0; iter < 1000; iter++){
        double dgds   = R*T*(-2.0*rlog((1.0-s)) + rlog((s)) + rlog((1.0+s))) + DG22p + W22p*(1.0-2.0*s);
        double d2gds2 = R*T*(2.0/(1.0-s) + 1.0/s + 1.0/(1.0+s)) - 2.0*W22p;
        double sNew   = s - dgds/d2gds2;
        sNew = (sNew > 1.0 - 2.220446049250313e-16) ? 1.0 - 2.220446049250313e-16 : sNew;
        sNew = (sNew < 2.220446049250313e-16) ? 2.220446049250313e-16 : sNew;
        int converged = (fabs(sNew - s) <= 10.0*2.220446049250313e-16);
        s = sNew;
        if (converged) break;
    }
    double S = -R*( 2.0*(1.0-s)*rlog((1.0-s)) + s*rlog((s)) + (1.0+s)*rlog((1.0+s)) - 2.0*rlog((2.0)) );
    double H = DG22p*s + W22p*s*(1.0-s);
    return H - T*S;
}

typedef struct { double gs, gr0s, gr1s, gr2s, gr0r1, gr0r2, gr1r2, w22p, w12, w13, w14; } GH_mel_gc;

static GH_mel_gc GH_mel_gcoef(void){
    const double DG22p=12000.0, W12=9000.0, W12p=13000.0, W13=16000.0, W14=0.0,
                 W22p=38354.0,  W23=9000.0, W2p3=13000.0, W24=9000.0,  W2p4=13000.0, W34=0.0;
    GH_mel_gc c;
    c.gs    = DG22p+W12p-W12;
    c.gr0s  = W22p-W12p+W12;
    c.gr1s  = W2p3-W23-W12p+W12;
    c.gr2s  = W2p4-W24-W12p+W12;
    c.gr0r1 = W23-W12-W13;
    c.gr0r2 = W24-W14-W12;
    c.gr1r2 = W34-W14-W13;
    c.w22p = W22p; c.w12 = W12; c.w13 = W13; c.w14 = W14;
    return c;
}

static double GH_mel_dgds(const GH_mel_gc *c, double T, double r0, double r1, double r2, double s){
    double xAl3T1 = r0-s, xSi4T1 = r2+s, xAl3T2 = 0.5*(r0+s), xSi4T2 = 1.0-0.5*(r0+s);
    return GH_R_J*T*(-rlog((xAl3T1))+rlog((xSi4T1))+rlog((xAl3T2))-rlog((xSi4T2)))
         + c->gs - 2.0*c->w22p*s + c->gr0s*r0 + c->gr1s*r1 + c->gr2s*r2;
}

double obj_gh_mel(unsigned n, const double *x, double *grad, void *SS_ref_db){
    SS_ref *d = (SS_ref *) SS_ref_db;

    double T    = d->T;
    double Rgas = GH_R_J;
    double *p   = d->p;
    double *gb  = d->gb_lvl;
    double *mu_Gex = d->mu_Gex;

    for (int i = 0; i < 4; i++){ p[i] = x[i]; }
    double r0 = p[1], r1 = p[2], r2 = p[3];
    double xMg2T1 = p[0] + d->d_em[0];

    GH_mel_gc c = GH_mel_gcoef();

    double eps_s = 1.0e-10;
    double s_lo = fmax(fmax(r0-(1.0-eps_s), eps_s-r2), 2.0*eps_s-r0);
    double s_hi = fmin(fmin(r0-eps_s, (1.0-eps_s)-r2), 2.0*(1.0-eps_s)-r0);

    double dgds_lo = GH_mel_dgds(&c, T, r0, r1, r2, s_lo);
    double dgds_hi = GH_mel_dgds(&c, T, r0, r1, r2, s_hi);

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
            double fm = GH_mel_dgds(&c, T, r0, r1, r2, mid);
            if ((fa > 0.0) == (fm > 0.0)){ a = mid; fa = fm; }
            else { b = mid; }
            if (b-a < 1.0e-14){ break; }
        }
        s = 0.5*(a+b);
    }

    double xCaOct = 1.0-r2, xNaOct = r2 + d->d_em[3], xFe2T1 = r1 + d->d_em[2];
    double xAl3T1 = r0-s, xSi4T1 = r2+s, xAl3T2 = 0.5*(r0+s), xSi4T2 = 1.0-0.5*(r0+s);

    double H = c.gs*s + c.w12*r0*(1.0-r0) - c.w22p*s*s + c.w13*r1*(1.0-r1)
             + c.gr0s*r0*s + c.gr0r1*r0*r1 + c.gr1s*r1*s
             + c.w14*r2*(1.0-r2) + c.gr0r2*r0*r2 + c.gr1r2*r1*r2
             + c.gr2s*r2*s;

    double S = -Rgas*( 2.0*xCaOct*rlog((xCaOct)) + 2.0*xNaOct*rlog((xNaOct))
                      + xMg2T1*rlog((xMg2T1)) + xFe2T1*rlog((xFe2T1))
                      + xAl3T1*rlog((xAl3T1)) + xSi4T1*rlog((xSi4T1))
                      + 2.0*xAl3T2*rlog((xAl3T2)) + 2.0*xSi4T2*rlog((xSi4T2)) );

    double Graw = H - T*S;

    double geh_G = GH_mel_mix_geh_order_G(T);
    double Gex = Graw - p[1]*geh_G;

    double dGdr0 = Rgas*T*(rlog((xAl3T1))+rlog((xAl3T2))-rlog((xSi4T2))+1.0)
                 + c.w12*(1.0-2.0*r0) + c.gr0s*s
                 + c.gr0r1*r1 + c.gr0r2*r2;
    double dGdr1 = Rgas*T*(rlog((xFe2T1))+1.0)
                 + c.w13*(1.0-2.0*r1) + c.gr0r1*r0
                 + c.gr1s*s + c.gr1r2*r2;
    double dGdr2 = Rgas*T*(-2.0*rlog((xCaOct))+2.0*rlog((xNaOct))+rlog((xSi4T1))+1.0)
                 + c.w14*(1.0-2.0*r2) + c.gr0r2*r0
                 + c.gr1r2*r1 + c.gr2s*s;

    double dGex[4];
    dGex[0] = Rgas*T*(rlog((xMg2T1))+1.0);
    dGex[1] = dGdr0 - geh_G;
    dGex[2] = dGdr1;
    dGex[3] = dGdr2;

    GH_mu_from_dGex(d, p, 1.0, 4, Gex, dGex);

    return GH_df_assemble(d, p, gb, mu_Gex, 4, grad);
}

double obj_gh_cum(unsigned n, const double *x, double *grad, void *SS_ref_db){
    SS_ref *d = (SS_ref *) SS_ref_db;

    double T    = d->T;
    double Rgas = GH_R_J;
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

    double eps_s = 1.0e-8;
    double s0 = 0.0, s1 = 0.0;
    for (int iter = 0; iter < 200; iter++){
        double xfe2m4  = p1+3.0*s0/7.0+2.0*s1/7.0;
        double xfe2m13 = p1-4.0*s0/7.0+2.0*s1/7.0;
        double xfe2m2  = p1+3.0*s0/7.0-5.0*s1/7.0;
        double xmg2m4 = 1.0-xfe2m4, xmg2m13 = 1.0-xfe2m13, xmg2m2 = 1.0-xfe2m2;

        double dgds0 = Rgas*T*( (6.0/7.0)*rlog((xfe2m4/xmg2m4)) - (12.0/7.0)*rlog((xfe2m13/xmg2m13)) + (6.0/7.0)*rlog((xfe2m2/xmg2m2)) )
                     + HS13 + HRS13*r + 2.0*HS13S13*s0 + HS13S2*s1;
        double dgds1 = Rgas*T*( (4.0/7.0)*rlog((xfe2m4/xmg2m4)) + (6.0/7.0)*rlog((xfe2m13/xmg2m13)) - (10.0/7.0)*rlog((xfe2m2/xmg2m2)) )
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

    double S = -Rgas*( 2.0*xfe2m4*rlog((xfe2m4)) + 2.0*xmg2m4*rlog((xmg2m4))
                      + 3.0*xfe2m13*rlog((xfe2m13)) + 3.0*xmg2m13*rlog((xmg2m13))
                      + 2.0*xfe2m2*rlog((xfe2m2)) + 2.0*xmg2m2*rlog((xmg2m2)) );

    double Gex = H - T*S;

    double dGdp1 = 4.0*H0c*(p0-p1) + 2.0*HRS13*s0 + 2.0*HRS2*s1
                 + Rgas*T*( 2.0*rlog((xfe2m4/xmg2m4)) + 3.0*rlog((xfe2m13/xmg2m13)) + 2.0*rlog((xfe2m2/xmg2m2)) );

    double dGex[2];
    dGex[0] = 0.0;
    dGex[1] = dGdp1;

    GH_mu_from_dGex(d, p, 1.0, 2, Gex, dGex);

    return GH_df_assemble(d, p, gb, mu_Gex, 2, grad);
}

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

static inline double GH_spn_pos(double v){ return (v <= 0.0) ? 2.220446049250313e-16 : v; }

static void GH_spn_site_fracs(double p0,double p2,double p3,double p4,double s0,double s1,double s2,double d_em2,
    double *xmg2tet,double *xfe2tet,double *xal3tet,double *xfe3tet,
    double *xmg2oct,double *xfe2oct,double *xal3oct,double *xfe3oct){
    double eps = 2.220446049250313e-16;
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
           + 0.5*Rgas*T*rlog(((xmg2tet*xfe2oct)/(xfe2tet*xmg2oct)));
    *dgds1 = gs2 + gx2s2*p3 + gx3s2*p0 + gx4s2*p4 + gx5s2*p2
           + gs1s2*s0 + 2.0*gs2s2*s1 + gs2s3*p0 + gs2s4*s2
           + Rgas*T*rlog(((xfe2tet*xal3oct)/(xal3tet*xfe2oct)));
    *dgds2 = hs4 - T*ss4 + gx2s4*p3 + gx3s4*p0 + gx4s4*p4 + gx5s4*p2
           + gs1s4*s0 + gs2s4*s1 + gs3s4*p0 + 2.0*gs4s4*s2
           + Rgas*T*rlog(((xfe2tet*xfe3oct)/(xfe3tet*xfe2oct)));
}

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
    double xmg2tet=GH_spn_pos((p3+s0)/2.0), xfe2tet=GH_spn_pos(p4-0.5*p3-0.5*s0+s1+p0+s2), xal3tet=GH_spn_pos(1.0-p0-p4-p2-s1), xfe3tet=GH_spn_pos(p2-s2+d_em2);
    double xmg2oct=GH_spn_pos((p3-s0)/4.0), xfe2oct=GH_spn_pos((2.0-p3+s0-2.0*s1-2.0*p0-2.0*s2)/4.0), xal3oct=GH_spn_pos((1.0-p0-p4-p2+s1)/2.0), xfe3oct=GH_spn_pos((p2+s2)/2.0+d_em2);
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

    double S = -Rgas*( xmg2tet*rlog((xmg2tet)) + xfe2tet*rlog((xfe2tet))
                      + xal3tet*rlog((xal3tet)) + xfe3tet*rlog((xfe3tet))
                      + 2.0*xmg2oct*rlog((xmg2oct)) + 2.0*xfe2oct*rlog((xfe2oct))
                      + 2.0*xal3oct*rlog((xal3oct)) + 2.0*xfe3oct*rlog((xfe3oct))
                      + 2.0*xcr3oct*rlog((xcr3oct)) + 2.0*xti4oct*rlog((xti4oct)) ) + gc[4]*s2;

    double Gvol = p4*p2*(WV1*p2+WV2*p4)*(Pval-1.0);
    return H - T*S + Gvol;
}

static void GH_spn_solve_s_from(double p0,double p2,double p3,double p4,double T,double Rgas,const double *gc,double d_em2,
                                 double s0init,double s1init,double s2init,double *s0o,double *s1o,double *s2o){
    const int maxIter = 200;
    double s0=s0init, s1=s1init, s2=s2init;
    double s0Old=2.0, s1Old=2.0, s2Old=2.0;
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

        double r0=-dgds0, r1=-dgds1, r2=-dgds2;
        double det  = d00*(d11*d22-d12*d12) - d01*(d01*d22-d12*d02) + d02*(d01*d12-d11*d02);
        double det0 =  r0*(d11*d22-d12*d12) - d01*( r1*d22-d12* r2) + d02*( r1*d12-d11* r2);
        double det1 = d00*( r1*d22- r2*d12) -  r0*(d01*d22-d12*d02) + d02*(d01* r2- r1*d02);
        double det2 = d00*(d11* r2-d12* r1) - d01*(d01* r2- r1*d02) +  r0*(d01*d12-d11*d02);
        double deltaS0 = det0/det, deltaS1 = det1/det, deltaS2 = det2/det;
        if (!isfinite(deltaS0) || !isfinite(deltaS1) || !isfinite(deltaS2)){ break; }

        double lambda = 1.0;
        double val, slope;

        val=xmg2tet; slope=deltaS0/2.0;
        if      (slope != 0.0 && val+lambda*slope < 0.0) lambda = -val/slope;
        else if (slope != 0.0 && val+lambda*slope > 1.0) lambda = (1.0-val)/slope;

        val=xfe2tet; slope=-deltaS0/2.0+deltaS1+deltaS2;
        if      (slope != 0.0 && val+lambda*slope < 0.0) lambda = -val/slope;
        else if (slope != 0.0 && val+lambda*slope > 1.0) lambda = (1.0-val)/slope;

        val=xal3tet; slope=-deltaS1;
        if      (slope != 0.0 && val+lambda*slope < 0.0) lambda = -val/slope;
        else if (slope != 0.0 && val+lambda*slope > 1.0) lambda = (1.0-val)/slope;

        val=xfe3tet; slope=-deltaS2;
        if      (slope != 0.0 && val+lambda*slope < 0.0) lambda = -val/slope;
        else if (slope != 0.0 && val+lambda*slope > 1.0) lambda = (1.0-val)/slope;

        val=xmg2oct; slope=-deltaS0/4.0;
        if      (slope != 0.0 && val+lambda*slope < 0.0) lambda = -val/slope;
        else if (slope != 0.0 && val+lambda*slope > 1.0) lambda = (1.0-val)/slope;

        val=xfe2oct; slope=deltaS0/4.0-deltaS1/2.0-deltaS2/2.0;
        if      (slope != 0.0 && val+lambda*slope < 0.0) lambda = -val/slope;
        else if (slope != 0.0 && val+lambda*slope > 1.0) lambda = (1.0-val)/slope;

        val=xal3oct; slope=deltaS1/2.0;
        if      (slope != 0.0 && val+lambda*slope < 0.0) lambda = -val/slope;
        else if (slope != 0.0 && val+lambda*slope > 1.0) lambda = (1.0-val)/slope;

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

static void GH_spn_gcoef(double gc[31], double gx[10]){
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

    gc[0]  = 0.5*(-WOCT+W24U-W14+0.5*HX+0.5*HEX+H24);
    gc[1]  = W11+H11;
    gc[2]  = W13P-W13+H33;
    gc[3]  = W15P-W15+H55;
    gc[4]  = S55;
    gc[5]  = 0.5*(WOCT-WTET);
    gc[6]  = 0.5*(WTET-W22+W11-WOCT+2.0*W2P4U-2.0*W1P4-HEX+2.0*H24);
    gc[7]  = 0.5*(WTET-WOCT+2.0*W3PU4U-2.0*W3P4-W2P3U-W23PU+W22+W1P3+W13P-W11-HEX+2.0*H24-2.0*H23);
    gc[8]  = 0.5*(WTET-WOCT-W11+W22+2.0*W4U5PU-2.0*W45P-W2P5U-W25PU+W1P5+W15P-HEX-2.0*H25+2.0*H24);
    gc[9]  = 0.5*(W2P3U-W22-W1P3+W11);
    gc[10] = W1P3-W13-W11;
    gc[11] = W33-W13P+W13;
    gc[12] = W35P-W35-W15P+W15;
    gc[13] = 0.5*(WTET-W24U+W14-0.5*HX+0.5*HEX-H24);
    gc[14] = -W11+W1P4-W14;
    gc[15] = W3P4-W34-W13P+W13;
    gc[16] = W45P-W45-W15P+W15;
    gc[17] = 0.5*(-W22+W11+W2P5U-W1P5);
    gc[18] = -W11+W1P5-W15;
    gc[19] = W3P5-W35-W13P+W13;
    gc[20] = W55-W15P+W15;
    gc[21] = 0.25*(-WTET-WOCT+HX);
    gc[22] = 0.5*(WTET-W22+W11+WOCT-HX);
    gc[23] = 0.5*(WTET+WOCT-W2P3U-W23PU+W22+W1P3+W13P-W11-HX);
    gc[24] = 0.5*(WTET+WOCT+W22-W11-W2P5U-W25PU+W1P5+W15P-HX);
    gc[25] = -W11;
    gc[26] = W1P3P-W1P3-W13P+W13;
    gc[27] = W1P5P-W1P5-W15P+W15;
    gc[28] = -W33;
    gc[29] = W3P5P-W3P5-W35P+W35;
    gc[30] = -W55;

    gx[0] = gx2x2; gx[1] = gx2x3; gx[2] = gx2x4; gx[3] = gx2x5; gx[4] = gx3x3;
    gx[5] = gx3x4; gx[6] = gx3x5; gx[7] = gx4x4; gx[8] = gx4x5; gx[9] = gx5x5;
}

static void GH_spn_mix_endmember_G(double T, double *CR_G, double *HC_G, double *MT_G, double *SP_G, double *UV_G){
    const double R = 8.3143;
    const double eps = 2.220446049250313e-16;

    double gc[31], gx[10];
    GH_spn_gcoef(gc, gx);
    const double gs1 = gc[0], gs2 = gc[1], gs3 = gc[2], hs4 = gc[3], S55 = gc[4];
    const double gx2x2 = gx[0], gx3x3 = gx[4], gx4x4 = gx[7], gx5x5 = gx[9];
    const double gx2s1 = gc[5], gx2s2 = gc[6], gx3s3 = gc[11], gx5s4 = gc[20];
    const double gs1s1 = gc[21], gs1s2 = gc[22], gs2s2 = gc[25], gs3s3 = gc[28], gs4s4 = gc[30];

    double s0=0.5, s1=0.9, s2=0.1;
    for (int iter = 0; iter < 1000; iter++){
        double dgds0 = gs1 + gs2/2.0 + gx2s1 + gx2s2/2.0 + 2.0*gs1s1*s0
                     + gs1s2*(0.5+s0) + gs2s2*(1.0+s0)/2.0
                     + R*T*(0.5*rlog((1.0+s0)) + 0.5*rlog((3.0+s0)) - rlog((1.0-s0)));
        double dgds1 = gs2 + 2.0*gs2s2*s1
                     + R*T*(rlog((s1)) + rlog((1.0+s1)) - 2.0*rlog((1.0-s1)));
        double dgds2 = hs4 - T*S55 + gx5s4 + 2.0*gs4s4*s2
                     + R*T*(rlog((s2)) - 2.0*rlog((1.0-s2)) + rlog((1.0+s2)));
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

    double SP_S = -R*(0.5*(1.0+s0)*rlog((1.0+s0))+(1.0-s0)*rlog((1.0-s0))+0.5*(3.0+s0)*rlog((3.0+s0))-5.0*rlog((2.0)));
    double SP_H = gs1*s0 + gs2*(1.0+s0)/2.0 + gx2x2 + gx2s1*s0 + gx2s2*(1.0+s0)/2.0
                + gs1s1*s0*s0 + gs1s2*s0*(1.0+s0)/2.0 + gs2s2*(1.0+s0)*(1.0+s0)/4.0;
    *SP_G = SP_H - T*SP_S;

    double HC_S = -R*(s1*rlog((s1))+2.0*(1.0-s1)*rlog((1.0-s1))+(1.0+s1)*rlog((1.0+s1))-2.0*rlog((2.0)));
    double HC_H = gs2*s1 + gs2s2*s1*s1;
    *HC_G = HC_H - T*HC_S;

    *CR_G = gs3 + gx3x3 + gx3s3 + gs3s3;

    *UV_G = gx4x4 - T*(R*2.0*rlog((2.0)));

    double MT_S = -R*(s2*rlog((s2))+2.0*(1.0-s2)*rlog((1.0-s2))+(1.0+s2)*rlog((1.0+s2))-2.0*rlog((2.0))) + S55*s2;
    double MT_H = hs4*s2 + gx5x5 + gx5s4*s2 + gs4s4*s2*s2;
    *MT_G = MT_H - T*MT_S;
}

double obj_gh_spn(unsigned n, const double *x, double *grad, void *SS_ref_db){
    SS_ref *d = (SS_ref *) SS_ref_db;

    double T    = d->T;
    double Rgas = GH_R_J;
    double *p   = d->p;
    double *gb  = d->gb_lvl;
    double *mu_Gex = d->mu_Gex;

    for (int i = 0; i < 5; i++){ p[i] = x[i]; }
    double p0 = p[0], p2 = p[2], p3 = p[3], p4 = p[4];

    const double WV1=-0.1250, WV2=0.1018;

    double gc[31], gx[10];
    GH_spn_gcoef(gc, gx);
    const double gx2x2 = gx[0], gx2x3 = gx[1], gx2x4 = gx[2], gx2x5 = gx[3], gx3x3 = gx[4],
                 gx3x4 = gx[5], gx3x5 = gx[6], gx4x4 = gx[7], gx4x5 = gx[8], gx5x5 = gx[9];

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
    if (!d->gh_multistart_order){
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
                                    WV1,WV2,d->P*1000.0,d->d_em[0],d->d_em[2],d->d_em[4]);
            if (G < bestG){ bestG = G; bs0 = cs0; bs1 = cs1; bs2 = cs2; }
        }
        s0 = bs0; s1 = bs1; s2 = bs2;
    }

    double xmg2tet=GH_spn_pos((p3+s0)/2.0), xfe2tet=GH_spn_pos(p4-0.5*p3-0.5*s0+s1+p0+s2), xal3tet=GH_spn_pos(1.0-p0-p4-p2-s1), xfe3tet=GH_spn_pos(p2-s2+d->d_em[2]);
    double xmg2oct=GH_spn_pos((p3-s0)/4.0), xfe2oct=GH_spn_pos((2.0-p3+s0-2.0*s1-2.0*p0-2.0*s2)/4.0), xal3oct=GH_spn_pos((1.0-p0-p4-p2+s1)/2.0), xfe3oct=GH_spn_pos((p2+s2)/2.0+d->d_em[2]);
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

    double S = -Rgas*( xmg2tet*rlog((xmg2tet)) + xfe2tet*rlog((xfe2tet))
                      + xal3tet*rlog((xal3tet)) + xfe3tet*rlog((xfe3tet))
                      + 2.0*xmg2oct*rlog((xmg2oct)) + 2.0*xfe2oct*rlog((xfe2oct))
                      + 2.0*xal3oct*rlog((xal3oct)) + 2.0*xfe3oct*rlog((xfe3oct))
                      + 2.0*xcr3oct*rlog((xcr3oct)) + 2.0*xti4oct*rlog((xti4oct)) ) + gc[4]*s2;

    double Gvol = p4*p2*(WV1*p2+WV2*p4)*((d->P*1000.0)-1.0);
    double Graw = H - T*S + Gvol;

    double CR_G, HC_G, MT_G, SP_G, UV_G;
    GH_spn_mix_endmember_G(T, &CR_G, &HC_G, &MT_G, &SP_G, &UV_G);
    double ENDMEMBERS = p0*CR_G + p[1]*HC_G + p2*MT_G + p3*SP_G + p4*UV_G;
    double Gex = Graw - ENDMEMBERS;

    double dGdp0 = gc[2] + gc[7]*p3 + gx2x3*p3 + 2.0*gx3x3*p0 + gx3x4*p4 + gx3x5*p2 + gc[9]*s0 + gc[10]*s1 + 2.0*gc[11]*p0 + gc[12]*s2
                 + gc[15]*p4 + gc[19]*p2 + gc[23]*s0 + gc[26]*s1 + 2.0*gc[28]*p0 + gc[29]*s2
                 + Rgas*T*(rlog((xfe2tet/xal3tet)) + 2.0*rlog((xcr3oct)) - rlog((xfe2oct)) - rlog((xal3oct)));

    double dGdp2 = gx5x5*2.0*p2 + gx2x5*p3 + gx3x5*p0 + gx4x5*p4 + gc[17]*s0 + gc[18]*s1 + gc[19]*p0 + gc[20]*s2
                 + Rgas*T*(rlog((xfe3tet/xal3tet)) + rlog((xfe3oct/xal3oct)))
                 + p4*(WV1*p2+WV2*p4)*((d->P*1000.0)-1.0) + p2*p4*WV1*((d->P*1000.0)-1.0);
    double dGdp3 = gx2x2*2.0*p3 + gx2x3*p0 + gx2x4*p4 + gx2x5*p2 + gc[5]*s0 + gc[6]*s1 + gc[7]*p0 + gc[8]*s2
                 + 0.5*Rgas*T*(rlog((xmg2tet/xfe2tet)) + rlog((xmg2oct/xfe2oct)));
    double dGdp4 = gx4x4*2.0*p4 + gx2x4*p3 + gx3x4*p0 + gx4x5*p2 + gc[13]*s0 + gc[14]*s1 + gc[15]*p0 + gc[16]*s2
                 + Rgas*T*(rlog((xfe2tet/xal3tet)) + rlog((xti4oct/xal3oct)))
                 + p2*(WV1*p2+WV2*p4)*((d->P*1000.0)-1.0) + p2*p4*WV2*((d->P*1000.0)-1.0);

    double dGex[5];
    dGex[0] = dGdp0 - CR_G;
    dGex[1] = -HC_G;
    dGex[2] = dGdp2 - MT_G;
    dGex[3] = dGdp3 - SP_G;
    dGex[4] = dGdp4 - UV_G;

    GH_mu_from_dGex(d, p, 1.0, 5, Gex, dGex);

    GH_sf_spn(p, d->sf);
    return GH_df_assemble(d, p, gb, mu_Gex, 5, grad);
}

typedef struct { GH_hsv d, w, dwhm, whm2, dwcrn, whm; double wcrn; } GH_rhm_spc;
static const GH_rhm_spc GH_rhm_SPC[3] = {
    { {17477.0, 0.0, 0.010758}, {3189.0, 0.0, 0.035089}, {-5626.63, 0.0, 0.013701}, {-833.14, 0.0, 0.0}, {0.0, 0.0, 0.013701}, {22535.6, 0.0, -0.11764}, 22535.6 },
    { {17477.0, 0.0, 0.010758}, {3189.0, 0.0, 0.035089}, {-5626.63, 0.0, 0.013701}, {-833.14, 0.0, 0.0}, {0.0, 0.0, 0.013701}, {22535.6, 0.0, -0.11764}, 22535.6 },
    { {17477.0, 0.0, 0.010758}, {3189.0, 0.0, 0.035089}, {-5626.63, 0.0, 0.013701}, {-833.14, 0.0, 0.0}, {0.0, 0.0, 0.013701}, {22535.6, 0.0, -0.11764}, 22535.6 }
};
static const GH_hsv GH_rhm_PAIR[3]  = { {2600.0, 0.0, 0.0}, {2200.0, 0.0, 0.0}, {2600.0, 0.0, 0.0} };
static const double GH_rhm_PAIRT[3] = { 88099.9, 30244.0, 2600.0 };
static const double GH_rhm_whmcrn = 69000.0;

typedef struct { double d[3], w[3], dwhm[3], dwcrn[3], whm[3], wcrn[3], a[3], wsum[3], wdif[3]; } GH_rhm_gp;

static void GH_rhm_gpar(double p, GH_rhm_gp *g){
    double dP = p - 1.0;
    for (int i = 0; i < 3; i++){
        const GH_rhm_spc *c = &GH_rhm_SPC[i];
        g->d[i]     = GH_hsv_G(c->d, 0.0, dP);
        g->w[i]     = GH_hsv_G(c->w, 0.0, dP);
        g->dwhm[i]  = GH_hsv_G(c->dwhm, 0.0, dP);
        g->dwcrn[i] = GH_hsv_G(c->dwcrn, 0.0, dP);
        g->whm[i]   = GH_hsv_G(c->whm, 0.0, dP);
        g->wcrn[i]  = c->wcrn;
        g->a[i]     = g->dwhm[i]/2.0 + GH_hsv_G(c->whm2, 0.0, dP)/4.0;
        g->wsum[i]  = GH_hsv_G(GH_rhm_PAIR[i], 0.0, dP) + GH_rhm_PAIRT[i];
        g->wdif[i]  = GH_rhm_PAIR[i].h - GH_rhm_PAIRT[i];
    }
}
static inline double GH_rhm_B(const GH_rhm_gp *g, int i, double sumr, double r3){
    return g->d[i] + g->a[i]*(1.0-sumr) - g->dwcrn[i]*r3/2.0;
}

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
static void GH_rhm_dgds(const double r[4], const double s[3], double t, const GH_rhm_gp *g,
                         double *dgds0, double *dgds1, double *dgds2){
    double xfe2a,xmg2a,xmn2a,xti4a,xal3a,xfe3a,xfe2b,xmg2b,xmn2b,xti4b,xal3b,xfe3b;
    GH_rhm_site_fracs(r,s,&xfe2a,&xmg2a,&xmn2a,&xti4a,&xal3a,&xfe3a,&xfe2b,&xmg2b,&xmn2b,&xti4b,&xal3b,&xfe3b);
    double R = GH_rhm_R;
    double r0=r[0],r1=r[1],r2=r[2],r3=r[3], s0=s[0],s1=s[1],s2=s[2];
    double sumr = r0+r1+r2+r3;

    *dgds0 = 0.5*R*t*( rlog((xfe2a)) - rlog((xti4a)) - rlog((xfe2b)) + rlog((xti4b)) )
           + g->wdif[0]*s1/2.0 + g->wdif[1]*s2/2.0
           - 2.0*GH_rhm_B(g,0,sumr,r3)*s0
           + 2.0*g->w[0]*s0*(r0*r0-2.0*s0*s0);

    *dgds1 = 0.5*R*t*( rlog((xmg2a)) - rlog((xti4a)) - rlog((xmg2b)) + rlog((xti4b)) )
           + g->wdif[0]*s0/2.0 + g->wdif[2]*s2/2.0
           - 2.0*GH_rhm_B(g,1,sumr,r3)*s1
           + 2.0*g->w[1]*s1*(r1*r1-2.0*s1*s1);

    *dgds2 = 0.5*R*t*( rlog((xmn2a)) - rlog((xti4a)) - rlog((xmn2b)) + rlog((xti4b)) )
           + g->wdif[1]*s0/2.0 + g->wdif[2]*s1/2.0
           - 2.0*GH_rhm_B(g,2,sumr,r3)*s2
           + 2.0*g->w[2]*s2*(r2*r2-2.0*s2*s2);
}
static void GH_rhm_d2gds2(const double r[4], const double s[3], double t, const GH_rhm_gp *g,
                           double *d00,double *d01,double *d02,double *d11,double *d12,double *d22){
    double xfe2a,xmg2a,xmn2a,xti4a,xal3a,xfe3a,xfe2b,xmg2b,xmn2b,xti4b,xal3b,xfe3b;
    GH_rhm_site_fracs(r,s,&xfe2a,&xmg2a,&xmn2a,&xti4a,&xal3a,&xfe3a,&xfe2b,&xmg2b,&xmn2b,&xti4b,&xal3b,&xfe3b);
    double R = GH_rhm_R;
    double r0=r[0],r1=r[1],r2=r[2],r3=r[3], s0=s[0],s1=s[1],s2=s[2];
    double sumr = r0+r1+r2+r3;

    *d00 = 0.25*R*t*( 1.0/xfe2a + 1.0/xti4a + 1.0/xfe2b + 1.0/xti4b )
         - 2.0*GH_rhm_B(g,0,sumr,r3) + 2.0*g->w[0]*(r0*r0-6.0*s0*s0);
    *d01 = 0.25*R*t*( 1.0/xti4a + 1.0/xti4b ) + g->wdif[0]/2.0;
    *d02 = 0.25*R*t*( 1.0/xti4a + 1.0/xti4b ) + g->wdif[1]/2.0;
    *d11 = 0.25*R*t*( 1.0/xmg2a + 1.0/xti4a + 1.0/xmg2b + 1.0/xti4b )
         - 2.0*GH_rhm_B(g,1,sumr,r3) + 2.0*g->w[1]*(r1*r1-6.0*s1*s1);
    *d12 = 0.25*R*t*( 1.0/xti4a + 1.0/xti4b ) + g->wdif[2]/2.0;
    *d22 = 0.25*R*t*( 1.0/xmn2a + 1.0/xti4a + 1.0/xmn2b + 1.0/xti4b )
         - 2.0*GH_rhm_B(g,2,sumr,r3) + 2.0*g->w[2]*(r2*r2-6.0*s2*s2);
}

static void GH_rhm_solve_s_from(const double r[4], double t, const GH_rhm_gp *g,
                                 double *s0o, double *s1o, double *s2o){
    const int maxIter = 100;
    double eps = 2.220446049250313e-16;
    double s[3] = { 0.9*r[0], 0.9*r[1], 0.9*r[2] };
    double sOld[3] = { 2.0, 2.0, 2.0 };
    int iter = 0;
    while ( (fabs(s[0]-sOld[0]) > 10.0*eps ||
             fabs(s[1]-sOld[1]) > 10.0*eps ||
             fabs(s[2]-sOld[2]) > 10.0*eps) && iter < maxIter ){
        double dgds0,dgds1,dgds2;
        GH_rhm_dgds(r, s, t, g, &dgds0, &dgds1, &dgds2);
        double d00,d01,d02,d11,d12,d22;
        GH_rhm_d2gds2(r, s, t, g, &d00,&d01,&d02,&d11,&d12,&d22);

        sOld[0]=s[0]; sOld[1]=s[1]; sOld[2]=s[2];

        double g0=-dgds0, g1=-dgds1, g2=-dgds2;
        double det  = d00*(d11*d22-d12*d12) - d01*(d01*d22-d12*d02) + d02*(d01*d12-d11*d02);
        double det0 =  g0*(d11*d22-d12*d12) - d01*( g1*d22-d12* g2) + d02*( g1*d12-d11* g2);
        double det1 = d00*( g1*d22- g2*d12) -  g0*(d01*d22-d12*d02) + d02*(d01* g2- g1*d02);
        double det2 = d00*(d11* g2-d12* g1) - d01*(d01* g2- g1*d02) +  g0*(d01*d12-d11*d02);
        double delta0 = det0/det, delta1 = det1/det, delta2 = det2/det;

        s[0] = sOld[0] + delta0;
        s[1] = sOld[1] + delta1;
        s[2] = sOld[2] + delta2;

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

    return -R*( xfe2a*rlog((xfe2a)) + xmg2a*rlog((xmg2a)) + xmn2a*rlog((xmn2a)) + xti4a*rlog((xti4a))
              + xfe2b*rlog((xfe2b)) + xmg2b*rlog((xmg2b)) + xmn2b*rlog((xmn2b)) + xti4b*rlog((xti4b))
              - 2.0*xfe2ID*rlog((xfe2ID)) - 2.0*xmg2ID*rlog((xmg2ID)) - 2.0*xmn2ID*rlog((xmn2ID)) - 2.0*xti4ID*rlog((xti4ID)) )
         - (1.0-GH_rhm_SROconst)*2.0*R*( xfe3ID*rlog((xfe3ID)) + xal3ID*rlog((xal3ID)) + xfe2ID*rlog((xfe2ID))
                                        + xmg2ID*rlog((xmg2ID)) + xmn2ID*rlog((xmn2ID)) + xti4ID*rlog((xti4ID)) )
         + GH_rhm_SROconst*2.0*R*rlog((2.0));
}
static double GH_rhm_H(const double r[4], const double s[3], const GH_rhm_gp *g){
    double r0=r[0],r1=r[1],r2=r[2],r3=r[3], s0=s[0],s1=s[1],s2=s[2];
    double sumr = r0+r1+r2+r3;
    return
          (g->whm[0]+g->dwhm[0]*(1.0-2.0*r0-r1-r2-r3))*r0*(1.0-sumr)
        + (g->whm[1]+g->dwhm[1]*(1.0-r0-2.0*r1-r2-r3))*r1*(1.0-sumr)
        + (g->whm[2]+g->dwhm[2]*(1.0-r0-r1-2.0*r2-r3))*r2*(1.0-sumr)
        +  GH_rhm_whmcrn*r3*(1.0-sumr)
        + g->wsum[0]*r0*r1/2.0 + g->wdif[0]*s0*s1/2.0
        + g->wsum[1]*r0*r2/2.0 + g->wdif[1]*s0*s2/2.0
        + g->wsum[2]*r1*r2/2.0 + g->wdif[2]*s1*s2/2.0
        + (g->wcrn[0]+g->dwcrn[0]*(r0-r3))*r0*r3
        + (g->wcrn[1]+g->dwcrn[1]*(r1-r3))*r1*r3
        + (g->wcrn[2]+g->dwcrn[2]*(r2-r3))*r2*r3
        + GH_rhm_B(g,0,sumr,r3)*(r0*r0-s0*s0)
        + GH_rhm_B(g,1,sumr,r3)*(r1*r1-s1*s1)
        + GH_rhm_B(g,2,sumr,r3)*(r2*r2-s2*s2)
        + g->w[0]*s0*s0*(r0*r0-s0*s0)
        + g->w[1]*s1*s1*(r1*r1-s1*s1)
        + g->w[2]*s2*s2*(r2*r2-s2*s2);
}
static void GH_rhm_dGdr(const double r[4], const double s[3], double t, const GH_rhm_gp *g,
                         double *dgdr0,double *dgdr1,double *dgdr2,double *dgdr3){
    double xfe2a,xmg2a,xmn2a,xti4a,xal3a,xfe3a,xfe2b,xmg2b,xmn2b,xti4b,xal3b,xfe3b;
    GH_rhm_site_fracs(r,s,&xfe2a,&xmg2a,&xmn2a,&xti4a,&xal3a,&xfe3a,&xfe2b,&xmg2b,&xmn2b,&xti4b,&xal3b,&xfe3b);
    double xfe2ID,xmg2ID,xmn2ID,xti4ID,xal3ID,xfe3ID;
    GH_rhm_ID_fracs(r,&xfe2ID,&xmg2ID,&xmn2ID,&xti4ID,&xal3ID,&xfe3ID);
    double R = GH_rhm_R;
    double r0=r[0],r1=r[1],r2=r[2],r3=r[3], s0=s[0],s1=s[1],s2=s[2];
    double sumr = r0+r1+r2+r3;
    double q0 = r0*r0-s0*s0, q1 = r1*r1-s1*s1, q2 = r2*r2-s2*s2;
    double u0 = g->whm[0]+g->dwhm[0]*(1.0-2.0*r0-r1-r2-r3);
    double u1 = g->whm[1]+g->dwhm[1]*(1.0-r0-2.0*r1-r2-r3);
    double u2 = g->whm[2]+g->dwhm[2]*(1.0-r0-r1-2.0*r2-r3);
    double hmsum = g->dwhm[0]*r0*(1.0-sumr) + g->dwhm[1]*r1*(1.0-sumr) + g->dwhm[2]*r2*(1.0-sumr);

    *dgdr0 = R*t*( 0.5*rlog((xfe2a)) + 0.5*rlog((xti4a)) + 0.5*rlog((xfe2b)) + 0.5*rlog((xti4b)) - rlog((xfe2ID)) - rlog((xti4ID)) )
           + (1.0-GH_rhm_SROconst)*2.0*R*t*(- rlog((xfe3ID)) + 0.5*rlog((xfe2ID)) + 0.5*rlog((xti4ID)) )
           + u0*(1.0-2.0*r0-r1-r2-r3) - u1*r1 - u2*r2
           - g->dwhm[0]*r0*(1.0-sumr) - hmsum
           - GH_rhm_whmcrn*r3
           + g->wsum[0]*r1/2.0 + g->wsum[1]*r2/2.0
           + (g->wcrn[0]+g->dwcrn[0]*(r0-r3))*r3 + g->dwcrn[0]*r0*r3
           + 2.0*GH_rhm_B(g,0,sumr,r3)*r0
           - g->a[0]*q0 - g->a[1]*q1 - g->a[2]*q2
           + 2.0*g->w[0]*s0*s0*r0;

    *dgdr1 = R*t*( 0.5*rlog((xmg2a)) + 0.5*rlog((xti4a)) + 0.5*rlog((xmg2b)) + 0.5*rlog((xti4b)) - rlog((xmg2ID)) - rlog((xti4ID)) )
           + (1.0-GH_rhm_SROconst)*2.0*R*t*(- rlog((xfe3ID)) + 0.5*rlog((xmg2ID)) + 0.5*rlog((xti4ID)) )
           - u0*r0 + u1*(1.0-r0-2.0*r1-r2-r3) - u2*r2
           - g->dwhm[1]*r1*(1.0-sumr) - hmsum
           - GH_rhm_whmcrn*r3
           + g->wsum[0]*r0/2.0 + g->wsum[2]*r2/2.0
           + (g->wcrn[1]+g->dwcrn[1]*(r1-r3))*r3 + g->dwcrn[1]*r1*r3
           - g->a[0]*q0
           + 2.0*GH_rhm_B(g,1,sumr,r3)*r1
           - g->a[1]*q1 - g->a[2]*q2
           + 2.0*g->w[1]*s1*s1*r1;

    *dgdr2 = R*t*( 0.5*rlog((xmn2a)) + 0.5*rlog((xti4a)) + 0.5*rlog((xmn2b)) + 0.5*rlog((xti4b)) - rlog((xmn2ID)) - rlog((xti4ID)) )
           + (1.0-GH_rhm_SROconst)*2.0*R*t*(- rlog((xfe3ID)) + 0.5*rlog((xmn2ID)) + 0.5*rlog((xti4ID)) )
           - u0*r0 - u1*r1 + u2*(1.0-r0-r1-2.0*r2-r3)
           - g->dwhm[2]*r2*(1.0-sumr) - hmsum
           - GH_rhm_whmcrn*r3
           + g->wsum[1]*r0/2.0 + g->wsum[2]*r1/2.0
           + (g->wcrn[2]+g->dwcrn[2]*(r2-r3))*r3 + g->dwcrn[2]*r2*r3
           - g->a[0]*q0 - g->a[1]*q1
           + 2.0*GH_rhm_B(g,2,sumr,r3)*r2
           - g->a[2]*q2
           + 2.0*g->w[2]*s2*s2*r2;

    *dgdr3 = (1.0-GH_rhm_SROconst)*2.0*R*t*(- rlog((xfe3ID)) + rlog((xal3ID)) )
           - u0*r0 - u1*r1 - u2*r2
           - hmsum
           + GH_rhm_whmcrn*(1.0-r0-r1-r2-2.0*r3)
           + (g->wcrn[0]+g->dwcrn[0]*(r0-r3))*r0 - g->dwcrn[0]*r0*r3
           + (g->wcrn[1]+g->dwcrn[1]*(r1-r3))*r1 - g->dwcrn[1]*r1*r3
           + (g->wcrn[2]+g->dwcrn[2]*(r2-r3))*r2 - g->dwcrn[2]*r2*r3
           - (g->a[0] + g->dwcrn[0]/2.0)*q0
           - (g->a[1] + g->dwcrn[1]/2.0)*q1
           - (g->a[2] + g->dwcrn[2]/2.0)*q2;
}

static double GH_rhm_mix_endmember_order_G(double T, double dG, double wG){
    double R = 8.3143;
    double eps = 2.220446049250313e-16;
    double s = 0.98;
    for (int iter = 0; iter < 1000; iter++){
        double dgds   = R*T*(rlog((1.0+s)) - rlog((1.0-s)))
                      - 2.0*dG*s + wG*(2.0*s - 4.0*s*s*s);
        double d2gds2 = R*T*(1.0/(1.0+s) + 1.0/(1.0-s))
                      - 2.0*dG + wG*(2.0 - 12.0*s*s);
        double sNew   = s - dgds/d2gds2;
        sNew = (sNew > 1.0 - 10.0*eps) ? 1.0 - 10.0*eps : sNew;
        sNew = (sNew < 0.0) ? 0.0 : sNew;
        int converged = (fabs(sNew - s) <= 10.0*eps);
        s = sNew;
        if (converged) break;
    }
    double S = -R*((1.0+s)*rlog((1.0+s)) + (1.0-s)*rlog((1.0-s)) - 2.0*rlog((2.0)));
    double H = dG*(1.0 - s*s) + wG*s*s*(1.0 - s*s);
    return H - T*S;
}
static double GH_rhm_mix_endmember_SRO_G(double T){
    const double SROconst = 0.0730205;
    double R = 8.3143;
    return -T*(SROconst*2.0*R*rlog((2.0)));
}

enum { IL, GK, PY };
enum { FE2, MG2, MN2, TI4, FE3 };

typedef struct {
    double h[3], whm[3], dwhm[3];
    double dis[3][3], ord[3][3], dw[3][3];
} GH_r91_par;

static const GH_r91_par GH_R91 = {
    .h    = {17000.0, 17000.0, 17000.0},
    .whm  = {5000.0, 5000.0, 5000.0},
    .dwhm = {-19000.0, -19000.0, -19000.0},
    .dis  = { [IL][GK] = 25000.0, [IL][PY] = 2200.0, [GK][PY] = 25000.0 },
    .ord  = { [IL][GK] = 5000.0,  [IL][PY] = 2200.0, [GK][PY] = 5000.0 },
    .dw   = { [IL][GK] = 12000.0, [GK][IL] = 8000.0, [IL][PY] = 17000.0, [PY][IL] = 17000.0,
              [GK][PY] = 8000.0,  [PY][GK] = 12000.0 },
};

typedef struct {
    double h[3], whm[3], dis[3][3];
    double ahm[3], bex[3][3], cord[3][3];
} GH_r91_coef;

typedef struct { double a[5], b[5]; } GH_r91_sf;

static void GH_r91_coef_make(const GH_r91_par *P, GH_r91_coef *c){
    for (int i = 0; i < 3; i++){
        c->h[i]   = P->h[i];
        c->whm[i] = P->whm[i];
        c->ahm[i] = P->dwhm[i] + P->h[i];
        for (int j = 0; j < 3; j++){
            c->dis[i][j]  = P->dis[i][j];
            c->bex[i][j]  = P->dw[i][j] - P->h[i];
            c->cord[i][j] = P->ord[i][j] - P->dis[i][j] - P->dw[i][j] - P->dw[j][i];
        }
    }
}

static void GH_r91_site_fracs(const double r[3], const double s[3], GH_r91_sf *f){
    f->a[FE2] = (r[0] + s[0])/2.0;
    f->a[MG2] = (r[1] + s[1])/2.0;
    f->a[MN2] = (r[2] + s[2])/2.0;
    f->a[TI4] = (r[0] - s[0] + r[1] - s[1] + r[2] - s[2])/2.0;
    f->a[FE3] = 1.0 - r[0] - r[1] - r[2];
    f->b[FE2] = (r[0] - s[0])/2.0;
    f->b[MG2] = (r[1] - s[1])/2.0;
    f->b[MN2] = (r[2] - s[2])/2.0;
    f->b[TI4] = (r[0] + s[0] + r[1] + s[1] + r[2] + s[2])/2.0;
    f->b[FE3] = 1.0 - r[0] - r[1] - r[2];
    for (int i = 0; i < 5; i++){
        if (f->a[i] <= 0.0) f->a[i] = DBL_EPSILON;
        if (f->b[i] <= 0.0) f->b[i] = DBL_EPSILON;
    }
}

static double GH_r91_entropy(const GH_r91_sf *f){
    double acc = 0.0;
    for (int i = 0; i < 5; i++){ acc += f->a[i]*log(f->a[i]); }
    for (int i = 0; i < 5; i++){ acc += f->b[i]*log(f->b[i]); }
    return -GH_R_J*acc;
}

static double GH_r91_enthalpy(const GH_r91_coef *c, const double r[3], const double s[3]){
    return c->h[IL]*(r[0] - s[0]*s[0]) + c->h[GK]*(r[1] - s[1]*s[1]) + c->h[PY]*(r[2] - s[2]*s[2])
            + c->whm[IL]*r[0]*(1.0 - r[0] - r[1] - r[2]) + c->dis[IL][GK]*r[0]*r[1] + c->dis[IL][PY]*r[0]*r[2]
            + c->whm[GK]*r[1]*(1.0 - r[0] - r[1] - r[2]) + c->dis[GK][PY]*r[1]*r[2]
            + c->whm[PY]*r[2]*(1.0 - r[0] - r[1] - r[2]) - 2.0*s[0]*s[0]*(1.0 - r[0] - r[1] - r[2])*c->ahm[IL]
            - 2.0*s[1]*s[1]*(1.0 - r[0] - r[1] - r[2])*c->ahm[GK]
            - 2.0*s[2]*s[2]*(1.0 - r[0] - r[1] - r[2])*c->ahm[PY] + 2.0*s[0]*s[0]*r[1]*c->bex[IL][GK]
            + 2.0*s[0]*s[0]*r[2]*c->bex[IL][PY] + 2.0*s[1]*s[1]*r[0]*c->bex[GK][IL]
            + 2.0*s[1]*s[1]*r[2]*c->bex[GK][PY] + 2.0*s[2]*s[2]*r[0]*c->bex[PY][IL]
            + 2.0*s[2]*s[2]*r[1]*c->bex[PY][GK] + s[0]*s[1]*c->cord[IL][GK] + s[0]*s[2]*c->cord[IL][PY]
            + s[1]*s[2]*c->cord[GK][PY];
}

static double GH_r91_gibbs(const GH_r91_coef *c, const double r[3], const double s[3], const GH_r91_sf *f, double t){
    return GH_r91_enthalpy(c, r, s) - t*GH_r91_entropy(f);
}

static void GH_r91_dg_dr(const GH_r91_coef *c, const double r[3], const double s[3], const GH_r91_sf *f, double t, double dr[3]){
    const double *xa = f->a, *xb = f->b;
    dr[0] = GH_R_J*t*(0.5*log(xa[FE2]) + 0.5*log(xa[TI4]) - log(xa[FE3]) + 0.5*log(xb[FE2]) + 0.5*log(xb[TI4]) - log(xb[FE3]))
            + c->h[IL] + c->whm[IL]*(1.0 - 2.0*r[0] - r[1] - r[2]) + c->dis[IL][GK]*r[1] + c->dis[IL][PY]*r[2]
            - c->whm[GK]*r[1] - c->whm[PY]*r[2] + 2.0*s[0]*s[0]*c->ahm[IL] + 2.0*s[1]*s[1]*c->ahm[GK]
            + 2.0*s[2]*s[2]*c->ahm[PY] + 2.0*s[1]*s[1]*c->bex[GK][IL] + 2.0*s[2]*s[2]*c->bex[PY][IL];
    dr[1] = GH_R_J*t*(0.5*log(xa[MG2]) + 0.5*log(xa[TI4]) - log(xa[FE3]) + 0.5*log(xb[MG2]) + 0.5*log(xb[TI4]) - log(xb[FE3]))
            + c->h[GK] - c->whm[IL]*r[0] + c->dis[IL][GK]*r[0] + c->whm[GK]*(1.0 - r[0] - 2.0*r[1] - r[2])
            + c->dis[GK][PY]*r[2] - c->whm[PY]*r[2] + 2.0*s[0]*s[0]*c->ahm[IL] + 2.0*s[1]*s[1]*c->ahm[GK]
            + 2.0*s[2]*s[2]*c->ahm[PY] + 2.0*s[0]*s[0]*c->bex[IL][GK] + 2.0*s[2]*s[2]*c->bex[PY][GK];
    dr[2] = GH_R_J*t*(0.5*log(xa[MN2]) + 0.5*log(xa[TI4]) - log(xa[FE3]) + 0.5*log(xb[MN2]) + 0.5*log(xb[TI4]) - log(xb[FE3]))
            + c->h[PY] - c->whm[IL]*r[0] + c->dis[IL][PY]*r[0] - c->whm[GK]*r[1] + c->dis[GK][PY]*r[1]
            + c->whm[PY]*(1.0 - r[0] - r[1] - 2.0*r[2]) + 2.0*s[0]*s[0]*c->ahm[IL] + 2.0*s[1]*s[1]*c->ahm[GK]
            + 2.0*s[2]*s[2]*c->ahm[PY] + 2.0*s[0]*s[0]*c->bex[IL][PY] + 2.0*s[1]*s[1]*c->bex[GK][PY];
}

static void GH_r91_dg_ds(const GH_r91_coef *c, const double r[3], const double s[3], const GH_r91_sf *f, double t, double dg[3]){
    const double *xa = f->a, *xb = f->b;
    dg[0] = 0.5*GH_R_J*t*(log(xa[FE2]) - log(xa[TI4]) - log(xb[FE2]) + log(xb[TI4])) - 2.0*c->h[IL]*s[0]
            - 4.0*s[0]*(1.0 - r[0] - r[1] - r[2])*c->ahm[IL] + 4.0*s[0]*r[1]*c->bex[IL][GK]
            + 4.0*s[0]*r[2]*c->bex[IL][PY] + s[1]*c->cord[IL][GK] + s[2]*c->cord[IL][PY];
    dg[1] = 0.5*GH_R_J*t*(log(xa[MG2]) - log(xa[TI4]) - log(xb[MG2]) + log(xb[TI4])) - 2.0*c->h[GK]*s[1]
            - 4.0*s[1]*(1.0 - r[0] - r[1] - r[2])*c->ahm[GK] + 4.0*s[1]*r[0]*c->bex[GK][IL]
            + 4.0*s[1]*r[2]*c->bex[GK][PY] + s[0]*c->cord[IL][GK] + s[2]*c->cord[GK][PY];
    dg[2] = 0.5*GH_R_J*t*(log(xa[MN2]) - log(xa[TI4]) - log(xb[MN2]) + log(xb[TI4])) - 2.0*c->h[PY]*s[2]
            - 4.0*s[2]*(1.0 - r[0] - r[1] - r[2])*c->ahm[PY] + 4.0*s[2]*r[0]*c->bex[PY][IL]
            + 4.0*s[2]*r[1]*c->bex[PY][GK] + s[0]*c->cord[IL][PY] + s[1]*c->cord[GK][PY];
}

static void GH_r91_d2g_ds2(const GH_r91_coef *c, const double r[3], const GH_r91_sf *f, double t, double H[3][3]){
    const double *xa = f->a, *xb = f->b;
    H[0][0] = 0.25*GH_R_J*t*(1.0/xa[FE2] + 1.0/xa[TI4] + 1.0/xb[FE2] + 1.0/xb[TI4]) - 2.0*c->h[IL]
            - 4.0*(1.0 - r[0] - r[1] - r[2])*c->ahm[IL] + 4.0*r[1]*c->bex[IL][GK] + 4.0*r[2]*c->bex[IL][PY];
    H[0][1] = 0.25*GH_R_J*t*(1.0/xa[TI4] + 1.0/xb[TI4]) + c->cord[IL][GK];
    H[0][2] = 0.25*GH_R_J*t*(1.0/xa[TI4] + 1.0/xb[TI4]) + c->cord[IL][PY];
    H[1][1] = 0.25*GH_R_J*t*(1.0/xa[MG2] + 1.0/xa[TI4] + 1.0/xb[MG2] + 1.0/xb[TI4]) - 2.0*c->h[GK]
            - 4.0*(1.0 - r[0] - r[1] - r[2])*c->ahm[GK] + 4.0*r[0]*c->bex[GK][IL] + 4.0*r[2]*c->bex[GK][PY];
    H[1][2] = 0.25*GH_R_J*t*(1.0/xa[TI4] + 1.0/xb[TI4]) + c->cord[GK][PY];
    H[2][2] = 0.25*GH_R_J*t*(1.0/xa[MN2] + 1.0/xa[TI4] + 1.0/xb[MN2] + 1.0/xb[TI4]) - 2.0*c->h[PY]
            - 4.0*(1.0 - r[0] - r[1] - r[2])*c->ahm[PY] + 4.0*r[0]*c->bex[PY][IL] + 4.0*r[1]*c->bex[PY][GK];
    H[1][0] = H[0][1]; H[2][0] = H[0][2]; H[2][1] = H[1][2];
}

static double GH_r91_pure_S(double s){
    return - GH_R_J*( (1.0 + s)*log(1.0 + s) + (1.0 - s)*log(1.0 - s) - 2.0*log(2.0) );
}

static double GH_r91_pure_G(double h, double s, double t){
    return h*(1.0 - s*s) - t*(GH_r91_pure_S(s));
}

static double GH_r91_pure_dg_ds(double h, double s, double t){
    return GH_R_J*t*( log(1.0 + s) - log(1.0 - s) ) - 2.0*h*s;
}

static double GH_r91_pure_d2g_ds2(double h, double s, double t){
    return GH_R_J*t*( 1.0/(1.0 + s) + 1.0/(1.0 - s) ) - 2.0*h;
}

static void GH_r91_pure_ends(const GH_r91_coef *c, double t, double ends[4]){
    double sOld[3], sNew[3];
    int iter = 0;
    for (int i = 0; i < 3; i++){ sOld[i] = 2.0; sNew[i] = 0.9; }
    while ((fabs(sNew[0]-sOld[0]) > 10.0*DBL_EPSILON || fabs(sNew[1]-sOld[1]) > 10.0*DBL_EPSILON ||
            fabs(sNew[2]-sOld[2]) > 10.0*DBL_EPSILON) && iter < 100){
        double s[3] = {sNew[0], sNew[1], sNew[2]};
        for (int i = 0; i < 3; i++){ sOld[i] = s[i]; }
        for (int i = 0; i < 3; i++){
            s[i] += -GH_r91_pure_dg_ds(c->h[i], sOld[i], t)/GH_r91_pure_d2g_ds2(c->h[i], sOld[i], t);
            if (s[i] > 1.0 - DBL_EPSILON) s[i] = 1.0 - DBL_EPSILON;
            if (s[i] < DBL_EPSILON) s[i] = DBL_EPSILON;
        }
        for (int i = 0; i < 3; i++){ sNew[i] = s[i]; }
        iter++;
    }
    ends[0] = GH_r91_pure_G(c->h[GK], sOld[GK], t);
    ends[1] = 0.0;
    ends[2] = GH_r91_pure_G(c->h[IL], sOld[IL], t);
    ends[3] = GH_r91_pure_G(c->h[PY], sOld[PY], t);
}

static void GH_r91_order(const GH_r91_coef *c, double t, const double r[3], double s_out[3], GH_r91_sf *sf){
    GH_r91_sf f = {{0.0}, {0.0}};
    double sOld[3], sNew[3];
    int iter = 0;
    for (int i = 0; i < 3; i++){ sOld[i] = 2.0; sNew[i] = 0.9*r[i]; }
    while ((fabs(sNew[0]-sOld[0]) > 10.0*DBL_EPSILON || fabs(sNew[1]-sOld[1]) > 10.0*DBL_EPSILON ||
            fabs(sNew[2]-sOld[2]) > 10.0*DBL_EPSILON) && iter < 100){
        double s[3] = {sNew[0], sNew[1], sNew[2]};
        double dgds[3], H[3][3], inv[3][3];

        GH_r91_site_fracs(r, s, &f);
        GH_r91_dg_ds(c, r, s, &f, t, dgds);
        GH_r91_d2g_ds2(c, r, &f, t, H);
        GH_invert3(H, inv);

        for (int i = 0; i < 3; i++){ sOld[i] = s[i]; }
        for (int i = 0; i < 3; i++){
            for (int j = 0; j < 3; j++){ s[i] += -inv[i][j]*dgds[j]; }
            if (s[i] > r[i] - DBL_EPSILON) s[i] = r[i] - DBL_EPSILON;
            if (s[i] < 0.0) s[i] = 0.0;
        }
        for (int i = 0; i < 3; i++){ sNew[i] = s[i]; }
        iter++;
    }
    for (int i = 0; i < 3; i++){ s_out[i] = sOld[i]; }
    *sf = f;
}

static double GH_r91_gmix(double t, const double r[3], double dx[3]){
    GH_r91_coef c;
    GH_r91_sf f;
    double s[3], ends[4];
    GH_r91_coef_make(&GH_R91, &c);
    GH_r91_order(&c, t, r, s, &f);
    GH_r91_pure_ends(&c, t, ends);
    double g = GH_r91_gibbs(&c, r, s, &f, t);
    g -= (r[1]*ends[0] + (1.0 - r[0] - r[1] - r[2])*ends[1] + r[0]*ends[2] + r[2]*ends[3]);
    if (dx){
        GH_r91_dg_dr(&c, r, s, &f, t, dx);
        dx[0] -= (ends[2] - ends[1]);
        dx[1] -= (ends[0] - ends[1]);
        dx[2] -= (ends[3] - ends[1]);
    }
    return g;
}

static double obj_gh_rhm91(unsigned n, const double *x, double *grad, void *SS_ref_db){
    SS_ref *d = (SS_ref *) SS_ref_db;

    double T    = d->T;
    double *p   = d->p;
    double *gb  = d->gb_lvl;
    double *mu_Gex = d->mu_Gex;

    for (int i = 0; i < 4; i++){ p[i] = x[i]; }
    double sum = p[0] + p[1] + p[2] + p[3];
    double r[3] = { p[2]/sum, p[0]/sum, p[3]/sum };

    double dx[3];
    double Gex = GH_r91_gmix(T, r, dx);

    double dGex[4];
    dGex[0] = dx[1];
    dGex[1] = 0.0;
    dGex[2] = dx[0];
    dGex[3] = dx[2];

    GH_mu_from_dGex(d, p, sum, 4, Gex, dGex);

    return GH_df_assemble(d, p, gb, mu_Gex, 4, grad);
}

double obj_gh_rhm(unsigned n, const double *x, double *grad, void *SS_ref_db){
    SS_ref *d = (SS_ref *) SS_ref_db;
    if (d->EM_database == 2){ return obj_gh_rhm91(n, x, grad, SS_ref_db); }

    double T    = d->T;
    double Pbar = d->P*1000.0;
    double *p   = d->p;
    double *gb  = d->gb_lvl;
    double *mu_Gex = d->mu_Gex;

    for (int i = 0; i < 5; i++){ p[i] = x[i]; }

    double r[4] = { p[2], p[0], p[3], p[4] };

    GH_rhm_gp gp;
    GH_rhm_gpar(Pbar, &gp);

    double s0, s1, s2;
    GH_rhm_solve_s_from(r, T, &gp, &s0, &s1, &s2);
    double s[3] = { s0, s1, s2 };

    double Graw = GH_rhm_H(r, s, &gp) - T*GH_rhm_S(r, s, T);

    double dgdr0, dgdr1, dgdr2, dgdr3;
    GH_rhm_dGdr(r, s, T, &gp, &dgdr0, &dgdr1, &dgdr2, &dgdr3);

    double corr_ord = GH_rhm_mix_endmember_order_G(T, gp.d[0], gp.w[0]);
    double corr_sro = GH_rhm_mix_endmember_SRO_G(T);
    double ENDMEMBERS = (p[0]+p[2]+p[3])*corr_ord + (p[1]+p[4])*corr_sro;
    double Gex = Graw - ENDMEMBERS;

    double dGex[5];
    dGex[0] = dgdr1 - corr_ord;
    dGex[1] = -corr_sro;
    dGex[2] = dgdr0 - corr_ord;
    dGex[3] = dgdr2 - corr_ord;
    dGex[4] = dgdr3 - corr_sro;

    GH_mu_from_dGex(d, p, 1.0, 5, Gex, dGex);

    return GH_df_assemble(d, p, gb, mu_Gex, 5, grad);
}

static const GH_hsv GH_nph_EX  = { -31744.63, -20.92, -1.05 };
static const GH_hsv GH_nph_X   = { -13893.18,  12.55,  0.0  };
static const GH_hsv GH_nph_G23 = {  73520.0,    0.0,   0.0  };
static const GH_hsv GH_nph_G24 = {  42560.0,    0.0,   0.0  };
static const GH_hsv GH_nph_WNAKLS = {  6861.76, 0.0, 0.33 };
static const GH_hsv GH_nph_WNAKSS = { 51002.96, 0.0, 0.54 };
static const double GH_nph_WVN=14644.0, GH_nph_WVK=8368.0;
static const double GH_nph_S4=-15.8765;
static const double GH_nph_R = 8.3143;
static const double GH_nph_PENALTY = 1.4901161193847656e-08;

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
static void GH_nph_gcoef(double t, double p, double c[10]){
    double dP = p - 1.0;
    double GEX = GH_hsv_G(GH_nph_EX, t, dP), GX = GH_hsv_G(GH_nph_X, t, dP);
    double WNAKLS = GH_hsv_G(GH_nph_WNAKLS, t, dP), WNAKSS = GH_hsv_G(GH_nph_WNAKSS, t, dP);
    double G23 = GH_hsv_G(GH_nph_G23, t, dP), G24 = GH_hsv_G(GH_nph_G24, t, dP);
    c[0] = 0.25*(2.0*GEX+GX+3.0*WNAKLS-WNAKSS);
    c[1] = GX+WNAKLS+WNAKSS;
    c[2] = (1.0/16.0)*(3.0*GX-9.0*WNAKLS-WNAKSS);
    c[3] = -0.5*(GX+3.0*WNAKLS-WNAKSS);
    c[4] = GH_nph_WVN;
    c[5] = 0.5*(2.0*G23+GEX-GX-2.0*WNAKLS-2.0*GH_nph_WVN+2.0*GH_nph_WVK);
    c[6] = 0.125*(GX-GEX-2.0*G23-6.0*WNAKLS-6.0*GH_nph_WVN+6.0*GH_nph_WVK);
    c[7] = (6.0*G23+3.0*GEX-15.0*GX-18.0*WNAKLS-4.0*WNAKSS+18.0*GH_nph_WVN-18.0*GH_nph_WVK)/18.0;
    c[8] = -GH_nph_WVN;
    c[9] = (18.0*GH_nph_WVN-18.0*GH_nph_WVK-18.0*WNAKLS+4.0*WNAKSS-9.0*GEX-3.0*GX+6.0*G23-36.0*G24)/24.0;
}
static double GH_nph_dgds0(const double c[10], double r0, double r1, double r2, double s0, double t){
    double xkls,xvcls,xnals,xkss,xcass,xnass;
    GH_nph_site_fracs(r0,r1,r2,s0,&xkls,&xvcls,&xnals,&xkss,&xcass,&xnass);
    double R = GH_nph_R;
    return R*t*(0.75*rlog((xkls)) - 0.75*rlog((xnals)) - 0.75*rlog((xkss)) + 0.75*rlog((xnass)))
         + c[0] + 2.0*c[2]*s0 + c[3]*r0 + c[6]*r1 + c[9]*r2;
}
static double GH_nph_d2gds0s0(const double c[10], double r0, double r1, double r2, double s0, double t){
    double xkls,xvcls,xnals,xkss,xcass,xnass;
    GH_nph_site_fracs(r0,r1,r2,s0,&xkls,&xvcls,&xnals,&xkss,&xcass,&xnass);
    double R = GH_nph_R;
    return R*t*(0.75*0.75/xkls + 0.75*0.75/xnals + 3.0*0.25*0.25/xkss + 3.0*0.25*0.25/xnass)
         + 2.0*c[2];
}

static double GH_nph_solve_s0(const double c[10], double r0, double r1, double r2, double t){
    const int MAX_ITER = 200;
    double eps = 2.220446049250313e-16;

    double sNew = -4.0*r0*(r1+2.0*r2/3.0)/(4.0-r1-2.0*r2);
    double xkls,xvcls,xnals,xkss,xcass,xnass;
    GH_nph_site_fracs(r0,r1,r2,sNew,&xkls,&xvcls,&xnals,&xkss,&xcass,&xnass);
    sNew = xkls - xkss;

    double sOld = 2.0;
    int iter = 0;
    while ((fabs(sNew-sOld) > 10.0*eps) && (iter < MAX_ITER)){
        double s0 = sNew;
        double dgds  = GH_nph_dgds0(c,r0,r1,r2,s0,t);
        double d2gds2 = GH_nph_d2gds0s0(c,r0,r1,r2,s0,t);
        sOld = s0;

        double sCorr = -dgds/d2gds2;
        double lambda = 1.0;
        double sTrial = sOld + lambda*sCorr;
        double xkls_,xvcls_,xnals_,xkss_,xcass_,xnass_;
        GH_nph_site_fracs(r0,r1,r2,sTrial,&xkls_,&xvcls_,&xnals_,&xkss_,&xcass_,&xnass_);

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
static double GH_nph_G(const double c[10], double r0, double r1, double r2, double s0, double t){
    double xkls,xvcls,xnals,xkss,xcass,xnass;
    GH_nph_site_fracs(r0,r1,r2,s0,&xkls,&xvcls,&xnals,&xkss,&xcass,&xnass);
    double R = GH_nph_R;
    double Gsite = R*t*(xkls*rlog((xkls)) + xvcls*rlog((xvcls)) + xnals*rlog((xnals))
                      - (1.0-3.0*xcass)*rlog((1.0-3.0*xcass))
                      + 3.0*xkss*rlog((xkss)) + 3.0*xcass*rlog((xcass)) + 3.0*xnass*rlog((xnass)));
    return c[0]*s0 + c[1]*r0*(1.0-r0) + c[2]*s0*s0 + c[3]*r0*s0 + c[4]*r1*(1.0-r1)
         + c[5]*r0*r1 + c[6]*r1*s0 + c[7]*r0*r2 + c[8]*r1*r2 + c[9]*r2*s0
         - t*GH_nph_S4*r2 + Gsite;
}
static void GH_nph_dGdr(const double c[10], double r0, double r1, double r2, double s0, double t,
                         double *dgdr0, double *dgdr1, double *dgdr2){
    double xkls,xvcls,xnals,xkss,xcass,xnass;
    GH_nph_site_fracs(r0,r1,r2,s0,&xkls,&xvcls,&xnals,&xkss,&xcass,&xnass);
    double R = GH_nph_R;
    double G4 = -t*GH_nph_S4;

    *dgdr0 = R*t*(rlog((xkls)) - rlog((xnals)) + 3.0*rlog((xkss)) - 3.0*rlog((xnass)))
           + c[1]*(1.0-2.0*r0) + c[3]*s0 + c[5]*r1 + c[7]*r2;

    *dgdr1 = R*t*(rlog((xvcls)) - rlog((xnals))) + c[4]*(1.0-2.0*r1)
           + c[5]*r0 + c[6]*s0 + c[8]*r2
           - GH_nph_PENALTY/(r1*r1);

    *dgdr2 = G4 + R*t*(rlog((xvcls)) - rlog((xnals)) + rlog((1.0-3.0*xcass)) + rlog((xcass)) - rlog((xnass)) + 1.0)
           + c[7]*r0 + c[8]*r1 + c[9]*s0;
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

    double c[10];
    GH_nph_gcoef(T, Pbar, c);
    double s0 = GH_nph_solve_s0(c, r0, r1, r2, T);

    double Gex = GH_nph_G(c, r0, r1, r2, s0, T) + GH_nph_PENALTY/r1;

    double dgdr0, dgdr1, dgdr2;
    GH_nph_dGdr(c, r0, r1, r2, s0, T, &dgdr0, &dgdr1, &dgdr2);

    double dGex[4];
    dGex[0] = 0.0;
    dGex[1] = dgdr0;
    dGex[2] = dgdr1;
    dGex[3] = dgdr2;

    GH_mu_from_dGex(d, p, 1.0, 4, Gex, dGex);

    GH_sf_nph(p, d->sf);
    return GH_df_assemble(d, p, gb, mu_Gex, 4, grad);
}

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

    const GH_hsv D[4] = {
        { 32145.67,  10.46, 0.0  },
        { -5648.4,    0.0, -0.04 },
        { 14644.0,  -18.7,  0.0  },
        { 35184.0,   -7.17, 0.0  }
    };
    const double WKALS=29288.0, WKVKALS=30334.0, WNVKALS=14646.0, WKCKALS=14644.0;
    const double S4=-15.8765;
    const double R = 8.3143;
    const double PENALTY = 1.4901161193847656e-08;

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

    double g1 = GH_hsv_G(D[0], T, Pbar-1.0), g2 = GH_hsv_G(D[1], T, Pbar-1.0);
    double g3 = GH_hsv_G(D[2], T, Pbar-1.0), g4 = GH_hsv_G(D[3], T, Pbar-1.0);
    double G4 = -T*S4;

    double Gex = g1*(1.0-sumr) + g2*r0 + g3*r1 + (g4+G4)*r2
               + 4.0*R*T*(xk*rlog((xk)) + (xvc-xca)*rlog((xvc-xca)) + xna*rlog((xna)) + xca*rlog((xca)))
               + WKALS*r0*(1.0-sumr) + WNVKALS*r1*(1.0-sumr) + WKVKALS*r0*r1 + WKCKALS*r0*r2
               + PENALTY/r1;

    double dgdr0 = g2 - g1 + R*T*4.0*(rlog((xk))-rlog((xna)))
                 + WKALS*(1.0-2.0*r0-r1-r2) - WNVKALS*r1 + WKVKALS*r1 + WKCKALS*r2;
    double dgdr1 = g3 - g1 + R*T*(rlog((xvc-xca))-rlog((xna)))
                 - WKALS*r0 + WNVKALS*(1.0-r0-2.0*r1-r2) + WKVKALS*r0
                 - PENALTY/(r1*r1);
    double dgdr2 = g4 - g1 + G4
                 + R*T*(rlog((xca))-2.0*rlog((xna))-1.0)
                 - WKALS*r0 - WNVKALS*r1 + WKCKALS*r0;

    double dGex[4];
    dGex[0] = 0.0;
    dGex[1] = dgdr0;
    dGex[2] = dgdr1;
    dGex[3] = dgdr2;

    GH_mu_from_dGex(d, p, 1.0, 4, Gex, dGex);

    GH_sf_nph(p, d->sf);
    return GH_df_assemble(d, p, gb, mu_Gex, 4, grad);
}

typedef struct { GH_hsv c0, cS1, cS2, cX2, cX3, cX4, cX5, cX6, cX7, cS1S1, cS1S2, cS2S2, cX2S1, cX2S2, cX2X2, cX2X3, cX2X4, cX2X5, cX2X6, cX2X7, cX3S1, cX3S2, cX3X3, cX3X4, cX3X5, cX3X6, cX3X7, cX4S1, cX4S2, cX4X4, cX4X5, cX4X6, cX4X7, cX5S1, cX5S2, cX5X5, cX5X6, cX5X7, cX6S1, cX6S2, cX6X6, cX6X7, cX7S1, cX7S2, cX7X7, cX2S2S2, cX2X2S2, cX2X2X7, cX2X3S2, cX2X3X7, cX2X4S2, cX2X4X7, cX2X5S2, cX2X5X7, cX2X6S2, cX2X6X7, cX2X7S2, cX2X7X7, cX3S2S2, cX3X3S2, cX3X3X7, cX3X4S2, cX3X4X7, cX3X5S2, cX3X5X7, cX3X6S2, cX3X6X7, cX3X7S2, cX3X7X7, cX4S2S2, cX4X4S2, cX4X4X7, cX4X5S2, cX4X5X7, cX4X6S2, cX4X6X7, cX4X7S2, cX4X7X7, cX5S2S2, cX5X5S2, cX5X5X7, cX5X6S2, cX5X6X7, cX5X7S2, cX5X7X7, cX6S2S2, cX6X6S2, cX6X6X7, cX6X7S2, cX6X7X7, cX7S2S2, cX7X7S2, cX7X7X7; } GH_pyx_par;
typedef struct { double c0, cS1, cS2, cX2, cX3, cX4, cX5, cX6, cX7, cS1S1, cS1S2, cS2S2, cX2S1, cX2S2, cX2X2, cX2X3, cX2X4, cX2X5, cX2X6, cX2X7, cX3S1, cX3S2, cX3X3, cX3X4, cX3X5, cX3X6, cX3X7, cX4S1, cX4S2, cX4X4, cX4X5, cX4X6, cX4X7, cX5S1, cX5S2, cX5X5, cX5X6, cX5X7, cX6S1, cX6S2, cX6X6, cX6X7, cX7S1, cX7S2, cX7X7, cX2S2S2, cX2X2S2, cX2X2X7, cX2X3S2, cX2X3X7, cX2X4S2, cX2X4X7, cX2X5S2, cX2X5X7, cX2X6S2, cX2X6X7, cX2X7S2, cX2X7X7, cX3S2S2, cX3X3S2, cX3X3X7, cX3X4S2, cX3X4X7, cX3X5S2, cX3X5X7, cX3X6S2, cX3X6X7, cX3X7S2, cX3X7X7, cX4S2S2, cX4X4S2, cX4X4X7, cX4X5S2, cX4X5X7, cX4X6S2, cX4X6X7, cX4X7S2, cX4X7X7, cX5S2S2, cX5X5S2, cX5X5X7, cX5X6S2, cX5X6X7, cX5X7S2, cX5X7X7, cX6S2S2, cX6X6S2, cX6X6X7, cX6X7S2, cX6X7X7, cX7S2S2, cX7X7S2, cX7X7X7; } GH_pyx_g;

static const GH_pyx_par GH_PYX_CPX = {
    .c0 = {0.0, 0.0, 0.0},
    .cS1 = {-2677.12, -0.0, 0.0},
    .cS2 = {-9811.9766408, 0.0, -0.0064949},
    .cX2 = {7029.12, 0.0, 0.0},
    .cX3 = {16318.0, 0.0, 0.0},
    .cX4 = {18828.0, 0.0, 0.0},
    .cX5 = {23597.12, 0.0, 0.0},
    .cX6 = {-13220.68, 0.0, 0.0},
    .cX7 = {11610.1033592, -3.97551128, -0.0106789},
    .cS1S1 = {-8917.5, 0.0, 0.0},
    .cS1S2 = {7804.93, 0.0, 0.0},
    .cS2S2 = {7460.595, 0.0, -0.089956},
    .cX2S1 = {-3819.645, 0.0, 0.0},
    .cX2S2 = {-3117.08, 0.0, 0.087864},
    .cX2X2 = {-7029.12, 0.0, -0.0},
    .cX2X3 = {-8948.56, 0.0, -0.0},
    .cX2X4 = {-18247.52, 0.0, -0.0},
    .cX2X5 = {-20147.715, 0.0, -0.0},
    .cX2X6 = {-5424.3875, 0.0, -0.0},
    .cX2X7 = {657.9867184, -5.8088564, -0.0862098},
    .cX3S1 = {-19744.44, 0.0, 0.0},
    .cX3S2 = {7537.4816796, -0.54009164, 0.00878155},
    .cX3X3 = {-16318.0, 0.0, 0.0},
    .cX3X4 = {-18357.0, 0.0, 0.0},
    .cX3X5 = {-5465.44, 0.0, 0.0},
    .cX3X6 = {8160.4, 0.0, 0.0},
    .cX3X7 = {22496.3516796, 0.35523636, -0.07179645},
    .cX4S1 = {3347.06, 0.0, 0.0},
    .cX4S2 = {7968.9216796, -0.54009164, 0.00878155},
    .cX4X4 = {-18828.0, 0.0, 0.0},
    .cX4X5 = {-15951.94, 0.0, 0.0},
    .cX4X6 = {23252.15, 0.0, 0.0},
    .cX4X7 = {27147.9216796, 0.35523636, -0.07179645},
    .cX5S1 = {-2065.5, 0.0, 0.0},
    .cX5S2 = {-854.3466408, -1.08018328, 0.0091951},
    .cX5X5 = {-9937.0, 0.0, 0.0},
    .cX5X6 = {-16226.25, 0.0, 0.0},
    .cX5X7 = {17668.5133592, -0.18485528, -0.0713829},
    .cX6S1 = {-7322.0, 0.0, 0.0},
    .cX6S2 = {7437.9366796, -0.54009164, 0.00459755},
    .cX6X6 = {8735.875, 0.0, 0.0},
    .cX6X7 = {-651.4583204, -0.09242764, -0.03569145},
    .cX7S1 = {4240.94, 0.0, 0.0},
    .cX7S2 = {5460.6166408, 2.9044282, 0.01904688266},
    .cX7X7 = {-17564.4583592, 5.7997562, 0.04624288266},
    .cX2S2S2 = {169.975, 0.0, -0.005753},
    .cX2X2S2 = {-339.95, 0.0, 0.011506},
    .cX2X2X7 = {339.95, 0.0, -0.011506},
    .cX2X3S2 = {5177.7, 0.0, -0.043932},
    .cX2X3X7 = {-679.9, 0.0, 0.023012},
    .cX2X4S2 = {5177.7, 0.0, -0.043932},
    .cX2X4X7 = {-679.9, 0.0, 0.023012},
    .cX2X5S2 = {5177.7, 0.0, -0.043932},
    .cX2X5X7 = {-679.9, 0.0, 0.023012},
    .cX2X6S2 = {2588.85, 0.0, -0.021966},
    .cX2X6X7 = {-339.95, 0.0, 0.011506},
    .cX2X7S2 = {10695.35, 0.0, -0.09937},
    .cX2X7X7 = {-1012.005, 3.64848984, 0.08420296532},
    .cX3S2S2 = {-5628.7875, 0.0, 0.0509925},
    .cX3X3S2 = {1464.4, 0.0, -0.016736},
    .cX3X3X7 = {-1464.4, 0.0, 0.016736},
    .cX3X4S2 = {2366.575, 0.0, -0.030857},
    .cX3X4X7 = {-3491.025, 0.0, 0.036087},
    .cX3X5S2 = {3491.025, 0.0, -0.036087},
    .cX3X5X7 = {-2366.575, 0.0, 0.030857},
    .cX3X6S2 = {1464.4, 0.0, -0.016736},
    .cX3X6X7 = {-1464.4, 0.0, 0.016736},
    .cX3X7S2 = {1234.28, 0.91212246, -0.02928800867},
    .cX3X7X7 = {-1923.3325, 0.91212246, 0.02013549133},
    .cX4S2S2 = {-5628.7875, 0.0, 0.0509925},
    .cX4X4S2 = {1464.4, 0.0, -0.016736},
    .cX4X4X7 = {-1464.4, 0.0, 0.016736},
    .cX4X5S2 = {3491.025, 0.0, -0.036087},
    .cX4X5X7 = {-2366.575, 0.0, 0.030857},
    .cX4X6S2 = {2026.625, 0.0, -0.019351},
    .cX4X6X7 = {-902.175, 0.0, 0.014121},
    .cX4X7S2 = {1234.28, 0.91212246, -0.02928800867},
    .cX4X7X7 = {-1923.3325, 0.91212246, 0.02013549133},
    .cX5S2S2 = {-5347.675, 0.0, 0.049685},
    .cX5X5S2 = {2588.85, 0.0, -0.021966},
    .cX5X5X7 = {-339.95, 0.0, 0.011506},
    .cX5X6S2 = {2588.85, 0.0, -0.021966},
    .cX5X6X7 = {-339.95, 0.0, 0.011506},
    .cX5X7S2 = {3697.61, 1.82424492, -0.03451801734},
    .cX5X7X7 = {258.885, 1.82424492, 0.01621298266},
    .cX6S2S2 = {-2673.8375, 0.0, 0.0248425},
    .cX6X6S2 = {787.76875, 0.0, -0.00614525},
    .cX6X6X7 = {55.56875, 0.0, 0.00222275},
    .cX6X7S2 = {1848.805, 0.91212246, -0.01725900867},
    .cX6X7X7 = {129.4425, 0.91212246, 0.00810649133},
    .cX7S2S2 = {-10695.35, 0.0, 0.09937},
    .cX7X7S2 = {-4246.76, -1.82424492, -0.01255198266},
    .cX7X7X7 = {590.99, -1.82424492, -0.04497798266}
};

static const GH_pyx_par GH_PYX_OPX = {
    .c0 = {2331.7432000000003, 1.3807200000000002, 0.075312},
    .cS1 = {-2677.12, 0.0, 0.0},
    .cS2 = {-7635.799999999999, 0.0, -0.09523640000000068},
    .cX2 = {8086.416799999999, -0.8284320000000002, 0.03710980000000141},
    .cX3 = {39105.4968, -3.38072, -0.126602},
    .cX4 = {41615.4968, -3.38072, -0.126602},
    .cX5 = {46384.6168, -3.38072, -0.126602},
    .cX6 = {-1826.9316, -1.69036, -0.063301},
    .cX7 = {11203.4968, -4.78469688, -0.14879160000000066},
    .cS1S1 = {-8917.5, 0.0, 0.0},
    .cS1S2 = {7804.93, 0.0, 0.0},
    .cS2S2 = {15637.7, 0.0, -0.14957800000000002},
    .cX2S1 = {-3819.6450000000004, 0.0, 0.0},
    .cX2S2 = {-10773.800000000001, 0.0, 0.142256},
    .cX2X2 = {-8368.0, 0.0, -0.014121},
    .cX2X3 = {-9618.0, 0.0, -0.0070605},
    .cX2X4 = {-18916.96, 0.0, -0.0070605},
    .cX2X5 = {-21486.595, 0.0, -0.014121},
    .cX2X6 = {-6093.8275, 0.0, -0.0070605},
    .cX2X7 = {7050.039999999999, -1.15955376, -0.01317580000000132},
    .cX3S1 = {-19744.440000000002, 0.0, 0.0},
    .cX3S2 = {5934.734999999998, -0.28988844, 0.07267179999999968},
    .cX3X3 = {-16318.000000000002, 0.0, 0.0},
    .cX3X4 = {-18357.0, 0.0, 0.0},
    .cX3X5 = {-5465.440000000004, 0.0, 0.0},
    .cX3X6 = {8160.4000000000015, 0.0, 0.0},
    .cX3X7 = {-2484.9950000000044, -0.28988844, -0.06121620000000034},
    .cX4S1 = {3347.0599999999986, 0.0, 0.0},
    .cX4S2 = {6366.1749999999965, -0.28988844, 0.07267179999999968},
    .cX4X4 = {-18828.0, 0.0, 0.0},
    .cX4X5 = {-15951.939999999999, 0.0, 0.0},
    .cX4X6 = {23252.149999999994, 0.0, 0.0},
    .cX4X7 = {2166.574999999997, -0.28988844, -0.06121620000000034},
    .cX5S1 = {-2065.5, 0.0, 0.0},
    .cX5S2 = {-4133.060000000004, -0.57977688, 0.10141159999999935},
    .cX5X5 = {-9937.0, 0.0, 0.0},
    .cX5X6 = {-16226.249999999996, 0.0, 0.0},
    .cX5X7 = {-8988.800000000007, -0.57977688, -0.03247640000000067},
    .cX6S1 = {-7322.0, 0.0, 0.0},
    .cX6S2 = {5798.579999999997, -0.28988844, 0.050705799999999676},
    .cX6X6 = {8735.874999999996, 0.0, 0.0},
    .cX6X7 = {-13980.114999999998, -0.28988844, -0.016238200000000334},
    .cX7S1 = {4240.940000000002, 0.0, 0.0},
    .cX7S2 = {2782.3600000000006, 0.57977688, 0.022016400000000658},
    .cX7X7 = {-28168.78, 0.57977688, -0.01250159999999934},
    .cX2S2S2 = {209.19999999999982, 0.0, -0.007322},
    .cX2X2S2 = {-418.39999999999964, 0.0, 0.014644},
    .cX2X2X7 = {418.39999999999964, 0.0, -0.014644},
    .cX2X3S2 = {10041.6, 0.0, -0.071128},
    .cX2X3X7 = {-836.7999999999993, 0.0, 0.029288},
    .cX2X4S2 = {10041.6, 0.0, -0.071128},
    .cX2X4X7 = {-836.7999999999993, 0.0, 0.029288},
    .cX2X5S2 = {10041.6, 0.0, -0.071128},
    .cX2X5X7 = {-836.7999999999993, 0.0, 0.029288},
    .cX2X6S2 = {5020.8, 0.0, -0.035564},
    .cX2X6X7 = {-418.39999999999964, 0.0, 0.014644},
    .cX2X7S2 = {20501.6, 0.0, -0.15690000000000004},
    .cX2X7X7 = {-2301.199999999997, 0.0, 0.080542},
    .cX3S2S2 = {-10826.1, 0.0, 0.08106500000000001},
    .cX3X3S2 = {2719.6, 0.0, -0.025104},
    .cX3X3X7 = {-2719.6, 0.0, 0.025104},
    .cX3X4S2 = {4288.599999999999, 0.0, -0.044978000000000004},
    .cX3X4X7 = {-6589.8, 0.0, 0.055438},
    .cX3X5S2 = {6589.8, 0.0, -0.055438},
    .cX3X5X7 = {-4288.599999999999, 0.0, 0.044978000000000004},
    .cX3X6S2 = {2719.6, 0.0, -0.025104},
    .cX3X6X7 = {-2719.6, 0.0, 0.025104},
    .cX3X7S2 = {1255.199999999999, 0.0, -0.043932},
    .cX3X7X7 = {-4236.3, 0.0, 0.025627},
    .cX4S2S2 = {-10826.1, 0.0, 0.08106500000000001},
    .cX4X4S2 = {2719.6, 0.0, -0.025104},
    .cX4X4X7 = {-2719.6, 0.0, 0.025104},
    .cX4X5S2 = {6589.8, 0.0, -0.055438},
    .cX4X5X7 = {-4288.599999999999, 0.0, 0.044978000000000004},
    .cX4X6S2 = {3870.2, 0.0, -0.030334},
    .cX4X6X7 = {-1568.9999999999998, 0.0, 0.019874000000000003},
    .cX4X7S2 = {1255.199999999999, 0.0, -0.043932},
    .cX4X7X7 = {-4236.3, 0.0, 0.025627},
    .cX5S2S2 = {-10250.8, 0.0, 0.07845000000000002},
    .cX5X5S2 = {5020.8, 0.0, -0.035564},
    .cX5X5X7 = {-418.39999999999964, 0.0, 0.014644},
    .cX5X6S2 = {5020.8, 0.0, -0.035564},
    .cX5X6X7 = {-418.39999999999964, 0.0, 0.014644},
    .cX5X7S2 = {5857.5999999999985, 0.0, -0.064852},
    .cX5X7X7 = {-209.19999999999982, 0.0, 0.007322},
    .cX6S2S2 = {-5125.4, 0.0, 0.03922500000000001},
    .cX6X6S2 = {1542.85, 0.0, -0.010198500000000001},
    .cX6X6X7 = {183.05000000000007, 0.0, 0.0023534999999999997},
    .cX6X7S2 = {2928.7999999999993, 0.0, -0.032426},
    .cX6X7X7 = {-104.59999999999991, 0.0, 0.003661},
    .cX7S2S2 = {-20501.6, 0.0, 0.15690000000000004},
    .cX7X7S2 = {-8368.0, 0.0, 0.012552000000000008},
    .cX7X7X7 = {1255.199999999999, 0.0, -0.043932}
};

static void GH_pyx_fold(const GH_pyx_par *P, double T, double Pv, double SS2, GH_pyx_g *g){
    g->c0 = GH_hsv_G(P->c0, T, Pv);
    g->cS1 = GH_hsv_G(P->cS1, T, Pv);
    g->cS2 = P->cS2.h - T*SS2 + Pv*P->cS2.v;
    g->cX2 = GH_hsv_G(P->cX2, T, Pv);
    g->cX3 = GH_hsv_G(P->cX3, T, Pv);
    g->cX4 = GH_hsv_G(P->cX4, T, Pv);
    g->cX5 = GH_hsv_G(P->cX5, T, Pv);
    g->cX6 = GH_hsv_G(P->cX6, T, Pv);
    g->cX7 = GH_hsv_G(P->cX7, T, Pv);
    g->cS1S1 = GH_hsv_G(P->cS1S1, T, Pv);
    g->cS1S2 = GH_hsv_G(P->cS1S2, T, Pv);
    g->cS2S2 = GH_hsv_G(P->cS2S2, T, Pv);
    g->cX2S1 = GH_hsv_G(P->cX2S1, T, Pv);
    g->cX2S2 = GH_hsv_G(P->cX2S2, T, Pv);
    g->cX2X2 = GH_hsv_G(P->cX2X2, T, Pv);
    g->cX2X3 = GH_hsv_G(P->cX2X3, T, Pv);
    g->cX2X4 = GH_hsv_G(P->cX2X4, T, Pv);
    g->cX2X5 = GH_hsv_G(P->cX2X5, T, Pv);
    g->cX2X6 = GH_hsv_G(P->cX2X6, T, Pv);
    g->cX2X7 = GH_hsv_G(P->cX2X7, T, Pv);
    g->cX3S1 = GH_hsv_G(P->cX3S1, T, Pv);
    g->cX3S2 = GH_hsv_G(P->cX3S2, T, Pv);
    g->cX3X3 = GH_hsv_G(P->cX3X3, T, Pv);
    g->cX3X4 = GH_hsv_G(P->cX3X4, T, Pv);
    g->cX3X5 = GH_hsv_G(P->cX3X5, T, Pv);
    g->cX3X6 = GH_hsv_G(P->cX3X6, T, Pv);
    g->cX3X7 = GH_hsv_G(P->cX3X7, T, Pv);
    g->cX4S1 = GH_hsv_G(P->cX4S1, T, Pv);
    g->cX4S2 = GH_hsv_G(P->cX4S2, T, Pv);
    g->cX4X4 = GH_hsv_G(P->cX4X4, T, Pv);
    g->cX4X5 = GH_hsv_G(P->cX4X5, T, Pv);
    g->cX4X6 = GH_hsv_G(P->cX4X6, T, Pv);
    g->cX4X7 = GH_hsv_G(P->cX4X7, T, Pv);
    g->cX5S1 = GH_hsv_G(P->cX5S1, T, Pv);
    g->cX5S2 = GH_hsv_G(P->cX5S2, T, Pv);
    g->cX5X5 = GH_hsv_G(P->cX5X5, T, Pv);
    g->cX5X6 = GH_hsv_G(P->cX5X6, T, Pv);
    g->cX5X7 = GH_hsv_G(P->cX5X7, T, Pv);
    g->cX6S1 = GH_hsv_G(P->cX6S1, T, Pv);
    g->cX6S2 = GH_hsv_G(P->cX6S2, T, Pv);
    g->cX6X6 = GH_hsv_G(P->cX6X6, T, Pv);
    g->cX6X7 = GH_hsv_G(P->cX6X7, T, Pv);
    g->cX7S1 = GH_hsv_G(P->cX7S1, T, Pv);
    g->cX7S2 = GH_hsv_G(P->cX7S2, T, Pv);
    g->cX7X7 = GH_hsv_G(P->cX7X7, T, Pv);
    g->cX2S2S2 = GH_hsv_G(P->cX2S2S2, T, Pv);
    g->cX2X2S2 = GH_hsv_G(P->cX2X2S2, T, Pv);
    g->cX2X2X7 = GH_hsv_G(P->cX2X2X7, T, Pv);
    g->cX2X3S2 = GH_hsv_G(P->cX2X3S2, T, Pv);
    g->cX2X3X7 = GH_hsv_G(P->cX2X3X7, T, Pv);
    g->cX2X4S2 = GH_hsv_G(P->cX2X4S2, T, Pv);
    g->cX2X4X7 = GH_hsv_G(P->cX2X4X7, T, Pv);
    g->cX2X5S2 = GH_hsv_G(P->cX2X5S2, T, Pv);
    g->cX2X5X7 = GH_hsv_G(P->cX2X5X7, T, Pv);
    g->cX2X6S2 = GH_hsv_G(P->cX2X6S2, T, Pv);
    g->cX2X6X7 = GH_hsv_G(P->cX2X6X7, T, Pv);
    g->cX2X7S2 = GH_hsv_G(P->cX2X7S2, T, Pv);
    g->cX2X7X7 = GH_hsv_G(P->cX2X7X7, T, Pv);
    g->cX3S2S2 = GH_hsv_G(P->cX3S2S2, T, Pv);
    g->cX3X3S2 = GH_hsv_G(P->cX3X3S2, T, Pv);
    g->cX3X3X7 = GH_hsv_G(P->cX3X3X7, T, Pv);
    g->cX3X4S2 = GH_hsv_G(P->cX3X4S2, T, Pv);
    g->cX3X4X7 = GH_hsv_G(P->cX3X4X7, T, Pv);
    g->cX3X5S2 = GH_hsv_G(P->cX3X5S2, T, Pv);
    g->cX3X5X7 = GH_hsv_G(P->cX3X5X7, T, Pv);
    g->cX3X6S2 = GH_hsv_G(P->cX3X6S2, T, Pv);
    g->cX3X6X7 = GH_hsv_G(P->cX3X6X7, T, Pv);
    g->cX3X7S2 = GH_hsv_G(P->cX3X7S2, T, Pv);
    g->cX3X7X7 = GH_hsv_G(P->cX3X7X7, T, Pv);
    g->cX4S2S2 = GH_hsv_G(P->cX4S2S2, T, Pv);
    g->cX4X4S2 = GH_hsv_G(P->cX4X4S2, T, Pv);
    g->cX4X4X7 = GH_hsv_G(P->cX4X4X7, T, Pv);
    g->cX4X5S2 = GH_hsv_G(P->cX4X5S2, T, Pv);
    g->cX4X5X7 = GH_hsv_G(P->cX4X5X7, T, Pv);
    g->cX4X6S2 = GH_hsv_G(P->cX4X6S2, T, Pv);
    g->cX4X6X7 = GH_hsv_G(P->cX4X6X7, T, Pv);
    g->cX4X7S2 = GH_hsv_G(P->cX4X7S2, T, Pv);
    g->cX4X7X7 = GH_hsv_G(P->cX4X7X7, T, Pv);
    g->cX5S2S2 = GH_hsv_G(P->cX5S2S2, T, Pv);
    g->cX5X5S2 = GH_hsv_G(P->cX5X5S2, T, Pv);
    g->cX5X5X7 = GH_hsv_G(P->cX5X5X7, T, Pv);
    g->cX5X6S2 = GH_hsv_G(P->cX5X6S2, T, Pv);
    g->cX5X6X7 = GH_hsv_G(P->cX5X6X7, T, Pv);
    g->cX5X7S2 = GH_hsv_G(P->cX5X7S2, T, Pv);
    g->cX5X7X7 = GH_hsv_G(P->cX5X7X7, T, Pv);
    g->cX6S2S2 = GH_hsv_G(P->cX6S2S2, T, Pv);
    g->cX6X6S2 = GH_hsv_G(P->cX6X6S2, T, Pv);
    g->cX6X6X7 = GH_hsv_G(P->cX6X6X7, T, Pv);
    g->cX6X7S2 = GH_hsv_G(P->cX6X7S2, T, Pv);
    g->cX6X7X7 = GH_hsv_G(P->cX6X7X7, T, Pv);
    g->cX7S2S2 = GH_hsv_G(P->cX7S2S2, T, Pv);
    g->cX7X7S2 = GH_hsv_G(P->cX7X7S2, T, Pv);
    g->cX7X7X7 = GH_hsv_G(P->cX7X7X7, T, Pv);
}

static void GH_pyx_vertex_G(const GH_pyx_g *g, double ends[7]){
    ends[0] = g->c0;
    ends[1] = g->c0 + g->cX7 - g->cS2 + g->cX7X7 - g->cX7S2 + g->cS2S2 + g->cX7X7X7 - g->cX7X7S2 + g->cX7S2S2;
    ends[2] = g->c0 + g->cX2 + g->cX2X2;
    ends[3] = g->c0 + g->cX3 + g->cX3X3;
    ends[4] = g->c0 + g->cX4 + g->cX4X4;
    ends[5] = 0.0;
    ends[6] = g->c0
            + 0.5*g->cX3
            - 0.5*g->cX4
            + 0.5*g->cX5
            +     g->cX6
            -     g->cS1
            + 0.25*g->cX3X3
            - 0.25*g->cX3X4
            + 0.25*g->cX3X5
            + 0.5*g->cX3X6
            - 0.5*g->cX3S1
            + 0.25*g->cX4X4
            - 0.25*g->cX4X5
            - 0.5*g->cX4X6
            + 0.5*g->cX4S1
            + 0.25*g->cX5X5
            + 0.5*g->cX5X6
            - 0.5*g->cX5S1
            +     g->cX6X6
            -     g->cX6S1
            +     g->cS1S1;
}
static void GH_pyx_site_fracs(double r0,double r1,double r2,double r3,double r4,double r5,double s0,double s1,
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

static void GH_pyx_composites(double xmg2m1,double xfe2m1,double xti4m1,double xna1m2,double xsi4tet,
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

static void GH_pyx_dgds(double r0,double r1,double r2,double r3,double r4,double r5,double s0,double s1,
                         double T,double Rgas,const GH_pyx_g *g,double *dgds0,double *dgds1){
    double xal3m1,xfe2m1,xfe3m1,xmg2m1,xti4m1,xca2m2,xfe2m2,xmg2m2,xna1m2,xal3tet,xfe3tet,xsi4tet;
    GH_pyx_site_fracs(r0,r1,r2,r3,r4,r5,s0,s1,&xal3m1,&xfe2m1,&xfe3m1,&xmg2m1,&xti4m1,&xca2m2,&xfe2m2,&xmg2m2,&xna1m2,&xal3tet,&xfe3tet,&xsi4tet);
    *dgds0 = 0.5*Rgas*T*rlog(((xfe3m1*xal3tet)/(xal3m1*xfe3tet))) +                g->cS1 +                g->cX2S1*r0 +                g->cX3S1*r1 +                g->cX4S1*r2 +                g->cX5S1*r3 +                g->cX6S1*r4 +                g->cX7S1*r5 +                g->cS1S1*s0*2.0 +                g->cS1S2*s1;
    *dgds1 = 0.5*Rgas*T*rlog(((xmg2m1*xfe2m2)/(xfe2m1*xmg2m2))) +                g->cS2 +                g->cX2S2*r0 +                g->cX3S2*r1 +                g->cX4S2*r2 +                g->cX5S2*r3 +                g->cX6S2*r4 +                g->cX7S2*r5 +                g->cS1S2*s0 +                g->cS2S2*s1*2.0 +                g->cX2X2S2*r0*r0 +                g->cX2X3S2*r0*r1 +                g->cX2X4S2*r0*r2 +                g->cX2X5S2*r0*r3 +                g->cX2X6S2*r0*r4 +                g->cX2X7S2*r0*r5 +                g->cX2S2S2*r0*s1*2.0 +                g->cX3X3S2*r1*r1 +                g->cX3X4S2*r1*r2 +                g->cX3X5S2*r1*r3 +                g->cX3X6S2*r1*r4 +                g->cX3X7S2*r1*r5 +                g->cX3S2S2*r1*s1*2.0 +                g->cX4X4S2*r2*r2 +                g->cX4X5S2*r2*r3 +                g->cX4X6S2*r2*r4 +                g->cX4X7S2*r2*r5 +                g->cX4S2S2*r2*s1*2.0 +                g->cX5X5S2*r3*r3 +                g->cX5X6S2*r3*r4 +                g->cX5X7S2*r3*r5 +                g->cX5S2S2*r3*s1*2.0 +                g->cX6X6S2*r4*r4 +                g->cX6X7S2*r4*r5 +                g->cX6S2S2*r4*s1*2.0 +                g->cX7X7S2*r5*r5 +                g->cX7S2S2*r5*s1*2.0;
}

static void GH_bound_track(double val0, const double dval[6], double slope, double eps,
                           double *lo, double *hi, double dlo[6], double dhi[6]){
    double a, b;
    if (slope > 0.0){ a = (eps-val0)/slope; b = (1.0-eps-val0)/slope; }
    else if (slope < 0.0){ a = (1.0-eps-val0)/slope; b = (eps-val0)/slope; }
    else { return; }
    if (a > *lo){ *lo = a; for (int j = 0; j < 6; j++){ dlo[j] = -dval[j]/slope; } }
    if (b < *hi){ *hi = b; for (int j = 0; j < 6; j++){ dhi[j] = -dval[j]/slope; } }
}

static void GH_pyx_bound_chain(const double r[6], double s0, double s1, double g0, double g1, double eps, double dGdr[6]){
    double lo, hi, dlo[6], dhi[6];
    lo = -1.0e6; hi = 1.0e6;
    for (int j = 0; j < 6; j++){ dlo[j] = 0.0; dhi[j] = 0.0; }
    { double d[6] = {0,0,0,0.5,0.25,0};     GH_bound_track((2.0*r[3]+r[4])/4.0, d, -0.5, eps, &lo, &hi, dlo, dhi); }
    { double d[6] = {0,0,0,0.5,0.25,0};     GH_bound_track((2.0*r[3]+r[4])/4.0, d,  0.5, eps, &lo, &hi, dlo, dhi); }
    { double d[6] = {0,0.5,0,0.25,-0.125,0}; GH_bound_track((4.0*r[1]+2.0*r[3]-r[4])/8.0, d, 0.25, eps, &lo, &hi, dlo, dhi); }
    { double d[6] = {0,0,0.5,0.25,-0.125,0}; GH_bound_track((4.0*r[2]+2.0*r[3]-r[4])/8.0, d, -0.25, eps, &lo, &hi, dlo, dhi); }
    double tol0 = 1.0e-12*(1.0 + fabs(s0));
    if (fabs(s0 - lo) <= tol0 && (g0 > 0.0 || fabs(s0 - hi) > tol0)){ for (int j = 0; j < 6; j++){ dGdr[j] += g0*dlo[j]; } }
    else if (fabs(s0 - hi) <= tol0){ for (int j = 0; j < 6; j++){ dGdr[j] += g0*dhi[j]; } }

    lo = -1.0e6; hi = 1.0e6;
    for (int j = 0; j < 6; j++){ dlo[j] = 0.0; dhi[j] = 0.0; }
    { double d[6] = {1.0,0,0,0,0,-0.5};          GH_bound_track((2.0*r[0]-r[5])/2.0, d, -0.5, eps, &lo, &hi, dlo, dhi); }
    { double d[6] = {-1.0,-0.5,-0.5,-1.0,-0.5,0.5}; GH_bound_track(1.0-r[0]-r[3]-0.5*(r[1]+r[2]+r[4]-r[5]), d, 0.5, eps, &lo, &hi, dlo, dhi); }
    { double d[6] = {0,0,0,0,0,0.5};             GH_bound_track(r[5]/2.0, d,  0.5, eps, &lo, &hi, dlo, dhi); }
    { double d[6] = {0,0,0,0,0,0.5};             GH_bound_track(r[5]/2.0, d, -0.5, eps, &lo, &hi, dlo, dhi); }
    double tol1 = 1.0e-12*(1.0 + fabs(s1));
    if (fabs(s1 - lo) <= tol1 && (g1 > 0.0 || fabs(s1 - hi) > tol1)){ for (int j = 0; j < 6; j++){ dGdr[j] += g1*dlo[j]; } }
    else if (fabs(s1 - hi) <= tol1){ for (int j = 0; j < 6; j++){ dGdr[j] += g1*dhi[j]; } }
}

static void GH_pyx_bounds_s0(double r1,double r2,double r3,double r4,double eps,double *lo,double *hi){
    *lo = -1.0e6; *hi = 1.0e6;
    double val, slope;
    val = (2.0*r3+r4)/4.0;      slope = -0.5;  GH_narrow_bound(val,slope,eps,lo,hi);
    val = (2.0*r3+r4)/4.0;      slope =  0.5;  GH_narrow_bound(val,slope,eps,lo,hi);
    val = (4.0*r1+2.0*r3-r4)/8.0; slope = 0.25; GH_narrow_bound(val,slope,eps,lo,hi);
    val = (4.0*r2+2.0*r3-r4)/8.0; slope = -0.25; GH_narrow_bound(val,slope,eps,lo,hi);
}
static void GH_pyx_bounds_s1(double r0,double r1,double r2,double r3,double r4,double r5,double eps,double *lo,double *hi){
    *lo = -1.0e6; *hi = 1.0e6;
    double val, slope;
    val = (2.0*r0-r5)/2.0;                        slope = -0.5; GH_narrow_bound(val,slope,eps,lo,hi);
    val = 1.0-r0-r3-0.5*(r1+r2+r4-r5);            slope =  0.5; GH_narrow_bound(val,slope,eps,lo,hi);
    val = r5/2.0;                                  slope =  0.5; GH_narrow_bound(val,slope,eps,lo,hi);
    val = r5/2.0;                                  slope = -0.5; GH_narrow_bound(val,slope,eps,lo,hi);
}

static void GH_pyx_d2gds(double r0,double r1,double r2,double r3,double r4,double r5,double s0,double s1,
                          double T,double Rgas,const GH_pyx_g *g,double *d00,double *d01,double *d11){
    double xal3m1,xfe2m1,xfe3m1,xmg2m1,xti4m1,xca2m2,xfe2m2,xmg2m2,xna1m2,xal3tet,xfe3tet,xsi4tet;
    GH_pyx_site_fracs(r0,r1,r2,r3,r4,r5,s0,s1,&xal3m1,&xfe2m1,&xfe3m1,&xmg2m1,&xti4m1,&xca2m2,&xfe2m2,&xmg2m2,&xna1m2,&xal3tet,&xfe3tet,&xsi4tet);
    *d00 = Rgas*T*(0.25/xfe3m1+0.25/xal3m1+0.125/xal3tet+0.125/xfe3tet) + 2.0*g->cS1S1;
    *d01 = g->cS1S2;
    *d11 = Rgas*T*(0.25/xmg2m1+0.25/xfe2m1+0.25/xfe2m2+0.25/xmg2m2) + 2.0*g->cS2S2
         + 2.0*r0*g->cX2S2S2 + 2.0*r1*g->cX3S2S2 + 2.0*r2*g->cX4S2S2
         + 2.0*r3*g->cX5S2S2 + 2.0*r4*g->cX6S2S2 + 2.0*r5*g->cX7S2S2;
}

static void GH_pyx_solve_s_from(double r0,double r1,double r2,double r3,double r4,double r5,double T,double Rgas,const GH_pyx_g *g,
                                 double s0init,double s1init,double *s0o,double *s1o){
    const int maxIter = 1000;
    double eps_s = 1.0e-8;
    double s0=s0init, s1=s1init;
    double s0Old=2.0, s1Old=2.0;
    int iter = 0;
    while ( (fabs(s0-s0Old) > 10.0*2.220446049250313e-16 ||
             fabs(s1-s1Old) > 10.0*2.220446049250313e-16) && iter < maxIter ){
        double dgds0, dgds1;
        GH_pyx_dgds(r0,r1,r2,r3,r4,r5,s0,s1,T,Rgas,g,&dgds0,&dgds1);
        double d00, d01, d11;
        GH_pyx_d2gds(r0,r1,r2,r3,r4,r5,s0,s1,T,Rgas,g,&d00,&d01,&d11);

        s0Old = s0; s1Old = s1;

        double det = d00*d11 - d01*d01;
        double deltaS0 = -( d11*dgds0 - d01*dgds1)/det;
        double deltaS1 = -(-d01*dgds0 + d00*dgds1)/det;

        double sN0 = s0Old + deltaS0, sN1 = s1Old + deltaS1;
        double lo, hi;
        GH_pyx_bounds_s0(r1,r2,r3,r4,eps_s,&lo,&hi);
        if (sN0 < lo){ sN0 = lo; } if (sN0 > hi){ sN0 = hi; }
        GH_pyx_bounds_s1(r0,r1,r2,r3,r4,r5,eps_s,&lo,&hi);
        if (sN1 < lo){ sN1 = lo; } if (sN1 > hi){ sN1 = hi; }
        s0 = sN0; s1 = sN1;
        iter++;
    }
    *s0o = s0; *s1o = s1;
}

static void GH_pyx_guess_s(double r0,double r1,double r2,double r3,double r4,double r5,double *s0g,double *s1g){
    double totAl=r1+r3, totCa=1.0-r4-r5, totFe2=r0, totFe3=r2+r3;
    double totMg=1.0-r0-0.5*(r1+r2)-r3-0.5*r4+r5, totNa=r4, totTi=0.5*(r1+r2);
    double totM2=totCa+totNa, totM1=totTi+(totMg+totFe2-(1.0-totM2));
    *s0g = (totFe3+totAl != 0.0) ? (1.0-totM1)*(totFe3-totAl)/(totFe3+totAl) : 0.0;
    *s1g = (totFe2+totMg != 0.0) ? (1.0-totM2)*(totFe2-totMg)/(totFe2+totMg) : 0.0;
}
static double GH_pyx_G_at(double r0,double r1,double r2,double r3,double r4,double r5,double s0,double s1,
                           double T,double Rgas,const GH_pyx_g *g){
    double xal3m1,xfe2m1,xfe3m1,xmg2m1,xti4m1,xca2m2,xfe2m2,xmg2m2,xna1m2,xal3tet,xfe3tet,xsi4tet;
    GH_pyx_site_fracs(r0,r1,r2,r3,r4,r5,s0,s1,&xal3m1,&xfe2m1,&xfe3m1,&xmg2m1,&xti4m1,&xca2m2,&xfe2m2,&xmg2m2,&xna1m2,&xal3tet,&xfe3tet,&xsi4tet);
    double xm1_MgFeNoTi,xm1_notTi,xm1_MgFe,xtet_notSi,xrest_notNa,xm1_notMgFe,xm2_notNa;
    GH_pyx_composites(xmg2m1,xfe2m1,xti4m1,xna1m2,xsi4tet,&xm1_MgFeNoTi,&xm1_notTi,&xm1_MgFe,&xtet_notSi,&xrest_notNa,&xm1_notMgFe,&xm2_notNa);
    double SIC = -Rgas*(xmg2m1*rlog((xmg2m1))                + xfe2m1*rlog((xfe2m1))                + xal3m1*rlog((xal3m1))                + xfe3m1*rlog((xfe3m1))                + xti4m1*rlog((xti4m1))                + (xm1_MgFeNoTi)*rlog((xm1_MgFeNoTi))                - (xm1_notTi)*rlog((xm1_notTi))                - (xm1_MgFe)*rlog((xm1_MgFe))                - 2.0*(xtet_notSi)*rlog((xtet_notSi))                + 2.0*xal3tet*rlog((xal3tet))                + 2.0*xfe3tet*rlog((xfe3tet))                + xca2m2*rlog((xca2m2))                + xna1m2*rlog((xna1m2))                + xmg2m2*rlog((xmg2m2))                + xfe2m2*rlog((xfe2m2))                + (xrest_notNa)*rlog((xrest_notNa))                - (xm1_notMgFe)*rlog((xm1_notMgFe))                - (xm2_notNa)*rlog((xm2_notNa)) );
    return -T*(SIC) + g->c0 +              g->cX2*r0 +              g->cX3*r1 +              g->cX4*r2 +              g->cX5*r3 +              g->cX6*r4 +              g->cX7*r5 +              g->cS1*s0 +              g->cS2*s1 +              g->cX2X2*r0*r0 +              g->cX2X3*r0*r1 +              g->cX2X4*r0*r2 +              g->cX2X5*r0*r3 +              g->cX2X6*r0*r4 +              g->cX2X7*r0*r5 +              g->cX2S1*r0*s0 +              g->cX2S2*r0*s1 +              g->cX3X3*r1*r1 +              g->cX3X4*r1*r2 +              g->cX3X5*r1*r3 +              g->cX3X6*r1*r4 +              g->cX3X7*r1*r5 +              g->cX3S1*r1*s0 +              g->cX3S2*r1*s1 +              g->cX4X4*r2*r2 +              g->cX4X5*r2*r3 +              g->cX4X6*r2*r4 +              g->cX4X7*r2*r5 +              g->cX4S1*r2*s0 +              g->cX4S2*r2*s1 +              g->cX5X5*r3*r3 +              g->cX5X6*r3*r4 +              g->cX5X7*r3*r5 +              g->cX5S1*r3*s0 +              g->cX5S2*r3*s1 +              g->cX6X6*r4*r4 +              g->cX6X7*r4*r5 +              g->cX6S1*r4*s0 +              g->cX6S2*r4*s1 +              g->cX7X7*r5*r5 +              g->cX7S1*r5*s0 +              g->cX7S2*r5*s1 +              g->cS1S1*s0*s0 +              g->cS1S2*s0*s1 +              g->cS2S2*s1*s1 +              g->cX2X2X7*r0*r0*r5 +              g->cX2X2S2*r0*r0*s1 +              g->cX2X3X7*r0*r1*r5 +              g->cX2X3S2*r0*r1*s1 +              g->cX2X4X7*r0*r2*r5 +              g->cX2X4S2*r0*r2*s1 +              g->cX2X5X7*r0*r3*r5 +              g->cX2X5S2*r0*r3*s1 +              g->cX2X6X7*r0*r4*r5 +              g->cX2X6S2*r0*r4*s1 +              g->cX2X7X7*r0*r5*r5 +              g->cX2X7S2*r0*r5*s1 +              g->cX2S2S2*r0*s1*s1 +              g->cX3X3X7*r1*r1*r5 +              g->cX3X3S2*r1*r1*s1 +              g->cX3X4X7*r1*r2*r5 +              g->cX3X4S2*r1*r2*s1 +              g->cX3X5X7*r1*r3*r5 +              g->cX3X5S2*r1*r3*s1 +              g->cX3X6X7*r1*r4*r5 +              g->cX3X6S2*r1*r4*s1 +              g->cX3X7X7*r1*r5*r5 +              g->cX3X7S2*r1*r5*s1 +              g->cX3S2S2*r1*s1*s1 +              g->cX4X4X7*r2*r2*r5 +              g->cX4X4S2*r2*r2*s1 +              g->cX4X5X7*r2*r3*r5 +              g->cX4X5S2*r2*r3*s1 +              g->cX4X6X7*r2*r4*r5 +              g->cX4X6S2*r2*r4*s1 +              g->cX4X7X7*r2*r5*r5 +              g->cX4X7S2*r2*r5*s1 +              g->cX4S2S2*r2*s1*s1 +              g->cX5X5X7*r3*r3*r5 +              g->cX5X5S2*r3*r3*s1 +              g->cX5X6X7*r3*r4*r5 +              g->cX5X6S2*r3*r4*s1 +              g->cX5X7X7*r3*r5*r5 +              g->cX5X7S2*r3*r5*s1 +              g->cX5S2S2*r3*s1*s1 +              g->cX6X6X7*r4*r4*r5 +              g->cX6X6S2*r4*r4*s1 +              g->cX6X7X7*r4*r5*r5 +              g->cX6X7S2*r4*r5*s1 +              g->cX6S2S2*r4*s1*s1 +              g->cX7X7X7*r5*r5*r5 +              g->cX7X7S2*r5*r5*s1 +              g->cX7S2S2*r5*s1*s1;
}

static double GH_pyx_solve_and_G(double r0,double r1,double r2,double r3,double r4,double r5,
                                  double T,double Rgas,const GH_pyx_g *g,double *s0o,double *s1o){
    double s0g,s1g; GH_pyx_guess_s(r0,r1,r2,r3,r4,r5,&s0g,&s1g);
    double cs0,cs1; GH_pyx_solve_s_from(r0,r1,r2,r3,r4,r5,T,Rgas,g,s0g,s1g,&cs0,&cs1);
    double G = GH_pyx_G_at(r0,r1,r2,r3,r4,r5,cs0,cs1,T,Rgas,g);
    double lo0,hi0,lo1,hi1;
    GH_pyx_bounds_s0(r1,r2,r3,r4,1.0e-8,&lo0,&hi0);
    GH_pyx_bounds_s1(r0,r1,r2,r3,r4,r5,1.0e-8,&lo1,&hi1);
    double gs0,gs1;
    GH_pyx_dgds(r0,r1,r2,r3,r4,r5,cs0,cs1,T,Rgas,g,&gs0,&gs1);
    int free0 = (cs0 > lo0 + 1.0e-6*(hi0-lo0)) && (cs0 < hi0 - 1.0e-6*(hi0-lo0));
    int free1 = (cs1 > lo1 + 1.0e-6*(hi1-lo1)) && (cs1 < hi1 - 1.0e-6*(hi1-lo1));
    if ((free0 && fabs(gs0) > 1.0) || (free1 && fabs(gs1) > 1.0)){
        const int N = 20;
        double bG = 1.0e300, b0 = s0g, b1 = s1g;
        for (int i = 0; i <= N; i++){
            for (int j = 0; j <= N; j++){
                double a = lo0 + (hi0-lo0)*i/N, b = lo1 + (hi1-lo1)*j/N;
                double Gt = GH_pyx_G_at(r0,r1,r2,r3,r4,r5,a,b,T,Rgas,g);
                if (Gt < bG){ bG = Gt; b0 = a; b1 = b; }
            }
        }
        double ts0,ts1; GH_pyx_solve_s_from(r0,r1,r2,r3,r4,r5,T,Rgas,g,b0,b1,&ts0,&ts1);
        double Gt = GH_pyx_G_at(r0,r1,r2,r3,r4,r5,ts0,ts1,T,Rgas,g);
        if (Gt < G){ G = Gt; cs0 = ts0; cs1 = ts1; }
    }
    *s0o = cs0; *s1o = cs1;
    return G;
}

static double GH_pyx_pure_ES_dgds(double s,double T,double Rgas,const GH_pyx_g *g){
    return Rgas*T*(rlog((1.0+s))-rlog((1.0-s))) + g->cS1 + g->cX5S1 + g->cS1S1*s*2.0;
}
static double GH_pyx_pure_ES_G_at_s(double s,double T,double Rgas,const GH_pyx_g *g){
    return Rgas*T*((1.0-s)*rlog((1.0-s)) + (1.0+s)*rlog((1.0+s)) - 2.0*rlog((2.0)))
         + g->c0 + g->cX5 + g->cS1*s + g->cX5X5 + g->cX5S1*s + g->cS1S1*s*s;
}
static double GH_pyx_pure_ES_G_compute(double T,double Rgas,const GH_pyx_g *g){
    double eps = 1.0e-8, lo = -1.0+eps, hi = 1.0-eps;
    double fa = GH_pyx_pure_ES_dgds(lo,T,Rgas,g), fb = GH_pyx_pure_ES_dgds(hi,T,Rgas,g);
    double a = lo, b = hi;
    if (fa > 0.0 && fb > 0.0){ a = b = lo; }
    else if (fa < 0.0 && fb < 0.0){ a = b = hi; }
    else {
        for (int it = 0; it < 80; it++){
            double mid = 0.5*(a+b), fm = GH_pyx_pure_ES_dgds(mid,T,Rgas,g);
            if ((fa > 0.0) == (fm > 0.0)){ a = mid; fa = fm; } else { b = mid; }
            if (b-a < 1.0e-14){ break; }
        }
    }
    return GH_pyx_pure_ES_G_at_s(0.5*(a+b),T,Rgas,g);
}

static double GH_pyx_pure_ES_G(double T,double Rgas,double Pv,const GH_pyx_g *g){
    static __thread double c_T, c_R, c_P, c_G;
    static __thread int    c_ok = 0;
    if (c_ok && c_T == T && c_R == Rgas && c_P == Pv){ return c_G; }
    c_G = GH_pyx_pure_ES_G_compute(T,Rgas,g);
    c_T = T; c_R = Rgas; c_P = Pv; c_ok = 1;
    return c_G;
}

static double GH_pyx_obj(unsigned n, const double *x, double *grad, void *SS_ref_db, const GH_pyx_par *P, double SS2x, int opx){
    SS_ref *d = (SS_ref *) SS_ref_db;

    double T    = d->T;
    double Rgas = GH_R_J;
    double Pv   = d->P*1000.0 - 1.0;
    double SS2  = (d->EM_database == 0) ? SS2x : 0.0;
    double SS2v = (d->EM_database == 0) ? -1.08018328 : 0.0;
    double *p   = d->p;
    double *gb  = d->gb_lvl;
    double *mu_Gex = d->mu_Gex;

    GH_pyx_g g, gv;
    GH_pyx_fold(P, T, Pv, SS2, &g);
    if (opx){ GH_pyx_fold(&GH_PYX_CPX, T, Pv, SS2v, &gv); }
    const GH_pyx_g *gvp = opx ? &gv : &g;

    for (int i = 0; i < 7; i++){ p[i] = x[i]; }
    double p1=p[1], p2=p[2], p3=p[3], p4=p[4], p5=p[5], p6=p[6];

    double r0=p2, r1=p3+0.5*p6, r2=p4-0.5*p6, r3=p5+0.5*p6, r4=p6, r5=p1;

    double s0, s1;
    double Graw = GH_pyx_solve_and_G(r0,r1,r2,r3,r4,r5,T,Rgas,&g,&s0,&s1);

    if (GH_nograd(grad, p, 7, 1.0)){
        double ends[7];
        GH_pyx_vertex_G(gvp, ends);
        ends[5] = GH_pyx_pure_ES_G(T,Rgas,Pv,gvp);
        double Emix = 0.0;
        for (int k = 0; k < 7; k++){ Emix += p[k]*ends[k]; }
        GH_sf_pyx(p, d->sf, opx);
        return GH_df_fast(d, p, gb, 7, Graw - Emix, 0.0);
    }

    double xal3m1,xfe2m1,xfe3m1,xmg2m1,xti4m1,xca2m2,xfe2m2,xmg2m2,xna1m2,xal3tet,xfe3tet,xsi4tet;
    GH_pyx_site_fracs(r0,r1,r2,r3,r4,r5,s0,s1,&xal3m1,&xfe2m1,&xfe3m1,&xmg2m1,&xti4m1,&xca2m2,&xfe2m2,&xmg2m2,&xna1m2,&xal3tet,&xfe3tet,&xsi4tet);
    double xm1_MgFeNoTi,xm1_notTi,xm1_MgFe,xtet_notSi,xrest_notNa,xm1_notMgFe,xm2_notNa;
    GH_pyx_composites(xmg2m1,xfe2m1,xti4m1,xna1m2,xsi4tet,&xm1_MgFeNoTi,&xm1_notTi,&xm1_MgFe,&xtet_notSi,&xrest_notNa,&xm1_notMgFe,&xm2_notNa);

    double dgdr0 = Rgas*T*(rlog((xfe2m1)) - rlog((xmg2m1))) +                g.cX2 +                g.cX2X2*r0*2.0 +                g.cX2X3*r1 +                g.cX2X4*r2 +                g.cX2X5*r3 +                g.cX2X6*r4 +                g.cX2X7*r5 +                g.cX2S1*s0 +                g.cX2S2*s1 +                g.cX2X2X7*r0*r5*2.0 +                g.cX2X2S2*r0*s1*2.0 +                g.cX2X3X7*r1*r5 +                g.cX2X3S2*r1*s1 +                g.cX2X4X7*r2*r5 +                g.cX2X4S2*r2*s1 +                g.cX2X5X7*r3*r5 +                g.cX2X5S2*r3*s1 +                g.cX2X6X7*r4*r5 +                g.cX2X6S2*r4*s1 +                g.cX2X7X7*r5*r5 +                g.cX2X7S2*r5*s1 +                g.cX2S2S2*s1*s1;
    double dgdr1 = Rgas*T*(0.5*rlog((xti4m1)) - 0.5*rlog((xmg2m1)) + 0.5*rlog((xm1_notTi))                  - rlog((xm1_MgFeNoTi)) + 0.5*rlog((xm1_MgFe))                  + 0.5*rlog((xrest_notNa)) - rlog((xtet_notSi))                  - 0.5*rlog((xm1_notMgFe)) + rlog((xal3tet))) +                g.cX3 +                g.cX2X3*r0 +                g.cX3X3*r1*2.0 +                g.cX3X4*r2 +                g.cX3X5*r3 +                g.cX3X6*r4 +                g.cX3X7*r5 +                g.cX3S1*s0 +                g.cX3S2*s1 +                g.cX2X3X7*r0*r5 +                g.cX2X3S2*r0*s1 +                g.cX3X3X7*r1*r5*2.0 +                g.cX3X3S2*r1*s1*2.0 +                g.cX3X4X7*r2*r5 +                g.cX3X4S2*r2*s1 +                g.cX3X5X7*r3*r5 +                g.cX3X5S2*r3*s1 +                g.cX3X6X7*r4*r5 +                g.cX3X6S2*r4*s1 +                g.cX3X7X7*r5*r5 +                g.cX3X7S2*r5*s1 +                g.cX3S2S2*s1*s1;
    double dgdr2 = Rgas*T*(0.5*rlog((xti4m1)) - 0.5*rlog((xmg2m1)) + 0.5*rlog((xm1_notTi))                  - rlog((xm1_MgFeNoTi)) + 0.5*rlog((xm1_MgFe))                  + 0.5*rlog((xrest_notNa)) - rlog((xtet_notSi))                  - 0.5*rlog((xm1_notMgFe)) + rlog((xfe3tet))) +                g.cX4 +                g.cX2X4*r0 +                g.cX3X4*r1 +                g.cX4X4*r2*2.0 +                g.cX4X5*r3 +                g.cX4X6*r4 +                g.cX4X7*r5 +                g.cX4S1*s0 +                g.cX4S2*s1 +                g.cX2X4X7*r0*r5 +                g.cX2X4S2*r0*s1 +                g.cX3X4X7*r1*r5 +                g.cX3X4S2*r1*s1 +                g.cX4X4X7*r2*r5*2.0 +                g.cX4X4S2*r2*s1*2.0 +                g.cX4X5X7*r3*r5 +                g.cX4X5S2*r3*s1 +                g.cX4X6X7*r4*r5 +                g.cX4X6S2*r4*s1 +                g.cX4X7X7*r5*r5 +                g.cX4X7S2*r5*s1 +                g.cX4S2S2*s1*s1;
    double dgdr3 = Rgas*T*(0.5*rlog((xal3m1)) + 0.5*rlog((xfe3m1)) - rlog((xmg2m1))                  - rlog((xm1_MgFeNoTi)) + rlog((xm1_MgFe))                  - rlog((xtet_notSi)) + 0.5*rlog((xal3tet)) + 0.5*rlog((xfe3tet))                  + rlog((xrest_notNa)) - rlog((xm1_notMgFe))) +                g.cX5 +                g.cX2X5*r0 +                g.cX3X5*r1 +                g.cX4X5*r2 +                g.cX5X5*r3*2.0 +                g.cX5X6*r4 +                g.cX5X7*r5 +                g.cX5S1*s0 +                g.cX5S2*s1 +                g.cX2X5X7*r0*r5 +                g.cX2X5S2*r0*s1 +                g.cX3X5X7*r1*r5 +                g.cX3X5S2*r1*s1 +                g.cX4X5X7*r2*r5 +                g.cX4X5S2*r2*s1 +                g.cX5X5X7*r3*r5*2.0 +                g.cX5X5S2*r3*s1*2.0 +                g.cX5X6X7*r4*r5 +                g.cX5X6S2*r4*s1 +                g.cX5X7X7*r5*r5 +                g.cX5X7S2*r5*s1 +                g.cX5S2S2*s1*s1;
    double dgdr4 = Rgas*T*(0.25*rlog((xal3m1)) + 0.25*rlog((xfe3m1)) - 0.5*rlog((xmg2m1))                  - 0.5*rlog((xm1_MgFeNoTi)) + 0.5*rlog((xm1_MgFe))                  + 0.5*rlog((xtet_notSi)) - 0.25*rlog((xal3tet))                  - 0.25*rlog((xfe3tet)) - rlog((xca2m2)) + rlog((xna1m2))                  - 0.5*rlog((xrest_notNa))                  - 0.5*rlog((xm1_notMgFe)) + rlog((xm2_notNa))) +                g.cX6 +                g.cX2X6*r0 +                g.cX3X6*r1 +                g.cX4X6*r2 +                g.cX5X6*r3 +                g.cX6X6*r4*2.0 +                g.cX6X7*r5 +                g.cX6S1*s0 +                g.cX6S2*s1 +                g.cX2X6X7*r0*r5 +                g.cX2X6S2*r0*s1 +                g.cX3X6X7*r1*r5 +                g.cX3X6S2*r1*s1 +                g.cX4X6X7*r2*r5 +                g.cX4X6S2*r2*s1 +                g.cX5X6X7*r3*r5 +                g.cX5X6S2*r3*s1 +                g.cX6X6X7*r4*r5*2.0 +                g.cX6X6S2*r4*s1*2.0 +                g.cX6X7X7*r5*r5 +                g.cX6X7S2*r5*s1 +                g.cX6S2S2*s1*s1;
    double dgdr5 = 0.5*Rgas*T*(rlog((xmg2m1)) - 2.0*rlog((xca2m2))                  + rlog((xmg2m2)) + rlog((xfe2m2/xfe2m1))) +                g.cX7 +                g.cX2X7*r0 +                g.cX3X7*r1 +                g.cX4X7*r2 +                g.cX5X7*r3 +                g.cX6X7*r4 +                g.cX7X7*r5*2.0 +                g.cX7S1*s0 +                g.cX7S2*s1 +                g.cX2X2X7*r0*r0 +                g.cX2X3X7*r0*r1 +                g.cX2X4X7*r0*r2 +                g.cX2X5X7*r0*r3 +                g.cX2X6X7*r0*r4 +                g.cX2X7X7*r0*r5*2.0 +                g.cX2X7S2*r0*s1 +                g.cX3X3X7*r1*r1 +                g.cX3X4X7*r1*r2 +                g.cX3X5X7*r1*r3 +                g.cX3X6X7*r1*r4 +                g.cX3X7X7*r1*r5*2.0 +                g.cX3X7S2*r1*s1 +                g.cX4X4X7*r2*r2 +                g.cX4X5X7*r2*r3 +                g.cX4X6X7*r2*r4 +                g.cX4X7X7*r2*r5*2.0 +                g.cX4X7S2*r2*s1 +                g.cX5X5X7*r3*r3 +                g.cX5X6X7*r3*r4 +                g.cX5X7X7*r3*r5*2.0 +                g.cX5X7S2*r3*s1 +                g.cX6X6X7*r4*r4 +                g.cX6X7X7*r4*r5*2.0 +                g.cX6X7S2*r4*s1 +                g.cX7X7X7*r5*r5*3.0 +                g.cX7X7S2*r5*s1*2.0 +                g.cX7S2S2*s1*s1;

    double es_end = GH_pyx_pure_ES_G(T,Rgas,Pv,gvp);
    double ends[7];
    GH_pyx_vertex_G(gvp, ends);
    ends[5] = es_end;

    double Emix = 0.0;
    for (int k = 0; k < 7; k++){ Emix += p[k]*ends[k]; }
    double Gex = Graw - Emix;

    {
        double gs0, gs1;
        GH_pyx_dgds(r0,r1,r2,r3,r4,r5,s0,s1,T,Rgas,&g,&gs0,&gs1);
        double rr[6] = {r0,r1,r2,r3,r4,r5};
        double dGdr[6] = {dgdr0,dgdr1,dgdr2,dgdr3,dgdr4,dgdr5};
        GH_pyx_bound_chain(rr, s0, s1, gs0, gs1, 1.0e-8, dGdr);
        dgdr0 = dGdr[0]; dgdr1 = dGdr[1]; dgdr2 = dGdr[2]; dgdr3 = dGdr[3]; dgdr4 = dGdr[4]; dgdr5 = dGdr[5];
    }
    double dGex[7];
    dGex[0] = 0.0;
    dGex[1] = dgdr5;
    dGex[2] = dgdr0;
    dGex[3] = dgdr1;
    dGex[4] = dgdr2;
    dGex[5] = dgdr3;
    dGex[6] = 0.5*dgdr1 - 0.5*dgdr2 + 0.5*dgdr3 + dgdr4;
    for (int k = 1; k < 7; k++){ dGex[k] -= ends[k] - ends[0]; }

    GH_mu_from_dGex(d, p, 1.0, 7, Gex, dGex);

    GH_sf_pyx(p, d->sf, opx);
    return GH_df_assemble(d, p, gb, mu_Gex, 7, grad);
}

double obj_gh_cpx(unsigned n, const double *x, double *grad, void *SS_ref_db){
    return GH_pyx_obj(n, x, grad, SS_ref_db, &GH_PYX_CPX, -1.08018328, 0);
}

double obj_gh_opx(unsigned n, const double *x, double *grad, void *SS_ref_db){
    return GH_pyx_obj(n, x, grad, SS_ref_db, &GH_PYX_OPX, -0.57977688, 1);
}

typedef struct {
    double g11p, g22p, g33p, w11p, w12, w12p, w13, w13p, w1p2, w1p2p, w1p3, w1p3p, w22p, w23, w23p, w2p3, w2p3p, w33p;
} GH_oox_par;

static const GH_oox_par GH_OOX = {
    .g11p = 22000.0,
    .g22p = 22000.0,
    .g33p = 22000.0,
    .w11p = 6000.0,
    .w12 = 26000.0,
    .w12p = 26000.0,
    .w13 = 45000.0,
    .w13p = 45000.0,
    .w1p2 = 26000.0,
    .w1p2p = 26000.0,
    .w1p3 = 45000.0,
    .w1p3p = 45000.0,
    .w22p = 6000.0,
    .w23 = 8400.0,
    .w23p = 8400.0,
    .w2p3 = 8400.0,
    .w2p3p = 8400.0,
    .w33p = 6000.0,
};

typedef struct {
    double g0, gx2, gx3, gs1, gs2, gs3, gx2x2, gx2x3, gx2s1, gx2s2, gx2s3, gx3x3, gx3s1, gx3s2, gx3s3, gs1s1, gs1s2, gs1s3, gs2s2, gs2s3, gs3s3;
} GH_oox_coef;

static void GH_oox_coef_make(const GH_oox_par *P, GH_oox_coef *c){
    c->g0 = P->g11p/2.0 + P->w11p/4.0;
    c->gx2 = (P->g22p - P->g11p - P->w11p)/2.0 + (P->w12 + P->w12p + P->w1p2 + P->w1p2p)/4.0;
    c->gx3 = (P->g33p - P->g11p - P->w11p)/2.0 + (P->w13 + P->w13p + P->w1p3 + P->w1p3p)/4.0;
    c->gs1 = -P->g11p/2.0;
    c->gs2 = -P->g22p/2.0 + (P->w12 - P->w12p + P->w1p2 - P->w1p2p)/4.0;
    c->gs3 = -P->g33p/2.0 + (P->w13 - P->w13p + P->w1p3 - P->w1p3p)/4.0;
    c->gx2x2 = (P->w11p - P->w12 - P->w12p - P->w1p2 - P->w1p2p + P->w22p)/4.0;
    c->gx2x3 = P->w11p/2.0
            + (-P->w12 - P->w12p - P->w13 - P->w13p - P->w1p2 - P->w1p2p - P->w1p3 - P->w1p3p + P->w23 + P->w23p + P->w2p3 + P->w2p3p)/4.0;
    c->gx2s1 = (P->w12 + P->w12p - P->w1p2 - P->w1p2p)/4.0;
    c->gx2s2 = (-P->w12 + P->w12p - P->w1p2 + P->w1p2p)/4.0;
    c->gx2s3 = (-P->w13 + P->w13p - P->w1p3 + P->w1p3p + P->w23 - P->w23p + P->w2p3 - P->w2p3p)/4.0;
    c->gx3x3 = (P->w11p - P->w13 - P->w13p - P->w1p3 - P->w1p3p + P->w33p)/4.0;
    c->gx3s1 = (P->w13 + P->w13p - P->w1p3 - P->w1p3p)/4.0;
    c->gx3s2 = (-P->w12 + P->w12p - P->w1p2 + P->w1p2p + P->w23 + P->w23p - P->w2p3 - P->w2p3p)/4.0;
    c->gx3s3 = (-P->w13 + P->w13p - P->w1p3 + P->w1p3p)/4.0;
    c->gs1s1 = -P->w11p/4.0;
    c->gs1s2 = (P->w12 - P->w12p - P->w1p2 + P->w1p2p)/4.0;
    c->gs1s3 = (P->w13 - P->w13p - P->w1p3 + P->w1p3p)/4.0;
    c->gs2s2 = -P->w22p/4.0;
    c->gs2s3 = (P->w23 - P->w23p - P->w2p3 + P->w2p3p)/4.0;
    c->gs3s3 = -P->w33p/4.0;
}

static void GH_oox_site_fracs(double r0, double r1, const double s[3], double x[8]){
    x[0] = (r0 + s[1])/2.0;
    x[1] = (r1 + s[2])/2.0;
    x[2] = (1.0 - s[0] - s[1] - s[2])/2.0;
    x[3] = (1.0 - r0 - r1 + s[0])/2.0;
    x[4] = (r0 - s[1])/4.0;
    x[5] = (r1 - s[2])/4.0;
    x[6] = (1.0 + 2.0*r0 + 2.0*r1 + s[0] + s[1] + s[2])/4.0;
    x[7] = (3.0 - 3.0*r0 - 3.0*r1 - s[0])/4.0;
    for (int i = 0; i < 8; i++){ if (x[i] <= 0.0){ x[i] = DBL_EPSILON; } }
}

static void GH_oox_order(const GH_oox_coef *c, double T, double Rgas, double r0, double r1, double s[3], double x[8]){
    double sOld[3] = {2.0, 2.0, 2.0};
    double sNew[3] = {0.9*(1.0-r0-r1), 0.9*r0, 0.9*r1};
    int iter = 0;
    while ((fabs(sNew[0]-sOld[0]) > 10.0*DBL_EPSILON || fabs(sNew[1]-sOld[1]) > 10.0*DBL_EPSILON ||
            fabs(sNew[2]-sOld[2]) > 10.0*DBL_EPSILON) && iter < 100){
        double sc[3] = {sNew[0], sNew[1], sNew[2]};
        GH_oox_site_fracs(r0, r1, sc, x);
        double fe2M1 = x[0], mg2M1 = x[1], ti4M1 = x[2], fe3M1 = x[3];
        double fe2M2 = x[4], mg2M2 = x[5], ti4M2 = x[6], fe3M2 = x[7];
        double RT = Rgas*T;
        double g[3], H[3][3];
        g[0] = RT*(-0.5*log(ti4M1) + 0.5*log(fe3M1) + 0.5*log(ti4M2) - 0.5*log(fe3M2))
             + c->gs1 + c->gx2s1*r0 + c->gx3s1*r1 + 2.0*c->gs1s1*sc[0] + c->gs1s2*sc[1] + c->gs1s3*sc[2];
        g[1] = RT*( 0.5*log(fe2M1) - 0.5*log(ti4M1) - 0.5*log(fe2M2) + 0.5*log(ti4M2))
             + c->gs2 + c->gx2s2*r0 + c->gx3s2*r1 + c->gs1s2*sc[0] + 2.0*c->gs2s2*sc[1] + c->gs2s3*sc[2];
        g[2] = RT*( 0.5*log(mg2M1) - 0.5*log(ti4M1) - 0.5*log(mg2M2) + 0.5*log(ti4M2))
             + c->gs3 + c->gx2s3*r0 + c->gx3s3*r1 + c->gs1s3*sc[0] + c->gs2s3*sc[1] + 2.0*c->gs3s3*sc[2];
        H[0][0] = RT*(0.25/ti4M1 + 0.25/fe3M1 + 0.125/ti4M2 + 0.125/fe3M2) + 2.0*c->gs1s1;
        H[0][1] = RT*(0.25/ti4M1 + 0.125/ti4M2) + c->gs1s2;
        H[0][2] = RT*(0.25/ti4M1 + 0.125/ti4M2) + c->gs1s3;
        H[1][1] = RT*(0.25/fe2M1 + 0.25/ti4M1 + 0.125/fe2M2 + 0.125/ti4M2) + 2.0*c->gs2s2;
        H[1][2] = RT*(0.25/ti4M1 + 0.125/ti4M2) + c->gs2s3;
        H[2][2] = RT*(0.25/mg2M1 + 0.25/ti4M1 + 0.125/mg2M2 + 0.125/ti4M2) + 2.0*c->gs3s3;
        H[1][0] = H[0][1]; H[2][0] = H[0][2]; H[2][1] = H[1][2];
        double inv[3][3];
        GH_invert3(H, inv);
        for (int i = 0; i < 3; i++){ sOld[i] = sc[i]; }
        for (int i = 0; i < 3; i++){
            sNew[i] = sc[i];
            for (int j = 0; j < 3; j++){ sNew[i] -= inv[i][j]*g[j]; }
        }
        iter++;
    }
    for (int i = 0; i < 3; i++){ s[i] = sOld[i]; }
    GH_oox_site_fracs(r0, r1, s, x);
}

static void GH_oox_pure_ends(const GH_oox_coef *c, double T, double Rgas, double ends[3]){
    double sOld[3] = {2.0, 2.0, 2.0}, sNew[3] = {0.9, 0.9, 0.9};
    int iter = 0;
    while ((fabs(sNew[0]-sOld[0]) > 10.0*DBL_EPSILON || fabs(sNew[1]-sOld[1]) > 10.0*DBL_EPSILON ||
            fabs(sNew[2]-sOld[2]) > 10.0*DBL_EPSILON) && iter < 100){
        double sc[3] = {sNew[0], sNew[1], sNew[2]};
        double RT = Rgas*T, g[3], h[3];
        g[0] = RT*(log((1.0+sc[0])/2.0)/2.0 - log((1.0-sc[0])/2.0)/2.0 - log((3.0-sc[0])/4.0)/2.0 + log((1.0+sc[0])/4.0)/2.0)
             + c->gs1 + 2.0*c->gs1s1*sc[0];
        g[1] = RT*(log((1.0+sc[1])/2.0)/2.0 - log((1.0-sc[1])/4.0)/2.0) + c->gs2 + c->gx2s2 + 2.0*c->gs2s2*sc[1];
        g[2] = RT*(log((1.0+sc[2])/2.0)/2.0 - log((1.0-sc[2])/4.0)/2.0) + c->gs3 + c->gx3s3 + 2.0*c->gs3s3*sc[2];
        h[0] = RT*(0.5/(1.0+sc[0]) + 0.5/(1.0-sc[0]) + 0.5/(3.0-sc[0]) + 0.5/(1.0+sc[0])) + 2.0*c->gs1s1;
        h[1] = RT*(0.5/(1.0+sc[1]) + 0.5/(1.0-sc[1])) + 2.0*c->gs2s2;
        h[2] = RT*(0.5/(1.0+sc[2]) + 0.5/(1.0-sc[2])) + 2.0*c->gs3s3;
        for (int i = 0; i < 3; i++){ sOld[i] = sc[i]; }
        for (int i = 0; i < 3; i++){
            double v = sc[i] - g[i]/h[i];
            if (v >  1.0 - DBL_EPSILON){ v =  1.0 - DBL_EPSILON; }
            if (v < -1.0 + DBL_EPSILON){ v = -1.0 + DBL_EPSILON; }
            sNew[i] = v;
        }
        iter++;
    }
    double s0 = sOld[0], s1 = sOld[1], s2 = sOld[2];
    double S_pb = -Rgas*(((1.0+s0)/2.0)*log((1.0+s0)/2.0) + ((1.0-s0)/2.0)*log((1.0-s0)/2.0)
                       + ((3.0-s0)/2.0)*log((3.0-s0)/4.0) + ((1.0+s0)/2.0)*log((1.0+s0)/4.0));
    double S_fe = -Rgas*(((1.0+s1)/2.0)*log((1.0+s1)/2.0) + ((1.0-s1)/2.0)*log((1.0-s1)/4.0));
    double S_mg = -Rgas*(((1.0+s2)/2.0)*log((1.0+s2)/2.0) + ((1.0-s2)/2.0)*log((1.0-s2)/4.0));
    ends[0] = c->g0 + c->gs1*s0 + c->gs1s1*s0*s0 - T*S_pb;
    ends[1] = c->gx2 + c->gs2*s1 + c->gx2x2 + c->gx2s2*s1 + c->gs2s2*s1*s1 - T*S_fe;
    ends[2] = c->gx3 + c->gs3*s2 + c->gx3x3 + c->gx3s3*s2 + c->gs3s3*s2*s2 - T*S_mg;
}

double obj_gh_oox(unsigned n, const double *x, double *grad, void *SS_ref_db){
    SS_ref *d = (SS_ref *) SS_ref_db;

    double T    = d->T;
    double Rgas = GH_R_J;
    double *p   = d->p;
    double *gb  = d->gb_lvl;
    double *mu_Gex = d->mu_Gex;

    for (int i = 0; i < 3; i++){ p[i] = x[i]; }
    double sum = p[0] + p[1] + p[2];
    double r0 = p[1]/sum, r1 = p[2]/sum;

    GH_oox_coef c;
    GH_oox_coef_make(&GH_OOX, &c);
    double s[3], xs[8];
    GH_oox_order(&c, T, Rgas, r0, r1, s, xs);
    double fe2M1 = xs[0], mg2M1 = xs[1], ti4M1 = xs[2], fe3M1 = xs[3];
    double fe2M2 = xs[4], mg2M2 = xs[5], ti4M2 = xs[6], fe3M2 = xs[7];

    double S = -Rgas*(fe2M1*log(fe2M1) + mg2M1*log(mg2M1) + ti4M1*log(ti4M1) + fe3M1*log(fe3M1)
                    + 2.0*fe2M2*log(fe2M2) + 2.0*mg2M2*log(mg2M2) + 2.0*ti4M2*log(ti4M2) + 2.0*fe3M2*log(fe3M2));
    double H = c.g0 + c.gx2*r0 + c.gx3*r1 + c.gs1*s[0] + c.gs2*s[1] + c.gs3*s[2]
             + c.gx2x2*r0*r0 + c.gx2x3*r0*r1 + c.gx2s1*r0*s[0]
             + c.gx2s2*r0*s[1] + c.gx2s3*r0*s[2] + c.gx3x3*r1*r1
             + c.gx3s1*r1*s[0] + c.gx3s2*r1*s[1] + c.gx3s3*r1*s[2]
             + c.gs1s1*s[0]*s[0] + c.gs1s2*s[0]*s[1] + c.gs1s3*s[0]*s[2]
             + c.gs2s2*s[1]*s[1] + c.gs2s3*s[1]*s[2] + c.gs3s3*s[2]*s[2];

    double ends[3];
    GH_oox_pure_ends(&c, T, Rgas, ends);
    double Gex = H - T*S - ((1.0-r0-r1)*ends[0] + r0*ends[1] + r1*ends[2]);

    double dGdr0 = Rgas*T*(0.5*log(fe2M1) - 0.5*log(fe3M1) + 0.5*log(fe2M2) + log(ti4M2) - 1.5*log(fe3M2))
                 + c.gx2 + 2.0*c.gx2x2*r0 + c.gx2x3*r1 + c.gx2s1*s[0] + c.gx2s2*s[1] + c.gx2s3*s[2];
    double dGdr1 = Rgas*T*(0.5*log(mg2M1) - 0.5*log(fe3M1) + 0.5*log(mg2M2) + log(ti4M2) - 1.5*log(fe3M2))
                 + c.gx3 + c.gx2x3*r0 + 2.0*c.gx3x3*r1 + c.gx3s1*s[0] + c.gx3s2*s[1] + c.gx3s3*s[2];

    double dGex[3];
    dGex[0] = 0.0;
    dGex[1] = dGdr0 - (ends[1] - ends[0]);
    dGex[2] = dGdr1 - (ends[2] - ends[0]);

    GH_mu_from_dGex(d, p, sum, 3, Gex, dGex);

    return GH_df_assemble(d, p, gb, mu_Gex, 3, grad);
}

typedef struct {
    double gx123_4, gx13_24, gx134_2, w123, w4, dz, gxca_13, gxca_2, wcamg, wcafe, dwcamg, dwcafe;
} GH_amp_par;

static const GH_amp_par GH_AMP[2] = {
    [0] = { .gx123_4 = 13.7354*1000.0, .gx13_24 = 5.25246*1000.0, .gx134_2 = 16.7461*1000.0, .w123 = 5.52380*1000.0, .w4 = 3.72280*1000.0, .dz = -0.210304*1000.0, .gxca_13 = 17.801*1000.0, .gxca_2 = -11.395*1000.0, .wcamg = 30.100*1000.0, .wcafe = 19.500*1000.0, .dwcamg = 0.0000*1000.0, .dwcafe = 0.0000*1000.0 },
    [1] = { .gx123_4 = 12.1027*1000.0, .gx13_24 = 3.29327*1000.0, .gx134_2 = 15.1135*1000.0, .w123 = 5.11564*1000.0, .w4 = 3.55954*1000.0, .dz = -0.781732*1000.0, .gxca_13 = 17.801*1000.0, .gxca_2 = -11.395*1000.0, .wcamg = 28.100*1000.0, .wcafe = 19.500*1000.0, .dwcamg = 0.0000*1000.0, .dwcafe = 0.0000*1000.0 },
};

static const double GH_AMP_GEX  = -25.8556*1000.0;
static const double GH_AMP_GDEX = 5.90595*1000.0;

static const GH_hsv GH_AMP_D[3] = {
    { (-12073132.08) - (-12067920.38), (535.2587) - (540.2587), (26.310) - (26.340) },
    { (-9627014.85) - (-9623550.0),    (720.0) - (725.0),       (27.810) - (27.840) },
    { (-12307863.0) - (-12307863.0),   (544.500) - (549.500),   (27.312) - (27.312) }
};

typedef struct {
    double g0, gr, gs13, gs2, grr, grs13, grs2, gs13s13, gs13s2, gs2s2, gp, grp, gps13, gps2, gpp, gps13s13, gpps13, gpps2, grps13, grps2, gps2s2, gps13s2, grrp, grpp, gppp;
} GH_amp_coef;

typedef struct { double xfe2m4, xmg2m4, xca2m4, xfe2m13, xmg2m13, xfe2m2, xmg2m2; } GH_amp_sf;

static void GH_amp_coef_make(int clino, double t, double p, GH_amp_coef *c){
    const GH_amp_par *P = &GH_AMP[clino];
    double g1 = clino ? 0.0 : GH_hsv_G(GH_AMP_D[0], t, p - 1.0);
    double g2 = clino ? 0.0 : GH_hsv_G(GH_AMP_D[1], t, p - 1.0);
    double g3 = clino ? 0.0 : GH_hsv_G(GH_AMP_D[2], t, p - 1.0);

    c->g0 = g1/2.0 + g2/2.0 + 5.0*P->w123/4.0 + P->w4/2.0 + P->gx123_4/4.0;
    c->gr = -g1/2.0 + g2/2.0;
    c->gs13 = 3.0*GH_AMP_GEX/5.0 + GH_AMP_GDEX;
    c->gs2 = 2.0*GH_AMP_GEX/5.0 - GH_AMP_GDEX;
    c->grr = -5.0*P->w123/4.0 - P->w4/2.0 - P->gx123_4/4.0;
    c->grs13 = 5.0*P->w123/14.0 - 6.0*P->w4/7.0 - 5.0*P->gx123_4/28.0 + P->dz/2.0;
    c->grs2 = 15.0*P->w123/14.0 - 4.0*P->w4/7.0 - P->gx123_4/28.0 - P->dz/2.0;
    c->gs13s13 = -125.0*P->w123/98.0 - 18.0*P->w4/49.0 - 43.0*P->gx123_4/196.0 + P->gx13_24/2.0 - P->dz/14.0;
    c->gs13s2 = 115.0*P->w123/49.0 - 24.0*P->w4/49.0 + 30.0*P->gx123_4/49.0 - P->gx13_24/2.0 - P->gx134_2/2.0
            - P->dz/7.0;
    c->gs2s2 = -145.0*P->w123/98.0 - 8.0*P->w4/49.0 - 37.0*P->gx123_4/196.0 + P->gx134_2/2.0 + 3.0*P->dz/14.0;
    c->gp = -g1 + g3 - P->w4 + P->wcafe + P->wcamg + P->gxca_13/2.0 + P->gxca_2 - P->dwcafe/2.0
            - P->dwcamg/2.0;
    c->grp = -P->w4 + P->wcafe - P->wcamg + P->gxca_13/2.0 + P->gxca_2 - P->dwcafe + P->dwcamg;
    c->gps13 = -6.0*P->w4/7.0 + 6.0*P->wcafe/7.0 - 6.0*P->wcamg/7.0 - 4.0*P->gxca_13/7.0 - P->gxca_2/7.0
            - 6.0*P->dwcafe/7.0 + 6.0*P->dwcamg/7.0;
    c->gps2 = -4.0*P->w4/7.0 + 4.0*P->wcafe/7.0 - 4.0*P->wcamg/7.0 + 2.0*P->gxca_13/7.0 - 3.0*P->gxca_2/7.0
            - 4.0*P->dwcafe/7.0 + 4.0*P->dwcamg/7.0;
    c->gpp = -2.0*P->wcamg + P->dwcafe + 3.0*P->dwcamg;
    c->gps13s13 = -18.0*P->dwcafe/49.0 - 18.0*P->dwcamg/49.0;
    c->gpps13 = 6.0*P->dwcafe/7.0 - 18.0*P->dwcamg/7.0;
    c->gpps2 = 4.0*P->dwcafe/7.0 - 12.0*P->dwcamg/7.0;
    c->grps13 = -6.0*P->dwcafe/7.0 - 6.0*P->dwcamg/7.0;
    c->grps2 = -4.0*P->dwcafe/7.0 - 4.0*P->dwcamg/7.0;
    c->gps2s2 = -8.0*P->dwcafe/49.0 - 8.0*P->dwcamg/49.0;
    c->gps13s2 = -24.0*P->dwcafe/49.0 - 24.0*P->dwcamg/49.0;
    c->grrp = -P->dwcafe/2.0 - P->dwcamg/2.0;
    c->grpp = P->dwcafe - 3.0*P->dwcamg;
    c->gppp = -4.0*P->dwcamg;
}

static void GH_amp_site_fracs(const double r[2], const double s[2], GH_amp_sf *f){
    f->xca2m4  = r[1];
    f->xfe2m4  = (1.0 + r[0] + 6.0*s[0]/7.0 +  4.0*s[1]/7.0)/2.0;
    f->xfe2m13 = (1.0 + r[0] - 8.0*s[0]/7.0 +  4.0*s[1]/7.0)/2.0;
    f->xfe2m2  = (1.0 + r[0] + 6.0*s[0]/7.0 - 10.0*s[1]/7.0)/2.0;
    f->xmg2m4  = 1.0 - f->xfe2m4  - f->xca2m4;
    f->xmg2m13 = 1.0 - f->xfe2m13;
    f->xmg2m2  = 1.0 - f->xfe2m2;
}

static void GH_amp_clamp_sf(GH_amp_sf *f){
    if (f->xfe2m4  <= DBL_EPSILON) f->xfe2m4  = DBL_EPSILON;
    if (f->xfe2m13 <= DBL_EPSILON) f->xfe2m13 = DBL_EPSILON;
    if (f->xfe2m2  <= DBL_EPSILON) f->xfe2m2  = DBL_EPSILON;
    if (f->xmg2m4  <= DBL_EPSILON) f->xmg2m4  = DBL_EPSILON;
    if (f->xmg2m13 <= DBL_EPSILON) f->xmg2m13 = DBL_EPSILON;
    if (f->xmg2m2  <= DBL_EPSILON) f->xmg2m2  = DBL_EPSILON;
}

static int GH_amp_sf_out_of_range(const GH_amp_sf *f){
    return f->xfe2m4  < 0.0 || f->xfe2m4  > 1.0 || f->xmg2m4  < 0.0 || f->xmg2m4  > 1.0
        || f->xfe2m13 < 0.0 || f->xfe2m13 > 1.0 || f->xmg2m13 < 0.0 || f->xmg2m13 > 1.0
        || f->xfe2m2  < 0.0 || f->xfe2m2  > 1.0 || f->xmg2m2  < 0.0 || f->xmg2m2  > 1.0;
}

static double GH_amp_tiny(double x){ return (x > 0.0) ? x : 1.0e-300; }

static GH_amp_sf GH_amp_floor_sf(const GH_amp_sf *f){
    GH_amp_sf g;
    g.xfe2m4  = GH_amp_tiny(f->xfe2m4);  g.xmg2m4  = GH_amp_tiny(f->xmg2m4);  g.xca2m4 = GH_amp_tiny(f->xca2m4);
    g.xfe2m13 = GH_amp_tiny(f->xfe2m13); g.xmg2m13 = GH_amp_tiny(f->xmg2m13);
    g.xfe2m2  = GH_amp_tiny(f->xfe2m2);  g.xmg2m2  = GH_amp_tiny(f->xmg2m2);
    return g;
}

static double GH_amp_entropy_ideal(const GH_amp_sf *f){
    return -GH_R_J*( 2.0*f->xfe2m4*log(f->xfe2m4) + 2.0*f->xmg2m4*log(f->xmg2m4) + 2.0*f->xca2m4*log(f->xca2m4)
                     + 3.0*f->xfe2m13*log(f->xfe2m13) + 3.0*f->xmg2m13*log(f->xmg2m13)
                     + 2.0*f->xfe2m2*log(f->xfe2m2) + 2.0*f->xmg2m2*log(f->xmg2m2) );
}

static double GH_amp_gibbs(const GH_amp_coef *c, const double r[2], const double s[2], const GH_amp_sf *f, double t){
    return -t*GH_amp_entropy_ideal(f) + c->g0 + c->gr*r[0] + c->gs13*s[0] + c->gs2*s[1] + c->grr*r[0]*r[0]
            + c->grs13*r[0]*s[0] + c->grs2*r[0]*s[1] + c->gs13s13*s[0]*s[0] + c->gs13s2*s[0]*s[1]
            + c->gs2s2*s[1]*s[1] + c->gp*r[1] + c->grp*r[0]*r[1] + c->gps13*r[1]*s[0] + c->gps2*r[1]*s[1]
            + c->gpp*r[1]*r[1] + c->gps13s13*r[1]*s[0]*s[0] + c->gpps13*r[1]*r[1]*s[0]
            + c->gpps2*r[1]*r[1]*s[1] + c->grps13*r[0]*r[1]*s[0] + c->grps2*r[0]*r[1]*s[1]
            + c->gps2s2*r[1]*s[1]*s[1] + c->gps13s2*r[1]*s[0]*s[1] + c->grrp*r[0]*r[0]*r[1]
            + c->grpp*r[0]*r[1]*r[1] + c->gppp*r[1]*r[1]*r[1];
}

static void GH_amp_dg_dr(const GH_amp_coef *c, const double r[2], const double s[2], const GH_amp_sf *f, double t, double dr[2]){
    dr[0] = GH_R_J*t*(log(f->xfe2m4) - log(f->xmg2m4) + 3.0*log(f->xfe2m13)/2.0 - 3.0*log(f->xmg2m13)/2.0 + log(f->xfe2m2) - log(f->xmg2m2))
            + c->gr + 2.0*c->grr*r[0] + c->grs13*s[0] + c->grs2*s[1] + c->grp*r[1] + c->grps13*r[1]*s[0]
            + c->grps2*r[1]*s[1] + 2.0*c->grrp*r[0]*r[1] + c->grpp*r[1]*r[1];
    dr[1] = GH_R_J*t*(-2.0*log(f->xmg2m4) + 2.0*log(f->xca2m4)) + c->gp + c->grp*r[0] + c->gps13*s[0]
            + c->gps2*s[1] + 2.0*c->gpp*r[1] + c->gps13s13*s[0]*s[0] + 2.0*c->gpps13*r[1]*s[0]
            + 2.0*c->gpps2*r[1]*s[1] + c->grps13*r[0]*s[0] + c->grps2*r[0]*s[1] + c->gps2s2*s[1]*s[1]
            + c->gps13s2*s[0]*s[1] + c->grrp*r[0]*r[0] + 2.0*c->grpp*r[0]*r[1] + 3.0*c->gppp*r[1]*r[1];
}

static void GH_amp_dg_ds(const GH_amp_coef *c, const double r[2], const double s[2], const GH_amp_sf *f, double t, double dg[2]){
    dg[0] = GH_R_J*t*(6.0*log(f->xfe2m4)/7.0 - 6.0*log(f->xmg2m4)/7.0 - 12.0*log(f->xfe2m13)/7.0 + 12.0*log(f->xmg2m13)/7.0 + 6.0*log(f->xfe2m2)/7.0 - 6.0*log(f->xmg2m2)/7.0)
            + c->gs13 + c->grs13*r[0] + 2.0*c->gs13s13*s[0] + c->gs13s2*s[1] + c->gps13*r[1]
            + 2.0*c->gps13s13*r[1]*s[0] + c->gpps13*r[1]*r[1] + c->grps13*r[0]*r[1] + c->gps13s2*r[1]*s[1];
    dg[1] = GH_R_J*t*(4.0*log(f->xfe2m4)/7.0 - 4.0*log(f->xmg2m4)/7.0 + 6.0*log(f->xfe2m13)/7.0 - 6.0*log(f->xmg2m13)/7.0 - 10.0*log(f->xfe2m2)/7.0 + 10.0*log(f->xmg2m2)/7.0)
            + c->gs2 + c->grs2*r[0] + c->gs13s2*s[0] + 2.0*c->gs2s2*s[1] + c->gps2*r[1] + c->gpps2*r[1]*r[1]
            + c->grps2*r[0]*r[1] + 2.0*c->gps2s2*r[1]*s[1] + c->gps13s2*r[1]*s[0];
}

static void GH_amp_d2g_ds2(const GH_amp_coef *c, const double r[2], const GH_amp_sf *f, double t, double H[2][2]){
    H[0][0] = GH_R_J*t*(6.0*3.0/(7.0*7.0*f->xfe2m4) + 6.0*3.0/(7.0*7.0*f->xmg2m4) + 12.0*4.0/(7.0*7.0*f->xfe2m13) + 12.0*4.0/(7.0*7.0*f->xmg2m13) + 6.0*3.0/(7.0*7.0*f->xfe2m2) + 6.0*3.0/(7.0*7.0*f->xmg2m2))
            + 2.0*c->gs13s13 + 2.0*c->gps13s13*r[1];
    H[0][1] = GH_R_J*t*(6.0*2.0/(7.0*7.0*f->xfe2m4) + 6.0*2.0/(7.0*7.0*f->xmg2m4) - 12.0*2.0/(7.0*7.0*f->xfe2m13) - 12.0*2.0/(7.0*7.0*f->xmg2m13) - 6.0*5.0/(7.0*7.0*f->xfe2m2) - 6.0*5.0/(7.0*7.0*f->xmg2m2))
            + c->gs13s2 + c->gps13s2*r[1];
    H[1][0] = H[0][1];
    H[1][1] = GH_R_J*t*(4.0*2.0/(7.0*7.0*f->xfe2m4) + 4.0*2.0/(7.0*7.0*f->xmg2m4) + 6.0*2.0/(7.0*7.0*f->xfe2m13) + 6.0*2.0/(7.0*7.0*f->xmg2m13) + 10.0*5.0/(7.0*7.0*f->xfe2m2) + 10.0*5.0/(7.0*7.0*f->xmg2m2))
            + 2.0*c->gs2s2 + 2.0*c->gps2s2*r[1];
}

static void GH_amp_order(const GH_amp_coef *c, double t, const double r[2], double s_out[2], GH_amp_sf *sf){
    GH_amp_sf f = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    double sOld[2] = {2.0, 2.0}, sNew[2];
    double xfe, xmg;
    int iter = 0;

    f.xca2m4 = r[1];
    xfe      = (1.0+r[0])/2.0;
    xmg      = (7.0-2.0*f.xca2m4-7.0*xfe)/7.0;
    sNew[0]  = ((1.0-f.xca2m4)*xfe - xfe)/(xfe+xmg);
    sNew[1]  = ((1.0-f.xca2m4)*xfe - xfe)/(xfe+xmg);

    if (xfe < 10.0*DBL_EPSILON){
        f.xfe2m4  = DBL_EPSILON;
        f.xfe2m13 = DBL_EPSILON;
        f.xfe2m2  = DBL_EPSILON;
        f.xmg2m4  = 1.0 - f.xca2m4;
        f.xmg2m13 = 1.0;
        f.xmg2m2  = 1.0;
        sOld[0] = sNew[0];
        sOld[1] = sNew[1];
    }

    while ((fabs(sNew[0]-sOld[0]) > 10.0*DBL_EPSILON || fabs(sNew[1]-sOld[1]) > 10.0*DBL_EPSILON) && iter < 200){
        double s[2], sCorr[2], lambda = 1.0, dgds[2], d2gds2[2][2];

        s[0] = sNew[0]; s[1] = sNew[1];

        GH_amp_site_fracs(r, s, &f);
        GH_amp_clamp_sf(&f);
        GH_amp_dg_ds(c, r, s, &f, t, dgds);
        GH_amp_d2g_ds2(c, r, &f, t, d2gds2);

        sOld[0] = s[0]; sOld[1] = s[1];

        double det = d2gds2[0][0]*d2gds2[1][1] - d2gds2[0][1]*d2gds2[1][0];
        double i00 =  d2gds2[1][1]/det, i01 = -d2gds2[0][1]/det;
        double i10 = -d2gds2[1][0]/det, i11 =  d2gds2[0][0]/det;
        sCorr[0] = -(i00*dgds[0] + i01*dgds[1]);
        sCorr[1] = -(i10*dgds[0] + i11*dgds[1]);
        s[0] = sOld[0] + lambda*sCorr[0];
        s[1] = sOld[1] + lambda*sCorr[1];

        GH_amp_site_fracs(r, s, &f);

        while (GH_amp_sf_out_of_range(&f) && lambda > DBL_EPSILON){
            lambda /= 2.0;
            s[0] = sOld[0] + lambda*sCorr[0];
            s[1] = sOld[1] + lambda*sCorr[1];
            GH_amp_site_fracs(r, s, &f);
        }

        sNew[0] = s[0]; sNew[1] = s[1];
        iter++;
    }

    s_out[0] = sOld[0]; s_out[1] = sOld[1];
    *sf = f;
}

static double GH_amp_gmix(double t, double p, const double r[2], double dr[2]){
    GH_amp_coef cC, cO;
    GH_amp_sf sfC, sfO, fC, fO;
    double sC[2], sO[2];
    GH_amp_coef_make(1, t, p, &cC);
    GH_amp_coef_make(0, t, p, &cO);
    GH_amp_order(&cC, t, r, sC, &sfC);
    fC = GH_amp_floor_sf(&sfC);
    double gC = GH_amp_gibbs(&cC, r, sC, &fC, t);
    GH_amp_order(&cO, t, r, sO, &sfO);
    fO = GH_amp_floor_sf(&sfO);
    double gO = GH_amp_gibbs(&cO, r, sO, &fO, t);
    if (gO < gC){
        if (dr){ GH_amp_dg_dr(&cO, r, sO, &fO, t, dr); }
        return gO;
    }
    if (dr){ GH_amp_dg_dr(&cC, r, sC, &fC, t, dr); }
    return gC;
}

double obj_gh_amp(unsigned n, const double *x, double *grad, void *SS_ref_db){
    SS_ref *d = (SS_ref *) SS_ref_db;

    double T    = d->T;
    double Pbar = d->P*1000.0;
    double *p   = d->p;
    double *gb  = d->gb_lvl;
    double *mu_Gex = d->mu_Gex;

    for (int i = 0; i < 3; i++){ p[i] = x[i]; }
    double sum = p[0] + p[1] + p[2];
    double r[2] = { 2.0*p[1]/sum - 1.0, p[2]/sum };

    double dr[2];
    double Gex = GH_amp_gmix(T, Pbar, r, dr);

    double dGex[3];
    dGex[0] = 0.0;
    dGex[1] = 2.0*dr[0];
    dGex[2] = dr[1];

    GH_mu_from_dGex(d, p, sum, 3, Gex, dGex);

    GH_sf_amp(p, d->sf);
    return GH_df_assemble(d, p, gb, mu_Gex, 3, grad);
}

void GH_SS_objective_init_function(    obj_type            *SS_objective,
                                        global_variable      gv                  ){
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
        else if (strcmp( gv.SS_list[iss], "oox") == 0 ){
            SS_objective[iss] = obj_gh_oox;
        }
        else if (strcmp( gv.SS_list[iss], "amp") == 0 ){
            SS_objective[iss] = obj_gh_amp;
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
        else if (strcmp( gv.SS_list[iss], "oox") == 0 ){
            PC_read[iss] = obj_gh_oox;
        }
        else if (strcmp( gv.SS_list[iss], "amp") == 0 ){
            PC_read[iss] = obj_gh_amp;
        }
    }
}

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
            strcmp( gv.SS_list[iss], "kls") == 0 ||
            strcmp( gv.SS_list[iss], "oox") == 0 ||
            strcmp( gv.SS_list[iss], "amp") == 0){
            P2X_read[iss] = p2x_gh_generic;
        }
    }
}

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
