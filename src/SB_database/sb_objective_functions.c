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
	List of objective functions used for non-linear minimization and to generate pseudocompounds
*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <complex.h> 

#include "nlopt.h"
#include "../MAGEMin.h"
#include "../toolkit.h"
#include "sb_objective_functions.h"

/**************************************************************************************/
/**************************************************************************************/
/****************SB11 DATABASE (Stixrude & Lithgow-Bertelloni, 2011)*******************/
/**************************************************************************************/
/**************************************************************************************/
/**
    Objective function for sb11_plg
*/
double obj_sb11_plg(unsigned n, const double *x, double *grad, void *SS_ref_db){
    SS_ref *d         = (SS_ref *) SS_ref_db;

    int n_em          = d->n_em;
    double P          = d->P;
    double T          = d->T;
    double R          = d->R;
    double *dfx       = d->dfx;
    double **N        = d->N;
    double *Vec1      = d->Vec1;
    double *Vec2      = d->Vec2;
    double *p         = d->p;
    double *gb        = d->gb_lvl;
    double *mu_Gex    = d->mu_Gex;

    for (int i = 0; i < n_em; i++){
        p[i]   = x[i];
    }

    double tmp = 0.0;
    double Gex = 0.0;
    for (int i = 0; i < n_em; i++){
        Gex = 0.0;
        int it    = 0;
        for (int j = 0; j < d->n_xeos; j++){
            tmp = (d->eye[i][j] - p[j]);
            for (int k = j+1; k < n_em; k++){
                Gex -= tmp*(d->eye[i][k] - p[k])*(d->W[it]);
                it += 1;
            }
        }
        mu_Gex[i] = Gex/1000.0;
    }

    d->sum_apep = 0.0;
    for (int i = 0; i < n_em; i++){
        d->sum_apep += d->ape[i]*p[i];
    }
    d->factor = d->fbc/d->sum_apep;

    double Sconfig    = R*T*(p[0]*log(p[0]) + p[1]*log(p[1]));

    d->df_raw = Sconfig;
    for (int i = 0; i < n_em; i++){
        d->df_raw += (mu_Gex[i] + gb[i])*p[i];
    }
    d->df = d->df_raw * d->factor;

    if (grad){
        grad[0] = (R*T*(1 + log(p[0])) + mu_Gex[0] + gb[0] )* d->factor;
        grad[1] = (R*T*(1 + log(p[1])) + mu_Gex[1] + gb[1] )* d->factor;
    }
    return d->df;
}

/**
    Objective function for sb11_sp
*/
double obj_sb11_sp(unsigned n, const double *x, double *grad, void *SS_ref_db){
    SS_ref *d         = (SS_ref *) SS_ref_db;

    int n_em          = d->n_em;
    double P          = d->P;
    double T          = d->T;
    double R          = d->R;
    double *dfx       = d->dfx;
    double **N        = d->N;
    double *Vec1      = d->Vec1;
    double *Vec2      = d->Vec2;
    double *p         = d->p;
    double *gb        = d->gb_lvl;
    double *mu_Gex    = d->mu_Gex;

    for (int i = 0; i < n_em; i++){
        p[i]   = x[i];
    }

    double tmp = 0.0;
    double Gex = 0.0;
    for (int i = 0; i < n_em; i++){
        Gex = 0.0;
        int it    = 0;
        for (int j = 0; j < d->n_xeos; j++){
            tmp = (d->eye[i][j] - p[j]);
            for (int k = j+1; k < n_em; k++){
                Gex -= tmp*(d->eye[i][k] - p[k])*(d->W[it]);
                it += 1;
            }
        }
        mu_Gex[i] = Gex/1000.0;
    }

    d->sum_apep = 0.0;
    for (int i = 0; i < n_em; i++){
        d->sum_apep += d->ape[i]*p[i];
    }
    d->factor = d->fbc/d->sum_apep;

    double Sconfig    = R*T*(8.0*(0.875*p[0] + 0.875*p[1])*log(0.875*p[0] + 0.875*p[1]) + 3.0*p[0]*log(0.75*p[0]) + p[0]*log(0.125*p[0]) + 4.0*(0.25*p[0] + 0.25*p[1])*log(0.25*p[0] + 0.25*p[1]) + p[1]*log(0.125*p[1]) + 3.0*p[1]*log(0.75*p[1]));

    d->df_raw = Sconfig;
    for (int i = 0; i < n_em; i++){
        d->df_raw += (mu_Gex[i] + gb[i])*p[i];
    }
    d->df = d->df_raw * d->factor;

    if (grad){
        grad[0] = (R*T*(8.0 + 7.0*log(0.875*p[0] + 0.875*p[1]) + 3.0*log(0.75*p[0]) + log(0.125*p[0]) + log(0.25*p[0] + 0.25*p[1]) + 0.125 / 0.125 + 2.25 / 0.75) + mu_Gex[0] + gb[0] )* d->factor;
        grad[1] = (R*T*(8.0 + 7.0*log(0.875*p[0] + 0.875*p[1]) + log(0.125*p[1]) + 3.0*log(0.75*p[1]) + log(0.25*p[0] + 0.25*p[1]) + 0.125 / 0.125 + 2.25 / 0.75) + mu_Gex[1] + gb[1] )* d->factor;
    }
    return d->df;
}

/**
    Objective function for sb11_ol
*/
double obj_sb11_ol(unsigned n, const double *x, double *grad, void *SS_ref_db){
    SS_ref *d         = (SS_ref *) SS_ref_db;

    int n_em          = d->n_em;
    double P          = d->P;
    double T          = d->T;
    double R          = d->R;
    double *dfx       = d->dfx;
    double **N        = d->N;
    double *Vec1      = d->Vec1;
    double *Vec2      = d->Vec2;
    double *p         = d->p;
    double *gb        = d->gb_lvl;
    double *mu_Gex    = d->mu_Gex;

    for (int i = 0; i < n_em; i++){
        p[i]   = x[i];
    }

    double tmp = 0.0;
    double Gex = 0.0;
    for (int i = 0; i < n_em; i++){
        Gex = 0.0;
        int it    = 0;
        for (int j = 0; j < d->n_xeos; j++){
            tmp = (d->eye[i][j] - p[j]);
            for (int k = j+1; k < n_em; k++){
                Gex -= tmp*(d->eye[i][k] - p[k])*(d->W[it]);
                it += 1;
            }
        }
        mu_Gex[i] = Gex/1000.0;
    }

    d->sum_apep = 0.0;
    for (int i = 0; i < n_em; i++){
        d->sum_apep += d->ape[i]*p[i];
    }
    d->factor = d->fbc/d->sum_apep;

    double Sconfig    = R*T*(2.0*p[0]*log(p[0]) + 2.0*p[1]*log(p[1]));

    d->df_raw = Sconfig;
    for (int i = 0; i < n_em; i++){
        d->df_raw += (mu_Gex[i] + gb[i])*p[i];
    }
    d->df = d->df_raw * d->factor;

    if (grad){
        grad[0] = (R*T*(2.0 + 2.0*log(p[0])) + mu_Gex[0] + gb[0] )* d->factor;
        grad[1] = (R*T*(2.0 + 2.0*log(p[1])) + mu_Gex[1] + gb[1] )* d->factor;
    }
    return d->df;
}

/**
    Objective function for sb11_wa
*/
double obj_sb11_wa(unsigned n, const double *x, double *grad, void *SS_ref_db){
    SS_ref *d         = (SS_ref *) SS_ref_db;

    int n_em          = d->n_em;
    double P          = d->P;
    double T          = d->T;
    double R          = d->R;
    double *dfx       = d->dfx;
    double **N        = d->N;
    double *Vec1      = d->Vec1;
    double *Vec2      = d->Vec2;
    double *p         = d->p;
    double *gb        = d->gb_lvl;
    double *mu_Gex    = d->mu_Gex;

    for (int i = 0; i < n_em; i++){
        p[i]   = x[i];
    }

    double tmp = 0.0;
    double Gex = 0.0;
    for (int i = 0; i < n_em; i++){
        Gex = 0.0;
        int it    = 0;
        for (int j = 0; j < d->n_xeos; j++){
            tmp = (d->eye[i][j] - p[j]);
            for (int k = j+1; k < n_em; k++){
                Gex -= tmp*(d->eye[i][k] - p[k])*(d->W[it]);
                it += 1;
            }
        }
        mu_Gex[i] = Gex/1000.0;
    }

    d->sum_apep = 0.0;
    for (int i = 0; i < n_em; i++){
        d->sum_apep += d->ape[i]*p[i];
    }
    d->factor = d->fbc/d->sum_apep;

    double Sconfig    = R*T*(2.0*p[0]*log(p[0]) + 2.0*p[1]*log(p[1]));

    d->df_raw = Sconfig;
    for (int i = 0; i < n_em; i++){
        d->df_raw += (mu_Gex[i] + gb[i])*p[i];
    }
    d->df = d->df_raw * d->factor;

    if (grad){
        grad[0] = (R*T*(2.0 + 2.0*log(p[0])) + mu_Gex[0] + gb[0] )* d->factor;
        grad[1] = (R*T*(2.0 + 2.0*log(p[1])) + mu_Gex[1] + gb[1] )* d->factor;
    }
    return d->df;
}

/**
    Objective function for sb11_ri
*/
double obj_sb11_ri(unsigned n, const double *x, double *grad, void *SS_ref_db){
    SS_ref *d         = (SS_ref *) SS_ref_db;

    int n_em          = d->n_em;
    double P          = d->P;
    double T          = d->T;
    double R          = d->R;
    double *dfx       = d->dfx;
    double **N        = d->N;
    double *Vec1      = d->Vec1;
    double *Vec2      = d->Vec2;
    double *p         = d->p;
    double *gb        = d->gb_lvl;
    double *mu_Gex    = d->mu_Gex;

    for (int i = 0; i < n_em; i++){
        p[i]   = x[i];
    }

    double tmp = 0.0;
    double Gex = 0.0;
    for (int i = 0; i < n_em; i++){
        Gex = 0.0;
        int it    = 0;
        for (int j = 0; j < d->n_xeos; j++){
            tmp = (d->eye[i][j] - p[j]);
            for (int k = j+1; k < n_em; k++){
                Gex -= tmp*(d->eye[i][k] - p[k])*(d->W[it]);
                it += 1;
            }
        }
        mu_Gex[i] = Gex/1000.0;
    }

    d->sum_apep = 0.0;
    for (int i = 0; i < n_em; i++){
        d->sum_apep += d->ape[i]*p[i];
    }
    d->factor = d->fbc/d->sum_apep;

    double Sconfig    = R*T*(2.0*p[0]*log(p[0]) + 2.0*p[1]*log(p[1]));

    d->df_raw = Sconfig;
    for (int i = 0; i < n_em; i++){
        d->df_raw += (mu_Gex[i] + gb[i])*p[i];
    }
    d->df = d->df_raw * d->factor;

    if (grad){
        grad[0] = (R*T*(2.0 + 2.0*log(p[0])) + mu_Gex[0] + gb[0] )* d->factor;
        grad[1] = (R*T*(2.0 + 2.0*log(p[1])) + mu_Gex[1] + gb[1] )* d->factor;
    }
    return d->df;
}

/**
    Objective function for sb11_opx
*/
double obj_sb11_opx(unsigned n, const double *x, double *grad, void *SS_ref_db){
    SS_ref *d         = (SS_ref *) SS_ref_db;

    int n_em          = d->n_em;
    double P          = d->P;
    double T          = d->T;
    double R          = d->R;
    double *dfx       = d->dfx;
    double **N        = d->N;
    double *Vec1      = d->Vec1;
    double *Vec2      = d->Vec2;
    double *p         = d->p;
    double *gb        = d->gb_lvl;
    double *mu_Gex    = d->mu_Gex;

    for (int i = 0; i < n_em; i++){
        p[i]   = x[i];
    }

    double tmp = 0.0;
    double Gex = 0.0;
    for (int i = 0; i < n_em; i++){
        Gex = 0.0;
        int it    = 0;
        for (int j = 0; j < d->n_xeos; j++){
            tmp = (d->eye[i][j] - p[j]);
            for (int k = j+1; k < n_em; k++){
                Gex -= tmp*(d->eye[i][k] - p[k])*(d->W[it]);
                it += 1;
            }
        }
        mu_Gex[i] = Gex/1000.0;
    }

    d->sum_apep = 0.0;
    for (int i = 0; i < n_em; i++){
        d->sum_apep += d->ape[i]*p[i];
    }
    d->factor = d->fbc/d->sum_apep;

    double Sconfig    = R*T*(p[0]*log(p[0]) + (p[0] + p[2])*log(p[0] + p[2]) + 2*p[1]*log(p[1]) + (p[2] + p[3])*log(p[2] + p[3]) + p[3]*log(p[3]));

    d->df_raw = Sconfig;
    for (int i = 0; i < n_em; i++){
        d->df_raw += (mu_Gex[i] + gb[i])*p[i];
    }
    d->df = d->df_raw * d->factor;

    if (grad){
        grad[0] = (R*T*(2 + log(p[0]) + log(p[0] + p[2])) + mu_Gex[0] + gb[0] )* d->factor;
        grad[1] = (R*T*(2 + 2*log(p[1])) + mu_Gex[1] + gb[1] )* d->factor;
        grad[2] = (R*T*(2 + log(p[0] + p[2]) + log(p[2] + p[3])) + mu_Gex[2] + gb[2] )* d->factor;
        grad[3] = (R*T*(2 + log(p[3]) + log(p[2] + p[3])) + mu_Gex[3] + gb[3] )* d->factor;
    }
    return d->df;
}

/**
    Objective function for sb11_cpx
*/
double obj_sb11_cpx(unsigned n, const double *x, double *grad, void *SS_ref_db){
    SS_ref *d         = (SS_ref *) SS_ref_db;

    int n_em          = d->n_em;
    double P          = d->P;
    double T          = d->T;
    double R          = d->R;
    double *dfx       = d->dfx;
    double **N        = d->N;
    double *Vec1      = d->Vec1;
    double *Vec2      = d->Vec2;
    double *p         = d->p;
    double *gb        = d->gb_lvl;
    double *mu_Gex    = d->mu_Gex;
    double *mat_phi   = d->mat_phi;

    for (int i = 0; i < n_em; i++){
        p[i]   = x[i];
    }

    d->sum_v = 0.0;
    for (int i = 0; i < n_em; i++){
        d->sum_v += p[i]*d->v[i];
    }
    for (int i = 0; i < n_em; i++){
        d->mat_phi[i] = (p[i]*d->v[i])/d->sum_v;
    }

    double tmp = 0.0;
    double Gex = 0.0;
    for (int i = 0; i < d->n_em; i++){
        Gex = 0.0;
        int it = 0;
        for (int j = 0; j < d->n_xeos; j++){
            tmp = (d->eye[i][j] - d->mat_phi[j]);
            for (int k = j+1; k < d->n_em; k++){
                Gex -= tmp*(d->eye[i][k] - d->mat_phi[k])*(d->W[it]*2.0*d->v[i]/(d->v[j]+d->v[k]));
                it += 1;
            }
        }
        mu_Gex[i] = Gex;
    }
    d->sum_apep = 0.0;
    for (int i = 0; i < n_em; i++){
        d->sum_apep += d->ape[i]*p[i];
    }
    d->factor = d->fbc/d->sum_apep;

    double Sconfig    = R*T*(p[0]*log(p[0]) + (p[0] + p[3] + p[4])*log(p[0] + p[3] + p[4]) + 2.0*(p[0] + p[1] + p[2] + 0.5*p[3] + p[4])*log(p[0] + p[1] + p[2] + 0.5*p[3] + p[4]) + (p[1] + p[3])*log(p[1] + p[3]) + p[1]*log(p[1]) + (p[2] + p[4])*log(p[2] + p[4]) + p[2]*log(p[2]) + p[3]*log(0.5*p[3]));

    d->df_raw = Sconfig;
    for (int i = 0; i < n_em; i++){
        d->df_raw += (mu_Gex[i] + gb[i])*p[i];
    }
    d->df = d->df_raw * d->factor;

    if (grad){
        grad[0] = (R*T*(4.0 + log(p[0]) + log(p[0] + p[3] + p[4]) + 2.0*log(p[0] + p[1] + p[2] + 0.5*p[3] + p[4])) + mu_Gex[0] + gb[0] )* d->factor;
        grad[1] = (R*T*(4.0 + log(p[1] + p[3]) + 2.0*log(p[0] + p[1] + p[2] + 0.5*p[3] + p[4]) + log(p[1])) + mu_Gex[1] + gb[1] )* d->factor;
        grad[2] = (R*T*(4.0 + log(p[2] + p[4]) + log(p[2]) + 2.0*log(p[0] + p[1] + p[2] + 0.5*p[3] + p[4])) + mu_Gex[2] + gb[2] )* d->factor;
        grad[3] = (R*T*(3 + log(0.5*p[3]) + 0.5 / 0.5 + log(p[0] + p[3] + p[4]) + log(p[1] + p[3]) + log(p[0] + p[1] + p[2] + 0.5*p[3] + p[4])) + mu_Gex[3] + gb[3] )* d->factor;
        grad[4] = (R*T*(4.0 + log(p[0] + p[3] + p[4]) + log(p[2] + p[4]) + 2.0*log(p[0] + p[1] + p[2] + 0.5*p[3] + p[4])) + mu_Gex[4] + gb[4] )* d->factor;
    }
    return d->df;
}

/**
    Objective function for sb11_hpcpx
*/
double obj_sb11_hpcpx(unsigned n, const double *x, double *grad, void *SS_ref_db){
    SS_ref *d         = (SS_ref *) SS_ref_db;

    int n_em          = d->n_em;
    double P          = d->P;
    double T          = d->T;
    double R          = d->R;
    double *dfx       = d->dfx;
    double **N        = d->N;
    double *Vec1      = d->Vec1;
    double *Vec2      = d->Vec2;
    double *p         = d->p;
    double *gb        = d->gb_lvl;
    double *mu_Gex    = d->mu_Gex;

    for (int i = 0; i < n_em; i++){
        p[i]   = x[i];
    }

    double tmp = 0.0;
    double Gex = 0.0;
    for (int i = 0; i < n_em; i++){
        Gex = 0.0;
        int it    = 0;
        for (int j = 0; j < d->n_xeos; j++){
            tmp = (d->eye[i][j] - p[j]);
            for (int k = j+1; k < n_em; k++){
                Gex -= tmp*(d->eye[i][k] - p[k])*(d->W[it]);
                it += 1;
            }
        }
        mu_Gex[i] = Gex/1000.0;
    }

    d->sum_apep = 0.0;
    for (int i = 0; i < n_em; i++){
        d->sum_apep += d->ape[i]*p[i];
    }
    d->factor = d->fbc/d->sum_apep;

    double Sconfig    = R*T*(2.0*p[0]*log(p[0]) + 2.0*p[1]*log(p[1]));

    d->df_raw = Sconfig;
    for (int i = 0; i < n_em; i++){
        d->df_raw += (mu_Gex[i] + gb[i])*p[i];
    }
    d->df = d->df_raw * d->factor;

    if (grad){
        grad[0] = (R*T*(2.0 + 2.0*log(p[0])) + mu_Gex[0] + gb[0] )* d->factor;
        grad[1] = (R*T*(2.0 + 2.0*log(p[1])) + mu_Gex[1] + gb[1] )* d->factor;
    }
    return d->df;
}

/**
    Objective function for sb11_ak
*/
double obj_sb11_ak(unsigned n, const double *x, double *grad, void *SS_ref_db){
    SS_ref *d         = (SS_ref *) SS_ref_db;

    int n_em          = d->n_em;
    double P          = d->P;
    double T          = d->T;
    double R          = d->R;
    double *dfx       = d->dfx;
    double **N        = d->N;
    double *Vec1      = d->Vec1;
    double *Vec2      = d->Vec2;
    double *p         = d->p;
    double *gb        = d->gb_lvl;
    double *mu_Gex    = d->mu_Gex;

    for (int i = 0; i < n_em; i++){
        p[i]   = x[i];
    }

    double tmp = 0.0;
    double Gex = 0.0;
    for (int i = 0; i < n_em; i++){
        Gex = 0.0;
        int it    = 0;
        for (int j = 0; j < d->n_xeos; j++){
            tmp = (d->eye[i][j] - p[j]);
            for (int k = j+1; k < n_em; k++){
                Gex -= tmp*(d->eye[i][k] - p[k])*(d->W[it]);
                it += 1;
            }
        }
        mu_Gex[i] = Gex/1000.0;
    }

    d->sum_apep = 0.0;
    for (int i = 0; i < n_em; i++){
        d->sum_apep += d->ape[i]*p[i];
    }
    d->factor = d->fbc/d->sum_apep;

    double Sconfig    = R*T*(2*p[0]*log(p[0]) + (p[1] + p[2])*log(p[1] + p[2]) + p[1]*log(p[1]) + p[2]*log(p[2]));

    d->df_raw = Sconfig;
    for (int i = 0; i < n_em; i++){
        d->df_raw += (mu_Gex[i] + gb[i])*p[i];
    }
    d->df = d->df_raw * d->factor;

    if (grad){
        grad[0] = (R*T*(2 + 2*log(p[0])) + mu_Gex[0] + gb[0] )* d->factor;
        grad[1] = (R*T*(2 + log(p[1] + p[2]) + log(p[1])) + mu_Gex[1] + gb[1] )* d->factor;
        grad[2] = (R*T*(2 + log(p[1] + p[2]) + log(p[2])) + mu_Gex[2] + gb[2] )* d->factor;
    }
    return d->df;
}

/**
    Objective function for sb11_gtmj
*/
double obj_sb11_gtmj(unsigned n, const double *x, double *grad, void *SS_ref_db){
    SS_ref *d         = (SS_ref *) SS_ref_db;

    int n_em          = d->n_em;
    double P          = d->P;
    double T          = d->T;
    double R          = d->R;
    double *dfx       = d->dfx;
    double **N        = d->N;
    double *Vec1      = d->Vec1;
    double *Vec2      = d->Vec2;
    double *p         = d->p;
    double *gb        = d->gb_lvl;
    double *mu_Gex    = d->mu_Gex;

    for (int i = 0; i < n_em; i++){
        p[i]   = x[i];
    }

    double tmp = 0.0;
    double Gex = 0.0;
    for (int i = 0; i < n_em; i++){
        Gex = 0.0;
        int it    = 0;
        for (int j = 0; j < d->n_xeos; j++){
            tmp = (d->eye[i][j] - p[j]);
            for (int k = j+1; k < n_em; k++){
                Gex -= tmp*(d->eye[i][k] - p[k])*(d->W[it]);
                it += 1;
            }
        }
        mu_Gex[i] = Gex/1000.0;
    }

    d->sum_apep = 0.0;
    for (int i = 0; i < n_em; i++){
        d->sum_apep += d->ape[i]*p[i];
    }
    d->factor = d->fbc/d->sum_apep;

    double Sconfig    = R*T*(3.0*p[0]*log(p[0]) + (p[0] + p[3] + p[4])*log(p[0] + p[3] + p[4]) + (p[0] + p[1] + p[3] + p[4])*log(p[0] + p[1] + p[3] + p[4]) + p[1]*log(0.3333333333333333*p[1]) + (p[1] + p[2])*log(p[1] + p[2]) + 2.0*p[1]*log(0.6666666666666666*p[1]) + p[2]*log(p[2]) + 3.0*(p[2] + p[3])*log(p[2] + p[3]) + 3.0*p[4]*log(p[4]));

    d->df_raw = Sconfig;
    for (int i = 0; i < n_em; i++){
        d->df_raw += (mu_Gex[i] + gb[i])*p[i];
    }
    d->df = d->df_raw * d->factor;

    if (grad){
        grad[0] = (R*T*(5.0 + 3.0*log(p[0]) + log(p[0] + p[3] + p[4]) + log(p[0] + p[1] + p[3] + p[4])) + mu_Gex[0] + gb[0] )* d->factor;
        grad[1] = (R*T*(2 + 0.3333333333333333 / 0.3333333333333333 + log(0.3333333333333333*p[1]) + log(p[1] + p[2]) + log(p[0] + p[1] + p[3] + p[4]) + 1.3333333333333333 / 0.6666666666666666 + 2.0*log(0.6666666666666666*p[1])) + mu_Gex[1] + gb[1] )* d->factor;
        grad[2] = (R*T*(5.0 + log(p[1] + p[2]) + log(p[2]) + 3.0*log(p[2] + p[3])) + mu_Gex[2] + gb[2] )* d->factor;
        grad[3] = (R*T*(5.0 + log(p[0] + p[3] + p[4]) + log(p[0] + p[1] + p[3] + p[4]) + 3.0*log(p[2] + p[3])) + mu_Gex[3] + gb[3] )* d->factor;
        grad[4] = (R*T*(5.0 + log(p[0] + p[3] + p[4]) + log(p[0] + p[1] + p[3] + p[4]) + 3.0*log(p[4])) + mu_Gex[4] + gb[4] )* d->factor;
    }
    return d->df;
}

/**
    Objective function for sb11_pv
*/
double obj_sb11_pv(unsigned n, const double *x, double *grad, void *SS_ref_db){
    SS_ref *d         = (SS_ref *) SS_ref_db;

    int n_em          = d->n_em;
    double P          = d->P;
    double T          = d->T;
    double R          = d->R;
    double *dfx       = d->dfx;
    double **N        = d->N;
    double *Vec1      = d->Vec1;
    double *Vec2      = d->Vec2;
    double *p         = d->p;
    double *gb        = d->gb_lvl;
    double *mu_Gex    = d->mu_Gex;
    double *mat_phi   = d->mat_phi;

    for (int i = 0; i < n_em; i++){
        p[i]   = x[i];
    }

    d->sum_v = 0.0;
    for (int i = 0; i < n_em; i++){
        d->sum_v += p[i]*d->v[i];
    }
    for (int i = 0; i < n_em; i++){
        d->mat_phi[i] = (p[i]*d->v[i])/d->sum_v;
    }

    double tmp = 0.0;
    double Gex = 0.0;
    for (int i = 0; i < d->n_em; i++){
        Gex = 0.0;
        int it = 0;
        for (int j = 0; j < d->n_xeos; j++){
            tmp = (d->eye[i][j] - d->mat_phi[j]);
            for (int k = j+1; k < d->n_em; k++){
                Gex -= tmp*(d->eye[i][k] - d->mat_phi[k])*(d->W[it]*2.0*d->v[i]/(d->v[j]+d->v[k]));
                it += 1;
            }
        }
        mu_Gex[i] = Gex;
    }
    d->sum_apep = 0.0;
    for (int i = 0; i < n_em; i++){
        d->sum_apep += d->ape[i]*p[i];
    }
    d->factor = d->fbc/d->sum_apep;

    double Sconfig    = R*T*(2*p[0]*log(p[0]) + (p[1] + p[2])*log(p[1] + p[2]) + p[1]*log(p[1]) + p[2]*log(p[2]));

    d->df_raw = Sconfig;
    for (int i = 0; i < n_em; i++){
        d->df_raw += (mu_Gex[i] + gb[i])*p[i];
    }
    d->df = d->df_raw * d->factor;

    if (grad){
        grad[0] = (R*T*(2 + 2*log(p[0])) + mu_Gex[0] + gb[0] )* d->factor;
        grad[1] = (R*T*(2 + log(p[1] + p[2]) + log(p[1])) + mu_Gex[1] + gb[1] )* d->factor;
        grad[2] = (R*T*(2 + log(p[1] + p[2]) + log(p[2])) + mu_Gex[2] + gb[2] )* d->factor;
    }
    return d->df;
}

/**
    Objective function for sb11_ppv
*/
double obj_sb11_ppv(unsigned n, const double *x, double *grad, void *SS_ref_db){
    SS_ref *d         = (SS_ref *) SS_ref_db;

    int n_em          = d->n_em;
    double P          = d->P;
    double T          = d->T;
    double R          = d->R;
    double *dfx       = d->dfx;
    double **N        = d->N;
    double *Vec1      = d->Vec1;
    double *Vec2      = d->Vec2;
    double *p         = d->p;
    double *gb        = d->gb_lvl;
    double *mu_Gex    = d->mu_Gex;

    for (int i = 0; i < n_em; i++){
        p[i]   = x[i];
    }

    double tmp = 0.0;
    double Gex = 0.0;
    for (int i = 0; i < n_em; i++){
        Gex = 0.0;
        int it    = 0;
        for (int j = 0; j < d->n_xeos; j++){
            tmp = (d->eye[i][j] - p[j]);
            for (int k = j+1; k < n_em; k++){
                Gex -= tmp*(d->eye[i][k] - p[k])*(d->W[it]);
                it += 1;
            }
        }
        mu_Gex[i] = Gex/1000.0;
    }

    d->sum_apep = 0.0;
    for (int i = 0; i < n_em; i++){
        d->sum_apep += d->ape[i]*p[i];
    }
    d->factor = d->fbc/d->sum_apep;

    double Sconfig    = R*T*(p[0]*log(p[0]) + (p[0] + p[1])*log(p[0] + p[1]) + p[1]*log(p[1]) + 2*p[2]*log(p[2]));

    d->df_raw = Sconfig;
    for (int i = 0; i < n_em; i++){
        d->df_raw += (mu_Gex[i] + gb[i])*p[i];
    }
    d->df = d->df_raw * d->factor;

    if (grad){
        grad[0] = (R*T*(2 + log(p[0]) + log(p[0] + p[1])) + mu_Gex[0] + gb[0] )* d->factor;
        grad[1] = (R*T*(2 + log(p[0] + p[1]) + log(p[1])) + mu_Gex[1] + gb[1] )* d->factor;
        grad[2] = (R*T*(2 + 2*log(p[2])) + mu_Gex[2] + gb[2] )* d->factor;
    }
    return d->df;
}

/**
    Objective function for sb11_mw
*/
double obj_sb11_mw(unsigned n, const double *x, double *grad, void *SS_ref_db){
    SS_ref *d         = (SS_ref *) SS_ref_db;

    int n_em          = d->n_em;
    double P          = d->P;
    double T          = d->T;
    double R          = d->R;
    double *dfx       = d->dfx;
    double **N        = d->N;
    double *Vec1      = d->Vec1;
    double *Vec2      = d->Vec2;
    double *p         = d->p;
    double *gb        = d->gb_lvl;
    double *mu_Gex    = d->mu_Gex;

    for (int i = 0; i < n_em; i++){
        p[i]   = x[i];
    }

    double tmp = 0.0;
    double Gex = 0.0;
    for (int i = 0; i < n_em; i++){
        Gex = 0.0;
        int it    = 0;
        for (int j = 0; j < d->n_xeos; j++){
            tmp = (d->eye[i][j] - p[j]);
            for (int k = j+1; k < n_em; k++){
                Gex -= tmp*(d->eye[i][k] - p[k])*(d->W[it]);
                it += 1;
            }
        }
        mu_Gex[i] = Gex/1000.0;
    }

    d->sum_apep = 0.0;
    for (int i = 0; i < n_em; i++){
        d->sum_apep += d->ape[i]*p[i];
    }
    d->factor = d->fbc/d->sum_apep;

    double Sconfig    = R*T*(p[0]*log(p[0]) + p[1]*log(p[1]));

    d->df_raw = Sconfig;
    for (int i = 0; i < n_em; i++){
        d->df_raw += (mu_Gex[i] + gb[i])*p[i];
    }
    d->df = d->df_raw * d->factor;

    if (grad){
        grad[0] = (R*T*(1 + log(p[0])) + mu_Gex[0] + gb[0] )* d->factor;
        grad[1] = (R*T*(1 + log(p[1])) + mu_Gex[1] + gb[1] )* d->factor;
    }
    return d->df;
}

/**
    Objective function for sb11_cf
*/
double obj_sb11_cf(unsigned n, const double *x, double *grad, void *SS_ref_db){
    SS_ref *d         = (SS_ref *) SS_ref_db;

    int n_em          = d->n_em;
    double P          = d->P;
    double T          = d->T;
    double R          = d->R;
    double *dfx       = d->dfx;
    double **N        = d->N;
    double *Vec1      = d->Vec1;
    double *Vec2      = d->Vec2;
    double *p         = d->p;
    double *gb        = d->gb_lvl;
    double *mu_Gex    = d->mu_Gex;

    for (int i = 0; i < n_em; i++){
        p[i]   = x[i];
    }

    double tmp = 0.0;
    double Gex = 0.0;
    for (int i = 0; i < n_em; i++){
        Gex = 0.0;
        int it    = 0;
        for (int j = 0; j < d->n_xeos; j++){
            tmp = (d->eye[i][j] - p[j]);
            for (int k = j+1; k < n_em; k++){
                Gex -= tmp*(d->eye[i][k] - p[k])*(d->W[it]);
                it += 1;
            }
        }
        mu_Gex[i] = Gex/1000.0;
    }

    d->sum_apep = 0.0;
    for (int i = 0; i < n_em; i++){
        d->sum_apep += d->ape[i]*p[i];
    }
    d->factor = d->fbc/d->sum_apep;

    double Sconfig    = R*T*(p[0]*log(p[0]) + (p[0] + p[2])*log(p[0] + p[2]) + 2*p[1]*log(p[1]) + p[2]*log(p[2]));

    d->df_raw = Sconfig;
    for (int i = 0; i < n_em; i++){
        d->df_raw += (mu_Gex[i] + gb[i])*p[i];
    }
    d->df = d->df_raw * d->factor;

    if (grad){
        grad[0] = (R*T*(2 + log(p[0]) + log(p[0] + p[2])) + mu_Gex[0] + gb[0] )* d->factor;
        grad[1] = (R*T*(2 + 2*log(p[1])) + mu_Gex[1] + gb[1] )* d->factor;
        grad[2] = (R*T*(2 + log(p[2]) + log(p[0] + p[2])) + mu_Gex[2] + gb[2] )* d->factor;
    }
    return d->df;
}


/**
    associate the array of pointer with the right solution phase
*/
void SB_sb11_objective_init_function(	obj_type 			*SS_objective,
                                        global_variable 	 gv				){	
    for (int iss = 0; iss < gv.len_ss; iss++){
        if      (strcmp( gv.SS_list[iss], "plg")  == 0 ){
            SS_objective[iss]  = obj_sb11_plg; 		}
        else if (strcmp( gv.SS_list[iss], "sp")  == 0 ){
            SS_objective[iss]  = obj_sb11_sp; 		}
        else if (strcmp( gv.SS_list[iss], "ol")  == 0 ){
            SS_objective[iss]  = obj_sb11_ol; 		}
        else if (strcmp( gv.SS_list[iss], "wa")  == 0 ){
            SS_objective[iss]  = obj_sb11_wa; 		}
        else if (strcmp( gv.SS_list[iss], "ri")  == 0 ){
            SS_objective[iss]  = obj_sb11_ri; 		}
        else if (strcmp( gv.SS_list[iss], "opx")  == 0 ){
            SS_objective[iss]  = obj_sb11_opx; 		}
        else if (strcmp( gv.SS_list[iss], "cpx")  == 0 ){
            SS_objective[iss]  = obj_sb11_cpx; 		}
        else if (strcmp( gv.SS_list[iss], "hpcpx")  == 0 ){
            SS_objective[iss]  = obj_sb11_hpcpx; 		}
        else if (strcmp( gv.SS_list[iss], "ak")  == 0 ){
            SS_objective[iss]  = obj_sb11_ak; 		}
        else if (strcmp( gv.SS_list[iss], "gtmj")  == 0 ){
            SS_objective[iss]  = obj_sb11_gtmj; 		}
        else if (strcmp( gv.SS_list[iss], "pv")  == 0 ){
            SS_objective[iss]  = obj_sb11_pv; 		}
        else if (strcmp( gv.SS_list[iss], "ppv")  == 0 ){
            SS_objective[iss]  = obj_sb11_ppv; 		}
        else if (strcmp( gv.SS_list[iss], "mw")  == 0 ){
            SS_objective[iss]  = obj_sb11_mw; 		}
        else if (strcmp( gv.SS_list[iss], "cf")  == 0 ){
            SS_objective[iss]  = obj_sb11_cf; 		}
        else{
            printf("\nsolid solution '%s' is not in the database, cannot be initialized\n", gv.SS_list[iss]);
        }	
    };
}

void SB_SS_objective_init_function(	obj_type 			*SS_objective,
									global_variable 	 gv				){

	if (gv.EM_database == 0){				// Stixrude 2011 //
		SB_sb11_objective_init_function(		SS_objective,
												gv							);
	}
	
}

/**
    associate the array of pointer with the right solution phase
*/
void SB_sb11_PC_init(	PC_type 			*PC_read,
                        global_variable 	 gv				){	
    for (int iss = 0; iss < gv.len_ss; iss++){
        if      (strcmp( gv.SS_list[iss], "plg")  == 0 ){
            PC_read[iss]   = obj_sb11_plg; 		}
        else if (strcmp( gv.SS_list[iss], "sp")  == 0 ){
            PC_read[iss]   = obj_sb11_sp; 		}
        else if (strcmp( gv.SS_list[iss], "ol")  == 0 ){
            PC_read[iss]   = obj_sb11_ol; 		}
        else if (strcmp( gv.SS_list[iss], "wa")  == 0 ){
            PC_read[iss]   = obj_sb11_wa; 		}
        else if (strcmp( gv.SS_list[iss], "ri")  == 0 ){
            PC_read[iss]   = obj_sb11_ri; 		}
        else if (strcmp( gv.SS_list[iss], "opx")  == 0 ){
            PC_read[iss]   = obj_sb11_opx; 		}
        else if (strcmp( gv.SS_list[iss], "cpx")  == 0 ){
            PC_read[iss]   = obj_sb11_cpx; 		}
        else if (strcmp( gv.SS_list[iss], "hpcpx")  == 0 ){
            PC_read[iss]   = obj_sb11_hpcpx; 		}
        else if (strcmp( gv.SS_list[iss], "ak")  == 0 ){
            PC_read[iss]   = obj_sb11_ak; 		}
        else if (strcmp( gv.SS_list[iss], "gtmj")  == 0 ){
            PC_read[iss]   = obj_sb11_gtmj; 		}
        else if (strcmp( gv.SS_list[iss], "pv")  == 0 ){
            PC_read[iss]   = obj_sb11_pv; 		}
        else if (strcmp( gv.SS_list[iss], "ppv")  == 0 ){
            PC_read[iss]   = obj_sb11_ppv; 		}
        else if (strcmp( gv.SS_list[iss], "mw")  == 0 ){
            PC_read[iss]   = obj_sb11_mw; 		}
        else if (strcmp( gv.SS_list[iss], "cf")  == 0 ){
            PC_read[iss]   = obj_sb11_cf; 		}
        else{
            printf("\nsolid solution '%s' is not in the database, cannot be initiated\n", gv.SS_list[iss]);
        }	
    };
}


void SB_PC_init(	                PC_type 			*PC_read,
									global_variable 	 gv				){

	if (gv.EM_database == 0){				// metapelite database //
		SB_sb11_PC_init(		PC_read,
								gv							);
	}
	
}



SS_ref SB_PC_function(	global_variable 	 gv,
                        PC_type             *PC_read,
						SS_ref 				 SS_ref_db, 
						bulk_info 	 		 z_b,
						int    			     ph_id				){

	double G0 = 0.0;

    G0 = (*PC_read[ph_id])(		         SS_ref_db.n_xeos, 
                                         SS_ref_db.iguess, 
                                         SS_ref_db.dfx, 
                                        &SS_ref_db		        );	

    SS_ref_db.df = G0;
	
	/** initialize composition */
	for (int j = 0; j < gv.len_ox; j++){
	   SS_ref_db.ss_comp[j] = 0.0;
	}
	
	// /* set mu = 0 for absent oxides */
	// for (int j = 0; j < SS_ref_db.n_em; j++){
	//    SS_ref_db.mu[j] *= SS_ref_db.z_em[j];
	// } 
	
	/* find solution phase composition*/
	for (int i = 0; i < SS_ref_db.n_em; i++){
	   for (int j = 0; j < gv.len_ox; j++){
		   SS_ref_db.ss_comp[j] += SS_ref_db.Comp[i][j]*SS_ref_db.p[i]*SS_ref_db.z_em[i];
	   } 
	}
	
	/* check if site fractions are satisfied */
	SS_ref_db.sf_ok = 1;
	// for (int i = 0; i < SS_ref_db.n_sf; i++){
	// 	if (SS_ref_db.sf[i] < gv.eps_sf_pc || isnan(SS_ref_db.sf[i]) == 1|| isinf(SS_ref_db.sf[i]) == 1){
	// 		SS_ref_db.sf_ok = 0;	
	// 		break;
	// 	}
	// }

    return SS_ref_db;
}

