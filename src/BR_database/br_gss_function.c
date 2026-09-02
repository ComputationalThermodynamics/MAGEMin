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
/* Source: Pourteau et al. (2014), Contrib Mineral Petrol. */
#include <stdio.h>
#include <string.h>

#include "br_gss_function.h"

static int BR_boiled_out(int len_ox, const double *Comp, double *bulk_rock){
    for (int j = 0; j < len_ox; j++){
        if (Comp[j] != 0.0 && bulk_rock[j] == 0.0){
            return 1;
        }
    }
    return 0;
}

SS_ref G_SS_br_ctd_function(SS_ref SS_ref_db, char *research_group, int EM_dataset, int len_ox, bulk_info z_b, double eps){
    strcpy(SS_ref_db.fName,"ctd_BR");
    int i;
    int n_em = SS_ref_db.n_em;

    char *EM_tmp[] = {"fctd","mctd"};
    for (i = 0; i < n_em; i++){
        strcpy(SS_ref_db.EM_list[i], EM_tmp[i]);
    };

    SS_ref_db.W[0] = 0.0;

    em_data fctd_eq = get_em_data(research_group, EM_dataset, len_ox, z_b, SS_ref_db.P, SS_ref_db.T, "fctd", "equilibrium");
    em_data mctd_eq = get_em_data(research_group, EM_dataset, len_ox, z_b, SS_ref_db.P, SS_ref_db.T, "mctd", "equilibrium");

    SS_ref_db.gbase[0] = fctd_eq.gb;
    SS_ref_db.gbase[1] = mctd_eq.gb;

    SS_ref_db.ElShearMod[0] = fctd_eq.ElShearMod;
    SS_ref_db.ElShearMod[1] = mctd_eq.ElShearMod;

    SS_ref_db.ElBulkMod[0] = fctd_eq.ElBulkMod;
    SS_ref_db.ElBulkMod[1] = mctd_eq.ElBulkMod;

    SS_ref_db.ElCp[0] = fctd_eq.ElCp;
    SS_ref_db.ElCp[1] = mctd_eq.ElCp;

    SS_ref_db.ElExpansivity[0] = fctd_eq.ElExpansivity;
    SS_ref_db.ElExpansivity[1] = mctd_eq.ElExpansivity;

    for (i = 0; i < len_ox; i++){
        SS_ref_db.Comp[0][i] = fctd_eq.C[i];
        SS_ref_db.Comp[1][i] = mctd_eq.C[i];
    }

    for (i = 0; i < n_em; i++){
        SS_ref_db.z_em[i] = 1.0;
    };

    for (i = 0; i < SS_ref_db.n_xeos; i++){
        SS_ref_db.bounds_ref[i][0] = 0.0+eps;
        SS_ref_db.bounds_ref[i][1] = 1.0-eps;
    }

    for (i = 0; i < n_em; i++){ SS_ref_db.d_em[i] = 0.0; }
    if (BR_boiled_out(len_ox, fctd_eq.C, z_b.bulk_rock)){ SS_ref_db.d_em[0] = 1.0; SS_ref_db.z_em[0] = 0.0; SS_ref_db.bounds_ref[0][0] = 0.0; SS_ref_db.bounds_ref[0][1] = 0.0; }
    if (BR_boiled_out(len_ox, mctd_eq.C, z_b.bulk_rock)){ SS_ref_db.d_em[1] = 1.0; SS_ref_db.z_em[1] = 0.0; SS_ref_db.bounds_ref[1][0] = 0.0; SS_ref_db.bounds_ref[1][1] = 0.0; }

    return SS_ref_db;
}

SS_ref G_SS_br_car_function(SS_ref SS_ref_db, char *research_group, int EM_dataset, int len_ox, bulk_info z_b, double eps){
    strcpy(SS_ref_db.fName,"car_BR");
    int i;
    int n_em = SS_ref_db.n_em;

    char *EM_tmp[] = {"fcar","mcar"};
    for (i = 0; i < n_em; i++){
        strcpy(SS_ref_db.EM_list[i], EM_tmp[i]);
    };

    SS_ref_db.W[0] = 0.0;

    em_data fcar_eq = get_em_data(research_group, EM_dataset, len_ox, z_b, SS_ref_db.P, SS_ref_db.T, "fcar", "equilibrium");
    em_data mcar_eq = get_em_data(research_group, EM_dataset, len_ox, z_b, SS_ref_db.P, SS_ref_db.T, "mcar", "equilibrium");

    SS_ref_db.gbase[0] = fcar_eq.gb;
    SS_ref_db.gbase[1] = mcar_eq.gb;

    SS_ref_db.ElShearMod[0] = fcar_eq.ElShearMod;
    SS_ref_db.ElShearMod[1] = mcar_eq.ElShearMod;

    SS_ref_db.ElBulkMod[0] = fcar_eq.ElBulkMod;
    SS_ref_db.ElBulkMod[1] = mcar_eq.ElBulkMod;

    SS_ref_db.ElCp[0] = fcar_eq.ElCp;
    SS_ref_db.ElCp[1] = mcar_eq.ElCp;

    SS_ref_db.ElExpansivity[0] = fcar_eq.ElExpansivity;
    SS_ref_db.ElExpansivity[1] = mcar_eq.ElExpansivity;

    for (i = 0; i < len_ox; i++){
        SS_ref_db.Comp[0][i] = fcar_eq.C[i];
        SS_ref_db.Comp[1][i] = mcar_eq.C[i];
    }

    for (i = 0; i < n_em; i++){
        SS_ref_db.z_em[i] = 1.0;
    };

    for (i = 0; i < SS_ref_db.n_xeos; i++){
        SS_ref_db.bounds_ref[i][0] = 0.0+eps;
        SS_ref_db.bounds_ref[i][1] = 1.0-eps;
    }

    for (i = 0; i < n_em; i++){ SS_ref_db.d_em[i] = 0.0; }
    if (BR_boiled_out(len_ox, fcar_eq.C, z_b.bulk_rock)){ SS_ref_db.d_em[0] = 1.0; SS_ref_db.z_em[0] = 0.0; SS_ref_db.bounds_ref[0][0] = 0.0; SS_ref_db.bounds_ref[0][1] = 0.0; }
    if (BR_boiled_out(len_ox, mcar_eq.C, z_b.bulk_rock)){ SS_ref_db.d_em[1] = 1.0; SS_ref_db.z_em[1] = 0.0; SS_ref_db.bounds_ref[1][0] = 0.0; SS_ref_db.bounds_ref[1][1] = 0.0; }

    return SS_ref_db;
}

SS_ref G_SS_br_chl_function(SS_ref SS_ref_db, char *research_group, int EM_dataset, int len_ox, bulk_info z_b, double eps){
    strcpy(SS_ref_db.fName,"chl_BR");
    int i;
    int n_em = SS_ref_db.n_em;

    /* dph, clin, feam, sud, ames */
    char *EM_tmp[] = {"dph","clin","feam","sud","ames"};
    for (i = 0; i < n_em; i++){
        strcpy(SS_ref_db.EM_list[i], EM_tmp[i]);
    };

    double T = SS_ref_db.T;
    double P = SS_ref_db.P;
    /* M1 pairs, order: FeMg, FeAl, Fev, MgAl, Mgv, Alv (WH,WS,WV; FeMg=0) */
    SS_ref_db.W[0] = 0.0;
    SS_ref_db.W[1] = 1200.0    - T*31.0   + P*0.7;
    SS_ref_db.W[2] = 2000.0    - T*(-15.0)+ P*0.4;
    SS_ref_db.W[3] = -9400.0   - T*(-30.0)+ P*(-0.2);
    SS_ref_db.W[4] = 10000.0   - T*(-25.0)+ P*0.9;
    SS_ref_db.W[5] = -10000.0  - T*(-30.0)+ P*0.9;

    em_data dph_eq  = get_em_data(research_group, EM_dataset, len_ox, z_b, SS_ref_db.P, SS_ref_db.T, "dph",  "equilibrium");
    em_data clin_eq = get_em_data(research_group, EM_dataset, len_ox, z_b, SS_ref_db.P, SS_ref_db.T, "clin", "equilibrium");
    em_data feam_eq = get_em_data(research_group, EM_dataset, len_ox, z_b, SS_ref_db.P, SS_ref_db.T, "feam", "equilibrium");
    em_data sud_eq  = get_em_data(research_group, EM_dataset, len_ox, z_b, SS_ref_db.P, SS_ref_db.T, "sud",  "equilibrium");
    em_data ames_eq = get_em_data(research_group, EM_dataset, len_ox, z_b, SS_ref_db.P, SS_ref_db.T, "ames", "equilibrium");
    em_data *em_eq[5] = { &dph_eq, &clin_eq, &feam_eq, &sud_eq, &ames_eq };

    for (i = 0; i < n_em; i++){
        SS_ref_db.gbase[i]         = em_eq[i]->gb;
        SS_ref_db.ElShearMod[i]    = em_eq[i]->ElShearMod;
        SS_ref_db.ElBulkMod[i]     = em_eq[i]->ElBulkMod;
        SS_ref_db.ElCp[i]          = em_eq[i]->ElCp;
        SS_ref_db.ElExpansivity[i] = em_eq[i]->ElExpansivity;
        for (int j = 0; j < len_ox; j++){
            SS_ref_db.Comp[i][j] = em_eq[i]->C[j];
        }
        SS_ref_db.z_em[i] = 1.0;
        SS_ref_db.d_em[i] = 0.0;
    }

    for (i = 0; i < SS_ref_db.n_xeos; i++){
        SS_ref_db.bounds_ref[i][0] = 0.0+eps;
        SS_ref_db.bounds_ref[i][1] = 1.0-eps;
    }

    for (i = 0; i < n_em; i++){
        if (BR_boiled_out(len_ox, em_eq[i]->C, z_b.bulk_rock)){
            SS_ref_db.d_em[i] = 1.0; SS_ref_db.z_em[i] = 0.0;
            SS_ref_db.bounds_ref[i][0] = 0.0; SS_ref_db.bounds_ref[i][1] = 0.0;
        }
    }

    return SS_ref_db;
}

SS_ref G_SS_br_mica_function(SS_ref SS_ref_db, char *research_group, int EM_dataset, int len_ox, bulk_info z_b, double eps){
    strcpy(SS_ref_db.fName,"mica_BR");
    int i;
    int n_em = SS_ref_db.n_em;

    /* mu, pa, cel, fcel, prlph */
    char *EM_tmp[] = {"mu","pa","cel","fcel","prlph"};
    for (i = 0; i < n_em; i++){
        strcpy(SS_ref_db.EM_list[i], EM_tmp[i]);
    };

    double T = SS_ref_db.T;
    double P = SS_ref_db.P;
    /* I-site (K,Na,v), subregular: KKNa,KNaNa,NaNav,Navv,KKv,Kvv,KNav(ternary) */
    SS_ref_db.W[0] = 12230.0  - T*(-5.0)   + P*0.67;
    SS_ref_db.W[1] = 19456.0  - T*(-1.65)  + P*(-0.46);
    SS_ref_db.W[2] = 40000.0  - T*5.0      + P*0.0;
    SS_ref_db.W[3] = 40000.0  - T*5.0      + P*0.0;
    SS_ref_db.W[4] = 35000.0  - T*25.0     + P*0.0;
    SS_ref_db.W[5] = 45000.0  - T*10.0     + P*0.0;
    SS_ref_db.W[6] = 95843.0  - T*19.13    + P*0.105;
    /* M2-site (Mg,Fe,Al), symmetric: MgAl, FeAl (FeMg=0) */
    SS_ref_db.W[7] = -30500.0 - T*15.0     + P*0.78;
    SS_ref_db.W[8] = -5500.0  - T*15.0     + P*0.65;

    em_data mu_eq    = get_em_data(research_group, EM_dataset, len_ox, z_b, SS_ref_db.P, SS_ref_db.T, "mu",    "equilibrium");
    em_data pa_eq    = get_em_data(research_group, EM_dataset, len_ox, z_b, SS_ref_db.P, SS_ref_db.T, "pa",    "equilibrium");
    em_data cel_eq   = get_em_data(research_group, EM_dataset, len_ox, z_b, SS_ref_db.P, SS_ref_db.T, "cel",   "equilibrium");
    em_data fcel_eq  = get_em_data(research_group, EM_dataset, len_ox, z_b, SS_ref_db.P, SS_ref_db.T, "fcel",  "equilibrium");
    em_data prlph_eq = get_em_data(research_group, EM_dataset, len_ox, z_b, SS_ref_db.P, SS_ref_db.T, "prlph", "equilibrium");
    em_data *em_eq[5] = { &mu_eq, &pa_eq, &cel_eq, &fcel_eq, &prlph_eq };

    for (i = 0; i < n_em; i++){
        SS_ref_db.gbase[i]         = em_eq[i]->gb;
        SS_ref_db.ElShearMod[i]    = em_eq[i]->ElShearMod;
        SS_ref_db.ElBulkMod[i]     = em_eq[i]->ElBulkMod;
        SS_ref_db.ElCp[i]          = em_eq[i]->ElCp;
        SS_ref_db.ElExpansivity[i] = em_eq[i]->ElExpansivity;
        for (int j = 0; j < len_ox; j++){
            SS_ref_db.Comp[i][j] = em_eq[i]->C[j];
        }
        SS_ref_db.z_em[i] = 1.0;
        SS_ref_db.d_em[i] = 0.0;
    }

    for (i = 0; i < SS_ref_db.n_xeos; i++){
        SS_ref_db.bounds_ref[i][0] = 0.0+eps;
        SS_ref_db.bounds_ref[i][1] = 1.0-eps;
    }

    for (i = 0; i < n_em; i++){
        if (BR_boiled_out(len_ox, em_eq[i]->C, z_b.bulk_rock)){
            SS_ref_db.d_em[i] = 1.0; SS_ref_db.z_em[i] = 0.0;
            SS_ref_db.bounds_ref[i][0] = 0.0; SS_ref_db.bounds_ref[i][1] = 0.0;
        }
    }

    return SS_ref_db;
}

/* Shared assembler for a 2-endmember binary (any site multiplicity, any W).
   Caller sets fName/EM_tmp/W[0] before calling. */
static SS_ref BR_assemble_binary(SS_ref SS_ref_db, char *research_group, int EM_dataset, int len_ox, bulk_info z_b, double eps, char *em0, char *em1){
    int i;
    int n_em = SS_ref_db.n_em;

    char *EM_tmp[] = {em0, em1};
    for (i = 0; i < n_em; i++){
        strcpy(SS_ref_db.EM_list[i], EM_tmp[i]);
    };

    em_data em0_eq = get_em_data(research_group, EM_dataset, len_ox, z_b, SS_ref_db.P, SS_ref_db.T, em0, "equilibrium");
    em_data em1_eq = get_em_data(research_group, EM_dataset, len_ox, z_b, SS_ref_db.P, SS_ref_db.T, em1, "equilibrium");
    em_data *em_eq[2] = { &em0_eq, &em1_eq };

    for (i = 0; i < n_em; i++){
        SS_ref_db.gbase[i]         = em_eq[i]->gb;
        SS_ref_db.ElShearMod[i]    = em_eq[i]->ElShearMod;
        SS_ref_db.ElBulkMod[i]     = em_eq[i]->ElBulkMod;
        SS_ref_db.ElCp[i]          = em_eq[i]->ElCp;
        SS_ref_db.ElExpansivity[i] = em_eq[i]->ElExpansivity;
        for (int j = 0; j < len_ox; j++){
            SS_ref_db.Comp[i][j] = em_eq[i]->C[j];
        }
        SS_ref_db.z_em[i] = 1.0;
        SS_ref_db.d_em[i] = 0.0;
    }

    for (i = 0; i < SS_ref_db.n_xeos; i++){
        SS_ref_db.bounds_ref[i][0] = 0.0+eps;
        SS_ref_db.bounds_ref[i][1] = 1.0-eps;
    }

    for (i = 0; i < n_em; i++){
        if (BR_boiled_out(len_ox, em_eq[i]->C, z_b.bulk_rock)){
            SS_ref_db.d_em[i] = 1.0; SS_ref_db.z_em[i] = 0.0;
            SS_ref_db.bounds_ref[i][0] = 0.0; SS_ref_db.bounds_ref[i][1] = 0.0;
        }
    }

    return SS_ref_db;
}

SS_ref G_SS_br_talc_function(SS_ref SS_ref_db, char *research_group, int EM_dataset, int len_ox, bulk_info z_b, double eps){
    strcpy(SS_ref_db.fName,"talc_BR");
    SS_ref_db.W[0] = 0.0;
    return BR_assemble_binary(SS_ref_db, research_group, EM_dataset, len_ox, z_b, eps, "ftlc", "mtlc");
}

SS_ref G_SS_br_ilm_function(SS_ref SS_ref_db, char *research_group, int EM_dataset, int len_ox, bulk_info z_b, double eps){
    strcpy(SS_ref_db.fName,"ilm_BR");
    SS_ref_db.W[0] = 0.0;
    return BR_assemble_binary(SS_ref_db, research_group, EM_dataset, len_ox, z_b, eps, "ilm", "gk");
}

SS_ref G_SS_br_bt_function(SS_ref SS_ref_db, char *research_group, int EM_dataset, int len_ox, bulk_info z_b, double eps){
    strcpy(SS_ref_db.fName,"bt_BR");
    SS_ref_db.W[0] = 0.0;
    return BR_assemble_binary(SS_ref_db, research_group, EM_dataset, len_ox, z_b, eps, "phl", "ann");
}

SS_ref G_SS_br_ol_function(SS_ref SS_ref_db, char *research_group, int EM_dataset, int len_ox, bulk_info z_b, double eps){
    strcpy(SS_ref_db.fName,"ol_BR");
    SS_ref_db.W[0] = 0.0;
    return BR_assemble_binary(SS_ref_db, research_group, EM_dataset, len_ox, z_b, eps, "fa", "fo");
}

SS_ref G_SS_br_ep_function(SS_ref SS_ref_db, char *research_group, int EM_dataset, int len_ox, bulk_info z_b, double eps){
    strcpy(SS_ref_db.fName,"ep_BR");
    SS_ref_db.W[0] = 0.0;
    return BR_assemble_binary(SS_ref_db, research_group, EM_dataset, len_ox, z_b, eps, "czo", "ep");
}

SS_ref G_SS_br_opx_function(SS_ref SS_ref_db, char *research_group, int EM_dataset, int len_ox, bulk_info z_b, double eps){
    strcpy(SS_ref_db.fName,"opx_BR");
    SS_ref_db.W[0] = 0.0;
    return BR_assemble_binary(SS_ref_db, research_group, EM_dataset, len_ox, z_b, eps, "en", "fs");
}

SS_ref G_SS_br_amph_function(SS_ref SS_ref_db, char *research_group, int EM_dataset, int len_ox, bulk_info z_b, double eps){
    strcpy(SS_ref_db.fName,"amph_BR");
    SS_ref_db.W[0] = 0.0;
    return BR_assemble_binary(SS_ref_db, research_group, EM_dataset, len_ox, z_b, eps, "parg", "fepa");
}

SS_ref G_SS_br_spl_function(SS_ref SS_ref_db, char *research_group, int EM_dataset, int len_ox, bulk_info z_b, double eps){
    strcpy(SS_ref_db.fName,"spl_BR");
    double T = SS_ref_db.T; double P = SS_ref_db.P;
    SS_ref_db.W[0] = -4661.1144 - T*0.0 + P*0.0;
    return BR_assemble_binary(SS_ref_db, research_group, EM_dataset, len_ox, z_b, eps, "spin", "herc");
}

SS_ref G_SS_br_stau_function(SS_ref SS_ref_db, char *research_group, int EM_dataset, int len_ox, bulk_info z_b, double eps){
    strcpy(SS_ref_db.fName,"stau_BR");
    double T = SS_ref_db.T;
    SS_ref_db.W[0] = -20000.0 - T*0.0;
    return BR_assemble_binary(SS_ref_db, research_group, EM_dataset, len_ox, z_b, eps, "fst", "mst");
}

SS_ref G_SS_br_crd_function(SS_ref SS_ref_db, char *research_group, int EM_dataset, int len_ox, bulk_info z_b, double eps){
    strcpy(SS_ref_db.fName,"crd_BR");
    SS_ref_db.W[0] = 0.0;
    return BR_assemble_binary(SS_ref_db, research_group, EM_dataset, len_ox, z_b, eps, "crd", "fcrd");
}

SS_ref G_SS_br_grt_function(SS_ref SS_ref_db, char *research_group, int EM_dataset, int len_ox, bulk_info z_b, double eps){
    strcpy(SS_ref_db.fName,"grt_BR");
    double P = SS_ref_db.P;
    /* W_112 (py-py-alm), W_122 (py-alm-alm); W_S=0 for both */
    SS_ref_db.W[0] = 230.0  + P*0.01;
    SS_ref_db.W[1] = 3720.0 + P*0.06;
    return BR_assemble_binary(SS_ref_db, research_group, EM_dataset, len_ox, z_b, eps, "py", "alm");
}

SS_ref G_SS_br_omph_function(SS_ref SS_ref_db, char *research_group, int EM_dataset, int len_ox, bulk_info z_b, double eps){
    strcpy(SS_ref_db.fName,"omph_BR");
    int i;
    int n_em = SS_ref_db.n_em;

    /* di, jd, hed */
    char *EM_tmp[] = {"di","jd","hed"};
    for (i = 0; i < n_em; i++){
        strcpy(SS_ref_db.EM_list[i], EM_tmp[i]);
    };

    double T = SS_ref_db.T;
    /* pairwise subregular terms, local pair order as printed in Pourt14:
       W[0..1] = Jd-Di (W_jd2di, W_jddi2), W[2..3] = Di-Hd (W_di2hd, W_dihd2),
       W[4..5] = Jd-Hd (W_jd2hd, W_jdhd2). All W_V = 0. */
    SS_ref_db.W[0] = 25587.5749  - T*15.9243;
    SS_ref_db.W[1] = 5920.2458   - T*10.7164;
    SS_ref_db.W[2] = 31448.6254  - T*28.8272;
    SS_ref_db.W[3] = -65167.1449 - T*(-45.4535);
    SS_ref_db.W[4] = 3670.0      - T*9.0;
    SS_ref_db.W[5] = 15090.0     - T*8.0;

    em_data di_eq  = get_em_data(research_group, EM_dataset, len_ox, z_b, SS_ref_db.P, SS_ref_db.T, "di",  "equilibrium");
    em_data jd_eq  = get_em_data(research_group, EM_dataset, len_ox, z_b, SS_ref_db.P, SS_ref_db.T, "jd",  "equilibrium");
    em_data hed_eq = get_em_data(research_group, EM_dataset, len_ox, z_b, SS_ref_db.P, SS_ref_db.T, "hed", "equilibrium");
    em_data *em_eq[3] = { &di_eq, &jd_eq, &hed_eq };

    for (i = 0; i < n_em; i++){
        SS_ref_db.gbase[i]         = em_eq[i]->gb;
        SS_ref_db.ElShearMod[i]    = em_eq[i]->ElShearMod;
        SS_ref_db.ElBulkMod[i]     = em_eq[i]->ElBulkMod;
        SS_ref_db.ElCp[i]          = em_eq[i]->ElCp;
        SS_ref_db.ElExpansivity[i] = em_eq[i]->ElExpansivity;
        for (int j = 0; j < len_ox; j++){
            SS_ref_db.Comp[i][j] = em_eq[i]->C[j];
        }
        SS_ref_db.z_em[i] = 1.0;
        SS_ref_db.d_em[i] = 0.0;
    }

    for (i = 0; i < SS_ref_db.n_xeos; i++){
        SS_ref_db.bounds_ref[i][0] = 0.0+eps;
        SS_ref_db.bounds_ref[i][1] = 1.0-eps;
    }

    for (i = 0; i < n_em; i++){
        if (BR_boiled_out(len_ox, em_eq[i]->C, z_b.bulk_rock)){
            SS_ref_db.d_em[i] = 1.0; SS_ref_db.z_em[i] = 0.0;
            SS_ref_db.bounds_ref[i][0] = 0.0; SS_ref_db.bounds_ref[i][1] = 0.0;
        }
    }

    return SS_ref_db;
}

SS_ref G_SS_br_amphx_function(SS_ref SS_ref_db, char *research_group, int EM_dataset, int len_ox, bulk_info z_b, double eps){
    strcpy(SS_ref_db.fName,"amphx_BR");
    int i;
    int n_em = SS_ref_db.n_em;

    /* tr, tsch, parg */
    char *EM_tmp[] = {"tr","tsch","parg"};
    for (i = 0; i < n_em; i++){
        strcpy(SS_ref_db.EM_list[i], EM_tmp[i]);
    };

    /* simple symmetric pairs on raw endmember fractions: Tr-Tsch, Tr-Parg, Tsch-Parg */
    SS_ref_db.W[0] = 21431.0;
    SS_ref_db.W[1] = 20459.5;
    SS_ref_db.W[2] = 20459.5;

    em_data tr_eq   = get_em_data(research_group, EM_dataset, len_ox, z_b, SS_ref_db.P, SS_ref_db.T, "tr",   "equilibrium");
    em_data tsch_eq = get_em_data(research_group, EM_dataset, len_ox, z_b, SS_ref_db.P, SS_ref_db.T, "tsch", "equilibrium");
    em_data parg_eq = get_em_data(research_group, EM_dataset, len_ox, z_b, SS_ref_db.P, SS_ref_db.T, "parg", "equilibrium");
    em_data *em_eq[3] = { &tr_eq, &tsch_eq, &parg_eq };

    for (i = 0; i < n_em; i++){
        SS_ref_db.gbase[i]         = em_eq[i]->gb;
        SS_ref_db.ElShearMod[i]    = em_eq[i]->ElShearMod;
        SS_ref_db.ElBulkMod[i]     = em_eq[i]->ElBulkMod;
        SS_ref_db.ElCp[i]          = em_eq[i]->ElCp;
        SS_ref_db.ElExpansivity[i] = em_eq[i]->ElExpansivity;
        for (int j = 0; j < len_ox; j++){
            SS_ref_db.Comp[i][j] = em_eq[i]->C[j];
        }
        SS_ref_db.z_em[i] = 1.0;
        SS_ref_db.d_em[i] = 0.0;
    }

    for (i = 0; i < SS_ref_db.n_xeos; i++){
        SS_ref_db.bounds_ref[i][0] = 0.0+eps;
        SS_ref_db.bounds_ref[i][1] = 1.0-eps;
    }

    for (i = 0; i < n_em; i++){
        if (BR_boiled_out(len_ox, em_eq[i]->C, z_b.bulk_rock)){
            SS_ref_db.d_em[i] = 1.0; SS_ref_db.z_em[i] = 0.0;
            SS_ref_db.bounds_ref[i][0] = 0.0; SS_ref_db.bounds_ref[i][1] = 0.0;
        }
    }

    return SS_ref_db;
}

SS_ref G_SS_br_fsp_function(SS_ref SS_ref_db, char *research_group, int EM_dataset, int len_ox, bulk_info z_b, double eps){
    strcpy(SS_ref_db.fName,"fsp_BR");
    int i;
    int n_em = SS_ref_db.n_em;

    /* ab, kfs, an */
    char *EM_tmp[] = {"ab","kfs","an"};
    for (i = 0; i < n_em; i++){
        strcpy(SS_ref_db.EM_list[i], EM_tmp[i]);
    };

    double T = SS_ref_db.T;
    double P = SS_ref_db.P;
    /* Fuhrman & Lindsley (1988): W1=ab2kfs, W2=abkfs2, W3=ab2an, W4=aban2,
       W5=an2kfs, W6=ankfs2, W7=ternary (ab-an-kfs) */
    SS_ref_db.W[0] = 27320.0   - T*10.3;
    SS_ref_db.W[1] = 18810.0   - T*10.3;
    SS_ref_db.W[2] = 8471.0;
    SS_ref_db.W[3] = 28226.0;
    SS_ref_db.W[4] = 47396.0;
    SS_ref_db.W[5] = 52468.0   + P*(-0.120);
    SS_ref_db.W[6] = 100045.5  - T*10.3    + P*(-0.76);

    em_data ab_eq  = get_em_data(research_group, EM_dataset, len_ox, z_b, SS_ref_db.P, SS_ref_db.T, "ab",  "equilibrium");
    em_data kfs_eq = get_em_data(research_group, EM_dataset, len_ox, z_b, SS_ref_db.P, SS_ref_db.T, "kfs", "equilibrium");
    em_data an_eq  = get_em_data(research_group, EM_dataset, len_ox, z_b, SS_ref_db.P, SS_ref_db.T, "an",  "equilibrium");
    em_data *em_eq[3] = { &ab_eq, &kfs_eq, &an_eq };

    for (i = 0; i < n_em; i++){
        SS_ref_db.gbase[i]         = em_eq[i]->gb;
        SS_ref_db.ElShearMod[i]    = em_eq[i]->ElShearMod;
        SS_ref_db.ElBulkMod[i]     = em_eq[i]->ElBulkMod;
        SS_ref_db.ElCp[i]          = em_eq[i]->ElCp;
        SS_ref_db.ElExpansivity[i] = em_eq[i]->ElExpansivity;
        for (int j = 0; j < len_ox; j++){
            SS_ref_db.Comp[i][j] = em_eq[i]->C[j];
        }
        SS_ref_db.z_em[i] = 1.0;
        SS_ref_db.d_em[i] = 0.0;
    }

    for (i = 0; i < SS_ref_db.n_xeos; i++){
        SS_ref_db.bounds_ref[i][0] = 0.0+eps;
        SS_ref_db.bounds_ref[i][1] = 1.0-eps;
    }

    for (i = 0; i < n_em; i++){
        if (BR_boiled_out(len_ox, em_eq[i]->C, z_b.bulk_rock)){
            SS_ref_db.d_em[i] = 1.0; SS_ref_db.z_em[i] = 0.0;
            SS_ref_db.bounds_ref[i][0] = 0.0; SS_ref_db.bounds_ref[i][1] = 0.0;
        }
    }

    return SS_ref_db;
}

SS_ref G_SS_br_EM_function( global_variable gv, SS_ref SS_ref_db, int EM_dataset, bulk_info z_b, char *name){
    double eps  = gv.bnd_val;
    double P    = SS_ref_db.P;
    double T    = SS_ref_db.T;

    SS_ref_db.ss_flags[0] = 1;

    for (int FD = 0; FD < gv.n_Diff; FD++){
        SS_ref_db.P = P + gv.gb_P_eps*gv.pdev[0][FD];
        SS_ref_db.T = T + gv.gb_T_eps*gv.pdev[1][FD];

        if (strcmp( name, "ctd") == 0 ){
            SS_ref_db = G_SS_br_ctd_function(SS_ref_db, gv.research_group, EM_dataset, gv.len_ox, z_b, eps);
        }
        else if (strcmp( name, "car") == 0 ){
            SS_ref_db = G_SS_br_car_function(SS_ref_db, gv.research_group, EM_dataset, gv.len_ox, z_b, eps);
        }
        else if (strcmp( name, "chl") == 0 ){
            SS_ref_db = G_SS_br_chl_function(SS_ref_db, gv.research_group, EM_dataset, gv.len_ox, z_b, eps);
        }
        else if (strcmp( name, "mica") == 0 ){
            SS_ref_db = G_SS_br_mica_function(SS_ref_db, gv.research_group, EM_dataset, gv.len_ox, z_b, eps);
        }
        else if (strcmp( name, "talc") == 0 ){
            SS_ref_db = G_SS_br_talc_function(SS_ref_db, gv.research_group, EM_dataset, gv.len_ox, z_b, eps);
        }
        else if (strcmp( name, "ilm") == 0 ){
            SS_ref_db = G_SS_br_ilm_function(SS_ref_db, gv.research_group, EM_dataset, gv.len_ox, z_b, eps);
        }
        else if (strcmp( name, "bt") == 0 ){
            SS_ref_db = G_SS_br_bt_function(SS_ref_db, gv.research_group, EM_dataset, gv.len_ox, z_b, eps);
        }
        else if (strcmp( name, "ol") == 0 ){
            SS_ref_db = G_SS_br_ol_function(SS_ref_db, gv.research_group, EM_dataset, gv.len_ox, z_b, eps);
        }
        else if (strcmp( name, "ep") == 0 ){
            SS_ref_db = G_SS_br_ep_function(SS_ref_db, gv.research_group, EM_dataset, gv.len_ox, z_b, eps);
        }
        else if (strcmp( name, "opx") == 0 ){
            SS_ref_db = G_SS_br_opx_function(SS_ref_db, gv.research_group, EM_dataset, gv.len_ox, z_b, eps);
        }
        else if (strcmp( name, "amph") == 0 ){
            SS_ref_db = G_SS_br_amph_function(SS_ref_db, gv.research_group, EM_dataset, gv.len_ox, z_b, eps);
        }
        else if (strcmp( name, "spl") == 0 ){
            SS_ref_db = G_SS_br_spl_function(SS_ref_db, gv.research_group, EM_dataset, gv.len_ox, z_b, eps);
        }
        else if (strcmp( name, "stau") == 0 ){
            SS_ref_db = G_SS_br_stau_function(SS_ref_db, gv.research_group, EM_dataset, gv.len_ox, z_b, eps);
        }
        else if (strcmp( name, "crd") == 0 ){
            SS_ref_db = G_SS_br_crd_function(SS_ref_db, gv.research_group, EM_dataset, gv.len_ox, z_b, eps);
        }
        else if (strcmp( name, "grt") == 0 ){
            SS_ref_db = G_SS_br_grt_function(SS_ref_db, gv.research_group, EM_dataset, gv.len_ox, z_b, eps);
        }
        else if (strcmp( name, "omph") == 0 ){
            SS_ref_db = G_SS_br_omph_function(SS_ref_db, gv.research_group, EM_dataset, gv.len_ox, z_b, eps);
        }
        else if (strcmp( name, "amphx") == 0 ){
            SS_ref_db = G_SS_br_amphx_function(SS_ref_db, gv.research_group, EM_dataset, gv.len_ox, z_b, eps);
        }
        else if (strcmp( name, "fsp") == 0 ){
            SS_ref_db = G_SS_br_fsp_function(SS_ref_db, gv.research_group, EM_dataset, gv.len_ox, z_b, eps);
        }
        else{
            printf("\nsolid solution '%s' is not in the 'br' database\n", name);
        }

        for (int j = 0; j < SS_ref_db.n_em; j++){
            SS_ref_db.mu_array[FD][j] = SS_ref_db.gbase[j];
        }
    }

    for (int j = 0; j < SS_ref_db.n_xeos; j++){
        SS_ref_db.bounds[j][0] = SS_ref_db.bounds_ref[j][0];
        SS_ref_db.bounds[j][1] = SS_ref_db.bounds_ref[j][1];
    }

    for (int i = 0; i < SS_ref_db.n_em; i++){
        SS_ref_db.ape[i] = 0.0;
        for (int j = 0; j < gv.len_ox; j++){
            SS_ref_db.ape[i] += SS_ref_db.Comp[i][j]*z_b.apo[j];
        }
    }

    SS_ref_db.fbc = z_b.fbc;

    return SS_ref_db;
}
