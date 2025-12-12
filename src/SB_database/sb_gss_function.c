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
Function to calculate the reference chemical potential of solid-solutions        

Stixrude thermodynamic database for the mantle minerals
*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <complex.h> 

#include "../MAGEMin.h"
#include "../initialize.h"
#include "../all_solution_phases.h"
#include "../simplex_levelling.h"
#include "../toolkit.h"
/**
    Solution phase data for sb11_plg
*/
SS_ref G_SS_sb11_plg_function(SS_ref SS_ref_db, char* research_group, int EM_dataset, int len_ox, bulk_info z_b, double eps){

    int i, j;
    int n_em = SS_ref_db.n_em;

    char   *EM_tmp[]            = {"an","ab",};
    for (int i = 0; i < SS_ref_db.n_em; i++){
        strcpy(SS_ref_db.EM_list[i],EM_tmp[i]);
    };

    // Site mixing composition;
    SS_ref_db.C[0][0] = 1.0;    SS_ref_db.C[0][1] = 0.0;
    SS_ref_db.C[1][0] = 0.0;    SS_ref_db.C[1][1] = 1.0;

    // pre-computed Nullspace;
    SS_ref_db.N[0][0] = -0.7071067811865475;
    SS_ref_db.N[1][0] = 0.7071067811865476;

    SS_ref_db.W[0] = 26000.0 ;


    em_data an            = get_em_data(    research_group, EM_dataset,
                                            len_ox,
                                            z_b,
                                            SS_ref_db.P,
                                            SS_ref_db.T,
                                            "an",
                                            "equilibrium");

    em_data ab            = get_em_data(    research_group, EM_dataset,
                                            len_ox,
                                            z_b,
                                            SS_ref_db.P,
                                            SS_ref_db.T,
                                            "ab",
                                            "equilibrium");

    SS_ref_db.gbase[0]        = an.gb;
    SS_ref_db.gbase[1]        = ab.gb;

    SS_ref_db.ElShearMod[0]        = an.ElShearMod;
    SS_ref_db.ElShearMod[1]        = ab.ElShearMod;

    SS_ref_db.ElBulkMod[0]        = an.ElBulkMod;
    SS_ref_db.ElBulkMod[1]        = ab.ElBulkMod;

    SS_ref_db.ElCp[0]        = an.ElCp;
    SS_ref_db.ElCp[1]        = ab.ElCp;

    SS_ref_db.ElExpansivity[0]        = an.ElExpansivity;
    SS_ref_db.ElExpansivity[1]        = ab.ElExpansivity;

    for (i = 0; i < len_ox; i++){
        SS_ref_db.Comp[0][i]    = an.C[i];
        SS_ref_db.Comp[1][i]    = ab.C[i];
    }

    for (i = 0; i < n_em; i++){
        SS_ref_db.z_em[i] = 1.0;
    };

    SS_ref_db.bounds_ref[0][0] = 0.0+eps;  SS_ref_db.bounds_ref[0][1] = 1.0-eps;
    SS_ref_db.bounds_ref[1][0] = 0.0+eps;  SS_ref_db.bounds_ref[1][1] = 1.0-eps;


    return SS_ref_db;
}

/**
    Solution phase data for sb11_sp
*/
SS_ref G_SS_sb11_sp_function(SS_ref SS_ref_db, char* research_group, int EM_dataset, int len_ox, bulk_info z_b, double eps){

    int i, j;
    int n_em = SS_ref_db.n_em;

    char   *EM_tmp[]            = {"sp","hc",};
    for (int i = 0; i < SS_ref_db.n_em; i++){
        strcpy(SS_ref_db.EM_list[i],EM_tmp[i]);
    };

    // Site mixing composition;
    SS_ref_db.C[0][0] = 0.25;    SS_ref_db.C[0][1] = 0.25;
    SS_ref_db.C[1][0] = 0.0;    SS_ref_db.C[1][1] = 0.75;
    SS_ref_db.C[2][0] = 0.75;    SS_ref_db.C[2][1] = 0.0;
    SS_ref_db.C[3][0] = 0.875;    SS_ref_db.C[3][1] = 0.875;
    SS_ref_db.C[4][0] = 0.0;    SS_ref_db.C[4][1] = 0.125;
    SS_ref_db.C[5][0] = 0.125;    SS_ref_db.C[5][1] = 0.0;

    // pre-computed Nullspace;
    SS_ref_db.N[0][0] = -0.7071067811865475;
    SS_ref_db.N[1][0] = 0.7071067811865476;

    SS_ref_db.W[0] = 5876.46 ;


    em_data sp            = get_em_data(    research_group, EM_dataset,
                                            len_ox,
                                            z_b,
                                            SS_ref_db.P,
                                            SS_ref_db.T,
                                            "sp",
                                            "equilibrium");

    em_data hc            = get_em_data(    research_group, EM_dataset,
                                            len_ox,
                                            z_b,
                                            SS_ref_db.P,
                                            SS_ref_db.T,
                                            "hc",
                                            "equilibrium");

    SS_ref_db.gbase[0]        = sp.gb;
    SS_ref_db.gbase[1]        = hc.gb;

    SS_ref_db.ElShearMod[0]        = sp.ElShearMod;
    SS_ref_db.ElShearMod[1]        = hc.ElShearMod;

    SS_ref_db.ElBulkMod[0]        = sp.ElBulkMod;
    SS_ref_db.ElBulkMod[1]        = hc.ElBulkMod;

    SS_ref_db.ElCp[0]        = sp.ElCp;
    SS_ref_db.ElCp[1]        = hc.ElCp;

    SS_ref_db.ElExpansivity[0]        = sp.ElExpansivity;
    SS_ref_db.ElExpansivity[1]        = hc.ElExpansivity;

    for (i = 0; i < len_ox; i++){
        SS_ref_db.Comp[0][i]    = sp.C[i];
        SS_ref_db.Comp[1][i]    = hc.C[i];
    }

    for (i = 0; i < n_em; i++){
        SS_ref_db.z_em[i] = 1.0;
    };

    SS_ref_db.bounds_ref[0][0] = 0.0+eps;  SS_ref_db.bounds_ref[0][1] = 1.0-eps;
    SS_ref_db.bounds_ref[1][0] = 0.0+eps;  SS_ref_db.bounds_ref[1][1] = 1.0-eps;


    return SS_ref_db;
}

/**
    Solution phase data for sb11_ol
*/
SS_ref G_SS_sb11_ol_function(SS_ref SS_ref_db, char* research_group, int EM_dataset, int len_ox, bulk_info z_b, double eps){

    int i, j;
    int n_em = SS_ref_db.n_em;

    char   *EM_tmp[]            = {"fo","fa",};
    for (int i = 0; i < SS_ref_db.n_em; i++){
        strcpy(SS_ref_db.EM_list[i],EM_tmp[i]);
    };

    // Site mixing composition;
    SS_ref_db.C[0][0] = 0.0;    SS_ref_db.C[0][1] = 1.0;
    SS_ref_db.C[1][0] = 1.0;    SS_ref_db.C[1][1] = 0.0;

    // pre-computed Nullspace;
    SS_ref_db.N[0][0] = -0.7071067811865475;
    SS_ref_db.N[1][0] = 0.7071067811865476;

    SS_ref_db.W[0] = 7813.22 ;


    em_data fo            = get_em_data(    research_group, EM_dataset,
                                            len_ox,
                                            z_b,
                                            SS_ref_db.P,
                                            SS_ref_db.T,
                                            "fo",
                                            "equilibrium");

    em_data fa            = get_em_data(    research_group, EM_dataset,
                                            len_ox,
                                            z_b,
                                            SS_ref_db.P,
                                            SS_ref_db.T,
                                            "fa",
                                            "equilibrium");

    SS_ref_db.gbase[0]        = fo.gb;
    SS_ref_db.gbase[1]        = fa.gb;

    SS_ref_db.ElShearMod[0]        = fo.ElShearMod;
    SS_ref_db.ElShearMod[1]        = fa.ElShearMod;

    SS_ref_db.ElBulkMod[0]        = fo.ElBulkMod;
    SS_ref_db.ElBulkMod[1]        = fa.ElBulkMod;

    SS_ref_db.ElCp[0]        = fo.ElCp;
    SS_ref_db.ElCp[1]        = fa.ElCp;

    SS_ref_db.ElExpansivity[0]        = fo.ElExpansivity;
    SS_ref_db.ElExpansivity[1]        = fa.ElExpansivity;

    for (i = 0; i < len_ox; i++){
        SS_ref_db.Comp[0][i]    = fo.C[i];
        SS_ref_db.Comp[1][i]    = fa.C[i];
    }

    for (i = 0; i < n_em; i++){
        SS_ref_db.z_em[i] = 1.0;
    };

    SS_ref_db.bounds_ref[0][0] = 0.0+eps;  SS_ref_db.bounds_ref[0][1] = 1.0-eps;
    SS_ref_db.bounds_ref[1][0] = 0.0+eps;  SS_ref_db.bounds_ref[1][1] = 1.0-eps;


    return SS_ref_db;
}

/**
    Solution phase data for sb11_wa
*/
SS_ref G_SS_sb11_wa_function(SS_ref SS_ref_db, char* research_group, int EM_dataset, int len_ox, bulk_info z_b, double eps){

    int i, j;
    int n_em = SS_ref_db.n_em;

    char   *EM_tmp[]            = {"mgwa","fewa",};
    for (int i = 0; i < SS_ref_db.n_em; i++){
        strcpy(SS_ref_db.EM_list[i],EM_tmp[i]);
    };

    // Site mixing composition;
    SS_ref_db.C[0][0] = 0.0;    SS_ref_db.C[0][1] = 1.0;
    SS_ref_db.C[1][0] = 1.0;    SS_ref_db.C[1][1] = 0.0;

    // pre-computed Nullspace;
    SS_ref_db.N[0][0] = -0.7071067811865475;
    SS_ref_db.N[1][0] = 0.7071067811865476;

    SS_ref_db.W[0] = 16747.18 ;


    em_data mgwa            = get_em_data(    research_group, EM_dataset,
                                            len_ox,
                                            z_b,
                                            SS_ref_db.P,
                                            SS_ref_db.T,
                                            "mgwa",
                                            "equilibrium");

    em_data fewa            = get_em_data(    research_group, EM_dataset,
                                            len_ox,
                                            z_b,
                                            SS_ref_db.P,
                                            SS_ref_db.T,
                                            "fewa",
                                            "equilibrium");

    SS_ref_db.gbase[0]        = mgwa.gb;
    SS_ref_db.gbase[1]        = fewa.gb;

    SS_ref_db.ElShearMod[0]        = mgwa.ElShearMod;
    SS_ref_db.ElShearMod[1]        = fewa.ElShearMod;

    SS_ref_db.ElBulkMod[0]        = mgwa.ElBulkMod;
    SS_ref_db.ElBulkMod[1]        = fewa.ElBulkMod;

    SS_ref_db.ElCp[0]        = mgwa.ElCp;
    SS_ref_db.ElCp[1]        = fewa.ElCp;

    SS_ref_db.ElExpansivity[0]        = mgwa.ElExpansivity;
    SS_ref_db.ElExpansivity[1]        = fewa.ElExpansivity;

    for (i = 0; i < len_ox; i++){
        SS_ref_db.Comp[0][i]    = mgwa.C[i];
        SS_ref_db.Comp[1][i]    = fewa.C[i];
    }

    for (i = 0; i < n_em; i++){
        SS_ref_db.z_em[i] = 1.0;
    };

    SS_ref_db.bounds_ref[0][0] = 0.0+eps;  SS_ref_db.bounds_ref[0][1] = 1.0-eps;
    SS_ref_db.bounds_ref[1][0] = 0.0+eps;  SS_ref_db.bounds_ref[1][1] = 1.0-eps;


    return SS_ref_db;
}

/**
    Solution phase data for sb11_ri
*/
SS_ref G_SS_sb11_ri_function(SS_ref SS_ref_db, char* research_group, int EM_dataset, int len_ox, bulk_info z_b, double eps){

    int i, j;
    int n_em = SS_ref_db.n_em;

    char   *EM_tmp[]            = {"mgri","feri",};
    for (int i = 0; i < SS_ref_db.n_em; i++){
        strcpy(SS_ref_db.EM_list[i],EM_tmp[i]);
    };

    // Site mixing composition;
    SS_ref_db.C[0][0] = 0.0;    SS_ref_db.C[0][1] = 1.0;
    SS_ref_db.C[1][0] = 1.0;    SS_ref_db.C[1][1] = 0.0;

    // pre-computed Nullspace;
    SS_ref_db.N[0][0] = -0.7071067811865475;
    SS_ref_db.N[1][0] = 0.7071067811865476;

    SS_ref_db.W[0] = 9340.84 ;


    em_data mgri            = get_em_data(    research_group, EM_dataset,
                                            len_ox,
                                            z_b,
                                            SS_ref_db.P,
                                            SS_ref_db.T,
                                            "mgri",
                                            "equilibrium");

    em_data feri            = get_em_data(    research_group, EM_dataset,
                                            len_ox,
                                            z_b,
                                            SS_ref_db.P,
                                            SS_ref_db.T,
                                            "feri",
                                            "equilibrium");

    SS_ref_db.gbase[0]        = mgri.gb;
    SS_ref_db.gbase[1]        = feri.gb;

    SS_ref_db.ElShearMod[0]        = mgri.ElShearMod;
    SS_ref_db.ElShearMod[1]        = feri.ElShearMod;

    SS_ref_db.ElBulkMod[0]        = mgri.ElBulkMod;
    SS_ref_db.ElBulkMod[1]        = feri.ElBulkMod;

    SS_ref_db.ElCp[0]        = mgri.ElCp;
    SS_ref_db.ElCp[1]        = feri.ElCp;

    SS_ref_db.ElExpansivity[0]        = mgri.ElExpansivity;
    SS_ref_db.ElExpansivity[1]        = feri.ElExpansivity;

    for (i = 0; i < len_ox; i++){
        SS_ref_db.Comp[0][i]    = mgri.C[i];
        SS_ref_db.Comp[1][i]    = feri.C[i];
    }

    for (i = 0; i < n_em; i++){
        SS_ref_db.z_em[i] = 1.0;
    };

    SS_ref_db.bounds_ref[0][0] = 0.0+eps;  SS_ref_db.bounds_ref[0][1] = 1.0-eps;
    SS_ref_db.bounds_ref[1][0] = 0.0+eps;  SS_ref_db.bounds_ref[1][1] = 1.0-eps;


    return SS_ref_db;
}

/**
    Solution phase data for sb11_opx
*/
SS_ref G_SS_sb11_opx_function(SS_ref SS_ref_db, char* research_group, int EM_dataset, int len_ox, bulk_info z_b, double eps){

    int i, j;
    int n_em = SS_ref_db.n_em;

    char   *EM_tmp[]            = {"en","fs","mgts","odi",};
    for (int i = 0; i < SS_ref_db.n_em; i++){
        strcpy(SS_ref_db.EM_list[i],EM_tmp[i]);
    };

    // Site mixing composition;
    SS_ref_db.C[0][0] = 0.0;    SS_ref_db.C[0][1] = 0.0;    SS_ref_db.C[0][2] = 0.0;    SS_ref_db.C[0][3] = 1.0;
    SS_ref_db.C[1][0] = 0.0;    SS_ref_db.C[1][1] = 1.0;    SS_ref_db.C[1][2] = 0.0;    SS_ref_db.C[1][3] = 0.0;
    SS_ref_db.C[2][0] = 1.0;    SS_ref_db.C[2][1] = 0.0;    SS_ref_db.C[2][2] = 1.0;    SS_ref_db.C[2][3] = 0.0;
    SS_ref_db.C[3][0] = 0.0;    SS_ref_db.C[3][1] = 0.0;    SS_ref_db.C[3][2] = 1.0;    SS_ref_db.C[3][3] = 0.0;
    SS_ref_db.C[4][0] = 0.0;    SS_ref_db.C[4][1] = 1.0;    SS_ref_db.C[4][2] = 0.0;    SS_ref_db.C[4][3] = 0.0;
    SS_ref_db.C[5][0] = 1.0;    SS_ref_db.C[5][1] = 0.0;    SS_ref_db.C[5][2] = 0.0;    SS_ref_db.C[5][3] = 1.0;

    // pre-computed Nullspace;
    SS_ref_db.N[0][0] = -0.5;    SS_ref_db.N[0][1] = -0.5;    SS_ref_db.N[0][2] = -0.5;
    SS_ref_db.N[1][0] = 0.8333333333333334;    SS_ref_db.N[1][1] = -0.16666666666666666;    SS_ref_db.N[1][2] = -0.16666666666666666;
    SS_ref_db.N[2][0] = -0.16666666666666666;    SS_ref_db.N[2][1] = 0.8333333333333334;    SS_ref_db.N[2][2] = -0.16666666666666666;
    SS_ref_db.N[3][0] = -0.16666666666666666;    SS_ref_db.N[3][1] = -0.16666666666666666;    SS_ref_db.N[3][2] = 0.8333333333333334;

    SS_ref_db.W[0] = 0.0 ;
    SS_ref_db.W[1] = 0.0 ;
    SS_ref_db.W[2] = 32113.52 ;
    SS_ref_db.W[3] = 0.0 ;
    SS_ref_db.W[4] = 0.0 ;
    SS_ref_db.W[5] = 48353.16 ;


    em_data en            = get_em_data(    research_group, EM_dataset,
                                            len_ox,
                                            z_b,
                                            SS_ref_db.P,
                                            SS_ref_db.T,
                                            "en",
                                            "equilibrium");

    em_data fs            = get_em_data(    research_group, EM_dataset,
                                            len_ox,
                                            z_b,
                                            SS_ref_db.P,
                                            SS_ref_db.T,
                                            "fs",
                                            "equilibrium");

    em_data mgts            = get_em_data(    research_group, EM_dataset,
                                            len_ox,
                                            z_b,
                                            SS_ref_db.P,
                                            SS_ref_db.T,
                                            "mgts",
                                            "equilibrium");

    em_data odi            = get_em_data(    research_group, EM_dataset,
                                            len_ox,
                                            z_b,
                                            SS_ref_db.P,
                                            SS_ref_db.T,
                                            "odi",
                                            "equilibrium");

    SS_ref_db.gbase[0]        = en.gb;
    SS_ref_db.gbase[1]        = fs.gb;
    SS_ref_db.gbase[2]        = mgts.gb;
    SS_ref_db.gbase[3]        = odi.gb;

    SS_ref_db.ElShearMod[0]        = en.ElShearMod;
    SS_ref_db.ElShearMod[1]        = fs.ElShearMod;
    SS_ref_db.ElShearMod[2]        = mgts.ElShearMod;
    SS_ref_db.ElShearMod[3]        = odi.ElShearMod;

    SS_ref_db.ElBulkMod[0]        = en.ElBulkMod;
    SS_ref_db.ElBulkMod[1]        = fs.ElBulkMod;
    SS_ref_db.ElBulkMod[2]        = mgts.ElBulkMod;
    SS_ref_db.ElBulkMod[3]        = odi.ElBulkMod;

    SS_ref_db.ElCp[0]        = en.ElCp;
    SS_ref_db.ElCp[1]        = fs.ElCp;
    SS_ref_db.ElCp[2]        = mgts.ElCp;
    SS_ref_db.ElCp[3]        = odi.ElCp;

    SS_ref_db.ElExpansivity[0]        = en.ElExpansivity;
    SS_ref_db.ElExpansivity[1]        = fs.ElExpansivity;
    SS_ref_db.ElExpansivity[2]        = mgts.ElExpansivity;
    SS_ref_db.ElExpansivity[3]        = odi.ElExpansivity;

    for (i = 0; i < len_ox; i++){
        SS_ref_db.Comp[0][i]    = en.C[i];
        SS_ref_db.Comp[1][i]    = fs.C[i];
        SS_ref_db.Comp[2][i]    = mgts.C[i];
        SS_ref_db.Comp[3][i]    = odi.C[i];
    }

    for (i = 0; i < n_em; i++){
        SS_ref_db.z_em[i] = 1.0;
    };

    SS_ref_db.bounds_ref[0][0] = 0.0+eps;  SS_ref_db.bounds_ref[0][1] = 1.0-eps;
    SS_ref_db.bounds_ref[1][0] = 0.0+eps;  SS_ref_db.bounds_ref[1][1] = 1.0-eps;
    SS_ref_db.bounds_ref[2][0] = 0.0+eps;  SS_ref_db.bounds_ref[2][1] = 1.0-eps;
    SS_ref_db.bounds_ref[3][0] = 0.0+eps;  SS_ref_db.bounds_ref[3][1] = 1.0-eps;


    return SS_ref_db;
}

/**
    Solution phase data for sb11_cpx
*/
SS_ref G_SS_sb11_cpx_function(SS_ref SS_ref_db, char* research_group, int EM_dataset, int len_ox, bulk_info z_b, double eps){

    int i, j;
    int n_em = SS_ref_db.n_em;

    char   *EM_tmp[]            = {"di","he","cen","cats","jd",};
    for (int i = 0; i < SS_ref_db.n_em; i++){
        strcpy(SS_ref_db.EM_list[i],EM_tmp[i]);
    };

    // Site mixing composition;
    SS_ref_db.C[0][0] = 1.0;    SS_ref_db.C[0][1] = 1.0;    SS_ref_db.C[0][2] = 0.0;    SS_ref_db.C[0][3] = 1.0;    SS_ref_db.C[0][4] = 0.0;
    SS_ref_db.C[1][0] = 0.0;    SS_ref_db.C[1][1] = 0.0;    SS_ref_db.C[1][2] = 1.0;    SS_ref_db.C[1][3] = 0.0;    SS_ref_db.C[1][4] = 0.0;
    SS_ref_db.C[2][0] = 0.0;    SS_ref_db.C[2][1] = 0.0;    SS_ref_db.C[2][2] = 0.0;    SS_ref_db.C[2][3] = 0.0;    SS_ref_db.C[2][4] = 1.0;
    SS_ref_db.C[3][0] = 0.0;    SS_ref_db.C[3][1] = 0.0;    SS_ref_db.C[3][2] = 0.0;    SS_ref_db.C[3][3] = 1.0;    SS_ref_db.C[3][4] = 1.0;
    SS_ref_db.C[4][0] = 0.0;    SS_ref_db.C[4][1] = 1.0;    SS_ref_db.C[4][2] = 0.0;    SS_ref_db.C[4][3] = 0.0;    SS_ref_db.C[4][4] = 0.0;
    SS_ref_db.C[5][0] = 1.0;    SS_ref_db.C[5][1] = 0.0;    SS_ref_db.C[5][2] = 1.0;    SS_ref_db.C[5][3] = 0.0;    SS_ref_db.C[5][4] = 0.0;
    SS_ref_db.C[6][0] = 1.0;    SS_ref_db.C[6][1] = 1.0;    SS_ref_db.C[6][2] = 1.0;    SS_ref_db.C[6][3] = 0.5;    SS_ref_db.C[6][4] = 1.0;
    SS_ref_db.C[7][0] = 0.0;    SS_ref_db.C[7][1] = 0.0;    SS_ref_db.C[7][2] = 0.0;    SS_ref_db.C[7][3] = 0.5;    SS_ref_db.C[7][4] = 0.0;

    // pre-computed Nullspace;
    SS_ref_db.N[0][0] = -0.44721359549995787;    SS_ref_db.N[0][1] = -0.44721359549995787;    SS_ref_db.N[0][2] = -0.44721359549995787;    SS_ref_db.N[0][3] = -0.44721359549995787;
    SS_ref_db.N[1][0] = 0.8618033988749895;    SS_ref_db.N[1][1] = -0.13819660112501048;    SS_ref_db.N[1][2] = -0.13819660112501048;    SS_ref_db.N[1][3] = -0.13819660112501048;
    SS_ref_db.N[2][0] = -0.13819660112501048;    SS_ref_db.N[2][1] = 0.8618033988749895;    SS_ref_db.N[2][2] = -0.13819660112501048;    SS_ref_db.N[2][3] = -0.13819660112501048;
    SS_ref_db.N[3][0] = -0.13819660112501048;    SS_ref_db.N[3][1] = -0.13819660112501048;    SS_ref_db.N[3][2] = 0.8618033988749895;    SS_ref_db.N[3][3] = -0.13819660112501048;
    SS_ref_db.N[4][0] = -0.13819660112501048;    SS_ref_db.N[4][1] = -0.13819660112501048;    SS_ref_db.N[4][2] = -0.13819660112501048;    SS_ref_db.N[4][3] = 0.8618033988749895;

    SS_ref_db.W[0] = 0.0 ;
    SS_ref_db.W[1] = 24740.0 ;
    SS_ref_db.W[2] = 26000.0 ;
    SS_ref_db.W[3] = 24300.0 ;
    SS_ref_db.W[4] = 24740.0 ;
    SS_ref_db.W[5] = 0.0 ;
    SS_ref_db.W[6] = 0.0 ;
    SS_ref_db.W[7] = 60531.36 ;
    SS_ref_db.W[8] = 0.0 ;
    SS_ref_db.W[9] = 10000.0 ;

    SS_ref_db.v[0] = 1.0;
    SS_ref_db.v[1] = 1.0;
    SS_ref_db.v[2] = 1.0;
    SS_ref_db.v[3] = 3.5;
    SS_ref_db.v[4] = 1.0;

    em_data di            = get_em_data(    research_group, EM_dataset,
                                            len_ox,
                                            z_b,
                                            SS_ref_db.P,
                                            SS_ref_db.T,
                                            "di",
                                            "equilibrium");

    em_data he            = get_em_data(    research_group, EM_dataset,
                                            len_ox,
                                            z_b,
                                            SS_ref_db.P,
                                            SS_ref_db.T,
                                            "he",
                                            "equilibrium");

    em_data cen            = get_em_data(    research_group, EM_dataset,
                                            len_ox,
                                            z_b,
                                            SS_ref_db.P,
                                            SS_ref_db.T,
                                            "cen",
                                            "equilibrium");

    em_data cats            = get_em_data(    research_group, EM_dataset,
                                            len_ox,
                                            z_b,
                                            SS_ref_db.P,
                                            SS_ref_db.T,
                                            "cats",
                                            "equilibrium");

    em_data jd            = get_em_data(    research_group, EM_dataset,
                                            len_ox,
                                            z_b,
                                            SS_ref_db.P,
                                            SS_ref_db.T,
                                            "jd",
                                            "equilibrium");

    SS_ref_db.gbase[0]        = di.gb;
    SS_ref_db.gbase[1]        = he.gb;
    SS_ref_db.gbase[2]        = cen.gb;
    SS_ref_db.gbase[3]        = cats.gb;
    SS_ref_db.gbase[4]        = jd.gb;

    SS_ref_db.ElShearMod[0]        = di.ElShearMod;
    SS_ref_db.ElShearMod[1]        = he.ElShearMod;
    SS_ref_db.ElShearMod[2]        = cen.ElShearMod;
    SS_ref_db.ElShearMod[3]        = cats.ElShearMod;
    SS_ref_db.ElShearMod[4]        = jd.ElShearMod;

    SS_ref_db.ElBulkMod[0]        = di.ElBulkMod;
    SS_ref_db.ElBulkMod[1]        = he.ElBulkMod;
    SS_ref_db.ElBulkMod[2]        = cen.ElBulkMod;
    SS_ref_db.ElBulkMod[3]        = cats.ElBulkMod;
    SS_ref_db.ElBulkMod[4]        = jd.ElBulkMod;

    SS_ref_db.ElCp[0]        = di.ElCp;
    SS_ref_db.ElCp[1]        = he.ElCp;
    SS_ref_db.ElCp[2]        = cen.ElCp;
    SS_ref_db.ElCp[3]        = cats.ElCp;
    SS_ref_db.ElCp[4]        = jd.ElCp;

    SS_ref_db.ElExpansivity[0]        = di.ElExpansivity;
    SS_ref_db.ElExpansivity[1]        = he.ElExpansivity;
    SS_ref_db.ElExpansivity[2]        = cen.ElExpansivity;
    SS_ref_db.ElExpansivity[3]        = cats.ElExpansivity;
    SS_ref_db.ElExpansivity[4]        = jd.ElExpansivity;

    for (i = 0; i < len_ox; i++){
        SS_ref_db.Comp[0][i]    = di.C[i];
        SS_ref_db.Comp[1][i]    = he.C[i];
        SS_ref_db.Comp[2][i]    = cen.C[i];
        SS_ref_db.Comp[3][i]    = cats.C[i];
        SS_ref_db.Comp[4][i]    = jd.C[i];
    }

    for (i = 0; i < n_em; i++){
        SS_ref_db.z_em[i] = 1.0;
    };

    SS_ref_db.bounds_ref[0][0] = 0.0+eps;  SS_ref_db.bounds_ref[0][1] = 1.0-eps;
    SS_ref_db.bounds_ref[1][0] = 0.0+eps;  SS_ref_db.bounds_ref[1][1] = 1.0-eps;
    SS_ref_db.bounds_ref[2][0] = 0.0+eps;  SS_ref_db.bounds_ref[2][1] = 1.0-eps;
    SS_ref_db.bounds_ref[3][0] = 0.0+eps;  SS_ref_db.bounds_ref[3][1] = 1.0-eps;
    SS_ref_db.bounds_ref[4][0] = 0.0+eps;  SS_ref_db.bounds_ref[4][1] = 1.0-eps;


    return SS_ref_db;
}

/**
    Solution phase data for sb11_hpcpx
*/
SS_ref G_SS_sb11_hpcpx_function(SS_ref SS_ref_db, char* research_group, int EM_dataset, int len_ox, bulk_info z_b, double eps){

    int i, j;
    int n_em = SS_ref_db.n_em;

    char   *EM_tmp[]            = {"hpcen","hpcfs",};
    for (int i = 0; i < SS_ref_db.n_em; i++){
        strcpy(SS_ref_db.EM_list[i],EM_tmp[i]);
    };

    // Site mixing composition;
    SS_ref_db.C[0][0] = 0.0;    SS_ref_db.C[0][1] = 1.0;
    SS_ref_db.C[1][0] = 1.0;    SS_ref_db.C[1][1] = 0.0;

    // pre-computed Nullspace;
    SS_ref_db.N[0][0] = -0.7071067811865475;
    SS_ref_db.N[1][0] = 0.7071067811865476;

    SS_ref_db.W[0] = 0.0 ;


    em_data hpcen            = get_em_data(    research_group, EM_dataset,
                                            len_ox,
                                            z_b,
                                            SS_ref_db.P,
                                            SS_ref_db.T,
                                            "hpcen",
                                            "equilibrium");

    em_data hpcfs            = get_em_data(    research_group, EM_dataset,
                                            len_ox,
                                            z_b,
                                            SS_ref_db.P,
                                            SS_ref_db.T,
                                            "hpcfs",
                                            "equilibrium");

    SS_ref_db.gbase[0]        = hpcen.gb;
    SS_ref_db.gbase[1]        = hpcfs.gb;

    SS_ref_db.ElShearMod[0]        = hpcen.ElShearMod;
    SS_ref_db.ElShearMod[1]        = hpcfs.ElShearMod;

    SS_ref_db.ElBulkMod[0]        = hpcen.ElBulkMod;
    SS_ref_db.ElBulkMod[1]        = hpcfs.ElBulkMod;

    SS_ref_db.ElCp[0]        = hpcen.ElCp;
    SS_ref_db.ElCp[1]        = hpcfs.ElCp;

    SS_ref_db.ElExpansivity[0]        = hpcen.ElExpansivity;
    SS_ref_db.ElExpansivity[1]        = hpcfs.ElExpansivity;

    for (i = 0; i < len_ox; i++){
        SS_ref_db.Comp[0][i]    = hpcen.C[i];
        SS_ref_db.Comp[1][i]    = hpcfs.C[i];
    }

    for (i = 0; i < n_em; i++){
        SS_ref_db.z_em[i] = 1.0;
    };

    SS_ref_db.bounds_ref[0][0] = 0.0+eps;  SS_ref_db.bounds_ref[0][1] = 1.0-eps;
    SS_ref_db.bounds_ref[1][0] = 0.0+eps;  SS_ref_db.bounds_ref[1][1] = 1.0-eps;


    return SS_ref_db;
}

/**
    Solution phase data for sb11_ak
*/
SS_ref G_SS_sb11_ak_function(SS_ref SS_ref_db, char* research_group, int EM_dataset, int len_ox, bulk_info z_b, double eps){

    int i, j;
    int n_em = SS_ref_db.n_em;

    char   *EM_tmp[]            = {"mgak","feak","co",};
    for (int i = 0; i < SS_ref_db.n_em; i++){
        strcpy(SS_ref_db.EM_list[i],EM_tmp[i]);
    };

    // Site mixing composition;
    SS_ref_db.C[0][0] = 0.0;    SS_ref_db.C[0][1] = 0.0;    SS_ref_db.C[0][2] = 1.0;
    SS_ref_db.C[1][0] = 0.0;    SS_ref_db.C[1][1] = 1.0;    SS_ref_db.C[1][2] = 0.0;
    SS_ref_db.C[2][0] = 1.0;    SS_ref_db.C[2][1] = 0.0;    SS_ref_db.C[2][2] = 0.0;
    SS_ref_db.C[3][0] = 1.0;    SS_ref_db.C[3][1] = 1.0;    SS_ref_db.C[3][2] = 0.0;
    SS_ref_db.C[4][0] = 0.0;    SS_ref_db.C[4][1] = 0.0;    SS_ref_db.C[4][2] = 1.0;

    // pre-computed Nullspace;
    SS_ref_db.N[0][0] = -0.5773502691896257;    SS_ref_db.N[0][1] = -0.5773502691896257;
    SS_ref_db.N[1][0] = 0.7886751345948129;    SS_ref_db.N[1][1] = -0.2113248654051871;
    SS_ref_db.N[2][0] = -0.2113248654051871;    SS_ref_db.N[2][1] = 0.7886751345948129;

    SS_ref_db.W[0] = 0.0 ;
    SS_ref_db.W[1] = 66000.0 ;
    SS_ref_db.W[2] = 0.0 ;


    em_data mgak            = get_em_data(    research_group, EM_dataset,
                                            len_ox,
                                            z_b,
                                            SS_ref_db.P,
                                            SS_ref_db.T,
                                            "mgak",
                                            "equilibrium");

    em_data feak            = get_em_data(    research_group, EM_dataset,
                                            len_ox,
                                            z_b,
                                            SS_ref_db.P,
                                            SS_ref_db.T,
                                            "feak",
                                            "equilibrium");

    em_data co            = get_em_data(    research_group, EM_dataset,
                                            len_ox,
                                            z_b,
                                            SS_ref_db.P,
                                            SS_ref_db.T,
                                            "co",
                                            "equilibrium");

    SS_ref_db.gbase[0]        = mgak.gb;
    SS_ref_db.gbase[1]        = feak.gb;
    SS_ref_db.gbase[2]        = co.gb;

    SS_ref_db.ElShearMod[0]        = mgak.ElShearMod;
    SS_ref_db.ElShearMod[1]        = feak.ElShearMod;
    SS_ref_db.ElShearMod[2]        = co.ElShearMod;

    SS_ref_db.ElBulkMod[0]        = mgak.ElBulkMod;
    SS_ref_db.ElBulkMod[1]        = feak.ElBulkMod;
    SS_ref_db.ElBulkMod[2]        = co.ElBulkMod;

    SS_ref_db.ElCp[0]        = mgak.ElCp;
    SS_ref_db.ElCp[1]        = feak.ElCp;
    SS_ref_db.ElCp[2]        = co.ElCp;

    SS_ref_db.ElExpansivity[0]        = mgak.ElExpansivity;
    SS_ref_db.ElExpansivity[1]        = feak.ElExpansivity;
    SS_ref_db.ElExpansivity[2]        = co.ElExpansivity;

    for (i = 0; i < len_ox; i++){
        SS_ref_db.Comp[0][i]    = mgak.C[i];
        SS_ref_db.Comp[1][i]    = feak.C[i];
        SS_ref_db.Comp[2][i]    = co.C[i];
    }

    for (i = 0; i < n_em; i++){
        SS_ref_db.z_em[i] = 1.0;
    };

    SS_ref_db.bounds_ref[0][0] = 0.0+eps;  SS_ref_db.bounds_ref[0][1] = 1.0-eps;
    SS_ref_db.bounds_ref[1][0] = 0.0+eps;  SS_ref_db.bounds_ref[1][1] = 1.0-eps;
    SS_ref_db.bounds_ref[2][0] = 0.0+eps;  SS_ref_db.bounds_ref[2][1] = 1.0-eps;


    return SS_ref_db;
}

/**
    Solution phase data for sb11_gtmj
*/
SS_ref G_SS_sb11_gtmj_function(SS_ref SS_ref_db, char* research_group, int EM_dataset, int len_ox, bulk_info z_b, double eps){

    int i, j;
    int n_em = SS_ref_db.n_em;

    char   *EM_tmp[]            = {"py","alm","gr","mgmj","jdmj",};
    for (int i = 0; i < SS_ref_db.n_em; i++){
        strcpy(SS_ref_db.EM_list[i],EM_tmp[i]);
    };

    // Site mixing composition;
    SS_ref_db.C[0][0] = 0.0;    SS_ref_db.C[0][1] = 0.0;    SS_ref_db.C[0][2] = 1.0;    SS_ref_db.C[0][3] = 0.0;    SS_ref_db.C[0][4] = 0.0;
    SS_ref_db.C[1][0] = 0.0;    SS_ref_db.C[1][1] = 0.0;    SS_ref_db.C[1][2] = 0.0;    SS_ref_db.C[1][3] = 0.0;    SS_ref_db.C[1][4] = 0.3333333333333333;
    SS_ref_db.C[2][0] = 0.0;    SS_ref_db.C[2][1] = 1.0;    SS_ref_db.C[2][2] = 0.0;    SS_ref_db.C[2][3] = 0.0;    SS_ref_db.C[2][4] = 0.0;
    SS_ref_db.C[3][0] = 1.0;    SS_ref_db.C[3][1] = 0.0;    SS_ref_db.C[3][2] = 0.0;    SS_ref_db.C[3][3] = 1.0;    SS_ref_db.C[3][4] = 0.0;
    SS_ref_db.C[4][0] = 0.0;    SS_ref_db.C[4][1] = 0.0;    SS_ref_db.C[4][2] = 0.0;    SS_ref_db.C[4][3] = 0.0;    SS_ref_db.C[4][4] = 0.6666666666666666;
    SS_ref_db.C[5][0] = 1.0;    SS_ref_db.C[5][1] = 1.0;    SS_ref_db.C[5][2] = 1.0;    SS_ref_db.C[5][3] = 0.0;    SS_ref_db.C[5][4] = 1.0;
    SS_ref_db.C[6][0] = 0.0;    SS_ref_db.C[6][1] = 0.0;    SS_ref_db.C[6][2] = 0.0;    SS_ref_db.C[6][3] = 1.0;    SS_ref_db.C[6][4] = 0.0;
    SS_ref_db.C[7][0] = 0.0;    SS_ref_db.C[7][1] = 0.0;    SS_ref_db.C[7][2] = 0.0;    SS_ref_db.C[7][3] = 1.0;    SS_ref_db.C[7][4] = 1.0;
    SS_ref_db.C[8][0] = 1.0;    SS_ref_db.C[8][1] = 1.0;    SS_ref_db.C[8][2] = 1.0;    SS_ref_db.C[8][3] = 0.0;    SS_ref_db.C[8][4] = 0.0;

    // pre-computed Nullspace;
    SS_ref_db.N[0][0] = -0.44721359549995787;    SS_ref_db.N[0][1] = -0.44721359549995787;    SS_ref_db.N[0][2] = -0.44721359549995787;    SS_ref_db.N[0][3] = -0.44721359549995787;
    SS_ref_db.N[1][0] = 0.8618033988749895;    SS_ref_db.N[1][1] = -0.13819660112501048;    SS_ref_db.N[1][2] = -0.13819660112501048;    SS_ref_db.N[1][3] = -0.13819660112501048;
    SS_ref_db.N[2][0] = -0.13819660112501048;    SS_ref_db.N[2][1] = 0.8618033988749895;    SS_ref_db.N[2][2] = -0.13819660112501048;    SS_ref_db.N[2][3] = -0.13819660112501048;
    SS_ref_db.N[3][0] = -0.13819660112501048;    SS_ref_db.N[3][1] = -0.13819660112501048;    SS_ref_db.N[3][2] = 0.8618033988749895;    SS_ref_db.N[3][3] = -0.13819660112501048;
    SS_ref_db.N[4][0] = -0.13819660112501048;    SS_ref_db.N[4][1] = -0.13819660112501048;    SS_ref_db.N[4][2] = -0.13819660112501048;    SS_ref_db.N[4][3] = 0.8618033988749895;

    SS_ref_db.W[0] = 0.0 ;
    SS_ref_db.W[1] = 30000.0 ;
    SS_ref_db.W[2] = 21202.78 ;
    SS_ref_db.W[3] = 0.0 ;
    SS_ref_db.W[4] = 0.0 ;
    SS_ref_db.W[5] = 0.0 ;
    SS_ref_db.W[6] = 0.0 ;
    SS_ref_db.W[7] = 57775.96 ;
    SS_ref_db.W[8] = 0.0 ;
    SS_ref_db.W[9] = 0.0 ;


    em_data py            = get_em_data(    research_group, EM_dataset,
                                            len_ox,
                                            z_b,
                                            SS_ref_db.P,
                                            SS_ref_db.T,
                                            "py",
                                            "equilibrium");

    em_data alm            = get_em_data(    research_group, EM_dataset,
                                            len_ox,
                                            z_b,
                                            SS_ref_db.P,
                                            SS_ref_db.T,
                                            "alm",
                                            "equilibrium");

    em_data gr            = get_em_data(    research_group, EM_dataset,
                                            len_ox,
                                            z_b,
                                            SS_ref_db.P,
                                            SS_ref_db.T,
                                            "gr",
                                            "equilibrium");

    em_data mgmj            = get_em_data(    research_group, EM_dataset,
                                            len_ox,
                                            z_b,
                                            SS_ref_db.P,
                                            SS_ref_db.T,
                                            "mgmj",
                                            "equilibrium");

    em_data jdmj            = get_em_data(    research_group, EM_dataset,
                                            len_ox,
                                            z_b,
                                            SS_ref_db.P,
                                            SS_ref_db.T,
                                            "jdmj",
                                            "equilibrium");

    SS_ref_db.gbase[0]        = py.gb;
    SS_ref_db.gbase[1]        = alm.gb;
    SS_ref_db.gbase[2]        = gr.gb;
    SS_ref_db.gbase[3]        = mgmj.gb;
    SS_ref_db.gbase[4]        = jdmj.gb;

    SS_ref_db.ElShearMod[0]        = py.ElShearMod;
    SS_ref_db.ElShearMod[1]        = alm.ElShearMod;
    SS_ref_db.ElShearMod[2]        = gr.ElShearMod;
    SS_ref_db.ElShearMod[3]        = mgmj.ElShearMod;
    SS_ref_db.ElShearMod[4]        = jdmj.ElShearMod;

    SS_ref_db.ElBulkMod[0]        = py.ElBulkMod;
    SS_ref_db.ElBulkMod[1]        = alm.ElBulkMod;
    SS_ref_db.ElBulkMod[2]        = gr.ElBulkMod;
    SS_ref_db.ElBulkMod[3]        = mgmj.ElBulkMod;
    SS_ref_db.ElBulkMod[4]        = jdmj.ElBulkMod;

    SS_ref_db.ElCp[0]        = py.ElCp;
    SS_ref_db.ElCp[1]        = alm.ElCp;
    SS_ref_db.ElCp[2]        = gr.ElCp;
    SS_ref_db.ElCp[3]        = mgmj.ElCp;
    SS_ref_db.ElCp[4]        = jdmj.ElCp;

    SS_ref_db.ElExpansivity[0]        = py.ElExpansivity;
    SS_ref_db.ElExpansivity[1]        = alm.ElExpansivity;
    SS_ref_db.ElExpansivity[2]        = gr.ElExpansivity;
    SS_ref_db.ElExpansivity[3]        = mgmj.ElExpansivity;
    SS_ref_db.ElExpansivity[4]        = jdmj.ElExpansivity;

    for (i = 0; i < len_ox; i++){
        SS_ref_db.Comp[0][i]    = py.C[i];
        SS_ref_db.Comp[1][i]    = alm.C[i];
        SS_ref_db.Comp[2][i]    = gr.C[i];
        SS_ref_db.Comp[3][i]    = mgmj.C[i];
        SS_ref_db.Comp[4][i]    = jdmj.C[i];
    }

    for (i = 0; i < n_em; i++){
        SS_ref_db.z_em[i] = 1.0;
    };

    SS_ref_db.bounds_ref[0][0] = 0.0+eps;  SS_ref_db.bounds_ref[0][1] = 1.0-eps;
    SS_ref_db.bounds_ref[1][0] = 0.0+eps;  SS_ref_db.bounds_ref[1][1] = 1.0-eps;
    SS_ref_db.bounds_ref[2][0] = 0.0+eps;  SS_ref_db.bounds_ref[2][1] = 1.0-eps;
    SS_ref_db.bounds_ref[3][0] = 0.0+eps;  SS_ref_db.bounds_ref[3][1] = 1.0-eps;
    SS_ref_db.bounds_ref[4][0] = 0.0+eps;  SS_ref_db.bounds_ref[4][1] = 1.0-eps;


    return SS_ref_db;
}

/**
    Solution phase data for sb11_pv
*/
SS_ref G_SS_sb11_pv_function(SS_ref SS_ref_db, char* research_group, int EM_dataset, int len_ox, bulk_info z_b, double eps){

    int i, j;
    int n_em = SS_ref_db.n_em;

    char   *EM_tmp[]            = {"mgpv","fepv","alpv",};
    for (int i = 0; i < SS_ref_db.n_em; i++){
        strcpy(SS_ref_db.EM_list[i],EM_tmp[i]);
    };

    // Site mixing composition;
    SS_ref_db.C[0][0] = 0.0;    SS_ref_db.C[0][1] = 0.0;    SS_ref_db.C[0][2] = 1.0;
    SS_ref_db.C[1][0] = 0.0;    SS_ref_db.C[1][1] = 1.0;    SS_ref_db.C[1][2] = 0.0;
    SS_ref_db.C[2][0] = 1.0;    SS_ref_db.C[2][1] = 0.0;    SS_ref_db.C[2][2] = 0.0;
    SS_ref_db.C[3][0] = 1.0;    SS_ref_db.C[3][1] = 1.0;    SS_ref_db.C[3][2] = 0.0;
    SS_ref_db.C[4][0] = 0.0;    SS_ref_db.C[4][1] = 0.0;    SS_ref_db.C[4][2] = 1.0;

    // pre-computed Nullspace;
    SS_ref_db.N[0][0] = -0.5773502691896257;    SS_ref_db.N[0][1] = -0.5773502691896257;
    SS_ref_db.N[1][0] = 0.7886751345948129;    SS_ref_db.N[1][1] = -0.2113248654051871;
    SS_ref_db.N[2][0] = -0.2113248654051871;    SS_ref_db.N[2][1] = 0.7886751345948129;

    SS_ref_db.W[0] = 0.0 ;
    SS_ref_db.W[1] = 116000.0 ;
    SS_ref_db.W[2] = 0.0 ;

    SS_ref_db.v[0] = 1.0;
    SS_ref_db.v[1] = 1.0;
    SS_ref_db.v[2] = 0.39;

    em_data mgpv            = get_em_data(    research_group, EM_dataset,
                                            len_ox,
                                            z_b,
                                            SS_ref_db.P,
                                            SS_ref_db.T,
                                            "mgpv",
                                            "equilibrium");

    em_data fepv            = get_em_data(    research_group, EM_dataset,
                                            len_ox,
                                            z_b,
                                            SS_ref_db.P,
                                            SS_ref_db.T,
                                            "fepv",
                                            "equilibrium");

    em_data alpv            = get_em_data(    research_group, EM_dataset,
                                            len_ox,
                                            z_b,
                                            SS_ref_db.P,
                                            SS_ref_db.T,
                                            "alpv",
                                            "equilibrium");

    SS_ref_db.gbase[0]        = mgpv.gb;
    SS_ref_db.gbase[1]        = fepv.gb;
    SS_ref_db.gbase[2]        = alpv.gb;

    SS_ref_db.ElShearMod[0]        = mgpv.ElShearMod;
    SS_ref_db.ElShearMod[1]        = fepv.ElShearMod;
    SS_ref_db.ElShearMod[2]        = alpv.ElShearMod;

    SS_ref_db.ElBulkMod[0]        = mgpv.ElBulkMod;
    SS_ref_db.ElBulkMod[1]        = fepv.ElBulkMod;
    SS_ref_db.ElBulkMod[2]        = alpv.ElBulkMod;

    SS_ref_db.ElCp[0]        = mgpv.ElCp;
    SS_ref_db.ElCp[1]        = fepv.ElCp;
    SS_ref_db.ElCp[2]        = alpv.ElCp;

    SS_ref_db.ElExpansivity[0]        = mgpv.ElExpansivity;
    SS_ref_db.ElExpansivity[1]        = fepv.ElExpansivity;
    SS_ref_db.ElExpansivity[2]        = alpv.ElExpansivity;

    for (i = 0; i < len_ox; i++){
        SS_ref_db.Comp[0][i]    = mgpv.C[i];
        SS_ref_db.Comp[1][i]    = fepv.C[i];
        SS_ref_db.Comp[2][i]    = alpv.C[i];
    }

    for (i = 0; i < n_em; i++){
        SS_ref_db.z_em[i] = 1.0;
    };

    SS_ref_db.bounds_ref[0][0] = 0.0+eps;  SS_ref_db.bounds_ref[0][1] = 1.0-eps;
    SS_ref_db.bounds_ref[1][0] = 0.0+eps;  SS_ref_db.bounds_ref[1][1] = 1.0-eps;
    SS_ref_db.bounds_ref[2][0] = 0.0+eps;  SS_ref_db.bounds_ref[2][1] = 1.0-eps;


    return SS_ref_db;
}

/**
    Solution phase data for sb11_ppv
*/
SS_ref G_SS_sb11_ppv_function(SS_ref SS_ref_db, char* research_group, int EM_dataset, int len_ox, bulk_info z_b, double eps){

    int i, j;
    int n_em = SS_ref_db.n_em;

    char   *EM_tmp[]            = {"mppv","fppv","appv",};
    for (int i = 0; i < SS_ref_db.n_em; i++){
        strcpy(SS_ref_db.EM_list[i],EM_tmp[i]);
    };

    // Site mixing composition;
    SS_ref_db.C[0][0] = 0.0;    SS_ref_db.C[0][1] = 0.0;    SS_ref_db.C[0][2] = 1.0;
    SS_ref_db.C[1][0] = 0.0;    SS_ref_db.C[1][1] = 1.0;    SS_ref_db.C[1][2] = 0.0;
    SS_ref_db.C[2][0] = 1.0;    SS_ref_db.C[2][1] = 0.0;    SS_ref_db.C[2][2] = 0.0;
    SS_ref_db.C[3][0] = 1.0;    SS_ref_db.C[3][1] = 1.0;    SS_ref_db.C[3][2] = 0.0;
    SS_ref_db.C[4][0] = 0.0;    SS_ref_db.C[4][1] = 0.0;    SS_ref_db.C[4][2] = 1.0;

    // pre-computed Nullspace;
    SS_ref_db.N[0][0] = -0.5773502691896257;    SS_ref_db.N[0][1] = -0.5773502691896257;
    SS_ref_db.N[1][0] = 0.7886751345948129;    SS_ref_db.N[1][1] = -0.2113248654051871;
    SS_ref_db.N[2][0] = -0.2113248654051871;    SS_ref_db.N[2][1] = 0.7886751345948129;

    SS_ref_db.W[0] = 0.0 ;
    SS_ref_db.W[1] = 60000.0 ;
    SS_ref_db.W[2] = 0.0 ;


    em_data mppv            = get_em_data(    research_group, EM_dataset,
                                            len_ox,
                                            z_b,
                                            SS_ref_db.P,
                                            SS_ref_db.T,
                                            "mppv",
                                            "equilibrium");

    em_data fppv            = get_em_data(    research_group, EM_dataset,
                                            len_ox,
                                            z_b,
                                            SS_ref_db.P,
                                            SS_ref_db.T,
                                            "fppv",
                                            "equilibrium");

    em_data appv            = get_em_data(    research_group, EM_dataset,
                                            len_ox,
                                            z_b,
                                            SS_ref_db.P,
                                            SS_ref_db.T,
                                            "appv",
                                            "equilibrium");

    SS_ref_db.gbase[0]        = mppv.gb;
    SS_ref_db.gbase[1]        = fppv.gb;
    SS_ref_db.gbase[2]        = appv.gb;

    SS_ref_db.ElShearMod[0]        = mppv.ElShearMod;
    SS_ref_db.ElShearMod[1]        = fppv.ElShearMod;
    SS_ref_db.ElShearMod[2]        = appv.ElShearMod;

    SS_ref_db.ElBulkMod[0]        = mppv.ElBulkMod;
    SS_ref_db.ElBulkMod[1]        = fppv.ElBulkMod;
    SS_ref_db.ElBulkMod[2]        = appv.ElBulkMod;

    SS_ref_db.ElCp[0]        = mppv.ElCp;
    SS_ref_db.ElCp[1]        = fppv.ElCp;
    SS_ref_db.ElCp[2]        = appv.ElCp;

    SS_ref_db.ElExpansivity[0]        = mppv.ElExpansivity;
    SS_ref_db.ElExpansivity[1]        = fppv.ElExpansivity;
    SS_ref_db.ElExpansivity[2]        = appv.ElExpansivity;

    for (i = 0; i < len_ox; i++){
        SS_ref_db.Comp[0][i]    = mppv.C[i];
        SS_ref_db.Comp[1][i]    = fppv.C[i];
        SS_ref_db.Comp[2][i]    = appv.C[i];
    }

    for (i = 0; i < n_em; i++){
        SS_ref_db.z_em[i] = 1.0;
    };

    SS_ref_db.bounds_ref[0][0] = 0.0+eps;  SS_ref_db.bounds_ref[0][1] = 1.0-eps;
    SS_ref_db.bounds_ref[1][0] = 0.0+eps;  SS_ref_db.bounds_ref[1][1] = 1.0-eps;
    SS_ref_db.bounds_ref[2][0] = 0.0+eps;  SS_ref_db.bounds_ref[2][1] = 1.0-eps;


    return SS_ref_db;
}

/**
    Solution phase data for sb11_mw
*/
SS_ref G_SS_sb11_mw_function(SS_ref SS_ref_db, char* research_group, int EM_dataset, int len_ox, bulk_info z_b, double eps){

    int i, j;
    int n_em = SS_ref_db.n_em;

    char   *EM_tmp[]            = {"pe","wu",};
    for (int i = 0; i < SS_ref_db.n_em; i++){
        strcpy(SS_ref_db.EM_list[i],EM_tmp[i]);
    };

    // Site mixing composition;
    SS_ref_db.C[0][0] = 0.0;    SS_ref_db.C[0][1] = 1.0;
    SS_ref_db.C[1][0] = 1.0;    SS_ref_db.C[1][1] = 0.0;

    // pre-computed Nullspace;
    SS_ref_db.N[0][0] = -0.7071067811865475;
    SS_ref_db.N[1][0] = 0.7071067811865476;

    SS_ref_db.W[0] = 13000.0 ;


    em_data pe            = get_em_data(    research_group, EM_dataset,
                                            len_ox,
                                            z_b,
                                            SS_ref_db.P,
                                            SS_ref_db.T,
                                            "pe",
                                            "equilibrium");

    em_data wu            = get_em_data(    research_group, EM_dataset,
                                            len_ox,
                                            z_b,
                                            SS_ref_db.P,
                                            SS_ref_db.T,
                                            "wu",
                                            "equilibrium");

    SS_ref_db.gbase[0]        = pe.gb;
    SS_ref_db.gbase[1]        = wu.gb;

    SS_ref_db.ElShearMod[0]        = pe.ElShearMod;
    SS_ref_db.ElShearMod[1]        = wu.ElShearMod;

    SS_ref_db.ElBulkMod[0]        = pe.ElBulkMod;
    SS_ref_db.ElBulkMod[1]        = wu.ElBulkMod;

    SS_ref_db.ElCp[0]        = pe.ElCp;
    SS_ref_db.ElCp[1]        = wu.ElCp;

    SS_ref_db.ElExpansivity[0]        = pe.ElExpansivity;
    SS_ref_db.ElExpansivity[1]        = wu.ElExpansivity;

    for (i = 0; i < len_ox; i++){
        SS_ref_db.Comp[0][i]    = pe.C[i];
        SS_ref_db.Comp[1][i]    = wu.C[i];
    }

    for (i = 0; i < n_em; i++){
        SS_ref_db.z_em[i] = 1.0;
    };

    SS_ref_db.bounds_ref[0][0] = 0.0+eps;  SS_ref_db.bounds_ref[0][1] = 1.0-eps;
    SS_ref_db.bounds_ref[1][0] = 0.0+eps;  SS_ref_db.bounds_ref[1][1] = 1.0-eps;


    return SS_ref_db;
}

/**
    Solution phase data for sb11_cf
*/
SS_ref G_SS_sb11_cf_function(SS_ref SS_ref_db, char* research_group, int EM_dataset, int len_ox, bulk_info z_b, double eps){

    int i, j;
    int n_em = SS_ref_db.n_em;

    char   *EM_tmp[]            = {"mgcf","fecf","nacf",};
    for (int i = 0; i < SS_ref_db.n_em; i++){
        strcpy(SS_ref_db.EM_list[i],EM_tmp[i]);
    };

    // Site mixing composition;
    SS_ref_db.C[0][0] = 0.0;    SS_ref_db.C[0][1] = 1.0;    SS_ref_db.C[0][2] = 0.0;
    SS_ref_db.C[1][0] = 1.0;    SS_ref_db.C[1][1] = 0.0;    SS_ref_db.C[1][2] = 0.0;
    SS_ref_db.C[2][0] = 0.0;    SS_ref_db.C[2][1] = 0.0;    SS_ref_db.C[2][2] = 1.0;
    SS_ref_db.C[3][0] = 0.0;    SS_ref_db.C[3][1] = 0.0;    SS_ref_db.C[3][2] = 1.0;
    SS_ref_db.C[4][0] = 1.0;    SS_ref_db.C[4][1] = 1.0;    SS_ref_db.C[4][2] = 0.0;

    // pre-computed Nullspace;
    SS_ref_db.N[0][0] = -0.5773502691896257;    SS_ref_db.N[0][1] = -0.5773502691896257;
    SS_ref_db.N[1][0] = 0.7886751345948129;    SS_ref_db.N[1][1] = -0.2113248654051871;
    SS_ref_db.N[2][0] = -0.2113248654051871;    SS_ref_db.N[2][1] = 0.7886751345948129;

    SS_ref_db.W[0] = 0.0 ;
    SS_ref_db.W[1] = 0.0 ;
    SS_ref_db.W[2] = 0.0 ;


    em_data mgcf            = get_em_data(    research_group, EM_dataset,
                                            len_ox,
                                            z_b,
                                            SS_ref_db.P,
                                            SS_ref_db.T,
                                            "mgcf",
                                            "equilibrium");

    em_data fecf            = get_em_data(    research_group, EM_dataset,
                                            len_ox,
                                            z_b,
                                            SS_ref_db.P,
                                            SS_ref_db.T,
                                            "fecf",
                                            "equilibrium");

    em_data nacf            = get_em_data(    research_group, EM_dataset,
                                            len_ox,
                                            z_b,
                                            SS_ref_db.P,
                                            SS_ref_db.T,
                                            "nacf",
                                            "equilibrium");

    SS_ref_db.gbase[0]        = mgcf.gb;
    SS_ref_db.gbase[1]        = fecf.gb;
    SS_ref_db.gbase[2]        = nacf.gb;

    SS_ref_db.ElShearMod[0]        = mgcf.ElShearMod;
    SS_ref_db.ElShearMod[1]        = fecf.ElShearMod;
    SS_ref_db.ElShearMod[2]        = nacf.ElShearMod;

    SS_ref_db.ElBulkMod[0]        = mgcf.ElBulkMod;
    SS_ref_db.ElBulkMod[1]        = fecf.ElBulkMod;
    SS_ref_db.ElBulkMod[2]        = nacf.ElBulkMod;

    SS_ref_db.ElCp[0]        = mgcf.ElCp;
    SS_ref_db.ElCp[1]        = fecf.ElCp;
    SS_ref_db.ElCp[2]        = nacf.ElCp;

    SS_ref_db.ElExpansivity[0]        = mgcf.ElExpansivity;
    SS_ref_db.ElExpansivity[1]        = fecf.ElExpansivity;
    SS_ref_db.ElExpansivity[2]        = nacf.ElExpansivity;

    for (i = 0; i < len_ox; i++){
        SS_ref_db.Comp[0][i]    = mgcf.C[i];
        SS_ref_db.Comp[1][i]    = fecf.C[i];
        SS_ref_db.Comp[2][i]    = nacf.C[i];
    }

    for (i = 0; i < n_em; i++){
        SS_ref_db.z_em[i] = 1.0;
    };

    SS_ref_db.bounds_ref[0][0] = 0.0+eps;  SS_ref_db.bounds_ref[0][1] = 1.0-eps;
    SS_ref_db.bounds_ref[1][0] = 0.0+eps;  SS_ref_db.bounds_ref[1][1] = 1.0-eps;
    SS_ref_db.bounds_ref[2][0] = 0.0+eps;  SS_ref_db.bounds_ref[2][1] = 1.0-eps;


    return SS_ref_db;
}

SS_ref G_SS_sb11_EM_function(       global_variable        gv,
                                    SS_ref            SS_ref_db,
                                    int            EM_dataset,
                                    bulk_info            z_b,
                                    char            *name){            

    double eps 		   	= gv.bnd_val;
    double P 			= SS_ref_db.P;
    double T 			= SS_ref_db.T;

    SS_ref_db.ss_flags[0]  = 1;

    /* Associate the right solid-solution data */
    for (int FD = 0; FD < gv.n_Diff; FD++){	

        if (FD == 8 || FD == 9){
            SS_ref_db.P = 1.+ gv.gb_P_eps*gv.pdev[0][FD];
            SS_ref_db.T = T + gv.gb_T_eps*gv.pdev[1][FD];
        }
        else{
            SS_ref_db.P = P + gv.gb_P_eps*gv.pdev[0][FD];
            SS_ref_db.T = T + gv.gb_T_eps*gv.pdev[1][FD];
        }

        if (strcmp( name, "plg") == 0 ){
  			if ( P > 500.0){
				SS_ref_db.ss_flags[0]  = 0;
			}
            SS_ref_db  = G_SS_sb11_plg_function(SS_ref_db, gv.research_group, EM_dataset, gv.len_ox, z_b, eps);	}
        else if (strcmp( name, "sp") == 0 ){
            SS_ref_db  = G_SS_sb11_sp_function(SS_ref_db, gv.research_group, EM_dataset, gv.len_ox, z_b, eps);	}
        else if (strcmp( name, "ol") == 0 ){
            SS_ref_db  = G_SS_sb11_ol_function(SS_ref_db, gv.research_group, EM_dataset, gv.len_ox, z_b, eps);	}
        else if (strcmp( name, "wa") == 0 ){
            SS_ref_db  = G_SS_sb11_wa_function(SS_ref_db, gv.research_group, EM_dataset, gv.len_ox, z_b, eps);	}
        else if (strcmp( name, "ri") == 0 ){
            SS_ref_db  = G_SS_sb11_ri_function(SS_ref_db, gv.research_group, EM_dataset, gv.len_ox, z_b, eps);	}
        else if (strcmp( name, "opx") == 0 ){
            SS_ref_db  = G_SS_sb11_opx_function(SS_ref_db, gv.research_group, EM_dataset, gv.len_ox, z_b, eps);	}
        else if (strcmp( name, "cpx") == 0 ){
            SS_ref_db  = G_SS_sb11_cpx_function(SS_ref_db, gv.research_group, EM_dataset, gv.len_ox, z_b, eps);	}
        else if (strcmp( name, "hpcpx") == 0 ){
            SS_ref_db  = G_SS_sb11_hpcpx_function(SS_ref_db, gv.research_group, EM_dataset, gv.len_ox, z_b, eps);	}
        else if (strcmp( name, "ak") == 0 ){
            SS_ref_db  = G_SS_sb11_ak_function(SS_ref_db, gv.research_group, EM_dataset, gv.len_ox, z_b, eps);	}
        else if (strcmp( name, "gtmj") == 0 ){
            SS_ref_db  = G_SS_sb11_gtmj_function(SS_ref_db, gv.research_group, EM_dataset, gv.len_ox, z_b, eps);	}
        else if (strcmp( name, "pv") == 0 ){
            SS_ref_db  = G_SS_sb11_pv_function(SS_ref_db, gv.research_group, EM_dataset, gv.len_ox, z_b, eps);	}
        else if (strcmp( name, "ppv") == 0 ){
            SS_ref_db  = G_SS_sb11_ppv_function(SS_ref_db, gv.research_group, EM_dataset, gv.len_ox, z_b, eps);	}
        else if (strcmp( name, "mw") == 0 ){
            SS_ref_db  = G_SS_sb11_mw_function(SS_ref_db, gv.research_group, EM_dataset, gv.len_ox, z_b, eps);	}
        else if (strcmp( name, "cf") == 0 ){
            SS_ref_db  = G_SS_sb11_cf_function(SS_ref_db, gv.research_group, EM_dataset, gv.len_ox, z_b, eps);	}
        else{
            printf("\nsolid solution '%s' is not in the database\n",name);	}

        for (int j = 0; j < SS_ref_db.n_em; j++){
            SS_ref_db.mu_array[FD][j] = SS_ref_db.gbase[j];
        }
    }

    for (int j = 0; j < SS_ref_db.n_em; j++){
        SS_ref_db.bounds[j][0] = SS_ref_db.bounds_ref[j][0];
        SS_ref_db.bounds[j][1] = SS_ref_db.bounds_ref[j][1];
    }

    /* Calculate the number of atoms in the bulk-rock composition */
    double fbc     = 0.0;
    for (int i = 0; i < gv.len_ox; i++){
        fbc += z_b.bulk_rock[i]*z_b.apo[i];
    }

    for (int i = 0; i < SS_ref_db.n_em; i++){
        SS_ref_db.ape[i] = 0.0;
        for (int j = 0; j < gv.len_ox; j++){
            SS_ref_db.ape[i] += SS_ref_db.Comp[i][j]*z_b.apo[j];
        }
    }

    SS_ref_db.fbc = z_b.fbc;

    if (gv.verbose == 1){
        printf(" %4s:",name);
        for (int j = 0; j < SS_ref_db.n_em; j++){
            printf(" %+12.5f",SS_ref_db.gbase[j]);
        }
        printf("\n");
        printf(" S   C   A   F   M   N\n");
        for (int i = 0; i < SS_ref_db.n_em; i++){
            for (int j = 0; j < gv.len_ox; j++){
                printf(" %.1f",SS_ref_db.Comp[i][j]);
            }
            printf("\n");
        }
        printf("\n");
    }

    return SS_ref_db;
};


/**
    Solution phase data for sb21_plg
*/
SS_ref G_SS_sb21_plg_function(SS_ref SS_ref_db, char* research_group, int EM_dataset, int len_ox, bulk_info z_b, double eps){

    int i, j;
    int n_em = SS_ref_db.n_em;

    char   *EM_tmp[]            = {"an","ab",};
    for (int i = 0; i < SS_ref_db.n_em; i++){
        strcpy(SS_ref_db.EM_list[i],EM_tmp[i]);
    };

    // Site mixing composition;
    SS_ref_db.C[0][0] = 1.0;    SS_ref_db.C[0][1] = 0.0;
    SS_ref_db.C[1][0] = 0.0;    SS_ref_db.C[1][1] = 1.0;

    // pre-computed Nullspace;
    SS_ref_db.N[0][0] = -0.7071067811865475;
    SS_ref_db.N[1][0] = 0.7071067811865476;

    SS_ref_db.W[0] = 13000.0 ;


    em_data an            = get_em_data(    research_group, EM_dataset,
                                            len_ox,
                                            z_b,
                                            SS_ref_db.P,
                                            SS_ref_db.T,
                                            "an",
                                            "equilibrium");

    em_data ab            = get_em_data(    research_group, EM_dataset,
                                            len_ox,
                                            z_b,
                                            SS_ref_db.P,
                                            SS_ref_db.T,
                                            "ab",
                                            "equilibrium");

    SS_ref_db.gbase[0]        = an.gb;
    SS_ref_db.gbase[1]        = ab.gb;

    SS_ref_db.ElShearMod[0]        = an.ElShearMod;
    SS_ref_db.ElShearMod[1]        = ab.ElShearMod;

    SS_ref_db.ElBulkMod[0]        = an.ElBulkMod;
    SS_ref_db.ElBulkMod[1]        = ab.ElBulkMod;

    SS_ref_db.ElCp[0]        = an.ElCp;
    SS_ref_db.ElCp[1]        = ab.ElCp;

    SS_ref_db.ElExpansivity[0]        = an.ElExpansivity;
    SS_ref_db.ElExpansivity[1]        = ab.ElExpansivity;

    for (i = 0; i < len_ox; i++){
        SS_ref_db.Comp[0][i]    = an.C[i];
        SS_ref_db.Comp[1][i]    = ab.C[i];
    }

    for (i = 0; i < n_em; i++){
        SS_ref_db.z_em[i] = 1.0;
    };

    SS_ref_db.bounds_ref[0][0] = 0.0+eps;  SS_ref_db.bounds_ref[0][1] = 1.0-eps;
    SS_ref_db.bounds_ref[1][0] = 0.0+eps;  SS_ref_db.bounds_ref[1][1] = 1.0-eps;


    return SS_ref_db;
}

/**
    Solution phase data for sb21_sp
*/
SS_ref G_SS_sb21_sp_function(SS_ref SS_ref_db, char* research_group, int EM_dataset, int len_ox, bulk_info z_b, double eps){

    int i, j;
    int n_em = SS_ref_db.n_em;

    char   *EM_tmp[]            = {"sp","hc",};
    for (int i = 0; i < SS_ref_db.n_em; i++){
        strcpy(SS_ref_db.EM_list[i],EM_tmp[i]);
    };

    // Site mixing composition;
    SS_ref_db.C[0][0] = 0.25;    SS_ref_db.C[0][1] = 0.25;
    SS_ref_db.C[1][0] = 0.0;    SS_ref_db.C[1][1] = 0.75;
    SS_ref_db.C[2][0] = 0.75;    SS_ref_db.C[2][1] = 0.0;
    SS_ref_db.C[3][0] = 0.875;    SS_ref_db.C[3][1] = 0.875;
    SS_ref_db.C[4][0] = 0.0;    SS_ref_db.C[4][1] = 0.125;
    SS_ref_db.C[5][0] = 0.125;    SS_ref_db.C[5][1] = 0.0;

    // pre-computed Nullspace;
    SS_ref_db.N[0][0] = -0.7071067811865475;
    SS_ref_db.N[1][0] = 0.7071067811865476;

    SS_ref_db.W[0] = -533.21 ;


    em_data sp            = get_em_data(    research_group, EM_dataset,
                                            len_ox,
                                            z_b,
                                            SS_ref_db.P,
                                            SS_ref_db.T,
                                            "sp",
                                            "equilibrium");

    em_data hc            = get_em_data(    research_group, EM_dataset,
                                            len_ox,
                                            z_b,
                                            SS_ref_db.P,
                                            SS_ref_db.T,
                                            "hc",
                                            "equilibrium");

    SS_ref_db.gbase[0]        = sp.gb;
    SS_ref_db.gbase[1]        = hc.gb;

    SS_ref_db.ElShearMod[0]        = sp.ElShearMod;
    SS_ref_db.ElShearMod[1]        = hc.ElShearMod;

    SS_ref_db.ElBulkMod[0]        = sp.ElBulkMod;
    SS_ref_db.ElBulkMod[1]        = hc.ElBulkMod;

    SS_ref_db.ElCp[0]        = sp.ElCp;
    SS_ref_db.ElCp[1]        = hc.ElCp;

    SS_ref_db.ElExpansivity[0]        = sp.ElExpansivity;
    SS_ref_db.ElExpansivity[1]        = hc.ElExpansivity;

    for (i = 0; i < len_ox; i++){
        SS_ref_db.Comp[0][i]    = sp.C[i];
        SS_ref_db.Comp[1][i]    = hc.C[i];
    }

    for (i = 0; i < n_em; i++){
        SS_ref_db.z_em[i] = 1.0;
    };

    SS_ref_db.bounds_ref[0][0] = 0.0+eps;  SS_ref_db.bounds_ref[0][1] = 1.0-eps;
    SS_ref_db.bounds_ref[1][0] = 0.0+eps;  SS_ref_db.bounds_ref[1][1] = 1.0-eps;


    return SS_ref_db;
}

/**
    Solution phase data for sb21_ol
*/
SS_ref G_SS_sb21_ol_function(SS_ref SS_ref_db, char* research_group, int EM_dataset, int len_ox, bulk_info z_b, double eps){

    int i, j;
    int n_em = SS_ref_db.n_em;

    char   *EM_tmp[]            = {"fo","fa",};
    for (int i = 0; i < SS_ref_db.n_em; i++){
        strcpy(SS_ref_db.EM_list[i],EM_tmp[i]);
    };

    // Site mixing composition;
    SS_ref_db.C[0][0] = 0.0;    SS_ref_db.C[0][1] = 1.0;
    SS_ref_db.C[1][0] = 1.0;    SS_ref_db.C[1][1] = 0.0;

    // pre-computed Nullspace;
    SS_ref_db.N[0][0] = -0.7071067811865475;
    SS_ref_db.N[1][0] = 0.7071067811865476;

    SS_ref_db.W[0] = 4694.66 ;


    em_data fo            = get_em_data(    research_group, EM_dataset,
                                            len_ox,
                                            z_b,
                                            SS_ref_db.P,
                                            SS_ref_db.T,
                                            "fo",
                                            "equilibrium");

    em_data fa            = get_em_data(    research_group, EM_dataset,
                                            len_ox,
                                            z_b,
                                            SS_ref_db.P,
                                            SS_ref_db.T,
                                            "fa",
                                            "equilibrium");

    SS_ref_db.gbase[0]        = fo.gb;
    SS_ref_db.gbase[1]        = fa.gb;

    SS_ref_db.ElShearMod[0]        = fo.ElShearMod;
    SS_ref_db.ElShearMod[1]        = fa.ElShearMod;

    SS_ref_db.ElBulkMod[0]        = fo.ElBulkMod;
    SS_ref_db.ElBulkMod[1]        = fa.ElBulkMod;

    SS_ref_db.ElCp[0]        = fo.ElCp;
    SS_ref_db.ElCp[1]        = fa.ElCp;

    SS_ref_db.ElExpansivity[0]        = fo.ElExpansivity;
    SS_ref_db.ElExpansivity[1]        = fa.ElExpansivity;

    for (i = 0; i < len_ox; i++){
        SS_ref_db.Comp[0][i]    = fo.C[i];
        SS_ref_db.Comp[1][i]    = fa.C[i];
    }

    for (i = 0; i < n_em; i++){
        SS_ref_db.z_em[i] = 1.0;
    };

    SS_ref_db.bounds_ref[0][0] = 0.0+eps;  SS_ref_db.bounds_ref[0][1] = 1.0-eps;
    SS_ref_db.bounds_ref[1][0] = 0.0+eps;  SS_ref_db.bounds_ref[1][1] = 1.0-eps;


    return SS_ref_db;
}

/**
    Solution phase data for sb21_wa
*/
SS_ref G_SS_sb21_wa_function(SS_ref SS_ref_db, char* research_group, int EM_dataset, int len_ox, bulk_info z_b, double eps){

    int i, j;
    int n_em = SS_ref_db.n_em;

    char   *EM_tmp[]            = {"mgwa","fewa",};
    for (int i = 0; i < SS_ref_db.n_em; i++){
        strcpy(SS_ref_db.EM_list[i],EM_tmp[i]);
    };

    // Site mixing composition;
    SS_ref_db.C[0][0] = 0.0;    SS_ref_db.C[0][1] = 1.0;
    SS_ref_db.C[1][0] = 1.0;    SS_ref_db.C[1][1] = 0.0;

    // pre-computed Nullspace;
    SS_ref_db.N[0][0] = -0.7071067811865475;
    SS_ref_db.N[1][0] = 0.7071067811865476;

    SS_ref_db.W[0] = 13202.38 ;


    em_data mgwa            = get_em_data(    research_group, EM_dataset,
                                            len_ox,
                                            z_b,
                                            SS_ref_db.P,
                                            SS_ref_db.T,
                                            "mgwa",
                                            "equilibrium");

    em_data fewa            = get_em_data(    research_group, EM_dataset,
                                            len_ox,
                                            z_b,
                                            SS_ref_db.P,
                                            SS_ref_db.T,
                                            "fewa",
                                            "equilibrium");

    SS_ref_db.gbase[0]        = mgwa.gb;
    SS_ref_db.gbase[1]        = fewa.gb;

    SS_ref_db.ElShearMod[0]        = mgwa.ElShearMod;
    SS_ref_db.ElShearMod[1]        = fewa.ElShearMod;

    SS_ref_db.ElBulkMod[0]        = mgwa.ElBulkMod;
    SS_ref_db.ElBulkMod[1]        = fewa.ElBulkMod;

    SS_ref_db.ElCp[0]        = mgwa.ElCp;
    SS_ref_db.ElCp[1]        = fewa.ElCp;

    SS_ref_db.ElExpansivity[0]        = mgwa.ElExpansivity;
    SS_ref_db.ElExpansivity[1]        = fewa.ElExpansivity;

    for (i = 0; i < len_ox; i++){
        SS_ref_db.Comp[0][i]    = mgwa.C[i];
        SS_ref_db.Comp[1][i]    = fewa.C[i];
    }

    for (i = 0; i < n_em; i++){
        SS_ref_db.z_em[i] = 1.0;
    };

    SS_ref_db.bounds_ref[0][0] = 0.0+eps;  SS_ref_db.bounds_ref[0][1] = 1.0-eps;
    SS_ref_db.bounds_ref[1][0] = 0.0+eps;  SS_ref_db.bounds_ref[1][1] = 1.0-eps;


    return SS_ref_db;
}

/**
    Solution phase data for sb21_ri
*/
SS_ref G_SS_sb21_ri_function(SS_ref SS_ref_db, char* research_group, int EM_dataset, int len_ox, bulk_info z_b, double eps){

    int i, j;
    int n_em = SS_ref_db.n_em;

    char   *EM_tmp[]            = {"mgri","feri",};
    for (int i = 0; i < SS_ref_db.n_em; i++){
        strcpy(SS_ref_db.EM_list[i],EM_tmp[i]);
    };

    // Site mixing composition;
    SS_ref_db.C[0][0] = 0.0;    SS_ref_db.C[0][1] = 1.0;
    SS_ref_db.C[1][0] = 1.0;    SS_ref_db.C[1][1] = 0.0;

    // pre-computed Nullspace;
    SS_ref_db.N[0][0] = -0.7071067811865475;
    SS_ref_db.N[1][0] = 0.7071067811865476;

    SS_ref_db.W[0] = 7600.74 ;


    em_data mgri            = get_em_data(    research_group, EM_dataset,
                                            len_ox,
                                            z_b,
                                            SS_ref_db.P,
                                            SS_ref_db.T,
                                            "mgri",
                                            "equilibrium");

    em_data feri            = get_em_data(    research_group, EM_dataset,
                                            len_ox,
                                            z_b,
                                            SS_ref_db.P,
                                            SS_ref_db.T,
                                            "feri",
                                            "equilibrium");

    SS_ref_db.gbase[0]        = mgri.gb;
    SS_ref_db.gbase[1]        = feri.gb;

    SS_ref_db.ElShearMod[0]        = mgri.ElShearMod;
    SS_ref_db.ElShearMod[1]        = feri.ElShearMod;

    SS_ref_db.ElBulkMod[0]        = mgri.ElBulkMod;
    SS_ref_db.ElBulkMod[1]        = feri.ElBulkMod;

    SS_ref_db.ElCp[0]        = mgri.ElCp;
    SS_ref_db.ElCp[1]        = feri.ElCp;

    SS_ref_db.ElExpansivity[0]        = mgri.ElExpansivity;
    SS_ref_db.ElExpansivity[1]        = feri.ElExpansivity;

    for (i = 0; i < len_ox; i++){
        SS_ref_db.Comp[0][i]    = mgri.C[i];
        SS_ref_db.Comp[1][i]    = feri.C[i];
    }

    for (i = 0; i < n_em; i++){
        SS_ref_db.z_em[i] = 1.0;
    };

    SS_ref_db.bounds_ref[0][0] = 0.0+eps;  SS_ref_db.bounds_ref[0][1] = 1.0-eps;
    SS_ref_db.bounds_ref[1][0] = 0.0+eps;  SS_ref_db.bounds_ref[1][1] = 1.0-eps;


    return SS_ref_db;
}

/**
    Solution phase data for sb21_opx
*/
SS_ref G_SS_sb21_opx_function(SS_ref SS_ref_db, char* research_group, int EM_dataset, int len_ox, bulk_info z_b, double eps){

    int i, j;
    int n_em = SS_ref_db.n_em;

    char   *EM_tmp[]            = {"en","fs","mgts","odi",};
    for (int i = 0; i < SS_ref_db.n_em; i++){
        strcpy(SS_ref_db.EM_list[i],EM_tmp[i]);
    };

    // Site mixing composition;
    SS_ref_db.C[0][0] = 0.0;    SS_ref_db.C[0][1] = 0.0;    SS_ref_db.C[0][2] = 0.0;    SS_ref_db.C[0][3] = 1.0;
    SS_ref_db.C[1][0] = 0.0;    SS_ref_db.C[1][1] = 1.0;    SS_ref_db.C[1][2] = 0.0;    SS_ref_db.C[1][3] = 0.0;
    SS_ref_db.C[2][0] = 1.0;    SS_ref_db.C[2][1] = 0.0;    SS_ref_db.C[2][2] = 1.0;    SS_ref_db.C[2][3] = 0.0;
    SS_ref_db.C[3][0] = 0.0;    SS_ref_db.C[3][1] = 0.0;    SS_ref_db.C[3][2] = 1.0;    SS_ref_db.C[3][3] = 0.0;
    SS_ref_db.C[4][0] = 0.0;    SS_ref_db.C[4][1] = 1.0;    SS_ref_db.C[4][2] = 0.0;    SS_ref_db.C[4][3] = 0.0;
    SS_ref_db.C[5][0] = 1.0;    SS_ref_db.C[5][1] = 0.0;    SS_ref_db.C[5][2] = 0.0;    SS_ref_db.C[5][3] = 1.0;

    // pre-computed Nullspace;
    SS_ref_db.N[0][0] = -0.5;    SS_ref_db.N[0][1] = -0.5;    SS_ref_db.N[0][2] = -0.5;
    SS_ref_db.N[1][0] = 0.8333333333333334;    SS_ref_db.N[1][1] = -0.16666666666666666;    SS_ref_db.N[1][2] = -0.16666666666666666;
    SS_ref_db.N[2][0] = -0.16666666666666666;    SS_ref_db.N[2][1] = 0.8333333333333334;    SS_ref_db.N[2][2] = -0.16666666666666666;
    SS_ref_db.N[3][0] = -0.16666666666666666;    SS_ref_db.N[3][1] = -0.16666666666666666;    SS_ref_db.N[3][2] = 0.8333333333333334;

    SS_ref_db.W[0] = 0.0 ;
    SS_ref_db.W[1] = 0.0 ;
    SS_ref_db.W[2] = 32217.44 ;
    SS_ref_db.W[3] = 0.0 ;
    SS_ref_db.W[4] = 32217.44 ;
    SS_ref_db.W[5] = 48370.41 ;


    em_data en            = get_em_data(    research_group, EM_dataset,
                                            len_ox,
                                            z_b,
                                            SS_ref_db.P,
                                            SS_ref_db.T,
                                            "en",
                                            "equilibrium");

    em_data fs            = get_em_data(    research_group, EM_dataset,
                                            len_ox,
                                            z_b,
                                            SS_ref_db.P,
                                            SS_ref_db.T,
                                            "fs",
                                            "equilibrium");

    em_data mgts            = get_em_data(    research_group, EM_dataset,
                                            len_ox,
                                            z_b,
                                            SS_ref_db.P,
                                            SS_ref_db.T,
                                            "mgts",
                                            "equilibrium");

    em_data odi            = get_em_data(    research_group, EM_dataset,
                                            len_ox,
                                            z_b,
                                            SS_ref_db.P,
                                            SS_ref_db.T,
                                            "odi",
                                            "equilibrium");

    SS_ref_db.gbase[0]        = en.gb;
    SS_ref_db.gbase[1]        = fs.gb;
    SS_ref_db.gbase[2]        = mgts.gb;
    SS_ref_db.gbase[3]        = odi.gb;

    SS_ref_db.ElShearMod[0]        = en.ElShearMod;
    SS_ref_db.ElShearMod[1]        = fs.ElShearMod;
    SS_ref_db.ElShearMod[2]        = mgts.ElShearMod;
    SS_ref_db.ElShearMod[3]        = odi.ElShearMod;

    SS_ref_db.ElBulkMod[0]        = en.ElBulkMod;
    SS_ref_db.ElBulkMod[1]        = fs.ElBulkMod;
    SS_ref_db.ElBulkMod[2]        = mgts.ElBulkMod;
    SS_ref_db.ElBulkMod[3]        = odi.ElBulkMod;

    SS_ref_db.ElCp[0]        = en.ElCp;
    SS_ref_db.ElCp[1]        = fs.ElCp;
    SS_ref_db.ElCp[2]        = mgts.ElCp;
    SS_ref_db.ElCp[3]        = odi.ElCp;

    SS_ref_db.ElExpansivity[0]        = en.ElExpansivity;
    SS_ref_db.ElExpansivity[1]        = fs.ElExpansivity;
    SS_ref_db.ElExpansivity[2]        = mgts.ElExpansivity;
    SS_ref_db.ElExpansivity[3]        = odi.ElExpansivity;

    for (i = 0; i < len_ox; i++){
        SS_ref_db.Comp[0][i]    = en.C[i];
        SS_ref_db.Comp[1][i]    = fs.C[i];
        SS_ref_db.Comp[2][i]    = mgts.C[i];
        SS_ref_db.Comp[3][i]    = odi.C[i];
    }

    for (i = 0; i < n_em; i++){
        SS_ref_db.z_em[i] = 1.0;
    };

    SS_ref_db.bounds_ref[0][0] = 0.0+eps;  SS_ref_db.bounds_ref[0][1] = 1.0-eps;
    SS_ref_db.bounds_ref[1][0] = 0.0+eps;  SS_ref_db.bounds_ref[1][1] = 1.0-eps;
    SS_ref_db.bounds_ref[2][0] = 0.0+eps;  SS_ref_db.bounds_ref[2][1] = 1.0-eps;
    SS_ref_db.bounds_ref[3][0] = 0.0+eps;  SS_ref_db.bounds_ref[3][1] = 1.0-eps;


    return SS_ref_db;
}

/**
    Solution phase data for sb21_cpx
*/
SS_ref G_SS_sb21_cpx_function(SS_ref SS_ref_db, char* research_group, int EM_dataset, int len_ox, bulk_info z_b, double eps){

    int i, j;
    int n_em = SS_ref_db.n_em;

    char   *EM_tmp[]            = {"di","he","cen","cats","jd",};
    for (int i = 0; i < SS_ref_db.n_em; i++){
        strcpy(SS_ref_db.EM_list[i],EM_tmp[i]);
    };

    // Site mixing composition;
    SS_ref_db.C[0][0] = 1.0;    SS_ref_db.C[0][1] = 1.0;    SS_ref_db.C[0][2] = 0.0;    SS_ref_db.C[0][3] = 1.0;    SS_ref_db.C[0][4] = 0.0;
    SS_ref_db.C[1][0] = 0.0;    SS_ref_db.C[1][1] = 0.0;    SS_ref_db.C[1][2] = 1.0;    SS_ref_db.C[1][3] = 0.0;    SS_ref_db.C[1][4] = 0.0;
    SS_ref_db.C[2][0] = 0.0;    SS_ref_db.C[2][1] = 0.0;    SS_ref_db.C[2][2] = 0.0;    SS_ref_db.C[2][3] = 0.0;    SS_ref_db.C[2][4] = 1.0;
    SS_ref_db.C[3][0] = 0.0;    SS_ref_db.C[3][1] = 0.0;    SS_ref_db.C[3][2] = 0.0;    SS_ref_db.C[3][3] = 1.0;    SS_ref_db.C[3][4] = 1.0;
    SS_ref_db.C[4][0] = 0.0;    SS_ref_db.C[4][1] = 1.0;    SS_ref_db.C[4][2] = 0.0;    SS_ref_db.C[4][3] = 0.0;    SS_ref_db.C[4][4] = 0.0;
    SS_ref_db.C[5][0] = 1.0;    SS_ref_db.C[5][1] = 0.0;    SS_ref_db.C[5][2] = 1.0;    SS_ref_db.C[5][3] = 0.0;    SS_ref_db.C[5][4] = 0.0;
    SS_ref_db.C[6][0] = 1.0;    SS_ref_db.C[6][1] = 1.0;    SS_ref_db.C[6][2] = 1.0;    SS_ref_db.C[6][3] = 0.5;    SS_ref_db.C[6][4] = 1.0;
    SS_ref_db.C[7][0] = 0.0;    SS_ref_db.C[7][1] = 0.0;    SS_ref_db.C[7][2] = 0.0;    SS_ref_db.C[7][3] = 0.5;    SS_ref_db.C[7][4] = 0.0;

    // pre-computed Nullspace;
    SS_ref_db.N[0][0] = -0.44721359549995787;    SS_ref_db.N[0][1] = -0.44721359549995787;    SS_ref_db.N[0][2] = -0.44721359549995787;    SS_ref_db.N[0][3] = -0.44721359549995787;
    SS_ref_db.N[1][0] = 0.8618033988749895;    SS_ref_db.N[1][1] = -0.13819660112501048;    SS_ref_db.N[1][2] = -0.13819660112501048;    SS_ref_db.N[1][3] = -0.13819660112501048;
    SS_ref_db.N[2][0] = -0.13819660112501048;    SS_ref_db.N[2][1] = 0.8618033988749895;    SS_ref_db.N[2][2] = -0.13819660112501048;    SS_ref_db.N[2][3] = -0.13819660112501048;
    SS_ref_db.N[3][0] = -0.13819660112501048;    SS_ref_db.N[3][1] = -0.13819660112501048;    SS_ref_db.N[3][2] = 0.8618033988749895;    SS_ref_db.N[3][3] = -0.13819660112501048;
    SS_ref_db.N[4][0] = -0.13819660112501048;    SS_ref_db.N[4][1] = -0.13819660112501048;    SS_ref_db.N[4][2] = -0.13819660112501048;    SS_ref_db.N[4][3] = 0.8618033988749895;

    SS_ref_db.W[0] = 0.0 ;
    SS_ref_db.W[1] = 24740.0 ;
    SS_ref_db.W[2] = 26000.0 ;
    SS_ref_db.W[3] = 24300.0 ;
    SS_ref_db.W[4] = 24740.0 ;
    SS_ref_db.W[5] = 26000.0 ;
    SS_ref_db.W[6] = 24300.0 ;
    SS_ref_db.W[7] = 60132.81 ;
    SS_ref_db.W[8] = 46046.07 ;
    SS_ref_db.W[9] = 10000.0 ;

    SS_ref_db.v[0] = 1.0;
    SS_ref_db.v[1] = 1.0;
    SS_ref_db.v[2] = 1.0;
    SS_ref_db.v[3] = 3.5;
    SS_ref_db.v[4] = 1.0;

    em_data di            = get_em_data(    research_group, EM_dataset,
                                            len_ox,
                                            z_b,
                                            SS_ref_db.P,
                                            SS_ref_db.T,
                                            "di",
                                            "equilibrium");

    em_data he            = get_em_data(    research_group, EM_dataset,
                                            len_ox,
                                            z_b,
                                            SS_ref_db.P,
                                            SS_ref_db.T,
                                            "he",
                                            "equilibrium");

    em_data cen            = get_em_data(    research_group, EM_dataset,
                                            len_ox,
                                            z_b,
                                            SS_ref_db.P,
                                            SS_ref_db.T,
                                            "cen",
                                            "equilibrium");

    em_data cats            = get_em_data(    research_group, EM_dataset,
                                            len_ox,
                                            z_b,
                                            SS_ref_db.P,
                                            SS_ref_db.T,
                                            "cats",
                                            "equilibrium");

    em_data jd            = get_em_data(    research_group, EM_dataset,
                                            len_ox,
                                            z_b,
                                            SS_ref_db.P,
                                            SS_ref_db.T,
                                            "jd",
                                            "equilibrium");

    SS_ref_db.gbase[0]        = di.gb;
    SS_ref_db.gbase[1]        = he.gb;
    SS_ref_db.gbase[2]        = cen.gb;
    SS_ref_db.gbase[3]        = cats.gb;
    SS_ref_db.gbase[4]        = jd.gb;

    SS_ref_db.ElShearMod[0]        = di.ElShearMod;
    SS_ref_db.ElShearMod[1]        = he.ElShearMod;
    SS_ref_db.ElShearMod[2]        = cen.ElShearMod;
    SS_ref_db.ElShearMod[3]        = cats.ElShearMod;
    SS_ref_db.ElShearMod[4]        = jd.ElShearMod;

    SS_ref_db.ElBulkMod[0]        = di.ElBulkMod;
    SS_ref_db.ElBulkMod[1]        = he.ElBulkMod;
    SS_ref_db.ElBulkMod[2]        = cen.ElBulkMod;
    SS_ref_db.ElBulkMod[3]        = cats.ElBulkMod;
    SS_ref_db.ElBulkMod[4]        = jd.ElBulkMod;

    SS_ref_db.ElCp[0]        = di.ElCp;
    SS_ref_db.ElCp[1]        = he.ElCp;
    SS_ref_db.ElCp[2]        = cen.ElCp;
    SS_ref_db.ElCp[3]        = cats.ElCp;
    SS_ref_db.ElCp[4]        = jd.ElCp;

    SS_ref_db.ElExpansivity[0]        = di.ElExpansivity;
    SS_ref_db.ElExpansivity[1]        = he.ElExpansivity;
    SS_ref_db.ElExpansivity[2]        = cen.ElExpansivity;
    SS_ref_db.ElExpansivity[3]        = cats.ElExpansivity;
    SS_ref_db.ElExpansivity[4]        = jd.ElExpansivity;

    for (i = 0; i < len_ox; i++){
        SS_ref_db.Comp[0][i]    = di.C[i];
        SS_ref_db.Comp[1][i]    = he.C[i];
        SS_ref_db.Comp[2][i]    = cen.C[i];
        SS_ref_db.Comp[3][i]    = cats.C[i];
        SS_ref_db.Comp[4][i]    = jd.C[i];
    }

    for (i = 0; i < n_em; i++){
        SS_ref_db.z_em[i] = 1.0;
    };

    SS_ref_db.bounds_ref[0][0] = 0.0+eps;  SS_ref_db.bounds_ref[0][1] = 1.0-eps;
    SS_ref_db.bounds_ref[1][0] = 0.0+eps;  SS_ref_db.bounds_ref[1][1] = 1.0-eps;
    SS_ref_db.bounds_ref[2][0] = 0.0+eps;  SS_ref_db.bounds_ref[2][1] = 1.0-eps;
    SS_ref_db.bounds_ref[3][0] = 0.0+eps;  SS_ref_db.bounds_ref[3][1] = 1.0-eps;
    SS_ref_db.bounds_ref[4][0] = 0.0+eps;  SS_ref_db.bounds_ref[4][1] = 1.0-eps;


    return SS_ref_db;
}

/**
    Solution phase data for sb21_hpcpx
*/
SS_ref G_SS_sb21_hpcpx_function(SS_ref SS_ref_db, char* research_group, int EM_dataset, int len_ox, bulk_info z_b, double eps){

    int i, j;
    int n_em = SS_ref_db.n_em;

    char   *EM_tmp[]            = {"hpcen","hpcfs",};
    for (int i = 0; i < SS_ref_db.n_em; i++){
        strcpy(SS_ref_db.EM_list[i],EM_tmp[i]);
    };

    // Site mixing composition;
    SS_ref_db.C[0][0] = 0.0;    SS_ref_db.C[0][1] = 1.0;
    SS_ref_db.C[1][0] = 1.0;    SS_ref_db.C[1][1] = 0.0;

    // pre-computed Nullspace;
    SS_ref_db.N[0][0] = -0.7071067811865475;
    SS_ref_db.N[1][0] = 0.7071067811865476;

    SS_ref_db.W[0] = 0.0 ;


    em_data hpcen            = get_em_data(    research_group, EM_dataset,
                                            len_ox,
                                            z_b,
                                            SS_ref_db.P,
                                            SS_ref_db.T,
                                            "hpcen",
                                            "equilibrium");

    em_data hpcfs            = get_em_data(    research_group, EM_dataset,
                                            len_ox,
                                            z_b,
                                            SS_ref_db.P,
                                            SS_ref_db.T,
                                            "hpcfs",
                                            "equilibrium");

    SS_ref_db.gbase[0]        = hpcen.gb;
    SS_ref_db.gbase[1]        = hpcfs.gb;

    SS_ref_db.ElShearMod[0]        = hpcen.ElShearMod;
    SS_ref_db.ElShearMod[1]        = hpcfs.ElShearMod;

    SS_ref_db.ElBulkMod[0]        = hpcen.ElBulkMod;
    SS_ref_db.ElBulkMod[1]        = hpcfs.ElBulkMod;

    SS_ref_db.ElCp[0]        = hpcen.ElCp;
    SS_ref_db.ElCp[1]        = hpcfs.ElCp;

    SS_ref_db.ElExpansivity[0]        = hpcen.ElExpansivity;
    SS_ref_db.ElExpansivity[1]        = hpcfs.ElExpansivity;

    for (i = 0; i < len_ox; i++){
        SS_ref_db.Comp[0][i]    = hpcen.C[i];
        SS_ref_db.Comp[1][i]    = hpcfs.C[i];
    }

    for (i = 0; i < n_em; i++){
        SS_ref_db.z_em[i] = 1.0;
    };

    SS_ref_db.bounds_ref[0][0] = 0.0+eps;  SS_ref_db.bounds_ref[0][1] = 1.0-eps;
    SS_ref_db.bounds_ref[1][0] = 0.0+eps;  SS_ref_db.bounds_ref[1][1] = 1.0-eps;


    return SS_ref_db;
}

/**
    Solution phase data for sb21_ak
*/
SS_ref G_SS_sb21_ak_function(SS_ref SS_ref_db, char* research_group, int EM_dataset, int len_ox, bulk_info z_b, double eps){

    int i, j;
    int n_em = SS_ref_db.n_em;

    char   *EM_tmp[]            = {"mgak","feak","co",};
    for (int i = 0; i < SS_ref_db.n_em; i++){
        strcpy(SS_ref_db.EM_list[i],EM_tmp[i]);
    };

    // Site mixing composition;
    SS_ref_db.C[0][0] = 0.0;    SS_ref_db.C[0][1] = 0.0;    SS_ref_db.C[0][2] = 1.0;
    SS_ref_db.C[1][0] = 0.0;    SS_ref_db.C[1][1] = 1.0;    SS_ref_db.C[1][2] = 0.0;
    SS_ref_db.C[2][0] = 1.0;    SS_ref_db.C[2][1] = 0.0;    SS_ref_db.C[2][2] = 0.0;
    SS_ref_db.C[3][0] = 1.0;    SS_ref_db.C[3][1] = 1.0;    SS_ref_db.C[3][2] = 0.0;
    SS_ref_db.C[4][0] = 0.0;    SS_ref_db.C[4][1] = 0.0;    SS_ref_db.C[4][2] = 1.0;

    // pre-computed Nullspace;
    SS_ref_db.N[0][0] = -0.5773502691896257;    SS_ref_db.N[0][1] = -0.5773502691896257;
    SS_ref_db.N[1][0] = 0.7886751345948129;    SS_ref_db.N[1][1] = -0.2113248654051871;
    SS_ref_db.N[2][0] = -0.2113248654051871;    SS_ref_db.N[2][1] = 0.7886751345948129;

    SS_ref_db.W[0] = 0.0 ;
    SS_ref_db.W[1] = 59348.69 ;
    SS_ref_db.W[2] = 59348.69 ;


    em_data mgak            = get_em_data(    research_group, EM_dataset,
                                            len_ox,
                                            z_b,
                                            SS_ref_db.P,
                                            SS_ref_db.T,
                                            "mgak",
                                            "equilibrium");

    em_data feak            = get_em_data(    research_group, EM_dataset,
                                            len_ox,
                                            z_b,
                                            SS_ref_db.P,
                                            SS_ref_db.T,
                                            "feak",
                                            "equilibrium");

    em_data co            = get_em_data(    research_group, EM_dataset,
                                            len_ox,
                                            z_b,
                                            SS_ref_db.P,
                                            SS_ref_db.T,
                                            "co",
                                            "equilibrium");

    SS_ref_db.gbase[0]        = mgak.gb;
    SS_ref_db.gbase[1]        = feak.gb;
    SS_ref_db.gbase[2]        = co.gb;

    SS_ref_db.ElShearMod[0]        = mgak.ElShearMod;
    SS_ref_db.ElShearMod[1]        = feak.ElShearMod;
    SS_ref_db.ElShearMod[2]        = co.ElShearMod;

    SS_ref_db.ElBulkMod[0]        = mgak.ElBulkMod;
    SS_ref_db.ElBulkMod[1]        = feak.ElBulkMod;
    SS_ref_db.ElBulkMod[2]        = co.ElBulkMod;

    SS_ref_db.ElCp[0]        = mgak.ElCp;
    SS_ref_db.ElCp[1]        = feak.ElCp;
    SS_ref_db.ElCp[2]        = co.ElCp;

    SS_ref_db.ElExpansivity[0]        = mgak.ElExpansivity;
    SS_ref_db.ElExpansivity[1]        = feak.ElExpansivity;
    SS_ref_db.ElExpansivity[2]        = co.ElExpansivity;

    for (i = 0; i < len_ox; i++){
        SS_ref_db.Comp[0][i]    = mgak.C[i];
        SS_ref_db.Comp[1][i]    = feak.C[i];
        SS_ref_db.Comp[2][i]    = co.C[i];
    }

    for (i = 0; i < n_em; i++){
        SS_ref_db.z_em[i] = 1.0;
    };

    SS_ref_db.bounds_ref[0][0] = 0.0+eps;  SS_ref_db.bounds_ref[0][1] = 1.0-eps;
    SS_ref_db.bounds_ref[1][0] = 0.0+eps;  SS_ref_db.bounds_ref[1][1] = 1.0-eps;
    SS_ref_db.bounds_ref[2][0] = 0.0+eps;  SS_ref_db.bounds_ref[2][1] = 1.0-eps;


    return SS_ref_db;
}

/**
    Solution phase data for sb21_gtmj
*/
SS_ref G_SS_sb21_gtmj_function(SS_ref SS_ref_db, char* research_group, int EM_dataset, int len_ox, bulk_info z_b, double eps){

    int i, j;
    int n_em = SS_ref_db.n_em;

    char   *EM_tmp[]            = {"py","alm","gr","mgmj","jdmj",};
    for (int i = 0; i < SS_ref_db.n_em; i++){
        strcpy(SS_ref_db.EM_list[i],EM_tmp[i]);
    };

    // Site mixing composition;
    SS_ref_db.C[0][0] = 0.0;    SS_ref_db.C[0][1] = 0.0;    SS_ref_db.C[0][2] = 1.0;    SS_ref_db.C[0][3] = 0.0;    SS_ref_db.C[0][4] = 0.0;
    SS_ref_db.C[1][0] = 0.0;    SS_ref_db.C[1][1] = 0.0;    SS_ref_db.C[1][2] = 0.0;    SS_ref_db.C[1][3] = 0.0;    SS_ref_db.C[1][4] = 0.3333333333333333;
    SS_ref_db.C[2][0] = 0.0;    SS_ref_db.C[2][1] = 1.0;    SS_ref_db.C[2][2] = 0.0;    SS_ref_db.C[2][3] = 0.0;    SS_ref_db.C[2][4] = 0.0;
    SS_ref_db.C[3][0] = 1.0;    SS_ref_db.C[3][1] = 0.0;    SS_ref_db.C[3][2] = 0.0;    SS_ref_db.C[3][3] = 1.0;    SS_ref_db.C[3][4] = 0.0;
    SS_ref_db.C[4][0] = 0.0;    SS_ref_db.C[4][1] = 0.0;    SS_ref_db.C[4][2] = 0.0;    SS_ref_db.C[4][3] = 0.0;    SS_ref_db.C[4][4] = 0.6666666666666666;
    SS_ref_db.C[5][0] = 1.0;    SS_ref_db.C[5][1] = 1.0;    SS_ref_db.C[5][2] = 1.0;    SS_ref_db.C[5][3] = 0.0;    SS_ref_db.C[5][4] = 1.0;
    SS_ref_db.C[6][0] = 0.0;    SS_ref_db.C[6][1] = 0.0;    SS_ref_db.C[6][2] = 0.0;    SS_ref_db.C[6][3] = 1.0;    SS_ref_db.C[6][4] = 0.0;
    SS_ref_db.C[7][0] = 0.0;    SS_ref_db.C[7][1] = 0.0;    SS_ref_db.C[7][2] = 0.0;    SS_ref_db.C[7][3] = 1.0;    SS_ref_db.C[7][4] = 1.0;
    SS_ref_db.C[8][0] = 1.0;    SS_ref_db.C[8][1] = 1.0;    SS_ref_db.C[8][2] = 1.0;    SS_ref_db.C[8][3] = 0.0;    SS_ref_db.C[8][4] = 0.0;

    // pre-computed Nullspace;
    SS_ref_db.N[0][0] = -0.44721359549995787;    SS_ref_db.N[0][1] = -0.44721359549995787;    SS_ref_db.N[0][2] = -0.44721359549995787;    SS_ref_db.N[0][3] = -0.44721359549995787;
    SS_ref_db.N[1][0] = 0.8618033988749895;    SS_ref_db.N[1][1] = -0.13819660112501048;    SS_ref_db.N[1][2] = -0.13819660112501048;    SS_ref_db.N[1][3] = -0.13819660112501048;
    SS_ref_db.N[2][0] = -0.13819660112501048;    SS_ref_db.N[2][1] = 0.8618033988749895;    SS_ref_db.N[2][2] = -0.13819660112501048;    SS_ref_db.N[2][3] = -0.13819660112501048;
    SS_ref_db.N[3][0] = -0.13819660112501048;    SS_ref_db.N[3][1] = -0.13819660112501048;    SS_ref_db.N[3][2] = 0.8618033988749895;    SS_ref_db.N[3][3] = -0.13819660112501048;
    SS_ref_db.N[4][0] = -0.13819660112501048;    SS_ref_db.N[4][1] = -0.13819660112501048;    SS_ref_db.N[4][2] = -0.13819660112501048;    SS_ref_db.N[4][3] = 0.8618033988749895;

    SS_ref_db.W[0] = 0.0 ;
    SS_ref_db.W[1] = 21117.58 + 1.03e-5 * SS_ref_db.P;
    SS_ref_db.W[2] = 22672.42 ;
    SS_ref_db.W[3] = 22672.42 ;
    SS_ref_db.W[4] = 21117.58 ;
    SS_ref_db.W[5] = 22672.42 ;
    SS_ref_db.W[6] = 22672.42 ;
    SS_ref_db.W[7] = 60718.2 + 1.03e-5 * SS_ref_db.P;
    SS_ref_db.W[8] = 60718.2 ;
    SS_ref_db.W[9] = 70879.14 ;


    em_data py            = get_em_data(    research_group, EM_dataset,
                                            len_ox,
                                            z_b,
                                            SS_ref_db.P,
                                            SS_ref_db.T,
                                            "py",
                                            "equilibrium");

    em_data alm            = get_em_data(    research_group, EM_dataset,
                                            len_ox,
                                            z_b,
                                            SS_ref_db.P,
                                            SS_ref_db.T,
                                            "alm",
                                            "equilibrium");

    em_data gr            = get_em_data(    research_group, EM_dataset,
                                            len_ox,
                                            z_b,
                                            SS_ref_db.P,
                                            SS_ref_db.T,
                                            "gr",
                                            "equilibrium");

    em_data mgmj            = get_em_data(    research_group, EM_dataset,
                                            len_ox,
                                            z_b,
                                            SS_ref_db.P,
                                            SS_ref_db.T,
                                            "mgmj",
                                            "equilibrium");

    em_data jdmj            = get_em_data(    research_group, EM_dataset,
                                            len_ox,
                                            z_b,
                                            SS_ref_db.P,
                                            SS_ref_db.T,
                                            "jdmj",
                                            "equilibrium");

    SS_ref_db.gbase[0]        = py.gb;
    SS_ref_db.gbase[1]        = alm.gb;
    SS_ref_db.gbase[2]        = gr.gb;
    SS_ref_db.gbase[3]        = mgmj.gb;
    SS_ref_db.gbase[4]        = jdmj.gb;

    SS_ref_db.ElShearMod[0]        = py.ElShearMod;
    SS_ref_db.ElShearMod[1]        = alm.ElShearMod;
    SS_ref_db.ElShearMod[2]        = gr.ElShearMod;
    SS_ref_db.ElShearMod[3]        = mgmj.ElShearMod;
    SS_ref_db.ElShearMod[4]        = jdmj.ElShearMod;

    SS_ref_db.ElBulkMod[0]        = py.ElBulkMod;
    SS_ref_db.ElBulkMod[1]        = alm.ElBulkMod;
    SS_ref_db.ElBulkMod[2]        = gr.ElBulkMod;
    SS_ref_db.ElBulkMod[3]        = mgmj.ElBulkMod;
    SS_ref_db.ElBulkMod[4]        = jdmj.ElBulkMod;

    SS_ref_db.ElCp[0]        = py.ElCp;
    SS_ref_db.ElCp[1]        = alm.ElCp;
    SS_ref_db.ElCp[2]        = gr.ElCp;
    SS_ref_db.ElCp[3]        = mgmj.ElCp;
    SS_ref_db.ElCp[4]        = jdmj.ElCp;

    SS_ref_db.ElExpansivity[0]        = py.ElExpansivity;
    SS_ref_db.ElExpansivity[1]        = alm.ElExpansivity;
    SS_ref_db.ElExpansivity[2]        = gr.ElExpansivity;
    SS_ref_db.ElExpansivity[3]        = mgmj.ElExpansivity;
    SS_ref_db.ElExpansivity[4]        = jdmj.ElExpansivity;

    for (i = 0; i < len_ox; i++){
        SS_ref_db.Comp[0][i]    = py.C[i];
        SS_ref_db.Comp[1][i]    = alm.C[i];
        SS_ref_db.Comp[2][i]    = gr.C[i];
        SS_ref_db.Comp[3][i]    = mgmj.C[i];
        SS_ref_db.Comp[4][i]    = jdmj.C[i];
    }

    for (i = 0; i < n_em; i++){
        SS_ref_db.z_em[i] = 1.0;
    };

    SS_ref_db.bounds_ref[0][0] = 0.0+eps;  SS_ref_db.bounds_ref[0][1] = 1.0-eps;
    SS_ref_db.bounds_ref[1][0] = 0.0+eps;  SS_ref_db.bounds_ref[1][1] = 1.0-eps;
    SS_ref_db.bounds_ref[2][0] = 0.0+eps;  SS_ref_db.bounds_ref[2][1] = 1.0-eps;
    SS_ref_db.bounds_ref[3][0] = 0.0+eps;  SS_ref_db.bounds_ref[3][1] = 1.0-eps;
    SS_ref_db.bounds_ref[4][0] = 0.0+eps;  SS_ref_db.bounds_ref[4][1] = 1.0-eps;


    return SS_ref_db;
}

/**
    Solution phase data for sb21_pv
*/
SS_ref G_SS_sb21_pv_function(SS_ref SS_ref_db, char* research_group, int EM_dataset, int len_ox, bulk_info z_b, double eps){

    int i, j;
    int n_em = SS_ref_db.n_em;

    char   *EM_tmp[]            = {"mgpv","fepv","alpv",};
    for (int i = 0; i < SS_ref_db.n_em; i++){
        strcpy(SS_ref_db.EM_list[i],EM_tmp[i]);
    };

    // Site mixing composition;
    SS_ref_db.C[0][0] = 0.0;    SS_ref_db.C[0][1] = 0.0;    SS_ref_db.C[0][2] = 1.0;
    SS_ref_db.C[1][0] = 0.0;    SS_ref_db.C[1][1] = 1.0;    SS_ref_db.C[1][2] = 0.0;
    SS_ref_db.C[2][0] = 1.0;    SS_ref_db.C[2][1] = 0.0;    SS_ref_db.C[2][2] = 0.0;
    SS_ref_db.C[3][0] = 1.0;    SS_ref_db.C[3][1] = 1.0;    SS_ref_db.C[3][2] = 0.0;
    SS_ref_db.C[4][0] = 0.0;    SS_ref_db.C[4][1] = 0.0;    SS_ref_db.C[4][2] = 1.0;

    // pre-computed Nullspace;
    SS_ref_db.N[0][0] = -0.5773502691896257;    SS_ref_db.N[0][1] = -0.5773502691896257;
    SS_ref_db.N[1][0] = 0.7886751345948129;    SS_ref_db.N[1][1] = -0.2113248654051871;
    SS_ref_db.N[2][0] = -0.2113248654051871;    SS_ref_db.N[2][1] = 0.7886751345948129;

    SS_ref_db.W[0] = -11396.17 ;
    SS_ref_db.W[1] = 34979.87 ;
    SS_ref_db.W[2] = 0.0 ;


    em_data mgpv            = get_em_data(    research_group, EM_dataset,
                                            len_ox,
                                            z_b,
                                            SS_ref_db.P,
                                            SS_ref_db.T,
                                            "mgpv",
                                            "equilibrium");

    em_data fepv            = get_em_data(    research_group, EM_dataset,
                                            len_ox,
                                            z_b,
                                            SS_ref_db.P,
                                            SS_ref_db.T,
                                            "fepv",
                                            "equilibrium");

    em_data alpv            = get_em_data(    research_group, EM_dataset,
                                            len_ox,
                                            z_b,
                                            SS_ref_db.P,
                                            SS_ref_db.T,
                                            "alpv",
                                            "equilibrium");

    SS_ref_db.gbase[0]        = mgpv.gb;
    SS_ref_db.gbase[1]        = fepv.gb;
    SS_ref_db.gbase[2]        = alpv.gb;

    SS_ref_db.ElShearMod[0]        = mgpv.ElShearMod;
    SS_ref_db.ElShearMod[1]        = fepv.ElShearMod;
    SS_ref_db.ElShearMod[2]        = alpv.ElShearMod;

    SS_ref_db.ElBulkMod[0]        = mgpv.ElBulkMod;
    SS_ref_db.ElBulkMod[1]        = fepv.ElBulkMod;
    SS_ref_db.ElBulkMod[2]        = alpv.ElBulkMod;

    SS_ref_db.ElCp[0]        = mgpv.ElCp;
    SS_ref_db.ElCp[1]        = fepv.ElCp;
    SS_ref_db.ElCp[2]        = alpv.ElCp;

    SS_ref_db.ElExpansivity[0]        = mgpv.ElExpansivity;
    SS_ref_db.ElExpansivity[1]        = fepv.ElExpansivity;
    SS_ref_db.ElExpansivity[2]        = alpv.ElExpansivity;

    for (i = 0; i < len_ox; i++){
        SS_ref_db.Comp[0][i]    = mgpv.C[i];
        SS_ref_db.Comp[1][i]    = fepv.C[i];
        SS_ref_db.Comp[2][i]    = alpv.C[i];
    }

    for (i = 0; i < n_em; i++){
        SS_ref_db.z_em[i] = 1.0;
    };

    SS_ref_db.bounds_ref[0][0] = 0.0+eps;  SS_ref_db.bounds_ref[0][1] = 1.0-eps;
    SS_ref_db.bounds_ref[1][0] = 0.0+eps;  SS_ref_db.bounds_ref[1][1] = 1.0-eps;
    SS_ref_db.bounds_ref[2][0] = 0.0+eps;  SS_ref_db.bounds_ref[2][1] = 1.0-eps;


    return SS_ref_db;
}

/**
    Solution phase data for sb21_ppv
*/
SS_ref G_SS_sb21_ppv_function(SS_ref SS_ref_db, char* research_group, int EM_dataset, int len_ox, bulk_info z_b, double eps){

    int i, j;
    int n_em = SS_ref_db.n_em;

    char   *EM_tmp[]            = {"mppv","fppv","appv",};
    for (int i = 0; i < SS_ref_db.n_em; i++){
        strcpy(SS_ref_db.EM_list[i],EM_tmp[i]);
    };

    // Site mixing composition;
    SS_ref_db.C[0][0] = 0.0;    SS_ref_db.C[0][1] = 0.0;    SS_ref_db.C[0][2] = 1.0;
    SS_ref_db.C[1][0] = 0.0;    SS_ref_db.C[1][1] = 1.0;    SS_ref_db.C[1][2] = 0.0;
    SS_ref_db.C[2][0] = 1.0;    SS_ref_db.C[2][1] = 0.0;    SS_ref_db.C[2][2] = 0.0;
    SS_ref_db.C[3][0] = 1.0;    SS_ref_db.C[3][1] = 1.0;    SS_ref_db.C[3][2] = 0.0;
    SS_ref_db.C[4][0] = 0.0;    SS_ref_db.C[4][1] = 0.0;    SS_ref_db.C[4][2] = 1.0;

    // pre-computed Nullspace;
    SS_ref_db.N[0][0] = -0.5773502691896257;    SS_ref_db.N[0][1] = -0.5773502691896257;
    SS_ref_db.N[1][0] = 0.7886751345948129;    SS_ref_db.N[1][1] = -0.2113248654051871;
    SS_ref_db.N[2][0] = -0.2113248654051871;    SS_ref_db.N[2][1] = 0.7886751345948129;

    SS_ref_db.W[0] = -10955.49 ;
    SS_ref_db.W[1] = 34979.87 ;
    SS_ref_db.W[2] = 34979.87 ;


    em_data mppv            = get_em_data(    research_group, EM_dataset,
                                            len_ox,
                                            z_b,
                                            SS_ref_db.P,
                                            SS_ref_db.T,
                                            "mppv",
                                            "equilibrium");

    em_data fppv            = get_em_data(    research_group, EM_dataset,
                                            len_ox,
                                            z_b,
                                            SS_ref_db.P,
                                            SS_ref_db.T,
                                            "fppv",
                                            "equilibrium");

    em_data appv            = get_em_data(    research_group, EM_dataset,
                                            len_ox,
                                            z_b,
                                            SS_ref_db.P,
                                            SS_ref_db.T,
                                            "appv",
                                            "equilibrium");

    SS_ref_db.gbase[0]        = mppv.gb;
    SS_ref_db.gbase[1]        = fppv.gb;
    SS_ref_db.gbase[2]        = appv.gb;

    SS_ref_db.ElShearMod[0]        = mppv.ElShearMod;
    SS_ref_db.ElShearMod[1]        = fppv.ElShearMod;
    SS_ref_db.ElShearMod[2]        = appv.ElShearMod;

    SS_ref_db.ElBulkMod[0]        = mppv.ElBulkMod;
    SS_ref_db.ElBulkMod[1]        = fppv.ElBulkMod;
    SS_ref_db.ElBulkMod[2]        = appv.ElBulkMod;

    SS_ref_db.ElCp[0]        = mppv.ElCp;
    SS_ref_db.ElCp[1]        = fppv.ElCp;
    SS_ref_db.ElCp[2]        = appv.ElCp;

    SS_ref_db.ElExpansivity[0]        = mppv.ElExpansivity;
    SS_ref_db.ElExpansivity[1]        = fppv.ElExpansivity;
    SS_ref_db.ElExpansivity[2]        = appv.ElExpansivity;

    for (i = 0; i < len_ox; i++){
        SS_ref_db.Comp[0][i]    = mppv.C[i];
        SS_ref_db.Comp[1][i]    = fppv.C[i];
        SS_ref_db.Comp[2][i]    = appv.C[i];
    }

    for (i = 0; i < n_em; i++){
        SS_ref_db.z_em[i] = 1.0;
    };

    SS_ref_db.bounds_ref[0][0] = 0.0+eps;  SS_ref_db.bounds_ref[0][1] = 1.0-eps;
    SS_ref_db.bounds_ref[1][0] = 0.0+eps;  SS_ref_db.bounds_ref[1][1] = 1.0-eps;
    SS_ref_db.bounds_ref[2][0] = 0.0+eps;  SS_ref_db.bounds_ref[2][1] = 1.0-eps;


    return SS_ref_db;
}

/**
    Solution phase data for sb21_cf
*/
SS_ref G_SS_sb21_cf_function(SS_ref SS_ref_db, char* research_group, int EM_dataset, int len_ox, bulk_info z_b, double eps){

    int i, j;
    int n_em = SS_ref_db.n_em;

    char   *EM_tmp[]            = {"mgcf","fecf","nacf",};
    for (int i = 0; i < SS_ref_db.n_em; i++){
        strcpy(SS_ref_db.EM_list[i],EM_tmp[i]);
    };

    // Site mixing composition;
    SS_ref_db.C[0][0] = 0.0;    SS_ref_db.C[0][1] = 1.0;    SS_ref_db.C[0][2] = 0.0;
    SS_ref_db.C[1][0] = 1.0;    SS_ref_db.C[1][1] = 0.0;    SS_ref_db.C[1][2] = 0.0;
    SS_ref_db.C[2][0] = 0.0;    SS_ref_db.C[2][1] = 0.0;    SS_ref_db.C[2][2] = 1.0;
    SS_ref_db.C[3][0] = 0.0;    SS_ref_db.C[3][1] = 0.0;    SS_ref_db.C[3][2] = 1.0;
    SS_ref_db.C[4][0] = 1.0;    SS_ref_db.C[4][1] = 1.0;    SS_ref_db.C[4][2] = 0.0;

    // pre-computed Nullspace;
    SS_ref_db.N[0][0] = -0.5773502691896257;    SS_ref_db.N[0][1] = -0.5773502691896257;
    SS_ref_db.N[1][0] = 0.7886751345948129;    SS_ref_db.N[1][1] = -0.2113248654051871;
    SS_ref_db.N[2][0] = -0.2113248654051871;    SS_ref_db.N[2][1] = 0.7886751345948129;

    SS_ref_db.W[0] = 0.0 ;
    SS_ref_db.W[1] = 60825.08 ;
    SS_ref_db.W[2] = 60825.08 ;

    SS_ref_db.v[0] = 1.0;
    SS_ref_db.v[1] = 1.0;
    SS_ref_db.v[2] = 4.4532;

    em_data mgcf            = get_em_data(    research_group, EM_dataset,
                                            len_ox,
                                            z_b,
                                            SS_ref_db.P,
                                            SS_ref_db.T,
                                            "mgcf",
                                            "equilibrium");

    em_data fecf            = get_em_data(    research_group, EM_dataset,
                                            len_ox,
                                            z_b,
                                            SS_ref_db.P,
                                            SS_ref_db.T,
                                            "fecf",
                                            "equilibrium");

    em_data nacf            = get_em_data(    research_group, EM_dataset,
                                            len_ox,
                                            z_b,
                                            SS_ref_db.P,
                                            SS_ref_db.T,
                                            "nacf",
                                            "equilibrium");

    SS_ref_db.gbase[0]        = mgcf.gb;
    SS_ref_db.gbase[1]        = fecf.gb;
    SS_ref_db.gbase[2]        = nacf.gb;

    SS_ref_db.ElShearMod[0]        = mgcf.ElShearMod;
    SS_ref_db.ElShearMod[1]        = fecf.ElShearMod;
    SS_ref_db.ElShearMod[2]        = nacf.ElShearMod;

    SS_ref_db.ElBulkMod[0]        = mgcf.ElBulkMod;
    SS_ref_db.ElBulkMod[1]        = fecf.ElBulkMod;
    SS_ref_db.ElBulkMod[2]        = nacf.ElBulkMod;

    SS_ref_db.ElCp[0]        = mgcf.ElCp;
    SS_ref_db.ElCp[1]        = fecf.ElCp;
    SS_ref_db.ElCp[2]        = nacf.ElCp;

    SS_ref_db.ElExpansivity[0]        = mgcf.ElExpansivity;
    SS_ref_db.ElExpansivity[1]        = fecf.ElExpansivity;
    SS_ref_db.ElExpansivity[2]        = nacf.ElExpansivity;

    for (i = 0; i < len_ox; i++){
        SS_ref_db.Comp[0][i]    = mgcf.C[i];
        SS_ref_db.Comp[1][i]    = fecf.C[i];
        SS_ref_db.Comp[2][i]    = nacf.C[i];
    }

    for (i = 0; i < n_em; i++){
        SS_ref_db.z_em[i] = 1.0;
    };

    SS_ref_db.bounds_ref[0][0] = 0.0+eps;  SS_ref_db.bounds_ref[0][1] = 1.0-eps;
    SS_ref_db.bounds_ref[1][0] = 0.0+eps;  SS_ref_db.bounds_ref[1][1] = 1.0-eps;
    SS_ref_db.bounds_ref[2][0] = 0.0+eps;  SS_ref_db.bounds_ref[2][1] = 1.0-eps;


    return SS_ref_db;
}

/**
    Solution phase data for sb21_mw
*/
SS_ref G_SS_sb21_mw_function(SS_ref SS_ref_db, char* research_group, int EM_dataset, int len_ox, bulk_info z_b, double eps){

    int i, j;
    int n_em = SS_ref_db.n_em;

    char   *EM_tmp[]            = {"pe","wu","anao",};
    for (int i = 0; i < SS_ref_db.n_em; i++){
        strcpy(SS_ref_db.EM_list[i],EM_tmp[i]);
    };

    // Site mixing composition;
    SS_ref_db.C[0][0] = 0.0;    SS_ref_db.C[0][1] = 1.0;    SS_ref_db.C[0][2] = 0.0;
    SS_ref_db.C[1][0] = 1.0;    SS_ref_db.C[1][1] = 0.0;    SS_ref_db.C[1][2] = 0.0;
    SS_ref_db.C[2][0] = 0.0;    SS_ref_db.C[2][1] = 0.0;    SS_ref_db.C[2][2] = 1.0;
    SS_ref_db.C[3][0] = 0.0;    SS_ref_db.C[3][1] = 0.0;    SS_ref_db.C[3][2] = 1.0;
    SS_ref_db.C[4][0] = 0.0;    SS_ref_db.C[4][1] = 1.0;    SS_ref_db.C[4][2] = 0.0;
    SS_ref_db.C[5][0] = 1.0;    SS_ref_db.C[5][1] = 0.0;    SS_ref_db.C[5][2] = 0.0;

    // pre-computed Nullspace;
    SS_ref_db.N[0][0] = -0.5773502691896257;    SS_ref_db.N[0][1] = -0.5773502691896257;
    SS_ref_db.N[1][0] = 0.7886751345948129;    SS_ref_db.N[1][1] = -0.2113248654051871;
    SS_ref_db.N[2][0] = -0.2113248654051871;    SS_ref_db.N[2][1] = 0.7886751345948129;

    SS_ref_db.W[0] = 44000.0 + 4.4e-6 * SS_ref_db.P;
    SS_ref_db.W[1] = 120000.0 ;
    SS_ref_db.W[2] = 120000.0 ;


    em_data pe            = get_em_data(    research_group, EM_dataset,
                                            len_ox,
                                            z_b,
                                            SS_ref_db.P,
                                            SS_ref_db.T,
                                            "pe",
                                            "equilibrium");

    em_data wu            = get_em_data(    research_group, EM_dataset,
                                            len_ox,
                                            z_b,
                                            SS_ref_db.P,
                                            SS_ref_db.T,
                                            "wu",
                                            "equilibrium");

    em_data anao            = get_em_data(    research_group, EM_dataset,
                                            len_ox,
                                            z_b,
                                            SS_ref_db.P,
                                            SS_ref_db.T,
                                            "anao",
                                            "equilibrium");

    SS_ref_db.gbase[0]        = pe.gb;
    SS_ref_db.gbase[1]        = wu.gb;
    SS_ref_db.gbase[2]        = anao.gb;

    SS_ref_db.ElShearMod[0]        = pe.ElShearMod;
    SS_ref_db.ElShearMod[1]        = wu.ElShearMod;
    SS_ref_db.ElShearMod[2]        = anao.ElShearMod;

    SS_ref_db.ElBulkMod[0]        = pe.ElBulkMod;
    SS_ref_db.ElBulkMod[1]        = wu.ElBulkMod;
    SS_ref_db.ElBulkMod[2]        = anao.ElBulkMod;

    SS_ref_db.ElCp[0]        = pe.ElCp;
    SS_ref_db.ElCp[1]        = wu.ElCp;
    SS_ref_db.ElCp[2]        = anao.ElCp;

    SS_ref_db.ElExpansivity[0]        = pe.ElExpansivity;
    SS_ref_db.ElExpansivity[1]        = wu.ElExpansivity;
    SS_ref_db.ElExpansivity[2]        = anao.ElExpansivity;

    for (i = 0; i < len_ox; i++){
        SS_ref_db.Comp[0][i]    = pe.C[i];
        SS_ref_db.Comp[1][i]    = wu.C[i];
        SS_ref_db.Comp[2][i]    = anao.C[i];
    }

    for (i = 0; i < n_em; i++){
        SS_ref_db.z_em[i] = 1.0;
    };

    SS_ref_db.bounds_ref[0][0] = 0.0+eps;  SS_ref_db.bounds_ref[0][1] = 1.0-eps;
    SS_ref_db.bounds_ref[1][0] = 0.0+eps;  SS_ref_db.bounds_ref[1][1] = 1.0-eps;
    SS_ref_db.bounds_ref[2][0] = 0.0+eps;  SS_ref_db.bounds_ref[2][1] = 1.0-eps;


    return SS_ref_db;
}

/**
    Solution phase data for sb21_nal
*/
SS_ref G_SS_sb21_nal_function(SS_ref SS_ref_db, char* research_group, int EM_dataset, int len_ox, bulk_info z_b, double eps){

    int i, j;
    int n_em = SS_ref_db.n_em;

    char   *EM_tmp[]            = {"mnal","fnal","nnal",};
    for (int i = 0; i < SS_ref_db.n_em; i++){
        strcpy(SS_ref_db.EM_list[i],EM_tmp[i]);
    };

    // Site mixing composition;
    SS_ref_db.C[0][0] = 1.0;    SS_ref_db.C[0][1] = 1.0;    SS_ref_db.C[0][2] = 1.0;
    SS_ref_db.C[1][0] = 0.0;    SS_ref_db.C[1][1] = 1.0;    SS_ref_db.C[1][2] = 0.0;
    SS_ref_db.C[2][0] = 1.0;    SS_ref_db.C[2][1] = 0.0;    SS_ref_db.C[2][2] = 0.0;
    SS_ref_db.C[3][0] = 0.0;    SS_ref_db.C[3][1] = 0.0;    SS_ref_db.C[3][2] = 1.0;
    SS_ref_db.C[4][0] = 0.16666666666666666;    SS_ref_db.C[4][1] = 0.16666666666666666;    SS_ref_db.C[4][2] = 0.5;
    SS_ref_db.C[5][0] = 0.8333333333333334;    SS_ref_db.C[5][1] = 0.8333333333333334;    SS_ref_db.C[5][2] = 0.5;

    // pre-computed Nullspace;
    SS_ref_db.N[0][0] = -0.5773502691896257;    SS_ref_db.N[0][1] = -0.5773502691896257;
    SS_ref_db.N[1][0] = 0.7886751345948129;    SS_ref_db.N[1][1] = -0.2113248654051871;
    SS_ref_db.N[2][0] = -0.2113248654051871;    SS_ref_db.N[2][1] = 0.7886751345948129;

    SS_ref_db.W[0] = 0.0 ;
    SS_ref_db.W[1] = -60781.47 ;
    SS_ref_db.W[2] = -60781.47 ;


    em_data mnal            = get_em_data(    research_group, EM_dataset,
                                            len_ox,
                                            z_b,
                                            SS_ref_db.P,
                                            SS_ref_db.T,
                                            "mnal",
                                            "equilibrium");

    em_data fnal            = get_em_data(    research_group, EM_dataset,
                                            len_ox,
                                            z_b,
                                            SS_ref_db.P,
                                            SS_ref_db.T,
                                            "fnal",
                                            "equilibrium");

    em_data nnal            = get_em_data(    research_group, EM_dataset,
                                            len_ox,
                                            z_b,
                                            SS_ref_db.P,
                                            SS_ref_db.T,
                                            "nnal",
                                            "equilibrium");

    SS_ref_db.gbase[0]        = mnal.gb;
    SS_ref_db.gbase[1]        = fnal.gb;
    SS_ref_db.gbase[2]        = nnal.gb;

    SS_ref_db.ElShearMod[0]        = mnal.ElShearMod;
    SS_ref_db.ElShearMod[1]        = fnal.ElShearMod;
    SS_ref_db.ElShearMod[2]        = nnal.ElShearMod;

    SS_ref_db.ElBulkMod[0]        = mnal.ElBulkMod;
    SS_ref_db.ElBulkMod[1]        = fnal.ElBulkMod;
    SS_ref_db.ElBulkMod[2]        = nnal.ElBulkMod;

    SS_ref_db.ElCp[0]        = mnal.ElCp;
    SS_ref_db.ElCp[1]        = fnal.ElCp;
    SS_ref_db.ElCp[2]        = nnal.ElCp;

    SS_ref_db.ElExpansivity[0]        = mnal.ElExpansivity;
    SS_ref_db.ElExpansivity[1]        = fnal.ElExpansivity;
    SS_ref_db.ElExpansivity[2]        = nnal.ElExpansivity;

    for (i = 0; i < len_ox; i++){
        SS_ref_db.Comp[0][i]    = mnal.C[i];
        SS_ref_db.Comp[1][i]    = fnal.C[i];
        SS_ref_db.Comp[2][i]    = nnal.C[i];
    }

    for (i = 0; i < n_em; i++){
        SS_ref_db.z_em[i] = 1.0;
    };

    SS_ref_db.bounds_ref[0][0] = 0.0+eps;  SS_ref_db.bounds_ref[0][1] = 1.0-eps;
    SS_ref_db.bounds_ref[1][0] = 0.0+eps;  SS_ref_db.bounds_ref[1][1] = 1.0-eps;
    SS_ref_db.bounds_ref[2][0] = 0.0+eps;  SS_ref_db.bounds_ref[2][1] = 1.0-eps;


    return SS_ref_db;
}

SS_ref G_SS_sb21_EM_function(       global_variable          gv,
                                    SS_ref                   SS_ref_db,
                                    int                      EM_dataset,
                                    bulk_info                z_b,
                                    char                    *name){            

    double eps                  = gv.bnd_val;
    double P                    = SS_ref_db.P;
    double T                    = SS_ref_db.T;

    SS_ref_db.ss_flags[0]  = 1;

    /* Associate the right solid-solution data */
    for (int FD = 0; FD < gv.n_Diff; FD++){

        if (FD == 8 || FD == 9){
            SS_ref_db.P = 1.+ gv.gb_P_eps*gv.pdev[0][FD];
            SS_ref_db.T = T + gv.gb_T_eps*gv.pdev[1][FD];
        }
        else{
            SS_ref_db.P = P + gv.gb_P_eps*gv.pdev[0][FD];
            SS_ref_db.T = T + gv.gb_T_eps*gv.pdev[1][FD];
        }

        if (strcmp( name, "plg") == 0 ){
  			if ( P > 500.0){
				SS_ref_db.ss_flags[0]  = 0;
			}
            SS_ref_db  = G_SS_sb21_plg_function(SS_ref_db, gv.research_group, EM_dataset, gv.len_ox, z_b, eps); }
        else if (strcmp( name, "sp") == 0 ){
            SS_ref_db  = G_SS_sb21_sp_function(SS_ref_db, gv.research_group, EM_dataset, gv.len_ox, z_b, eps);  }
        else if (strcmp( name, "ol") == 0 ){
            SS_ref_db  = G_SS_sb21_ol_function(SS_ref_db, gv.research_group, EM_dataset, gv.len_ox, z_b, eps);  }
        else if (strcmp( name, "wa") == 0 ){
            SS_ref_db  = G_SS_sb21_wa_function(SS_ref_db, gv.research_group, EM_dataset, gv.len_ox, z_b, eps);  }
        else if (strcmp( name, "ri") == 0 ){
            SS_ref_db  = G_SS_sb21_ri_function(SS_ref_db, gv.research_group, EM_dataset, gv.len_ox, z_b, eps);  }
        else if (strcmp( name, "opx") == 0 ){
            SS_ref_db  = G_SS_sb21_opx_function(SS_ref_db, gv.research_group, EM_dataset, gv.len_ox, z_b, eps); }
        else if (strcmp( name, "cpx") == 0 ){
            SS_ref_db  = G_SS_sb21_cpx_function(SS_ref_db, gv.research_group, EM_dataset, gv.len_ox, z_b, eps); }
        else if (strcmp( name, "hpcpx") == 0 ){
            SS_ref_db  = G_SS_sb21_hpcpx_function(SS_ref_db, gv.research_group, EM_dataset, gv.len_ox, z_b, eps);       }
        else if (strcmp( name, "ak") == 0 ){
            SS_ref_db  = G_SS_sb21_ak_function(SS_ref_db, gv.research_group, EM_dataset, gv.len_ox, z_b, eps);  }
        else if (strcmp( name, "gtmj") == 0 ){
            SS_ref_db  = G_SS_sb21_gtmj_function(SS_ref_db, gv.research_group, EM_dataset, gv.len_ox, z_b, eps);        }
        else if (strcmp( name, "pv") == 0 ){
            SS_ref_db  = G_SS_sb21_pv_function(SS_ref_db, gv.research_group, EM_dataset, gv.len_ox, z_b, eps);  }
        else if (strcmp( name, "ppv") == 0 ){
            SS_ref_db  = G_SS_sb21_ppv_function(SS_ref_db, gv.research_group, EM_dataset, gv.len_ox, z_b, eps); }
        else if (strcmp( name, "cf") == 0 ){
            SS_ref_db  = G_SS_sb21_cf_function(SS_ref_db, gv.research_group, EM_dataset, gv.len_ox, z_b, eps);  }
        else if (strcmp( name, "mw") == 0 ){
            SS_ref_db  = G_SS_sb21_mw_function(SS_ref_db, gv.research_group, EM_dataset, gv.len_ox, z_b, eps);  }
        else if (strcmp( name, "nal") == 0 ){
            SS_ref_db  = G_SS_sb21_nal_function(SS_ref_db, gv.research_group, EM_dataset, gv.len_ox, z_b, eps); }
        else{
            printf("\nsolid solution '%s' is not in the database\n",name);      }

        for (int j = 0; j < SS_ref_db.n_em; j++){
            SS_ref_db.mu_array[FD][j] = SS_ref_db.gbase[j];
        }
    }

    for (int j = 0; j < SS_ref_db.n_em; j++){
        SS_ref_db.bounds[j][0] = SS_ref_db.bounds_ref[j][0];
        SS_ref_db.bounds[j][1] = SS_ref_db.bounds_ref[j][1];
    }

    /* Calculate the number of atoms in the bulk-rock composition */
    double fbc     = 0.0;
    for (int i = 0; i < gv.len_ox; i++){
        fbc += z_b.bulk_rock[i]*z_b.apo[i];
    }

    for (int i = 0; i < SS_ref_db.n_em; i++){
        SS_ref_db.ape[i] = 0.0;
        for (int j = 0; j < gv.len_ox; j++){
            SS_ref_db.ape[i] += SS_ref_db.Comp[i][j]*z_b.apo[j];
        }
    }

    SS_ref_db.fbc = z_b.fbc;

    if (gv.verbose == 1){
        printf(" %4s:",name);
        for (int j = 0; j < SS_ref_db.n_em; j++){
            printf(" %+12.5f",SS_ref_db.gbase[j]);
        }
        printf("\n");
        printf(" S   C   A   F   M   N\n");
        for (int i = 0; i < SS_ref_db.n_em; i++){
            for (int j = 0; j < gv.len_ox; j++){
                printf(" %.1f",SS_ref_db.Comp[i][j]);
            }
            printf("\n");
        }
        printf("\n");
    }

    return SS_ref_db;
};