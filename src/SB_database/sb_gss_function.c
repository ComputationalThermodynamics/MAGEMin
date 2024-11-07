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

    char   *EM_tmp[] 		= {"an","ab",};
    for (int i = 0; i < SS_ref_db.n_em; i++){
        strcpy(SS_ref_db.EM_list[i],EM_tmp[i]);
    };

    SS_ref_db.W[0] = 26000.0;


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

    for (i = 0; i < len_ox; i++){
        SS_ref_db.Comp[0][i] 	= an.C[i];
        SS_ref_db.Comp[1][i] 	= ab.C[i];
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

    char   *EM_tmp[] 		= {"sp","hc",};
    for (int i = 0; i < SS_ref_db.n_em; i++){
        strcpy(SS_ref_db.EM_list[i],EM_tmp[i]);
    };

    SS_ref_db.W[0] = 5000.0;


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

    for (i = 0; i < len_ox; i++){
        SS_ref_db.Comp[0][i] 	= sp.C[i];
        SS_ref_db.Comp[1][i] 	= hc.C[i];
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

    char   *EM_tmp[] 		= {"fa","fo",};
    for (int i = 0; i < SS_ref_db.n_em; i++){
        strcpy(SS_ref_db.EM_list[i],EM_tmp[i]);
    };

    SS_ref_db.W[0] = 7600.0;


    em_data fa            = get_em_data(    research_group, EM_dataset,
                                            len_ox,
                                            z_b,
                                            SS_ref_db.P,
                                            SS_ref_db.T,
                                            "fa",
                                            "equilibrium");

    em_data fo            = get_em_data(    research_group, EM_dataset,
                                            len_ox,
                                            z_b,
                                            SS_ref_db.P,
                                            SS_ref_db.T,
                                            "fo",
                                            "equilibrium");

    SS_ref_db.gbase[0]        = fa.gb;
    SS_ref_db.gbase[1]        = fo.gb;

    SS_ref_db.ElShearMod[0]        = fa.ElShearMod;
    SS_ref_db.ElShearMod[1]        = fo.ElShearMod;

    for (i = 0; i < len_ox; i++){
        SS_ref_db.Comp[0][i] 	= fa.C[i];
        SS_ref_db.Comp[1][i] 	= fo.C[i];
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

    char   *EM_tmp[] 		= {"fewa","mgwa",};
    for (int i = 0; i < SS_ref_db.n_em; i++){
        strcpy(SS_ref_db.EM_list[i],EM_tmp[i]);
    };

    SS_ref_db.W[0] = 16500.0;


    em_data fewa            = get_em_data(    research_group, EM_dataset,
                                            len_ox,
                                            z_b,
                                            SS_ref_db.P,
                                            SS_ref_db.T,
                                            "fewa",
                                            "equilibrium");

    em_data mgwa            = get_em_data(    research_group, EM_dataset,
                                            len_ox,
                                            z_b,
                                            SS_ref_db.P,
                                            SS_ref_db.T,
                                            "mgwa",
                                            "equilibrium");

    SS_ref_db.gbase[0]        = fewa.gb;
    SS_ref_db.gbase[1]        = mgwa.gb;

    SS_ref_db.ElShearMod[0]        = fewa.ElShearMod;
    SS_ref_db.ElShearMod[1]        = mgwa.ElShearMod;

    for (i = 0; i < len_ox; i++){
        SS_ref_db.Comp[0][i] 	= fewa.C[i];
        SS_ref_db.Comp[1][i] 	= mgwa.C[i];
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

    char   *EM_tmp[] 		= {"mgri","feri",};
    for (int i = 0; i < SS_ref_db.n_em; i++){
        strcpy(SS_ref_db.EM_list[i],EM_tmp[i]);
    };

    SS_ref_db.W[0] = 9100.0;


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

    for (i = 0; i < len_ox; i++){
        SS_ref_db.Comp[0][i] 	= mgri.C[i];
        SS_ref_db.Comp[1][i] 	= feri.C[i];
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

    char   *EM_tmp[] 		= {"mgts","fs","en","odi",};
    for (int i = 0; i < SS_ref_db.n_em; i++){
        strcpy(SS_ref_db.EM_list[i],EM_tmp[i]);
    };

    SS_ref_db.W[0] = 0.0;
    SS_ref_db.W[1] = 0.0;
    SS_ref_db.W[2] = 32100.0;
    SS_ref_db.W[3] = 0.0;
    SS_ref_db.W[4] = 0.0;
    SS_ref_db.W[5] = 48000.0;


    em_data mgts            = get_em_data(    research_group, EM_dataset,
                                            len_ox,
                                            z_b,
                                            SS_ref_db.P,
                                            SS_ref_db.T,
                                            "mgts",
                                            "equilibrium");

    em_data fs            = get_em_data(    research_group, EM_dataset,
                                            len_ox,
                                            z_b,
                                            SS_ref_db.P,
                                            SS_ref_db.T,
                                            "fs",
                                            "equilibrium");

    em_data en            = get_em_data(    research_group, EM_dataset,
                                            len_ox,
                                            z_b,
                                            SS_ref_db.P,
                                            SS_ref_db.T,
                                            "en",
                                            "equilibrium");

    em_data odi            = get_em_data(    research_group, EM_dataset,
                                            len_ox,
                                            z_b,
                                            SS_ref_db.P,
                                            SS_ref_db.T,
                                            "odi",
                                            "equilibrium");

    SS_ref_db.gbase[0]        = mgts.gb;
    SS_ref_db.gbase[1]        = fs.gb;
    SS_ref_db.gbase[2]        = en.gb;
    SS_ref_db.gbase[3]        = odi.gb;

    SS_ref_db.ElShearMod[0]        = mgts.ElShearMod;
    SS_ref_db.ElShearMod[1]        = fs.ElShearMod;
    SS_ref_db.ElShearMod[2]        = en.ElShearMod;
    SS_ref_db.ElShearMod[3]        = odi.ElShearMod;

    for (i = 0; i < len_ox; i++){
        SS_ref_db.Comp[0][i] 	= mgts.C[i];
        SS_ref_db.Comp[1][i] 	= fs.C[i];
        SS_ref_db.Comp[2][i] 	= en.C[i];
        SS_ref_db.Comp[3][i] 	= odi.C[i];
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

    char   *EM_tmp[] 		= {"he","jd","cen","cats","di",};
    for (int i = 0; i < SS_ref_db.n_em; i++){
        strcpy(SS_ref_db.EM_list[i],EM_tmp[i]);
    };

    SS_ref_db.W[0] = 24700.0;
    SS_ref_db.W[1] = 0.0;
    SS_ref_db.W[2] = 26000.0;
    SS_ref_db.W[3] = 24700.0;
    SS_ref_db.W[4] = 24300.0;
    SS_ref_db.W[5] = 0.0;
    SS_ref_db.W[6] = 60600.0;
    SS_ref_db.W[7] = 10000.0;
    SS_ref_db.W[8] = 0.0;
    SS_ref_db.W[9] = 0.0;

    SS_ref_db.v[0] = 1.0;
    SS_ref_db.v[1] = 1.0;
    SS_ref_db.v[2] = 1.0;
    SS_ref_db.v[3] = 3.5;
    SS_ref_db.v[4] = 1.0;

    em_data he            = get_em_data(    research_group, EM_dataset,
                                            len_ox,
                                            z_b,
                                            SS_ref_db.P,
                                            SS_ref_db.T,
                                            "he",
                                            "equilibrium");

    em_data jd            = get_em_data(    research_group, EM_dataset,
                                            len_ox,
                                            z_b,
                                            SS_ref_db.P,
                                            SS_ref_db.T,
                                            "jd",
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

    em_data di            = get_em_data(    research_group, EM_dataset,
                                            len_ox,
                                            z_b,
                                            SS_ref_db.P,
                                            SS_ref_db.T,
                                            "di",
                                            "equilibrium");

    SS_ref_db.gbase[0]        = he.gb;
    SS_ref_db.gbase[1]        = jd.gb;
    SS_ref_db.gbase[2]        = cen.gb;
    SS_ref_db.gbase[3]        = cats.gb;
    SS_ref_db.gbase[4]        = di.gb;

    SS_ref_db.ElShearMod[0]        = he.ElShearMod;
    SS_ref_db.ElShearMod[1]        = jd.ElShearMod;
    SS_ref_db.ElShearMod[2]        = cen.ElShearMod;
    SS_ref_db.ElShearMod[3]        = cats.ElShearMod;
    SS_ref_db.ElShearMod[4]        = di.ElShearMod;

    for (i = 0; i < len_ox; i++){
        SS_ref_db.Comp[0][i] 	= he.C[i];
        SS_ref_db.Comp[1][i] 	= jd.C[i];
        SS_ref_db.Comp[2][i] 	= cen.C[i];
        SS_ref_db.Comp[3][i] 	= cats.C[i];
        SS_ref_db.Comp[4][i] 	= di.C[i];
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

    char   *EM_tmp[] 		= {"hpcen","hpcfs",};
    for (int i = 0; i < SS_ref_db.n_em; i++){
        strcpy(SS_ref_db.EM_list[i],EM_tmp[i]);
    };

    SS_ref_db.W[0] = 0.0;


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

    for (i = 0; i < len_ox; i++){
        SS_ref_db.Comp[0][i] 	= hpcen.C[i];
        SS_ref_db.Comp[1][i] 	= hpcfs.C[i];
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

    char   *EM_tmp[] 		= {"co","mgak","feak",};
    for (int i = 0; i < SS_ref_db.n_em; i++){
        strcpy(SS_ref_db.EM_list[i],EM_tmp[i]);
    };

    SS_ref_db.W[0] = 66000.0;
    SS_ref_db.W[1] = 0.0;
    SS_ref_db.W[2] = 0.0;


    em_data co            = get_em_data(    research_group, EM_dataset,
                                            len_ox,
                                            z_b,
                                            SS_ref_db.P,
                                            SS_ref_db.T,
                                            "co",
                                            "equilibrium");

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

    SS_ref_db.gbase[0]        = co.gb;
    SS_ref_db.gbase[1]        = mgak.gb;
    SS_ref_db.gbase[2]        = feak.gb;

    SS_ref_db.ElShearMod[0]        = co.ElShearMod;
    SS_ref_db.ElShearMod[1]        = mgak.ElShearMod;
    SS_ref_db.ElShearMod[2]        = feak.ElShearMod;

    for (i = 0; i < len_ox; i++){
        SS_ref_db.Comp[0][i] 	= co.C[i];
        SS_ref_db.Comp[1][i] 	= mgak.C[i];
        SS_ref_db.Comp[2][i] 	= feak.C[i];
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

    char   *EM_tmp[] 		= {"alm","jdmj","mgmj","py","gr",};
    for (int i = 0; i < SS_ref_db.n_em; i++){
        strcpy(SS_ref_db.EM_list[i],EM_tmp[i]);
    };

    SS_ref_db.W[0] = 30000.0;
    SS_ref_db.W[1] = 0.0;
    SS_ref_db.W[2] = 0.0;
    SS_ref_db.W[3] = 0.0;
    SS_ref_db.W[4] = 0.0;
    SS_ref_db.W[5] = 0.0;
    SS_ref_db.W[6] = 21300.0;
    SS_ref_db.W[7] = 0.0;
    SS_ref_db.W[8] = 0.0;
    SS_ref_db.W[9] = 58000.0;


    em_data alm            = get_em_data(    research_group, EM_dataset,
                                            len_ox,
                                            z_b,
                                            SS_ref_db.P,
                                            SS_ref_db.T,
                                            "alm",
                                            "equilibrium");

    em_data jdmj            = get_em_data(    research_group, EM_dataset,
                                            len_ox,
                                            z_b,
                                            SS_ref_db.P,
                                            SS_ref_db.T,
                                            "jdmj",
                                            "equilibrium");

    em_data mgmj            = get_em_data(    research_group, EM_dataset,
                                            len_ox,
                                            z_b,
                                            SS_ref_db.P,
                                            SS_ref_db.T,
                                            "mgmj",
                                            "equilibrium");

    em_data py            = get_em_data(    research_group, EM_dataset,
                                            len_ox,
                                            z_b,
                                            SS_ref_db.P,
                                            SS_ref_db.T,
                                            "py",
                                            "equilibrium");

    em_data gr            = get_em_data(    research_group, EM_dataset,
                                            len_ox,
                                            z_b,
                                            SS_ref_db.P,
                                            SS_ref_db.T,
                                            "gr",
                                            "equilibrium");

    SS_ref_db.gbase[0]        = alm.gb;
    SS_ref_db.gbase[1]        = jdmj.gb;
    SS_ref_db.gbase[2]        = mgmj.gb;
    SS_ref_db.gbase[3]        = py.gb;
    SS_ref_db.gbase[4]        = gr.gb;

    SS_ref_db.ElShearMod[0]        = alm.ElShearMod;
    SS_ref_db.ElShearMod[1]        = jdmj.ElShearMod;
    SS_ref_db.ElShearMod[2]        = mgmj.ElShearMod;
    SS_ref_db.ElShearMod[3]        = py.ElShearMod;
    SS_ref_db.ElShearMod[4]        = gr.ElShearMod;

    for (i = 0; i < len_ox; i++){
        SS_ref_db.Comp[0][i] 	= alm.C[i];
        SS_ref_db.Comp[1][i] 	= jdmj.C[i];
        SS_ref_db.Comp[2][i] 	= mgmj.C[i];
        SS_ref_db.Comp[3][i] 	= py.C[i];
        SS_ref_db.Comp[4][i] 	= gr.C[i];
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

    char   *EM_tmp[] 		= {"alpv","fepv","mgpv",};
    for (int i = 0; i < SS_ref_db.n_em; i++){
        strcpy(SS_ref_db.EM_list[i],EM_tmp[i]);
    };

    SS_ref_db.W[0] = 0.0;
    SS_ref_db.W[1] = 116000.0;
    SS_ref_db.W[2] = 0.0;

    SS_ref_db.v[0] = 1.0;
    SS_ref_db.v[1] = 1.0;
    SS_ref_db.v[2] = 0.39;

    em_data alpv            = get_em_data(    research_group, EM_dataset,
                                            len_ox,
                                            z_b,
                                            SS_ref_db.P,
                                            SS_ref_db.T,
                                            "alpv",
                                            "equilibrium");

    em_data fepv            = get_em_data(    research_group, EM_dataset,
                                            len_ox,
                                            z_b,
                                            SS_ref_db.P,
                                            SS_ref_db.T,
                                            "fepv",
                                            "equilibrium");

    em_data mgpv            = get_em_data(    research_group, EM_dataset,
                                            len_ox,
                                            z_b,
                                            SS_ref_db.P,
                                            SS_ref_db.T,
                                            "mgpv",
                                            "equilibrium");

    SS_ref_db.gbase[0]        = alpv.gb;
    SS_ref_db.gbase[1]        = fepv.gb;
    SS_ref_db.gbase[2]        = mgpv.gb;

    SS_ref_db.ElShearMod[0]        = alpv.ElShearMod;
    SS_ref_db.ElShearMod[1]        = fepv.ElShearMod;
    SS_ref_db.ElShearMod[2]        = mgpv.ElShearMod;

    for (i = 0; i < len_ox; i++){
        SS_ref_db.Comp[0][i] 	= alpv.C[i];
        SS_ref_db.Comp[1][i] 	= fepv.C[i];
        SS_ref_db.Comp[2][i] 	= mgpv.C[i];
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

    char   *EM_tmp[] 		= {"feppv","mgppv","alppv",};
    for (int i = 0; i < SS_ref_db.n_em; i++){
        strcpy(SS_ref_db.EM_list[i],EM_tmp[i]);
    };

    SS_ref_db.W[0] = 0.0;
    SS_ref_db.W[1] = 0.0;
    SS_ref_db.W[2] = 60000.0;


    em_data feppv            = get_em_data(    research_group, EM_dataset,
                                            len_ox,
                                            z_b,
                                            SS_ref_db.P,
                                            SS_ref_db.T,
                                            "feppv",
                                            "equilibrium");

    em_data mgppv            = get_em_data(    research_group, EM_dataset,
                                            len_ox,
                                            z_b,
                                            SS_ref_db.P,
                                            SS_ref_db.T,
                                            "mgppv",
                                            "equilibrium");

    em_data alppv            = get_em_data(    research_group, EM_dataset,
                                            len_ox,
                                            z_b,
                                            SS_ref_db.P,
                                            SS_ref_db.T,
                                            "alppv",
                                            "equilibrium");

    SS_ref_db.gbase[0]        = feppv.gb;
    SS_ref_db.gbase[1]        = mgppv.gb;
    SS_ref_db.gbase[2]        = alppv.gb;

    SS_ref_db.ElShearMod[0]        = feppv.ElShearMod;
    SS_ref_db.ElShearMod[1]        = mgppv.ElShearMod;
    SS_ref_db.ElShearMod[2]        = alppv.ElShearMod;

    for (i = 0; i < len_ox; i++){
        SS_ref_db.Comp[0][i] 	= feppv.C[i];
        SS_ref_db.Comp[1][i] 	= mgppv.C[i];
        SS_ref_db.Comp[2][i] 	= alppv.C[i];
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

    char   *EM_tmp[] 		= {"wu","pe",};
    for (int i = 0; i < SS_ref_db.n_em; i++){
        strcpy(SS_ref_db.EM_list[i],EM_tmp[i]);
    };

    SS_ref_db.W[0] = 13000.0;


    em_data wu            = get_em_data(    research_group, EM_dataset,
                                            len_ox,
                                            z_b,
                                            SS_ref_db.P,
                                            SS_ref_db.T,
                                            "wu",
                                            "equilibrium");

    em_data pe            = get_em_data(    research_group, EM_dataset,
                                            len_ox,
                                            z_b,
                                            SS_ref_db.P,
                                            SS_ref_db.T,
                                            "pe",
                                            "equilibrium");

    SS_ref_db.gbase[0]        = wu.gb;
    SS_ref_db.gbase[1]        = pe.gb;

    SS_ref_db.ElShearMod[0]        = wu.ElShearMod;
    SS_ref_db.ElShearMod[1]        = pe.ElShearMod;

    for (i = 0; i < len_ox; i++){
        SS_ref_db.Comp[0][i] 	= wu.C[i];
        SS_ref_db.Comp[1][i] 	= pe.C[i];
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

    char   *EM_tmp[] 		= {"mgcf","nacf","fecf",};
    for (int i = 0; i < SS_ref_db.n_em; i++){
        strcpy(SS_ref_db.EM_list[i],EM_tmp[i]);
    };

    SS_ref_db.W[0] = 0.0;
    SS_ref_db.W[1] = 0.0;
    SS_ref_db.W[2] = 0.0;


    em_data mgcf            = get_em_data(    research_group, EM_dataset,
                                            len_ox,
                                            z_b,
                                            SS_ref_db.P,
                                            SS_ref_db.T,
                                            "mgcf",
                                            "equilibrium");

    em_data nacf            = get_em_data(    research_group, EM_dataset,
                                            len_ox,
                                            z_b,
                                            SS_ref_db.P,
                                            SS_ref_db.T,
                                            "nacf",
                                            "equilibrium");

    em_data fecf            = get_em_data(    research_group, EM_dataset,
                                            len_ox,
                                            z_b,
                                            SS_ref_db.P,
                                            SS_ref_db.T,
                                            "fecf",
                                            "equilibrium");

    SS_ref_db.gbase[0]        = mgcf.gb;
    SS_ref_db.gbase[1]        = nacf.gb;
    SS_ref_db.gbase[2]        = fecf.gb;

    SS_ref_db.ElShearMod[0]        = mgcf.ElShearMod;
    SS_ref_db.ElShearMod[1]        = nacf.ElShearMod;
    SS_ref_db.ElShearMod[2]        = fecf.ElShearMod;

    for (i = 0; i < len_ox; i++){
        SS_ref_db.Comp[0][i] 	= mgcf.C[i];
        SS_ref_db.Comp[1][i] 	= nacf.C[i];
        SS_ref_db.Comp[2][i] 	= fecf.C[i];
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
