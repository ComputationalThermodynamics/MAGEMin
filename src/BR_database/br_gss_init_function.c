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
#include <string.h>

#include "br_gss_init_function.h"

SS_ref G_SS_br_ctd_init_function(SS_ref SS_ref_db, global_variable gv){
    (void) gv;
    SS_ref_db.is_liq    = 0;
    SS_ref_db.override  = 0;
    SS_ref_db.symmetry  = 1;
    SS_ref_db.n_cat     = 0;
    SS_ref_db.n_xeos    = 2;
    SS_ref_db.n_em      = 2;
    SS_ref_db.n_sf      = 2;
    SS_ref_db.n_w       = 1;

    return SS_ref_db;
}

SS_ref G_SS_br_car_init_function(SS_ref SS_ref_db, global_variable gv){
    (void) gv;
    SS_ref_db.is_liq    = 0;
    SS_ref_db.override  = 0;
    SS_ref_db.symmetry  = 1;
    SS_ref_db.n_cat     = 0;
    SS_ref_db.n_xeos    = 2;
    SS_ref_db.n_em      = 2;
    SS_ref_db.n_sf      = 2;
    SS_ref_db.n_w       = 1;

    return SS_ref_db;
}

SS_ref G_SS_br_chl_init_function(SS_ref SS_ref_db, global_variable gv){
    (void) gv;
    SS_ref_db.is_liq    = 0;
    SS_ref_db.override  = 0;
    SS_ref_db.symmetry  = 1;
    SS_ref_db.n_cat     = 0;
    SS_ref_db.n_xeos    = 5;
    SS_ref_db.n_em      = 5;
    SS_ref_db.n_sf      = 9;
    SS_ref_db.n_w       = 6;

    return SS_ref_db;
}

SS_ref G_SS_br_mica_init_function(SS_ref SS_ref_db, global_variable gv){
    (void) gv;
    SS_ref_db.is_liq    = 0;
    SS_ref_db.override  = 0;
    SS_ref_db.symmetry  = 1;
    SS_ref_db.n_cat     = 0;
    SS_ref_db.n_xeos    = 5;
    SS_ref_db.n_em      = 5;
    SS_ref_db.n_sf      = 8;
    SS_ref_db.n_w       = 9;

    return SS_ref_db;
}

/* Ideal or single-symmetric-W binaries (M(mult):species site), (SITE)/(IDEAL) tag only */
SS_ref G_SS_br_talc_init_function(SS_ref SS_ref_db, global_variable gv){
    (void) gv;
    SS_ref_db.is_liq    = 0;
    SS_ref_db.override  = 0;
    SS_ref_db.symmetry  = 1;
    SS_ref_db.n_cat     = 0;
    SS_ref_db.n_xeos    = 2;
    SS_ref_db.n_em      = 2;
    SS_ref_db.n_sf      = 2;
    SS_ref_db.n_w       = 1;
    return SS_ref_db;
}

SS_ref G_SS_br_ilm_init_function(SS_ref SS_ref_db, global_variable gv){
    (void) gv;
    SS_ref_db.is_liq    = 0;
    SS_ref_db.override  = 0;
    SS_ref_db.symmetry  = 1;
    SS_ref_db.n_cat     = 0;
    SS_ref_db.n_xeos    = 2;
    SS_ref_db.n_em      = 2;
    SS_ref_db.n_sf      = 2;
    SS_ref_db.n_w       = 1;
    return SS_ref_db;
}

SS_ref G_SS_br_bt_init_function(SS_ref SS_ref_db, global_variable gv){
    (void) gv;
    SS_ref_db.is_liq    = 0;
    SS_ref_db.override  = 0;
    SS_ref_db.symmetry  = 1;
    SS_ref_db.n_cat     = 0;
    SS_ref_db.n_xeos    = 2;
    SS_ref_db.n_em      = 2;
    SS_ref_db.n_sf      = 2;
    SS_ref_db.n_w       = 1;
    return SS_ref_db;
}

SS_ref G_SS_br_ol_init_function(SS_ref SS_ref_db, global_variable gv){
    (void) gv;
    SS_ref_db.is_liq    = 0;
    SS_ref_db.override  = 0;
    SS_ref_db.symmetry  = 1;
    SS_ref_db.n_cat     = 0;
    SS_ref_db.n_xeos    = 2;
    SS_ref_db.n_em      = 2;
    SS_ref_db.n_sf      = 2;
    SS_ref_db.n_w       = 1;
    return SS_ref_db;
}

SS_ref G_SS_br_ep_init_function(SS_ref SS_ref_db, global_variable gv){
    (void) gv;
    SS_ref_db.is_liq    = 0;
    SS_ref_db.override  = 0;
    SS_ref_db.symmetry  = 1;
    SS_ref_db.n_cat     = 0;
    SS_ref_db.n_xeos    = 2;
    SS_ref_db.n_em      = 2;
    SS_ref_db.n_sf      = 2;
    SS_ref_db.n_w       = 1;
    return SS_ref_db;
}

SS_ref G_SS_br_opx_init_function(SS_ref SS_ref_db, global_variable gv){
    (void) gv;
    SS_ref_db.is_liq    = 0;
    SS_ref_db.override  = 0;
    SS_ref_db.symmetry  = 1;
    SS_ref_db.n_cat     = 0;
    SS_ref_db.n_xeos    = 2;
    SS_ref_db.n_em      = 2;
    SS_ref_db.n_sf      = 2;
    SS_ref_db.n_w       = 1;
    return SS_ref_db;
}

SS_ref G_SS_br_amph_init_function(SS_ref SS_ref_db, global_variable gv){
    (void) gv;
    SS_ref_db.is_liq    = 0;
    SS_ref_db.override  = 0;
    SS_ref_db.symmetry  = 1;
    SS_ref_db.n_cat     = 0;
    SS_ref_db.n_xeos    = 2;
    SS_ref_db.n_em      = 2;
    SS_ref_db.n_sf      = 2;
    SS_ref_db.n_w       = 1;
    return SS_ref_db;
}

SS_ref G_SS_br_spl_init_function(SS_ref SS_ref_db, global_variable gv){
    (void) gv;
    SS_ref_db.is_liq    = 0;
    SS_ref_db.override  = 0;
    SS_ref_db.symmetry  = 1;
    SS_ref_db.n_cat     = 0;
    SS_ref_db.n_xeos    = 2;
    SS_ref_db.n_em      = 2;
    SS_ref_db.n_sf      = 2;
    SS_ref_db.n_w       = 1;
    return SS_ref_db;
}

SS_ref G_SS_br_stau_init_function(SS_ref SS_ref_db, global_variable gv){
    (void) gv;
    SS_ref_db.is_liq    = 0;
    SS_ref_db.override  = 0;
    SS_ref_db.symmetry  = 1;
    SS_ref_db.n_cat     = 0;
    SS_ref_db.n_xeos    = 2;
    SS_ref_db.n_em      = 2;
    SS_ref_db.n_sf      = 2;
    SS_ref_db.n_w       = 1;
    return SS_ref_db;
}

SS_ref G_SS_br_crd_init_function(SS_ref SS_ref_db, global_variable gv){
    (void) gv;
    SS_ref_db.is_liq    = 0;
    SS_ref_db.override  = 0;
    SS_ref_db.symmetry  = 1;
    SS_ref_db.n_cat     = 0;
    SS_ref_db.n_xeos    = 2;
    SS_ref_db.n_em      = 2;
    SS_ref_db.n_sf      = 2;
    SS_ref_db.n_w       = 1;
    return SS_ref_db;
}

/* Asymmetric-subregular binary (MARGULES,SITE), 2 W terms */
SS_ref G_SS_br_grt_init_function(SS_ref SS_ref_db, global_variable gv){
    (void) gv;
    SS_ref_db.is_liq    = 0;
    SS_ref_db.override  = 0;
    SS_ref_db.symmetry  = 1;
    SS_ref_db.n_cat     = 0;
    SS_ref_db.n_xeos    = 2;
    SS_ref_db.n_em      = 2;
    SS_ref_db.n_sf      = 2;
    SS_ref_db.n_w       = 2;
    return SS_ref_db;
}

/* Ternary ideal + subregular Margules, no site-splitting (MARGULES,IDEAL) */
SS_ref G_SS_br_omph_init_function(SS_ref SS_ref_db, global_variable gv){
    (void) gv;
    SS_ref_db.is_liq    = 0;
    SS_ref_db.override  = 0;
    SS_ref_db.symmetry  = 1;
    SS_ref_db.n_cat     = 0;
    SS_ref_db.n_xeos    = 3;
    SS_ref_db.n_em      = 3;
    SS_ref_db.n_sf      = 3;
    SS_ref_db.n_w       = 6;
    return SS_ref_db;
}

/* Reciprocal 2-site entropy (A(1):v,Na - M1(2):Mg,Al) + symmetric ternary
   Margules on raw endmember fractions (IDEAL,MARGULES) */
SS_ref G_SS_br_amphx_init_function(SS_ref SS_ref_db, global_variable gv){
    (void) gv;
    SS_ref_db.is_liq    = 0;
    SS_ref_db.override  = 0;
    SS_ref_db.symmetry  = 1;
    SS_ref_db.n_cat     = 0;
    SS_ref_db.n_xeos    = 3;
    SS_ref_db.n_em      = 3;
    SS_ref_db.n_sf      = 4;
    SS_ref_db.n_w       = 3;
    return SS_ref_db;
}

/* Feldspar: (EXT,MARGULES), ternary ideal entropy (no site given, mult=1)
   + 3-pair subregular + ternary Margules (Fuhrman & Lindsley 1988 - the
   only one of Pourt14's 3 competing FSP calibrations under a correctly-
   spelled MARGULES PARAMETER header; Benisek 2004 and Ghiorso 1984 are
   both deliberately deactivated, same misspelled-keyword convention as
   CPX/P2n). Al-Si order-disorder (D1/D2) intentionally dropped - fixed-G0
   endmembers, same simplification decided for this dataset. */
SS_ref G_SS_br_fsp_init_function(SS_ref SS_ref_db, global_variable gv){
    (void) gv;
    SS_ref_db.is_liq    = 0;
    SS_ref_db.override  = 0;
    SS_ref_db.symmetry  = 1;
    SS_ref_db.n_cat     = 0;
    SS_ref_db.n_xeos    = 3;
    SS_ref_db.n_em      = 3;
    SS_ref_db.n_sf      = 3;
    SS_ref_db.n_w       = 7;
    return SS_ref_db;
}

void BR_SS_init(SS_init_type *SS_init, global_variable gv){
    for (int iss = 0; iss < gv.len_ss; iss++){
        if (strcmp( gv.SS_list[iss], "ctd") == 0 ){
            SS_init[iss] = G_SS_br_ctd_init_function;
        }
        else if (strcmp( gv.SS_list[iss], "car") == 0 ){
            SS_init[iss] = G_SS_br_car_init_function;
        }
        else if (strcmp( gv.SS_list[iss], "chl") == 0 ){
            SS_init[iss] = G_SS_br_chl_init_function;
        }
        else if (strcmp( gv.SS_list[iss], "mica") == 0 ){
            SS_init[iss] = G_SS_br_mica_init_function;
        }
        else if (strcmp( gv.SS_list[iss], "talc") == 0 ){
            SS_init[iss] = G_SS_br_talc_init_function;
        }
        else if (strcmp( gv.SS_list[iss], "ilm") == 0 ){
            SS_init[iss] = G_SS_br_ilm_init_function;
        }
        else if (strcmp( gv.SS_list[iss], "bt") == 0 ){
            SS_init[iss] = G_SS_br_bt_init_function;
        }
        else if (strcmp( gv.SS_list[iss], "ol") == 0 ){
            SS_init[iss] = G_SS_br_ol_init_function;
        }
        else if (strcmp( gv.SS_list[iss], "ep") == 0 ){
            SS_init[iss] = G_SS_br_ep_init_function;
        }
        else if (strcmp( gv.SS_list[iss], "opx") == 0 ){
            SS_init[iss] = G_SS_br_opx_init_function;
        }
        else if (strcmp( gv.SS_list[iss], "amph") == 0 ){
            SS_init[iss] = G_SS_br_amph_init_function;
        }
        else if (strcmp( gv.SS_list[iss], "spl") == 0 ){
            SS_init[iss] = G_SS_br_spl_init_function;
        }
        else if (strcmp( gv.SS_list[iss], "stau") == 0 ){
            SS_init[iss] = G_SS_br_stau_init_function;
        }
        else if (strcmp( gv.SS_list[iss], "crd") == 0 ){
            SS_init[iss] = G_SS_br_crd_init_function;
        }
        else if (strcmp( gv.SS_list[iss], "grt") == 0 ){
            SS_init[iss] = G_SS_br_grt_init_function;
        }
        else if (strcmp( gv.SS_list[iss], "omph") == 0 ){
            SS_init[iss] = G_SS_br_omph_init_function;
        }
        else if (strcmp( gv.SS_list[iss], "amphx") == 0 ){
            SS_init[iss] = G_SS_br_amphx_init_function;
        }
        else if (strcmp( gv.SS_list[iss], "fsp") == 0 ){
            SS_init[iss] = G_SS_br_fsp_init_function;
        }
    }
}
