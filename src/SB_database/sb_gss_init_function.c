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
  Function to allocate memory for solid-solutions        
*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <complex.h> 

#include "../MAGEMin.h"
#include "../all_solution_phases.h"
#include "sb_gss_init_function.h"
/**
    allocate memory for sb11_plg
*/
SS_ref G_SS_sb11_plg_init_function(SS_ref SS_ref_db,  global_variable gv){

    SS_ref_db.is_liq    = 0;
    SS_ref_db.override  = 0;
    SS_ref_db.symmetry  = 1;
    SS_ref_db.n_cat     = 2;
    SS_ref_db.n_xeos    = 2;
    SS_ref_db.n_em      = 2;
    SS_ref_db.n_sf      = 1;
    SS_ref_db.n_w       = 1;

     return SS_ref_db;
}

/**
    allocate memory for sb11_sp
*/
SS_ref G_SS_sb11_sp_init_function(SS_ref SS_ref_db,  global_variable gv){

    SS_ref_db.is_liq    = 0;
    SS_ref_db.override  = 0;
    SS_ref_db.symmetry  = 1;
    SS_ref_db.n_cat     = 6;
    SS_ref_db.n_xeos    = 2;
    SS_ref_db.n_em      = 2;
    SS_ref_db.n_sf      = 2;
    SS_ref_db.n_w       = 1;

     return SS_ref_db;
}

/**
    allocate memory for sb11_ol
*/
SS_ref G_SS_sb11_ol_init_function(SS_ref SS_ref_db,  global_variable gv){

    SS_ref_db.is_liq    = 0;
    SS_ref_db.override  = 0;
    SS_ref_db.symmetry  = 1;
    SS_ref_db.n_cat     = 2;
    SS_ref_db.n_xeos    = 2;
    SS_ref_db.n_em      = 2;
    SS_ref_db.n_sf      = 1;
    SS_ref_db.n_w       = 1;

     return SS_ref_db;
}

/**
    allocate memory for sb11_wa
*/
SS_ref G_SS_sb11_wa_init_function(SS_ref SS_ref_db,  global_variable gv){

    SS_ref_db.is_liq    = 0;
    SS_ref_db.override  = 0;
    SS_ref_db.symmetry  = 1;
    SS_ref_db.n_cat     = 2;
    SS_ref_db.n_xeos    = 2;
    SS_ref_db.n_em      = 2;
    SS_ref_db.n_sf      = 1;
    SS_ref_db.n_w       = 1;

     return SS_ref_db;
}

/**
    allocate memory for sb11_ri
*/
SS_ref G_SS_sb11_ri_init_function(SS_ref SS_ref_db,  global_variable gv){

    SS_ref_db.is_liq    = 0;
    SS_ref_db.override  = 0;
    SS_ref_db.symmetry  = 1;
    SS_ref_db.n_cat     = 2;
    SS_ref_db.n_xeos    = 2;
    SS_ref_db.n_em      = 2;
    SS_ref_db.n_sf      = 1;
    SS_ref_db.n_w       = 1;

     return SS_ref_db;
}

/**
    allocate memory for sb11_opx
*/
SS_ref G_SS_sb11_opx_init_function(SS_ref SS_ref_db,  global_variable gv){

    SS_ref_db.is_liq    = 0;
    SS_ref_db.override  = 0;
    SS_ref_db.symmetry  = 1;
    SS_ref_db.n_cat     = 6;
    SS_ref_db.n_xeos    = 4;
    SS_ref_db.n_em      = 4;
    SS_ref_db.n_sf      = 2;
    SS_ref_db.n_w       = 6;

     return SS_ref_db;
}

/**
    allocate memory for sb11_cpx
*/
SS_ref G_SS_sb11_cpx_init_function(SS_ref SS_ref_db,  global_variable gv){

    SS_ref_db.is_liq    = 0;
    SS_ref_db.override  = 0;
    SS_ref_db.symmetry  = 0;
    SS_ref_db.n_cat     = 8;
    SS_ref_db.n_xeos    = 5;
    SS_ref_db.n_em      = 5;
    SS_ref_db.n_sf      = 3;
    SS_ref_db.n_w       = 10;
    SS_ref_db.n_v       = 5;

     return SS_ref_db;
}

/**
    allocate memory for sb11_hpcpx
*/
SS_ref G_SS_sb11_hpcpx_init_function(SS_ref SS_ref_db,  global_variable gv){

    SS_ref_db.is_liq    = 0;
    SS_ref_db.override  = 0;
    SS_ref_db.symmetry  = 1;
    SS_ref_db.n_cat     = 2;
    SS_ref_db.n_xeos    = 2;
    SS_ref_db.n_em      = 2;
    SS_ref_db.n_sf      = 1;
    SS_ref_db.n_w       = 1;

     return SS_ref_db;
}

/**
    allocate memory for sb11_ak
*/
SS_ref G_SS_sb11_ak_init_function(SS_ref SS_ref_db,  global_variable gv){

    SS_ref_db.is_liq    = 0;
    SS_ref_db.override  = 0;
    SS_ref_db.symmetry  = 1;
    SS_ref_db.n_cat     = 5;
    SS_ref_db.n_xeos    = 3;
    SS_ref_db.n_em      = 3;
    SS_ref_db.n_sf      = 2;
    SS_ref_db.n_w       = 3;

     return SS_ref_db;
}

/**
    allocate memory for sb11_gtmj
*/
SS_ref G_SS_sb11_gtmj_init_function(SS_ref SS_ref_db,  global_variable gv){

    SS_ref_db.is_liq    = 0;
    SS_ref_db.override  = 0;
    SS_ref_db.symmetry  = 1;
    SS_ref_db.n_cat     = 9;
    SS_ref_db.n_xeos    = 5;
    SS_ref_db.n_em      = 5;
    SS_ref_db.n_sf      = 3;
    SS_ref_db.n_w       = 10;

     return SS_ref_db;
}

/**
    allocate memory for sb11_pv
*/
SS_ref G_SS_sb11_pv_init_function(SS_ref SS_ref_db,  global_variable gv){

    SS_ref_db.is_liq    = 0;
    SS_ref_db.override  = 0;
    SS_ref_db.symmetry  = 0;
    SS_ref_db.n_cat     = 5;
    SS_ref_db.n_xeos    = 3;
    SS_ref_db.n_em      = 3;
    SS_ref_db.n_sf      = 2;
    SS_ref_db.n_w       = 3;
    SS_ref_db.n_v       = 3;

     return SS_ref_db;
}

/**
    allocate memory for sb11_ppv
*/
SS_ref G_SS_sb11_ppv_init_function(SS_ref SS_ref_db,  global_variable gv){

    SS_ref_db.is_liq    = 0;
    SS_ref_db.override  = 0;
    SS_ref_db.symmetry  = 1;
    SS_ref_db.n_cat     = 5;
    SS_ref_db.n_xeos    = 3;
    SS_ref_db.n_em      = 3;
    SS_ref_db.n_sf      = 2;
    SS_ref_db.n_w       = 3;

     return SS_ref_db;
}

/**
    allocate memory for sb11_mw
*/
SS_ref G_SS_sb11_mw_init_function(SS_ref SS_ref_db,  global_variable gv){

    SS_ref_db.is_liq    = 0;
    SS_ref_db.override  = 0;
    SS_ref_db.symmetry  = 1;
    SS_ref_db.n_cat     = 2;
    SS_ref_db.n_xeos    = 2;
    SS_ref_db.n_em      = 2;
    SS_ref_db.n_sf      = 1;
    SS_ref_db.n_w       = 1;

     return SS_ref_db;
}

/**
    allocate memory for sb11_cf
*/
SS_ref G_SS_sb11_cf_init_function(SS_ref SS_ref_db,  global_variable gv){

    SS_ref_db.is_liq    = 0;
    SS_ref_db.override  = 0;
    SS_ref_db.symmetry  = 1;
    SS_ref_db.n_cat     = 5;
    SS_ref_db.n_xeos    = 3;
    SS_ref_db.n_em      = 3;
    SS_ref_db.n_sf      = 2;
    SS_ref_db.n_w       = 3;

     return SS_ref_db;
}

void SB_SS_init_sb11(       SS_init_type        *SS_init,
                            global_variable      gv                             ){

    for (int iss = 0; iss < gv.len_ss; iss++){
        if      (strcmp( gv.SS_list[iss], "plg")  == 0 ){
            SS_init[iss]  = G_SS_sb11_plg_init_function;                }
        else if (strcmp( gv.SS_list[iss], "sp")  == 0 ){
            SS_init[iss]  = G_SS_sb11_sp_init_function;                 }
        else if (strcmp( gv.SS_list[iss], "ol")  == 0 ){
            SS_init[iss]  = G_SS_sb11_ol_init_function;                 }
        else if (strcmp( gv.SS_list[iss], "wa")  == 0 ){
            SS_init[iss]  = G_SS_sb11_wa_init_function;                 }
        else if (strcmp( gv.SS_list[iss], "ri")  == 0 ){
            SS_init[iss]  = G_SS_sb11_ri_init_function;                 }
        else if (strcmp( gv.SS_list[iss], "opx")  == 0 ){
            SS_init[iss]  = G_SS_sb11_opx_init_function;                }
        else if (strcmp( gv.SS_list[iss], "cpx")  == 0 ){
            SS_init[iss]  = G_SS_sb11_cpx_init_function;                }
        else if (strcmp( gv.SS_list[iss], "hpcpx")  == 0 ){
            SS_init[iss]  = G_SS_sb11_hpcpx_init_function;              }
        else if (strcmp( gv.SS_list[iss], "ak")  == 0 ){
            SS_init[iss]  = G_SS_sb11_ak_init_function;                 }
        else if (strcmp( gv.SS_list[iss], "gtmj")  == 0 ){
            SS_init[iss]  = G_SS_sb11_gtmj_init_function;               }
        else if (strcmp( gv.SS_list[iss], "pv")  == 0 ){
            SS_init[iss]  = G_SS_sb11_pv_init_function;                 }
        else if (strcmp( gv.SS_list[iss], "ppv")  == 0 ){
            SS_init[iss]  = G_SS_sb11_ppv_init_function;                }
        else if (strcmp( gv.SS_list[iss], "mw")  == 0 ){
            SS_init[iss]  = G_SS_sb11_mw_init_function;                 }
        else if (strcmp( gv.SS_list[iss], "cf")  == 0 ){
            SS_init[iss]  = G_SS_sb11_cf_init_function;                 }
        else{
            printf("\nsolid solution '%s' is not in the database, cannot be initiated\n", gv.SS_list[iss]);
        }
    }
}

/**
    allocate memory for sb21_plg
*/
SS_ref G_SS_sb21_plg_init_function(SS_ref SS_ref_db,  global_variable gv){

    SS_ref_db.is_liq    = 0;
    SS_ref_db.override  = 0;
    SS_ref_db.symmetry  = 1;
    SS_ref_db.n_cat     = 2;
    SS_ref_db.n_xeos    = 2;
    SS_ref_db.n_em      = 2;
    SS_ref_db.n_sf      = 1;
    SS_ref_db.n_w       = 1;

     return SS_ref_db;
}

/**
    allocate memory for sb21_sp
*/
SS_ref G_SS_sb21_sp_init_function(SS_ref SS_ref_db,  global_variable gv){

    SS_ref_db.is_liq    = 0;
    SS_ref_db.override  = 0;
    SS_ref_db.symmetry  = 1;
    SS_ref_db.n_cat     = 6;
    SS_ref_db.n_xeos    = 2;
    SS_ref_db.n_em      = 2;
    SS_ref_db.n_sf      = 2;
    SS_ref_db.n_w       = 1;

     return SS_ref_db;
}

/**
    allocate memory for sb21_ol
*/
SS_ref G_SS_sb21_ol_init_function(SS_ref SS_ref_db,  global_variable gv){

    SS_ref_db.is_liq    = 0;
    SS_ref_db.override  = 0;
    SS_ref_db.symmetry  = 1;
    SS_ref_db.n_cat     = 2;
    SS_ref_db.n_xeos    = 2;
    SS_ref_db.n_em      = 2;
    SS_ref_db.n_sf      = 1;
    SS_ref_db.n_w       = 1;

     return SS_ref_db;
}

/**
    allocate memory for sb21_wa
*/
SS_ref G_SS_sb21_wa_init_function(SS_ref SS_ref_db,  global_variable gv){

    SS_ref_db.is_liq    = 0;
    SS_ref_db.override  = 0;
    SS_ref_db.symmetry  = 1;
    SS_ref_db.n_cat     = 2;
    SS_ref_db.n_xeos    = 2;
    SS_ref_db.n_em      = 2;
    SS_ref_db.n_sf      = 1;
    SS_ref_db.n_w       = 1;

     return SS_ref_db;
}

/**
    allocate memory for sb21_ri
*/
SS_ref G_SS_sb21_ri_init_function(SS_ref SS_ref_db,  global_variable gv){

    SS_ref_db.is_liq    = 0;
    SS_ref_db.override  = 0;
    SS_ref_db.symmetry  = 1;
    SS_ref_db.n_cat     = 2;
    SS_ref_db.n_xeos    = 2;
    SS_ref_db.n_em      = 2;
    SS_ref_db.n_sf      = 1;
    SS_ref_db.n_w       = 1;

     return SS_ref_db;
}

/**
    allocate memory for sb21_opx
*/
SS_ref G_SS_sb21_opx_init_function(SS_ref SS_ref_db,  global_variable gv){

    SS_ref_db.is_liq    = 0;
    SS_ref_db.override  = 0;
    SS_ref_db.symmetry  = 1;
    SS_ref_db.n_cat     = 6;
    SS_ref_db.n_xeos    = 4;
    SS_ref_db.n_em      = 4;
    SS_ref_db.n_sf      = 2;
    SS_ref_db.n_w       = 6;

     return SS_ref_db;
}

/**
    allocate memory for sb21_cpx
*/
SS_ref G_SS_sb21_cpx_init_function(SS_ref SS_ref_db,  global_variable gv){

    SS_ref_db.is_liq    = 0;
    SS_ref_db.override  = 0;
    SS_ref_db.symmetry  = 0;
    SS_ref_db.n_cat     = 8;
    SS_ref_db.n_xeos    = 5;
    SS_ref_db.n_em      = 5;
    SS_ref_db.n_sf      = 3;
    SS_ref_db.n_w       = 10;
    SS_ref_db.n_v       = 5;

     return SS_ref_db;
}

/**
    allocate memory for sb21_hpcpx
*/
SS_ref G_SS_sb21_hpcpx_init_function(SS_ref SS_ref_db,  global_variable gv){

    SS_ref_db.is_liq    = 0;
    SS_ref_db.override  = 0;
    SS_ref_db.symmetry  = 1;
    SS_ref_db.n_cat     = 2;
    SS_ref_db.n_xeos    = 2;
    SS_ref_db.n_em      = 2;
    SS_ref_db.n_sf      = 1;
    SS_ref_db.n_w       = 1;

     return SS_ref_db;
}

/**
    allocate memory for sb21_ak
*/
SS_ref G_SS_sb21_ak_init_function(SS_ref SS_ref_db,  global_variable gv){

    SS_ref_db.is_liq    = 0;
    SS_ref_db.override  = 0;
    SS_ref_db.symmetry  = 1;
    SS_ref_db.n_cat     = 5;
    SS_ref_db.n_xeos    = 3;
    SS_ref_db.n_em      = 3;
    SS_ref_db.n_sf      = 2;
    SS_ref_db.n_w       = 3;

     return SS_ref_db;
}

/**
    allocate memory for sb21_gtmj
*/
SS_ref G_SS_sb21_gtmj_init_function(SS_ref SS_ref_db,  global_variable gv){

    SS_ref_db.is_liq    = 0;
    SS_ref_db.override  = 0;
    SS_ref_db.symmetry  = 1;
    SS_ref_db.n_cat     = 9;
    SS_ref_db.n_xeos    = 5;
    SS_ref_db.n_em      = 5;
    SS_ref_db.n_sf      = 3;
    SS_ref_db.n_w       = 10;

     return SS_ref_db;
}

/**
    allocate memory for sb21_pv
*/
SS_ref G_SS_sb21_pv_init_function(SS_ref SS_ref_db,  global_variable gv){

    SS_ref_db.is_liq    = 0;
    SS_ref_db.override  = 0;
    SS_ref_db.symmetry  = 1;
    SS_ref_db.n_cat     = 5;
    SS_ref_db.n_xeos    = 3;
    SS_ref_db.n_em      = 3;
    SS_ref_db.n_sf      = 2;
    SS_ref_db.n_w       = 3;

     return SS_ref_db;
}

/**
    allocate memory for sb21_ppv
*/
SS_ref G_SS_sb21_ppv_init_function(SS_ref SS_ref_db,  global_variable gv){

    SS_ref_db.is_liq    = 0;
    SS_ref_db.override  = 0;
    SS_ref_db.symmetry  = 1;
    SS_ref_db.n_cat     = 5;
    SS_ref_db.n_xeos    = 3;
    SS_ref_db.n_em      = 3;
    SS_ref_db.n_sf      = 2;
    SS_ref_db.n_w       = 3;

     return SS_ref_db;
}

/**
    allocate memory for sb21_cf
*/
SS_ref G_SS_sb21_cf_init_function(SS_ref SS_ref_db,  global_variable gv){

    SS_ref_db.is_liq    = 0;
    SS_ref_db.override  = 0;
    SS_ref_db.symmetry  = 0;
    SS_ref_db.n_cat     = 5;
    SS_ref_db.n_xeos    = 3;
    SS_ref_db.n_em      = 3;
    SS_ref_db.n_sf      = 2;
    SS_ref_db.n_w       = 3;
    SS_ref_db.n_v       = 3;

     return SS_ref_db;
}

/**
    allocate memory for sb21_mw
*/
SS_ref G_SS_sb21_mw_init_function(SS_ref SS_ref_db,  global_variable gv){

    SS_ref_db.is_liq    = 0;
    SS_ref_db.override  = 0;
    SS_ref_db.symmetry  = 1;
    SS_ref_db.n_cat     = 6;
    SS_ref_db.n_xeos    = 3;
    SS_ref_db.n_em      = 3;
    SS_ref_db.n_sf      = 2;
    SS_ref_db.n_w       = 3;

     return SS_ref_db;
}

/**
    allocate memory for sb21_nal
*/
SS_ref G_SS_sb21_nal_init_function(SS_ref SS_ref_db,  global_variable gv){

    SS_ref_db.is_liq    = 0;
    SS_ref_db.override  = 0;
    SS_ref_db.symmetry  = 1;
    SS_ref_db.n_cat     = 6;
    SS_ref_db.n_xeos    = 3;
    SS_ref_db.n_em      = 3;
    SS_ref_db.n_sf      = 3;
    SS_ref_db.n_w       = 3;

     return SS_ref_db;
}

void SB_SS_init_sb21(          SS_init_type                *SS_init,
                            global_variable      gv                             ){

    for (int iss = 0; iss < gv.len_ss; iss++){
        if      (strcmp( gv.SS_list[iss], "plg")  == 0 ){
            SS_init[iss]  = G_SS_sb21_plg_init_function;                }
        else if (strcmp( gv.SS_list[iss], "sp")  == 0 ){
            SS_init[iss]  = G_SS_sb21_sp_init_function;                 }
        else if (strcmp( gv.SS_list[iss], "ol")  == 0 ){
            SS_init[iss]  = G_SS_sb21_ol_init_function;                 }
        else if (strcmp( gv.SS_list[iss], "wa")  == 0 ){
            SS_init[iss]  = G_SS_sb21_wa_init_function;                 }
        else if (strcmp( gv.SS_list[iss], "ri")  == 0 ){
            SS_init[iss]  = G_SS_sb21_ri_init_function;                 }
        else if (strcmp( gv.SS_list[iss], "opx")  == 0 ){
            SS_init[iss]  = G_SS_sb21_opx_init_function;                }
        else if (strcmp( gv.SS_list[iss], "cpx")  == 0 ){
            SS_init[iss]  = G_SS_sb21_cpx_init_function;                }
        else if (strcmp( gv.SS_list[iss], "hpcpx")  == 0 ){
            SS_init[iss]  = G_SS_sb21_hpcpx_init_function;              }
        else if (strcmp( gv.SS_list[iss], "ak")  == 0 ){
            SS_init[iss]  = G_SS_sb21_ak_init_function;                 }
        else if (strcmp( gv.SS_list[iss], "gtmj")  == 0 ){
            SS_init[iss]  = G_SS_sb21_gtmj_init_function;               }
        else if (strcmp( gv.SS_list[iss], "pv")  == 0 ){
            SS_init[iss]  = G_SS_sb21_pv_init_function;                 }
        else if (strcmp( gv.SS_list[iss], "ppv")  == 0 ){
            SS_init[iss]  = G_SS_sb21_ppv_init_function;                }
        else if (strcmp( gv.SS_list[iss], "cf")  == 0 ){
            SS_init[iss]  = G_SS_sb21_cf_init_function;                 }
        else if (strcmp( gv.SS_list[iss], "mw")  == 0 ){
            SS_init[iss]  = G_SS_sb21_mw_init_function;                 }
        else if (strcmp( gv.SS_list[iss], "nal")  == 0 ){
            SS_init[iss]  = G_SS_sb21_nal_init_function;                }
        else{
            printf("\nsolid solution '%s' is not in the database, cannot be initiated\n", gv.SS_list[iss]);
        }
    }
}

void SB_SS_init(	        	    SS_init_type 		*SS_init,
									global_variable 	 gv				){

	if (gv.EM_database == 0){				// database //
		SB_SS_init_sb11(	 				SS_init,
											gv							);
	}
	if (gv.EM_database == 1){				// database //
		SB_SS_init_sb21(	 				SS_init,
											gv							);
	}
}