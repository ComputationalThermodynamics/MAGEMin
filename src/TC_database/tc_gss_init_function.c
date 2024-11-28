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
#include "tc_gss_init_function.h"


/**************************************************************************************/
/* Import text file extracted from Miron et al. (2017) database:                      */
/* aq17-thermofun.json from db.thermohub.org, v. 17.07.2021 16:49:29                  */
/*--------------------------------------------------------------------------          */
/* Miron, G. D., Wagner, T., Kulik, D. A., & Lothenbach, B. (2017). An                */
/* internally consistent thermodynamic dataset for aqueous species in the             */
/* system Ca-Mg-Na-K-Al-Si-OHC-Cl to 800 C and 5 kbar. American Journal               */
/* of Science, 317(7), 755-806.                                                       */
/* DOI: https://doi.org/10.2475/07.2017.01                                            */
/**************************************************************************************/


SS_ref G_SS_aq17_init_function(SS_ref SS_ref_db,  global_variable gv){
    
    SS_ref_db.n_cat     = 0;
    SS_ref_db.is_liq    = 0;
    SS_ref_db.symmetry  = 0;
    SS_ref_db.n_sf      = 25;
    SS_ref_db.n_em      = 25;//44;
    // SS_ref_db.n_v       = 15;
    // SS_ref_db.n_w       = 105;
    SS_ref_db.n_xeos    = 25;//44;
    
    return SS_ref_db;
}

/**************************************************************************************/
/**************************************************************************************/
/**********************METABASITE DATABASE (Gree et al., 2016)*************************/
/**************************************************************************************/
/**************************************************************************************/

/**
    allocate memory for L
*/
SS_ref G_SS_mb_liq_init_function(SS_ref SS_ref_db,  global_variable gv){
    
    SS_ref_db.n_cat     = 0;
    SS_ref_db.is_liq    = 1;
    SS_ref_db.symmetry  = 1;
    SS_ref_db.n_sf      = 11;
    SS_ref_db.n_em      = 9;
    SS_ref_db.n_w       = 36;
    SS_ref_db.n_xeos    = 8;
    
    return SS_ref_db;
}

/**
    allocate memory for hb
*/
SS_ref G_SS_mb_hb_init_function(SS_ref SS_ref_db,  global_variable gv){
    
    SS_ref_db.n_cat     = 0;
    SS_ref_db.is_liq    = 0;
    SS_ref_db.symmetry  = 0;
    SS_ref_db.n_sf      = 18;
    SS_ref_db.n_em      = 11;
    SS_ref_db.n_v       = 11;
    SS_ref_db.n_w       = 55;
    SS_ref_db.n_xeos    = 10;
    
    return SS_ref_db;
}

/**
    allocate memory for aug
*/
SS_ref G_SS_mb_aug_init_function(SS_ref SS_ref_db,  global_variable gv){
    
    SS_ref_db.n_cat     = 0;
    SS_ref_db.is_liq    = 0;
    SS_ref_db.symmetry  = 0;
    SS_ref_db.n_sf      = 12;
    SS_ref_db.n_em      = 8;
    SS_ref_db.n_v       = 8;
    SS_ref_db.n_w       = 28;
    SS_ref_db.n_xeos    = 7;
    
    return SS_ref_db;
}

/**
    allocate memory for dio
*/
SS_ref G_SS_mb_dio_init_function(SS_ref SS_ref_db,  global_variable gv){
    
    SS_ref_db.n_cat     = 0;
    SS_ref_db.is_liq    = 0;
    SS_ref_db.symmetry  = 1;
    SS_ref_db.n_sf      = 12;
    SS_ref_db.n_em      = 7;
    SS_ref_db.n_w       = 21;
    SS_ref_db.n_xeos    = 6;
    
    return SS_ref_db;
}

/**
    allocate memory for opx
*/
SS_ref G_SS_mb_opx_init_function(SS_ref SS_ref_db,  global_variable gv){
    
    SS_ref_db.n_cat     = 0;
    SS_ref_db.is_liq    = 0;
    SS_ref_db.symmetry  = 0;
    SS_ref_db.n_sf      = 9;
    SS_ref_db.n_em      = 6;
    SS_ref_db.n_v       = 6;
    SS_ref_db.n_w       = 15;
    SS_ref_db.n_xeos    = 5;
    
    return SS_ref_db;
}

/**
    allocate memory for g
*/
SS_ref G_SS_mb_g_init_function(SS_ref SS_ref_db,  global_variable gv){
    
    SS_ref_db.n_cat     = 0;
    SS_ref_db.is_liq    = 0;
    SS_ref_db.symmetry  = 0;
    SS_ref_db.n_sf      = 5;
    SS_ref_db.n_em      = 4;
    SS_ref_db.n_v       = 4;
    SS_ref_db.n_w       = 6;
    SS_ref_db.n_xeos    = 3;
    
    return SS_ref_db;
}

/**
    allocate memory for ol
*/
SS_ref G_SS_mb_ol_init_function(SS_ref SS_ref_db,  global_variable gv){
    
    SS_ref_db.n_cat     = 0;
    SS_ref_db.is_liq    = 0;
    SS_ref_db.symmetry  = 1;
    SS_ref_db.n_sf      = 2;
    SS_ref_db.n_em      = 2;
    SS_ref_db.n_w       = 1;
    SS_ref_db.n_xeos    = 1;
    
    return SS_ref_db;
}

/**
    allocate memory for fsp
*/
SS_ref G_SS_mb_fsp_init_function(SS_ref SS_ref_db,  global_variable gv){
    
    SS_ref_db.n_cat     = 0;
    SS_ref_db.is_liq    = 0;
    SS_ref_db.symmetry  = 0;
    SS_ref_db.n_sf      = 5;
    SS_ref_db.n_em      = 3;
    SS_ref_db.n_v       = 3;
    SS_ref_db.n_w       = 3;
    SS_ref_db.n_xeos    = 2;
    
    return SS_ref_db;
}

/**
    allocate memory for abc
*/
SS_ref G_SS_mb_abc_init_function(SS_ref SS_ref_db,  global_variable gv){
    
    SS_ref_db.n_cat     = 0;
    SS_ref_db.is_liq    = 0;
    SS_ref_db.symmetry  = 0;
    SS_ref_db.n_sf      = 2;
    SS_ref_db.n_em      = 2;
    SS_ref_db.n_v       = 2;
    SS_ref_db.n_w       = 1;
    SS_ref_db.n_xeos    = 1;
    
    return SS_ref_db;
}

/**
    allocate memory for k4tr
*/
SS_ref G_SS_mb_k4tr_init_function(SS_ref SS_ref_db,  global_variable gv){
    
    SS_ref_db.n_cat     = 0;
    SS_ref_db.is_liq    = 0;
    SS_ref_db.symmetry  = 0;
    SS_ref_db.n_sf      = 5;
    SS_ref_db.n_em      = 3;
    SS_ref_db.n_v       = 3;
    SS_ref_db.n_w       = 3;
    SS_ref_db.n_xeos    = 2;
    
    return SS_ref_db;
}
/**
    allocate memory for spn
*/
SS_ref G_SS_mb_spn_init_function(SS_ref SS_ref_db, global_variable gv){
    
    SS_ref_db.n_cat     = 0;
    SS_ref_db.is_liq    = 0;
    SS_ref_db.symmetry  = 1;
    SS_ref_db.n_sf      = 4;
    SS_ref_db.n_em      = 3;
    SS_ref_db.n_w       = 3;
    SS_ref_db.n_xeos    = 2;
    
    return SS_ref_db;
}
/**
    allocate memory for sp
*/
SS_ref G_SS_mb_sp_init_function(SS_ref SS_ref_db,  global_variable gv){
    
    SS_ref_db.n_cat     = 0;
    SS_ref_db.is_liq    = 0;
    SS_ref_db.symmetry  = 1;
    SS_ref_db.n_sf      = 5;
    SS_ref_db.n_em      = 4;
    SS_ref_db.n_w       = 6;
    SS_ref_db.n_xeos    = 3;
    
    return SS_ref_db;
}

/**
    allocate memory for ilm
*/
SS_ref G_SS_mb_ilm_init_function(SS_ref SS_ref_db,  global_variable gv){
    
    SS_ref_db.n_cat     = 0;
    SS_ref_db.is_liq    = 0;
    SS_ref_db.symmetry  = 1;
    SS_ref_db.n_sf      = 6;
    SS_ref_db.n_em      = 3;
    SS_ref_db.n_w       = 3;
    SS_ref_db.n_xeos    = 2;
    
    return SS_ref_db;
}

/**
    allocate memory for ilmm
*/
SS_ref G_SS_mb_ilmm_init_function(SS_ref SS_ref_db,  global_variable gv){
    
    SS_ref_db.n_cat     = 0;
    SS_ref_db.is_liq    = 0;
    SS_ref_db.symmetry  = 1;
    SS_ref_db.n_sf      = 7;
    SS_ref_db.n_em      = 4;
    SS_ref_db.n_w       = 6;
    SS_ref_db.n_xeos    = 3;
    
    return SS_ref_db;
}

/**
    allocate memory for ep
*/
SS_ref G_SS_mb_ep_init_function(SS_ref SS_ref_db,  global_variable gv){
    
    SS_ref_db.n_cat     = 0;
    SS_ref_db.is_liq    = 0;
    SS_ref_db.symmetry  = 1;
    SS_ref_db.n_sf      = 4;
    SS_ref_db.n_em      = 3;
    SS_ref_db.n_w       = 3;
    SS_ref_db.n_xeos    = 2;
    
    return SS_ref_db;
}

/**
    allocate memory for bi
*/
SS_ref G_SS_mb_bi_init_function(SS_ref SS_ref_db,  global_variable gv){
    
    SS_ref_db.n_cat     = 0;
    SS_ref_db.is_liq    = 0;
    SS_ref_db.symmetry  = 1;
    SS_ref_db.n_sf      = 11;
    SS_ref_db.n_em      = 6;
    SS_ref_db.n_w       = 15;
    SS_ref_db.n_xeos    = 5;
    
    return SS_ref_db;
}

/**
    allocate memory for mu
*/
SS_ref G_SS_mb_mu_init_function(SS_ref SS_ref_db,  global_variable gv){
    
    SS_ref_db.n_cat     = 0;
    SS_ref_db.is_liq    = 0;
    SS_ref_db.symmetry  = 0;
    SS_ref_db.n_sf      = 10;
    SS_ref_db.n_em      = 6;
    SS_ref_db.n_v       = 6;
    SS_ref_db.n_w       = 15;
    SS_ref_db.n_xeos    = 5;
    
    return SS_ref_db;
}

/**
    allocate memory for chl
*/
SS_ref G_SS_mb_chl_init_function(SS_ref SS_ref_db,  global_variable gv){
    
    SS_ref_db.n_cat     = 0;
    SS_ref_db.is_liq    = 0;
    SS_ref_db.symmetry  = 1;
    SS_ref_db.n_sf      = 11;
    SS_ref_db.n_em      = 7;
    SS_ref_db.n_w       = 21;
    SS_ref_db.n_xeos    = 6;
    
    return SS_ref_db;
}


/**************************************************************************************/
/**************************************************************************************/
/*********************METAPELITE DATABASE (White et al., 2014)*************************/
/**************************************************************************************/
/**************************************************************************************/

/**
    allocate memory for liq_mp
*/

SS_ref G_SS_mp_liq_init_function(SS_ref SS_ref_db,  global_variable gv){
    
    SS_ref_db.n_cat     = 0;
    SS_ref_db.is_liq    = 1;
    SS_ref_db.symmetry  = 1;
    SS_ref_db.n_sf      = 10;
    SS_ref_db.n_em      = 8;
    SS_ref_db.n_w       = 28;
    SS_ref_db.n_xeos    = 7;
    
    return SS_ref_db;
}

/**
    allocate memory for fsp_mp
*/
SS_ref G_SS_mp_fsp_init_function(SS_ref SS_ref_db,  global_variable gv){
    
    SS_ref_db.n_cat     = 0;
    SS_ref_db.is_liq    = 0;
    SS_ref_db.symmetry  = 0;
    SS_ref_db.n_sf      = 5;
    SS_ref_db.n_em      = 3;
	SS_ref_db.n_v       = 3;
    SS_ref_db.n_w       = 3;
    SS_ref_db.n_xeos    = 2;
    
    return SS_ref_db;
}

/**
    allocate memory for bi_mp
*/
SS_ref G_SS_mp_bi_init_function(SS_ref SS_ref_db,  global_variable gv){
    
    SS_ref_db.n_cat     = 0;
    SS_ref_db.is_liq    = 0;
    SS_ref_db.symmetry  = 1;
    SS_ref_db.n_sf      = 13;
    SS_ref_db.n_em      = 7;
    SS_ref_db.n_w       = 21;
    SS_ref_db.n_xeos    = 6;
    
    return SS_ref_db;
}

/**
    allocate memory for g_mp
*/
SS_ref G_SS_mp_g_init_function(SS_ref SS_ref_db,  global_variable gv){
    
    SS_ref_db.n_cat     = 0;
    SS_ref_db.is_liq    = 0;
    SS_ref_db.symmetry  = 0;
    SS_ref_db.n_sf      = 6;
    SS_ref_db.n_em      = 5;
    SS_ref_db.n_v       = 5;
    SS_ref_db.n_w       = 10;
    SS_ref_db.n_xeos    = 4;
    
    return SS_ref_db;
}

/**
    allocate memory for ep_mp
*/
SS_ref G_SS_mp_ep_init_function(SS_ref SS_ref_db,  global_variable gv){
    
    SS_ref_db.n_cat     = 0;
    SS_ref_db.is_liq    = 0;
    SS_ref_db.symmetry  = 1;
    SS_ref_db.n_sf      = 4;
    SS_ref_db.n_em      = 3;
    SS_ref_db.n_w       = 3;
    SS_ref_db.n_xeos    = 2;
    
    return SS_ref_db;
}

/**
    allocate memory for ma_mp
*/
SS_ref G_SS_mp_ma_init_function(SS_ref SS_ref_db,  global_variable gv){
    
    SS_ref_db.n_cat     = 0;
    SS_ref_db.is_liq    = 0;
    SS_ref_db.symmetry  = 0;
    SS_ref_db.n_sf      = 10;
    SS_ref_db.n_em      = 6;
    SS_ref_db.n_v       = 6;
    SS_ref_db.n_w       = 15;
    SS_ref_db.n_xeos    = 5;
    
    return SS_ref_db;
}

/**
    allocate memory for mu_mp
*/
SS_ref G_SS_mp_mu_init_function(SS_ref SS_ref_db,  global_variable gv){
    
    SS_ref_db.n_cat     = 0;
    SS_ref_db.is_liq    = 0;
    SS_ref_db.symmetry  = 0;
    SS_ref_db.n_sf      = 10;
    SS_ref_db.n_em      = 6;
    SS_ref_db.n_v       = 6;
    SS_ref_db.n_w       = 15;
    SS_ref_db.n_xeos    = 5;
    
    return SS_ref_db;
}

/**
    allocate memory for opx_mp
*/
SS_ref G_SS_mp_opx_init_function(SS_ref SS_ref_db,  global_variable gv){
    
    SS_ref_db.n_cat     = 0;
    SS_ref_db.is_liq    = 0;
    SS_ref_db.symmetry  = 0;
    SS_ref_db.n_sf      = 11;
    SS_ref_db.n_em      = 7;
    SS_ref_db.n_v       = 7;
    SS_ref_db.n_w       = 21;
    SS_ref_db.n_xeos    = 6;
    
    return SS_ref_db;
}

/**
    allocate memory for sa_mp
*/
SS_ref G_SS_mp_sa_init_function(SS_ref SS_ref_db,  global_variable gv){
    
    SS_ref_db.n_cat     = 0;
    SS_ref_db.is_liq    = 0;
    SS_ref_db.symmetry  = 1;
    SS_ref_db.n_sf      = 8;
    SS_ref_db.n_em      = 5;
    SS_ref_db.n_w       = 10;
    SS_ref_db.n_xeos    = 4;
    
    return SS_ref_db;
}

/**
    allocate memory for cd_mp
*/
SS_ref G_SS_mp_cd_init_function(SS_ref SS_ref_db,  global_variable gv){
    
    SS_ref_db.n_cat     = 0;
    SS_ref_db.is_liq    = 0;
    SS_ref_db.symmetry  = 1;
    SS_ref_db.n_sf      = 5;
    SS_ref_db.n_em      = 4;
    SS_ref_db.n_w       = 6;
    SS_ref_db.n_xeos    = 3;
    
    return SS_ref_db;
}

/**
    allocate memory for st_mp
*/
SS_ref G_SS_mp_st_init_function(SS_ref SS_ref_db,  global_variable gv){
    
    SS_ref_db.n_cat     = 0;
    SS_ref_db.is_liq    = 0;
    SS_ref_db.symmetry  = 1;
    SS_ref_db.n_sf      = 7;
    SS_ref_db.n_em      = 5;
    SS_ref_db.n_w       = 10;
    SS_ref_db.n_xeos    = 4;
    
    return SS_ref_db;
}

/**
    allocate memory for chl_mp
*/
SS_ref G_SS_mp_chl_init_function(SS_ref SS_ref_db,  global_variable gv){
    
    SS_ref_db.n_cat     = 0;
    SS_ref_db.is_liq    = 0;
    SS_ref_db.symmetry  = 1;
    SS_ref_db.n_sf      = 12;
    SS_ref_db.n_em      = 8;
    SS_ref_db.n_w       = 28;
    SS_ref_db.n_xeos    = 7;
    
    return SS_ref_db;
}

/**
    allocate memory for ctd_mp
*/
SS_ref G_SS_mp_ctd_init_function(SS_ref SS_ref_db,  global_variable gv){
    
    SS_ref_db.n_cat     = 0;
    SS_ref_db.is_liq    = 0;
    SS_ref_db.symmetry  = 1;
    SS_ref_db.n_sf      = 5;
    SS_ref_db.n_em      = 4;
    SS_ref_db.n_w       = 6;
    SS_ref_db.n_xeos    = 3;
    
    return SS_ref_db;
}

/**
    allocate memory for sp_mp
*/
SS_ref G_SS_mp_sp_init_function(SS_ref SS_ref_db,  global_variable gv){
    
    SS_ref_db.n_cat     = 0;
    SS_ref_db.is_liq    = 0;
    SS_ref_db.symmetry  = 1;
    SS_ref_db.n_sf      = 5;
    SS_ref_db.n_em      = 4;
    SS_ref_db.n_w       = 6;
    SS_ref_db.n_xeos    = 3;
    
    return SS_ref_db;
}

/**
    allocate memory for ilm
*/
SS_ref G_SS_mp_ilm_init_function(SS_ref SS_ref_db,  global_variable gv){
    
    SS_ref_db.n_cat     = 0;
    SS_ref_db.is_liq    = 0;
    SS_ref_db.symmetry  = 1;
    SS_ref_db.n_sf      = 6;
    SS_ref_db.n_em      = 3;
    SS_ref_db.n_w       = 3;
    SS_ref_db.n_xeos    = 2;
    
    return SS_ref_db;
}
/**
    allocate memory for ilmm_mp
*/
SS_ref G_SS_mp_ilmm_init_function(SS_ref SS_ref_db,  global_variable gv){
    
    SS_ref_db.n_cat     = 0;
    SS_ref_db.is_liq    = 0;
    SS_ref_db.symmetry  = 1;
    SS_ref_db.n_sf      = 7;
    SS_ref_db.n_em      = 5;
    SS_ref_db.n_w       = 10;
    SS_ref_db.n_xeos    = 4;
    
    return SS_ref_db;
}

/**
    allocate memory for mt_mp
*/
SS_ref G_SS_mp_mt_init_function(SS_ref SS_ref_db,  global_variable gv){
    
    SS_ref_db.n_cat     = 0;
    SS_ref_db.is_liq    = 0;
    SS_ref_db.symmetry  = 1;
    SS_ref_db.n_sf      = 5;
    SS_ref_db.n_em      = 3;
    SS_ref_db.n_w       = 3;
    SS_ref_db.n_xeos    = 2;
    
    return SS_ref_db;
}



/**************************************************************************************/
/**************************************************************************************/
/*********************IGNEOUS DATABASE (Holland et al., 2018)**************************/
/**************************************************************************************/
/**************************************************************************************/
/**
    allocate memory for fper
*/
SS_ref G_SS_ig_fper_init_function(SS_ref SS_ref_db,  global_variable gv){
    
    SS_ref_db.n_cat     = 0;
    SS_ref_db.is_liq    = 0;
    SS_ref_db.symmetry  = 1;
    SS_ref_db.n_sf      = 2;
    SS_ref_db.n_em      = 2;
    SS_ref_db.n_w       = 1;
    SS_ref_db.n_xeos    = 1;
    
    return SS_ref_db;
}

/** 
  allocate memory for biotite
*/
SS_ref G_SS_ig_bi_init_function(SS_ref SS_ref_db,  global_variable gv){		

	SS_ref_db.n_cat     = 0;
    SS_ref_db.is_liq    = 0;	
	SS_ref_db.symmetry  = 1;	
	SS_ref_db.n_sf      = 11;
	SS_ref_db.n_em      = 6;
	SS_ref_db.n_w       = 15;
	SS_ref_db.n_xeos    = 5;

	return SS_ref_db;
}

/** 
  allocate memory for clinopyroxene
*/
SS_ref G_SS_ig_cpx_init_function(SS_ref SS_ref_db,  global_variable gv){

	SS_ref_db.n_cat     = 0;
    SS_ref_db.is_liq    = 0;	
	SS_ref_db.symmetry  = 0;						  					  
	SS_ref_db.n_sf      = 13;
	SS_ref_db.n_em      = 10;
	SS_ref_db.n_v       = 10;
	SS_ref_db.n_w       = 45;
	SS_ref_db.n_xeos    = 9;	

	return SS_ref_db;
}

/** 
  allocate memory for cordierite
*/
SS_ref G_SS_ig_cd_init_function(SS_ref SS_ref_db,  global_variable gv){

	SS_ref_db.n_cat     = 0;
    SS_ref_db.is_liq    = 0;		
	SS_ref_db.symmetry  = 1;				  					  
	SS_ref_db.n_sf      = 4;
	SS_ref_db.n_em      = 3;
	SS_ref_db.n_w       = 3;
	SS_ref_db.n_xeos    = 2;

	return SS_ref_db;	
}

/** 
  allocate memory for epidote
*/
SS_ref G_SS_ig_ep_init_function(SS_ref SS_ref_db,  global_variable gv){

	SS_ref_db.n_cat     = 0;
    SS_ref_db.is_liq     = 0;	
	SS_ref_db.symmetry   = 1;						  					  
	SS_ref_db.n_sf       = 4;
	SS_ref_db.n_em       = 3;
	SS_ref_db.n_w        = 3;
	SS_ref_db.n_xeos     = 2;
  
	return SS_ref_db;	
}

/** 
  allocate memory for fluid
*/
SS_ref G_SS_ig_fl_init_function(SS_ref SS_ref_db,  global_variable gv){
    
    SS_ref_db.n_cat     = 0;
    SS_ref_db.is_liq    = 0;
    SS_ref_db.symmetry  = 1;
    SS_ref_db.n_sf      = 12;
    SS_ref_db.n_em      = 11;
    SS_ref_db.n_w       = 55;
    SS_ref_db.n_xeos    = 10;
    
    return SS_ref_db;
}


/** 
  allocate memory for garnet
*/
SS_ref G_SS_ig_g_init_function(SS_ref SS_ref_db,  global_variable gv){

	SS_ref_db.n_cat     = 0;
    SS_ref_db.is_liq      = 0;	
	SS_ref_db.symmetry    = 0;						  					  
	SS_ref_db.n_sf        = 8;
	SS_ref_db.n_em        = 6;
	SS_ref_db.n_v         = 6;
	SS_ref_db.n_w         = 15;
	SS_ref_db.n_xeos      = 5;

	return SS_ref_db;
}

/** 
  allocate memory for hornblende
*/
SS_ref G_SS_ig_hb_init_function(SS_ref SS_ref_db,  global_variable gv){

	SS_ref_db.n_cat     = 0;
    SS_ref_db.is_liq     = 0;	
	SS_ref_db.symmetry   = 0;					  					  
	SS_ref_db.n_sf       = 18;
	SS_ref_db.n_em       = 11;
	SS_ref_db.n_v        = 11;
	SS_ref_db.n_w        = 55;
	SS_ref_db.n_xeos     = 10;
		  
	return SS_ref_db;	
}


/**
    allocate memory for ilm
*/
SS_ref G_SS_ig_ilm_init_function(SS_ref SS_ref_db,  global_variable gv){
    
	SS_ref_db.n_cat     = 0;
    SS_ref_db.is_liq    = 0;	
	SS_ref_db.symmetry  = 1;							  					  
	SS_ref_db.n_sf      = 6;
	SS_ref_db.n_em      = 3;
	SS_ref_db.n_w       = 3;
	SS_ref_db.n_xeos    = 2; 
    
    return SS_ref_db;
}



/**
    allocate memory for liqHw
*/
SS_ref G_SS_ig_liq_init_function(SS_ref SS_ref_db,  global_variable gv){
    
    SS_ref_db.n_cat     = 0;
    SS_ref_db.is_liq    = 1;
    SS_ref_db.symmetry  = 0;
    SS_ref_db.n_sf      = 18;
    SS_ref_db.n_em      = 12;
    SS_ref_db.n_v       = 12;
    SS_ref_db.n_w       = 66;
    SS_ref_db.n_xeos    = 11;
    
    return SS_ref_db;
}

/** 
  allocate memory for muscovite
*/
SS_ref G_SS_ig_mu_init_function(SS_ref SS_ref_db,  global_variable gv){		

	SS_ref_db.n_cat     = 0;
    SS_ref_db.is_liq     = 0;	
	SS_ref_db.symmetry   = 0;						  					  
	SS_ref_db.n_sf       = 10;
	SS_ref_db.n_em       = 6;
	SS_ref_db.n_v        = 6;
	SS_ref_db.n_w        = 15;
	SS_ref_db.n_xeos     = 5;

	return SS_ref_db;
}

/** 
  allocate memory for olivine
*/
SS_ref G_SS_ig_ol_init_function(SS_ref SS_ref_db,  global_variable gv){		

	SS_ref_db.n_cat     = 0;
    SS_ref_db.is_liq     = 0;	
	SS_ref_db.symmetry   = 1;						  					  
	SS_ref_db.n_sf       = 5;
	SS_ref_db.n_em       = 4;
	SS_ref_db.n_w        = 6;
	SS_ref_db.n_xeos     = 3;

	return SS_ref_db;
}

/** 
  allocate memory for orthopyroxene
*/
SS_ref G_SS_ig_opx_init_function(SS_ref SS_ref_db,  global_variable gv){		

	SS_ref_db.n_cat     = 0;
    SS_ref_db.is_liq    = 0;	
	SS_ref_db.symmetry  = 0;					  					  
	SS_ref_db.n_sf      = 12;
	SS_ref_db.n_em      = 9;
	SS_ref_db.n_v       = 9;
	SS_ref_db.n_w       = 36;
	SS_ref_db.n_xeos    = 8;

	return SS_ref_db;
}

/** 
  allocate memory for plagioclase
*/
SS_ref G_SS_ig_fsp_init_function(SS_ref SS_ref_db,  global_variable gv){

	SS_ref_db.n_cat     = 0;
    SS_ref_db.is_liq     = 0;
	SS_ref_db.symmetry   = 0;
	SS_ref_db.n_sf       = 5;
	SS_ref_db.n_em       = 3;
	SS_ref_db.n_v        = 3;
	SS_ref_db.n_w        = 3;
	SS_ref_db.n_xeos     = 2;

	return SS_ref_db;										  
}

/** 
  allocate memory for spn
*/
SS_ref G_SS_ig_spn_init_function(SS_ref SS_ref_db,  global_variable gv){		

    SS_ref_db.n_cat     = 0;
    SS_ref_db.is_liq    = 0;
    SS_ref_db.symmetry  = 0;
    SS_ref_db.n_sf      = 10;
    SS_ref_db.n_em      = 8;
    SS_ref_db.n_v       = 8;
    SS_ref_db.n_w       = 28;
    SS_ref_db.n_xeos    = 7;

	return SS_ref_db;
}


/**************************************************************************************/
/**************************************************************************************/
/******************IGNEOUS ALKALINE DRY DATABASE (Weller et al., 2021)*****************/
/**************************************************************************************/
/**************************************************************************************/

/**
    allocate memory for liq
*/
SS_ref G_SS_igad_liq_init_function(SS_ref SS_ref_db,  global_variable gv){
    
    SS_ref_db.n_cat     = 0;
    SS_ref_db.is_liq    = 0;
    SS_ref_db.symmetry  = 0;
    SS_ref_db.n_sf      = 18;
    SS_ref_db.n_em      = 14;
    SS_ref_db.n_v       = 14;
    SS_ref_db.n_w       = 91;
    SS_ref_db.n_xeos    = 13;
    
    return SS_ref_db;
}


/**
    allocate memory for fsp
*/
SS_ref G_SS_igad_fsp_init_function(SS_ref SS_ref_db,  global_variable gv){
    
    SS_ref_db.n_cat     = 0;
    SS_ref_db.is_liq    = 0;
    SS_ref_db.symmetry  = 0;
    SS_ref_db.n_sf      = 5;
    SS_ref_db.n_em      = 3;
    SS_ref_db.n_v       = 3;
    SS_ref_db.n_w       = 3;
    SS_ref_db.n_xeos    = 2;
    
    return SS_ref_db;
}

/**
    allocate memory for spn
*/
SS_ref G_SS_igad_spn_init_function(SS_ref SS_ref_db,  global_variable gv){
    
    SS_ref_db.n_cat     = 0;
    SS_ref_db.is_liq    = 0;
    SS_ref_db.symmetry  = 0;
    SS_ref_db.n_sf      = 10;
    SS_ref_db.n_em      = 8;
    SS_ref_db.n_v       = 8;
    SS_ref_db.n_w       = 28;
    SS_ref_db.n_xeos    = 7;
    
    return SS_ref_db;
}

/**
    allocate memory for g
*/
SS_ref G_SS_igad_g_init_function(SS_ref SS_ref_db,  global_variable gv){
    
    SS_ref_db.n_cat     = 0;
    SS_ref_db.is_liq    = 0;
    SS_ref_db.symmetry  = 0;
    SS_ref_db.n_sf      = 8;
    SS_ref_db.n_em      = 6;
    SS_ref_db.n_v       = 6;
    SS_ref_db.n_w       = 15;
    SS_ref_db.n_xeos    = 5;
    
    return SS_ref_db;
}

/**
    allocate memory for ol
*/
SS_ref G_SS_igad_ol_init_function(SS_ref SS_ref_db,  global_variable gv){
    
    SS_ref_db.n_cat     = 0;
    SS_ref_db.is_liq    = 0;
    SS_ref_db.symmetry  = 1;
    SS_ref_db.n_sf      = 5;
    SS_ref_db.n_em      = 4;
    SS_ref_db.n_w       = 6;
    SS_ref_db.n_xeos    = 3;
    
    return SS_ref_db;
}

/**
    allocate memory for opx
*/
SS_ref G_SS_igad_opx_init_function(SS_ref SS_ref_db,  global_variable gv){
    
    SS_ref_db.n_cat     = 0;
    SS_ref_db.is_liq    = 0;
    SS_ref_db.symmetry  = 0;
    SS_ref_db.n_sf      = 12;
    SS_ref_db.n_em      = 9;
    SS_ref_db.n_v       = 9;
    SS_ref_db.n_w       = 36;
    SS_ref_db.n_xeos    = 8;
    
    return SS_ref_db;
}

/**
    allocate memory for cpx
*/
SS_ref G_SS_igad_cpx_init_function(SS_ref SS_ref_db,  global_variable gv){
    
    SS_ref_db.n_cat     = 0;
    SS_ref_db.is_liq    = 0;
    SS_ref_db.symmetry  = 0;
    SS_ref_db.n_sf      = 13;
    SS_ref_db.n_em      = 10;
    SS_ref_db.n_v       = 10;
    SS_ref_db.n_w       = 45;
    SS_ref_db.n_xeos    = 9;
    
    return SS_ref_db;
}

/**
    allocate memory for ilm
*/
SS_ref G_SS_igad_ilm_init_function(SS_ref SS_ref_db,  global_variable gv){
    
    SS_ref_db.n_cat     = 0;
    SS_ref_db.is_liq    = 0;
    SS_ref_db.symmetry  = 1;
    SS_ref_db.n_sf      = 8;
    SS_ref_db.n_em      = 5;
    SS_ref_db.n_w       = 10;
    SS_ref_db.n_xeos    = 4;
    
    return SS_ref_db;
}

/**
    allocate memory for ness
*/
SS_ref G_SS_igad_ness_init_function(SS_ref SS_ref_db,  global_variable gv){
    
    SS_ref_db.n_cat     = 0;
    SS_ref_db.is_liq    = 0;
    SS_ref_db.symmetry  = 0;
    SS_ref_db.n_sf      = 9;
    SS_ref_db.n_em      = 6;
    SS_ref_db.n_v       = 6;
    SS_ref_db.n_w       = 15;
    SS_ref_db.n_xeos    = 5;
    
    return SS_ref_db;
}

/**
    allocate memory for lct
*/
SS_ref G_SS_igad_lct_init_function(SS_ref SS_ref_db,  global_variable gv){
    
    SS_ref_db.n_cat     = 0;
    SS_ref_db.is_liq    = 0;
    SS_ref_db.symmetry  = 0;
    SS_ref_db.n_sf      = 2;
    SS_ref_db.n_em      = 2;
    SS_ref_db.n_v       = 2;
    SS_ref_db.n_w       = 1;
    SS_ref_db.n_xeos    = 1;
    
    return SS_ref_db;
}

/**
    allocate memory for kals
*/
SS_ref G_SS_igad_kals_init_function(SS_ref SS_ref_db,  global_variable gv){
    
    SS_ref_db.n_cat     = 0;
    SS_ref_db.is_liq    = 0;
    SS_ref_db.symmetry  = 0;
    SS_ref_db.n_sf      = 2;
    SS_ref_db.n_em      = 2;
    SS_ref_db.n_v       = 2;
    SS_ref_db.n_w       = 1;
    SS_ref_db.n_xeos    = 1;
    
    return SS_ref_db;
}

/**
    allocate memory for mel
*/
SS_ref G_SS_igad_mel_init_function(SS_ref SS_ref_db,  global_variable gv){
    
    SS_ref_db.n_cat     = 0;
    SS_ref_db.is_liq    = 0;
    SS_ref_db.symmetry  = 1;
    SS_ref_db.n_sf      = 8;
    SS_ref_db.n_em      = 5;
    SS_ref_db.n_w       = 10;
    SS_ref_db.n_xeos    = 4;
    
    return SS_ref_db;
}


/**************************************************************************************/
/**************************************************************************************/
/*********************Evan&Frost DATABASE (Evans&Frost , 2021)*************************/
/**************************************************************************************/
/**************************************************************************************/

/**
    allocate memory for fluid
*/
SS_ref G_SS_um_fluid_init_function(SS_ref SS_ref_db,  global_variable gv){
    
    SS_ref_db.n_cat     = 0;
    SS_ref_db.is_liq    = 0;
    SS_ref_db.symmetry  = -1;
    SS_ref_db.n_sf      = 2;
    SS_ref_db.n_em      = 2;
    SS_ref_db.n_xeos    = 1;
    
    return SS_ref_db;
}

/**
    allocate memory for ol
*/
SS_ref G_SS_um_ol_init_function(SS_ref SS_ref_db,  global_variable gv){
    
    SS_ref_db.n_cat     = 0;
    SS_ref_db.is_liq    = 0;
    SS_ref_db.symmetry  = 1;
    SS_ref_db.n_sf      = 2;
    SS_ref_db.n_em      = 2;
    SS_ref_db.n_w       = 1;
    SS_ref_db.n_xeos    = 1;
    
    return SS_ref_db;
}

/**
    allocate memory for br
*/
SS_ref G_SS_um_br_init_function(SS_ref SS_ref_db,  global_variable gv){
    
    SS_ref_db.n_cat     = 0;
    SS_ref_db.is_liq    = 0;
    SS_ref_db.symmetry  = -1;
    SS_ref_db.n_sf      = 2;
    SS_ref_db.n_em      = 2;
    SS_ref_db.n_xeos    = 1;
    
    return SS_ref_db;
}

/**
    allocate memory for ch
*/
SS_ref G_SS_um_ch_init_function(SS_ref SS_ref_db,  global_variable gv){
    
    SS_ref_db.n_cat     = 0;
    SS_ref_db.is_liq    = 0;
    SS_ref_db.symmetry  = 1;
    SS_ref_db.n_sf      = 2;
    SS_ref_db.n_em      = 2;
    SS_ref_db.n_w       = 1;
    SS_ref_db.n_xeos    = 1;
    
    return SS_ref_db;
}

/**
    allocate memory for atg
*/
SS_ref G_SS_um_atg_init_function(SS_ref SS_ref_db,  global_variable gv){
    
    SS_ref_db.n_cat     = 0;
    SS_ref_db.is_liq    = 0;
    SS_ref_db.symmetry  = 1;
    SS_ref_db.n_sf      = 8;
    SS_ref_db.n_em      = 5;
    SS_ref_db.n_w       = 10;
    SS_ref_db.n_xeos    = 4;
    
    return SS_ref_db;
}

/**
    allocate memory for g
*/
SS_ref G_SS_um_g_init_function(SS_ref SS_ref_db,  global_variable gv){
    
    SS_ref_db.n_cat     = 0;
    SS_ref_db.is_liq    = 0;
    SS_ref_db.symmetry  = 1;
    SS_ref_db.n_sf      = 2;
    SS_ref_db.n_em      = 2;
    SS_ref_db.n_w       = 1;
    SS_ref_db.n_xeos    = 1;
    
    return SS_ref_db;
}

/**
    allocate memory for ta
*/
SS_ref G_SS_um_ta_init_function(SS_ref SS_ref_db,  global_variable gv){
    
    SS_ref_db.n_cat     = 0;
    SS_ref_db.is_liq    = 0;
    SS_ref_db.symmetry  = 1;
    SS_ref_db.n_sf      = 9;
    SS_ref_db.n_em      = 6;
    SS_ref_db.n_w       = 15;
    SS_ref_db.n_xeos    = 5;
    
    return SS_ref_db;
}

/**
    allocate memory for chl
*/
SS_ref G_SS_um_chl_init_function(SS_ref SS_ref_db,  global_variable gv){
    
    SS_ref_db.n_cat     = 0;
    SS_ref_db.is_liq    = 0;
    SS_ref_db.symmetry  = 1;
    SS_ref_db.n_sf      = 11;
    SS_ref_db.n_em      = 7;
    SS_ref_db.n_w       = 21;
    SS_ref_db.n_xeos    = 6;
    
    return SS_ref_db;
}

/**
    allocate memory for anth
*/
SS_ref G_SS_um_anth_init_function(SS_ref SS_ref_db,  global_variable gv){
    
    SS_ref_db.n_cat     = 0;
    SS_ref_db.is_liq    = 0;
    SS_ref_db.symmetry  = 0;
    SS_ref_db.n_sf      = 9;
    SS_ref_db.n_em      = 5;
    SS_ref_db.n_v       = 5;
    SS_ref_db.n_w       = 10;
    SS_ref_db.n_xeos    = 4;
    
    return SS_ref_db;
}

/**
    allocate memory for spi
*/
SS_ref G_SS_um_spi_init_function(SS_ref SS_ref_db,  global_variable gv){
    
    SS_ref_db.n_cat     = 0;
    SS_ref_db.is_liq    = 0;
    SS_ref_db.symmetry  = 1;
    SS_ref_db.n_sf      = 4;
    SS_ref_db.n_em      = 3;
    SS_ref_db.n_w       = 3;
    SS_ref_db.n_xeos    = 2;
    
    return SS_ref_db;
}

/**
    allocate memory for opx
*/
SS_ref G_SS_um_opx_init_function(SS_ref SS_ref_db,  global_variable gv){
    
    SS_ref_db.n_cat     = 0;
    SS_ref_db.is_liq    = 0;
    SS_ref_db.symmetry  = 1;
    SS_ref_db.n_sf      = 8;
    SS_ref_db.n_em      = 5;
    SS_ref_db.n_w       = 10;
    SS_ref_db.n_xeos    = 4;
    
    return SS_ref_db;
}

/**
    allocate memory for po
*/
SS_ref G_SS_um_po_init_function(SS_ref SS_ref_db,  global_variable gv){
    
    SS_ref_db.n_cat     = 0;
    SS_ref_db.is_liq    = 0;
    SS_ref_db.symmetry  = 1;
    SS_ref_db.n_sf      = 2;
    SS_ref_db.n_em      = 2;
    SS_ref_db.n_w       = 1;
    SS_ref_db.n_xeos    = 1;
    
    return SS_ref_db;
}


/**
    allocate memory for pl4tr
*/
SS_ref G_SS_ume_pl4tr_init_function(SS_ref SS_ref_db, global_variable gv){
    
    SS_ref_db.n_cat     = 0;
    SS_ref_db.is_liq    = 0;
    SS_ref_db.symmetry  = 1;
    SS_ref_db.n_sf      = 4;
    SS_ref_db.n_em      = 2;
    SS_ref_db.n_w       = 1;
    SS_ref_db.n_xeos    = 1;
    
    return SS_ref_db;
}

/**
    allocate memory for hb
*/
SS_ref G_SS_ume_hb_init_function(SS_ref SS_ref_db, global_variable gv){
    
    SS_ref_db.n_cat     = 0;
    SS_ref_db.is_liq    = 0;
    SS_ref_db.symmetry  = 0;
    SS_ref_db.n_sf      = 14;
    SS_ref_db.n_em      = 9;
    SS_ref_db.n_v       = 9;
    SS_ref_db.n_w       = 36;
    SS_ref_db.n_xeos    = 8;
    
    return SS_ref_db;
}

/**
    allocate memory for aug
*/
SS_ref G_SS_ume_aug_init_function(SS_ref SS_ref_db, global_variable gv){
    
    SS_ref_db.n_cat     = 0;
    SS_ref_db.is_liq    = 0;
    SS_ref_db.symmetry  = 0;
    SS_ref_db.n_sf      = 12;
    SS_ref_db.n_em      = 8;
    SS_ref_db.n_v       = 8;
    SS_ref_db.n_w       = 28;
    SS_ref_db.n_xeos    = 7;
    
    return SS_ref_db;
}

/**************************************************************************************/
/**************************************************************************************/
/***********************Mantle DATABASE (Holland et al., 2013)*************************/
/**************************************************************************************/
/**************************************************************************************/


/**
    allocate memory for g
*/
SS_ref G_SS_mtl_g_init_function(SS_ref SS_ref_db,  global_variable gv){
    
    SS_ref_db.n_cat     = 0;
    SS_ref_db.is_liq    = 0;
    SS_ref_db.symmetry  = 1;
    SS_ref_db.n_sf      = 8;
    SS_ref_db.n_em      = 6;
    SS_ref_db.n_w       = 15;
    SS_ref_db.n_xeos    = 5;
    
    return SS_ref_db;
}

/**
    allocate memory for fp
*/
SS_ref G_SS_mtl_fp_init_function(SS_ref SS_ref_db,  global_variable gv){
    
    SS_ref_db.n_cat     = 0;
    SS_ref_db.is_liq    = 0;
    SS_ref_db.symmetry  = 1;
    SS_ref_db.n_sf      = 2;
    SS_ref_db.n_em      = 2;
    SS_ref_db.n_w       = 1;
    SS_ref_db.n_xeos    = 1;
    
    return SS_ref_db;
}

/**
    allocate memory for mpv
*/
SS_ref G_SS_mtl_mpv_init_function(SS_ref SS_ref_db,  global_variable gv){
    
    SS_ref_db.n_cat     = 0;
    SS_ref_db.is_liq    = 0;
    SS_ref_db.symmetry  = 1;
    SS_ref_db.n_sf      = 7;
    SS_ref_db.n_em      = 5;
    SS_ref_db.n_w       = 10;
    SS_ref_db.n_xeos    = 4;
    
    return SS_ref_db;
}
/**
    allocate memory for cpv
*/
SS_ref G_SS_mtl_cpv_init_function(SS_ref SS_ref_db,  global_variable gv){
    
    SS_ref_db.n_cat     = 0;
    SS_ref_db.is_liq    = 0;
    SS_ref_db.symmetry  = 1;
    SS_ref_db.n_sf      = 7;
    SS_ref_db.n_em      = 5;
    SS_ref_db.n_w       = 10;
    SS_ref_db.n_xeos    = 4;
    
    return SS_ref_db;
}

/**
    allocate memory for crn
*/
SS_ref G_SS_mtl_crn_init_function(SS_ref SS_ref_db,  global_variable gv){
    
    SS_ref_db.n_cat     = 0;
    SS_ref_db.is_liq    = 0;
    SS_ref_db.symmetry  = 1;
    SS_ref_db.n_sf      = 5;
    SS_ref_db.n_em      = 3;
    SS_ref_db.n_w       = 3;
    SS_ref_db.n_xeos    = 2;
    
    return SS_ref_db;
}

/**
    allocate memory for cf
*/
SS_ref G_SS_mtl_cf_init_function(SS_ref SS_ref_db,  global_variable gv){
    
    SS_ref_db.n_cat     = 0;
    SS_ref_db.is_liq    = 0;
    SS_ref_db.symmetry  = 1;
    SS_ref_db.n_sf      = 8;
    SS_ref_db.n_em      = 6;
    SS_ref_db.n_w       = 15;
    SS_ref_db.n_xeos    = 5;
    
    return SS_ref_db;
}

/**
    allocate memory for nal
*/
SS_ref G_SS_mtl_nal_init_function(SS_ref SS_ref_db,  global_variable gv){
    
    SS_ref_db.n_cat     = 0;
    SS_ref_db.is_liq    = 0;
    SS_ref_db.symmetry  = 1;
    SS_ref_db.n_sf      = 10;
    SS_ref_db.n_em      = 7;
    SS_ref_db.n_w       = 21;
    SS_ref_db.n_xeos    = 6;
    
    return SS_ref_db;
}

/**
    allocate memory for aki
*/
SS_ref G_SS_mtl_aki_init_function(SS_ref SS_ref_db,  global_variable gv){
    
    SS_ref_db.n_cat     = 0;
    SS_ref_db.is_liq    = 0;
    SS_ref_db.symmetry  = 1;
    SS_ref_db.n_sf      = 5;
    SS_ref_db.n_em      = 3;
    SS_ref_db.n_w       = 3;
    SS_ref_db.n_xeos    = 2;
    
    return SS_ref_db;
}

/**
    allocate memory for ol
*/
SS_ref G_SS_mtl_ol_init_function(SS_ref SS_ref_db,  global_variable gv){
    
    SS_ref_db.n_cat     = 0;
    SS_ref_db.is_liq    = 0;
    SS_ref_db.symmetry  = 1;
    SS_ref_db.n_sf      = 2;
    SS_ref_db.n_em      = 2;
    SS_ref_db.n_w       = 1;
    SS_ref_db.n_xeos    = 1;
    
    return SS_ref_db;
}

/**
    allocate memory for wad
*/
SS_ref G_SS_mtl_wad_init_function(SS_ref SS_ref_db,  global_variable gv){
    
    SS_ref_db.n_cat     = 0;
    SS_ref_db.is_liq    = 0;
    SS_ref_db.symmetry  = 1;
    SS_ref_db.n_sf      = 2;
    SS_ref_db.n_em      = 2;
    SS_ref_db.n_w       = 1;
    SS_ref_db.n_xeos    = 1;
    
    return SS_ref_db;
}

/**
    allocate memory for ring
*/
SS_ref G_SS_mtl_ring_init_function(SS_ref SS_ref_db,  global_variable gv){
    
    SS_ref_db.n_cat     = 0;
    SS_ref_db.is_liq    = 0;
    SS_ref_db.symmetry  = 1;
    SS_ref_db.n_sf      = 2;
    SS_ref_db.n_em      = 2;
    SS_ref_db.n_w       = 1;
    SS_ref_db.n_xeos    = 1;
    
    return SS_ref_db;
}

/**
    allocate memory for cpx
*/
SS_ref G_SS_mtl_cpx_init_function(SS_ref SS_ref_db,  global_variable gv){
    
    SS_ref_db.n_cat     = 0;
    SS_ref_db.is_liq    = 0;
    SS_ref_db.symmetry  = 0;
    SS_ref_db.n_sf      = 9;
    SS_ref_db.n_em      = 6;
    SS_ref_db.n_v       = 6;
    SS_ref_db.n_w       = 15;
    SS_ref_db.n_xeos    = 5;
    
    return SS_ref_db;
}

/**
    allocate memory for opx
*/
SS_ref G_SS_mtl_opx_init_function(SS_ref SS_ref_db,  global_variable gv){
    
    SS_ref_db.n_cat     = 0;
    SS_ref_db.is_liq    = 0;
    SS_ref_db.symmetry  = 0;
    SS_ref_db.n_sf      = 8;
    SS_ref_db.n_em      = 5;
    SS_ref_db.n_v       = 5;
    SS_ref_db.n_w       = 10;
    SS_ref_db.n_xeos    = 4;
    
    return SS_ref_db;
}

/**
    allocate memory for hpx
*/
SS_ref G_SS_mtl_hpx_init_function(SS_ref SS_ref_db,  global_variable gv){
    
    SS_ref_db.n_cat     = 0;
    SS_ref_db.is_liq    = 0;
    SS_ref_db.symmetry  = 0;
    SS_ref_db.n_sf      = 8;
    SS_ref_db.n_em      = 5;
    SS_ref_db.n_v       = 5;
    SS_ref_db.n_w       = 10;
    SS_ref_db.n_xeos    = 4;
    
    return SS_ref_db;
}


/**************************************************************************************/
/**************************************************************************************/
/*   Metapelite ext DB (White et al., 2014; Green et al., 2016; Evans & Forst, 2021)  */
/**************************************************************************************/
/**************************************************************************************/
/**
    allocate memory for g
*/
SS_ref G_SS_mpe_g_init_function(SS_ref SS_ref_db, global_variable gv){
    
    SS_ref_db.n_cat     = 0;
    SS_ref_db.is_liq    = 0;
    SS_ref_db.symmetry  = 0;
    SS_ref_db.n_sf      = 6;
    SS_ref_db.n_em      = 5;
    SS_ref_db.n_v       = 5;
    SS_ref_db.n_w       = 10;
    SS_ref_db.n_xeos    = 4;
    
    return SS_ref_db;
}

/**
    allocate memory for L
*/
SS_ref G_SS_mpe_liq_init_function(SS_ref SS_ref_db, global_variable gv){
    
    SS_ref_db.n_cat     = 0;
    SS_ref_db.is_liq    = 0;
    SS_ref_db.symmetry  = 1;
    SS_ref_db.n_sf      = 10;
    SS_ref_db.n_em      = 8;
    SS_ref_db.n_w       = 28;
    SS_ref_db.n_xeos    = 7;
    
    return SS_ref_db;
}

/**
    allocate memory for pl4tr
*/
SS_ref G_SS_mpe_fsp_init_function(SS_ref SS_ref_db, global_variable gv){
    
    SS_ref_db.n_cat     = 0;
    SS_ref_db.is_liq    = 0;
    SS_ref_db.symmetry  = 0;
    SS_ref_db.n_sf      = 5;
    SS_ref_db.n_em      = 3;
    SS_ref_db.n_v       = 3;
    SS_ref_db.n_w       = 3;
    SS_ref_db.n_xeos    = 2;
    
    return SS_ref_db;
}

/**
    allocate memory for ksp
*/
SS_ref G_SS_mpe_ksp_init_function(SS_ref SS_ref_db, global_variable gv){
    
    SS_ref_db.n_cat     = 0;
    SS_ref_db.is_liq    = 0;
    SS_ref_db.symmetry  = 0;
    SS_ref_db.n_sf      = 3;
    SS_ref_db.n_em      = 3;
    SS_ref_db.n_v       = 3;
    SS_ref_db.n_w       = 3;
    SS_ref_db.n_xeos    = 2;
    
    return SS_ref_db;
}

/**
    allocate memory for ep
*/
SS_ref G_SS_mpe_ep_init_function(SS_ref SS_ref_db, global_variable gv){
    
    SS_ref_db.n_cat     = 0;
    SS_ref_db.is_liq    = 0;
    SS_ref_db.symmetry  = 1;
    SS_ref_db.n_sf      = 4;
    SS_ref_db.n_em      = 3;
    SS_ref_db.n_w       = 3;
    SS_ref_db.n_xeos    = 2;
    
    return SS_ref_db;
}

/**
    allocate memory for ma
*/
SS_ref G_SS_mpe_ma_init_function(SS_ref SS_ref_db, global_variable gv){
    
    SS_ref_db.n_cat     = 0;
    SS_ref_db.is_liq    = 0;
    SS_ref_db.symmetry  = 0;
    SS_ref_db.n_sf      = 10;
    SS_ref_db.n_em      = 6;
    SS_ref_db.n_v       = 6;
    SS_ref_db.n_w       = 15;
    SS_ref_db.n_xeos    = 5;
    
    return SS_ref_db;
}

/**
    allocate memory for mu
*/
SS_ref G_SS_mpe_mu_init_function(SS_ref SS_ref_db, global_variable gv){
    
    SS_ref_db.n_cat     = 0;
    SS_ref_db.is_liq    = 0;
    SS_ref_db.symmetry  = 0;
    SS_ref_db.n_sf      = 10;
    SS_ref_db.n_em      = 6;
    SS_ref_db.n_v       = 6;
    SS_ref_db.n_w       = 15;
    SS_ref_db.n_xeos    = 5;
    
    return SS_ref_db;
}
/**
    allocate memory for bi_mpe
*/
SS_ref G_SS_mpe_bi_init_function(SS_ref SS_ref_db,  global_variable gv){
    
    SS_ref_db.n_cat     = 0;
    SS_ref_db.is_liq    = 0;
    SS_ref_db.symmetry  = 1;
    SS_ref_db.n_sf      = 13;
    SS_ref_db.n_em      = 7;
    SS_ref_db.n_w       = 21;
    SS_ref_db.n_xeos    = 6;
    
    return SS_ref_db;
}

/**
    allocate memory for opx
*/
SS_ref G_SS_mpe_opx_init_function(SS_ref SS_ref_db, global_variable gv){
    
    SS_ref_db.n_cat     = 0;
    SS_ref_db.is_liq    = 0;
    SS_ref_db.symmetry  = 0;
    SS_ref_db.n_sf      = 10;
    SS_ref_db.n_em      = 7;
    SS_ref_db.n_v       = 7;
    SS_ref_db.n_w       = 21;
    SS_ref_db.n_xeos    = 6;
    
    return SS_ref_db;
}

/**
    allocate memory for sa
*/
SS_ref G_SS_mpe_sa_init_function(SS_ref SS_ref_db, global_variable gv){
    
    SS_ref_db.n_cat     = 0;
    SS_ref_db.is_liq    = 0;
    SS_ref_db.symmetry  = 1;
    SS_ref_db.n_sf      = 8;
    SS_ref_db.n_em      = 5;
    SS_ref_db.n_w       = 10;
    SS_ref_db.n_xeos    = 4;
    
    return SS_ref_db;
}

/**
    allocate memory for cd
*/
SS_ref G_SS_mpe_cd_init_function(SS_ref SS_ref_db, global_variable gv){
    
    SS_ref_db.n_cat     = 0;
    SS_ref_db.is_liq    = 0;
    SS_ref_db.symmetry  = 1;
    SS_ref_db.n_sf      = 5;
    SS_ref_db.n_em      = 4;
    SS_ref_db.n_w       = 6;
    SS_ref_db.n_xeos    = 3;
    
    return SS_ref_db;
}

/**
    allocate memory for st
*/
SS_ref G_SS_mpe_st_init_function(SS_ref SS_ref_db, global_variable gv){
    
    SS_ref_db.n_cat     = 0;
    SS_ref_db.is_liq    = 0;
    SS_ref_db.symmetry  = 1;
    SS_ref_db.n_sf      = 7;
    SS_ref_db.n_em      = 5;
    SS_ref_db.n_w       = 10;
    SS_ref_db.n_xeos    = 4;
    
    return SS_ref_db;
}

/**
    allocate memory for chl
*/
SS_ref G_SS_mpe_chl_init_function(SS_ref SS_ref_db, global_variable gv){
    
    SS_ref_db.n_cat     = 0;
    SS_ref_db.is_liq    = 0;
    SS_ref_db.symmetry  = 1;
    SS_ref_db.n_sf      = 12;
    SS_ref_db.n_em      = 8;
    SS_ref_db.n_w       = 28;
    SS_ref_db.n_xeos    = 7;
    
    return SS_ref_db;
}

/**
    allocate memory for ctd
*/
SS_ref G_SS_mpe_ctd_init_function(SS_ref SS_ref_db, global_variable gv){
    
    SS_ref_db.n_cat     = 0;
    SS_ref_db.is_liq    = 0;
    SS_ref_db.symmetry  = 1;
    SS_ref_db.n_sf      = 5;
    SS_ref_db.n_em      = 4;
    SS_ref_db.n_w       = 6;
    SS_ref_db.n_xeos    = 3;
    
    return SS_ref_db;
}

/**
    allocate memory for sp
*/
SS_ref G_SS_mpe_sp_init_function(SS_ref SS_ref_db, global_variable gv){
    
    SS_ref_db.n_cat     = 0;
    SS_ref_db.is_liq    = 0;
    SS_ref_db.symmetry  = 1;
    SS_ref_db.n_sf      = 5;
    SS_ref_db.n_em      = 4;
    SS_ref_db.n_w       = 6;
    SS_ref_db.n_xeos    = 3;
    
    return SS_ref_db;
}

/**
    allocate memory for ilmm
*/
SS_ref G_SS_mpe_ilmm_init_function(SS_ref SS_ref_db, global_variable gv){
    
    SS_ref_db.n_cat     = 0;
    SS_ref_db.is_liq    = 0;
    SS_ref_db.symmetry  = 1;
    SS_ref_db.n_sf      = 7;
    SS_ref_db.n_em      = 5;
    SS_ref_db.n_w       = 10;
    SS_ref_db.n_xeos    = 4;
    
    return SS_ref_db;
}

/**
    allocate memory for ilm
*/
SS_ref G_SS_mpe_ilm_init_function(SS_ref SS_ref_db, global_variable gv){
    
    SS_ref_db.n_cat     = 0;
    SS_ref_db.is_liq    = 0;
    SS_ref_db.symmetry  = 1;
    SS_ref_db.n_sf      = 6;
    SS_ref_db.n_em      = 3;
    SS_ref_db.n_w       = 3;
    SS_ref_db.n_xeos    = 2;
    
    return SS_ref_db;
}

/**
    allocate memory for mt1
*/
SS_ref G_SS_mpe_mt_init_function(SS_ref SS_ref_db, global_variable gv){
    
    SS_ref_db.n_cat     = 0;
    SS_ref_db.is_liq    = 0;
    SS_ref_db.symmetry  = 1;
    SS_ref_db.n_sf      = 5;
    SS_ref_db.n_em      = 3;
    SS_ref_db.n_w       = 3;
    SS_ref_db.n_xeos    = 2;
    
    return SS_ref_db;
}

/**
    allocate memory for fl
*/
SS_ref G_SS_mpe_fl_init_function(SS_ref SS_ref_db, global_variable gv){
    
    SS_ref_db.n_cat     = 0;
    SS_ref_db.is_liq    = 0;
    SS_ref_db.symmetry  = 1;
    SS_ref_db.n_sf      = 2;
    SS_ref_db.n_em      = 2;
    SS_ref_db.n_v       = 0;
    SS_ref_db.n_w       = 0;
    SS_ref_db.n_xeos    = 1;
    
    return SS_ref_db;
}

/**
    allocate memory for occm
*/
SS_ref G_SS_mpe_occm_init_function(SS_ref SS_ref_db, global_variable gv){
    
    SS_ref_db.n_cat     = 0;
    SS_ref_db.is_liq    = 0;
    SS_ref_db.symmetry  = 0;
    SS_ref_db.n_sf      = 9;
    SS_ref_db.n_em      = 5;
    SS_ref_db.n_v       = 5;
    SS_ref_db.n_w       = 10;
    SS_ref_db.n_xeos    = 4;
    
    return SS_ref_db;
}

/**
    allocate memory for po
*/
SS_ref G_SS_mpe_po_init_function(SS_ref SS_ref_db, global_variable gv){
    
    SS_ref_db.n_cat     = 0;
    SS_ref_db.is_liq    = 0;
    SS_ref_db.symmetry  = 1;
    SS_ref_db.n_sf      = 2;
    SS_ref_db.n_em      = 2;
    SS_ref_db.n_w       = 1;
    SS_ref_db.n_xeos    = 1;
    
    return SS_ref_db;
}

/**
    allocate memory for hb
*/
SS_ref G_SS_mpe_hb_init_function(SS_ref SS_ref_db, global_variable gv){
    
    SS_ref_db.n_cat     = 0;
    SS_ref_db.is_liq    = 0;
    SS_ref_db.symmetry  = 0;
    SS_ref_db.n_sf      = 18;
    SS_ref_db.n_em      = 11;
    SS_ref_db.n_v       = 11;
    SS_ref_db.n_w       = 55;
    SS_ref_db.n_xeos    = 10;
    
    return SS_ref_db;
}

/**
    allocate memory for aug
*/
SS_ref G_SS_mpe_aug_init_function(SS_ref SS_ref_db, global_variable gv){
    
    SS_ref_db.n_cat     = 0;
    SS_ref_db.is_liq    = 0;
    SS_ref_db.symmetry  = 0;
    SS_ref_db.n_sf      = 12;
    SS_ref_db.n_em      = 8;
    SS_ref_db.n_v       = 8;
    SS_ref_db.n_w       = 28;
    SS_ref_db.n_xeos    = 7;
    
    return SS_ref_db;
}

/**
    allocate memory for dio
*/
SS_ref G_SS_mpe_dio_init_function(SS_ref SS_ref_db, global_variable gv){
    
    SS_ref_db.n_cat     = 0;
    SS_ref_db.is_liq    = 0;
    SS_ref_db.symmetry  = 1;
    SS_ref_db.n_sf      = 12;
    SS_ref_db.n_em      = 7;
    SS_ref_db.n_w       = 21;
    SS_ref_db.n_xeos    = 6;
    
    return SS_ref_db;
}

/**************************************************************************************/
/**************************************************************************************/
/**************************************************************************************/
/**************************************************************************************/

void TC_SS_init_mp(	                SS_init_type 		*SS_init,
									global_variable 	 gv				){

	for (int iss = 0; iss < gv.len_ss; iss++){

		if      (strcmp( gv.SS_list[iss], "liq")   == 0 ){
			SS_init[iss]  = G_SS_mp_liq_init_function; 		}
		else if (strcmp( gv.SS_list[iss], "fsp") == 0){
			SS_init[iss]  = G_SS_mp_fsp_init_function; 		}
		else if (strcmp( gv.SS_list[iss], "bi")    == 0){
			SS_init[iss]  = G_SS_mp_bi_init_function; 		}
		else if (strcmp( gv.SS_list[iss], "g")     == 0){
			SS_init[iss]  = G_SS_mp_g_init_function; 		}
		else if (strcmp( gv.SS_list[iss], "ep")    == 0){
			SS_init[iss]  = G_SS_mp_ep_init_function; 		}
		else if (strcmp( gv.SS_list[iss], "ma")    == 0){
			SS_init[iss]  = G_SS_mp_ma_init_function; 		}
		else if (strcmp( gv.SS_list[iss], "mu")    == 0){
			SS_init[iss]  = G_SS_mp_mu_init_function; 		}
		else if (strcmp( gv.SS_list[iss], "opx")   == 0){
			SS_init[iss]  = G_SS_mp_opx_init_function; 		}
		else if (strcmp( gv.SS_list[iss], "sa")    == 0){
			SS_init[iss]  = G_SS_mp_sa_init_function; 		}
		else if (strcmp( gv.SS_list[iss], "cd")    == 0){
			SS_init[iss]  = G_SS_mp_cd_init_function; 		}
		else if (strcmp( gv.SS_list[iss], "st")    == 0){
			SS_init[iss]  = G_SS_mp_st_init_function; 		}
		else if (strcmp( gv.SS_list[iss], "chl")   == 0){
			SS_init[iss]  = G_SS_mp_chl_init_function; 		}
		else if (strcmp( gv.SS_list[iss], "ctd")   == 0){
			SS_init[iss]  = G_SS_mp_ctd_init_function; 		}
		else if (strcmp( gv.SS_list[iss], "sp")    == 0){
			SS_init[iss]  = G_SS_mp_sp_init_function; 		}
		else if (strcmp( gv.SS_list[iss], "ilm")   == 0){
			SS_init[iss]  = G_SS_mp_ilm_init_function; 		}
		else if (strcmp( gv.SS_list[iss], "ilmm")   == 0){
			SS_init[iss]  = G_SS_mp_ilmm_init_function; 	}
		else if (strcmp( gv.SS_list[iss], "mt")    == 0){
			SS_init[iss]  = G_SS_mp_mt_init_function; 		}
		else if (strcmp( gv.SS_list[iss], "aq17")  == 0){
			SS_init[iss]  = G_SS_aq17_init_function; 	    }
		else{
			printf("\nsolid solution '%s' is not in the database, cannot be initiated\n", gv.SS_list[iss]);	
		}	
	};	
}

void TC_SS_init_mb(	                SS_init_type 		*SS_init,
									global_variable 	 gv				){
					 
	for (int iss = 0; iss < gv.len_ss; iss++){
        if (strcmp( gv.SS_list[iss], "liq")  == 0){
            SS_init[iss]  = G_SS_mb_liq_init_function;        }
        else if (strcmp( gv.SS_list[iss], "hb")  == 0){
            SS_init[iss]  = G_SS_mb_hb_init_function;         }
        else if (strcmp( gv.SS_list[iss], "aug")  == 0){
            SS_init[iss]  = G_SS_mb_aug_init_function;        }
        else if (strcmp( gv.SS_list[iss], "dio")  == 0){
            SS_init[iss]  = G_SS_mb_dio_init_function;        }
        else if (strcmp( gv.SS_list[iss], "opx")  == 0){
            SS_init[iss]  = G_SS_mb_opx_init_function;        }
        else if (strcmp( gv.SS_list[iss], "g")  == 0){
            SS_init[iss]  = G_SS_mb_g_init_function;          }
        else if (strcmp( gv.SS_list[iss], "ol")  == 0){
            SS_init[iss]  = G_SS_mb_ol_init_function;         }
        else if (strcmp( gv.SS_list[iss], "fsp")  == 0){
            SS_init[iss]  = G_SS_mb_fsp_init_function;        }
        else if (strcmp( gv.SS_list[iss], "abc")  == 0){
            SS_init[iss]  = G_SS_mb_abc_init_function;        }
        else if (strcmp( gv.SS_list[iss], "k4tr")  == 0){
            SS_init[iss]  = G_SS_mb_k4tr_init_function;       }
        else if (strcmp( gv.SS_list[iss], "sp")  == 0){
            SS_init[iss]  = G_SS_mb_sp_init_function;         }
        else if (strcmp( gv.SS_list[iss], "spn")  == 0){
            SS_init[iss]  = G_SS_mb_spn_init_function;        }
        else if (strcmp( gv.SS_list[iss], "ilm")  == 0){
            SS_init[iss]  = G_SS_mb_ilm_init_function;        }
        else if (strcmp( gv.SS_list[iss], "ilmm")  == 0){
            SS_init[iss]  = G_SS_mb_ilmm_init_function;       }
        else if (strcmp( gv.SS_list[iss], "ep")  == 0){
            SS_init[iss]  = G_SS_mb_ep_init_function;         }
        else if (strcmp( gv.SS_list[iss], "bi")  == 0){
            SS_init[iss]  = G_SS_mb_bi_init_function;         }
        else if (strcmp( gv.SS_list[iss], "mu")  == 0){
            SS_init[iss]  = G_SS_mb_mu_init_function;         }
        else if (strcmp( gv.SS_list[iss], "chl")  == 0){
            SS_init[iss]  = G_SS_mb_chl_init_function;        }
		else{
			printf("\nsolid solution '%s' is not in the database, cannot be initiated\n", gv.SS_list[iss]);	
		}	
	};	
}


void TC_SS_init_ig(	                SS_init_type 		*SS_init,
									global_variable 	 gv				){
					 
	for (int iss = 0; iss < gv.len_ss; iss++){

		if      (strcmp( gv.SS_list[iss], "bi")  == 0 ){
			SS_init[iss]  = G_SS_ig_bi_init_function; 		}
		else if (strcmp( gv.SS_list[iss], "fper")  == 0){
			SS_init[iss]  = G_SS_ig_fper_init_function; 		}
		else if (strcmp( gv.SS_list[iss], "cd")  == 0){
			SS_init[iss]  = G_SS_ig_cd_init_function; 		}
		else if (strcmp( gv.SS_list[iss], "cpx") == 0){
			SS_init[iss]  = G_SS_ig_cpx_init_function; 		}
		else if (strcmp( gv.SS_list[iss], "ep")  == 0){
			SS_init[iss]  = G_SS_ig_ep_init_function; 		}
		else if (strcmp( gv.SS_list[iss], "fl")  == 0){
			SS_init[iss]  = G_SS_ig_fl_init_function; 		}
		else if (strcmp( gv.SS_list[iss], "g")   == 0){
			SS_init[iss]  = G_SS_ig_g_init_function; 		    }
		else if (strcmp( gv.SS_list[iss], "hb")  == 0){
			SS_init[iss]  = G_SS_ig_hb_init_function; 		}
		else if (strcmp( gv.SS_list[iss], "ilm") == 0){
			SS_init[iss]  = G_SS_ig_ilm_init_function; 		}
		else if (strcmp( gv.SS_list[iss], "liq") == 0){
			SS_init[iss]  = G_SS_ig_liq_init_function; 		}
		else if (strcmp( gv.SS_list[iss], "mu")  == 0){
			SS_init[iss]  = G_SS_ig_mu_init_function; 		}
		else if (strcmp( gv.SS_list[iss], "ol")  == 0){
			SS_init[iss]  = G_SS_ig_ol_init_function; 		}
		else if (strcmp( gv.SS_list[iss], "opx") == 0){
			SS_init[iss]  = G_SS_ig_opx_init_function; 		}
		else if (strcmp( gv.SS_list[iss], "fsp") == 0){
			SS_init[iss]  = G_SS_ig_fsp_init_function; 		}
		else if (strcmp( gv.SS_list[iss], "spn") == 0){
			SS_init[iss]  = G_SS_ig_spn_init_function; 		}
		else{
			printf("\nsolid solution '%s' is not in the database, cannot be initiated\n", gv.SS_list[iss]);	
		}	
	};	
}


void TC_SS_init_igad(	            SS_init_type 		*SS_init,
									global_variable 	 gv				){
					 
	for (int iss = 0; iss < gv.len_ss; iss++){

        if (strcmp( gv.SS_list[iss], "cpx") == 0){
			SS_init[iss]  = G_SS_igad_cpx_init_function; 		}
		else if (strcmp( gv.SS_list[iss], "g")   == 0){
			SS_init[iss]  = G_SS_igad_g_init_function; 		    }
		else if (strcmp( gv.SS_list[iss], "ilm") == 0){
			SS_init[iss]  = G_SS_igad_ilm_init_function; 		}
		else if (strcmp( gv.SS_list[iss], "liq") == 0){
			SS_init[iss]  = G_SS_igad_liq_init_function; 		}
		else if (strcmp( gv.SS_list[iss], "ol")  == 0){
			SS_init[iss]  = G_SS_igad_ol_init_function; 		}
		else if (strcmp( gv.SS_list[iss], "opx") == 0){
			SS_init[iss]  = G_SS_igad_opx_init_function; 		}
		else if (strcmp( gv.SS_list[iss], "fsp") == 0){
			SS_init[iss]  = G_SS_igad_fsp_init_function; 		}
		else if (strcmp( gv.SS_list[iss], "spn") == 0){
			SS_init[iss]  = G_SS_igad_spn_init_function; 		}
		else if (strcmp( gv.SS_list[iss], "ness") == 0){
			SS_init[iss]  = G_SS_igad_ness_init_function; 		}
		else if (strcmp( gv.SS_list[iss], "lct") == 0){
			SS_init[iss]  = G_SS_igad_lct_init_function; 		}
		else if (strcmp( gv.SS_list[iss], "kals") == 0){
			SS_init[iss]  = G_SS_igad_kals_init_function; 		}
		else if (strcmp( gv.SS_list[iss], "mel") == 0){
			SS_init[iss]  = G_SS_igad_mel_init_function; 		}
		else{
			printf("\nsolid solution '%s' is not in the database, cannot be initiated\n", gv.SS_list[iss]);	
		}	
	};	
}

void TC_SS_init_um(	                SS_init_type 		*SS_init,
									global_variable 	 gv				){
						 
	for (int iss = 0; iss < gv.len_ss; iss++){

		if      (strcmp( gv.SS_list[iss], "fl")  == 0 ){
			SS_init[iss]  = G_SS_um_fluid_init_function; 		}
		else if (strcmp( gv.SS_list[iss], "ol")  == 0){
			SS_init[iss]  = G_SS_um_ol_init_function; 		}
		else if (strcmp( gv.SS_list[iss], "br") == 0){
			SS_init[iss]  = G_SS_um_br_init_function; 		}
		else if (strcmp( gv.SS_list[iss], "ch")  == 0){
			SS_init[iss]  = G_SS_um_ch_init_function; 		}
		else if (strcmp( gv.SS_list[iss], "atg")  == 0){
			SS_init[iss]  = G_SS_um_atg_init_function; 		}
		else if (strcmp( gv.SS_list[iss], "g")   == 0){
			SS_init[iss]  = G_SS_um_g_init_function; 		    }
		else if (strcmp( gv.SS_list[iss], "ta")  == 0){
			SS_init[iss]  = G_SS_um_ta_init_function; 		}
		else if (strcmp( gv.SS_list[iss], "chl") == 0){
			SS_init[iss]  = G_SS_um_chl_init_function; 		}
		else if (strcmp( gv.SS_list[iss], "anth") == 0){
			SS_init[iss]  = G_SS_um_anth_init_function; 		}
		else if (strcmp( gv.SS_list[iss], "spi")  == 0){
			SS_init[iss]  = G_SS_um_spi_init_function; 		}
		else if (strcmp( gv.SS_list[iss], "opx") == 0){
			SS_init[iss]  = G_SS_um_opx_init_function; 		}
		else if (strcmp( gv.SS_list[iss], "po") == 0){
			SS_init[iss]  = G_SS_um_po_init_function; 		}
		else{
			printf("\nsolid solution '%s' is not in the database, cannot be initiated\n", gv.SS_list[iss]);	
		}	
	};						 
}


void TC_SS_init_um_ext(	            SS_init_type 		*SS_init,
									global_variable 	 gv				){
						 
	for (int iss = 0; iss < gv.len_ss; iss++){

		if      (strcmp( gv.SS_list[iss], "fl")  == 0 ){
			SS_init[iss]  = G_SS_um_fluid_init_function; 		}
		else if (strcmp( gv.SS_list[iss], "ol")  == 0){
			SS_init[iss]  = G_SS_um_ol_init_function; 		}
		else if (strcmp( gv.SS_list[iss], "br") == 0){
			SS_init[iss]  = G_SS_um_br_init_function; 		}
		else if (strcmp( gv.SS_list[iss], "ch")  == 0){
			SS_init[iss]  = G_SS_um_ch_init_function; 		}
		else if (strcmp( gv.SS_list[iss], "atg")  == 0){
			SS_init[iss]  = G_SS_um_atg_init_function; 		}
		else if (strcmp( gv.SS_list[iss], "g")   == 0){
			SS_init[iss]  = G_SS_um_g_init_function; 	}
		else if (strcmp( gv.SS_list[iss], "ta")  == 0){
			SS_init[iss]  = G_SS_um_ta_init_function; 		}
		else if (strcmp( gv.SS_list[iss], "chl") == 0){
			SS_init[iss]  = G_SS_um_chl_init_function; 		}
		else if (strcmp( gv.SS_list[iss], "anth") == 0){
			SS_init[iss]  = G_SS_um_anth_init_function; 	}
		else if (strcmp( gv.SS_list[iss], "spi")  == 0){
			SS_init[iss]  = G_SS_um_spi_init_function; 		}
		else if (strcmp( gv.SS_list[iss], "opx") == 0){
			SS_init[iss]  = G_SS_um_opx_init_function; 		}
		else if (strcmp( gv.SS_list[iss], "po") == 0){
			SS_init[iss]  = G_SS_um_po_init_function; 		}
		else if (strcmp( gv.SS_list[iss], "pl4tr")  == 0){
			SS_init[iss]  = G_SS_ume_pl4tr_init_function; 	}
		else if (strcmp( gv.SS_list[iss], "hb") == 0){
			SS_init[iss]  = G_SS_ume_hb_init_function; 		}
		else if (strcmp( gv.SS_list[iss], "aug") == 0){
			SS_init[iss]  = G_SS_ume_aug_init_function; 	}
		else{
			printf("\nsolid solution '%s' is not in the database, cannot be initiated\n", gv.SS_list[iss]);	
		}	
	};						 
}

void TC_SS_init_mtl(	            SS_init_type 		*SS_init,
									global_variable 	 gv				){
						 
	for (int iss = 0; iss < gv.len_ss; iss++){

		if      (strcmp( gv.SS_list[iss], "g")  == 0 ){
			SS_init[iss]  = G_SS_mtl_g_init_function; 		}
		else if (strcmp( gv.SS_list[iss], "fp")  == 0){
			SS_init[iss]  = G_SS_mtl_fp_init_function; 		}
		else if (strcmp( gv.SS_list[iss], "mpv") == 0){
			SS_init[iss]  = G_SS_mtl_mpv_init_function; 	}
		else if (strcmp( gv.SS_list[iss], "cpv") == 0){
			SS_init[iss]  = G_SS_mtl_cpv_init_function; 	}
		else if (strcmp( gv.SS_list[iss], "crn")  == 0){
			SS_init[iss]  = G_SS_mtl_crn_init_function; 	}
		else if (strcmp( gv.SS_list[iss], "cf")  == 0){
			SS_init[iss]  = G_SS_mtl_cf_init_function; 		}
		else if (strcmp( gv.SS_list[iss], "nal")   == 0){
			SS_init[iss]  = G_SS_mtl_nal_init_function; 	}
		else if (strcmp( gv.SS_list[iss], "aki")  == 0){
			SS_init[iss]  = G_SS_mtl_aki_init_function; 	}
		else if (strcmp( gv.SS_list[iss], "ol") == 0){
			SS_init[iss]  = G_SS_mtl_ol_init_function; 		}
		else if (strcmp( gv.SS_list[iss], "wad") == 0){
			SS_init[iss]  = G_SS_mtl_wad_init_function; 	}
		else if (strcmp( gv.SS_list[iss], "ring")  == 0){
			SS_init[iss]  = G_SS_mtl_ring_init_function; 	}
		else if (strcmp( gv.SS_list[iss], "cpx") == 0){
			SS_init[iss]  = G_SS_mtl_cpx_init_function; 	}
		else if (strcmp( gv.SS_list[iss], "opx") == 0){
			SS_init[iss]  = G_SS_mtl_opx_init_function; 	}
		else if (strcmp( gv.SS_list[iss], "hpx")  == 0){
			SS_init[iss]  = G_SS_mtl_hpx_init_function; 	}
		else{
			printf("\nsolid solution '%s' is not in the database, cannot be initiated\n", gv.SS_list[iss]);	
		}	
	};						 
}


void TC_SS_init_mp_ext(	            SS_init_type 		*SS_init,
									global_variable 	 gv				){

	for (int iss = 0; iss < gv.len_ss; iss++){

		if      (strcmp( gv.SS_list[iss], "liq")   == 0 ){
			SS_init[iss]  = G_SS_mpe_liq_init_function; 	}
		else if (strcmp( gv.SS_list[iss], "fsp") == 0){
			SS_init[iss]  = G_SS_mpe_fsp_init_function; 	}
		else if (strcmp( gv.SS_list[iss], "bi")    == 0){
			SS_init[iss]  = G_SS_mpe_bi_init_function; 		}
		else if (strcmp( gv.SS_list[iss], "g")     == 0){
			SS_init[iss]  = G_SS_mpe_g_init_function; 		}
		else if (strcmp( gv.SS_list[iss], "ep")    == 0){
			SS_init[iss]  = G_SS_mpe_ep_init_function; 		}
		else if (strcmp( gv.SS_list[iss], "ma")    == 0){
			SS_init[iss]  = G_SS_mpe_ma_init_function; 		}
		else if (strcmp( gv.SS_list[iss], "mu")    == 0){
			SS_init[iss]  = G_SS_mpe_mu_init_function; 		}
		else if (strcmp( gv.SS_list[iss], "opx")   == 0){
			SS_init[iss]  = G_SS_mpe_opx_init_function;     }
		else if (strcmp( gv.SS_list[iss], "sa")    == 0){
			SS_init[iss]  = G_SS_mpe_sa_init_function; 		}
		else if (strcmp( gv.SS_list[iss], "cd")    == 0){
			SS_init[iss]  = G_SS_mpe_cd_init_function; 		}
		else if (strcmp( gv.SS_list[iss], "st")    == 0){
			SS_init[iss]  = G_SS_mpe_st_init_function; 		}
		else if (strcmp( gv.SS_list[iss], "chl")   == 0){
			SS_init[iss]  = G_SS_mpe_chl_init_function; 	}
		else if (strcmp( gv.SS_list[iss], "ctd")   == 0){
			SS_init[iss]  = G_SS_mpe_ctd_init_function; 	}
		else if (strcmp( gv.SS_list[iss], "sp")    == 0){
			SS_init[iss]  = G_SS_mpe_sp_init_function; 		}
		else if (strcmp( gv.SS_list[iss], "ilm")   == 0){
			SS_init[iss]  = G_SS_mpe_ilm_init_function; 	}
		else if (strcmp( gv.SS_list[iss], "ilmm")   == 0){
			SS_init[iss]  = G_SS_mpe_ilmm_init_function; 	}
		else if (strcmp( gv.SS_list[iss], "mt")    == 0){
			SS_init[iss]  = G_SS_mpe_mt_init_function; 		}
		else if (strcmp( gv.SS_list[iss], "occm")   == 0){
			SS_init[iss]  = G_SS_mpe_occm_init_function; 	}
		else if (strcmp( gv.SS_list[iss], "fl")   == 0){
			SS_init[iss]  = G_SS_mpe_fl_init_function; 	}
		else if (strcmp( gv.SS_list[iss], "po")    == 0){
			SS_init[iss]  = G_SS_mpe_po_init_function; 		}
		else if (strcmp( gv.SS_list[iss], "dio")   == 0){
			SS_init[iss]  = G_SS_mpe_dio_init_function; 	}
		else if (strcmp( gv.SS_list[iss], "aug")   == 0){
			SS_init[iss]  = G_SS_mpe_aug_init_function; 	}
		else if (strcmp( gv.SS_list[iss], "hb")    == 0){
			SS_init[iss]  = G_SS_mpe_hb_init_function; 		}
		else{
			printf("\nsolid solution '%s' is not in the database, cannot be initiated\n", gv.SS_list[iss]);	
		}	
	};	
}

void TC_SS_init(	        	    SS_init_type 		*SS_init,
									global_variable 	 gv				){

	if (gv.EM_database == 0){				// metapelite database //
		TC_SS_init_mp(	 				    SS_init,
											gv							);
	}
	if (gv.EM_database == 1){				// metabasite database //
		TC_SS_init_mb(	 				    SS_init,
											gv							);
	}
	else if (gv.EM_database == 2){			// igneous database //
		TC_SS_init_ig(	 				    SS_init,
											gv							);
	}
	else if (gv.EM_database == 3){			// igneous alkaline database //
		TC_SS_init_igad(	 				SS_init,
											gv							);
	}
	else if (gv.EM_database == 4){			// ultramafic database //
		TC_SS_init_um(	 				    SS_init,
											gv							);
	}
	else if (gv.EM_database == 5){			// ultramafic database //
		TC_SS_init_um_ext(	 				SS_init,
											gv							);
	}
	else if (gv.EM_database == 6){			// mantle database //
		TC_SS_init_mtl(	 				    SS_init,
											gv							);
	}
    else if (gv.EM_database == 7){			// metapelite ext database //
		TC_SS_init_mp_ext(	 				SS_init,
											gv							);
	}
}