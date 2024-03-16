/**

Function to allocate memory for solid-solutions        

Holland et al., 2018 - Melting of peridotites through to granites                 
Igneous dataset to use with tc-ds633.txt                                          
"bi","cpx","cd","ep","fl","g","hb","ilm","liq","mu", "ol", "opx","fsp","spn" 
 
PP & SS_flags
-------------
Declare flags needed for leveling and pge algorithms

+-------+-------+------+------+-------+-------+------+
| SS/PP |  IN	| ACT  | HLD  |  RMV  |  CYC  | REIN |
+=======+=======+======+======+=======+=======+======+
| [0]   | 0/1   |  0/1 | 0/1  | 0/1   | 0/n   | 0/1  |
+-------+-------+------+------+-------+-------+------+
| [1]   | 0/1   |  0/1 | 0/1  | 0/1   | 0/n   | 0/1  |
+-------+-------+------+------+-------+-------+------+
| [2]   | 0/1   |  0/1 | 0/1  | 0/1   | 0/n   | 0/1  |
+-------+-------+------+------+-------+-------+------+
|  ...  | 0/1   |  0/1 | 0/1  | 0/1   | 0/n   | 0/1  |
+-------+-------+------+------+-------+-------+------+
| [m]	| 0/1   |  0/1 | 0/1  | 0/1   | 0/n   | 0/1  |
+-------+-------+------+------+-------+-------+------+
													 
- IN:   allowed phase (satisfying bulk rock constraints)
- ACT:  considered phase (part of the active set of phases)
- HLD:  on hold (not in the active set but still scanned at every iteration)
- RMV:  removed (not considered anymore)
- REIN: phase reintroduced
- m: number of PP or SS
- n: number of cycles

*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <complex.h> 

#include "MAGEMin.h"
#include "gss_init_function.h"

/** 
  allocate memory to store considered phases

  note (01/05/2023):
  ------------------

  n (array dim) was originally equal to gv.len_ox + 1; However this has been changed following update
  on liq model (Green et al., 2023) and the use of 4 order variables summing up the total
  number of compositional variables to 14 for a chemical system of 11 oxides
  This means that the max number of the endmembers can be >= n-oxide + 1
*/
csd_phase_set CP_INIT_function(csd_phase_set cp, global_variable gv){
	
	int n 			= gv.max_ss_size_cp;
	
	/* initialize fractions flags and cycle arrays with zeros */
	cp.ss_flags 	= malloc (gv.n_flags  * sizeof(int));
	cp.name   		= malloc (20 * sizeof(char)		    );
	cp.p_em   		= malloc (n  * sizeof(double) 		);
	cp.xi_em  		= malloc (n  * sizeof(double) 		);
	cp.dguess 		= malloc (n  * sizeof(double) 		);
	cp.xeos   		= malloc (n  * sizeof(double) 		);
	cp.xeos_0   	= malloc (n  * sizeof(double) 		);
	cp.xeos_1   	= malloc (n  * sizeof(double) 		);
	cp.xeos_r   	= malloc (n  * sizeof(double) 		);
	cp.delta_mu		= malloc (n  * sizeof(double) 		);
	cp.dfx   		= malloc (n  * sizeof(double) 		);
	cp.mu    		= malloc (n  * sizeof(double) 		);
	cp.gbase    	= malloc (n  * sizeof(double) 		);
	cp.ss_comp		= malloc (n  * sizeof(double) 		);
	cp.sf			= malloc ((n*2)  * sizeof(double) 	);
	
	cp.phase_density  		= 0.0;
	cp.phase_cp				= 0.0;
	cp.phase_expansivity	= 0.0;
	cp.phase_entropy		= 0.0;
	cp.phase_enthalpy		= 0.0;
		
	return cp;
}

/** 
  allocate memory to store considered phases
*/
stb_system SP_INIT_function(stb_system sp, global_variable gv){

	sp.MAGEMin_ver   		= malloc(50  		* sizeof(char)				);
	sp.oxides 	     		= malloc(gv.len_ox  * sizeof(char*)				);
	
	for (int i = 0; i < gv.len_ox; i++){
		sp.oxides[i] 		= malloc(20 * sizeof(char));	
	}
	sp.bulk 				= malloc(gv.len_ox  * sizeof(double)			);	
	sp.gamma 				= malloc(gv.len_ox  * sizeof(double)			);	
	sp.bulk_S 				= malloc(gv.len_ox 	* sizeof(double)			);	
	sp.bulk_M 				= malloc(gv.len_ox 	* sizeof(double)			);	
	sp.bulk_F 				= malloc(gv.len_ox 	* sizeof(double)			);	

	sp.bulk_wt 				= malloc(gv.len_ox  * sizeof(double)			);	
	sp.bulk_S_wt 			= malloc(gv.len_ox 	* sizeof(double)			);	
	sp.bulk_M_wt 			= malloc(gv.len_ox 	* sizeof(double)			);	
	sp.bulk_F_wt 			= malloc(gv.len_ox 	* sizeof(double)			);	
	sp.ph 	     			= malloc(gv.len_ox  * sizeof(char*)				);
	sp.ph_frac 	     		= malloc(gv.len_ox  * sizeof(double)			);
	sp.ph_frac_wt     		= malloc(gv.len_ox  * sizeof(double)			);
	sp.ph_frac_vol     		= malloc(gv.len_ox  * sizeof(double)			);
	for (int i = 0; i < gv.len_ox; i++){
		sp.ph[i] 			= malloc(20 * sizeof(char));	
	}
	sp.ph_type 				= malloc(gv.len_ox 	* sizeof(int)				);	
	sp.ph_id 				= malloc(gv.len_ox 	* sizeof(int)				);	
	sp.PP 		 			= malloc(gv.len_ox  * sizeof(stb_PP_phase)		); 
	sp.SS 		 			= malloc(gv.len_ox  * sizeof(stb_SS_phase)		); 
	sp.mSS 		 			= malloc(gv.max_n_mSS  * sizeof(mstb_SS_phase)	); 

	for (int n = 0; n< gv.len_ox; n++){
		sp.PP[n].Comp 			= malloc(gv.len_ox 	* sizeof(double)		);
		sp.SS[n].Comp 			= malloc(gv.len_ox 	* sizeof(double)		);
		sp.PP[n].Comp_wt 		= malloc(gv.len_ox 	* sizeof(double)		);
		sp.SS[n].Comp_wt 		= malloc(gv.len_ox 	* sizeof(double)		);
		sp.SS[n].compVariables	= malloc(gv.len_ox*3   * sizeof(double)	    );
		sp.SS[n].emFrac			= malloc((gv.len_ox*3) * sizeof(double)		);
		sp.SS[n].emFrac_wt		= malloc((gv.len_ox*3) * sizeof(double)		);
		sp.SS[n].emChemPot		= malloc((gv.len_ox*3) * sizeof(double)		);
		sp.SS[n].compVariablesNames	= malloc(gv.len_ox*3 * sizeof(char*)	);
		sp.SS[n].emNames 	    = malloc((gv.len_ox*3) * sizeof(char*)		);
		sp.SS[n].emComp 	    = malloc((gv.len_ox*3) * sizeof(double*)	);
		sp.SS[n].emComp_wt 	    = malloc((gv.len_ox*3) * sizeof(double*)	);

		for (int i = 0; i < gv.len_ox*3; i++){
            sp.SS[n].compVariablesNames[i]		= malloc(20 * sizeof(char)	);
			sp.SS[n].emNames[i]		= malloc(20 * sizeof(char)				);
			sp.SS[n].emComp[i]		= malloc(gv.len_ox * sizeof(double)		);		
			sp.SS[n].emComp_wt[i]	= malloc(gv.len_ox * sizeof(double)		);		
		}
	}
    
    /** allocate memory for metastable phases len_ox * 2 to be safe?        */
	for (int n = 0; n< gv.max_n_mSS; n++){
        sp.mSS[n].ph_name	    = malloc(20 * sizeof(char)	                );
        sp.mSS[n].ph_type	    = malloc(20 * sizeof(char)	                );
        sp.mSS[n].info	        = malloc(20 * sizeof(char)	                );
        sp.mSS[n].comp_Ppc 	    = malloc((gv.len_ox) 	* sizeof(double)	);  
        sp.mSS[n].p_Ppc 	    = malloc((gv.len_ox*2) 	* sizeof(double)	);  
        sp.mSS[n].mu_Ppc 	    = malloc((gv.len_ox*2) 	* sizeof(double)	);  
        sp.mSS[n].xeos_Ppc 	    = malloc((gv.len_ox*2) 	* sizeof(double)	);  
	}

	return sp;
}





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


SS_ref G_SS_aq17_init_function(SS_ref SS_ref_db, int EM_database, global_variable gv){
    
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
SS_ref G_SS_mb_liq_init_function(SS_ref SS_ref_db, int EM_database, global_variable gv){
    
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
SS_ref G_SS_mb_hb_init_function(SS_ref SS_ref_db, int EM_database, global_variable gv){
    
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
SS_ref G_SS_mb_aug_init_function(SS_ref SS_ref_db, int EM_database, global_variable gv){
    
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
SS_ref G_SS_mb_dio_init_function(SS_ref SS_ref_db, int EM_database, global_variable gv){
    
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
SS_ref G_SS_mb_opx_init_function(SS_ref SS_ref_db, int EM_database, global_variable gv){
    
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
SS_ref G_SS_mb_g_init_function(SS_ref SS_ref_db, int EM_database, global_variable gv){
    
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
SS_ref G_SS_mb_ol_init_function(SS_ref SS_ref_db, int EM_database, global_variable gv){
    
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
SS_ref G_SS_mb_fsp_init_function(SS_ref SS_ref_db, int EM_database, global_variable gv){
    
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
SS_ref G_SS_mb_abc_init_function(SS_ref SS_ref_db, int EM_database, global_variable gv){
    
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
SS_ref G_SS_mb_k4tr_init_function(SS_ref SS_ref_db, int EM_database, global_variable gv){
    
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
    allocate memory for sp
*/
SS_ref G_SS_mb_sp_init_function(SS_ref SS_ref_db, int EM_database, global_variable gv){
    
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
SS_ref G_SS_mb_ilm_init_function(SS_ref SS_ref_db, int EM_database, global_variable gv){
    
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
SS_ref G_SS_mb_ilmm_init_function(SS_ref SS_ref_db, int EM_database, global_variable gv){
    
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
SS_ref G_SS_mb_ep_init_function(SS_ref SS_ref_db, int EM_database, global_variable gv){
    
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
SS_ref G_SS_mb_bi_init_function(SS_ref SS_ref_db, int EM_database, global_variable gv){
    
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
SS_ref G_SS_mb_mu_init_function(SS_ref SS_ref_db, int EM_database, global_variable gv){
    
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
SS_ref G_SS_mb_chl_init_function(SS_ref SS_ref_db, int EM_database, global_variable gv){
    
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

SS_ref G_SS_mp_liq_init_function(SS_ref SS_ref_db, int EM_database, global_variable gv){
    
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
SS_ref G_SS_mp_fsp_init_function(SS_ref SS_ref_db, int EM_database, global_variable gv){
    
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
SS_ref G_SS_mp_bi_init_function(SS_ref SS_ref_db, int EM_database, global_variable gv){
    
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
SS_ref G_SS_mp_g_init_function(SS_ref SS_ref_db, int EM_database, global_variable gv){
    
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
SS_ref G_SS_mp_ep_init_function(SS_ref SS_ref_db, int EM_database, global_variable gv){
    
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
SS_ref G_SS_mp_ma_init_function(SS_ref SS_ref_db, int EM_database, global_variable gv){
    
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
SS_ref G_SS_mp_mu_init_function(SS_ref SS_ref_db, int EM_database, global_variable gv){
    
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
SS_ref G_SS_mp_opx_init_function(SS_ref SS_ref_db, int EM_database, global_variable gv){
    
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
SS_ref G_SS_mp_sa_init_function(SS_ref SS_ref_db, int EM_database, global_variable gv){
    
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
SS_ref G_SS_mp_cd_init_function(SS_ref SS_ref_db, int EM_database, global_variable gv){
    
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
SS_ref G_SS_mp_st_init_function(SS_ref SS_ref_db, int EM_database, global_variable gv){
    
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
SS_ref G_SS_mp_chl_init_function(SS_ref SS_ref_db, int EM_database, global_variable gv){
    
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
SS_ref G_SS_mp_ctd_init_function(SS_ref SS_ref_db, int EM_database, global_variable gv){
    
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
SS_ref G_SS_mp_sp_init_function(SS_ref SS_ref_db, int EM_database, global_variable gv){
    
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
SS_ref G_SS_mp_ilm_init_function(SS_ref SS_ref_db, int EM_database, global_variable gv){
    
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
SS_ref G_SS_mp_ilmm_init_function(SS_ref SS_ref_db, int EM_database, global_variable gv){
    
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
SS_ref G_SS_mp_mt_init_function(SS_ref SS_ref_db, int EM_database, global_variable gv){
    
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
SS_ref G_SS_ig_fper_init_function(SS_ref SS_ref_db, int EM_database, global_variable gv){
    
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
SS_ref G_SS_ig_bi_init_function(SS_ref SS_ref_db, int EM_database, global_variable gv){		

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
SS_ref G_SS_ig_cpx_init_function(SS_ref SS_ref_db, int EM_database, global_variable gv){

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
SS_ref G_SS_ig_cd_init_function(SS_ref SS_ref_db, int EM_database, global_variable gv){

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
SS_ref G_SS_ig_ep_init_function(SS_ref SS_ref_db, int EM_database, global_variable gv){

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
SS_ref G_SS_ig_fl_init_function(SS_ref SS_ref_db, int EM_database, global_variable gv){
    
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
SS_ref G_SS_ig_g_init_function(SS_ref SS_ref_db, int EM_database, global_variable gv){

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
SS_ref G_SS_ig_hb_init_function(SS_ref SS_ref_db, int EM_database, global_variable gv){

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
SS_ref G_SS_ig_ilm_init_function(SS_ref SS_ref_db, int EM_database, global_variable gv){
    
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
SS_ref G_SS_ig_liq_init_function(SS_ref SS_ref_db, int EM_database, global_variable gv){
    
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
SS_ref G_SS_ig_mu_init_function(SS_ref SS_ref_db, int EM_database, global_variable gv){		

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
SS_ref G_SS_ig_ol_init_function(SS_ref SS_ref_db, int EM_database, global_variable gv){		

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
SS_ref G_SS_ig_opx_init_function(SS_ref SS_ref_db, int EM_database, global_variable gv){		

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
SS_ref G_SS_ig_fsp_init_function(SS_ref SS_ref_db, int EM_database, global_variable gv){

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
SS_ref G_SS_ig_spn_init_function(SS_ref SS_ref_db, int EM_database, global_variable gv){		

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
/*********************Evan&Frost DATABASE (Evans&Frost , 2021)*************************/
/**************************************************************************************/
/**************************************************************************************/

/**
    allocate memory for fluid
*/
SS_ref G_SS_um_fluid_init_function(SS_ref SS_ref_db, int EM_database, global_variable gv){
    
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
SS_ref G_SS_um_ol_init_function(SS_ref SS_ref_db, int EM_database, global_variable gv){
    
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
SS_ref G_SS_um_br_init_function(SS_ref SS_ref_db, int EM_database, global_variable gv){
    
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
SS_ref G_SS_um_ch_init_function(SS_ref SS_ref_db, int EM_database, global_variable gv){
    
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
SS_ref G_SS_um_atg_init_function(SS_ref SS_ref_db, int EM_database, global_variable gv){
    
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
SS_ref G_SS_um_g_init_function(SS_ref SS_ref_db, int EM_database, global_variable gv){
    
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
SS_ref G_SS_um_ta_init_function(SS_ref SS_ref_db, int EM_database, global_variable gv){
    
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
SS_ref G_SS_um_chl_init_function(SS_ref SS_ref_db, int EM_database, global_variable gv){
    
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
SS_ref G_SS_um_anth_init_function(SS_ref SS_ref_db, int EM_database, global_variable gv){
    
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
SS_ref G_SS_um_spi_init_function(SS_ref SS_ref_db, int EM_database, global_variable gv){
    
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
SS_ref G_SS_um_opx_init_function(SS_ref SS_ref_db, int EM_database, global_variable gv){
    
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
SS_ref G_SS_um_po_init_function(SS_ref SS_ref_db, int EM_database, global_variable gv){
    
    SS_ref_db.is_liq    = 0;
    SS_ref_db.symmetry  = 1;
    SS_ref_db.n_sf      = 2;
    SS_ref_db.n_em      = 2;
    SS_ref_db.n_w       = 1;
    SS_ref_db.n_xeos    = 1;
    
    return SS_ref_db;
}

/**
  attributes the right solution phase to the solution phase array
*/
SS_ref G_SS_init_EM_function(		int			 		 ph_id,
									SS_ref 		 		 SS_ref_db, 
									int 		 		 EM_database, 
									char 				*name, 
									global_variable 	 gv					){
	if (EM_database == 0) {	 //"bi","cd","cpx","ep","fl","g","hb","ilm","liq","mu","ol","opx","fsp","spn"	
		if      (strcmp( name, "liq")  == 0 ){
			SS_ref_db  = G_SS_mp_liq_init_function(SS_ref_db, EM_database, gv); 	}
		else if (strcmp( name, "fsp")  == 0){
			SS_ref_db  = G_SS_mp_fsp_init_function(SS_ref_db, EM_database, gv); 	}
		else if (strcmp( name, "bi")  == 0){
			SS_ref_db  = G_SS_mp_bi_init_function(SS_ref_db, EM_database, gv); 		}
		else if (strcmp( name, "g")  == 0){
			SS_ref_db  = G_SS_mp_g_init_function(SS_ref_db, EM_database, gv); 		}
		else if (strcmp( name, "ep")  == 0){
			SS_ref_db  = G_SS_mp_ep_init_function(SS_ref_db, EM_database, gv); 		}
		else if (strcmp( name, "ma")  == 0){
			SS_ref_db  = G_SS_mp_ma_init_function(SS_ref_db, EM_database, gv); 		}
		else if (strcmp( name, "mu")  == 0){
			SS_ref_db  = G_SS_mp_mu_init_function(SS_ref_db, EM_database, gv); 		}
		else if (strcmp( name, "opx")  == 0){
			SS_ref_db  = G_SS_mp_opx_init_function(SS_ref_db, EM_database, gv); 	}
		else if (strcmp( name, "sa")  == 0){
			SS_ref_db  = G_SS_mp_sa_init_function(SS_ref_db, EM_database, gv); 		}
		else if (strcmp( name, "cd")  == 0){
			SS_ref_db  = G_SS_mp_cd_init_function(SS_ref_db, EM_database, gv); 		}
		else if (strcmp( name, "st")  == 0){
			SS_ref_db  = G_SS_mp_st_init_function(SS_ref_db, EM_database, gv); 		}
		else if (strcmp( name, "chl")  == 0){
			SS_ref_db  = G_SS_mp_chl_init_function(SS_ref_db, EM_database, gv); 	}
		else if (strcmp( name, "ctd")  == 0){
			SS_ref_db  = G_SS_mp_ctd_init_function(SS_ref_db, EM_database, gv); 	}
		else if (strcmp( name, "sp")  == 0){
			SS_ref_db  = G_SS_mp_sp_init_function(SS_ref_db, EM_database, gv); 		}
		else if (strcmp( name, "ilm")  == 0){
			SS_ref_db  = G_SS_mp_ilm_init_function(SS_ref_db, EM_database, gv); 	}
		else if (strcmp( name, "ilmm")  == 0){
			SS_ref_db  = G_SS_mp_ilmm_init_function(SS_ref_db, EM_database, gv); 	}
		else if (strcmp( name, "mt")  == 0){
			SS_ref_db  = G_SS_mp_mt_init_function(SS_ref_db, EM_database, gv); 		}
		else if (strcmp( name, "aq17")  == 0){
			SS_ref_db  = G_SS_aq17_init_function(SS_ref_db, EM_database, gv); 		}

		else{
			printf("\nsolid solution '%s' is not in the database, cannot be initiated\n", name);	
		}	
	}	
    else if (EM_database == 1){
      if (strcmp( name, "liq") == 0 ){
         SS_ref_db  = G_SS_mb_liq_init_function(SS_ref_db, EM_database, gv); }
      else if (strcmp( name, "hb") == 0 ){
         SS_ref_db  = G_SS_mb_hb_init_function(SS_ref_db, EM_database, gv); }
      else if (strcmp( name, "aug") == 0 ){
         SS_ref_db  = G_SS_mb_aug_init_function(SS_ref_db, EM_database, gv); }
      else if (strcmp( name, "dio") == 0 ){
         SS_ref_db  = G_SS_mb_dio_init_function(SS_ref_db, EM_database, gv); }
      else if (strcmp( name, "opx") == 0 ){
         SS_ref_db  = G_SS_mb_opx_init_function(SS_ref_db, EM_database, gv); }
      else if (strcmp( name, "g") == 0 ){
         SS_ref_db  = G_SS_mb_g_init_function(SS_ref_db, EM_database, gv); }
      else if (strcmp( name, "ol") == 0 ){
         SS_ref_db  = G_SS_mb_ol_init_function(SS_ref_db, EM_database, gv); }
      else if (strcmp( name, "fsp") == 0 ){
         SS_ref_db  = G_SS_mb_fsp_init_function(SS_ref_db, EM_database, gv); }
      else if (strcmp( name, "abc") == 0 ){
         SS_ref_db  = G_SS_mb_abc_init_function(SS_ref_db, EM_database, gv); }
      else if (strcmp( name, "k4tr") == 0 ){
         SS_ref_db  = G_SS_mb_k4tr_init_function(SS_ref_db, EM_database, gv); }
      else if (strcmp( name, "sp") == 0 ){
         SS_ref_db  = G_SS_mb_sp_init_function(SS_ref_db, EM_database, gv); }
      else if (strcmp( name, "ilm") == 0 ){
         SS_ref_db  = G_SS_mb_ilm_init_function(SS_ref_db, EM_database, gv); }
      else if (strcmp( name, "ilmm") == 0 ){
         SS_ref_db  = G_SS_mb_ilmm_init_function(SS_ref_db, EM_database, gv); }
      else if (strcmp( name, "ep") == 0 ){
         SS_ref_db  = G_SS_mb_ep_init_function(SS_ref_db, EM_database, gv); }
      else if (strcmp( name, "bi") == 0 ){
         SS_ref_db  = G_SS_mb_bi_init_function(SS_ref_db, EM_database, gv); }
      else if (strcmp( name, "mu") == 0 ){
         SS_ref_db  = G_SS_mb_mu_init_function(SS_ref_db, EM_database, gv); }
      else if (strcmp( name, "chl") == 0 ){
         SS_ref_db  = G_SS_mb_chl_init_function(SS_ref_db, EM_database, gv); }
      else{
         printf("\nsolid solution '%s' is not in the database\n",name);
      }
   }						  
	else if (EM_database == 2) {	 //"bi","cd","cpx","ep","fl","g","hb","ilm","liq","mu","ol","opx","fsp","spn"	
		if      (strcmp( name, "bi")  == 0 ){
			SS_ref_db  = G_SS_ig_bi_init_function(SS_ref_db, EM_database, gv); 		}
		else if (strcmp( name, "fper") == 0 ){
			SS_ref_db  = G_SS_ig_fper_init_function(SS_ref_db, EM_database, gv); }
		else if (strcmp( name, "cd")  == 0){
			SS_ref_db  = G_SS_ig_cd_init_function(SS_ref_db, EM_database, gv); 		}
		else if (strcmp( name, "cpx") == 0){
			SS_ref_db  = G_SS_ig_cpx_init_function(SS_ref_db, EM_database, gv); 	}	
		else if (strcmp( name, "ep")  == 0){
			SS_ref_db  = G_SS_ig_ep_init_function(SS_ref_db, EM_database, gv); 		}
		else if (strcmp( name, "fl")  == 0){
			SS_ref_db  = G_SS_ig_fl_init_function(SS_ref_db, EM_database, gv); 		}		
		else if (strcmp( name, "g")   == 0){
			SS_ref_db  = G_SS_ig_g_init_function(SS_ref_db, EM_database, gv); 		}
		else if (strcmp( name, "hb")  == 0){
			SS_ref_db  = G_SS_ig_hb_init_function(SS_ref_db, EM_database, gv); 		}	
		else if (strcmp( name, "ilm") == 0){
			SS_ref_db  = G_SS_ig_ilm_init_function(SS_ref_db, EM_database, gv); 	}
		else if (strcmp( name, "liq") == 0){
			SS_ref_db  = G_SS_ig_liq_init_function(SS_ref_db, EM_database, gv); 	}
		else if (strcmp( name, "mu")  == 0){
			SS_ref_db  = G_SS_ig_mu_init_function(SS_ref_db, EM_database, gv); 		}	
		else if (strcmp( name, "ol")  == 0){
			SS_ref_db  = G_SS_ig_ol_init_function(SS_ref_db, EM_database, gv); 		}
		else if (strcmp( name, "opx") == 0){
			SS_ref_db  = G_SS_ig_opx_init_function(SS_ref_db, EM_database, gv); 	}
		else if (strcmp( name, "fsp") == 0){
			SS_ref_db  = G_SS_ig_fsp_init_function(SS_ref_db, EM_database, gv); 	}	
		else if (strcmp( name, "spn") == 0){
			SS_ref_db  = G_SS_ig_spn_init_function(SS_ref_db, EM_database, gv); 	}
		else{
			printf("\nsolid solution '%s' is not in the database, cannot be initiated\n", name);	
		}	
	}
	else if (EM_database == 4) {
		if      (strcmp( name, "fl")  == 0 ){
			SS_ref_db  = G_SS_um_fluid_init_function(SS_ref_db, EM_database, gv); 	}
		else if (strcmp( name, "ol")  == 0){
			SS_ref_db  = G_SS_um_ol_init_function(SS_ref_db, EM_database, gv); 		}
		else if (strcmp( name, "br") == 0){
			SS_ref_db  = G_SS_um_br_init_function(SS_ref_db, EM_database, gv); 	    }	
		else if (strcmp( name, "ch")  == 0){
			SS_ref_db  = G_SS_um_ch_init_function(SS_ref_db, EM_database, gv); 		}
		else if (strcmp( name, "atg")  == 0){
			SS_ref_db  = G_SS_um_atg_init_function(SS_ref_db, EM_database, gv); 	}		
		else if (strcmp( name, "g")   == 0){
			SS_ref_db  = G_SS_um_g_init_function(SS_ref_db, EM_database, gv); 		}
		else if (strcmp( name, "ta")  == 0){
			SS_ref_db  = G_SS_um_ta_init_function(SS_ref_db, EM_database, gv); 		}	
		else if (strcmp( name, "chl") == 0){
			SS_ref_db  = G_SS_um_chl_init_function(SS_ref_db, EM_database, gv); 	}
		else if (strcmp( name, "anth") == 0){
			SS_ref_db  = G_SS_um_anth_init_function(SS_ref_db, EM_database, gv); 	}
		else if (strcmp( name, "spi")  == 0){
			SS_ref_db  = G_SS_um_spi_init_function(SS_ref_db, EM_database, gv); 	}	
		else if (strcmp( name, "opx")  == 0){
			SS_ref_db  = G_SS_um_opx_init_function(SS_ref_db, EM_database, gv); 	}
		else if (strcmp( name, "po") == 0){
			SS_ref_db  = G_SS_um_po_init_function(SS_ref_db, EM_database, gv); 	    }
		else{
			printf("\nsolid solution '%s' is not in the database, cannot be initiated\n", name);	
		}	
	}

	/**
		Allocate memory for solution phase models and pseudocompound storage (memory is initialized in the reset function)
	*/
	int n_em   = SS_ref_db.n_em;
	int n_xeos = SS_ref_db.n_xeos;
	int n_sf   = SS_ref_db.n_sf;
	int sym    = SS_ref_db.symmetry;

    SS_ref_db.orderVar       = 0;
	
	SS_ref_db.EM_list 		 = malloc ((n_em) * sizeof (char*)	);
	for (int i = 0; i < n_em; i++){ 
		SS_ref_db.EM_list[i] = malloc(20 * sizeof(char)		);		
	}
	SS_ref_db.CV_list 		 = malloc ((n_xeos) * sizeof (char*)	);
	for (int i = 0; i < n_xeos; i++){ 
		SS_ref_db.CV_list[i] = malloc(20 * sizeof(char)		);		
	}
	if (sym == 0){
		SS_ref_db.W   		= malloc (SS_ref_db.n_w * sizeof (double) ); 
		SS_ref_db.v   		= malloc (SS_ref_db.n_v * sizeof (double) ); 
	}
	else if (sym == 1){
		SS_ref_db.W   		= malloc (SS_ref_db.n_w * sizeof (double) ); 
	}
	
	/* initialize fractions flags and cycle arrays with zeros */
	SS_ref_db.ss_flags      = malloc (gv.n_flags  * sizeof(int));
    SS_ref_db.solvus_id     = malloc ((gv.len_ss*4)   * sizeof (int)  	);

	/* dynamic memory allocation of data to send to NLopt */
	SS_ref_db.bounds 	= malloc (n_xeos * sizeof (double*)  ); 
	for (int i = 0; i < n_xeos; i++){
		SS_ref_db.bounds[i] = malloc (2 * sizeof (double) );
	}
	
	SS_ref_db.eye 			= malloc (n_em 			* sizeof (double*)); 
	SS_ref_db.dp_dx 		= malloc (n_em 			* sizeof (double*)); 
	SS_ref_db.Comp 			= malloc (n_em 			* sizeof (double*)); 
	for (int i = 0; i < n_em; i++){
		SS_ref_db.eye[i] 	= malloc (n_em	 		* sizeof (double) );
		SS_ref_db.Comp[i] 	= malloc (gv.len_ox 	* sizeof (double) );
		SS_ref_db.dp_dx[i] 	= malloc (n_xeos 		* sizeof (double) );
	}
	
	SS_ref_db.gbase   		= malloc (n_em   	 	* sizeof (double) ); 
	SS_ref_db.gb_lvl  		= malloc (n_em   	 	* sizeof (double) ); 
	SS_ref_db.z_em    		= malloc (n_em   	 	* sizeof (double) ); 
	SS_ref_db.d_em    		= malloc (n_em   	 	* sizeof (double) );
	SS_ref_db.density 		= malloc (n_em   	 	* sizeof (double) ); 
	SS_ref_db.dguess 		= malloc (n_xeos 		* sizeof (double) );
	SS_ref_db.iguess  		= malloc (n_xeos   	  	* sizeof (double) );
	SS_ref_db.mguess  		= malloc (n_xeos   	  	* sizeof (double) );
	SS_ref_db.idOrderVar  	= malloc (n_xeos   	  	* sizeof (double) );
    
	SS_ref_db.p       		= malloc (n_em       	* sizeof (double) ); 
	SS_ref_db.ElShearMod    = malloc (n_em       	* sizeof (double) ); 
	SS_ref_db.ape      		= malloc (n_em       	* sizeof (double) ); 
	SS_ref_db.mat_phi 		= malloc (n_em       	* sizeof (double) ); 
	SS_ref_db.mu_Gex  		= malloc (n_em       	* sizeof (double) ); 
	SS_ref_db.sf      		= malloc (n_sf       	* sizeof (double) ); 
	SS_ref_db.mu      		= malloc (n_em       	* sizeof (double) ); 
	SS_ref_db.dfx    		= malloc (n_xeos     	* sizeof (double) ); 
	SS_ref_db.ss_comp		= malloc (gv.len_ox  	* sizeof (double) ); 
	SS_ref_db.ElEntropy		= malloc (gv.len_ox  	* sizeof (double) );     
	SS_ref_db.xi_em   		= malloc (n_em   	 	* sizeof (double) ); 
	SS_ref_db.xeos    		= malloc (n_xeos     	* sizeof (double) ); 	
	SS_ref_db.xeos_sf_ok 	= malloc ((n_xeos) 		* sizeof (double) );

	/* memory allocation to store all gbase */
	SS_ref_db.mu_array = malloc ((gv.n_Diff) * sizeof (double*) ); 
	for (int i = 0; i < (gv.n_Diff); i++){
		SS_ref_db.mu_array[i] = malloc (n_em * sizeof (double) );
	}	
	
	/* dynamic memory allocation of data to send to NLopt */
	SS_ref_db.bounds_ref = malloc ((n_xeos) * sizeof (double*) ); 
	for (int i = 0; i < (n_xeos); i++){
		SS_ref_db.bounds_ref[i] = malloc (2 * sizeof (double) );
	}
	
	/* dynamic memory allocation of data to send to NLopt */
	SS_ref_db.ub   		= malloc ((n_xeos) * sizeof (double) ); 
	SS_ref_db.lb   		= malloc ((n_xeos) * sizeof (double) ); 
	SS_ref_db.tol_sf   	= malloc ((n_sf) * sizeof (double) ); 
	for (int j = 0; j < n_sf; j++){
		SS_ref_db.tol_sf[j] = gv.ineq_res;
	}
	
	/**
		Allocate memory for levelling pseudocompounds 
	*/
	// SS_ref_db.n_pc   	= gv.n_pc;
	SS_ref_db.n_pc   	= gv.n_SS_PC[ph_id];
    SS_ref_db.tot_pc   	= malloc (1 * sizeof (double) ); 
    SS_ref_db.id_pc   	= malloc (1 * sizeof (double) ); 
	SS_ref_db.G_pc   	= malloc ((SS_ref_db.n_pc) * sizeof (double) ); 
	SS_ref_db.DF_pc 	= malloc ((SS_ref_db.n_pc) * sizeof (double) ); 
	SS_ref_db.factor_pc = malloc ((SS_ref_db.n_pc) * sizeof (double) ); 
	SS_ref_db.info  	= malloc ((SS_ref_db.n_pc) * sizeof (int) 	 ); 
	SS_ref_db.p_pc 		= malloc ((SS_ref_db.n_pc) * sizeof (double*)); 
	// SS_ref_db.mu_pc 	= malloc ((SS_ref_db.n_pc) * sizeof (double*)); 
	
	for (int i = 0; i < (SS_ref_db.n_pc); i++){
		SS_ref_db.p_pc[i] 	 = malloc ((n_em) * sizeof (double) 	);
		// SS_ref_db.mu_pc[i] 	 = malloc ((n_em) * sizeof (double) 	);
	}
	SS_ref_db.comp_pc = malloc ((SS_ref_db.n_pc) * sizeof (double*) ); 
	for (int i = 0; i < (SS_ref_db.n_pc); i++){
		SS_ref_db.comp_pc[i] = malloc (gv.len_ox * sizeof (double) 	);
	}
	SS_ref_db.xeos_pc = malloc ((SS_ref_db.n_pc) * sizeof (double*) ); 
	for (int i = 0; i < (SS_ref_db.n_pc); i++){
		SS_ref_db.xeos_pc[i] = malloc ((n_xeos)  * sizeof (double) 	);
	}	

	/**
		Allocate memory for PGE pseudocompounds 
	*/
	SS_ref_db.n_Ppc   	= gv.n_Ppc;								/** maximum number of pseudocompounds to store */
	SS_ref_db.G_Ppc   	= malloc ((SS_ref_db.n_Ppc) * sizeof (double) ); 
	SS_ref_db.DF_Ppc 	= malloc ((SS_ref_db.n_Ppc) * sizeof (double) ); 
	SS_ref_db.info_Ppc 	= malloc ((SS_ref_db.n_Ppc) * sizeof (int) 	 ); 
	SS_ref_db.p_Ppc 	= malloc ((SS_ref_db.n_Ppc) * sizeof (double*)); 
	SS_ref_db.mu_Ppc 	= malloc ((SS_ref_db.n_Ppc) * sizeof (double*)); 
	
	for (int i = 0; i < (SS_ref_db.n_Ppc); i++){
		SS_ref_db.p_Ppc[i] 	 = malloc ((n_em) * sizeof (double) 		);
		SS_ref_db.mu_Ppc[i]  = malloc ((n_em) * sizeof (double) 		);
	}
	SS_ref_db.comp_Ppc = malloc ((SS_ref_db.n_Ppc) * sizeof (double*) 	); 
	for (int i = 0; i < (SS_ref_db.n_Ppc); i++){
		SS_ref_db.comp_Ppc[i] = malloc (gv.len_ox * sizeof (double) 	);
	}
	SS_ref_db.xeos_Ppc = malloc ((SS_ref_db.n_Ppc) * sizeof (double*) 	); 
	for (int i = 0; i < (SS_ref_db.n_Ppc); i++){
		SS_ref_db.xeos_Ppc[i] = malloc ((n_xeos)  * sizeof (double) 	);
	}	

	/* initiliazes eye matrix as there is no need to redo it afterward */
	for (int i = 0; i < n_em; i++){
		for (int j = 0; j < n_em; j++){
			SS_ref_db.eye[i][j] = 0.0;
		}	
	}
	
	/* initialize eye matrix */
	for (int j = 0; j < n_em; j++){
		SS_ref_db.eye[j][j]= 1.0;
	}

	for (int j = 0; j < n_xeos; j++){
		SS_ref_db.idOrderVar[j]= 1.0;
	}
    


	return SS_ref_db;
};

