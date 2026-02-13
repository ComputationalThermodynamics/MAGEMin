/*@ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 **
 **   Project      : MAGEMin
 **   License      : GNU GENERAL PUBLIC LICENSE Version 3, 29 June 2007
 **   Developers   : Nicolas Riel, Boris Kaus
 **   Contributors : Nickolas B. Moccetti, Dominguez, H., Assunção J., Green E., Berlie N., and Rummel L.
 **   Organization : Institute of Geosciences, Johannes-Gutenberg University, Mainz
 **   Contact      : nriel[at]uni-mainz.de, kaus[at]uni-mainz.de
 **
 ** ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ @*/
#ifndef __NLOPT_OPT_FUNCTION_H_
#define __NLOPT_OPT_FUNCTION_H_

#include "../MAGEMin.h"

typedef SS_ref (*NLopt_type) (		global_variable 	 gv,
									SS_ref 				 SS_ref_db		);

void TC_mp_NLopt_opt_init(	        NLopt_type 			*NLopt_opt,
									global_variable 	 gv				);
void TC_mb_NLopt_opt_init(	        NLopt_type 			*NLopt_opt,
									global_variable 	 gv				);
void TC_mb_ext_NLopt_opt_init(	    NLopt_type 			*NLopt_opt,
									global_variable 	 gv				);
void TC_ig_NLopt_opt_init(	        NLopt_type 			*NLopt_opt,
									global_variable 	 gv				);
void TC_igad_NLopt_opt_init(	    NLopt_type 			*NLopt_opt,
									global_variable 	 gv				);
void TC_um_NLopt_opt_init(	        NLopt_type 			*NLopt_opt,
									global_variable 	 gv				);
void TC_um_ext_NLopt_opt_init(	    NLopt_type 			*NLopt_opt,
									global_variable 	 gv				);
void TC_mtl_NLopt_opt_init(	    	NLopt_type 			*NLopt_opt,
									global_variable 	 gv				);
void TC_NLopt_opt_init(	        	NLopt_type 			*NLopt_opt,
									global_variable 	 gv				);
void TC_mpe_NLopt_opt_init(	        NLopt_type 			*NLopt_opt,
									global_variable 	 gv				);

SS_ref NLopt_opt_ig_spl_function(global_variable gv, SS_ref SS_ref_db);
SS_ref NLopt_opt_ig_cpx_function(global_variable gv, SS_ref SS_ref_db);
SS_ref NLopt_opt_ig_amp_function(global_variable gv, SS_ref SS_ref_db);

#endif

