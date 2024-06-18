/*@ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 **
 **   Project      : MAGEMin
 **   License      : GNU GENERAL PUBLIC LICENSE Version 3, 29 June 2007
 **   Developers   : Nicolas Riel, Boris Kaus
 **   Contributors : Dominguez, H., Green E., Berlie N., and Rummel L.
 **   Organization : Institute of Geosciences, Johannes-Gutenberg University, Mainz
 **   Contact      : nriel[at]uni-mainz.de, kaus[at]uni-mainz.de
 **
 ** ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ @*/
#ifndef __NLOPT_OPT_FUNCTION_H_
#define __NLOPT_OPT_FUNCTION_H_

#include "../MAGEMin.h"

typedef void (*sf_type) (		unsigned  		 m,
								double 			*result,
								unsigned 		 n,
								const double 	*x,
								double 			*grad,
								void 			*data		);

	
global_variable NLopt_global_opt_function(	bulk_info 			 z_b,
											global_variable 	 gv, 
											PP_ref 				*PP_ref_db,
											SS_ref 				*SS_ref_db,
	 										csd_phase_set  		*cp			);

typedef SS_ref (*NLopt_type) (		global_variable 	 gv,
									SS_ref 				 SS_ref_db		);

void TC_mp_NLopt_opt_init(	        NLopt_type 			*NLopt_opt,
									global_variable 	 gv				);
void TC_mb_NLopt_opt_init(	        NLopt_type 			*NLopt_opt,
									global_variable 	 gv				);
void TC_ig_NLopt_opt_init(	        NLopt_type 			*NLopt_opt,
									global_variable 	 gv				);
void TC_um_NLopt_opt_init(	        NLopt_type 			*NLopt_opt,
									global_variable 	 gv				);

void TC_NLopt_opt_init(	        	NLopt_type 			*NLopt_opt,
									global_variable 	 gv				);

#endif

