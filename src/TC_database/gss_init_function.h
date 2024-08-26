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
#ifndef __GSS_INIT_FUNCTION_H_
#define __GSS_INIT_FUNCTION_H_

#include "../MAGEMin.h"

typedef SS_ref (*SS_init_type) (	SS_ref 				 SS_ref_db,
									global_variable 	 gv				);
									
void TC_SS_init_mp(	                SS_init_type 		*SS_init,
									global_variable 	 gv				);				
void TC_SS_init_mb(	                SS_init_type 		*SS_init,
									global_variable 	 gv				);
void TC_SS_init_ig(	                SS_init_type 		*SS_init,
									global_variable 	 gv				);
void TC_SS_init_um(	                SS_init_type 		*SS_init,
									global_variable 	 gv				);
void TC_SS_init(	        	    SS_init_type 		*SS_init,
									global_variable 	 gv				);

SS_ref G_SS_init_EM_function(		SS_init_type   *SS_init,
									int				ph_id,
									SS_ref 			SS_ref_db,
									char 		   *name, 
									global_variable gv				);

csd_phase_set CP_INIT_function(		csd_phase_set 	cp, 
									global_variable gv				);
									
stb_system SP_INIT_function(		stb_system 		sp,
									global_variable gv				);

#endif
