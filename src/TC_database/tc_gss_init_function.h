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
#ifndef __TC_GSS_INIT_FUNCTION_H_
#define __TC_GSS_INIT_FUNCTION_H_

#include "../MAGEMin.h"
#include "../initialize.h"
								
void TC_SS_init_mp(	                SS_init_type 		*SS_init,
									global_variable 	 gv				);				
void TC_SS_init_mb(	                SS_init_type 		*SS_init,
									global_variable 	 gv				);
void TC_SS_init_ig(	                SS_init_type 		*SS_init,
									global_variable 	 gv				);
void TC_SS_init_igad(	            SS_init_type 		*SS_init,
									global_variable 	 gv				);
void TC_SS_init_um(	                SS_init_type 		*SS_init,
									global_variable 	 gv				);
void TC_SS_init_um_ext(	            SS_init_type 		*SS_init,
									global_variable 	 gv				);
void TC_SS_init_mp_ext(	            SS_init_type 		*SS_init,
									global_variable 	 gv				);
void TC_SS_init_mb_ext(	            SS_init_type 		*SS_init,
									global_variable 	 gv				);
void TC_SS_init_mtl(	            SS_init_type 		*SS_init,
									global_variable 	 gv				);
void TC_SS_init(	        	    SS_init_type 		*SS_init,
									global_variable 	 gv				);

#endif
