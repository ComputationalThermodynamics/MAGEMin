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
#ifndef __SB_GSS_FUNCTION_H_
#define __SB_GSS_FUNCTION_H_

#include "../MAGEMin.h"

SS_ref G_SS_sb11_EM_function(	global_variable  gv, 
							    SS_ref 			 SS_ref_db,
                                int 			 EM_dataset,
                                bulk_info 		 z_b,
                                char 			*name					);

SS_ref G_SS_sb21_EM_function(	global_variable  gv, 
                                SS_ref 			 SS_ref_db,
                                int 			 EM_dataset,
                                bulk_info 		 z_b,
                                char 			*name					);
                                
SS_ref G_SS_sb24_EM_function(	global_variable  gv, 
                                SS_ref 			 SS_ref_db,
                                int 			 EM_dataset,
                                bulk_info 		 z_b,
                                char 			*name					);						
#endif
