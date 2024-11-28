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
#ifndef __TC_GSS_FUNCTION_H_
#define __TC_GSS_FUNCTION_H_

#include "../MAGEMin.h"

SS_ref G_SS_mp_EM_function(		global_variable  gv, 
								SS_ref 			 SS_ref_db,
								int 			 EM_dataset,
								bulk_info 		 z_b,
								char 			*name					);
							
SS_ref G_SS_mb_EM_function(		global_variable  gv, 
								SS_ref 			 SS_ref_db,
								int 			 EM_dataset,
								bulk_info 		 z_b,
								char 			*name					);

SS_ref G_SS_ig_EM_function(		global_variable  gv, 
								SS_ref 			 SS_ref_db,
								int 			 EM_dataset,
								bulk_info 		 z_b,
								char 			*name					);
							
SS_ref G_SS_igad_EM_function(	global_variable  gv, 
								SS_ref 			 SS_ref_db,
								int 			 EM_dataset,
								bulk_info 		 z_b,
								char 			*name					);

SS_ref G_SS_um_EM_function(		global_variable  gv, 
								SS_ref 			 SS_ref_db,
								int 			 EM_dataset,
								bulk_info 		 z_b,
								char 			*name					);

SS_ref G_SS_um_ext_EM_function(	global_variable  gv, 
								SS_ref 			 SS_ref_db,
								int 			 EM_dataset,
								bulk_info 		 z_b,
								char 			*name					);

SS_ref G_SS_mtl_EM_function(	global_variable  gv, 
								SS_ref 			 SS_ref_db,
								int 			 EM_dataset,
								bulk_info 		 z_b,
								char 			*name					);

SS_ref G_SS_mpe_EM_function(	global_variable  gv, 
								SS_ref 			 SS_ref_db,
								int 			 EM_dataset,
								bulk_info 		 z_b,
								char 			*name					);
#endif
