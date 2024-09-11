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
#ifndef __INITIALIZE_H_
#define __INITIALIZE_H_

	#include <string.h>
	#include "uthash.h"
	#include "MAGEMin.h"
	#include "all_endmembers.h"
	#include "all_init_database.h"

	typedef SS_ref (*SS_init_type) (	SS_ref 				SS_ref_db,
										global_variable 	gv					);

	char** get_EM_DB_names_tc(			global_variable 	gv					);
	char** get_EM_DB_names_sb(			global_variable 	gv					);
	char** get_FS_DB_names(				global_variable 	gv					);

	global_variable global_variable_alloc( bulk_info  	    *z_b 				);

	csd_phase_set CP_INIT_function(		csd_phase_set 		cp, 
										global_variable 	gv					);
										
	stb_system SP_INIT_function(		stb_system 			sp,
										global_variable 	gv					);
								
	SS_ref G_SS_init_EM_function(		SS_init_type   	   *SS_init,
										int					ph_id,
										SS_ref 				SS_ref_db,
										char 		   	   *name, 
										global_variable 	gv					);

	global_variable reset_gv(			global_variable 	 gv,
										bulk_info 	 		 z_b,
										PP_ref 				*PP_ref_db,
										SS_ref 				*SS_ref_db			);

	void reset_sp(						global_variable 	 gv,
										stb_system  		*sp					);

	bulk_info reset_z_b_bulk(			global_variable 	 gv,
										bulk_info 	 		 z_b				);

	void reset_cp(						global_variable 	 gv,
										bulk_info 	 		z_b,
										csd_phase_set  		*cp					);

	void reset_SS(						global_variable 	 gv,
										bulk_info 	 		 z_b,
										SS_ref 				*SS_ref_db			);

	void init_simplex_A( 				simplex_data 		*splx_data,
										global_variable 	 gv					);

	void init_simplex_B_em(				simplex_data 		 *splx_data,
										global_variable 	 gv					);

	void reset_simplex_A( 				simplex_data 		*splx_data,
										bulk_info 	 		 z_b,
										global_variable 	 gv					);

	void reset_simplex_B_em(			simplex_data 		*splx_data,
										global_variable 	 gv					);


#endif
