/*@ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 **
 **   Project      : MAGEMin
 **   License      : GNU GENERAL PUBLIC LICENSE Version 3, 29 June 2007
 **   Developers   : Nicolas Riel, Boris Kaus
 **   Contributors : Moccetti, N. B., Dominguez, H., Assunção J., Green E., Dolejš, D., Berlie N., and Rummel L.
 **   Organization : Institute of Geosciences, Johannes-Gutenberg University, Mainz
 **   Contact      : nriel[at]uni-mainz.de, kaus[at]uni-mainz.de
 **
 ** ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ @*/
#ifndef __SS_MIN_FUNCTION_H_
#define __SS_MIN_FUNCTION_H_

#include "MAGEMin.h"
#include "all_solution_phases.h"

SS_ref SS_UPDATE_function(				global_variable 	 gv,
										SS_ref 				 SS_ref_db,
										bulk_info 	 		 z_b,
										char    			*name			);

/** Convert endmember fractions (p, length SS_ref_db.n_em, in gv.SS_list[ph_id]'s
    own EM_list order) to compositional variables (xeos), for one solution
    phase. Writes into both SS_ref_db.iguess[] (what PC_convert_function
    actually reads) and SS_ref_db.xeos[] (mirrored for callers inspecting
    either field). research_group must be "tc", "gh" or "br" ("sb" has no
    p->xeos mapping in MAGEMin at all). */
SS_ref P2X_convert_function(			global_variable 	 gv,
										SS_ref 				 SS_ref_db,
										int 				 ph_id,
										double 				*p,
										int 				 n_p			);

/** Given a solution phase's compositional variables (xeos, e.g. from
    P2X_convert_function) and its endmembers' reference energies at the
    desired P,T (already computed via ComputeG0_point on z_b/DB.SS_ref_db),
    compute the phase's Gibbs energy (.df, absolute/non-Gamma-rotated) and
    bulk composition (.ss_comp), matching the same non_rot_hyperplane +
    PC_function + SS_UPDATE_function pipeline used internally during a
    normal minimization. Does not compute density/volume/modulus - those
    need the separate compute_density_volume_modulus() step. */
SS_ref PC_convert_function(			global_variable 	 gv,
										SS_ref 				 SS_ref_db,
										bulk_info 	 		 z_b,
										int 				 ph_id			);


csd_phase_set CP_UPDATE_function(		global_variable 	 gv,
										SS_ref 				 SS_ref_db,
										csd_phase_set  		 cp, 
										bulk_info 	 		 z_b			);		

global_variable split_cp(				global_variable 	 gv,
										SS_ref 			    *SS_ref_db,
										csd_phase_set  		*cp				);

void init_PGE_from_LP(					global_variable 	 gv,
										PC_type 			*PC_read,
										obj_type 			*SS_objective,
										bulk_info 	 		 z_b,
										SS_ref 			    *SS_ref_db,
										csd_phase_set  		*cp				);

void ss_min_PGE(						global_variable 	 gv,
										PC_type				*PC_read,
										obj_type 			*SS_objective,
										NLopt_type			*NLopt_opt,
										bulk_info 	 		 z_b,
										SS_ref 				*SS_ref_db,
										csd_phase_set  		*cp				);
										
void ss_min_LP(							global_variable 	 gv,
										PC_type				*PC_read,
										
										obj_type 			*SS_objective,
										NLopt_type			*NLopt_opt,
										bulk_info 	 		 z_b,
										SS_ref 				*SS_ref_db,
										csd_phase_set  		*cp				);

void copy_to_Ppc(						int 				 pc_check,
										int 				 add,
										int 				 ph_id,
										global_variable 	 gv,

										obj_type 			*SS_objective,
										SS_ref 			    *SS_ref_db		);

void copy_to_cp(						int 				 i, 
										int 				 ph_id,
										global_variable 	 gv,
										SS_ref 			    *SS_ref_db,
										csd_phase_set  		*cp				);

global_variable init_ss_db(				int 				 EM_database,
										bulk_info 	 		 z_b,
										global_variable 	 gv,
										SS_ref 				*SS_ref_db		);

global_variable init_ss_db_sb(			int 				 EM_database,
										bulk_info 	 		 z_b,
										global_variable 	 gv,
										SS_ref 				*SS_ref_db		);

global_variable init_ss_db_gh(			int 				 EM_database,
										bulk_info 	 		 z_b,
										global_variable 	 gv,
										SS_ref 				*SS_ref_db		);

global_variable init_ss_db_br(			int 				 EM_database,
										bulk_info 	 		 z_b,
										global_variable 	 gv,
										SS_ref 				*SS_ref_db		);

#endif
