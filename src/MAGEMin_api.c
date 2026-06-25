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
/*@
 **   Minimal C API to call MAGEMin as a library from external C/C++ code.
 @*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "toolkit.h"
#include "io_function.h"
#include "gem_function.h"

#include "simplex_levelling.h"
#include "initialize.h"
#include "ss_min_function.h"
#include "pp_min_function.h"
#include "dump_function.h"
#include "PGE_function.h"
#include "phase_update_function.h"
#include "all_solution_phases.h"
#include "MAGEMin.h"

#include "MAGEMin_api.h"

MAGEMin_Handle *MAGEMin_Init(	const char *database,
								int         verbose			){

	if (database != NULL && strlen(database) > 4){
		printf(" MAGEMin_Init error: database acronym '%s' is too long\n",database);
		return NULL;
	}

	MAGEMin_Handle *h = malloc(sizeof(MAGEMin_Handle));

	h->gv = global_variable_alloc(		&h->z_b					);

	h->gv.verbose = verbose;
	if (database != NULL){
		strcpy(h->gv.db,database);
	}

	h->gv = SetupDatabase(				h->gv,
									   &h->z_b					);

	h->gv = global_variable_init(		h->gv,
									   &h->z_b					);

	h->DB = InitializeDatabases(		h->gv,
										h->gv.EM_database		);

	init_simplex_A(					&h->splx_data,
										 h->gv					);
	init_simplex_B_em(					&h->splx_data,
										 h->gv					);

	return h;
}

int MAGEMin_NOxides(			MAGEMin_Handle *h				){
	return h->gv.len_ox;
}

char **MAGEMin_OxideNames(		MAGEMin_Handle *h				){
	return h->gv.ox;
}

stb_system *MAGEMin_ComputeEquilibrium(	MAGEMin_Handle *h,
											double          P,
											double          T,
											const double   *bulk,
											const char     *sys_in		){

	if (strcmp(sys_in,"mol") != 0 && strcmp(sys_in,"wt") != 0){
		printf(" MAGEMin_ComputeEquilibrium error: sys_in must be \"mol\" or \"wt\"\n");
		return NULL;
	}
	strcpy(h->gv.sys_in,sys_in);

	h->z_b.P = P;
	h->z_b.T = T + 273.15;

	for (int i = 0; i < h->gv.len_ox; i++){
		h->gv.arg_bulk[i] = bulk[i];
	}

	/* dummy input_data: never dereferenced because gv.File stays "none",
	   which skips the file-based branch inside retrieve_bulk_PT */
	io_data dummy_input_data[1];
	memset(dummy_input_data,0,sizeof(dummy_input_data));

	h->z_b = retrieve_bulk_PT(			h->gv,
										dummy_input_data,
										0,
										h->z_b					);

	h->gv = reset_gv(					h->gv,
										h->z_b,
										h->DB.PP_ref_db,
										h->DB.SS_ref_db			);

	h->z_b = reset_z_b_bulk(			h->gv,
										h->z_b					);

	reset_simplex_A(				   &h->splx_data,
										h->z_b,
										h->gv					);
	reset_simplex_B_em(				   &h->splx_data,
										h->gv					);

	reset_cp(							h->gv,
										h->z_b,
										h->DB.cp				);

	reset_SS(							h->gv,
										h->z_b,
										h->DB.SS_ref_db			);

	reset_sp(							h->gv,
										h->DB.sp				);

	h->gv = ComputeG0_point(			h->gv.EM_database,
										h->z_b,
										h->gv,
										h->DB.PP_ref_db,
										h->DB.SS_ref_db			);

	io_data dummy_point;
	memset(&dummy_point,0,sizeof(dummy_point));

	h->gv = ComputeEquilibrium_Point(	h->gv.EM_database,
										dummy_point,
										h->z_b,
										h->gv,

									   &h->splx_data,
										h->DB.PP_ref_db,
										h->DB.SS_ref_db,
										h->DB.cp				);

	h->gv = ComputePostProcessing(		h->z_b,
										h->gv,
										h->DB.PP_ref_db,
										h->DB.SS_ref_db,
										h->DB.cp				);

	fill_output_struct(					h->gv,
									   &h->splx_data,
										h->z_b,

										h->DB.PP_ref_db,
										h->DB.SS_ref_db,
										h->DB.cp,
										h->DB.sp				);

	return &h->DB.sp[0];
}

void MAGEMin_Free(				MAGEMin_Handle *h				){
	if (h == NULL) return;

	FreeDatabases(						h->gv,
										h->DB,
										h->z_b					);

	free(h);
}
