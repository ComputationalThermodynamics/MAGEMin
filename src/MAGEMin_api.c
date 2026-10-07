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

/* research group (CLI --rg) a database acronym belongs to */
static const char *research_group_of(	const char *database		){
	if (database == NULL){ return "tc"; }
	if (	strcmp(database,"sb11") == 0 ||
			strcmp(database,"sb21") == 0 ||
			strcmp(database,"sb24") == 0		){ return "sb"; }
	if (	strcmp(database,"xMELTS") == 0 ||
			strcmp(database,"pMELTS") == 0 ||
			strcmp(database,"rMELTS") == 0		){ return "gh"; }
	if (	strcmp(database,"po") == 0			){ return "br"; }
	return "tc";
}

MAGEMin_Handle *MAGEMin_Init(	const char *database,
								int         verbose			){

	if (database != NULL && strlen(database) >= len_gv_db){
		printf(" MAGEMin_Init error: database acronym '%s' is too long\n",database);
		return NULL;
	}

	MAGEMin_Handle *h = malloc(sizeof(MAGEMin_Handle));

	h->gv = global_variable_alloc(		&h->z_b					);

	h->gv.verbose = verbose;
	if (database != NULL){
		snprintf(h->gv.db, len_gv_db, "%s", database);
	}
	snprintf(h->gv.research_group, len_gv_research_group, "%s", research_group_of(database));

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

	h->ss_suppressed = calloc(h->gv.len_ss, sizeof(int));
	h->pp_suppressed = calloc(h->gv.len_pp, sizeof(int));
	h->n_suppressed  = 0;
	if (	(h->gv.len_ss > 0 && h->ss_suppressed == NULL) ||
			(h->gv.len_pp > 0 && h->pp_suppressed == NULL)		){
		printf(" MAGEMin_Init error: out of memory\n");
		MAGEMin_Free(h);
		return NULL;
	}

	return h;
}

int MAGEMin_NOxides(			MAGEMin_Handle *h				){
	return h->gv.len_ox;
}

char **MAGEMin_OxideNames(		MAGEMin_Handle *h				){
	return h->gv.ox;
}

int MAGEMin_SetSolver(			MAGEMin_Handle *h,
								int             solver			){

	if (h == NULL || solver < 0 || solver > 3){
		printf(" MAGEMin_SetSolver error: solver must be 0, 1, 2 or 3\n");
		return -1;
	}

	/* SetupDatabase forces the legacy solver for these research groups */
	if (strcmp(h->gv.research_group,"tc") != 0){
		return solver == 0 ? 0 : 1;
	}

	h->gv.solver = solver;
	return 0;
}

int MAGEMin_SetBuffer(			MAGEMin_Handle *h,
								const char     *buffer,
								double          buffer_n		){

	if (h == NULL){ return -1; }

	if (buffer == NULL || strcmp(buffer,"none") == 0){
		snprintf(h->gv.buffer, len_gv_buffer, "%s", "none");
		h->gv.buffer_n = 0.0;
		return 0;
	}

	/* the CLI --buffer options (toolkit.c) plus "iw"; any other pure-phase
	   name would be accepted by MAGEMin but silently ignored */
	static const char *buffers[] = {	"O2","qfm","mw","qif","nno","hm","iw","cco",
										"aH2O","aO2","aMgO","aFeO","aAl2O3","aTiO2"	};
	int n_buffers = sizeof(buffers) / sizeof(buffers[0]);

	int found = 0;
	for (int k = 0; k < n_buffers; k++){
		if (strcmp(buffer,buffers[k]) == 0){ found = 1; break; }
	}
	if (found == 1){
		found = 0;
		for (int i = 0; i < h->gv.len_pp; i++){
			if (strcmp(buffer,h->gv.PP_list[i]) == 0){ found = 1; break; }
		}
	}
	if (found == 0 || strlen(buffer) >= len_gv_buffer){
		printf(" MAGEMin_SetBuffer error: buffer '%s' is not available for database '%s'\n",buffer,h->gv.db);
		return -1;
	}

	snprintf(h->gv.buffer, len_gv_buffer, "%s", buffer);
	h->gv.buffer_n = buffer_n;
	return 0;
}

/* index of name in list[0..n-1], or -1 */
static int find_name(			const char     *name,
								char          **list,
								int             n				){
	for (int i = 0; i < n; i++){
		if (strcmp(name,list[i]) == 0){ return i; }
	}
	return -1;
}

int MAGEMin_SetSuppressedPhases(	MAGEMin_Handle  *h,
									const char     **names,
									int              n				){

	if (h == NULL){ return -1; }
	if (names == NULL){ n = 0; }

	/* validate everything first, so a bad name leaves the previous set intact */
	for (int k = 0; k < n; k++){
		if (	names[k] == NULL ||
				(	find_name(names[k], h->gv.SS_list, h->gv.len_ss) < 0 &&
					find_name(names[k], h->gv.PP_list, h->gv.len_pp) < 0	)	){
			printf(" MAGEMin_SetSuppressedPhases error: unknown phase '%s' for database '%s'\n",
					names[k] == NULL ? "(null)" : names[k], h->gv.db);
			return -1;
		}
	}

	for (int i = 0; i < h->gv.len_ss; i++){ h->ss_suppressed[i] = 0; }
	for (int i = 0; i < h->gv.len_pp; i++){ h->pp_suppressed[i] = 0; }
	h->n_suppressed = 0;

	for (int k = 0; k < n; k++){
		int i = find_name(names[k], h->gv.SS_list, h->gv.len_ss);
		if (i >= 0){ h->ss_suppressed[i] = 1; }
		else       { h->pp_suppressed[find_name(names[k], h->gv.PP_list, h->gv.len_pp)] = 1; }
	}
	for (int i = 0; i < h->gv.len_ss; i++){ h->n_suppressed += h->ss_suppressed[i]; }
	for (int i = 0; i < h->gv.len_pp; i++){ h->n_suppressed += h->pp_suppressed[i]; }

	return 0;
}

int MAGEMin_NSolutionPhases(	MAGEMin_Handle *h				){
	return h->gv.len_ss;
}

char **MAGEMin_SolutionPhaseNames(	MAGEMin_Handle *h			){
	return h->gv.SS_list;
}

int MAGEMin_NPurePhases(		MAGEMin_Handle *h				){
	return h->gv.len_pp;
}

char **MAGEMin_PurePhaseNames(	MAGEMin_Handle *h				){
	return h->gv.PP_list;
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
	snprintf(h->gv.sys_in, len_gv_sys_in, "%s", sys_in);

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

	/* mbCpx/mbIlm/mpSp/mpIlm select which variant of a near-degenerate model
	   pair gets pseudocompounds at all; 2 generates both, so suppressing one
	   variant leaves the other reachable. Assigned on every call so a reused
	   handle never carries it over from an earlier suppressing call. */
	int both_variants = h->n_suppressed > 0 ? 2 : 0;
	h->gv.mbCpx = both_variants;
	h->gv.mbIlm = both_variants;
	h->gv.mpSp  = both_variants;
	h->gv.mpIlm = both_variants;

	h->gv = ComputeG0_point(			h->gv.EM_database,
										h->z_b,
										h->gv,
										h->DB.PP_ref_db,
										h->DB.SS_ref_db			);

	/* phase suppression: reset_gv cleared the flags and ComputeG0_point just
	   re-activated them, so this is the only place where it takes effect
	   before ComputeEquilibrium_Point reads them */
	for (int i = 0; i < h->gv.len_ss; i++){
		if (h->ss_suppressed[i] == 1){
			for (int f = 0; f < h->gv.n_flags; f++){ h->DB.SS_ref_db[i].ss_flags[f] = 0; }
		}
	}
	for (int i = 0; i < h->gv.len_pp; i++){
		if (h->pp_suppressed[i] == 1){
			for (int f = 0; f < h->gv.n_flags; f++){ h->gv.pp_flags[i][f] = 0; }
		}
	}

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

	free(h->ss_suppressed);
	free(h->pp_suppressed);

	FreeDatabases(						h->gv,
										h->DB,
										h->z_b,
									   &h->splx_data			);

	free(h);
}
