/**
Mineral Assemblage Gibbs Energy Minimization		  
--------------------------------------------

Contributors: 

- Main developers: Riel N., Kaus. B.
- Database translation and debugging: Green E., Berlie N., and Rummel L. 
     
Contacts: nriel[at]uni-mainz.de, kaus[at]uni-mainz.de 		 
                                                                                                  
MAGEMin is written as a parallel C library callable from any petrological/geodynamic tool. For a given set of pressure, temperature and bulk-rock composition MAGEMin uses a combination of linear programming, extended Partitioning Gibbs free Energy and gradient-based local minimization to compute the most stable mineral assemblage     
      
Available thermodynamic dataset                       
================================
 
Igneous thermodynamic dataset
*****************************
                    
- Holland et al., 2018 (see http://hpxeosandthermocalc.org)
- K2O-Na2O-CaO-FeO-MgO-Al2O3-SiO2-H2O-TiO2-O-Cr2O3 chemical system
- Equations of state for
	- Pure stoichiometric phases quartz (q), cristobalite (crst), tridymite (trd), coesite (coe), stishovite (stv), kyanite (ky), sillimanite (sill), andalusite (and), rutile (ru) and sphene (sph). 
	- Solution phases spinel (spn), biotite (bi), cordierite (cd), clinopyroxene (cpx), orthopyroxene (opx), epidote (ep), garnet (g), hornblende (hb), ilmenite (ilm), silicate melt (liq), muscovite (mu), olivine (ol), ternary feldspar (fsp), and aqueous fluid (fl).
                                                    
Imported libraries                       
==================

- LAPACKE (C version of LAPACK)                         
- NLopt  (https://nlopt.readthedocs.io/)                
- uthash (https://troydhanson.github.io/uthash/)        
- ketopt (https://github.com/attractivechaos/klib/blob/master/ketopt.h) 
   
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <complex.h> 

#include "uthash.h"
#include "ketopt.h"
#include "nlopt.h"                  // requires specifying this in the makefile
#include "mpi.h"
#include "Endmembers_tc-ds62.h"
#include "Endmembers_tc-ds633.h"
#include "Endmembers_tc-ds634.h"
#include "Endmembers_M2017.h"
#include "toolkit.h"
#include "io_function.h"
#include "gem_function.h"
#include "gss_init_function.h"
#include "gss_function.h"
#include "objective_functions.h"
#include "NLopt_opt_function.h"
#include "simplex_levelling.h"
#include "initialize.h"
#include "ss_min_function.h"
#include "pp_min_function.h"
#include "dump_function.h"
#include "PGE_function.h"
#include "phase_update_function.h"
#include "MAGEMin.h"
#include "simplex_levelling.h"

/** 
  Main routine
*/
int main(		int    argc, 
				char **argv
){
	int 	rank;
	
	MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	
	/* call the main MAGEMin routine */
	runMAGEMin(argc, argv);
   
	MPI_Finalize();
	return 0;
}

/** 
  call MAGEMin
*/
int runMAGEMin(			int    argc, 
						char **argv
){
	int 	i,j,k;
	int 	rank, numprocs;

	double 	time_taken;

	/* Select the endmember database */
   	
	clock_t t = clock(),u = clock(); 

	/*
	  initialize MPI communicators 
	*/
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	/*
		initiliaze structures
	*/
	global_variable gv;
	bulk_info 		z_b;	
	Databases 		DB;

	/*
	  allocate global variables 
	*/	
	gv = global_variable_alloc(		&z_b );
	

	/** 
		Read command-line arguments and set default parameters
	*/
	gv = ReadCommandLineOptions(	 gv,
									&z_b,
									 argc,
									 argv			);


	/*
	  initialize global structure to store shared variables (e.g. Gamma, SS and PP list, ...) 
	*/
	gv = global_variable_init( 		 gv,
									&z_b 			);


	/* 
	  Allocate both pure and solid-solution databases 
	*/
	DB = InitializeDatabases(		 gv,
									 gv.EM_database	);
	

	/*
	  initialize simplex (levelling stage using pseudocompounds) 
	*/
	simplex_data 					 splx_data;
	init_simplex_A(			   		&splx_data,
									 gv				);
	init_simplex_B_em(				&splx_data,
								 	 gv				);

	/*
	  initialize output	
	*/
	dump_init( gv );	

	/* 
	  get data from input file 
	*/
	io_data input_data[gv.n_points];
	if (strcmp( gv.File, "none") != 0){	
		read_in_data(gv, input_data, gv.n_points);			
	}

	/* 
	  get bulk rock composition parsed from args 
	*/
	if 		(gv.EM_database == 0){
		gv = get_bulk_metapelite( 	gv );
	}
	else if (gv.EM_database == 1){
		gv = get_bulk_metabasite( 	gv );
	}
	else if (gv.EM_database == 2){
		gv = get_bulk_igneous( 		gv );
	}
	else if (gv.EM_database == 4){
		gv = get_bulk_ultramafic( 	gv );
	}
	else{
		printf(" Wrong database...\n");
	}
	
	/****************************************************************************************/
	/**                               LAUNCH MINIMIZATION ROUTINE                          **/
	/****************************************************************************************/
	if (rank==0 && gv.verbose != -1){
    	printf("\nRunning MAGEMin %5s on %d cores {\n", gv.version, numprocs);
    	printf("═══════════════════════════════════════════════\n");
	}
	for (int point = 0; point < gv.n_points; point++){
        if ( (point % numprocs != rank)) continue;   	/** this ensures that, in parallel, not every point is computed by every processor (instead only every numprocs point). Only applied to Mode==0 */

		t              = clock();										/** reset loop timer 				*/
		gv.numPoint    = point; 										/** the number of the current point */

		z_b = retrieve_bulk_PT(				gv,
											input_data,
											point,
											z_b							);

		/* reset global variables flags 											*/
		gv = reset_gv(						gv,
											z_b,
											DB.PP_ref_db,
											DB.SS_ref_db				);

		/** reset bulk rock information (needed for parallel point calculation) 	*/
		z_b = reset_z_b_bulk(				gv,				
											z_b							);	
	
		/** reset simplex memory 													*/
		reset_simplex_A(			   	   &splx_data,
											z_b,
											gv							);
											
		reset_simplex_B_em(				   &splx_data,
											gv							);
				
		/** reset considered phases structure 										*/
		reset_cp(							gv,												
											z_b,
											DB.cp						);	
											
		/** reset pure and solution phases											*/
		reset_SS(							gv,												
											z_b,
											DB.SS_ref_db				);	
		
		/** reset stable phases														*/
		reset_sp(							gv,
											DB.sp						);

		/* Perform calculation for a single point 									*/	
		gv = ComputeG0_point(				gv.EM_database, 
											z_b,											/** bulk rock informations 			*/
											gv,												/** global variables (e.g. Gamma) 	*/

											DB.PP_ref_db,									/** pure phase database 			*/
											DB.SS_ref_db				);

		gv = ComputeEquilibrium_Point(		gv.EM_database, 
											input_data[point],
											z_b,											/** bulk rock informations 			*/
											gv,												/** global variables (e.g. Gamma) 	*/
											
										   &splx_data,
											DB.PP_ref_db,									/** pure phase database 			*/
											DB.SS_ref_db,									/** solid solution database 		*/
											DB.cp						);

		gv = ComputePostProcessing(			z_b,											/** bulk rock informations 			*/
											gv,												/** global variables (e.g. Gamma) 	*/
											DB.PP_ref_db,									/** pure phase database 			*/
											DB.SS_ref_db,									/** solid solution database 		*/
											DB.cp						);

		/* Fill structure holding stable phase equilibrium informations 			*/
		fill_output_struct(					gv,												/** global variables (e.g. Gamma) 	*/
										   &splx_data,
										   
											z_b,											/** bulk-rock informations 			*/
											DB.PP_ref_db,									/** pure phase database 			*/
											DB.SS_ref_db,									/** solution phase database 		*/
											DB.cp,
											DB.sp						);

		t 			= clock() - t; 
		time_taken 	= ((double)t)/CLOCKS_PER_SEC;
		gv.tot_time = time_taken*1000.0;

		/* Dump final results to files 												*/
		save_results_function(				gv,												/** global variables (e.g. Gamma) 	*/
											z_b,											/** bulk-rock informations 			*/
											DB.PP_ref_db,									/** pure phase database 			*/
											DB.SS_ref_db,									/** solution phase database 		*/
											DB.cp,
											DB.sp						);

		/* Print output to screen 													*/									/* in seconds 	 					*/
		PrintOutput(gv, rank, point, DB, time_taken, z_b);									/* print output on screen 			*/
	}
	/* end of loop over points */

	/* wait for all cores to be finished */
	MPI_Barrier(MPI_COMM_WORLD);		

	/* now merge the parallel output files into one*/
	mergeParallelFiles(gv);

	if (gv.output_matlab >= 1){
		mergeParallel_matlab(gv);
	}
	/* free memory allocated to solution and pure phases */
	FreeDatabases(gv, DB, z_b);

	/* print the time */
	u = clock() - u; 
	
	if (gv.verbose != -1){
		time_taken = ((double)u)/(CLOCKS_PER_SEC); 				/** in seconds */
		if (rank==0){
			printf("___________________________________\n");
			printf("MAGEMin comp time: %+3f ms }\n", time_taken*1000.0);
		}
	}
    return 0;
}

/** 
  Compute stable equilibrium at given P/T/C point
*/
	global_variable ComputePostProcessing(			bulk_info 	 		 z_b,
													global_variable 	 gv,
													PP_ref  			*PP_ref_db,
													SS_ref  			*SS_ref_db,
													csd_phase_set  		*cp					){
	

	gv = compute_phase_mol_fraction(	gv,
										PP_ref_db,
										SS_ref_db,
										cp				);


	gv = compute_activities (			gv.EM_database,	
										gv,
										PP_ref_db,
										z_b				);

	gv = compute_density_volume_modulus(gv.EM_database,
										z_b,										/** bulk rock informations 			*/
										gv,											/** global variables (e.g. Gamma) 	*/
										PP_ref_db,									/** pure phase database 			*/
										SS_ref_db,									/** solid solution database 		*/
										cp			);


	return gv;
}

/** 
  Compute the Gibbs energy of reference for PT conditions, for pure phases and solution phases
*/
	global_variable ComputeG0_point( 	int 				 EM_database,
										bulk_info 	 		 z_b,
										global_variable 	 gv,
										PP_ref  			*PP_ref_db,
										SS_ref  			*SS_ref_db				){

	/* initialize endmember database for given P-T point */
	gv = init_em_db(		EM_database,
							z_b,											/** bulk rock informations 			*/
							gv,												/** global variables (e.g. Gamma) 	*/
							PP_ref_db						);

	/* Calculate solution phase data at given P-T conditions (G0 based on G0 of endmembers) */
	gv = init_ss_db(		EM_database,
							z_b,
							gv,
							SS_ref_db						);

	return gv;
}



/** 
  Compute to levelling only (for testing)
*/
	global_variable ComputeLevellingOnly( 		int 				 EM_database,
												io_data 			 input_data,
												bulk_info 	 		 z_b,
												global_variable 	 gv,

												simplex_data	    *splx_data,
												PP_ref  			*PP_ref_db,
												SS_ref  			*SS_ref_db,
												csd_phase_set  		*cp						){

	/** pointer array to objective functions 								*/
	obj_type 								SS_objective[gv.len_ss];	

	if (EM_database == 0){				// metapelite database //
		SS_mp_objective_init_function(			SS_objective,
												gv							);
	}
	if (EM_database == 1){				// metabasite database //
		SS_mb_objective_init_function(			SS_objective,
												gv							);
	}
	else if (EM_database == 2){			// igneous database //
		SS_ig_objective_init_function(			SS_objective,
												gv							);
	}
	else if (EM_database == 4){			// ultramafic database //
		SS_um_objective_init_function(			SS_objective,
												gv							);
	}
	
	/****************************************************************************************/
	/**                                   LEVELLING                                        **/
	/****************************************************************************************/	
	// leveling mode = 0 is default, without using initial guess
	// leveling mode = 1 uses initial guess

	gv = Levelling(			z_b,										/** bulk rock informations 			*/
							gv,											/** global variables (e.g. Gamma) 	*/

							SS_objective,
							splx_data,
							PP_ref_db,									/** pure phase database 			*/
							SS_ref_db,									/** solution phase database 		*/
							cp							);

	return gv;

	}


/** 
  Compute stable equilibrium at given Pressure, Temperature and bulk-rock composition
*/
	global_variable ComputeEquilibrium_Point( 	int 				 EM_database,
												io_data 			 input_data,
												bulk_info 	 		 z_b,
												global_variable 	 gv,

												simplex_data	    *splx_data,
												PP_ref  			*PP_ref_db,
												SS_ref  			*SS_ref_db,
												csd_phase_set  		*cp						){

	/** pointer array to objective functions 								*/
	obj_type 								SS_objective[gv.len_ss];	

	if (EM_database == 0){				// metapelite database //
		SS_mp_objective_init_function(			SS_objective,
												gv							);
	}
	if (EM_database == 1){				// metabasite database //
		SS_mb_objective_init_function(			SS_objective,
												gv							);
	}
	else if (EM_database == 2){			// igneous database //
		SS_ig_objective_init_function(			SS_objective,
												gv							);
	}
	else if (EM_database == 4){			// ultramafic database //
		SS_um_objective_init_function(			SS_objective,
												gv							);
	}
	

	/****************************************************************************************/
	/**                                   LEVELLING                                        **/
	/****************************************************************************************/	
	// leveling mode = 0 is default, without using initial guess
	// leveling mode = 1 uses initial guess
	if (gv.leveling_mode == 0){
		gv = Levelling(			z_b,										/** bulk rock informations 			*/
								gv,											/** global variables (e.g. Gamma) 	*/

								SS_objective,
								splx_data,
								PP_ref_db,									/** pure phase database 			*/
								SS_ref_db,									/** solution phase database 		*/
								cp							);

	}
	else if (gv.leveling_mode == 1){
		gv = Initial_guess(		z_b,										/** bulk rock informations 			*/
								gv,											/** global variables (e.g. Gamma) 	*/

								splx_data,
								PP_ref_db,									/** pure phase database 			*/
								SS_ref_db,									/** solution phase database 		*/
								cp							);
	}


	/****************************************************************************************/
	/**                            PARTITIONING GIBBS ENERGY                               **/
	/****************************************************************************************/
	
	if (gv.solver == 0){ 		/* Legacy solver only */

		gv.div 		= 0;
		gv.status 	= -1;

		/* initialize legacy solver using results of levelling phase */
		for (int i = 0; i < gv.len_ox; i++){
			gv.gam_tot[i] = gv.gam_tot_0[i];
		}	

		gv = LP(				z_b,									/** bulk rock informations 			*/
								gv,										/** global variables (e.g. Gamma) 	*/

								SS_objective,
								splx_data,
								PP_ref_db,								/** pure phase database 			*/
								SS_ref_db,								/** solution phase database 		*/
								cp						);

	}
	else if (gv.solver == 1){	/* PGE + Legacy solver */

		if (z_b.T > gv.solver_switch_T){
			gv 		= PGE(			z_b,									/** bulk rock constraint 			*/ 
									gv,										/** global variables (e.g. Gamma) 	*/

									SS_objective,
									splx_data,
									PP_ref_db,								/** pure phase database 			*/
									SS_ref_db,								/** solution phase database 		*/
									cp							);

		}

		/**
			Launch legacy solver (LP, Theriak-Domino like algorithm)
		*/ 
		if ((gv.div == 1 || z_b.T <= gv.solver_switch_T ) && gv.solver == 1){
			if (gv.verbose == 1){
				if (gv.div == 1){	
					printf("\n[PGE failed (residual: %+.4f, PT: [%+.5f,%+.5f])-> legacy solver...]\n",gv.BR_norm,z_b.P,z_b.T-273.15);
				}
				if (z_b.T <= gv.solver_switch_T){	
					printf("\n Low Temperature conditions (T < %+5f) -> legacy solver...\n",gv.solver_switch_T);
				}
			}
			gv.div 		= 0;
			gv.status 	= -1;

			/* initialize legacy solver using results of levelling phase */
			for (int i = 0; i < gv.len_ox; i++){
				gv.gam_tot[i] = gv.gam_tot_0[i];
			}	

			gv = LP(				z_b,									/** bulk rock informations 			*/
									gv,										/** global variables (e.g. Gamma) 	*/

									SS_objective,
									splx_data,
									PP_ref_db,								/** pure phase database 			*/
									SS_ref_db,								/** solution phase database 		*/
									cp						);
		}
		else{	// here we compute the LP initial guess
			gv = run_LP_ig(					z_b,
											splx_data,
											gv,
													
											PP_ref_db,
											SS_ref_db			);
		}
	}
	else if (gv.solver == 2){
		for (int i = 0; i < gv.len_ss; i++){
			gv.n_ss_array[i] = 0;
		}

		double 	ig_liq = 0.0;
		int 	n_liq  = 0;
		int 	n_ss   = 0;

		for (int i = 0; i < gv.len_cp; i++){ 
			if (cp[i].ss_flags[1] == 1){
				gv.n_ss_array[cp[i].id] += 1;
				if (strcmp( gv.SS_list[cp[i].id], "liq")  == 0){
					ig_liq += cp[i].ss_n;
					n_liq  += 1;
				}
			}
		}

		for (int i = 0; i < gv.len_ss; i++){
			if (gv.n_ss_array[i] > n_ss){
				n_ss = gv.n_ss_array[i];
			}
		}

		gv.div 		= 0;
		gv.status 	= 0;

		/* initialize legacy solver using results of levelling phase */
		for (int i = 0; i < gv.len_ox; i++){
			gv.gam_tot[i] = gv.gam_tot_0[i];
		}	

		if (ig_liq > 0.25 || n_liq > 2 || n_ss > 6){
			gv 		= PGE(			z_b,									/** bulk rock constraint 			*/ 
									gv,										/** global variables (e.g. Gamma) 	*/

									SS_objective,
									splx_data,
									PP_ref_db,								/** pure phase database 			*/
									SS_ref_db,								/** solution phase database 		*/
									cp							);

			if (gv.div == 1){
				if (gv.verbose == 1){
					printf("\n[PGE failed (residual: %+.4f, PT: [%+.5f,%+.5f])-> legacy solver...]\n",gv.BR_norm,z_b.P,z_b.T-273.15);
				}
				gv.div 		= 0;
				gv.status 	= -1;
				
				/* initialize legacy solver using results of levelling phase */
				for (int i = 0; i < gv.len_ox; i++){
					gv.gam_tot[i] = gv.gam_tot_0[i];
				}	

				gv = LP(				z_b,									/** bulk rock informations 			*/
										gv,										/** global variables (e.g. Gamma) 	*/

										SS_objective,
										splx_data,
										PP_ref_db,								/** pure phase database 			*/
										SS_ref_db,								/** solution phase database 		*/
										cp						);
			}
			else{	// here we compute the LP initial guess
				gv = run_LP_ig(					z_b,
												splx_data,
												gv,
														
												PP_ref_db,
												SS_ref_db			);
			}
		}
		else{

			gv.div 		= 0;
			gv.status 	= 0;

			/* initialize legacy solver using results of levelling phase */
			for (int i = 0; i < gv.len_ox; i++){
				gv.gam_tot[i] = gv.gam_tot_0[i];
			}	

			gv = LP(				z_b,									/** bulk rock informations 			*/
									gv,										/** global variables (e.g. Gamma) 	*/

									SS_objective,
									splx_data,
									PP_ref_db,								/** pure phase database 			*/
									SS_ref_db,								/** solution phase database 		*/
									cp						);
		}
	}
	else {
		printf("  Wrong solver option: should be 0 (legacy) or 1 (PGE & legacy)");
	}

	if (gv.verbose == 1){
		gv = check_PC_driving_force( 	z_b,							/** bulk rock constraint 			*/ 
										gv,								/** global variables (e.g. Gamma) 	*/

										PP_ref_db,						/** pure phase database 			*/ 
										SS_ref_db,
										cp				); 	
		printf("\n\n\n");									
		printf("╔═════════════════════════════════════════════════════════════╗\n");
		printf("║                       COMPUTATION SUMMARY                   ║\n");
		printf("╚═════════════════════════════════════════════════════════════╝\n\n");
		printf(" Total number of computed pseudocompounds during iterations\n");
		printf("══════════════════════════════════════════════════════════════\n");
		for (int i = 0; i < gv.len_ss; i++){
			printf("%5s n_Ppc : %5d\n",gv.SS_list[i],SS_ref_db[i].id_Ppc);
		}
		printf("\n");
	

		printf(" Alg | ite  | duration   |  MASS norm | Gamma norm | Gibbs sys\n");
		printf("══════════════════════════════════════════════════════════════\n");

		for (int i = 0; i < gv.global_ite; i++){	
			if (gv.Alg[i] == 0){
				printf(" LP  | %4d | %+10f | %+10f | %+10f | %+10f\n",i,gv.ite_time[i],gv.PGE_mass_norm[i],gv.gamma_norm[i],gv.gibbs_ev[i]);
			}
			if (gv.Alg[i] == 1){
				printf(" PGE | %4d | %+10f | %+10f | %+10f | %+10f\n",i,gv.ite_time[i],gv.PGE_mass_norm[i],gv.gamma_norm[i],gv.gibbs_ev[i]);
			}	
			if (gv.Alg[i+1] - gv.Alg[i] == 1){
				printf("-----------------------------------------------------------------\n");
				printf("                          SWITCH FROM LP TO PGE                  \n");
				printf("-----------------------------------------------------------------\n");
			}
			if (gv.Alg[i+1] - gv.Alg[i] == -1 && i < gv.global_ite - 1){
				printf("-----------------------------------------------------------------\n");
				printf("                          SWITCH FROM PGE TO LP                  \n");
				printf("-----------------------------------------------------------------\n");
			}					
		}
		printf("\n");
	}

	return gv;
}

/** 
  	Get command line options
*/
global_variable ReadCommandLineOptions(	global_variable 	 gv,
										bulk_info 			*z_b,		
										int 				 argc, 
										char 			   **argv
){
	int i;
	static ko_longopt_t longopts[] = {
		{ "Verb", 		ko_optional_argument, 301 },
		{ "db", 		ko_optional_argument, 302 },
		{ "File", 		ko_optional_argument, 303 },
		{ "n_points ",	ko_optional_argument, 304 },
		{ "test",  		ko_optional_argument, 305 },
		{ "Temp", 		ko_optional_argument, 306 },
		{ "Pres",  		ko_optional_argument, 307 },
		{ "Phase",  	ko_optional_argument, 308 },
		{ "Gam",  		ko_optional_argument, 310 },
		{ "Bulk", 		ko_optional_argument, 311 },
        { "maxeval",    ko_optional_argument, 313 },
        { "version",    ko_optional_argument, 314 },
        { "help",    	ko_optional_argument, 315 },
        { "solver",    	ko_optional_argument, 316 },
        { "out_matlab", ko_optional_argument, 318 },
        { "sys_in",    	ko_optional_argument, 317 },
        { "buffer",    	ko_optional_argument, 319 },
        { "buffer_n",   ko_optional_argument, 320 },
     	{ "limitCaOpx", ko_optional_argument, 321 },
     	{ "CaOpxLim",   ko_optional_argument, 326 },
     	{ "fluidSpec",  ko_optional_argument, 322 },
     	{ "mbCpx",  	ko_optional_argument, 323 },
    	{ NULL, 0, 0 }
	};
	ketopt_t opt = KETOPT_INIT;

	int    c;
	while ((c = ketopt(&opt, argc, argv, 1, "", longopts)) >= 0) {
		if 		(c == 314){ printf("MAGEMin %20s\n",gv.version ); exit(0); }	
		else if (c == 315){ print_help( gv ); 					  exit(0); }	
 
        else if	(c == 301){ gv.verbose  		= atoi(opt.arg);				if (gv.verbose == 1){		printf("--verbose     : verbose              = %i \n", 	 	   		gv.verbose			);}} 	
		else if (c == 321){ gv.limitCaOpx   	= atoi(opt.arg);				if (gv.verbose == 1){		printf("--limitCaOpx  : limitCaOpx           = %i \n", 	 	   		gv.limitCaOpx		);}} 	
		else if (c == 326){ gv.CaOpxLim	   		= strtold(opt.arg,NULL); 		if (gv.verbose == 1){		printf("--CaOpxLim    : CaOpxLim             = %f \n", 	 	   		gv.CaOpxLim			);}} 	
		else if (c == 322){ gv.fluidSpec    	= atoi(opt.arg);				if (gv.verbose == 1){		printf("--fluidSpec   : fluidSpec            = %i \n", 	 	   		gv.fluidSpec		);}} 	
		else if (c == 323){ gv.mbCpx   			= atoi(opt.arg);				if (gv.verbose == 1){		printf("--mbCpx       : mbCpx                = %i \n", 	 	   		gv.mbCpx			);}} 	
		else if (c == 316){ gv.solver   		= atoi(opt.arg);				if (gv.verbose == 1){		printf("--solver      : solver               = %i \n", 	 	   		gv.solver			);}}																		
		else if (c == 318){ gv.output_matlab   	= atoi(opt.arg); 				if (gv.verbose == 1){		printf("--out_matlab  : out_matlab           = %i \n", 	 	   		gv.output_matlab	);}}																		
		else if (c == 304){ gv.n_points 		= atoi(opt.arg); 	 			if (gv.verbose == 1){		printf("--n_points    : n_points             = %i \n", 	 	   		gv.n_points			);}}
		else if (c == 305){ gv.test  			= atoi(opt.arg); 				if (gv.verbose == 1){		printf("--test        : Test                 = %i \n", 	 	  		gv.test				);}}
		else if (c == 320){ gv.buffer_n 		= strtold(opt.arg,NULL); 		if (gv.verbose == 1){		printf("--buffer_n    : buffer_n             = %f \n", 				gv.buffer_n			);}}
		else if (c == 306){ z_b->T				= strtold(opt.arg,NULL)+273.15;	if (gv.verbose == 1){		printf("--Temp        : Temperature          = %f C \n",            z_b->T-273.15000	);}}
		else if (c == 307){ z_b->P 				= strtold(opt.arg,NULL); 		if (gv.verbose == 1){		printf("--Pres        : Pressure             = %f kbar \n", 		z_b->P				);}}

		else if (c == 319){ strcpy(gv.buffer,opt.arg);							if (gv.verbose == 1){		printf("--buffer      : buffer               = %s \n", 	 	   		gv.buffer			);}} 	
		else if (c == 308){ strcpy(gv.Phase,opt.arg);		 					if (gv.verbose == 1){		printf("--Phase       : Phase name           = %s \n", 	   			gv.Phase			);}}
		else if (c == 303){ strcpy(gv.File,opt.arg);		 					if (gv.verbose == 1){		printf("--File        : File                 = %s \n", 	 	   		gv.File				);}}
		else if (c == 302){ strcpy(gv.db,opt.arg);		 						if (gv.verbose == 1){		printf("--db          : db                   = %s \n", 	 	   		gv.db				);}}
		else if (c == 317){ strcpy(gv.sys_in,opt.arg);		 					if (gv.verbose == 1){		printf("--sys_in      : sys_in               = %s \n", 	 	   		gv.sys_in			);}}

		else if (c == 310){
			char *p = strtok(opt.arg,",");
			size_t i = 0;
			while(p && i<11) {
					gv.arg_gamma[i++] = atof(p);
					p = strtok(NULL, ",");
			}
			if (gv.verbose == 1){
				printf("--Gam  	      : Gamma           = ");
				for (int j = 0; j < 11; j++){
					printf("%g ", gv.arg_gamma[j]);	
				} 
				printf(" dG \n");
			}
		 }
		else if (c == 311){
			char *p  = strtok(opt.arg,",");
			size_t i = 0;
			while(p && i<11) {
					gv.arg_bulk[i++] = atof(p);
					p = strtok(NULL, ",");
			}
			if (gv.verbose == 1){
				printf("--Bulk  	 : Bulk         = ");
				for (int j = 0; j < 11; j++){
					printf("%g ", gv.arg_bulk[j]);	
				} 
				printf(" \n");
			}
		 }
	}

	/* set-up database acronym here*/
	if 		(strcmp(gv.db, "mp") 	== 0){
		gv.EM_database = 0;
	}
	else if (strcmp(gv.db, "mb") 	== 0){
		gv.EM_database = 1;
	}
	else if (strcmp(gv.db, "ig") 	== 0){
		gv.EM_database = 2;
	}
	else if (strcmp(gv.db, "um") 	== 0){
		gv.EM_database = 4;
	}
	else {
		printf(" No or wrong database acronym has been provided, using default (Igneous [ig])\n");
		gv.EM_database = 2;
	}
	return gv;
} 

	
/** 
  Initiatizes the endmember and solid solution databases and adds them to a single struct
**/
Databases InitializeDatabases(	global_variable gv, 
								int 			EM_database
){
	Databases 	DB;
	int 		i;

	/* Allocate pure-phase database (to get gbase, comp and factor) 				*/
	DB.PP_ref_db = malloc ((gv.len_pp) 		* sizeof(PP_ref)); 
	
	/** 
		Allocate memory for each solution phase according to their specificities (n_em, sf etc) 
	*/
	DB.SS_ref_db = malloc ((gv.len_ss) 		* sizeof(SS_ref)); 
	for (i = 0; i < gv.len_ss; i++){
		DB.SS_ref_db[i] = G_SS_init_EM_function(		i,	
														DB.SS_ref_db[i], 
														EM_database, 
														gv.SS_list[i], 
														gv						);
	}

	/* Allocate memory of the considered set of phases 								*/
	DB.cp 		 = malloc ((gv.max_n_cp) 	* sizeof(csd_phase_set)); 
	for (i = 0; i < gv.max_n_cp; i++){
		DB.cp[i] = CP_INIT_function(		DB.cp[i], 
											gv									);
	}

	/* Allocate memory for the stable phase equilibrium 							*/
	DB.sp 		 = malloc (1 	* sizeof(stb_system)); 
	DB.sp[0] 	 = SP_INIT_function(		DB.sp[0], gv						);

	/* Endmember names */
	DB.EM_names  =	get_EM_DB_names(		gv									);

	/* Endmember names */
	DB.FS_names  =	get_FS_DB_names(		gv									);

	/* Create endmember Hashtable */
	struct EM2id *p_s, *tmp_p;
	struct EM_db EM_return;
	int n_em_db = gv.n_em_db;
    for (int i = 0; i < n_em_db; ++i) {
		char EM_name[20];
        p_s = (struct EM2id *)malloc(sizeof *p_s);
        strcpy(p_s->EM_tag, DB.EM_names[i]);
        p_s->id = i;
        HASH_ADD_STR( EM, EM_tag, p_s );
    }

	/* Create pure-phase hashtable */
	struct PP2id *pp_s, *tmp_pp;
    for (int i = 0; i < sizeof(gv.PP_list); ++i) {
        pp_s = (struct PP2id *)malloc(sizeof *pp_s);
        strcpy(pp_s->PP_tag, gv.PP_list[i]);
        pp_s->id = i;
        HASH_ADD_STR( PP, PP_tag, pp_s );
    }

	/* Create fluid species Hashtable */
	struct FS2id *fs_s, *tmp_fs;
	struct FS_db FS_return;
	int n_fs_db = gv.n_fs_db;
    for (int i = 0; i < n_fs_db; ++i) {
		char EM_name[20];
        fs_s = (struct FS2id *)malloc(sizeof *fs_s);
        strcpy(fs_s->FS_tag, DB.FS_names[i]);
        fs_s->id = i;
        HASH_ADD_STR( FS, FS_tag, fs_s );
    }

	return DB;
}

/** 
  Free the memory associated with the databases
**/
void FreeDatabases(		global_variable gv, 
						Databases 		DB,
						bulk_info 	 	z_b			){
	int i, j, n_xeos, n_em, n_ox, n_pc, n_Ppc, n_cp, sym, ndif, pp, ss;

	/*  ==================== SP ==============================  */
	n_ox = gv.len_ox;

	for ( i = 0; i < n_ox; i++){
		if  (DB.sp[0].oxides[i]				!=NULL)  free( DB.sp[0].oxides[i] 			);	
		if  (DB.sp[0].ph[i]					!=NULL)  free( DB.sp[0].ph[i] 				);	
	}

	for ( i = 0; i < n_ox; i++){
		if  (DB.sp[0].PP[i].Comp			!=NULL)  free( DB.sp[0].PP[i].Comp 			);	
		if  (DB.sp[0].PP[i].Comp_wt			!=NULL)  free( DB.sp[0].PP[i].Comp_wt 		);	
	}

	for ( i = 0; i < n_ox; i++){
		if  (DB.sp[0].SS[i].Comp			!=NULL)  free( DB.sp[0].SS[i].Comp 			);	
		if  (DB.sp[0].SS[i].Comp_wt			!=NULL)  free( DB.sp[0].SS[i].Comp_wt 		);	
		if  (DB.sp[0].SS[i].compVariables	!=NULL)  free( DB.sp[0].SS[i].compVariables );	
		if  (DB.sp[0].SS[i].emFrac			!=NULL)  free( DB.sp[0].SS[i].emFrac 		);	
		if  (DB.sp[0].SS[i].emFrac_wt		!=NULL)  free( DB.sp[0].SS[i].emFrac_wt 	);	
		if  (DB.sp[0].SS[i].emChemPot		!=NULL)  free( DB.sp[0].SS[i].emChemPot 	);	
		for ( j = 0; j < n_ox*3; j++){
			if  (DB.sp[0].SS[i].compVariablesNames[j]	!=NULL)  free( DB.sp[0].SS[i].compVariablesNames[j] 	);	
			if  (DB.sp[0].SS[i].emNames[j]				!=NULL)  free( DB.sp[0].SS[i].emNames[j] 				);	
			if  (DB.sp[0].SS[i].emComp[j]				!=NULL)  free( DB.sp[0].SS[i].emComp[j] 				);	
			if  (DB.sp[0].SS[i].emComp_wt[j]			!=NULL)  free( DB.sp[0].SS[i].emComp_wt[j] 				);	
		}
		if  (DB.sp[0].SS[i].compVariablesNames	!=NULL)  free( DB.sp[0].SS[i].compVariablesNames );	
		if  (DB.sp[0].SS[i].emNames				!=NULL)  free( DB.sp[0].SS[i].emNames 			);	
		if  (DB.sp[0].SS[i].emComp				!=NULL)  free( DB.sp[0].SS[i].emComp 			);	
		if  (DB.sp[0].SS[i].emComp_wt			!=NULL)  free( DB.sp[0].SS[i].emComp_wt 		);	
	}

	/* free metastable assemblage */
	for ( i = 0; i < gv.max_n_mSS; i++){
		if  (DB.sp[0].mSS[i].comp_Ppc		!=NULL)  free( DB.sp[0].mSS[i].comp_Ppc		);
		if  (DB.sp[0].mSS[i].p_Ppc			!=NULL)  free( DB.sp[0].mSS[i].p_Ppc		);
		if  (DB.sp[0].mSS[i].mu_Ppc			!=NULL)  free( DB.sp[0].mSS[i].mu_Ppc		);
		if  (DB.sp[0].mSS[i].xeos_Ppc		!=NULL)  free( DB.sp[0].mSS[i].xeos_Ppc		);
		if  (DB.sp[0].mSS[i].ph_name		!=NULL)  free( DB.sp[0].mSS[i].ph_name		);
		if  (DB.sp[0].mSS[i].ph_type		!=NULL)  free( DB.sp[0].mSS[i].ph_type		);
	}

	free(DB.sp[0].PP);
	free(DB.sp[0].SS);
	free(DB.sp[0].mSS);

	free(DB.sp[0].oxides);
	free(DB.sp[0].ph);

	free(DB.sp[0].bulk);
	free(DB.sp[0].gamma);
	free(DB.sp[0].bulk_S);
	free(DB.sp[0].bulk_M);
	free(DB.sp[0].bulk_F);

	free(DB.sp[0].bulk_wt);
	free(DB.sp[0].bulk_S_wt);
	free(DB.sp[0].bulk_M_wt);
	free(DB.sp[0].bulk_F_wt);
	free(DB.sp[0].ph_frac);
	free(DB.sp[0].ph_frac_wt);
	free(DB.sp[0].ph_frac_vol);

	free(DB.sp[0].ph_id);
	free(DB.sp[0].ph_type);
	free(DB.sp[0].MAGEMin_ver);


	/*  ==================== CP ==============================  */
	n_cp = gv.max_n_cp;
	for (i = 0; i < n_cp; i++) {
		free(DB.cp[i].ss_flags);
		free(DB.cp[i].name);
		free(DB.cp[i].p_em);
		free(DB.cp[i].xi_em);
		free(DB.cp[i].dguess);
		free(DB.cp[i].xeos);
		free(DB.cp[i].xeos_0);
		free(DB.cp[i].xeos_1);
		free(DB.cp[i].xeos_r);
		free(DB.cp[i].delta_mu);
		free(DB.cp[i].dfx);
		free(DB.cp[i].mu);
		free(DB.cp[i].gbase);
		free(DB.cp[i].ss_comp);
		free(DB.cp[i].sf);
	}

	int n_em_db = gv.n_em_db;
	for (i = 0; i < n_em_db; i++) {
		free(DB.EM_names[i]);
	}
	int n_fs_db = gv.n_fs_db;
	for (i = 0; i < n_fs_db; i++) {
		free(DB.FS_names[i]);
	}

	/*  ==================== SS_ref_db ==============================  */
	ndif 	= gv.n_Diff;
	n_Ppc  	= gv.n_Ppc;

	for (i = 0; i < gv.len_ss; i++) {
		// printf("SS being freed %s\n",gv.SS_list[i]);
		n_pc 	= gv.n_SS_PC[i];
		n_em 	= DB.SS_ref_db[i].n_em;
		n_xeos 	= DB.SS_ref_db[i].n_xeos;

		free(DB.SS_ref_db[i].ss_flags);
		free(DB.SS_ref_db[i].solvus_id);

		free(DB.SS_ref_db[i].gbase);
		free(DB.SS_ref_db[i].gb_lvl);
		free(DB.SS_ref_db[i].z_em);
		free(DB.SS_ref_db[i].d_em);
		free(DB.SS_ref_db[i].density);
		free(DB.SS_ref_db[i].dguess);
		free(DB.SS_ref_db[i].iguess);
		free(DB.SS_ref_db[i].mguess);
		free(DB.SS_ref_db[i].idOrderVar);

		free(DB.SS_ref_db[i].p);
		free(DB.SS_ref_db[i].ElShearMod);
		free(DB.SS_ref_db[i].ape);
		free(DB.SS_ref_db[i].mat_phi);
		free(DB.SS_ref_db[i].mu_Gex);	
		free(DB.SS_ref_db[i].sf);
		free(DB.SS_ref_db[i].mu);
		free(DB.SS_ref_db[i].dfx);
		free(DB.SS_ref_db[i].ss_comp);
		free(DB.SS_ref_db[i].ElEntropy);
		free(DB.SS_ref_db[i].xi_em);
		free(DB.SS_ref_db[i].xeos);

		free(DB.SS_ref_db[i].ub);
		free(DB.SS_ref_db[i].lb);

		sym    = DB.SS_ref_db[i].symmetry;
		if (sym == 0){
			free(DB.SS_ref_db[i].W);
			free(DB.SS_ref_db[i].v);
		}
		else if (sym == 1){
			free(DB.SS_ref_db[i].W);
		}

		for (j = 0; j < ndif; j++) {	free(DB.SS_ref_db[i].mu_array[j]);}	free(DB.SS_ref_db[i].mu_array);

		free(DB.SS_ref_db[i].G_pc);
		free(DB.SS_ref_db[i].DF_pc);
		free(DB.SS_ref_db[i].tot_pc);
		free(DB.SS_ref_db[i].id_pc);
		free(DB.SS_ref_db[i].factor_pc);
		free(DB.SS_ref_db[i].info);

		n_em = DB.SS_ref_db[i].n_em;
		for (j = 0; j < n_pc; j++) {
			free(DB.SS_ref_db[i].p_pc[j]);
			// free(DB.SS_ref_db[i].mu_pc[j]);
			free(DB.SS_ref_db[i].comp_pc[j]);
			free(DB.SS_ref_db[i].xeos_pc[j]);			
		}
		free(DB.SS_ref_db[i].p_pc);
		// free(DB.SS_ref_db[i].mu_pc);
		free(DB.SS_ref_db[i].comp_pc);
		free(DB.SS_ref_db[i].xeos_pc);

		/** free Ppc */
		for (j = 0; j < n_em; j++){ 
			free(DB.SS_ref_db[i].EM_list[j]);
			free(DB.SS_ref_db[i].eye[j]);
			free(DB.SS_ref_db[i].Comp[j]);
			free(DB.SS_ref_db[i].dp_dx[j]);
		}
		free(DB.SS_ref_db[i].EM_list);
		free(DB.SS_ref_db[i].eye);
		free(DB.SS_ref_db[i].Comp);
		free(DB.SS_ref_db[i].dp_dx);

		for (j = 0; j < n_xeos; j++){ 
			free(DB.SS_ref_db[i].CV_list[j]);
			free(DB.SS_ref_db[i].bounds[j]);
			free(DB.SS_ref_db[i].bounds_ref[j]);
		}

		free(DB.SS_ref_db[i].CV_list);
		free(DB.SS_ref_db[i].bounds);
		free(DB.SS_ref_db[i].bounds_ref);

		free(DB.SS_ref_db[i].G_Ppc);
		free(DB.SS_ref_db[i].DF_Ppc);
		free(DB.SS_ref_db[i].info_Ppc);

		for (j = 0; j < n_Ppc; j++) {
			free(DB.SS_ref_db[i].p_Ppc[j]);
			free(DB.SS_ref_db[i].mu_Ppc[j]);
			free(DB.SS_ref_db[i].comp_Ppc[j]);
			free(DB.SS_ref_db[i].xeos_Ppc[j]);			
		}
		free(DB.SS_ref_db[i].p_Ppc);
		free(DB.SS_ref_db[i].mu_Ppc);
		free(DB.SS_ref_db[i].comp_Ppc);
		free(DB.SS_ref_db[i].xeos_Ppc);
	}

	/*  ==================== GV ==============================  */
	/* It seems like unwrapped pointers using Julia creates conflicts when de-allocating memory */
	// free(gv.outpath);
	// free(gv.version);
	// free(gv.File);
	// free(gv.db);
	// free(gv.Phase);
	// free(gv.sys_in);
	// free(gv.buffer);	
	// free(gv.arg_bulk);
	// free(gv.arg_gamma);
	// free(gv.bulk_rock);

	n_ox 	= gv.len_ox;
	pp 		= gv.len_pp;
	ss 		= gv.len_ss;

	for (j = 0; j < n_ox; j++) {free(gv.ox[j]);} 		free(gv.ox);	
	for (j = 0; j < pp; j++) {	free(gv.PP_list[j]);} 	free(gv.PP_list);
	for (j = 0; j < ss; j++) {	free(gv.SS_list[j]);} 	free(gv.SS_list);
	for (j = 0; j < 2; j++) {	free(gv.pdev[j]);}		free(gv.pdev);
	for (j = 0; j < pp; j++) {	free(gv.pp_flags[j]);}	free(gv.pp_flags);
	for (j = 0; j < n_ox; j++) {free(gv.A[j]);}			free(gv.A);
	for (j = 0; j < n_ox; j++) {free(gv.A2[j]);}		free(gv.A2);

	free(gv.n_SS_PC);
	free(gv.n_min);
	free(gv.verifyPC);
	free(gv.n_ss_ph);
	free(gv.SS_PC_stp);

	free(gv.PGE_mass_norm);
	free(gv.Alg);
	free(gv.gamma_norm);
	free(gv.gibbs_ev);
	free(gv.ite_time);

	free(gv.V_cor);
	free(gv.dGamma);
	free(gv.gam_tot);
	free(gv.gam_tot_0);
	free(gv.delta_gam_tot);
	free(gv.mass_residual);

	free(gv.ipiv);
	free(gv.work);
	free(gv.n_solvi);

	free(gv.pp_n);
	free(gv.pp_n_mol);
	free(gv.pp_xi);
	free(gv.delta_pp_n);
	free(gv.delta_pp_xi);

	free(gv.A_PGE);
	free(gv.A0_PGE);
	free(gv.b_PGE);
	free(gv.cp_id);
	free(gv.pp_id);
	free(gv.dn_cp);
	free(gv.dn_pp);
	free(gv.b);
	free(gv.b1);
	free(gv.pc_id);
	free(gv.tmp1);
	free(gv.tmp2);
	free(gv.tmp3);
	free(gv.n_ss_array);
	/* ================ z_b ============= */
	free(z_b.apo);
	free(z_b.masspo);
	free(z_b.ElEntropy);
	free(z_b.id);
	free(z_b.bulk_rock);
	free(z_b.bulk_rock_cat);
	free(z_b.nzEl_array);
	free(z_b.zEl_array);

	/* ================ final free ============= */
	free(DB.EM_names);
	free(DB.FS_names);
	free(DB.PP_ref_db);
	free(DB.SS_ref_db);
	free(DB.sp);
	free(DB.cp);
}

/** 
  This prints output on screen
**/
void PrintOutput(	global_variable 	gv,
					int 				rank,
					int 				l,
					Databases 			DB,
					double 				time_taken,
					bulk_info 			z_b				){
						
	int i;
	if (gv.verbose !=-1){
		printf(" Status             : %12i ",gv.status);
		if (gv.verbose == 1){PrintStatus(gv.status);}
		printf("\n");
		printf(" Mass residual      : %+12.5e\n",gv.BR_norm);
    	printf(" Rank               : %12i \n",rank);
    	printf(" Point              : %12i \n",l);
    	printf(" Temperature        : %+12.5f\t [C] \n",   z_b.T - 273.15);
		printf(" Pressure           : %+12.5f\t [kbar]\n", z_b.P);
       
 		if (gv.verbose == 1){
 			printf("\n______________________________\n");
			printf("| Comp. Time: %.6f (ms) |\n", time_taken*1000);
			printf("| Min.  Time: %.6f (ms) |", gv.tot_min_time);
            printf("\n══════════════════════════════\n");
        }
    }

	if (gv.verbose != -1){
		printf("\n");
		printf(" SOL = [G: %.3f] (%i iterations, %.2f ms)\n",gv.G_system,gv.global_ite,time_taken*1000.0);
		printf(" GAM = [");
		for (i = 0; i < z_b.nzEl_val-1; i++){
			printf("%+8f,",gv.gam_tot[z_b.nzEl_array[i]]);
		}
		printf("%+8f",gv.gam_tot[z_b.nzEl_val-1]);
		printf("]\n\n");

		printf(" Phase : ");
		for (int i = 0; i < gv.len_cp; i++){
			if (DB.cp[i].ss_flags[1] == 1){
				printf(" %7s ", DB.cp[i].name);
			}
		}
		for (int i = 0; i < gv.len_pp; i++){
			if (gv.pp_flags[i][1] == 1){
				printf(" %7s ", gv.PP_list[i]);
			}
		}
		printf("\n");
		printf(" Mode  : ");
		for (int i = 0; i < gv.len_cp; i++){
			if (DB.cp[i].ss_flags[1] == 1){
				printf(" %.5f ", DB.cp[i].ss_n);
			}
		}
		for (int i = 0; i < gv.len_pp; i++){
			if (gv.pp_flags[i][1] == 1){
				printf(" %.5f ", gv.pp_n[i]);
			}
		}
		printf("\n");	
	}
}

/** 
  This converts the solver status code to human-readable text and prints it to screen
**/
void PrintStatus( int status )
{
	if (status == 0){printf("\t [success]");}
	if (status == 1){printf("\t [success, under-relaxed]");}
	if (status == 2){printf("\t [success, heavily under-relaxed]");}
	if (status == 3){printf("\t [failure, reached maximum iterations]");}
	if (status == 4){printf("\t [failure, terminated due to slow convergence or divergence]");}
}
