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
	- Solution phases spinel (spn), biotite (bi), cordierite (cd), clinopyroxene (cpx), orthopyroxene (opx), epidote (ep), garnet (g), hornblende (hb), ilmenite (ilm), silicate melt (liq), muscovite (mu), olivine (ol), ternary feldspar (pl4T), and aqueous fluid (fl).
                                                    
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
#include "Endmembers_tc-ds634.h"
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

// #define n_em_db 291

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
	if (gv.EM_database == 0){
		gv = get_bulk_metapelite( gv );
	}
	else if (gv.EM_database == 2){
		gv = get_bulk_igneous( gv );
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
		gv = ComputeEquilibrium_Point(		gv.EM_database, 
											input_data[point],
											z_b,											/** bulk rock informations 			*/
											gv,												/** global variables (e.g. Gamma) 	*/
											
										   &splx_data,
											DB.PP_ref_db,									/** pure phase database 			*/
											DB.SS_ref_db,									/** solid solution database 		*/
											DB.cp						);

		/* Perform calculation for a single point 									*/	
		gv = ComputePostProcessing(			gv.EM_database,
											z_b,											/** bulk rock informations 			*/
											gv,												/** global variables (e.g. Gamma) 	*/
											DB.PP_ref_db,									/** pure phase database 			*/
											DB.SS_ref_db,									/** solid solution database 		*/
											DB.cp						);
					
		/* Fill structure holding stable phase equilibrium informations 			*/
		fill_output_struct(					gv,												/** global variables (e.g. Gamma) 	*/
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

	if (gv.output_matlab == 1){
		mergeParallel_matlab(gv);
	}
	/* free memory allocated to solution and pure phases */
	FreeDatabases(gv, DB);

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
	global_variable ComputePostProcessing(			int 				 EM_database,
													bulk_info 	 		 z_b,
													global_variable 	 gv,
													PP_ref  			*PP_ref_db,
													SS_ref  			*SS_ref_db,
													csd_phase_set  		*cp					){
								
	PP_ref PP_db;					
								
	double muC, muN, muE, muW, muNN, muNNN, muNE, muNW, G, muN0, muC0;
	
	double P 			  = z_b.P;					/** PC function uses the z_b structure this is why the Pressure is saved here */
	double T 			  = z_b.T;					/** PC function uses the z_b structure this is why the Pressure is saved here */
	double sum_volume     = 0.0;
	double sum_volume_sol = 0.0;
	double dGdTPP, dGdTMP, dG2dT2, dGdP, dGdP_P0, dG2dP2, dG2dP2_N, dGdP_N;
	double mut, mut_N;
	double phase_isoTbulkModulus_P1;

	double density[gv.len_ox];
	int not_only_liq = 0;
	int ss;

	gv = compute_phase_mol_fraction(	gv,
										PP_ref_db,
										SS_ref_db,
										cp				);


	/** calculate oxygen fugacity: mu_O2 = G0_O2 + RTlog(fO2) */
	/* get O2 Gibbs energy of reference */
	double G0_O = 0.0;
	for (int i = 0; i < gv.len_pp; i++){
		if	(strcmp( gv.PP_list[i], "O2") == 0){
			G0_O = PP_ref_db[i].gbase*PP_ref_db[i].factor;
			break;
		}
	}
	/* get chemical potential of Oxygen (index)*/
	int O_ix = -1;
	for (int i = 0; i < gv.len_ox; i++){
		if	(strcmp( gv.ox[i], "O") == 0){
			O_ix = i;
			break;
		}
	}
	if (O_ix != -1){
		gv.system_fO2 = exp( (gv.gam_tot[O_ix]*2.0-G0_O) / (z_b.R*z_b.T));
	}
	else {
		if (gv.verbose == 1){
			printf("Oxygen fugacity could not be calculated, is O2 endmember included? Is pressure = 0.0?\n");
		}
	}

	/** calculate mass, volume and densities */
	for (int i = 0; i < gv.len_cp; i++){
		if (cp[i].ss_flags[1] == 1){

			if (strcmp( cp[i].name, "liq") != 0){
				not_only_liq = 1;
			}

			ss = cp[i].id;	
			
			for (int k = 0; k < cp[i].n_xeos; k++) {
				SS_ref_db[ss].iguess[k] = cp[i].xeos[k];
			}
										
			/** calculate Molar Mass of solution phase */
			cp[i].mass = 0.0;
			for (int k = 0; k < gv.len_ox; k++){
				cp[i].mass	+= cp[i].ss_comp[k]*z_b.masspo[k];
			}

			/** calculate cp of solution phase */
			cp[i].phase_cp 	 		 = 0.0;
			cp[i].volume 	 		 = 0.0;
			cp[i].phase_expansivity	 = 0.0;
			cp[i].phase_bulkModulus  = 0.0;
			cp[i].phase_shearModulus = 0.0;
			cp[i].phase_entropy  	 = 0.0;
			cp[i].phase_enthalpy 	 = 0.0;
			cp[i].phase_isoTbulkModulus = 0.0;
			cp[i].volume_P0 		 = 0.0;
			cp[i].thetaExp 			 = 0.0;
			phase_isoTbulkModulus_P1 = 0.0;

			for (int j = 0; j < cp[i].n_em; j++){ 
				if (SS_ref_db[ss].z_em[j] == 1.0){
					dG2dT2 					 = (SS_ref_db[ss].mu_array[0][j]-2.0*SS_ref_db[ss].mu_array[6][j]+SS_ref_db[ss].mu_array[1][j])/(gv.gb_T_eps*gv.gb_T_eps);
					dG2dP2 					 = (SS_ref_db[ss].mu_array[4][j]-2.0*SS_ref_db[ss].mu_array[5][j]+SS_ref_db[ss].mu_array[6][j])/(gv.gb_P_eps*gv.gb_P_eps);
					dG2dP2_N				 = (SS_ref_db[ss].mu_array[7][j]-2.0*SS_ref_db[ss].mu_array[4][j]+SS_ref_db[ss].mu_array[5][j])/(gv.gb_P_eps*gv.gb_P_eps);
					dGdTPP 					 = (SS_ref_db[ss].mu_array[2][j]-SS_ref_db[ss].mu_array[3][j])/(2.0*gv.gb_T_eps);
					dGdTMP 					 = (SS_ref_db[ss].mu_array[0][j]-SS_ref_db[ss].mu_array[1][j])/(2.0*gv.gb_T_eps);
					dGdP					 = (SS_ref_db[ss].mu_array[5][j]-SS_ref_db[ss].mu_array[6][j])/(gv.gb_P_eps);
					dGdP_N 					 = (SS_ref_db[ss].mu_array[4][j]-SS_ref_db[ss].mu_array[5][j])/(gv.gb_P_eps);
					dGdP_P0 				 = (SS_ref_db[ss].mu_array[8][j]-SS_ref_db[ss].mu_array[9][j])/(gv.gb_P_eps);
					/* heat capacity 	*/
					cp[i].phase_cp    		+= -T*(dG2dT2)*cp[i].p_em[j];
					
					/* volume 			*/
					cp[i].volume    		+= (dGdP)*cp[i].p_em[j];

					/* volume 			*/
					cp[i].volume_P0    		+= (dGdP_P0)*cp[i].p_em[j];

					/* entropy   		*/
					cp[i].phase_entropy 	+= -(dGdTMP)*cp[i].p_em[j];

					/* expansivity 		*/
					cp[i].phase_expansivity += (1.0/(dGdP)*((dGdTPP-dGdTMP)/(gv.gb_P_eps)))*cp[i].p_em[j];
					
					/* bulk modulus	*/
					cp[i].phase_bulkModulus += -dGdP/( dG2dP2 + pow(((dGdTPP-dGdTMP)/(gv.gb_P_eps)),2.0)/dG2dT2 ) * cp[i].p_em[j];

					/* iso bulk modulus	*/
					cp[i].phase_isoTbulkModulus += -dGdP/( dG2dP2 ) 	* cp[i].p_em[j];
					phase_isoTbulkModulus_P1	+= -dGdP_N/( dG2dP2_N ) * cp[i].p_em[j];
							
					/* shear modulus	*/
					cp[i].phase_shearModulus += SS_ref_db[ss].ElShearMod[j] * cp[i].p_em[j];
				}
			}	

			G = 0.0;
			for (int j = 0; j < gv.len_ox; j++){
				G += cp[i].ss_comp[j]*gv.gam_tot[j];
			}

			/* enthalpy   		*/
			cp[i].phase_enthalpy = cp[i].phase_entropy*T + G;
	
			/** calculate density from volume */
			cp[i].phase_density  = (cp[i].mass*1000.0)/(cp[i].volume*10.0);

			mut 				 = (3.0*cp[i].phase_isoTbulkModulus - 6.0*cp[i].phase_isoTbulkModulus*gv.poisson_ratio) / (2. + 2.0*gv.poisson_ratio)/10.0;
			mut_N 				 = (3.0*phase_isoTbulkModulus_P1 	- 6.0*phase_isoTbulkModulus_P1*gv.poisson_ratio) 	/ (2. + 2.0*gv.poisson_ratio)/10.0;

			cp[i].thetaExp 		 = (mut_N - mut)/gv.gb_P_eps - (cp[i].phase_bulkModulus*cp[i].phase_expansivity)/(cp[i].phase_cp*cp[i].phase_density);

			// printf(" %s: volR: %+10f mut: %+10f, mut_N: %+10f, thetaExp: %+10f\n",cp[i].name,cp[i].volume/cp[i].volume_P0,mut,mut_N,cp[i].thetaExp);

			if (strcmp( cp[i].name, "liq") == 0){
				gv.melt_density   	= cp[i].phase_density;
				gv.melt_fraction  	= cp[i].ss_n_mol;
				gv.melt_bulkModulus = cp[i].phase_bulkModulus/10.0;
			}

			/** get sum of volume*fraction*factor to calculate vol% from mol% */
			sum_volume += cp[i].volume*cp[i].ss_n_mol*cp[i].factor;

			if (strcmp( cp[i].name, "liq") != 0 && strcmp( cp[i].name, "fl") != 0){
				sum_volume_sol 		+= cp[i].volume*cp[i].ss_n_mol*cp[i].factor;
				gv.solid_fraction 	+= cp[i].ss_n_mol;
			}

		}
	}

	for (int i = 0; i < gv.len_pp; i++){
		/* if pure phase is active or on hold (PP cannot be removed from consideration */
		if (gv.pp_flags[i][1] == 1){

			/* calculate phase volume as V = dG/dP */
			PP_db    	 = G_EM_function(EM_database, gv.len_ox,z_b.bulk_rock, z_b.apo, 1., z_b.T, gv.PP_list[i], "equilibrium");
			muC0	 	 = PP_db.gbase;

			PP_db    	 = G_EM_function(EM_database, gv.len_ox,z_b.bulk_rock, z_b.apo, 1. + gv.gb_P_eps, z_b.T, gv.PP_list[i], "equilibrium");
			muN0 	 	 = PP_db.gbase;

			PP_db    	 = G_EM_function(EM_database, gv.len_ox,z_b.bulk_rock, z_b.apo, z_b.P, z_b.T, gv.PP_list[i], "equilibrium");
			muC	 	 	 = PP_db.gbase;

			PP_db    	 = G_EM_function(EM_database, gv.len_ox,z_b.bulk_rock, z_b.apo, z_b.P + gv.gb_P_eps, z_b.T, gv.PP_list[i], "equilibrium");
			muN 	 	 = PP_db.gbase;

			PP_db    	 = G_EM_function(EM_database, gv.len_ox,z_b.bulk_rock, z_b.apo, z_b.P + gv.gb_P_eps*2.0, z_b.T, gv.PP_list[i], "equilibrium");
			muNN 	 	 = PP_db.gbase;

			PP_db    	 = G_EM_function(EM_database, gv.len_ox,z_b.bulk_rock, z_b.apo, z_b.P + gv.gb_P_eps*3.0, z_b.T, gv.PP_list[i], "equilibrium");
			muNNN 	 	 = PP_db.gbase;

			PP_db    	 = G_EM_function(EM_database, gv.len_ox,z_b.bulk_rock, z_b.apo, z_b.P, z_b.T + gv.gb_T_eps, gv.PP_list[i], "equilibrium");
			muE	 		 = PP_db.gbase;

			PP_db    	 = G_EM_function(EM_database, gv.len_ox,z_b.bulk_rock, z_b.apo, z_b.P, z_b.T - gv.gb_T_eps, gv.PP_list[i], "equilibrium");			
			muW	 		 = PP_db.gbase;

			PP_db    	 = G_EM_function(EM_database, gv.len_ox,z_b.bulk_rock, z_b.apo, z_b.P + gv.gb_P_eps, z_b.T + gv.gb_T_eps, gv.PP_list[i], "equilibrium");
			muNE	 	 = PP_db.gbase;

			PP_db    	 = G_EM_function(EM_database, gv.len_ox,z_b.bulk_rock, z_b.apo, z_b.P + gv.gb_P_eps, z_b.T - gv.gb_T_eps, gv.PP_list[i], "equilibrium");			
			muNW	 	 = PP_db.gbase;

			/* Calculate mass per pure phase */
			PP_ref_db[i].mass = 0.0;
			for (int j = 0; j< gv.len_ox; j++){
				PP_ref_db[i].mass += PP_ref_db[i].Comp[j]*z_b.masspo[j];
			}

			dG2dT2 		= (muE-2.0*muC+muW)		/(gv.gb_T_eps*gv.gb_T_eps);
			dG2dP2 		= (muNN-2.0*muN+muC)	/(gv.gb_P_eps*gv.gb_P_eps);
			dG2dP2_N	= (muNNN-2.0*muNN+muN)	/(gv.gb_P_eps*gv.gb_P_eps);
			dGdTPP 		= (muNE-muNW)			/(2.0*gv.gb_T_eps);
			dGdTMP 		= (muE-muW)				/(2.0*gv.gb_T_eps);
			dGdP		= (muN-muC)				/(gv.gb_P_eps);
			dGdP_N		= (muNN-muN)			/(gv.gb_P_eps);
			dGdP_P0 	= (muN0-muC0)			/(gv.gb_P_eps);
			
			/* Calculate volume  per pure phase */
			PP_ref_db[i].volume  	   		= dGdP; 

			/* Calculate volume  per pure phase */
			PP_ref_db[i].volume_P0  	    = dGdP_P0; 
			
			/* Calculate density per pure phase */
			PP_ref_db[i].phase_density 		= (1000.0*PP_ref_db[i].mass)/(PP_ref_db[i].volume*10.0);
			
			/* calculate cp of pure phase */
			PP_ref_db[i].phase_cp 			= -T*(dG2dT2);
			
			/* expansivity 		*/
			PP_ref_db[i].phase_expansivity 	= 1.0/(dGdP)*((dGdTPP-dGdTMP)/(gv.gb_P_eps));
			
			/* entropy 		*/
			PP_ref_db[i].phase_entropy 		= -dGdTMP;
			
			/* enthalpy   		*/
			PP_ref_db[i].phase_enthalpy 	= PP_ref_db[i].phase_entropy*T + PP_ref_db[i].gbase;
	
			/* shear modulus	*/
			PP_ref_db[i].phase_bulkModulus	= -dGdP/( dG2dP2 + pow(((dGdTPP-dGdTMP)/(gv.gb_P_eps)),2.0)/dG2dT2 );
	
			/* shear modulus	*/
			PP_ref_db[i].phase_isoTbulkModulus	= -dGdP/( dG2dP2  );
			phase_isoTbulkModulus_P1			= -dGdP_N/( dG2dP2_N  );
	

			mut 				 = (3.0*PP_ref_db[i].phase_isoTbulkModulus - 6.0*PP_ref_db[i].phase_isoTbulkModulus*gv.poisson_ratio) / (2. + 2.0*gv.poisson_ratio)/10.0;
			mut_N 				 = (3.0*phase_isoTbulkModulus_P1 	- 6.0*phase_isoTbulkModulus_P1*gv.poisson_ratio) 	/ (2. + 2.0*gv.poisson_ratio)/10.0;

			PP_ref_db[i].thetaExp= (mut_N - mut)/gv.gb_P_eps - (PP_ref_db[i].phase_bulkModulus*PP_ref_db[i].phase_expansivity)/(PP_ref_db[i].phase_cp*PP_ref_db[i].phase_density);

			// printf(" %s: volR: %+10f mut: %+10f, mut_N: %+10f, thetaExp: %+10f\n",gv.PP_list[i],PP_ref_db[i].volume/PP_ref_db[i].volume_P0,mut,mut_N,PP_ref_db[i].thetaExp);

			/** get sum of volume*fraction*factor to calculate vol% from mol% */
			sum_volume 			+= PP_ref_db[i].volume*gv.pp_n_mol[i]*PP_ref_db[i].factor;
			sum_volume_sol 		+= PP_ref_db[i].volume*gv.pp_n_mol[i]*PP_ref_db[i].factor;
			gv.solid_fraction 	+= gv.pp_n_mol[i];
		}
	}

	/* calculate the bulk and shear modulus of the aggregate using the Voigt-Reuss-Hill averaging scheme with a weighting factor of 0.5 */
	double s1 = 0.0; double b1 = 0.0;
	double s2 = 0.0; double b2 = 0.0;
	double s1S = 0.0; double b1S = 0.0;
	double s2S = 0.0; double b2S = 0.0;

	for (int i = 0; i < gv.len_cp; i++){
		if (cp[i].ss_flags[1] == 1){
			s1 +=  cp[i].volume*cp[i].ss_n_mol*cp[i].factor/sum_volume *  (cp[i].phase_shearModulus/10.0);
			s2 += (cp[i].volume*cp[i].ss_n_mol*cp[i].factor/sum_volume) / (cp[i].phase_shearModulus/10.0);
			b1 +=  cp[i].volume*cp[i].ss_n_mol*cp[i].factor/sum_volume *  (cp[i].phase_bulkModulus /10.0);
			b2 += (cp[i].volume*cp[i].ss_n_mol*cp[i].factor/sum_volume) / (cp[i].phase_bulkModulus /10.0);
			if (strcmp( cp[i].name, "liq") != 0 && strcmp( cp[i].name, "fl") != 0){
				s1S +=  cp[i].volume*cp[i].ss_n_mol*cp[i].factor/sum_volume_sol *  (cp[i].phase_shearModulus/10.0);
				s2S += (cp[i].volume*cp[i].ss_n_mol*cp[i].factor/sum_volume_sol) / (cp[i].phase_shearModulus/10.0);
				b1S +=  cp[i].volume*cp[i].ss_n_mol*cp[i].factor/sum_volume_sol *  (cp[i].phase_bulkModulus /10.0);
				b2S += (cp[i].volume*cp[i].ss_n_mol*cp[i].factor/sum_volume_sol) / (cp[i].phase_bulkModulus /10.0);
			}

		}
	}
	for (int i = 0; i < gv.len_pp; i++){
		if (gv.pp_flags[i][1] == 1){
			s1 +=  PP_ref_db[i].volume*gv.pp_n_mol[i]*PP_ref_db[i].factor/sum_volume *  (PP_ref_db[i].phase_shearModulus/10.0);
			s2 += (PP_ref_db[i].volume*gv.pp_n_mol[i]*PP_ref_db[i].factor/sum_volume) / (PP_ref_db[i].phase_shearModulus/10.0);
			b1 +=  PP_ref_db[i].volume*gv.pp_n_mol[i]*PP_ref_db[i].factor/sum_volume *  (PP_ref_db[i].phase_bulkModulus /10.0);
			b2 += (PP_ref_db[i].volume*gv.pp_n_mol[i]*PP_ref_db[i].factor/sum_volume) / (PP_ref_db[i].phase_bulkModulus /10.0);

			s1S +=  PP_ref_db[i].volume*gv.pp_n_mol[i]*PP_ref_db[i].factor/sum_volume_sol *  (PP_ref_db[i].phase_shearModulus/10.0);
			s2S += (PP_ref_db[i].volume*gv.pp_n_mol[i]*PP_ref_db[i].factor/sum_volume_sol) / (PP_ref_db[i].phase_shearModulus/10.0);
			b1S +=  PP_ref_db[i].volume*gv.pp_n_mol[i]*PP_ref_db[i].factor/sum_volume_sol *  (PP_ref_db[i].phase_bulkModulus /10.0);
			b2S += (PP_ref_db[i].volume*gv.pp_n_mol[i]*PP_ref_db[i].factor/sum_volume_sol) / (PP_ref_db[i].phase_bulkModulus /10.0);
		}
	}

	// Voight-Reuss-Hill averaging
	gv.system_shearModulus 	= 0.50 * s1 + 0.50 * (1.0/(s2));
	gv.system_bulkModulus  	= 0.50 * b1 + 0.50 * (1.0/(b2));

	gv.solid_shearModulus 	= 0.50 * s1S + 0.50 * (1.0/(s2S));
	gv.solid_bulkModulus  	= 0.50 * b1S + 0.50 * (1.0/(b2S));

	/* calculate density of the system */
	for (int i = 0; i < gv.len_cp; i++){
		if (cp[i].ss_flags[1] == 1){
			gv.system_density += cp[i].phase_density*((cp[i].volume*cp[i].ss_n_mol*cp[i].factor)/sum_volume);
			gv.system_entropy += cp[i].phase_entropy*cp[i].ss_n_mol*cp[i].factor;
			if (strcmp( cp[i].name, "liq") != 0 && strcmp( cp[i].name, "fl") != 0){
				gv.solid_density += cp[i].phase_density*((cp[i].volume*cp[i].ss_n_mol*cp[i].factor)/sum_volume_sol);
			}
		}
	}
	for (int i = 0; i < gv.len_pp; i++){
		if (gv.pp_flags[i][1] == 1){
			gv.system_density += PP_ref_db[i].phase_density*((PP_ref_db[i].volume*gv.pp_n_mol[i]*PP_ref_db[i].factor)/sum_volume);
			gv.system_entropy += PP_ref_db[i].phase_entropy*gv.pp_n_mol[i]*PP_ref_db[i].factor;			
			gv.solid_density  += PP_ref_db[i].phase_density*((PP_ref_db[i].volume*gv.pp_n_mol[i]*PP_ref_db[i].factor)/sum_volume_sol);
		}
	}


	gv.system_volume = sum_volume;
	G = 0.0;
	for (int j = 0; j < gv.len_ox; j++){
		G += z_b.bulk_rock[j]*gv.gam_tot[j];
	}
	gv.system_enthalpy = gv.system_entropy*T + G;

	gv.system_Vp 	= sqrt((gv.system_bulkModulus +4.0/3.0*gv.system_shearModulus)/(gv.system_density/1e3));
	gv.system_Vs 	= sqrt(gv.system_shearModulus/(gv.system_density/1e3));

	gv.solid_Vp 	= sqrt((gv.solid_bulkModulus +4.0/3.0*gv.solid_shearModulus)/(gv.solid_density/1e3));
	gv.solid_Vs 	= sqrt(gv.solid_shearModulus/(gv.solid_density/1e3));

	// gv.V_cor[0] 	= gv.solid_Vp;
	// gv.V_cor[1] 	= gv.solid_Vs;

	// if (gv.calc_seismic_cor == 1){


	// 	gv.solid_Vs 	= anelastic_correction( 0,
	// 											gv.solid_Vs,
	// 											z_b.P,
	// 											z_b.T 		);

	// 	gv.V_cor[0] 	= gv.solid_Vp;
	// 	gv.V_cor[1] 	= gv.solid_Vs;

	// 	gv = wave_melt_correction(  	gv,
	// 									z_b,
	// 									0.1				);

	// }


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

	if (EM_database == 0){			// Igneous database //
		SS_mp_objective_init_function(			SS_objective,
												gv							);
	}
	else if (EM_database == 2){			// Igneous database //
		SS_ig_objective_init_function(			SS_objective,
												gv							);
	}


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


	/****************************************************************************************/
	/**                                   LEVELLING                                        **/
	/****************************************************************************************/	
	gv = Levelling(			z_b,										/** bulk rock informations 			*/
							gv,											/** global variables (e.g. Gamma) 	*/

							SS_objective,
							splx_data,
							PP_ref_db,									/** pure phase database 			*/
							SS_ref_db,									/** solution phase database 		*/
							cp							);


	/****************************************************************************************/
	/**                            PARTITIONING GIBBS ENERGY                               **/
	/****************************************************************************************/
	
	if (gv.solver == 0){ 		/* Legacy solver only */

		gv.div 		= 0;
		gv.status 	= 0;

		gv = init_LP(			z_b,
								splx_data,
								gv,
										
								PP_ref_db,
								SS_ref_db,
								cp						);	

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
		// if (gv.div == 1  && gv.solver == 1){	
			printf("\n[PGE failed -> legacy solver...]\n");
			gv.div 		= 0;
			gv.status 	= 0;

			gv = init_LP(			z_b,
									splx_data,
									gv,
											
									PP_ref_db,
									SS_ref_db,
									cp						);	

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
		printf("╔════════════════════════════════════════════════╗\n");
		printf("║               COMPUTATION SUMMARY              ║\n");
		printf("╚════════════════════════════════════════════════╝\n\n");
		printf(" Alg | ite  | duration   |  MASS norm | Gamma norm\n");
		printf("══════════════════════════════════════════════════\n");

		for (int i = 0; i < gv.global_ite; i++){	
			if (gv.Alg[i] == 0){
				printf(" LP  | %4d | %+10f | %+10f | %+10f\n",i,gv.ite_time[i],gv.PGE_mass_norm[i],gv.gamma_norm[i]);
			}
			if (gv.Alg[i] == 1){
				printf(" PGE | %4d | %+10f | %+10f | %+10f\n",i,gv.ite_time[i],gv.PGE_mass_norm[i],gv.gamma_norm[i]);
			}	
			if (gv.Alg[i+1] - gv.Alg[i] == 1){
				printf("--------------------------------------------------\n");
				printf("               SWITCH FROM LP TO PGE              \n");
				printf("--------------------------------------------------\n");
			}
			if (gv.Alg[i+1] - gv.Alg[i] == -1 && i < gv.global_ite - 1){
				printf("--------------------------------------------------\n");
				printf("               SWITCH FROM PGE TO LP              \n");
				printf("--------------------------------------------------\n");
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
		
    	{ NULL, 0, 0 }
	};
	ketopt_t opt = KETOPT_INIT;

	int    c;
	while ((c = ketopt(&opt, argc, argv, 1, "", longopts)) >= 0) {
		if 		(c == 314){ printf("MAGEMin %20s\n",gv.version ); exit(0); }	
		else if (c == 315){ print_help( gv ); 					  exit(0); }	
        else if	(c == 301){ gv.verbose  = atoi(opt.arg					); }
		else if (c == 316){ gv.solver   = atoi(opt.arg);		if (gv.verbose == 1){		printf("--solver      : solver               = %i \n", 	 	   		gv.solver			);}}																		
		else if (c == 318){ gv.output_matlab   = atoi(opt.arg); if (gv.verbose == 1){		printf("--out_matlab  : out_matlab           = %i \n", 	 	   		gv.output_matlab	);}}																		
		else if (c == 303){ strcpy(gv.File,opt.arg);		 	if (gv.verbose == 1){		printf("--File        : File                 = %s \n", 	 	   		gv.File				);}}
		else if (c == 302){ strcpy(gv.db,opt.arg);		 		if (gv.verbose == 1){		printf("--db          : db                   = %s \n", 	 	   		gv.db				);}}
		else if (c == 317){ strcpy(gv.sys_in,opt.arg);		 	if (gv.verbose == 1){		printf("--sys_in      : sys_in               = %s \n", 	 	   		gv.sys_in			);}}
		else if (c == 304){ gv.n_points = atoi(opt.arg); 	 	if (gv.verbose == 1){		printf("--n_points    : n_points             = %i \n", 	 	   		gv.n_points			);}}
		else if (c == 305){ gv.test  	= atoi(opt.arg); 		if (gv.verbose == 1){		printf("--test        : Test                 = %i \n", 	 	  		gv.test				);}}
		else if (c == 306){ z_b->T=strtof(opt.arg,NULL)+273.15; if (gv.verbose == 1){		printf("--Temp        : Temperature          = %f C \n",            z_b->T-273.15		);}}
		else if (c == 307){ z_b->P = strtof(opt.arg,NULL); 		if (gv.verbose == 1){		printf("--Pres        : Pressure             = %f kbar \n", 		z_b->P				);}}
		else if (c == 308){ strcpy(gv.Phase,opt.arg);		 	if (gv.verbose == 1){		printf("--Phase       : Phase name           = %s \n", 	   			gv.Phase			);}}
		else if (c == 313){ gv.maxeval  = strtof(opt.arg,NULL); if (gv.verbose == 1){
            if (gv.maxeval==0){     printf("--maxeval     : Max. # of local iter.    = infinite  \n"		); }
            else{                   printf("--maxeval     : Max. # of local iter.    = %i  \n", gv.maxeval		);}
            }
        }
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
	if (strcmp(gv.db, "mp") == 0){
		gv.EM_database = 0;
	}
	else if (strcmp(gv.db, "ig") == 0){
		gv.EM_database = 2;
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
	
	/* Allocate solid-solution reference database (to get gbase, comp and factor) 	*/
	DB.SS_ref_db = malloc ((gv.len_ss) 		* sizeof(SS_ref)); 
	
	/* Allocate memory of the considered set of phases 								*/
	DB.cp 		 = malloc ((gv.max_n_cp) 	* sizeof(csd_phase_set)); 
	
	/* Allocate memory of the considered set of phases 								*/
	DB.sp 		 = malloc (1 	* sizeof(stb_system)); 

	/** 
		Allocate memory for each solution phase according to their specificities (n_em, sf etc) 
	*/
	for (i = 0; i < gv.len_ss; i++){
		DB.SS_ref_db[i] = G_SS_init_EM_function(		i,	
														DB.SS_ref_db[i], 
														EM_database, 
														gv.SS_list[i], 
														gv						);
	}

	/* Allocate memory of the considered set of phases 								*/
	for (i = 0; i < gv.max_n_cp; i++){
		DB.cp[i] = CP_INIT_function(		DB.cp[i], 
											gv									);
	}
	
	/* Allocate memory for the stable phase equilibrium 							*/
	DB.sp[0] 	 = SP_INIT_function(		DB.sp[0], gv						);

	/* Endmember names */
	DB.EM_names  =	get_EM_DB_names(		gv									);

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

	return DB;
}

/** 
  Free the memory associated with the databases
**/
void FreeDatabases(		global_variable gv, 
						Databases 		DB	){

	int n_em_db = gv.n_em_db;
	for (int i = 0; i < n_em_db; i++) {
		free(DB.EM_names[i]);
	}
	free(DB.EM_names);
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
