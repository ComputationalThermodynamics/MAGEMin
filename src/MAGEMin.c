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
#include "Endmembers_tc-ds634.h"
#include "toolkit.h"
#include "io_function.h"
#include "gem_function.h"
#include "gss_init_function.h"
#include "gss_function.h"
#include "objective_functions.h"
#include "NLopt_opt_function.h"
#include "simplex_levelling.h"
#include "Initialize.h"
#include "ss_min_function.h"
#include "pp_min_function.h"
#include "dump_function.h"
#include "PGE_function.h"
#include "phase_update_function.h"
#include "MAGEMin.h"
#include "simplex_levelling.h"

#define n_em_db 291

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
	int i,j, k;
	int rank, numprocs;
	int EM_database;
	double time_taken;

	Databases DB;
	
	/* Select the endmember database */
   	EM_database = _tc_ds634_;
   	
	clock_t t,u; 
	t = clock();
	u = clock();
	
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	/* initialize global structure to store shared variables (e.g. Gamma, SS and PP list, ...) */
	global_variable gv;
	gv = global_variable_init();

	/* Allocate both pure and solid-solution databases */
	DB = InitializeDatabases(gv, EM_database);
	
	/** Declare bulk info structure */
	bulk_info z_b;

	/* Default conditions */
	double *bulk_rock = malloc ((gv.len_ox) * sizeof (double) ); 
	double 	Pres;
	double 	Temp;
	
	double 	P 			=  0.0;
	double 	T 			=  0.0;
	
	int 	test 		=  0;
	int 	Verb 		= -1;
	int     Mode 		=  0;
	int     n_points 	=  1;
	
	int 	maxeval		= -1;
	
	int		get_version;
	int		get_help;
	
	double	Gam[11],  Bulk[11], InitEM_Prop[15];
	char    File[50], Phase[50];

	/** 
	Read command-line arguments and set default parameters
	*/
	gv = ReadCommandLineOptions(	 gv,
									 argc, 
									 argv,  
									&Mode, 
									&Verb, 
									&test, 
									&n_points, 
									&Pres, 
									&Temp,
									 Bulk, 
									 Gam, 
									 File, 
									 Phase, 
									&maxeval,
									&get_version,
									&get_help		); 
									

	gv.verbose 	= Verb;
	gv.Mode 	= Mode;

    if (maxeval>-1){
        gv.maxeval = maxeval;   // otherwise we use default. Note that 0 = no limit
    }

	/* initial dumping logs and output */
	dump_init(gv);

	if (rank==0 && gv.verbose != -1){
    	printf("\nRunning MAGEMin %5s on %d cores {\n", gv.version, numprocs);
    	printf("═══════════════════════════════════\n");
	}

	/****************************************************************************************/
	/**              READ INPUT FILE FOR MULTIPLE P-T CONDITIONS, IF IT EXISTS             **/
	/****************************************************************************************/
	/* 	allocate space to store input data, P,T,X, gam, x-eos etc. */
	io_data input_data[n_points]; 

	/* get data from input file */
	if (strcmp( File, "none") != 0){	
		read_in_data(gv, input_data, File, n_points);			
	}
	
	/****************************************************************************************/
	/**                       DEFINE SOME TEST BULK-ROCK COMPOSITIONS                      **/
	/****************************************************************************************/
	/* get bulk rock composition parsed from args 							*/
	get_bulk(								bulk_rock,
											test,
											gv.len_ox 					);

	/* Override P,T & bulk with command-line options, if there is no file parsed: */
	if (Pres    > 0.0){ P = Pres 										;}
	if (Temp    > 0.0){ T = Temp + 273.15								;}
	
	if (Bulk[0] > 0.0) {
		for (i = 0; i < gv.len_ox; i++){ bulk_rock[i] = Bulk[i];}
	}
	
	/** Normalize composition to sum to 1. 									*/
	norm_array(								bulk_rock,
											gv.len_ox					);						
								
	/** Get zeros in bulk P and T 											*/				
	z_b = initialize_bulk_infos(			P, 
											T							);	

	/** allocate simplex data memory outside the MPI loop 					*/
	simplex_data 							splx_data;

	init_simplex_A(			   		   	   &splx_data,
											gv							);
										
	init_simplex_B_em(				   	   &splx_data,
											gv							);

	/** pointer array to objective functions 								*/
	obj_type 								SS_objective[gv.len_ss];	
	
	SS_objective_init_function(				SS_objective,
											gv							);
								
		
	/****************************************************************************************/
	/**                               LAUNCH MINIMIZATION ROUTINE                          **/
	/****************************************************************************************/
	for (int sgleP = 0; sgleP < n_points; sgleP++){
        if ((Mode==0) && (sgleP % numprocs != rank)) continue;   	/** this ensures that, in parallel, not every point is computed by every processor (instead only every numprocs point). Only applied to Mode==0 */

		t              = clock();									/** reset loop timer 				*/
		gv.numPoint    = sgleP; 									/** the number of the current point */

		/* If we read input from file: */
		if (strcmp( File, "none") != 0){						
			z_b.P = input_data[sgleP].P;
			z_b.T = input_data[sgleP].T + 273.15;					/** K to C 		*/

			for (int i = 0; i < gv.len_ox; i++){
				gv.gam_tot[i] = input_data[sgleP].in_gam[i];					
			}
			// if (input_data[sgleP].bulk_rock[i] > 0.0){
			// 	for (int i = 0; i < gv.len_ox; i++){
			// 		bulk_rock[i] = input_data[sgleP].bulk_rock[i];					
			// 	}	
			// }
		}
		

		/* reset global variables flags 											*/
		gv = reset_gv(						gv,
											z_b,
											DB.PP_ref_db,
											DB.SS_ref_db				);

		/** reset bulk rock information (needed for parallel point calculation) 	*/
		z_b = reset_z_b(					gv,				
											bulk_rock,								
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
		gv = ComputeEquilibrium_Point(		EM_database, 
											input_data[sgleP],
											Mode,
											z_b,											/** bulk rock informations 			*/
											gv,												/** global variables (e.g. Gamma) 	*/

											SS_objective,
										   &splx_data,
											DB.PP_ref_db,									/** pure phase database 			*/
											DB.SS_ref_db,									/** solid solution database 		*/
											DB.cp						);

		/* Perform calculation for a single point 									*/	
		gv = ComputePostProcessing(			EM_database,
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

		/* Dump final results to files 												*/
		dump_results_function(				gv,												/** global variables (e.g. Gamma) 	*/
											z_b,											/** bulk-rock informations 			*/
											DB.PP_ref_db,									/** pure phase database 			*/
											DB.SS_ref_db,									/** solution phase database 		*/
											DB.cp						);


		/* Print output to screen 													*/
		t 			= clock() - t; 
		time_taken 	= ((double)t)/CLOCKS_PER_SEC; 											/* in seconds 	 					*/
		PrintOutput(gv, rank, sgleP, DB, time_taken, z_b);									/* print output on screen 			*/
	}
	/* end of loop over points */

	// 2BM
	// for (int i = 0; i < splx_data.n_Ox; i++){
	// printf("%+10f ",splx_data.gamma_tot[i]);
	// printf("\n");
	// }

	/* wait for all cores to be finished */
	MPI_Barrier(MPI_COMM_WORLD);		

	/* now merge the parallel output files into one*/
	mergeParallelFiles(gv);

	if (gv.save_residual_evolution == 1){
		mergeParallel_residual_Files(gv);
	}

	/* free memory allocated to solution and pure phases */
	FreeDatabases(gv, DB);

	/** deallocate memory */
	destroy_simplex_A(&splx_data);
	destroy_simplex_B(&splx_data);

	/* print the time */
	u = clock() - u; 
	
	if (gv.verbose != -1){
		time_taken = ((double)u)/(CLOCKS_PER_SEC); 				/** in seconds */
		if (rank==0){
			printf("__________________________________\n");
			printf("MAGEMin comp time: %+3f ms }\n", time_taken*1000.);
		}
	}
    return 0;
}

/** 
  Compute stable equilibrium at given P/T/C point
*/
	global_variable ComputePostProcessing(			int 				 EM_database,
													bulk_info 	 z_b,
													global_variable 	 gv,
													PP_ref  			*PP_ref_db,
													SS_ref  			*SS_ref_db,
													csd_phase_set  		*cp					){
								
	PP_ref PP_db;					
								
	double muC, muN, muE, muW, muNN, muNE, muNW;
	
	double P 			  = z_b.P;					/** PC function uses the z_b structure this is why the Pressure is saved here */
	double T 			  = z_b.T;					/** PC function uses the z_b structure this is why the Pressure is saved here */
	double sum_volume     = 0.0;
	double sum_volume_sol = 0.0;
	double dGdTPP, dGdTMP, dG2dT2, dGdP, dG2dP2;

	double density[gv.len_ox];

	int ss;
	/** calculate mass, volume and densities */
	for (int i = 0; i < gv.len_cp; i++){
		if (cp[i].ss_flags[1] == 1){
			ss = cp[i].id;	
			
			for (int k = 0; k < cp[i].n_xeos; k++) {
				SS_ref_db[ss].iguess[k] = cp[i].xeos[k];
			}
			
			/* Associate the right solid-solution data */
			for (int FD = 0; FD < 7; FD++){

				z_b.P 		= P + gv.gb_P_eps*gv.numDiff[0][FD];
				z_b.T 		= T + gv.gb_T_eps*gv.numDiff[1][FD];
						
				SS_ref_db[ss] = raw_hyperplane(		gv, 
													SS_ref_db[ss],
													SS_ref_db[ss].mu_array[FD]	);
				
				SS_ref_db[ss] = PC_function(		gv,
													SS_ref_db[ss], 
													z_b,
													gv.SS_list[ss] 				);
													
				for (int j = 0; j < SS_ref_db[ss].n_em; j++){ 
					SS_ref_db[ss].mu_array[FD][j] = SS_ref_db[ss].mu[j];
				}
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

			for (int j = 0; j < cp[i].n_em; j++){ 
				if (SS_ref_db[ss].z_em[j] == 1.0){
					dG2dT2 					 = (SS_ref_db[ss].mu_array[0][j]-2.0*SS_ref_db[ss].mu_array[6][j]+SS_ref_db[ss].mu_array[1][j])/(gv.gb_T_eps*gv.gb_T_eps);
					dG2dP2 					 = (SS_ref_db[ss].mu_array[4][j]-2.0*SS_ref_db[ss].mu_array[5][j]+SS_ref_db[ss].mu_array[6][j])/(gv.gb_P_eps*gv.gb_P_eps);
					dGdTPP 					 = (SS_ref_db[ss].mu_array[2][j]-SS_ref_db[ss].mu_array[3][j])/(2.0*gv.gb_T_eps);
					dGdTMP 					 = (SS_ref_db[ss].mu_array[0][j]-SS_ref_db[ss].mu_array[1][j])/(2.0*gv.gb_T_eps);
					dGdP					 = (SS_ref_db[ss].mu_array[5][j]-SS_ref_db[ss].mu_array[6][j])/(gv.gb_P_eps);

					/* heat capacity 	*/
					cp[i].phase_cp    		+= -T*(dG2dT2)*cp[i].p_em[j];
					
					/* volume 			*/
					cp[i].volume    		+= (dGdP)*cp[i].p_em[j];
					
					/* expansivity 		*/
					cp[i].phase_expansivity += (1.0/(dGdP)*((dGdTPP-dGdTMP)/(gv.gb_P_eps)))*cp[i].p_em[j];
					
					/* bulk modulus	*/
					cp[i].phase_bulkModulus += -dGdP/( dG2dP2 + pow(((dGdTPP-dGdTMP)/(gv.gb_P_eps)),2.0)/dG2dT2 ) * cp[i].p_em[j];
					
					/* shear modulus	*/
					cp[i].phase_shearModulus += SS_ref_db[ss].ElShearMod[j] * cp[i].p_em[j];
				}
			}	
			
			/** calculate density from volume */
			cp[i].phase_density = (cp[i].mass*1000.0)/(cp[i].volume*10.0);

			if (strcmp( cp[i].name, "liq") == 0){
				gv.melt_density   	= cp[i].phase_density;
				gv.melt_fraction  	= cp[i].ss_n;
				gv.melt_bulkModulus = cp[i].phase_bulkModulus/10.0;
			}

			/** get sum of volume*fraction*factor to calculate vol% from mol% */
			sum_volume += cp[i].volume*cp[i].ss_n*cp[i].factor;

			if (strcmp( cp[i].name, "liq") != 0 && strcmp( cp[i].name, "fl") != 0){
				sum_volume_sol += cp[i].volume*cp[i].ss_n*cp[i].factor;
			}

		}
	}

	for (int i = 0; i < gv.len_pp; i++){

		/* if pure phase is active or on hold (PP cannot be removed from consideration */
		if (gv.pp_flags[i][1] == 1){

			/* calculate phase volume as V = dG/dP */
			PP_db    	 = G_EM_function(EM_database, z_b.bulk_rock, z_b.P, z_b.T, gv.PP_list[i], "equilibrium");
			muC	 	 	 = PP_db.gbase;
			PP_db    	 = G_EM_function(EM_database, z_b.bulk_rock, z_b.P + gv.gb_P_eps, z_b.T, gv.PP_list[i], "equilibrium");
			muN 	 	 = PP_db.gbase;
			PP_db    	 = G_EM_function(EM_database, z_b.bulk_rock, z_b.P + gv.gb_P_eps*2.0, z_b.T, gv.PP_list[i], "equilibrium");
			muNN 	 	 = PP_db.gbase;
		
			PP_db    	 = G_EM_function(EM_database, z_b.bulk_rock, z_b.P, z_b.T + gv.gb_T_eps, gv.PP_list[i], "equilibrium");
			muE	 		 = PP_db.gbase;
			PP_db    	 = G_EM_function(EM_database, z_b.bulk_rock, z_b.P, z_b.T - gv.gb_T_eps, gv.PP_list[i], "equilibrium");			
			muW	 		 = PP_db.gbase;

			PP_db    	 = G_EM_function(EM_database, z_b.bulk_rock, z_b.P + gv.gb_P_eps, z_b.T + gv.gb_T_eps, gv.PP_list[i], "equilibrium");
			muNE	 	 = PP_db.gbase;
			PP_db    	 = G_EM_function(EM_database, z_b.bulk_rock, z_b.P + gv.gb_P_eps, z_b.T - gv.gb_T_eps, gv.PP_list[i], "equilibrium");			
			muNW	 	 = PP_db.gbase;

			/* Calculate mass per pure phase */
			PP_ref_db[i].mass = 0.0;
			for (int j = 0; j< nEl; j++){
				PP_ref_db[i].mass += PP_ref_db[i].Comp[j]*z_b.masspo[j];
			}

			dG2dT2 		= (muE-2.0*muC+muW)/(gv.gb_T_eps*gv.gb_T_eps);
			dG2dP2 		= (muNN-2.0*muN+muC)/(gv.gb_P_eps*gv.gb_P_eps);
			dGdTPP 		= (muNE-muNW)/(2.0*gv.gb_T_eps);
			dGdTMP 		= (muE-muW)/(2.0*gv.gb_T_eps);
			dGdP		= (muN-muC)/(gv.gb_P_eps);

			/* Calculate volume  per pure phase */
			PP_ref_db[i].volume  	   		= dGdP; 
			
			/* Calculate density per pure phase */
			PP_ref_db[i].phase_density 		= (1000.0*PP_ref_db[i].mass)/(PP_ref_db[i].volume*10.0);
			
			/* calculate cp of pure phase */
			PP_ref_db[i].phase_cp 			= -T*(dG2dT2);
			
			/* expansivity 		*/
			PP_ref_db[i].phase_expansivity 	= 1.0/(dGdP)*((dGdTPP-dGdTMP)/(gv.gb_P_eps));
			
			/* shear modulus	*/
			PP_ref_db[i].phase_bulkModulus	= -dGdP/( dG2dP2 + pow(((dGdTPP-dGdTMP)/(gv.gb_P_eps)),2.0)/dG2dT2 );
	
			/** get sum of volume*fraction*factor to calculate vol% from mol% */
			sum_volume 		+= PP_ref_db[i].volume*gv.pp_n[i]*PP_ref_db[i].factor;
			sum_volume_sol 	+= PP_ref_db[i].volume*gv.pp_n[i]*PP_ref_db[i].factor;
		}
	}

	/* calculate the bulk and shear modulus of the aggregate using the Voigt-Reuss-Hill averaging scheme with a weighting factor of 0.5 */
	double s1 = 0.0; double b1 = 0.0;
	double s2 = 0.0; double b2 = 0.0;
	double s1S = 0.0; double b1S = 0.0;
	double s2S = 0.0; double b2S = 0.0;

	for (int i = 0; i < gv.len_cp; i++){
		if (cp[i].ss_flags[1] == 1){
			s1 +=  cp[i].volume*cp[i].ss_n*cp[i].factor/sum_volume *  (cp[i].phase_shearModulus/10.0);
			s2 += (cp[i].volume*cp[i].ss_n*cp[i].factor/sum_volume) / (cp[i].phase_shearModulus/10.0);
			b1 +=  cp[i].volume*cp[i].ss_n*cp[i].factor/sum_volume *  (cp[i].phase_bulkModulus /10.0);
			b2 += (cp[i].volume*cp[i].ss_n*cp[i].factor/sum_volume) / (cp[i].phase_bulkModulus /10.0);
			if (strcmp( cp[i].name, "liq") != 0 && strcmp( cp[i].name, "fl") != 0){
				s1S +=  cp[i].volume*cp[i].ss_n*cp[i].factor/sum_volume_sol *  (cp[i].phase_shearModulus/10.0);
				s2S += (cp[i].volume*cp[i].ss_n*cp[i].factor/sum_volume_sol) / (cp[i].phase_shearModulus/10.0);
				b1S +=  cp[i].volume*cp[i].ss_n*cp[i].factor/sum_volume_sol *  (cp[i].phase_bulkModulus /10.0);
				b2S += (cp[i].volume*cp[i].ss_n*cp[i].factor/sum_volume_sol) / (cp[i].phase_bulkModulus /10.0);
			}

		}
	}
	for (int i = 0; i < gv.len_pp; i++){
		if (gv.pp_flags[i][1] == 1){
			s1 +=  PP_ref_db[i].volume*gv.pp_n[i]*PP_ref_db[i].factor/sum_volume *  (PP_ref_db[i].phase_shearModulus/10.0);
			s2 += (PP_ref_db[i].volume*gv.pp_n[i]*PP_ref_db[i].factor/sum_volume) / (PP_ref_db[i].phase_shearModulus/10.0);
			b1 +=  PP_ref_db[i].volume*gv.pp_n[i]*PP_ref_db[i].factor/sum_volume *  (PP_ref_db[i].phase_bulkModulus /10.0);
			b2 += (PP_ref_db[i].volume*gv.pp_n[i]*PP_ref_db[i].factor/sum_volume) / (PP_ref_db[i].phase_bulkModulus /10.0);

			s1S +=  PP_ref_db[i].volume*gv.pp_n[i]*PP_ref_db[i].factor/sum_volume_sol *  (PP_ref_db[i].phase_shearModulus/10.0);
			s2S += (PP_ref_db[i].volume*gv.pp_n[i]*PP_ref_db[i].factor/sum_volume_sol) / (PP_ref_db[i].phase_shearModulus/10.0);
			b1S +=  PP_ref_db[i].volume*gv.pp_n[i]*PP_ref_db[i].factor/sum_volume_sol *  (PP_ref_db[i].phase_bulkModulus /10.0);
			b2S += (PP_ref_db[i].volume*gv.pp_n[i]*PP_ref_db[i].factor/sum_volume_sol) / (PP_ref_db[i].phase_bulkModulus /10.0);
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
			gv.system_density += cp[i].phase_density*((cp[i].volume*cp[i].ss_n*cp[i].factor)/sum_volume);
			if (strcmp( cp[i].name, "liq") != 0 && strcmp( cp[i].name, "fl") != 0){
				gv.solid_density += cp[i].phase_density*((cp[i].volume*cp[i].ss_n*cp[i].factor)/sum_volume_sol);
			}
		}
	}
	for (int i = 0; i < gv.len_pp; i++){
		if (gv.pp_flags[i][1] == 1){
			gv.system_density += PP_ref_db[i].phase_density*((PP_ref_db[i].volume*gv.pp_n[i]*PP_ref_db[i].factor)/sum_volume);
			gv.solid_density  += PP_ref_db[i].phase_density*((PP_ref_db[i].volume*gv.pp_n[i]*PP_ref_db[i].factor)/sum_volume_sol);
		}
	}

	gv.system_Vp 	= sqrt((gv.system_bulkModulus +4.0/3.0*gv.system_shearModulus)/(gv.system_density/1e3));
	gv.system_Vs 	= sqrt(gv.system_shearModulus/(gv.system_density/1e3));

	gv.solid_Vp 	= sqrt((gv.solid_bulkModulus +4.0/3.0*gv.solid_shearModulus)/(gv.solid_density/1e3));
	gv.solid_Vs 	= sqrt(gv.solid_shearModulus/(gv.solid_density/1e3));

	if (gv.verbose != -1){
		printf("\nSystem information\n");
		printf("═══════════════════\n");

		printf(" Mass residual      : %+12.5e\n\n",gv.BR_norm);

		printf(" Volume             : %+12.5f\n",sum_volume);
		printf(" Density            : %+12.5f\n\n",gv.system_density);

		printf(" Shear modulus      : %+12.5f\t [GPa]\n",gv.system_shearModulus);
		printf(" Vp           (VRH) : %+12.5f\t [km/s]\n",gv.system_Vp);
		printf(" Vs           (VRH) : %+12.5f\t [km/s]\n",gv.system_Vs);
		printf(" Vp/Vs        (VRH) : %+12.5f\t [km/s]\n\n",gv.system_Vp/gv.system_Vs);
	}

	return gv;
}

/** 
  Compute stable equilibrium at given P/T/C point
*/
global_variable ComputeEquilibrium_Point( 		int 				 EM_database,
												io_data 			 input_data,
												int 				 Mode,
												bulk_info 	 z_b,
												global_variable 	 gv,

												obj_type 			*SS_objective,
												simplex_data	    *splx_data,
												PP_ref  			*PP_ref_db,
												SS_ref  			*SS_ref_db,
												csd_phase_set  		*cp						){

	/* initialize endmember database for given P-T point */
	gv = init_em_db(		EM_database,
							z_b,										/** bulk rock informations */
							gv,											/** global variables (e.g. Gamma) */
							PP_ref_db				);

	/* Calculate solution phase data at given P-T conditions (G0 based on G0 of endmembers) */
	gv = init_ss_db(		EM_database,
							z_b,
							gv,
							SS_ref_db						);


	/* if Mode is 0, perform normal minimization with PGE */
	if (Mode == 0){
		/****************************************************************************************/
		/**                                   LEVELLING                                        **/
		/****************************************************************************************/	
		gv = Levelling(			z_b,									/** bulk rock informations */
								gv,										/** global variables (e.g. Gamma) */

								SS_objective,
							    splx_data,
								PP_ref_db,								/** pure phase database */
								SS_ref_db,								/** solution phase database */
								cp					);
		
		/****************************************************************************************/
		/**                                   MAIN LOOP (PGE)                                  **/
		/****************************************************************************************/
		gv = PGE(				z_b,									/** bulk rock constraint */ 
								gv,										/** global variables (e.g. Gamma) */

								SS_objective,
							    splx_data,
								PP_ref_db,								/** pure phase database */
								SS_ref_db,								/** solution phase database */
								cp					);

		if (0 == 1){
			for (int i = 0; i < gv.len_ss; i++){
				printf("#PGE_PC %s: %d\n", gv.SS_list[i], SS_ref_db[i].tot_Ppc);
				for (int k = 0; k < SS_ref_db[i].tot_Ppc; k++){
					for (int j = 0; j < gv.len_ox; j++){
						printf(" %g",SS_ref_db[i].comp_Ppc[k][j]);
					}
					printf("\n");
				}
			}
		}

	}
	/* if Mode = 1, spit out Gibbs energy and reference values with given compositional variables */
	else if (Mode == 1){
		printf("\n");
		printf("  Spit out Solution model informations for given input\n");
		printf("  ════════════════════════════════════════════════════\n");
		int id_cp = 0;
		for (int i = 0; i < input_data.n_phase; i++){
			/* simple function to get the array index of the given solution phase */
			int ss = get_phase_id(		gv,
										input_data.phase_names[i]	);
			
			for (int j = 0; j < SS_ref_db[ss].n_xeos; j++){
				SS_ref_db[ss].iguess[j] =  input_data.phase_xeos[i][j];
			}
			SS_ref_db[ss] = raw_hyperplane(		gv, 
												SS_ref_db[ss],
												SS_ref_db[ss].gbase		);
			
			SS_ref_db[ss] = PC_function(	gv,
											SS_ref_db[ss], 
											z_b,
											gv.SS_list[ss] 				);
											
			strcpy(cp[id_cp].name,gv.SS_list[ss]);				/* get phase name */	
			
			cp[id_cp].split 		= 0;							
			cp[id_cp].id 			= ss;						/* get phase id */
			cp[id_cp].n_xeos		= SS_ref_db[ss].n_xeos;		/* get number of compositional variables */
			cp[id_cp].n_em			= SS_ref_db[ss].n_em;		/* get number of endmembers */
			cp[id_cp].n_sf			= SS_ref_db[ss].n_sf;		/* get number of site fractions */
			
			cp[id_cp].df			= SS_ref_db[ss].df_raw;
			cp[id_cp].factor		= SS_ref_db[ss].factor;	
			
			cp[id_cp].ss_flags[0] 	= 1;							/* set flags */
			cp[id_cp].ss_flags[1] 	= 1;
			cp[id_cp].ss_flags[2] 	= 0;
			
			cp[id_cp].ss_n          = 1.0;			/* get initial phase fraction */
			
			for (int ii = 0; ii < cp[id_cp].n_xeos; ii++){
				cp[id_cp].xeos[ii]		= SS_ref_db[ss].iguess[ii]; 
				cp[id_cp].dfx[ii]		= SS_ref_db[ss].dfx[ii]; 
			}
			
			for (int ii = 0; ii < cp[id_cp].n_em; ii++){
				cp[id_cp].p_em[ii]		= SS_ref_db[ss].p[ii];
				cp[id_cp].xi_em[ii]		= SS_ref_db[ss].xi_em[ii];
				cp[id_cp].mu[ii]		= SS_ref_db[ss].mu[ii];
				cp[id_cp].gbase[ii]		= SS_ref_db[ss].gbase[ii];
			}
			for (int ii = 0; ii < SS_ref_db[ss].n_em; ii++){
				for (int jj = 0; jj < SS_ref_db[ss].n_xeos; jj++){
					cp[id_cp].dpdx[ii][jj] = SS_ref_db[ss].dp_dx[ii][jj];
				}
			}
			for (int ii = 0; ii < gv.len_ox; ii++){
				cp[id_cp].ss_comp[ii]	= SS_ref_db[ss].ss_comp[ii];
			}
			for (int ii = 0; ii < cp[id_cp].n_sf; ii++){
				cp[id_cp].sf[ii]		= SS_ref_db[ss].sf[ii];
			}	
			
			gv.id_solvi[ss][gv.n_solvi[ss]] = id_cp;
			gv.n_solvi[ss] 	   	   += 1;
			id_cp 				   += 1;
			gv.len_cp 			   += 1;
			gv.n_cp_phase 		   += 1;
			gv.n_phase             += 1;
		
			if (gv.verbose ==1){
				printf("   -> reading in %4s %+10f|",gv.SS_list[ss],SS_ref_db[ss].df);
				for (int j = 0; j < SS_ref_db[ss].n_xeos; j++){
					printf(" %+12.5f", input_data.phase_xeos[i][j]);
				}
				printf("\n");
			}
			
		}
	}
	/* if Mode = 2, perform search of local minima for given solution phase */
	else if (Mode == 2){
		printf("function has been deleted\n");
	}
	/* if Mode = 3, perform first stage levelling only */
	else if (Mode == 3){
		/* when Mode = 3, only first stage of levelling is activated */
		gv = Levelling(			z_b,									/** bulk rock informations */
								gv,										/** global variables (e.g. Gamma) */

								SS_objective,
							    splx_data,
								PP_ref_db,								/** pure phase database */
								SS_ref_db,								/** solution phase database */
								cp						);
	}

	return gv;
}

/** 
  	Get command line options
*/
global_variable ReadCommandLineOptions(	global_variable 	 gv,		
										int 				 argc, 
										char 			   **argv, 	
										int 				*Mode_out, 
										int 				*Verb_out, 
										int 				*test_out, 
										int 				*n_points_out, 
										double 				*P, 
										double 				*T, 
										double 				 Bulk[11], 
										double 				 Gam[11], 
										char 				 File[50], 
										char 				 Phase[50], 
										int 				*maxeval_out,
										int 				*get_version_out,
										int					*get_help		
){
	int i;
	static ko_longopt_t longopts[] = {
		{ "Verb", 		ko_optional_argument, 301 },
		{ "Mode", 		ko_optional_argument, 302 },
		{ "File", 		ko_optional_argument, 303 },
		{ "n_points ",	ko_optional_argument, 304 },
		{ "test",  		ko_optional_argument, 305 },
		{ "Temp", 		ko_optional_argument, 306 },
		{ "Pres",  		ko_optional_argument, 307 },
		{ "Phase",  	ko_optional_argument, 308 },
		{ "n_pc", 		ko_optional_argument, 309 },
		{ "Gam",  		ko_optional_argument, 310 },
		{ "Bulk", 		ko_optional_argument, 311 },
        { "maxeval",    ko_optional_argument, 313 },
        { "version",    ko_optional_argument, 314 },
        { "help",    	ko_optional_argument, 315 },
    	{ NULL, 0, 0 }
	};
	ketopt_t opt = KETOPT_INIT;
	
	int    c;
	int    Mode     =  0;
	int    Verb     =  gv.verbose;
	int    test     =  0;
	int    n_points =  1;
	int    n_pc     =  2;		/** number of pseudocompounds for Mode 2 */
	int    maxeval  = -1;

	double Temp , Pres;
	Temp   = 1100.0;
	Pres   = 12.0;
	 
	for (i = 0; i < nEl; i++) {
		Bulk[i] = 0.0;
		Gam[i]  = 0.0;
	}


	strcpy(File,"none"); // Filename to be read to have multiple P-T-bulk conditions to solve

	while ((c = ketopt(&opt, argc, argv, 1, "", longopts)) >= 0) {
		if 		(c == 314){ printf("MAGEMin %20s\n",gv.version ); exit(0); }	
		else if (c == 315){ print_help( gv ); 					  exit(0); }	
        else if	(c == 301){ Verb     = atoi(opt.arg	);}
		else if (c == 302){ Mode     = atoi(opt.arg);			if (Verb == 1){		printf("--Mode        : Mode                     = %i \n", 	 	   		Mode		);}}																		
		else if (c == 303){ strcpy(File,opt.arg);		 		if (Verb == 1){		printf("--File        : File                     = %s \n", 	 	   		File		);}}
		else if (c == 304){ n_points = atoi(opt.arg); 	 		if (Verb == 1){		printf("--n_points    : n_points                 = %i \n", 	 	   		n_points	);}}
		else if (c == 305){ test     = atoi(opt.arg); 		 	if (Verb == 1){		printf("--test        : Test                     = %i \n", 	 	  		test		);}}
		else if (c == 306){ Temp     = strtof(opt.arg,NULL); 	if (Verb == 1){		printf("--Temp        : Temperature              = %f C \n",             Temp		);}}
		else if (c == 307){ Pres     = strtof(opt.arg,NULL); 	if (Verb == 1){		printf("--Pres        : Pressure                 = %f kbar \n", 			Pres		);}}
		else if (c == 308){ strcpy(Phase,opt.arg);		 		if (Verb == 1){		printf("--Phase       : Phase name               = %s \n", 	   			Phase		);}}
		else if (c == 313){ maxeval  = strtof(opt.arg,NULL); 	if (Verb == 1){
            if (maxeval==0){        printf("--maxeval     : Max. # of local iter.    = infinite  \n"		); }
            else{                   printf("--maxeval     : Max. # of local iter.    = %i  \n", maxeval		);}
            }
        }
		else if (c == 310){
			char *p = strtok(opt.arg,",");
			size_t i = 0;
			while(p && i<11) {
					Gam[i++] = atof(p);
					p = strtok(NULL, ",");
			}
			if (Verb == 1){
				printf("--Gam  	      : Gamma           = ");
				for (int j = 0; j < 11; j++){
					printf("%g ", Gam[j]);	
				} 
				printf(" dG \n");
			}
		 }
		else if (c == 311){
			char *p  = strtok(opt.arg,",");
			size_t i = 0;
			while(p && i<11) {
					Bulk[i++] = atof(p);
					p = strtok(NULL, ",");
			}
			if (Verb == 1){
				printf("--Bulk  	 : Bulk         = ");
				for (int j = 0; j < 11; j++){
					printf("%g ", Bulk[j]);	
				} 
				printf(" \n");
			}
		 }
	}

	/** Output */
	*Verb_out 		= 	Verb;
	*Mode_out 		= 	Mode;
	*test_out	 	= 	test;
	*P        		= 	Pres;
	*T        		= 	Temp;
	*n_points_out 	= 	n_points;
    *maxeval_out    =   maxeval;
	
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

	/* Allocate memory for each solution phase according to their specificities (n_em, sf etc) */
	for (i = 0; i < gv.len_ss; i++){
		DB.SS_ref_db[i] = G_SS_INIT_EM_function(		DB.SS_ref_db[i], 
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
	DB.sp[0] = SP_INIT_function(		DB.sp[0], gv							);

	/* Endmember names */
	DB.EM_names  =	get_EM_DB_names(EM_database);

	/* Create endmember Hashtable */
	struct EM2id *p_s, *tmp_p;
	struct EM_db EM_return;
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

	CP_destroy(			gv, 
						DB.cp				);
	
	free(DB.cp);
	
	for (int i = 0; i < n_em_db; i++) {
		free(DB.EM_names[i]);
	}
	
	free(DB.EM_names);
	free(DB.PP_ref_db);
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
	if (gv.Mode==0 && gv.verbose !=-1){
		printf(" Status             : %12i ",gv.status);
		if (gv.verbose == 1){PrintStatus(gv.status);}
		printf("\n");

    	printf(" Rank               : %12i \n",rank);
    	printf(" Point              : %12i \n",l);
    	printf(" Temperature        : %+12.5f\t [C] \n",   z_b.T - 273.15);
		printf(" Pressure           : %+12.5f\t [kbar]\n", z_b.P);
       
 		if (gv.verbose == 1){
 			printf("\n______________________________\n");
			printf("| Comp. Time: %.6f (ms) |", time_taken*1000);
            printf("\n══════════════════════════════\n");
        }
    }

	if ((gv.verbose != -1) &&  (gv.Mode==0)){
		printf("\nSOLUTION: [G = %.3f] (%i iterations, %.2f ms)\n",gv.G_system,gv.global_ite,time_taken*1000.0);
		printf("[");
		for (i = 0; i < z_b.nzEl_val; i++){
			printf("%+8f,",gv.gam_tot[z_b.nzEl_array[i]]);
		}
		printf("]\n");
		for (int i = 0; i < gv.len_cp; i++){
			if (DB.cp[i].ss_flags[1] == 1){
				printf("%4s \t %.5f \n", DB.cp[i].name, DB.cp[i].ss_n);
			}
		}
		for (int i = 0; i < gv.len_pp; i++){
			if (gv.pp_flags[i][1] == 1){
				printf("%4s \t %.5f \n", gv.PP_list[i], gv.pp_n[i]);
			}
		}	
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
