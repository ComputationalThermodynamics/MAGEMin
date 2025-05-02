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
#include "uthash.h"
#include "MAGEMin.h"
#include "initialize.h"
#include "toolkit.h"
#include <stdio.h>

/* Function to allocate the memory of the data to be used/saved during PGE iterations */
global_variable global_variable_init( 	global_variable  	 gv,
										bulk_info 			*z_b 	){


	if (strcmp(gv.research_group, "tc") 	== 0 ){
	/* here we initialize MAGEMin using the THERMOCALC formalism */
		gv 	=	global_variable_TC_init( 	gv,
											z_b 	);	
	}
	else if (strcmp(gv.research_group, "sb") 	== 0 ){
	/* here we initialize MAGEMin using Stixrude formalism */
		gv 	=	global_variable_SB_init( 	gv,
											z_b 	);
	}
	else{
		printf(" wrong group, fix group name\n");
	}

	return gv;
}


/**
    Function to retrieve the endmember names from the database 
    Note the size of the array is n_em_db+1, required for the hashtable              
*/
char** get_EM_DB_names_tc(global_variable gv) {
    EM_db EM_return;
    int i, n_em_db;
    n_em_db = gv.n_em_db;
    char ** names = malloc((n_em_db+1) * sizeof(char*));
    for ( i = 0; i < n_em_db; i++){
        names[i] = malloc(20 * sizeof(char));
    }
    for ( i = 0; i < n_em_db; i++){	
        EM_return = Access_EM_DB(i, gv.EM_dataset);
        strcpy(names[i],EM_return.Name);
    }
    return names;
}

char** get_EM_DB_names_sb(global_variable gv) {

    EM_db_sb EM_return;
    int i, n_em_db;
    n_em_db = gv.n_em_db;
    char ** names = malloc((n_em_db+1) * sizeof(char*));
    for ( i = 0; i < n_em_db; i++){
        names[i] = malloc(20 * sizeof(char));
    }
    for ( i = 0; i < n_em_db; i++){	
        EM_return = Access_SB_EM_DB(i, gv.EM_dataset);
        strcpy(names[i],EM_return.Name);
    }
    return names;
}

/**
    Function to retrieve the species names from the database 
    Note the size of the array is n_em_db+1, required for the hashtable              
*/
char** get_FS_DB_names(global_variable gv) {
    FS_db FS_return;
    int i, n_fs_db;
    n_fs_db = gv.n_fs_db;
    char ** names = malloc((n_fs_db+1) * sizeof(char*));
    for ( i = 0; i < n_fs_db; i++){
        names[i] = malloc(20 * sizeof(char));
    }
    for ( i = 0; i < n_fs_db; i++){	
        FS_return = Access_FS_DB(i);
        strcpy(names[i],FS_return.Name);
    }
    return names;
}


/** 
	set default parameters necessary to initialize the system 
*/
global_variable global_variable_alloc( bulk_info  *z_b ){
	global_variable gv;

	int i,j,k;

	/** 
		allocate data necessary to initialize the system 
	*/
	/* system parameters 		*/
	gv.maxlen_ox 		= 15;
	gv.outpath 			= malloc (100 	* sizeof(char)			);
	gv.version 			= malloc (50  	* sizeof(char)			);
	gv.File 			= malloc (50 	* sizeof(char)			);
	gv.research_group 	= malloc (5 	* sizeof(char)			);
	gv.db 				= malloc (5 	* sizeof(char)			);
	gv.sys_in 			= malloc (5 	* sizeof(char)			);
	gv.buffer 			= malloc (10 	* sizeof(char)			);

	gv.arg_bulk 		= malloc (gv.maxlen_ox * sizeof(double)	);
	gv.arg_gamma 		= malloc (gv.maxlen_ox * sizeof(double)	);

	for (i = 0; i < gv.maxlen_ox; i++) {
		gv.arg_bulk[i]  = 0.0;
		gv.arg_gamma[i] = 0.0;
	}

	strcpy(gv.outpath,"./output/");				/** define the outpath to save logs and final results file	 						*/
	strcpy(gv.version,"1.7.5 [02/04/2025]");	/** MAGEMin version 																*/

	/* generate parameters        		*/
	strcpy(gv.buffer,"none");
	gv.EM_dataset 		= -1;
	gv.max_n_mSS		= 128;					/** maximum number of metastable pseudocompounds 									*/
	gv.max_n_cp 		= 256;					/** number of considered solution phases 											*/	
	gv.max_ss_size_cp   = 24;					/** maximum size for a solution phase saved in the cp structure                     */
	gv.buffer_n 		= 0.0;					/** factor for QFM buffer 															*/
	gv.limitCaOpx       = 0;					/** limit Ca-bearing  orthopyroxene (add-hoc correction) 							*/
	gv.CaOpxLim         = 1.0;					/** limit Ca-bearing  orthopyroxene (add-hoc correction) 							*/

	/* Phase selection 					*/
	gv.mbCpx 			= 0;					/** 0: omphacite LT, 1: augite HT													*/
	gv.mbIlm 			= 0;					/** 0: Ilmm, 1: Ilm 																*/
	gv.mpSp 			= 0;					/** 0: Sp LT, 1: Mt1													*/
	gv.mpIlm 			= 0;					/** 0: Ilmm, 1: Ilm 																*/

	// gv.calc_seismic_cor = 1;					/** compute seismic velocity corrections (melt and anelastic)						*/
	// gv.melt_pressure 	= 0.0;				/** [kbar] pressure shift in case of modelling melt pressure 						*/

	/* fluid speciation parameters 	    */
	gv.fluidSpec        = 0;					/** by default the fluid speciation option is deactivated 							*/
	gv.n_fs_db 			= 44; 					/** number of fluid species for the database 										*/

	/* residual tolerance 				*/
	gv.br_max_tol       = 1.0e-5;				/** value under which the solution is accepted to satisfy the mass constraint 		*/

	/* pc composite parameters */
	gv.pc_composite_dist= 2.5e-3;				/** parameter setting the distance for the pseudocompounds created around a minimized point 
													this parameter has a big impact on performances, it is advised to not change it */
	
	/* Magic PGE under-relaxing numbers */
	gv.relax_PGE_val    = 128.0;				/** restricting factor 																*/
	gv.PC_check_val1	= 1.0e-2;				/** br norm under which PC are tested for potential candidate to be added 			*/
	gv.PC_check_val2	= 1.0e-4;				/** br norm under which PC are tested for potential candidate to be added 			*/
	gv.PC_min_dist 		= 1.0;					/** factor multiplying the diagonal of the hyperbox of xeos s-tep 					*/
	gv.mSS_df_max_add 	= 0.4;					/** driving force under which a metastable solution phase is added to the assemblage */
	gv.mSS_df_min_add   = 1e-6;					/** driving force under which a metastable solution phase is added to the assemblage */

	/* levelling parameters 			*/
	gv.em2ss_shift		= 2e-7;					/** small value to shift x-eos of pure endmember from bounds after levelling 		*/
	gv.bnd_filter_pc    = 10.0;					/** value of driving force the pseudocompound is considered 						*/
	gv.bnd_filter_pge   = 2.5;					/** value of driving force the pseudocompound is considered 						*/
	gv.max_G_pc         = 2.5;					/** dG under which PC is considered after their generation		 					*/
	gv.eps_sf_pc		= 1e-12;				/** Minimum value of site fraction under which PC is rejected, 
													don't put it too high as it will conflict with bounds of x-eos					*/

	/* PGE LP pseudocompounds parameters */
	gv.launch_PGE 		= 0;
	gv.n_pc 			= 8192;
	gv.n_Ppc			= 8192;
	gv.max_LP_ite 		= 256;
	gv.save_Ppc_val     = 0.0; 					/** During PGE iterations, if the driving force is < save_Ppc_val, then the 
													pseudocompound is added to the Ppc list 										*/

	/* local minimizer options 	*/
	gv.bnd_val          = 1.0e-7;				/** boundary value for x-eos 										 				*/
	gv.box_size_mode_PGE= 0.25;					/** box edge size of the compositional variables used during PGE local minimization */
	gv.maxeval   		= 1024;					/** max number of evaluation of the obj function for mode 1 (PGE)					*/
	gv.maxgmTime        = 0.1; 					/** set a maximum minimization time for the local minimizer (sec)					*/
	gv.box_size_mode_LP	= 1.0;					/** box edge size of the compositional variables used during PGE local minimization */

	/* set of parameters to record the evolution of the norm of the mass constraint */
	gv.it_1             = 128;                  /** first critical iteration                                                        */
	gv.ur_1             = 4.;                   /** under relaxing factor on mass constraint if iteration is bigger than it_1       */
	gv.it_2             = 160;                  /** second critical iteration                                                       */
	gv.ur_2             = 8.;                   /** under relaxing factor on mass constraint if iteration is bigger than it_2       */
	gv.it_3             = 192;                  /** third critical iteration                                                        */
	gv.ur_3             = 16.;                  /** under relaxing factor on mass constraint if iteration is bigger than it_3       */
	gv.it_f             = 256;                  /** gives back failure when the number of iteration is bigger than it_f             */

	/* phase update options 			*/
	gv.min_df 			= -1e-6;					/** value under which a phase in hold is reintroduced */
	gv.re_in_df 		= -1e-6;
	/* numerical derivatives P,T steps (same value as TC) */
	gv.gb_P_eps			= 2e-3;					/** small value to calculate V using finite difference: V = dG/dP;					*/
	gv.gb_T_eps			= 2e-3;					/** small value to calculate V using finite difference: V = dG/dP;					*/
	gv.poisson_ratio 	= 0.3;					/** poisson ratio to compute elastic shear modulus 									*/

	/* initialize other values 			*/
	gv.mean_sum_xi		= 1.0;
	gv.sigma_sum_xi		= 1.0;
	gv.alpha        	= gv.max_fac;			/** active under-relaxing factor 													*/
	gv.tot_min_time 	= 0.0;
	gv.tot_time 		= 0.0;

	/* set default parameters (overwritten later from args)*/
	gv.EM_database  	=  0; 					
	gv.n_points 		=  1;
	gv.solver   		=  2;					/* 0 -> Legacy, 1 = PGE, Hybrid PGE/LP */
	gv.leveling_mode	=  0;
	gv.verbose 			=  0;
	gv.output_matlab 	=  0;
	gv.test     		= -1;

	/* default PT conditions for test */
	z_b->P 				= 12.0;		
	z_b->T 				= 1100.0 + 273.15;		
	z_b->R 				= 0.0083144;

	strcpy(gv.File,				"none"); 	/** Filename to be read to have multiple P-T-bulk conditions to solve 	*/
	strcpy(gv.sys_in,			"mol"); 	/** system unit 														*/
	strcpy(gv.db,				"ig"); 		/** database															*/
	strcpy(gv.research_group,	"tc"); 		/** Research group, THERMOCALC(tc) or  Stixrude-Lithgow-Bertelloni(sb)	*/

	return gv;
}


/** 
  allocate memory to store considered phases

  note (01/05/2023):
  ------------------

  n (array dim) was originally equal to gv.len_ox + 1; However this has been changed following update
  on liq model (Green et al., 2023) and the use of 4 order variables summing up the total
  number of compositional variables to 14 for a chemical system of 11 oxides
  This means that the max number of the endmembers can be >= n-oxide + 1
*/
csd_phase_set CP_INIT_function(csd_phase_set cp, global_variable gv){
	
	int n 			= gv.max_ss_size_cp;
	
	/* initialize fractions flags and cycle arrays with zeros */
	cp.ss_flags 	= malloc (gv.n_flags  * sizeof(int));
	cp.name   		= malloc (20 * sizeof(char)		    );
	cp.p_em   		= malloc (n  * sizeof(double) 		);
	cp.xi_em  		= malloc (n  * sizeof(double) 		);
	cp.dguess 		= malloc (n  * sizeof(double) 		);
	cp.xeos   		= malloc (n  * sizeof(double) 		);
	cp.xeos_0   	= malloc (n  * sizeof(double) 		);
	cp.xeos_1   	= malloc (n  * sizeof(double) 		);
	cp.xeos_r   	= malloc (n  * sizeof(double) 		);
	cp.delta_mu		= malloc (n  * sizeof(double) 		);
	cp.dfx   		= malloc (n  * sizeof(double) 		);
	cp.mu    		= malloc (n  * sizeof(double) 		);
	cp.gbase    	= malloc (n  * sizeof(double) 		);
	cp.ss_comp		= malloc (n  * sizeof(double) 		);
	cp.ss_comp_mol  = malloc (n  * sizeof(double) 		);
	cp.ss_comp_wt	= malloc (n  * sizeof(double) 		);
	cp.sf			= malloc ((n*2)  * sizeof(double) 	);
	
	cp.phase_density  		= 0.0;
	cp.phase_cp				= 0.0;
	cp.phase_expansivity	= 0.0;
	cp.phase_entropy		= 0.0;
	cp.phase_enthalpy		= 0.0;
		
	return cp;
}

/** 
  allocate memory to store considered phases
*/
stb_system SP_INIT_function(stb_system sp, global_variable gv){

	sp.MAGEMin_ver   		= malloc(50  		* sizeof(char)				);
	sp.dataset   		    = malloc(50  		* sizeof(char)				);
	sp.database   			= malloc(50  		* sizeof(char)				);
	sp.oxides 	     		= malloc(gv.len_ox  * sizeof(char*)				);
	sp.elements 	     	= malloc(gv.len_ox  * sizeof(char*)				);
		
	for (int i = 0; i < gv.len_ox; i++){
		sp.oxides[i] 		= malloc(20 * sizeof(char));	
		sp.elements[i] 		= malloc(20 * sizeof(char));	
	}
	sp.buffer 				= malloc(50  		* sizeof(char)				);
	sp.bulk 				= malloc(gv.len_ox  * sizeof(double)			);	
	sp.gamma 				= malloc(gv.len_ox  * sizeof(double)			);	
	sp.bulk_S 				= malloc(gv.len_ox 	* sizeof(double)			);	
	sp.bulk_M 				= malloc(gv.len_ox 	* sizeof(double)			);	
	sp.bulk_F 				= malloc(gv.len_ox 	* sizeof(double)			);	

	sp.bulk_wt 				= malloc(gv.len_ox  * sizeof(double)			);	
	sp.bulk_S_wt 			= malloc(gv.len_ox 	* sizeof(double)			);	
	sp.bulk_M_wt 			= malloc(gv.len_ox 	* sizeof(double)			);	
	sp.bulk_F_wt 			= malloc(gv.len_ox 	* sizeof(double)			);	
	sp.ph 	     			= malloc(gv.len_ox  * sizeof(char*)				);
	sp.sol_name 	     	= malloc(gv.len_ox  * sizeof(char*)				);
	sp.ph_frac 	     		= malloc(gv.len_ox  * sizeof(double)			);
	sp.ph_frac_wt     		= malloc(gv.len_ox  * sizeof(double)			);
	sp.ph_frac_1at    		= malloc(gv.len_ox  * sizeof(double)			);
	sp.ph_frac_vol     		= malloc(gv.len_ox  * sizeof(double)			);
	for (int i = 0; i < gv.len_ox; i++){
		sp.ph[i] 			= malloc(20 * sizeof(char));	
		sp.sol_name[i] 		= malloc(20 * sizeof(char));	
	}
	sp.ph_type 				= malloc(gv.len_ox 	* sizeof(int)				);	
	sp.ph_id 				= malloc(gv.len_ox 	* sizeof(int)				);	
	sp.ph_id_db 			= malloc(gv.len_ox 	* sizeof(int)				);	
	sp.PP 		 			= malloc(gv.len_ox  * sizeof(stb_PP_phase)		); 
	sp.SS 		 			= malloc(gv.len_ox  * sizeof(stb_SS_phase)		); 
	sp.mSS 		 			= malloc(gv.max_n_mSS  * sizeof(mstb_SS_phase)	); 

	for (int n = 0; n< gv.len_ox; n++){
		sp.PP[n].Comp 			= malloc(gv.len_ox 	* sizeof(double)		);
		sp.SS[n].Comp 			= malloc(gv.len_ox 	* sizeof(double)		);
		sp.PP[n].Comp_wt 		= malloc(gv.len_ox 	* sizeof(double)		);
		sp.SS[n].Comp_wt 		= malloc(gv.len_ox 	* sizeof(double)		);
		sp.PP[n].Comp_apfu		= malloc(gv.len_ox 	* sizeof(double)		);
		sp.SS[n].Comp_apfu		= malloc(gv.len_ox 	* sizeof(double)		);        
		sp.SS[n].compVariables	= malloc(gv.len_ox*3   * sizeof(double)	    );
        sp.SS[n].siteFractions	= malloc(gv.len_ox*3   * sizeof(double)	    );
		sp.SS[n].emFrac			= malloc((gv.len_ox*3) * sizeof(double)		);
		sp.SS[n].emFrac_wt		= malloc((gv.len_ox*3) * sizeof(double)		);
		sp.SS[n].emChemPot		= malloc((gv.len_ox*3) * sizeof(double)		);
		sp.SS[n].compVariablesNames	= malloc(gv.len_ox*3 * sizeof(char*)	);
		sp.SS[n].siteFractionsNames	= malloc(gv.len_ox*3 * sizeof(char*)	);
		sp.SS[n].emNames 	    = malloc((gv.len_ox*3) * sizeof(char*)		);
		sp.SS[n].emComp 	    = malloc((gv.len_ox*3) * sizeof(double*)	);
		sp.SS[n].emComp_wt 	    = malloc((gv.len_ox*3) * sizeof(double*)	);
        sp.SS[n].emComp_apfu    = malloc((gv.len_ox*3) * sizeof(double*)	);

		for (int i = 0; i < gv.len_ox*3; i++){
            sp.SS[n].compVariablesNames[i]		= malloc(20 * sizeof(char)	);
            sp.SS[n].siteFractionsNames[i]		= malloc(20 * sizeof(char)	);
			sp.SS[n].emNames[i]		= malloc(20 * sizeof(char)				);
			sp.SS[n].emComp[i]		= malloc(gv.len_ox * sizeof(double)		);		
			sp.SS[n].emComp_wt[i]	= malloc(gv.len_ox * sizeof(double)		);		
            sp.SS[n].emComp_apfu[i]	= malloc(gv.len_ox * sizeof(double)		);		
		}
	}
    
    /** allocate memory for metastable phases len_ox * 2 to be safe?        */
	for (int n = 0; n< gv.max_n_mSS; n++){
        sp.mSS[n].ph_name	    = malloc(20 * sizeof(char)	                );
        sp.mSS[n].ph_type	    = malloc(20 * sizeof(char)	                );
        sp.mSS[n].info	        = malloc(20 * sizeof(char)	                );
        sp.mSS[n].comp_Ppc 	    = malloc((gv.len_ox) 	* sizeof(double)	);  
        sp.mSS[n].p_Ppc 	    = malloc((gv.len_ox*2) 	* sizeof(double)	);  
        sp.mSS[n].mu_Ppc 	    = malloc((gv.len_ox*2) 	* sizeof(double)	);  
        sp.mSS[n].xeos_Ppc 	    = malloc((gv.len_ox*2) 	* sizeof(double)	);  
	}

	return sp;
}



/**
  attributes the right solution phase to the solution phase array
*/
SS_ref G_SS_init_EM_function(		SS_init_type		*SS_init,
									int			 		 ph_id,
									SS_ref 		 		 SS_ref_db,
									char 				*name, 
									global_variable 	 gv					){

	/* Retrieve the right data in the right place 	*/
	SS_ref_db  = (*SS_init[ph_id])(			SS_ref_db,
											gv							);
	/**
		Allocate memory for solution phase models and pseudocompound storage (memory is initialized in the reset function)
	*/
	int n_em   = SS_ref_db.n_em;
	int n_xeos = SS_ref_db.n_xeos;
	int n_sf   = SS_ref_db.n_sf;
	int sym    = SS_ref_db.symmetry;
	int n_cat  = SS_ref_db.n_cat;

    SS_ref_db.orderVar       = 0;
	SS_ref_db.fName 		 = malloc(20 * sizeof(char)		);		
	SS_ref_db.EM_list 		 = malloc ((n_em) * sizeof (char*)	);
	for (int i = 0; i < n_em; i++){ 
		SS_ref_db.EM_list[i] = malloc(20 * sizeof(char)		);		
	}
	SS_ref_db.CV_list 		 = malloc ((n_xeos) * sizeof (char*)	);
	for (int i = 0; i < n_xeos; i++){ 
		SS_ref_db.CV_list[i] = malloc(20 * sizeof(char)		);		
	}
	SS_ref_db.SF_list 		 = malloc ((n_sf) * sizeof (char*)	);
	for (int i = 0; i < n_sf; i++){ 
		SS_ref_db.SF_list[i] = malloc(20 * sizeof(char)		);		
	}
	if (sym == 0){
		SS_ref_db.W   		= malloc (SS_ref_db.n_w * sizeof (double) ); 
		SS_ref_db.v   		= malloc (SS_ref_db.n_v * sizeof (double) ); 
	}
	else if (sym == 1){
		SS_ref_db.W   		= malloc (SS_ref_db.n_w * sizeof (double) ); 
	}
	
	/* initialize fractions flags and cycle arrays with zeros */
	SS_ref_db.ss_flags      = malloc (gv.n_flags  * sizeof(int));
    SS_ref_db.solvus_id     = malloc ((gv.len_ss*4)   * sizeof (int)  	);

	/* dynamic memory allocation of data to send to NLopt */
	SS_ref_db.bounds 	= malloc (n_xeos * sizeof (double*)  ); 
	for (int i = 0; i < n_xeos; i++){
		SS_ref_db.bounds[i] = malloc (2 * sizeof (double) );
	}
	
	SS_ref_db.eye 			= malloc (n_em 			* sizeof (double*)); 
	SS_ref_db.dp_dx 		= malloc (n_em 			* sizeof (double*)); 
	SS_ref_db.Comp 			= malloc (n_em 			* sizeof (double*)); 
	for (int i = 0; i < n_em; i++){
		SS_ref_db.eye[i] 	= malloc (n_em	 		* sizeof (double) );
		SS_ref_db.Comp[i] 	= malloc (gv.len_ox 	* sizeof (double) );
		SS_ref_db.dp_dx[i] 	= malloc (n_xeos 		* sizeof (double) );
	}
	
	SS_ref_db.gbase   		= malloc (n_em   	 	* sizeof (double) ); 
	SS_ref_db.gb_lvl  		= malloc (n_em   	 	* sizeof (double) ); 
	SS_ref_db.z_em    		= malloc (n_em   	 	* sizeof (double) ); 
	SS_ref_db.d_em    		= malloc (n_em   	 	* sizeof (double) );
	SS_ref_db.density 		= malloc (n_em   	 	* sizeof (double) ); 
	SS_ref_db.dguess 		= malloc (n_xeos 		* sizeof (double) );
	SS_ref_db.iguess  		= malloc (n_xeos   	  	* sizeof (double) );
	SS_ref_db.mguess  		= malloc (n_xeos   	  	* sizeof (double) );
	SS_ref_db.idOrderVar  	= malloc (n_xeos   	  	* sizeof (double) );

	SS_ref_db.p       		= malloc (n_em       	* sizeof (double) ); 
	SS_ref_db.ElShearMod    = malloc (n_em       	* sizeof (double) ); 
	SS_ref_db.ElBulkMod     = malloc (n_em       	* sizeof (double) ); 
	SS_ref_db.ElCp    		= malloc (n_em       	* sizeof (double) ); 
	SS_ref_db.ElExpansivity = malloc (n_em       	* sizeof (double) ); 
	SS_ref_db.ape      		= malloc (n_em       	* sizeof (double) ); 
	SS_ref_db.mat_phi 		= malloc (n_em       	* sizeof (double) ); 
	SS_ref_db.mu_Gex  		= malloc (n_em       	* sizeof (double) ); 
	SS_ref_db.sf      		= malloc (n_sf       	* sizeof (double) ); 
	SS_ref_db.mu      		= malloc (n_em       	* sizeof (double) ); 
	SS_ref_db.dfx    		= malloc (n_xeos     	* sizeof (double) ); 
	SS_ref_db.ss_comp		= malloc (gv.len_ox  	* sizeof (double) ); 
	SS_ref_db.ElEntropy		= malloc (gv.len_ox  	* sizeof (double) );     
	SS_ref_db.xi_em   		= malloc (n_em   	 	* sizeof (double) ); 
	SS_ref_db.xeos    		= malloc (n_xeos     	* sizeof (double) ); 	

	/* memory allocation to store all gbase */
	SS_ref_db.mu_array = malloc ((gv.n_Diff) * sizeof (double*) ); 
	for (int i = 0; i < (gv.n_Diff); i++){
		SS_ref_db.mu_array[i] = malloc (n_em * sizeof (double) );
	}	
	
	/* dynamic memory allocation of data to send to NLopt */
	SS_ref_db.bounds_ref = malloc ((n_xeos) * sizeof (double*) ); 
	for (int i = 0; i < (n_xeos); i++){
		SS_ref_db.bounds_ref[i] = malloc (2 * sizeof (double) );
	}
	
	/* allocate memory when using Stixrude database */
	if (n_cat > 0){
		/* dynamic memory allocation of data to send to NLopt */
		SS_ref_db.C = malloc ((n_cat) * sizeof (double*) ); 
		for (int i = 0; i < (n_cat); i++){
			SS_ref_db.C[i] = malloc (n_em * sizeof (double) );
		}
		SS_ref_db.N = malloc ((n_em) * sizeof (double*) ); 
		for (int i = 0; i < (n_em); i++){
			SS_ref_db.N[i] = malloc ((n_em-1) * sizeof (double) );
		}
		SS_ref_db.Vec1 = malloc ((n_em-1) * sizeof (double) );
		SS_ref_db.Vec2 = malloc ((n_em) * sizeof (double) );
	}

	/* dynamic memory allocation of data to send to NLopt */
	SS_ref_db.ub   		= malloc ((n_xeos) * sizeof (double) ); 
	SS_ref_db.lb   		= malloc ((n_xeos) * sizeof (double) ); 

	/**
		Allocate memory for levelling pseudocompounds 
	*/
	// SS_ref_db.n_pc   	= gv.n_pc;
	SS_ref_db.n_pc   	= gv.n_SS_PC[ph_id];
    SS_ref_db.tot_pc   	= malloc (1 * sizeof (double) ); 
    SS_ref_db.id_pc   	= malloc (1 * sizeof (double) ); 
	SS_ref_db.G_pc   	= malloc ((SS_ref_db.n_pc) * sizeof (double) ); 
	SS_ref_db.DF_pc 	= malloc ((SS_ref_db.n_pc) * sizeof (double) ); 
	SS_ref_db.factor_pc = malloc ((SS_ref_db.n_pc) * sizeof (double) ); 
	SS_ref_db.info  	= malloc ((SS_ref_db.n_pc) * sizeof (int) 	 ); 
	SS_ref_db.p_pc 		= malloc ((SS_ref_db.n_pc) * sizeof (double*)); 
	// SS_ref_db.mu_pc 	= malloc ((SS_ref_db.n_pc) * sizeof (double*)); 
	
	for (int i = 0; i < (SS_ref_db.n_pc); i++){
		SS_ref_db.p_pc[i] 	 = malloc ((n_em) * sizeof (double) 	);
		// SS_ref_db.mu_pc[i] 	 = malloc ((n_em) * sizeof (double) 	);
	}
	SS_ref_db.comp_pc = malloc ((SS_ref_db.n_pc) * sizeof (double*) ); 
	for (int i = 0; i < (SS_ref_db.n_pc); i++){
		SS_ref_db.comp_pc[i] = malloc (gv.len_ox * sizeof (double) 	);
	}
	SS_ref_db.xeos_pc = malloc ((SS_ref_db.n_pc) * sizeof (double*) ); 
	for (int i = 0; i < (SS_ref_db.n_pc); i++){
		SS_ref_db.xeos_pc[i] = malloc ((n_xeos)  * sizeof (double) 	);
	}	

	/**
		Allocate memory for PGE pseudocompounds 
	*/
	SS_ref_db.n_Ppc   	= gv.n_Ppc;								/** maximum number of pseudocompounds to store */
	SS_ref_db.G_Ppc   	= malloc ((SS_ref_db.n_Ppc) * sizeof (double) ); 
	SS_ref_db.DF_Ppc 	= malloc ((SS_ref_db.n_Ppc) * sizeof (double) ); 
	SS_ref_db.info_Ppc 	= malloc ((SS_ref_db.n_Ppc) * sizeof (int) 	 ); 
	SS_ref_db.p_Ppc 	= malloc ((SS_ref_db.n_Ppc) * sizeof (double*)); 
	SS_ref_db.mu_Ppc 	= malloc ((SS_ref_db.n_Ppc) * sizeof (double*)); 
	
	for (int i = 0; i < (SS_ref_db.n_Ppc); i++){
		SS_ref_db.p_Ppc[i] 	 = malloc ((n_em) * sizeof (double) 		);
		SS_ref_db.mu_Ppc[i]  = malloc ((n_em) * sizeof (double) 		);
	}
	SS_ref_db.comp_Ppc = malloc ((SS_ref_db.n_Ppc) * sizeof (double*) 	); 
	for (int i = 0; i < (SS_ref_db.n_Ppc); i++){
		SS_ref_db.comp_Ppc[i] = malloc (gv.len_ox * sizeof (double) 	);
	}
	SS_ref_db.xeos_Ppc = malloc ((SS_ref_db.n_Ppc) * sizeof (double*) 	); 
	for (int i = 0; i < (SS_ref_db.n_Ppc); i++){
		SS_ref_db.xeos_Ppc[i] = malloc ((n_xeos)  * sizeof (double) 	);
	}	

	/* initiliazes eye matrix as there is no need to redo it afterward */
	for (int i = 0; i < n_em; i++){
		for (int j = 0; j < n_em; j++){
			SS_ref_db.eye[i][j] = 0.0;
		}	
	}
	
	/* initialize eye matrix */
	for (int j = 0; j < n_em; j++){
		SS_ref_db.eye[j][j]= 1.0;
	}

	for (int j = 0; j < n_xeos; j++){
		SS_ref_db.idOrderVar[j]= 1.0;
	}
    


	return SS_ref_db;
};


/**
	structure to transfer composition, oversized on purpose to accomodate for database with higher oxide number 
*/
typedef struct get_datas{
	double comp[15];			
} get_data;

void init_data(int len_ox,	void *comp_array ){
	get_data *d  = (get_data *) comp_array;
	for (int i = 0; i < len_ox; i++){
		d->comp[i] = 0.0;
	}
}

void init_pp(int len_ox, void *PP_db ){
	PP_ref *d  = (PP_ref *) PP_db;
	for (int i = 0; i < len_ox; i++){
		d->Comp[i] = 0.0;
	}
}


/** 
  function to easely get gb and comp in order to define solid solutions
*/
em_data get_em_data(	char        *research_group,
                        int 		 EM_dataset, 
						int          len_ox,
						bulk_info 	 z_b,
                        double       P,
                        double       T,
						char 		*name, 
						char 		*state		){

	em_data data; 
	PP_ref PP_db   		= G_EM_function(research_group, EM_dataset, len_ox, z_b.id, z_b.bulk_rock, z_b.apo, P, T, name, state);
   	data.ElShearMod  	= PP_db.phase_shearModulus;
	data.ElBulkMod  	= PP_db.phase_bulkModulus;
	data.ElCp  		 	= PP_db.phase_cp;
	data.ElExpansivity 	= PP_db.phase_expansivity;
   	data.gb  			= PP_db.gbase;

	for (int i = 0; i < len_ox; i++){
		data.C[i] = PP_db.Comp[i];
	}
	return data;
}



/**
  reset global variable for parallel calculations 
*/
global_variable reset_gv(					global_variable 	 gv,
											bulk_info 	 		 z_b,
											PP_ref 				*PP_ref_db,
											SS_ref 				*SS_ref_db
){

	int i,j,k;
	for (k = 0; k < gv.n_flags; k++){
		for (i = 0; i < gv.len_pp; i++){
			gv.pp_flags[i][k]   = 0;
		}
		for (int i = 0; i < gv.len_ss; i++){ 
			SS_ref_db[i].ss_flags[k]   = 0;
		}
	}
	
	// gv.H2O_id 	= -1;
	// gv.CaO_id 	= -1;
	// gv.Na2O_id 	= -1;
	// gv.FeO_id 	= -1;
	// gv.MgO_id 	= -1;
	// gv.K2O_id 	= -1;
	// gv.O_id 	= -1;
	// gv.MnO_id 	= -1;

	/* reset pure phases fractions and xi */
	for (int i = 0; i < gv.len_pp; i++){		
		gv.pp_n[i] 		  = 0.0;
		gv.pp_n_mol[i]	  = 0.0;
		gv.pp_n_wt[i]	  = 0.0;
		gv.delta_pp_n[i]  = 0.0;
		gv.pp_xi[i] 	  = 0.0;
		gv.delta_pp_xi[i] = 0.0;
	}
	
	/* reset pure phases */
	char liq_tail[] = "L";
	for (int i = 0; i < gv.len_pp; i++){
		if ( EndsWithTail(gv.PP_list[i], liq_tail) == 1 ){
			if (z_b.T < gv.min_melt_T){
				gv.pp_flags[i][0] = 0;
				gv.pp_flags[i][1] = 0;
				gv.pp_flags[i][2] = 0;
				gv.pp_flags[i][3] = 1;
				gv.pp_flags[i][4] = 0;
			}
			else{
				gv.pp_flags[i][0] = 1;
				gv.pp_flags[i][1] = 0;
				gv.pp_flags[i][2] = 0;
				gv.pp_flags[i][3] = 0;
				gv.pp_flags[i][4] = 0;
			}
		}
		else{

			if(strcmp( gv.PP_list[i], "O2") == 0){
				gv.pp_flags[i][0] = 0;
				gv.pp_flags[i][1] = 0;
				gv.pp_flags[i][2] = 0;
				gv.pp_flags[i][3] = 1;
				gv.pp_flags[i][4] = 0;
			}
			else{
				gv.pp_flags[i][0] = 1;
				gv.pp_flags[i][1] = 0;
				gv.pp_flags[i][2] = 1;
				gv.pp_flags[i][3] = 0;
				gv.pp_flags[i][4] = 0;
			}

		}
	}
	gv.leveling_mode	  = 0;
	gv.tot_time 	  	  = 0.0;
	gv.tot_min_time 	  = 0.0;
	gv.melt_fraction	  = 0.;
	gv.solid_fraction	  = 0.;
	gv.melt_density       = 0.;
	gv.melt_bulkModulus   = 0.;
	gv.launch_PGE		  = 0;

	gv.solid_density      = 0.;
	gv.solid_bulkModulus  = 0.;
	gv.solid_shearModulus = 0.;
	gv.solid_Vp 		  = 0.;
	gv.solid_Vs 		  = 0.;

	// gv.melt_pressure 	  = 0.;
	gv.system_fO2 		  = 0.;
	gv.system_deltaQFM	  = 0.;
	gv.system_aH2O	  	  = 0.;
	gv.system_aSiO2	  	  = 0.;
	gv.system_aTiO2	  	  = 0.;
	gv.system_aAl2O3  	  = 0.;
	gv.system_aMgO  	  = 0.;
	gv.system_aFeO  	  = 0.;
		
	gv.system_density     = 0.;
	gv.system_entropy     = 0.;
	gv.system_enthalpy    = 0.;
	gv.system_volume   	  = 0.;
	gv.system_cp    	  = 0.;
	gv.system_expansivity = 0.;
	gv.system_bulkModulus = 0.;
	gv.system_shearModulus= 0.;
	gv.system_Vp 		  = 0.;
	gv.system_Vs 		  = 0.;
	gv.V_cor[0]			  = 0.;
	gv.V_cor[1]			  = 0.;
	gv.PC_checked		  = 0;
	gv.check_PC1		  = 0;
	gv.check_PC2		  = 0;
	gv.len_cp 		  	  = 0;
	gv.ph_change  	      = 0;
	gv.BR_norm            = 1.0;					/** start with 1.0 																	*/
	gv.G_system           = 0.0;	
	gv.div				  = 0;
	gv.status			  = 0;
	gv.n_phase            = 0;                  /** reset the number of phases to start with */
	gv.global_ite		  = 0;					/** reset the number of global iteration to zero */
	gv.n_cp_phase         = 0;					/** reset the number of ss phases to start with */
	gv.n_pp_phase         = 0;					/** reset the number of pp phases to start with */
	gv.alpha          	  = gv.max_fac;

	/* reset iteration record */
	for (i = 0; i < gv.it_f; i++){	
		gv.Alg[i] 				= 0;
		gv.PGE_mass_norm[i] 	= 0.0;
		gv.gamma_norm[i] 		= 0.0;	
		gv.gibbs_ev[i] 			= 0.0;			
		gv.ite_time[i] 			= 0.0;
	}

	/* reset norm and residuals */
    for (i = 0; i < gv.len_ox; i++){	
        gv.mass_residual[i] = 0.0;
		gv.dGamma[i]     	= 0.0;
        gv.gam_tot[i]     	= 0.0;
        gv.gam_tot_0[i]   	= 0.0;
        gv.delta_gam_tot[i] = 0.0;
		gv.mass_residual[i] = 0.0;	
    }

    for (i = 0; i < gv.len_ss; i++){	
        gv.n_solvi[i] = 0;
		gv.n_ss_array[i] = 0.0;
    }

	for (i = 0; i < (gv.len_ox); i++){ 
		gv.pc_id[i] = -1;
		gv.b[i] 	= 0.0;
		gv.b1[i] 	= 0.0;
		gv.tmp1[i] 	= 0.0;
		gv.tmp2[i] 	= 0.0;
		gv.tmp3[i] 	= 0.0;
		for (j = 0; j < gv.len_ox; j++){
			gv.A[i][j] = 0.0;
			gv.A2[i][j] = 0.0;
		}
	}

	return gv;
}

/**
  reset stable phases entries
*/
void reset_sp(						global_variable 	 gv,
									stb_system  		*sp
){
	sp[0].X  							= 1.0;
	sp[0].M_sys							= 0.0;
	sp[0].aH2O	  	  					= 0.0;
	sp[0].aSiO2	  						= 0.0;
	sp[0].aTiO2  						= 0.0;
	sp[0].aAl2O3  						= 0.0;
	sp[0].aMgO  						= 0.0;
	sp[0].aFeO  						= 0.0;

	sp[0].alpha  						= 0.0;
	sp[0].cp  							= 0.0;
	sp[0].s_cp  						= 0.0;
	sp[0].cp_wt  						= 0.0;
	sp[0].V  							= 0.0;

	sp[0].frac_S_wt						= 0.0;
	sp[0].frac_M_wt						= 0.0;
	sp[0].frac_F_wt						= 0.0;

	sp[0].frac_S_vol					= 0.0;
	sp[0].frac_M_vol					= 0.0;
	sp[0].frac_F_vol					= 0.0;

	sp[0].frac_S						= 0.0;
	sp[0].frac_M						= 0.0;
	sp[0].frac_F						= 0.0;

	/* reset system */
	for (int i = 0; i < gv.len_ox; i++){
		strcpy(sp[0].ph[i],"");	
		strcpy(sp[0].sol_name[i],"");	
		sp[0].bulk[i] 					= 0.0;
		sp[0].bulk_wt[i] 				= 0.0;
		sp[0].gamma[i] 					= 0.0;
		sp[0].bulk_S[i] 				= 0.0;
		sp[0].bulk_M[i] 				= 0.0;
		sp[0].bulk_F[i] 				= 0.0;
		sp[0].bulk_S_wt[i] 				= 0.0;
		sp[0].bulk_M_wt[i] 				= 0.0;
		sp[0].bulk_F_wt[i] 				= 0.0;
		sp[0].ph_type[i] 				= -1;
		sp[0].ph_id[i] 					=  0;
		sp[0].ph_frac[i] 				=  0.0;
		sp[0].ph_frac_wt[i] 			=  0.0;
		sp[0].ph_frac_1at[i] 			=  0.0;
		sp[0].ph_frac_vol[i] 			=  0.0;
	}

	/* reset phases */
	for (int n = 0; n < gv.len_ox; n++){
		for (int i = 0; i < gv.len_ox; i++){
			sp[0].PP[n].Comp[i] 			= 0.0;
			sp[0].SS[n].Comp[i] 			= 0.0;
			sp[0].PP[n].Comp_wt[i] 			= 0.0;
			sp[0].SS[n].Comp_wt[i] 			= 0.0;
			sp[0].PP[n].Comp_apfu[i] 			= 0.0;
			sp[0].SS[n].Comp_apfu[i] 			= 0.0;
		}
		for (int i = 0; i < gv.len_ox*3; i++){
			sp[0].SS[n].compVariables[i] 	= 0.0;
			sp[0].SS[n].siteFractions[i] 	= 0.0;

			strcpy(sp[0].SS[n].compVariablesNames[i],"");	
			strcpy(sp[0].SS[n].siteFractionsNames[i],"");	
			strcpy(sp[0].SS[n].emNames[i],"");	

			sp[0].SS[n].emFrac[i] 			= 0.0;
			sp[0].SS[n].emFrac_wt[i] 		= 0.0;
			sp[0].SS[n].emChemPot[i] 		= 0.0;

			for (int j = 0; j < gv.len_ox; j++){
				sp[0].SS[n].emComp[i][j]	= 0.0;
				sp[0].SS[n].emComp_wt[i][j]	= 0.0;
				sp[0].SS[n].emComp_apfu[i][j]	= 0.0;
			}
		}
	}

	/* reset metastable phases */
	for (int n = 0; n < gv.max_n_mSS; n++){
		strcpy(sp[0].mSS[n].ph_name,"");
		strcpy(sp[0].mSS[n].ph_type,"");
		strcpy(sp[0].mSS[n].info,"");

		for (int i = 0; i < gv.len_ox; i++){
			sp[0].mSS[n].comp_Ppc[i] 		= 0.0;
		}
		for (int i = 0; i < gv.len_ox*2; i++){
			sp[0].mSS[n].p_Ppc[i] 			= 0.0;
			sp[0].mSS[n].mu_Ppc[i] 			= 0.0;
			sp[0].mSS[n].xeos_Ppc[i] 		= 0.0;
		}
	}

}

/**
  reset bulk rock composition informations (needed if the bulk-rock is not constant during parallel computation)
*/
bulk_info reset_z_b_bulk(			global_variable 	 gv,
									bulk_info 	 		 z_b
){
	int i, j, k;

	int sum = 0;
	for (i = 0; i < gv.len_ox; i++) {
		z_b.zEl_array[i] = 0.0;
		z_b.bulk_rock[i] = gv.bulk_rock[i];
		// if (gv.bulk_rock[i] > 0.0){
		if (gv.bulk_rock[i] != 0.0 ){ //|| gv.O_id == i
			sum += 1;
		}
	}

	/** calculate fbc to be used for normalization factor of solution phases */
	z_b.fbc			= 0.0; 
	for (i = 0; i < gv.len_ox; i++){
		z_b.fbc += z_b.bulk_rock[i]*z_b.apo[i];
	}
	
	z_b.nzEl_val = sum;						/** store number of non zero values */
	z_b.zEl_val  = gv.len_ox - sum;			/** store number of zero values 	*/
	
	if (z_b.zEl_val > 0){
		j = 0; k = 0;
		for (i = 0; i < gv.len_ox; i++){
			if (gv.bulk_rock[i] == 0.){
				z_b.zEl_array[j] = i;
				j += 1;
			}
			else{
				z_b.nzEl_array[k] = i;
				k += 1;
			}
		}
	}
	else {
		for (i = 0; i < gv.len_ox; i++){
			z_b.nzEl_array[i] = i;
		}
	}

	for ( i = 0; i < z_b.nzEl_val; i++){
		z_b.bulk_rock_cat[i] = z_b.bulk_rock[z_b.nzEl_array[i]];
	}
	for ( i = z_b.nzEl_val; i < gv.len_ox; i++){
		z_b.bulk_rock_cat[i] = 0.0;
	}
	
	return z_b;
};


/**
  reset considered phases entries
*/
void reset_cp(						global_variable 	 gv,
									bulk_info 	 		z_b,
									csd_phase_set  		*cp
){
	int n 			= 16;
	for (int i = 0; i < gv.max_n_cp; i++){		
		strcpy(cp[i].name,"");						/* get phase name */	
		cp[i].in_iter			=  0;
		cp[i].split				=  0;
		cp[i].id 				= -1;				/* get phaseid */
		cp[i].n_xeos			=  0;				/* get number of compositional variables */
		cp[i].n_em				=  0;				/* get number of endmembers */
		cp[i].n_sf				=  0;			
		cp[i].df 				=  0.0;
		cp[i].factor 			=  0.0;
		cp[i].factor_norm		=  0.0;
		for (int ii = 0; ii < gv.n_flags; ii++){
			cp[i].ss_flags[ii] 	= 0;
		}

		cp[i].ss_n        		= 0.0;				/* get initial phase fraction */
		cp[i].ss_n_mol      	= 0.0;				/* get initial phase fraction */
		cp[i].ss_n_wt       	= 0.0;				/* get initial phase fraction */
		cp[i].delta_ss_n    	= 0.0;				/* get initial phase fraction */
		
		for (int ii = 0; ii < n; ii++){
			cp[i].p_em[ii]      = 0.0;
			cp[i].xi_em[ii]     = 0.0;
			cp[i].dguess[ii]    = 0.0;
			cp[i].xeos[ii]      = 0.0;
			cp[i].xeos_0[ii]    = 0.0;
			cp[i].xeos_1[ii]    = 0.0;
			cp[i].xeos_r[ii]    = 0.0;
			cp[i].delta_mu[ii]  = 0.0;
			cp[i].dfx[ii]       = 0.0;
			cp[i].mu[ii]        = 0.0;
			cp[i].gbase[ii]     = 0.0;
			cp[i].ss_comp[ii]   = 0.0;
			cp[i].ss_comp_mol[ii]  = 0.0;
			cp[i].ss_comp_wt[ii]   = 0.0;
		}
		 
		for (int ii = 0; ii < n*2; ii++){
			cp[i].sf[ii]    	= 0.0;
		}
		cp[i].mass 				= 0.0;
		cp[i].volume 			= 0.0;
		cp[i].phase_density 	= 0.0;
		cp[i].phase_cp 			= 0.0;
	}

};

/**
  reset compositional variables (xeos) when something goes wrong during minimization 
*/
void reset_SS(						global_variable 	 gv,
									bulk_info 	 		 z_b,
									SS_ref 				*SS_ref_db
){
	/* reset solution phases */
	for (int iss = 0; iss < gv.len_ss; iss++){

		for (int j = 0; j < gv.n_flags; j++){	
			SS_ref_db[iss].ss_flags[j]   = 0;
		}

		SS_ref_db[iss].tot_pc[0] = 0;
		SS_ref_db[iss].id_pc[0]  = 0;
		for (int j = 0; j < gv.len_ss*4; j++){
			SS_ref_db[iss].solvus_id[j] = -1;	
		}

		/* reset levelling pseudocompounds */
		for (int i = 0; i < (SS_ref_db[iss].tot_pc[0] ); i++){
			SS_ref_db[iss].info[i]   = 0;
			SS_ref_db[iss].G_pc[i]   = 0.0;
			SS_ref_db[iss].DF_pc[i]  = 0.0;
			for (int j = 0; j < gv.len_ox; j++){
				SS_ref_db[iss].comp_pc[i][j]  = 0.0;
			}
			for (int j = 0; j < SS_ref_db[iss].n_em; j++){
				SS_ref_db[iss].p_pc[i][j]  = 0.0;	
			}
			for (int j = 0; j < (SS_ref_db[iss].n_xeos); j++){
				SS_ref_db[iss].xeos_pc[i][j]  = 0.0;
			}
			SS_ref_db[iss].factor_pc[i] = 0.0;
		}

		/* reset LP part of PGE (algo 2.0) */
		SS_ref_db[iss].tot_Ppc 	= 0;
		SS_ref_db[iss].id_Ppc  	= 0;
		for (int i = 0; i < (SS_ref_db[iss].n_Ppc); i++){
			SS_ref_db[iss].info_Ppc[i]   = 0;
			SS_ref_db[iss].G_Ppc[i]      = 0.0;
			SS_ref_db[iss].DF_Ppc[i]     = 0.0;
			for (int j = 0; j < gv.len_ox; j++){
				SS_ref_db[iss].comp_Ppc[i][j]  = 0.0;
			}
			for (int j = 0; j < SS_ref_db[iss].n_em; j++){
				SS_ref_db[iss].p_Ppc[i][j]  = 0;	
				SS_ref_db[iss].mu_Ppc[i][j] = 0;	
			}
			for (int j = 0; j < (SS_ref_db[iss].n_xeos); j++){
				SS_ref_db[iss].xeos_Ppc[i][j]  = 0.0;
			}
		}

		/* reset solution phase model parameters */
		for (int j = 0; j < SS_ref_db[iss].n_em; j++){
			SS_ref_db[iss].gb_lvl[j]     = 0.0;
			SS_ref_db[iss].gbase[j]      = 0.0;
			SS_ref_db[iss].xi_em[j]      = 0.0;
			SS_ref_db[iss].d_em[j]       = 0.0;
			SS_ref_db[iss].z_em[j]       = 1.0;
			SS_ref_db[iss].mu[j] 	     = 0.0;
		}
		SS_ref_db[iss].sum_xi		     = 0.0;
		SS_ref_db[iss].df			     = 0.0;
		SS_ref_db[iss].df_raw		     = 0.0;

		for (int k = 0; k < SS_ref_db[iss].n_xeos; k++) {					/** initialize initial guess using default thermocalc guess	*/ 
			SS_ref_db[iss].iguess[k]     = 0.0;
			SS_ref_db[iss].dguess[k]     = 0.0;
			SS_ref_db[iss].mguess[k]     = 0.0;
			SS_ref_db[iss].xeos[k]       = 0.0;
			SS_ref_db[iss].bounds[k][0]  = SS_ref_db[iss].bounds_ref[k][0];
			SS_ref_db[iss].bounds[k][1]  = SS_ref_db[iss].bounds_ref[k][1];
		}

		for (int j = 0; j < SS_ref_db[iss].n_em; j++){
			SS_ref_db[iss].p[j]     = 0.0;
			SS_ref_db[iss].ape[j]   = 0.0;
		}
	}

};

/**
  function to allocate memory for simplex linear programming (A)
*/	
void init_simplex_A( 	simplex_data 		*splx_data,
						global_variable 	 gv
){
	simplex_data *d  = (simplex_data *) splx_data;

	/* prescribe tolerance parameters */
	d->dG_B_tol	   = gv.re_in_df;
	d->min_F_tol   = 1e6;
	
	/* allocate reference assemblage memory */
	d->A           = malloc ((gv.len_ox*gv.len_ox)  * sizeof(double));
	d->Alu         = malloc ((gv.len_ox*gv.len_ox)  * sizeof(double));
	d->A1  		   = malloc ((gv.len_ox*gv.len_ox)  * sizeof(double));
	
	d->ph_id_A     = malloc (gv.len_ox * sizeof(int*));
    for (int i = 0; i < gv.len_ox; i++){
		d->ph_id_A[i] = malloc ((gv.len_ox*4)  * sizeof(int));
	}
	
	d->pivot   	   = malloc ((gv.len_ox) * sizeof(int));
	d->g0_A   	   = malloc ((gv.len_ox) * sizeof(double));
	d->dG_A   	   = malloc ((gv.len_ox) * sizeof(double));
	d->n_vec   	   = malloc ((gv.len_ox) * sizeof(double));
	d->stage   	   = malloc ((gv.len_ox) * sizeof(int));

	d->gamma_ps	   = malloc ((gv.len_ox) * sizeof(double));
	d->gamma_ss	   = malloc ((gv.len_ox) * sizeof(double));
	d->gamma_tot   = malloc ((gv.len_ox) * sizeof(double));
	d->gamma_delta = malloc ((gv.len_ox) * sizeof(double));

	/* initialize arrays */
	int k;
    for (int i = 0; i < gv.len_ox; i++){
		d->gamma_tot[i] 	= 0.0;
		d->gamma_delta[i] 	= 0.0;
		d->pivot[i]	   = 0;
		d->g0_A[i]     = 0.0;
		d->dG_A[i]     = 0.0;
		d->n_vec[i]    = 0.0;
		d->gamma_ps[i] = 0.0;
		d->gamma_ss[i] = 0.0;
		
		for (int j = 0; j < gv.len_ox; j++){
			k = i + j*gv.len_ox;
			d->A[k]  = 0.0;
			d->A1[k] = 0.0;
		}
		
		for (int j = 0; j < 4; j++){
			d->ph_id_A[i][j] = 0;
		}
	}

};

/**
  function to allocate memory for simplex linear programming (B)
*/	
void init_simplex_B_em(				simplex_data 		 *splx_data,
									global_variable 	 gv

){
	simplex_data *d  = (simplex_data *) splx_data;
	
	d->ph_id_B    = malloc (3  * sizeof(int));
	d->B   	 	  = malloc ((gv.len_ox) * sizeof(double));
	d->B1   	  = malloc ((gv.len_ox) * sizeof(double));

	/** initialize arrays */
	for (int j = 0; j < 3; j++){
		d->ph_id_B[j] = 0;
	}

	for (int j = 0; j < gv.len_ox; j++){
		d->B[j]   = 0.0;	
		d->B1[j]  = 0.0;	
	}
};

/**
  function to allocate memory for simplex linear programming (A)
*/	
void reset_simplex_A( 	simplex_data 		*splx_data,
 						bulk_info 	 		 z_b,
						global_variable 	 gv
){
	simplex_data *d  = (simplex_data *) splx_data;
	
	/* allocate reference assemblage memory */
	d->ph2swp      = -1;
	d->n_swp       =  0;
	d->swp         =  0;
	d->n_Ox        =  z_b.nzEl_val;

	/* initialize arrays */
	int k;
    for (int i = 0; i < gv.len_ox; i++){
		d->gamma_tot[i] 	= 0.0;
		d->gamma_delta[i] 	= 0.0;
		d->pivot[i]	  = 0;
		d->g0_A[i]     = 0.0;
		d->dG_A[i]     = 0.0;
		d->n_vec[i]    = 0.0;
		d->stage[i]    = 0;
		d->gamma_ps[i] = 0.0;
		d->gamma_ss[i] = 0.0;
		
		for (int j = 0; j < gv.len_ox; j++){
			k = i + j*gv.len_ox;
			d->A[k]  = 0.0;
			d->Alu[k] = 0.0;
			d->A1[k] = 0.0;
		}
		
		for (int j = 0; j < 4; j++){
			d->ph_id_A[i][j] = 0;
		}
	}


};

/**
  function to allocate memory for simplex linear programming (B)
*/	
void reset_simplex_B_em(			simplex_data 		*splx_data,
									global_variable 	 gv

){
	simplex_data *d  = (simplex_data *) splx_data;

	/** initialize arrays */
	for (int j = 0; j < 3; j++){
		d->ph_id_B[j] = 0;
	}

	for (int j = 0; j < gv.len_ox; j++){
		d->B[j]   = 0.0;	
		d->B1[j]  = 0.0;	
	}
};
