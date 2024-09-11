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
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <complex.h> 

#include "../uthash.h"
#include "../MAGEMin.h"
#include "../initialize.h"
#include "../all_endmembers.h"
#include "SB_init_database.h"

oxide_data oxide_info_sb = {
	15,						/* number of endmembers */
	{"SiO2"	,"Al2O3","CaO"	,"MgO"	,"FeO"	,"K2O"	,"Na2O"	,"TiO2"	,"O"	,"MnO"	,"Cr2O3","H2O"	,"CO2"	,"S"	,"Cl"		},
	{60.08  ,101.96 ,56.08  ,40.30  ,71.85  ,94.2   ,61.98  ,79.88  ,16.0   ,70.94	,151.99 ,18.015	,44.01	, 32.06	,35.453		},
	{3.0	,5.0	,2.0	,2.0	,2.0	,3.0	,3.0	,3.0	,1.0	,2.0 	,5.0	,3.0	,3.0	, 1.0	,1.0		},
	{66.7736,108.653,42.9947,40.3262,38.7162,69.1514,61.1729,70.3246,30.5827,40.1891,106.9795,69.5449,62.8768,9.5557,33.2556	},
	{2,3,1,1,1,1,1,2,1,1,3,1,2,0,0}
	// for the standard molar entropy the values are already normalized by the reference temperature = 298.15K (25°C) and expressed in kJ
};
// SiO2:0 CaO:1 Al2O3:2 FeO:3 MgO:4 Na2O:5 */
stx11_dataset stx11_db = {
	2011,							/* Endmember default dataset number */
	6,							/* number of oxides */			
	10,							/* number of pure phases */
	14,							/* number of solution phases */
	{"SiO2"	,"CaO"	,"Al2O3","FeO"	,"MgO"	,"Na2O"																							},
	{"neph"	,"ky"	,"seif"	,"st"	,"coe"	,"qtz"	,"capv"	,"aMgO"	,"aFeO"	,"aAl2O3"														},
	{"plg"	,"sp"	,"ol"	,"wa"	,"ri"	,"opx"	,"cpx"	,"hpcpx","ak"	,"gtmj"	,"pv"	,"ppv"	,"mw"	,"cf"							},
	
	{1		,1		,1		,1		,1		,1		,1		,1		,1 		,1 		,1 		,1 		,1		,1								}, // allow solvus?
	{1521	,3554	,121	,4124	,210	,2450	,5498	,420	,3088	,381	,3412	,231	,11		,2376							}, // # of pseudocompound
	{0.249	,0.124	,0.098	,0.249	,0.049	,0.145	,0.33	,0.0499	,0.198	,0.098	,0.249	,0.049	,1.0 	,0.198							}, // discretization step in endmember fraction

	6.0, 						/** max dG under which a phase is considered to be reintroduced  					*/
	673.15,						/** max temperature above which PGE solver is active 								*/
	773.15,						/** minimum temperature above which melt is considered 								*/

	8,							/** number of inner PGE iterations, this has to be made mass or dG dependent 		*/
	0.025,						/** maximum mol% phase change during one PGE iteration in wt% 						*/
	2.5,						/** maximum delta_G of reference change during PGE 									*/
	1.0,						/** maximum update factor during PGE under-relax < 0.0, over-relax > 0.0 	 		*/

	2e-1,						/** merge instances of solution phase if norm < val 								*/
	1e-4,						/** fraction of solution phase when re-introduced 									*/
	1e-6						/** objective function tolerance 				 									*/
};

/* Function to allocate the memory of the data to be used/saved during PGE iterations */
global_variable global_variable_SB_init( 	global_variable  	 gv,
											bulk_info 			*z_b 	){
	int i, j;

	if (gv.EM_database == 0){
		stx11_dataset db 	= stx11_db;
		gv.len_pp   		= db.n_pp;		
		gv.len_ss  			= db.n_ss;
		gv.len_ox  			= db.n_ox;

		gv.PC_df_add		= db.PC_df_add;					/** min value of df under which the PC is added 									*/
		gv.solver_switch_T  = db.solver_switch_T;
		gv.min_melt_T       = db.min_melt_T;				/** minimum temperature above which melt is considered 								*/

		gv.inner_PGE_ite    = db.inner_PGE_ite;				/** number of inner PGE iterations, this has to be made mass or dG dependent 		*/
		gv.max_n_phase  	= db.max_n_phase;				/** maximum mol% phase change during one PGE iteration in wt% 						*/
		gv.max_g_phase  	= db.max_g_phase;				/** maximum delta_G of reference change during PGE 									*/
		gv.max_fac          = db.max_fac;					/** maximum update factor during PGE under-relax < 0.0, over-relax > 0.0 	 		*/

		gv.merge_value		= db.merge_value;				/** merge instances of solution phase if norm < val 								*/
		gv.re_in_n          = db.re_in_n;					/** fraction of phase when being reintroduce.  										*/
		gv.obj_tol 			= db.obj_tol;

		gv.ox 				= malloc (gv.len_ox * sizeof(char*)		);
		for (i = 0; i < gv.len_ox; i++){
			gv.ox[i] 		= malloc(20 * sizeof(char));	
			strcpy(gv.ox[i],db.ox[i]);
		}

		gv.PP_list 			= malloc (gv.len_pp * sizeof(char*)		);
		for (i = 0; i < (gv.len_pp); i++){	
			gv.PP_list[i] 	= malloc(20 * sizeof(char));
			strcpy(gv.PP_list[i],db.PP[i]);
		}

		gv.SS_list 			= malloc ((gv.len_ss) * sizeof (char*)	);
		gv.n_SS_PC     		= malloc ((gv.len_ss) * sizeof (int) 	);
		gv.verifyPC  		= malloc ((gv.len_ss) * sizeof (int) 	);
		gv.SS_PC_stp     	= malloc ((gv.len_ss) * sizeof (double) );
		for (i = 0; i < gv.len_ss; i++){ 
			gv.SS_list[i] 	= malloc(20 * sizeof(char)				);
			strcpy(gv.SS_list[i],db.SS[i]);
			gv.verifyPC[i]  = db.verifyPC[i]; 
			gv.n_SS_PC[i] 	= db.n_SS_PC[i]; 
			gv.SS_PC_stp[i] = db.SS_PC_stp[i]; 	
		}
	}

	/**
	 ALLOCATE MEMORY OF OTHER GLOBAL VARIABLES
	*/
	gv.n_min     		= malloc ((gv.len_ss) * sizeof (int) 	);
	gv.n_ss_ph  		= malloc ((gv.len_ss) * sizeof (int) 	);
	gv.bulk_rock 		= malloc (gv.len_ox * sizeof(double)	);
	gv.PGE_mass_norm  	= malloc (gv.it_f*2 * sizeof (double) 	); 
	gv.Alg  			= malloc (gv.it_f*2 * sizeof (int) 		); 
	gv.gamma_norm  		= malloc (gv.it_f*2 * sizeof (double) 	); 
	gv.gibbs_ev  		= malloc (gv.it_f*2 * sizeof (double) 	); 
	gv.ite_time  		= malloc (gv.it_f*2 * sizeof (double) 	); 
	
	/* store values for numerical differentiation */
	/* The last entries MUST be [0-1][end] = 0.0  */
	gv.n_Diff = 11;
	gv.pdev = malloc (2 * sizeof(double*));			
	for (i = 0; i < 2; i++){
		gv.pdev[i] = malloc (gv.n_Diff * sizeof(double));
	}
	gv.pdev[0][0]  =  0.0;	gv.pdev[1][0]  =  1.0;	
	gv.pdev[0][1]  =  0.0;	gv.pdev[1][1]  = -1.0;	
	gv.pdev[0][2]  =  1.0;	gv.pdev[1][2]  =  1.0;	
	gv.pdev[0][3]  =  1.0;	gv.pdev[1][3]  = -1.0;	
	gv.pdev[0][4]  =  2.0;	gv.pdev[1][4]  =  0.0;
	gv.pdev[0][5]  =  1.0;	gv.pdev[1][5]  =  0.0;	
	gv.pdev[0][6]  =  0.0;	gv.pdev[1][6]  =  0.0;
	gv.pdev[0][7]  =  3.0;	gv.pdev[1][7]  =  0.0;
	gv.pdev[0][8]  =  1.0;	gv.pdev[1][8]  =  0.0;
	gv.pdev[0][9]  =  0.0;	gv.pdev[1][9]  =  0.0;
	gv.pdev[0][10] =  0.0;	gv.pdev[1][10] =  0.0;

	gv.V_cor = malloc (2 * sizeof(double));

	/* declare size of chemical potential (gamma) vector */
	gv.dGamma 			= malloc (gv.len_ox * sizeof(double)	);
	gv.gam_tot  		= malloc (gv.len_ox * sizeof (double) 	); 
	gv.gam_tot_0		= malloc (gv.len_ox * sizeof (double) 	); 
	gv.delta_gam_tot  	= malloc (gv.len_ox * sizeof (double) 	); 
	gv.mass_residual 	= malloc (gv.len_ox * sizeof(double)	);	

	gv.lwork 			= 64;
	gv.ipiv     		= malloc ((gv.len_ox*3) * sizeof (int) 	);
	gv.work     		= malloc ((gv.len_ox*gv.lwork) * sizeof (double) 	);
	gv.n_solvi			= malloc ((gv.len_ss) * sizeof (int) 	);

	/* size of the flag array */
	gv.n_flags     = 5;

	/* allocate memory for pure and solution phase fractions */
	gv.pp_n    			= malloc (gv.len_pp * sizeof(double)	);									/** pure phase fraction vector */
	gv.pp_n_mol 		= malloc (gv.len_pp * sizeof(double)	);									/** pure phase fraction vector */
	gv.pp_xi    		= malloc (gv.len_pp * sizeof(double)	);									/** pure phase fraction vector */
	gv.delta_pp_n 		= malloc (gv.len_pp * sizeof(double)	);									/** pure phase fraction vector */
	gv.delta_pp_xi 		= malloc (gv.len_pp * sizeof(double)	);									/** pure phase fraction vector */
	gv.pp_flags 		= malloc (gv.len_pp * sizeof(int*)		);

	for (i = 0; i < (gv.len_pp); i++){	
		gv.pp_flags[i]   	= malloc (gv.n_flags  * sizeof(int));
	}
		
	/**
		PGE Matrix and RHS
	*/
	/* PGE method matrix and gradient arrays */
	gv.A_PGE  = malloc ((gv.len_ox*gv.len_ox*9) 	* sizeof(double));			
	gv.A0_PGE = malloc ((gv.len_ox*gv.len_ox*9) 	* sizeof(double));			
	gv.b_PGE  = malloc ((gv.len_ox*3) 				* sizeof(double));			

	gv.cp_id  = malloc ((gv.len_ox) 				* sizeof(int)	);			
	gv.pp_id  = malloc ((gv.len_ox) 				* sizeof(int)	);			

	gv.dn_cp  = malloc ((gv.len_ox) 				* sizeof(double));			
	gv.dn_pp  = malloc ((gv.len_ox) 				* sizeof(double));			

	/* stoechiometry matrix */
	gv.A  = malloc ((gv.len_ox) * sizeof(double*));		
	gv.A2 = malloc ((gv.len_ox) * sizeof(double*));			
	for (i = 0; i < (gv.len_ox); i++){
		gv.A[i]  = malloc ((gv.len_ox) * sizeof(double));
		gv.A2[i] = malloc ((gv.len_ox) * sizeof(double));
	
	gv.pc_id= malloc (gv.len_ox * sizeof(int));}
	gv.b 	= malloc (gv.len_ox * sizeof(double));	
	gv.b1 	= malloc (gv.len_ox * sizeof(double));	
	gv.tmp1 = malloc (gv.len_ox * sizeof(double));	
	gv.tmp2 = malloc (gv.len_ox * sizeof(double));	
	gv.tmp3 = malloc (gv.len_ox * sizeof(double));	
	gv.n_ss_array = malloc (gv.len_ss* sizeof(double));	
	/** 
		allocate oxides informations 						
	*/	
	z_b->apo     		= malloc (gv.len_ox * sizeof (double) ); 
	z_b->masspo     	= malloc (gv.len_ox * sizeof (double) );
	z_b->ElEntropy     	= malloc (gv.len_ox * sizeof (double) );
	z_b->id     		= malloc (gv.len_ox * sizeof (int) 	  );

	/**
		retrieve the right set of oxide and their informations 
	*/
	oxide_data ox_in 	= oxide_info_sb;
	for (i = 0; i < gv.len_ox; i++){
		for (j = 0; j < ox_in.n_ox; j++){
			if (strcmp( gv.ox[i], ox_in.oxName[j]) == 0){
				if (strcmp( gv.ox[i], "Al2O3") 	== 0){
					gv.Al2O3_id = i;
				}
				else if (strcmp( gv.ox[i], "TiO2") 	== 0){
					gv.TiO2_id = i;
				}
				else if (strcmp( gv.ox[i], "O") 	== 0){
					gv.O_id = i;
				}											
				z_b->apo[i]     	= ox_in.atPerOx[j];
				z_b->masspo[i]  	= ox_in.oxMass[j];
				z_b->ElEntropy[i]   = ox_in.ElEntropy[j];
				z_b->id[i]  		= j;
				break;
			}
		}
	}

	z_b->bulk_rock_cat  = malloc (gv.len_ox * sizeof (double) ); 
	z_b->bulk_rock  	= malloc (gv.len_ox * sizeof (double) ); 
	z_b->nzEl_array 	= malloc (gv.len_ox * sizeof (int) ); 
	z_b->zEl_array 		= malloc (gv.len_ox * sizeof (int) ); 

	/* sets end-member dataset information */
	if (gv.EM_dataset == 2011){
			gv.n_em_db 			= 46;
	}
	// else if (gv.EM_dataset == 22){
	// 		gv.n_em_db 			= 289;
	// }
		

	return gv;
}

/* Provide a list of test bulk-rock composition for the metapelite database (White et al., 2014)*/
global_variable get_bulk_stx11( global_variable gv) {
 	if (gv.test != -1){
		if (gv.verbose == 1){
			printf("\n");
			printf("   - Minimization using in-built bulk-rock  : test %2d\n",gv.test);	
		}							
	}
	else{
		gv.test = 0;
		if (gv.verbose == 1){
			printf("\n");
			printf("   - No predefined bulk provided -> user custom bulk (if none provided, will run default KLB1)\n");	
		}	
	}
	if (gv.test == 0){ //KLB1
		/* SiO2 Al2O3 CaO MgO FeO K2O Na2O TiO2 O Cr2O3 H2O */
		/* Bulk rock composition of Peridotite from Holland et al., 2018, given by E. Green */
		gv.bulk_rock[0]  = 38.494 ;		/** SiO2 	*/
		gv.bulk_rock[1]  = 2.824;		/** CaO  	*/
		gv.bulk_rock[2]  = 1.776;		/** Al2O2 	*/
		gv.bulk_rock[3]  = 50.566;		/** MgO 	*/
		gv.bulk_rock[4]  = 5.886;		/** FeO 	*/
		gv.bulk_rock[5]  = 0.250;		/** Na2O 	*/	
	}	
	else{
		printf("Unknown test %i - please specify a different test! \n", gv.test);
	 	exit(EXIT_FAILURE);
	}
	return gv;
}
