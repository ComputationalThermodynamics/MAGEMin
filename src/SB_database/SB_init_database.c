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
	16,						/* number of endmembers */
    {"SiO2"	,"CaO"	,"Al2O3","FeO"	,"MgO"	,"Na2O"	,"K2O"	,"TiO2"	,"O"	,"MnO"	,"Cr2O3","H2O"	,"CO2"	,"S"	,"Cl", "Fe"},
    {"Si"	,"Ca"	,"Al"	,"Fe"	,"Mg"	,"Na"	,"K"	,"Ti"	,"O"	,"Mn"	,"Cr"	,"H"	,"C"	,"S"	,"Cl", "Fe"},
    {60.08  ,56.08  ,101.96 ,71.85  ,40.30  ,61.98  ,94.2   ,79.88  ,16.0   ,70.94	,151.99 ,18.015	,44.01	, 32.06	,35.453, 55.845},
    {3.0	,2.0	,5.0	,2.0	,2.0	,3.0	,3.0	,3.0	,1.0	,2.0 	,5.0	,3.0	,3.0	, 1.0	,1.0,   1.0},
    {66.7736,42.9947,108.653,38.7162,40.3262,61.1729,69.1514,70.3246,30.5827,40.1891,106.9795,69.5449,62.8768,9.5557,33.2556,1.0},
    {2,1,3,1,1,1,1,2,1,1,3,1,2,0,0,0},
    {1,1,2,1,1,2,2,1,1,1,2,2,1,1,1,1}
	// for the standard molar entropy the values are already normalized by the reference temperature = 298.15K (25°C) and expressed in kJ
};

// SiO2:0 CaO:1 Al2O3:2 FeO:3 MgO:4 Na2O:5 */
stx11_dataset stx11_db = {
	2011,						/* Endmember default dataset number */
	6,							/* number of oxides */			
	10,							/* number of pure phases */
	14,							/* number of solution phases */
	{"SiO2"	,"CaO"	,"Al2O3","FeO"	,"MgO"	,"Na2O"																							},
	{"neph"	,"ky"	,"st"	,"coe"	,"qtz"	,"capv"	,"co" 	,"aMgO"	,"aFeO"	,"aAl2O3"																},
	{"plg"	,"sp"	,"ol"	,"wa"	,"ri"	,"opx"	,"cpx"	,"hpcpx","ak"	,"gtmj"	,"pv"	,"ppv"	,"mw"	,"cf"							},
	
	{1		,1		,1		,1		,1		,1		,1		,1		,1 		,1 		,1 		,1 		,1		,1								}, // allow solvus?
	{11 	,11 	,11 	,11 	,11 	,286	,1001	,11 	,66 	,1001	,66 	,66 	,11		,66 							}, // # of pseudocompound
	{0.050	,0.050	,0.050	,0.050	,0.050	,0.050	,0.050	,0.050	,0.050	,0.050	,0.050	,0.050	,0.050 	,0.050							}, // discretization step in endmember fraction

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

// SiO2:0 CaO:1 Al2O3:2 FeO:3 MgO:4 Na2O:5 */
stx21_dataset stx21_db = {
	2021,						/* Endmember default dataset number */
	6,							/* number of oxides */			
	10,							/* number of pure phases */
	15,							/* number of solution phases */
	{"SiO2"	,"CaO"	,"Al2O3","FeO"	,"MgO"	,"Na2O"																							},
	{"neph"	,"ky"	,"st"	,"coe"	,"qtz"	,"capv"	,"co" 	,"aMgO"	,"aFeO"	,"aAl2O3"														},
	{"plg"	,"sp"	,"ol"	,"wa"	,"ri"	,"opx"	,"cpx"	,"hpcpx","ak"	,"gtmj"	,"pv"	,"ppv"	,"cf"	,"mw"	,"nal"					},
	
	{1		,1		,1		,1		,1		,1		,1		,1		,1 		,1 		,1 		,1 		,1		,1		,1						}, // allow solvus?
	{11 	,11 	,11 	,11 	,11 	,286	,1001	,11 	,66 	,1001	,66 	,66 	,66		,66 	,66						}, // # of pseudocompound
	{0.050	,0.050	,0.050	,0.050	,0.050	,0.050	,0.050	,0.050	,0.050	,0.050	,0.050	,0.050	,0.050 	,0.050	,0.050					}, // discretization step in endmember fraction

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

// SiO2:0 CaO:1 Al2O3:2 Fe:3 MgO:4 Na2O:5 Cr2O3:6 O:7 */
stx24_dataset stx24_db = {
	2024,						/* Endmember default dataset number */
	8,							/* number of oxides */			
	17,							/* number of pure phases */
	15,							/* number of solution phases */
	{"SiO2"	,"CaO"	,"Al2O3", "MgO"	,"Na2O" , "O"   ,"Cr2O3", "Fe"																		},
	{"neph"	,"ky"	,"st"	,"coe"	,"qtz"	,"capv"	, "O2" 	,  "fea",  "fee",  "feg", "apbo", "wo" , "lppv" , "pwo", "aMgO" ,"aFeO" ,"aAl2O3"					},
	{"plg"	,"sp"	,"ol"	,"wa"	,"ri"	,"opx"	,"cpx"	,"hpcpx", "ak"	,"gtmj"	,"pv"	 ,"ppv"	 ,"cf"	 ,"mw"	  ,"nal"					},
	
	{1		,1		,1		,1		,1		,1		,1		,1		,1 		,1 		,1 		,1 		,1		,1		,1						}, // allow solvus?
	{11 	,286 	,11 	,11 	,11 	,286	,3003	,11 	,1001 	,3003	,3003 	,1001 	,1001	,1001 	,66						}, // # of pseudocompound
	{0.050	,0.050	,0.050	,0.050	,0.050	,0.050	,0.050	,0.050	,0.050	,0.125	,0.125	,0.050	,0.050 	,0.050	,0.050					}, // discretization step in endmember fraction

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
		gv.EM_dataset 		= db.ds_version;	
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
		gv.act_PP     		= malloc ((gv.len_pp) * sizeof (int) 	);
		for (i = 0; i < gv.len_pp; i++){ gv.act_PP[i] = 1; 			};
	}
	else if (gv.EM_database == 1){
		stx21_dataset db 	= stx21_db;
		gv.EM_dataset 		= db.ds_version;	
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
		gv.act_PP     		= malloc ((gv.len_pp) * sizeof (int) 	);
		for (i = 0; i < gv.len_pp; i++){ gv.act_PP[i] = 1; 			};
	}
	else if (gv.EM_database == 2){
		stx24_dataset db 	= stx24_db;
		gv.EM_dataset 		= db.ds_version;	
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
		gv.act_PP     		= malloc ((gv.len_pp) * sizeof (int) 	);
		for (i = 0; i < gv.len_pp; i++){ gv.act_PP[i] = 1; 			};
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
	gv.n_Diff = 8;
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
	gv.pdev[0][6]  =  3.0;	gv.pdev[1][6]  =  0.0;
	gv.pdev[0][7]  =  0.0;	gv.pdev[1][7]  =  0.0;

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
	gv.pp_n_wt  		= malloc (gv.len_pp * sizeof(double)	);									/** pure phase fraction vector */
	gv.pp_n_vol  		= malloc (gv.len_pp * sizeof(double)	);
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
	z_b->opo     		= malloc (gv.len_ox * sizeof (double) );
	z_b->cpo     		= malloc (gv.len_ox * sizeof (double) );
	z_b->ElEntropy     	= malloc (gv.len_ox * sizeof (double) );
	z_b->id     		= malloc (gv.len_ox * sizeof (int) 	  );
	z_b->elName     	= malloc (gv.len_ox * sizeof (char*) );
	for (i = 0; i < (gv.len_ox); i++){	
		z_b->elName[i] 	= malloc(20 * sizeof(char));
	}
	/**
		retrieve the right set of oxide and their informations 
	*/
	gv.Al2O3_id = -1;
	gv.TiO2_id  = -1;
	gv.CaO_id 	= -1;
	gv.Na2O_id 	= -1;
	gv.FeO_id 	= -1;
	gv.Fe_id 	= -1;
	gv.MgO_id 	= -1;
	gv.Cr2O3_id = -1;
	gv.O_id 	= -1;
	oxide_data ox_in 	= oxide_info_sb;
	for (i = 0; i < gv.len_ox; i++){
		for (j = 0; j < ox_in.n_ox; j++){
			if (strcmp( gv.ox[i], ox_in.oxName[j]) == 0){
				if (strcmp( gv.ox[i], "Al2O3") 	== 0){
					gv.Al2O3_id = i;
				}
				// else if (strcmp( gv.ox[i], "TiO2") 	== 0){
				// 	gv.TiO2_id = i;
				// }
				else if (strcmp( gv.ox[i], "O") 	== 0){
					gv.O_id = i;
				}	
				else if (strcmp( gv.ox[i], "CaO") 	== 0){
					gv.CaO_id = i;
				}
				else if (strcmp( gv.ox[i], "Na2O") 	== 0){
					gv.Na2O_id = i;
				}
				else if (strcmp( gv.ox[i], "MgO") 	== 0){
					gv.MgO_id = i;
				}	
				else if (strcmp( gv.ox[i], "FeO") 	== 0){
					gv.FeO_id = i;
				}			
				else if (strcmp( gv.ox[i], "Fe") 	== 0){
					gv.Fe_id = i;
				}					
				else if (strcmp( gv.ox[i], "Cr2O3") == 0){
					gv.Cr2O3_id = i;
				}		
				z_b->apo[i]     	= ox_in.atPerOx[j];
				z_b->masspo[i]  	= ox_in.oxMass[j];
				z_b->opo[i]  		= ox_in.OPerOx[j];
				z_b->cpo[i]  		= ox_in.catPerOx[j];
				z_b->ElEntropy[i]   = ox_in.ElEntropy[j];
				strcpy(z_b->elName[i],ox_in.elName[j]);
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
			gv.n_em_db 		= 46;
	}
	else if (gv.EM_dataset == 2021){
		gv.n_em_db 			= 68;
	}
	else if (gv.EM_dataset == 2024){
		gv.n_em_db 			= 75;
	}

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
		/* SiO2 CaO  Al2O3 FeO MgO  Na2O  */
		/* Bulk rock composition of Peridotite from Holland et al., 2018, given by E. Green */
		gv.bulk_rock[0]  = 38.41;		/** SiO2 	*/
		gv.bulk_rock[1]  = 3.18;		/** CaO  	*/
		gv.bulk_rock[2]  = 1.8;			/** Al2O2 	*/
		gv.bulk_rock[3]  = 5.85;		/** FeO 	*/
		gv.bulk_rock[4]  = 50.49;		/** MgO 	*/
		gv.bulk_rock[5]  = 0.250;		/** Na2O 	*/	
	}	
	else if (gv.test == 1){ //Pyrolite
		/* SiO2 CaO  Al2O3 FeO MgO  Na2O  */
		gv.bulk_rock[0]  = 38.89 ;		/** SiO2 	*/
		gv.bulk_rock[1]  = 3.1;			/** CaO  	*/
		gv.bulk_rock[2]  = 2.2;			/** Al2O2 	*/
		gv.bulk_rock[3]  = 5.8;			/** FeO 	*/
		gv.bulk_rock[4]  = 50.0;		/** MgO 	*/
		gv.bulk_rock[5]  = 0.01;		/** Na2O 	*/
	}  
	else if (gv.test == 2){ //Harzburgite
		/* SiO2 CaO  Al2O3 FeO MgO  Na2O  */
		gv.bulk_rock[0]  = 36.39 ;		/** SiO2 	*/
		gv.bulk_rock[1]  = 0.9;			/** CaO  	*/
		gv.bulk_rock[2]  = 0.7;			/** Al2O2 	*/
		gv.bulk_rock[3]  = 5.4;			/** FeO 	*/
		gv.bulk_rock[4]  = 56.6;		/** MgO 	*/
		gv.bulk_rock[5]  = 0.01;		/** Na2O 	*/
	}  
	else{
		printf("Unknown test %i - please specify a different test! \n", gv.test);
	 	exit(EXIT_FAILURE);
	}
	return gv;
}

/* Provide a list of test bulk-rock composition for the metapelite database (White et al., 2014)*/
global_variable get_bulk_stx21( global_variable gv) {
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
	   /* SiO2 CaO  Al2O3 FeO MgO  Na2O  */
	   /* Bulk rock composition of Peridotite from Holland et al., 2018, given by E. Green */
	   gv.bulk_rock[0]  = 38.41;		/** SiO2 	*/
	   gv.bulk_rock[1]  = 3.18;		/** CaO  	*/
	   gv.bulk_rock[2]  = 1.8;			/** Al2O2 	*/
	   gv.bulk_rock[3]  = 5.85;		/** FeO 	*/
	   gv.bulk_rock[4]  = 50.49;		/** MgO 	*/
	   gv.bulk_rock[5]  = 0.250;		/** Na2O 	*/	
   }	
   else if (gv.test == 1){ //Pyrolite
	   /* SiO2 CaO  Al2O3 FeO MgO  Na2O  */
	   gv.bulk_rock[0]  = 38.89 ;		/** SiO2 	*/
	   gv.bulk_rock[1]  = 3.1;			/** CaO  	*/
	   gv.bulk_rock[2]  = 2.2;			/** Al2O2 	*/
	   gv.bulk_rock[3]  = 5.8;			/** FeO 	*/
	   gv.bulk_rock[4]  = 50.0;		/** MgO 	*/
	   gv.bulk_rock[5]  = 0.01;		/** Na2O 	*/
   }  
   else if (gv.test == 1){ //Harzburgite
	   /* SiO2 CaO  Al2O3 FeO MgO  Na2O  */
	   gv.bulk_rock[0]  = 36.39 ;		/** SiO2 	*/
	   gv.bulk_rock[1]  = 0.9;			/** CaO  	*/
	   gv.bulk_rock[2]  = 0.7;			/** Al2O2 	*/
	   gv.bulk_rock[3]  = 5.4;			/** FeO 	*/
	   gv.bulk_rock[4]  = 56.6;		/** MgO 	*/
	   gv.bulk_rock[5]  = 0.01;		/** Na2O 	*/
   }  
   else{
	   printf("Unknown test %i - please specify a different test! \n", gv.test);
		exit(EXIT_FAILURE);
   }
   return gv;
}

/* Provide a list of test bulk-rock composition for the metapelite database (White et al., 2014)*/
global_variable get_bulk_stx24( global_variable gv) {
	if (gv.test != -1){
	   if (gv.verbose == 1){
		   printf("\n");
		   printf("   - Minimization using in-built bulk-rock  : test %2d\n",gv.test);	
	   }							
   }
   else{
	   gv.test = 1;
	   if (gv.verbose == 1){
		   printf("\n");
		   printf("   - No predefined bulk provided -> user custom bulk (if none provided, will run default Pyrolite)\n");	
	   }	
   }

   if (gv.test == 0){ //Basalt
	   /* SiO2 CaO  Al2O3 Fe MgO Na2O Cr2O3 O  */
	   /* Bulk rock composition of Peridotite from Holland et al., 2018, given by E. Green */
	   gv.bulk_rock[0]  = 48.30;		/** SiO2 	*/
	   gv.bulk_rock[1]  = 12.87;		/** CaO  	*/
	   gv.bulk_rock[2]  = 9.48;		/** Al2O3 	*/
	   gv.bulk_rock[3]  = 13.96;		/** MgO 	*/
	   gv.bulk_rock[4]  = 2.07;			/** Na2O 	*/
	   gv.bulk_rock[5]  = 6.84;			/** O 		*/
	   gv.bulk_rock[6]  = 0.03;		/** Cr2O3 	*/
	   gv.bulk_rock[7]  = 6.46;			/** Fe 		*/
   }	
   else if (gv.test == 1){ //Pyrolite
	   /* SiO2 CaO  Al2O3 Fe MgO Na2O Cr2O3 O  */
	   gv.bulk_rock[0]  = 38.83;		/** SiO2 	*/
	   gv.bulk_rock[1]  = 2.94;			/** CaO  	*/
	   gv.bulk_rock[2]  = 2.03;			/** Al2O3 	*/
	   gv.bulk_rock[3]  = 50.02;		/** MgO 	*/
	   gv.bulk_rock[4]  = 0.11;			/** Na2O 	*/
	   gv.bulk_rock[5]  = 5.87;			/** O 		*/
	   gv.bulk_rock[6]  = 0.19;			/** Cr2O3 	*/
	   gv.bulk_rock[7]  = 5.81;			/** Fe 		*/
   }  
   else if (gv.test == 2){ //Harzburgite
	   /* SiO2 CaO  Al2O3 Fe MgO Na2O Cr2O3 O  */
	   gv.bulk_rock[0]  = 34.01;		/** SiO2 	*/
	   gv.bulk_rock[1]  = 0.76;			/** CaO  	*/
	   gv.bulk_rock[2]  = 0.46;			/** Al2O3 	*/
	   gv.bulk_rock[3]  = 53.62;		/** MgO 	*/
	   gv.bulk_rock[4]  = 0.008;		/** Na2O 	*/
	   gv.bulk_rock[5]  = 5.54;		/** O 		*/
	   gv.bulk_rock[6]  = 0.09;		/** Cr2O3 	*/
	   gv.bulk_rock[7]  = 5.52;		/** Fe 		*/
   }  
   else{
	   printf("Unknown test %i - please specify a different test! \n", gv.test);
		exit(EXIT_FAILURE);
   }
  
   return gv;
}

