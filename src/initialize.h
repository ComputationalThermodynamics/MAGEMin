#ifndef __INITIALIZE_H_
#define __INITIALIZE_H_

#include "uthash.h"
#include "Endmembers_tc-ds634.h"

enum {
	_tc_ds634_
};

/* Select required thermodynamic database */
struct EM_db Access_EM_DB(int id, int EM_database) {
	struct EM_db Entry_EM;

	//if (EM_database == _tc_ds633_){	
	Entry_EM = arr_em_db_tc_ds634[id]; 
	//}
	
	return Entry_EM;
}

/*---------------------------------------------------------------------------*/ 
/*  Hashtable for endmember in thermodynamic database                        */
struct EM2id {
    char EM_tag[20];           /* key (string is WITHIN the structure)       */
    int id;                    /* id of the key (array index)                */
    UT_hash_handle hh;         /* makes this structure hashable              */
};
struct EM2id *EM = NULL;

int find_EM_id(char *EM_tag) {
	struct EM2id *p_s;
	HASH_FIND_STR ( EM, EM_tag, p_s );
	return(p_s->id);	
}

/*  Hashtable for Pure-phases in the pure-phases list                        */
struct PP2id {
    char PP_tag[20];           /* key (string is WITHIN the structure)       */
    int id;                    /* id of the key (array index)                */
    UT_hash_handle hh;         /* makes this structure hashable              */
};
struct PP2id *PP = NULL;

int find_PP_id(char *PP_tag) {
	struct PP2id *pp_s;
	HASH_FIND_STR ( PP, PP_tag, pp_s );
	return(pp_s->id);	
}

/*  Function to retrieve the endmember names from the database               */
/*  Note the size of the array is n_em_db+1, required for the hashtable      */
char** get_EM_DB_names(int EM_database) {
	struct EM_db EM_return;
	int i;
	char ** names = malloc((n_em_db+1) * sizeof(char*));
	for ( i = 0; i < n_em_db; i++){
		names[i] = malloc(20 * sizeof(char));
	}
	for ( i = 0; i < n_em_db; i++){	
		EM_return = Access_EM_DB(i, EM_database);
		strcpy(names[i],EM_return.Name);
	}
	return names;
}

/* get position of zeros and non-zeros values in the bulk */
bulk_info initialize_bulk_infos(			double  P, 
											double  T			){
	bulk_info z_b;
	
	int i, j, k;
	z_b.P 			= P;
	z_b.T 			= T;
	z_b.R 			= 0.0083144;
	
	z_b.apo     	= malloc (nEl * sizeof (double) ); 
	z_b.apo[0]  	= 3.0;
	z_b.apo[1]  	= 5.0;
	z_b.apo[2]  	= 2.0;
	z_b.apo[3]  	= 2.0;
	z_b.apo[4]  	= 2.0;
	z_b.apo[5]  	= 3.0;
	z_b.apo[6]  	= 3.0;
	z_b.apo[7]  	= 3.0;
	z_b.apo[8]  	= 1.0;
	z_b.apo[9]  	= 5.0;
	z_b.apo[10] 	= 3.0;
	
	z_b.masspo     	= malloc (nEl * sizeof (double) );
	z_b.masspo[0]  	= 60.08;
	z_b.masspo[1]  	= 101.96;
	z_b.masspo[2]  	= 56.08;
	z_b.masspo[3]  	= 40.30;
	z_b.masspo[4]  	= 71.85;
	z_b.masspo[5]  	= 94.2;
	z_b.masspo[6]  	= 61.98;
	z_b.masspo[7]  	= 79.88;
	z_b.masspo[8]  	= 16.0;
	z_b.masspo[9]  	= 151.99;
	z_b.masspo[10] 	= 18.015;

	z_b.bulk_rock_cat  	= malloc (nEl * sizeof (double) ); 
	z_b.bulk_rock  		= malloc (nEl * sizeof (double) ); 

	return z_b;
}

/* Function to allocate the memory of the data to be used/saved during PGE iterations */
global_variable global_variable_init(){
	global_variable gv;
	
	gv.outpath 			= malloc(100 * sizeof(char));
	gv.version 			= malloc(50  * sizeof(char));
	gv.len_pp      		= 10;	
	
	/* Control center... */
	gv.save_residual_evolution = 1;				/** verbose needs to be set to 0 to save the residual evolution 					*/

	/* oxides and solution phases */
	char   *ox_tmp[] 		= {"SiO2"	,"Al2O3","CaO"	,"MgO"	,"FeO"	,"K2O"	,"Na2O"	,"TiO2"	,"O"	,"Cr2O3","H2O"								};
	char   *PP_tmp[] 		= {"q"		,"crst"	,"trd"	,"coe"	,"stv"	,"ky"	,"sill"	,"and"	,"ru"	,"sph"										};
	char   *SS_tmp[]     	= {"spn"	,"bi"	,"cd"	,"cpx"	,"ep"	,"g"	,"hb"	,"ilm"	,"liq"	,"mu"	,"ol"	,"opx"	,"pl4T"	,"fl"		};
	/* next entry is a flag to check for wrong local minimum/solvus when getting close to solution */
	int     verifyPC_tmp[]	= {1		,1		,1		,1		,1		,1		,1		,1		,1 		,1 		,1 		,1 		,1 		,1			};
	int 	n_SS_PC_tmp[]   = {1521		,1645	,121	,4384	,110	,1224	,4950	,420	,3099	,2376	,222	,1735	,231	,1			};
	double 	SS_PC_stp_tmp[] = {0.249	,0.124	,0.098	,0.249	,0.049	,0.199	,0.249	,0.0499	,0.198	,0.198	,0.098	,0.249	,0.049	,1.0 		};

	/* system parameters */
	strcpy(gv.outpath,"./output/");				/** define the outpath to save logs and final results file	 						*/
	strcpy(gv.version,"1.2.0 [26/06/2022]");	/** MAGEMin version 																*/

	gv.len_ox           = 11;					/** number of components in the system 												*/
	gv.max_n_cp 		= 128;					/** number of considered solution phases 											*/									
	gv.verbose          = 0;					/** verbose: -1, no verbose; 0, light verbose; 1, full verbose 										*/

	/* residual tolerance */
	gv.br_max_tol       = 1.0e-5;				/** value under which the solution is accepted to satisfy the mass constraint 		*/
	
	/* under-relaxing factors */
	gv.relax_PGE_val    = 128.0;				/** restricting factor 																*/
	gv.PC_check_val1	= 1.0e-2;				/** br norm under which PC are tested for potential candidate to be added 			*/
	gv.PC_check_val2	= 1.0e-4;				/** br norm under which PC are tested for potential candidate to be added 			*/
	gv.PC_min_dist 		= 1.0;					/** factor multiplying the diagonal of the hyperbox of xeos step 					*/
	gv.PC_df_add		= 4.0;					/** min value of df under which the PC is added 									*/

	/* levelling parameters */
	gv.em2ss_shift		= 1e-6;					/** small value to shift x-eos of pure endmember from bounds after levelling 		*/
	gv.bnd_filter_pc    = 10.0;					/** value of driving force the pseudocompound is considered 						*/
	gv.n_pc				= 5000;
	gv.max_G_pc         = 5.0;					/** dG under which PC is considered after their generation		 					*/
	gv.eps_sf_pc		= 1e-10;				/** Minimum value of site fraction under which PC is rejected, 
													don't put it too high as it will conflict with bounds of x-eos					*/

	/* PGE LP pseudocompounds parameters */
	gv.n_Ppc			= 2048;

	/* solvus tolerance */
	gv.merge_value      = 1e-1;					/** max norm distance between two instances of a solution phase						*/	
	
	/* local minimizer options */
	gv.bnd_val          = 1.0e-10;				/** boundary value for x-eos 										 				*/
	gv.obj_tol			= 1e-7;
	gv.ineq_res  	 	= 0.0;
	gv.box_size_mode_1	= 0.25;					/** edge size of the xeos hyperdensity used during PGE local minimization 			*/
	gv.maxeval_mode_1   = 1024;					/** max number of evaluation of the obj function for mode 1 (PGE)					*/

	/* Partitioning Gibbs Energy */
	gv.xi_em_cor   		= 0.99;	
	gv.outter_PGE_ite   = 1;					/** minimum number of outter PGE iterations, before a solution can be accepted 		*/
	gv.inner_PGE_ite    = 8;					/** number of inner PGE iterations, this has to be made mass or dG dependent 		*/
	gv.max_n_phase  	= 0.025;				/** maximum mol% phase change during one PGE iteration in wt% 						*/
	gv.max_g_phase  	= 2.5;					/** maximum delta_G of reference change during PGE 									*/
	gv.max_fac          = 1.0;					/** maximum update factor during PGE under-relax < 0.0, over-relax > 0.0 	 		*/

	/* set of parameters to record the evolution of the norm of the mass constraint                                                 */
	gv.it_1             = 128;                  /** first critical iteration                                                        */
	gv.ur_1             = 4.;                   /** under relaxing factor on mass constraint if iteration is bigger than it_1       */
	gv.it_2             = 160;                  /** second critical iteration                                                       */
	gv.ur_2             = 8.;                   /** under relaxing factor on mass constraint if iteration is bigger than it_2       */
	gv.it_3             = 192;                  /** third critical iteration                                                        */
	gv.ur_3             = 16.;                  /** under relaxing factor on mass constraint if iteration is bigger than it_3       */
	gv.it_f             = 256;                  /** gives back failure when the number of iteration is bigger than it_f             */

	/* phase update options */
	gv.re_in_n          = 1e-3;					/** fraction of phase when being reintroduce.  										*/

	/* density calculation */
	gv.gb_P_eps			= 2e-3;					/** small value to calculate V using finite difference: V = dG/dP;					*/
	gv.gb_T_eps			= 2e-3;					/** small value to calculate V using finite difference: V = dG/dP;					*/

	/* initialize other values */
	gv.mean_sum_xi		= 1.0;
	gv.sigma_sum_xi		= 1.0;
	gv.alpha        	= gv.max_fac;				/** active under-relaxing factor 													*/
	gv.len_ss          	= (int)(sizeof(n_SS_PC_tmp) / sizeof(n_SS_PC_tmp[0] ));					/** number of solution phases taken into accounnt									*/
	gv.maxeval		    = gv.maxeval_mode_1;

	/* declare chemical system */
	gv.PGE_mass_norm  	= malloc (gv.it_f * sizeof (double) 	); 
	gv.Alg  			= malloc (gv.it_f * sizeof (int) 		); 
	gv.gamma_norm  		= malloc (gv.it_f * sizeof (double) 	); 
	gv.ite_time  		= malloc (gv.it_f * sizeof (double) 	); 
	
	/* store values for numerical differentiation */
	gv.n_Diff = 7;
	gv.numDiff = malloc (2 * sizeof(double*));			
    for (int i = 0; i < 2; i++){
		gv.numDiff[i] = malloc (gv.n_Diff * sizeof(double));
	}
	gv.numDiff[0][0] = 0.0;	gv.numDiff[0][1] = 0.0;		gv.numDiff[0][2] = 1.0;	gv.numDiff[0][3] = 1.0;		gv.numDiff[0][4] = 2.0;	gv.numDiff[0][5] = 1.0;	gv.numDiff[0][6] = 0.0;
	gv.numDiff[1][0] = 1.0;	gv.numDiff[1][1] = -1.0;	gv.numDiff[1][2] = 1.0;	gv.numDiff[1][3] = -1.0;	gv.numDiff[1][4] = 0.0;	gv.numDiff[1][5] = 0.0;	gv.numDiff[1][6] = 0.0;
	gv.V_cor = malloc (2 * sizeof(double));
	/* declare size of chemical potential (gamma) vector */
	gv.dGamma 			= malloc (gv.len_ox * sizeof(double)	);
	gv.ox 				= malloc (gv.len_ox * sizeof(char*)		);
	gv.gam_tot  		= malloc (gv.len_ox * sizeof (double) 	); 
	gv.gam_tot_0		= malloc (gv.len_ox * sizeof (double) 	); 
	gv.delta_gam_tot  	= malloc (gv.len_ox * sizeof (double) 	); 
	gv.mass_residual 	= malloc (gv.len_ox * sizeof(double)	);	

	for (int i = 0; i < gv.len_ox; i++){
		gv.ox[i] 			= malloc(20 * sizeof(char));	
		strcpy(gv.ox[i],ox_tmp[i]);	
	}

	gv.n_SS_PC     		= malloc ((gv.len_ss) * sizeof (int) 	);
	gv.verifyPC  		= malloc ((gv.len_ss) * sizeof (int) 	);
	gv.SS_PC_stp     	= malloc ((gv.len_ss) * sizeof (double) );
	gv.SS_list 			= malloc ((gv.len_ss) * sizeof (char*)	);
	gv.n_solvi			= malloc ((gv.len_ss) * sizeof (int) 	);
    gv.id_solvi 		= malloc ((gv.len_ss) * sizeof (int*)	);
    
	for (int i = 0; i < (gv.len_ss); i++){	
		gv.id_solvi[i]   	= malloc (gv.max_n_cp  * sizeof(int));
	}
	
	for (int i = 0; i < gv.len_ss; i++){ 
		gv.verifyPC[i]      = verifyPC_tmp[i]; 
		gv.n_SS_PC[i] 		= n_SS_PC_tmp[i]; 
		gv.SS_PC_stp[i] 	= SS_PC_stp_tmp[i]; 
		gv.SS_list[i] 		= malloc(20 * sizeof(char)		);
		strcpy(gv.SS_list[i],SS_tmp[i]);			
	}
	
	/* size of the flag array */
    gv.n_flags     = 6;

	/* allocate memory for pure and solution phase fractions */
	gv.pp_n    			= malloc (gv.len_pp * sizeof(double)	);									/** pure phase fraction vector */
	gv.pp_n_0 			= malloc (gv.len_pp * sizeof(double)	);									/** pure phase fraction vector */
	gv.pp_xi    		= malloc (gv.len_pp * sizeof(double)	);									/** pure phase fraction vector */
	gv.delta_pp_n 		= malloc (gv.len_pp * sizeof(double)	);									/** pure phase fraction vector */
	gv.delta_pp_xi 		= malloc (gv.len_pp * sizeof(double)	);									/** pure phase fraction vector */
	gv.PP_list 			= malloc (gv.len_pp * sizeof(char*)		);
    gv.pp_flags 		= malloc (gv.len_pp * sizeof(int*)		);

	for (int i = 0; i < (gv.len_pp); i++){	
		gv.PP_list[i] 		= malloc(20 * sizeof(char));
		strcpy(gv.PP_list[i],PP_tmp[i]);
		gv.pp_flags[i]   	= malloc (gv.n_flags  * sizeof(int));
	}
		
	/**
		PGE Matrix and RHS
	*/
	/* PGE method matrix and gradient arrays */
	gv.A_PGE  = malloc ((gv.len_ox*gv.len_ox*4) 	* sizeof(double));			
	gv.A0_PGE = malloc ((gv.len_ox*gv.len_ox*4) 	* sizeof(double));			
	gv.b_PGE  = malloc ((gv.len_ox*gv.len_ox) 		* sizeof(double));			

	gv.cp_id  = malloc ((gv.len_ox) 				* sizeof(int)	);			
	gv.pp_id  = malloc ((gv.len_ox) 				* sizeof(int)	);			

	gv.dn_cp  = malloc ((gv.len_ox) 				* sizeof(double));			
	gv.dn_pp  = malloc ((gv.len_ox) 				* sizeof(double));			

	/* stoechiometry matrix */
	gv.A = malloc ((gv.len_ox) * sizeof(double*));			
    for (int i = 0; i < (gv.len_ox); i++){
		gv.A[i] = malloc ((gv.len_ox) * sizeof(double));
	}
	
	/* bulk rock vector */
	gv.b = malloc (gv.len_ox * sizeof(double));	
	
	return gv;
}

/* Get benchmark bulk rock composition given by Holland et al., 2018*/
void get_bulk(double *bulk_rock, int test, int n_El) {
	/*
	We implement the BR compositions, based on an email from Eleanor Green: 

	1) The bulk compositions that we used can be obtained from the relevant cited papers, applying the Fe2+:Fe3+ ratio that we provide on page 883. 
	I think all the dry bulk compositions will be in this list:
		
		% ---------------------------------------------------------------------
		% SiO2    Al2O3  CaO    MgO   FeO    K2O   Na2O   TiO2    O    Cr2O3    mole%  
		% ---------------------------------------------------------------------
		38.494  1.776  2.824 50.566 5.886  0.01  0.250  0.10  0.096  0.109   % 1  Kilbourne hole Recalc 2013
		52.14   9.68  13.10  12.27  8.15   0.02  3.07   1.40   0.12   0.02   % 2  G2  Pertermann & Hirschmann 03
		51.93   9.18  11.80  13.00  9.56   0.09  2.72   1.28   0.15   0.04   % 3  Fujii & Kushiro
		50.72   9.16  15.21  16.25  7.06   0.01  1.47   0.39   0.35   0.01   % 4  RE46-5 Yang et al 1996
		45.25   8.89  12.22  24.68  6.45   0.03  1.39   0.67   0.11   0.02   % 5  MIX1G Hirschmann et al 2003
		53.03   9.41  13.17  12.46  7.84   0.15  3.04   1.10   0.39   0.02   % 6  GS104-2-1 Tormey et al 1987
		53.21   9.41  12.21  12.21  8.65   0.09  2.90   1.21   0.69   0.02   % 7  N-MORB Gale et al 2013
		52.12   9.20  13.43  15.66  6.97   0.11  2.08   0.64   0.33   0.01   % 8  ARP 74 Fujii & Bougault 1983
		81.93   8.65   1.16   0.16  0.57   3.19  4.21   0.12   0.04   0.01   % 9  Granite Stern et al 1975
		80.64   9.68   1.75   0     0      3.96  3.96   0      0      0      % 10 Granite R1 Whitney 1975
		60.00  25.00  25.00   0     0       0     0     0      0      0      % 11 anorthite + quartz
		% ---------------------------------------------------------------------
 	*/

	if (test == 0){
		/* SiO2 Al2O3 CaO MgO FeO K2O Na2O TiO2 O Cr2O3 H2O */
		/* Bulk rock composition of Peridotite from Holland et al., 2018, given by E. Green */
		bulk_rock[0]  = 38.494 ;	/** SiO2 */
		bulk_rock[1]  = 1.776;		/** Al2O2 */
		bulk_rock[2]  = 2.824;		/** CaO  */
		bulk_rock[3]  = 50.566;		/** MgO */
		bulk_rock[4]  = 5.886;		/** FeO */
		bulk_rock[5]  = 0.01;		/** K2O	 */
		bulk_rock[6]  = 0.250;		/** Na2O */
		bulk_rock[7]  = 0.10;		/** TiO2 */
		bulk_rock[8]  = 0.096;		/** O */
		bulk_rock[9]  = 0.109;		/** Cr2O3 */
		bulk_rock[10] =	0.0;	
	}
	
	else if (test == 1){
		/* SiO2 Al2O3 CaO MgO FeO K2O Na2O TiO2 O Cr2O3 H2O */
		/* Bulk rock composition of RE46 - Icelandic basalt -Yang et al., 1996, given by E. Green */
		/*   50.72   9.16  15.21  16.25  7.06   0.01  1.47   0.39   0.35   0.01  */
		bulk_rock[0] = 50.72;	
		bulk_rock[1] = 9.16;	
		bulk_rock[2] = 15.21;	
		bulk_rock[3] = 16.25;	
		bulk_rock[4] = 7.06;	
		bulk_rock[5] = 0.01;	
		bulk_rock[6]  = 1.47;
		bulk_rock[7]  = 0.39;
		bulk_rock[8]  = 0.35;
		bulk_rock[9]  = 0.01;
		bulk_rock[10] =	0.0;
	}
	else if (test == 2){
		/* SiO2 Al2O3 CaO MgO FeO K2O Na2O TiO2 O Cr2O3 H2O */
		/* N_MORB, Gale et al., 2013, given by E. Green */
		bulk_rock[0] = 53.21;	
		bulk_rock[1] = 9.41;	
		bulk_rock[2] = 12.21;	
		bulk_rock[3] = 12.21;	
		bulk_rock[4] = 8.65;	
		bulk_rock[5] = 0.09;	
		bulk_rock[6]  = 2.90;
		bulk_rock[7]  = 1.21;
		bulk_rock[8]  = 0.69;
		bulk_rock[9]  = 0.02;
		bulk_rock[10] =	0.0;
	}
	else if (test == 3){
		/* SiO2 Al2O3 CaO MgO FeO K2O Na2O TiO2 O Cr2O3 H2O */
		/* MIX1-G, Hirschmann et al., 2003, given by E. Green */
		bulk_rock[0] = 45.25;	
		bulk_rock[1] = 8.89;	
		bulk_rock[2] = 12.22;	
		bulk_rock[3] = 24.68;	
		bulk_rock[4] = 6.45;	
		bulk_rock[5] = 0.03;	
		bulk_rock[6]  = 1.39;
		bulk_rock[7]  = 0.67;
		bulk_rock[8]  = 0.11;
		bulk_rock[9]  = 0.02;
		bulk_rock[10] =	0.0;
	}
	else if (test == 4){
		/* SiO2 Al2O3 CaO MgO FeO K2O Na2O TiO2 O Cr2O3 H2O */
		/* High Al basalt Baker 1983 */
		bulk_rock[0] = 54.40;	
		bulk_rock[1] = 12.96;	
		bulk_rock[2] = 11.31;	
		bulk_rock[3] = 7.68;	
		bulk_rock[4] = 8.63;	
		bulk_rock[5] = 0.54;	
		bulk_rock[6]  = 3.93;
		bulk_rock[7]  = 0.79;
		bulk_rock[8]  = 0.41;
		bulk_rock[9]  = 0.01;
		bulk_rock[10] =	0.0;
	}	
	else if (test == 5){
		/* SiO2 Al2O3 CaO MgO FeO K2O Na2O TiO2 O Cr2O3 H2O */
		/* Tonalite 101 */
		bulk_rock[0] = 66.01;	
		bulk_rock[1] = 11.98;	
		bulk_rock[2] = 7.06;	
		bulk_rock[3] = 4.16;	
		bulk_rock[4] = 5.30;	
		bulk_rock[5] = 1.57;	
		bulk_rock[6]  = 4.12;
		bulk_rock[7]  = 0.66;
		bulk_rock[8]  = 0.97;
		bulk_rock[9]  = 0.01;
		bulk_rock[10] =	50.0;
	}		
	else if (test == 6){
		/* SiO2 Al2O3 CaO MgO FeO K2O Na2O TiO2 O Cr2O3 H2O */
		/* Bulk rock composition of test 8 */
		bulk_rock[0] = 50.0810;	
		bulk_rock[1] = 8.6901;	
		bulk_rock[2] = 11.6698;	
		bulk_rock[3] = 12.1438;	
		bulk_rock[4] = 7.7832;	
		bulk_rock[5] = 0.2150;
		bulk_rock[6]  = 2.4978;
		bulk_rock[7]  = 1.0059;
		bulk_rock[8]  = 0.4670;
		bulk_rock[9]  = 0.0100;
		bulk_rock[10] =	5.4364;
	}
	else if (test == 7){
		/* SiO2 Al2O3 CaO MgO FeO K2O Na2O TiO2 O Cr2O3 H2O */
		/* Kl3 */
		bulk_rock[0] = 63.242;	
		bulk_rock[1] = 7.267;	
		bulk_rock[2] = 1.815;	
		bulk_rock[3] = 0.751;	
		bulk_rock[4] = 1.311;	
		bulk_rock[5] = 2.774;
		bulk_rock[6]  = 2.622;
		bulk_rock[7]  = 0.144;
		bulk_rock[8]  = 0.5;
		bulk_rock[9]  = 0.01;
		bulk_rock[10] =	19.986;
	}
	else if (test == 8){
		/* SiO2 Al2O3 CaO MgO FeO K2O Na2O TiO2 O Cr2O3 H2O */
		/* Kl3 */
		bulk_rock[0] = 54.65;	
		bulk_rock[1] = 9.04;	
		bulk_rock[2] = 10.69;	
		bulk_rock[3] = 13.9;	
		bulk_rock[4] = 7.63;	
		bulk_rock[5] = 0.51;
		bulk_rock[6]  = 2.62;
		bulk_rock[7]  = 0.71;
		bulk_rock[8]  = 0.4;
		bulk_rock[9]  = 0.1;
		bulk_rock[10] =	0.0;
	}
	else{
		printf("Unknown test %i - please specify a different test! \n", test);
	 	exit(EXIT_FAILURE);
	}
}


/**
  reset global variable for parallel calculations 
*/
global_variable reset_gv(					global_variable 	 gv,
											bulk_info 	 z_b,
											PP_ref 				*PP_ref_db,
											SS_ref 				*SS_ref_db
){
	gv.solver 			  = 0;

	int i,j,k;
	for (k = 0; k < gv.n_flags; k++){
		for (i = 0; i < gv.len_pp; i++){
			gv.pp_flags[i][k]   = 0;
		}
		for (int i = 0; i < gv.len_ss; i++){ 
			SS_ref_db[i].ss_flags[k]   = 0;
		}
	}
	
	/* reset pure phases fractions and xi */
	for (int i = 0; i < gv.len_pp; i++){		
		gv.pp_n[i] 		  = 0.0;
		gv.pp_n_0[i]	  = 0.0;
		gv.delta_pp_n[i]  = 0.0;
		gv.pp_xi[i] 	  = 0.0;
		gv.delta_pp_xi[i] = 0.0;
	}
	
	/* reset pure phases */
	char liq_tail[] = "L";
	for (int i = 0; i < gv.len_pp; i++){
		if ( EndsWithTail(gv.PP_list[i], liq_tail) == 1 ){
			if (z_b.T < 773.0){
				gv.pp_flags[i][0] = 0;
				gv.pp_flags[i][1] = 0;
				gv.pp_flags[i][2] = 0;
				gv.pp_flags[i][3] = 1;
			}
			else{
				gv.pp_flags[i][0] = 1;
				gv.pp_flags[i][1] = 0;
				gv.pp_flags[i][2] = 0;
				gv.pp_flags[i][3] = 0;
			}
		}
		else{
			gv.pp_flags[i][0] = 1;
			gv.pp_flags[i][1] = 0;
			gv.pp_flags[i][2] = 1;
			gv.pp_flags[i][3] = 0;
		}
	}

	gv.LP_PGE_switch	  = 0;
	gv.melt_fraction	  = 0.;
	gv.melt_density       = 0.;
	gv.melt_bulkModulus   = 0.;

	gv.solid_density      = 0.;
	gv.solid_bulkModulus  = 0.;
	gv.solid_shearModulus = 0.;
	gv.solid_Vp 		  = 0.;
	gv.solid_Vs 		  = 0.;

	gv.system_density     = 0.;
	gv.system_bulkModulus = 0.;
	gv.system_shearModulus= 0.;
	gv.system_Vp 		  = 0.;
	gv.system_Vs 		  = 0.;
	gv.V_cor[0]			  = 0.;
	gv.V_cor[1]			  = 0.;
	gv.check_PC1		  = 0;
	gv.check_PC2		  = 0;
	gv.maxeval		      = gv.maxeval_mode_1;
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
		gv.ite_time[i] 			= 0.0;
	}

	/* reset norm and residuals */
    for (i = 0; i < gv.len_ox; i++){	
        gv.mass_residual[i] = 0.0;
        gv.gam_tot[i]     	= 0.0;
        gv.gam_tot_0[i]   	= 0.0;
        gv.delta_gam_tot[i] = 0.0;
		gv.mass_residual[i] = 0.0;	
    }

    for (i = 0; i < gv.len_ss; i++){	
        gv.n_solvi[i] = 0;
		for (k = 0; k < gv.max_n_cp; k++){	
			gv.id_solvi[i][k] = 0;
		} 
    }

	for (i = 0; i < (gv.len_ox); i++){ gv.b[i] = 0.0;
		for (j = 0; j < gv.len_ox; j++){
			gv.A[i][j] = 0.0;
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

	sp[0].frac_S_wt						= 0.0;
	sp[0].frac_M_wt						= 0.0;
	sp[0].frac_F_wt						= 0.0;


	/* reset system */
	for (int i = 0; i < gv.len_ox; i++){
		strcpy(sp[0].ph[i],"");	
		sp[0].bulk[i] 					= 0.0;
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
	}

	/* reset phases */
	for (int n = 0; n < gv.len_ox; n++){
		for (int i = 0; i < gv.len_ox; i++){
			sp[0].PP[n].Comp[i] 			= 0.0;
			sp[0].SS[n].Comp[i] 			= 0.0;
			sp[0].SS[n].compVariables[i] 	= 0.0;
		}
		for (int i = 0; i < gv.len_ox+1; i++){
			strcpy(sp[0].SS[n].emNames[i],"");	
			
			sp[0].SS[n].emFrac[i] 			= 0.0;
			sp[0].SS[n].emChemPot[i] 		= 0.0;

			for (int j = 0; j < gv.len_ox; j++){
				sp[0].SS[n].emComp[i][j]	= 0.0;
			}
		}
	}

}

/**
  reset bulk rock composition informations (needed if the bulk-rock is not constant during parallel computation)
*/
bulk_info reset_z_b_bulk(			global_variable 	 gv,
									double 				*bulk,
									bulk_info 	 		 z_b
){
	int i, j, k;

	int sum = 0;
	for (i = 0; i < nEl; i++) {
		z_b.bulk_rock[i] = bulk[i];
		if (bulk[i] > 0.0){
			sum += 1;
		}
	}

	/** calculate fbc to be used for normalization factor of liq */
	z_b.fbc			= 0.0; 
	for (i = 0; i < nEl; i++){
		z_b.fbc += z_b.bulk_rock[i]*z_b.apo[i];
	}
	
	z_b.nzEl_val = sum;					/** store number of non zero values */
	z_b.zEl_val  = nEl - sum;			/** store number of zero values */
	
	z_b.nzEl_array  = malloc (z_b.nzEl_val * sizeof (int) ); 
	if (z_b.zEl_val > 0){
		z_b.zEl_array   = malloc (z_b.zEl_val * sizeof (int) ); 
		j = 0; k = 0;
		for (i = 0; i < nEl; i++){
			if (bulk[i] == 0.){
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
		for (i = 0; i < nEl; i++){
			z_b.nzEl_array[i] = i;
		}
	}

	for ( i = 0; i < z_b.nzEl_val; i++){
		z_b.bulk_rock_cat[i] = z_b.bulk_rock[z_b.nzEl_array[i]];
	}
	for ( i = z_b.nzEl_val; i < nEl; i++){
		z_b.bulk_rock_cat[i] = 0.0;
	}
	
	return z_b;
};


/**
  reset considered phases entries
*/
void reset_cp(						global_variable 	 gv,
									bulk_info 	 z_b,
									csd_phase_set  		*cp
){
	
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
		
		for (int ii = 0; ii < gv.n_flags; ii++){
			cp[i].ss_flags[ii] 	= 0;
		}

		cp[i].ss_n        		= 0.0;				/* get initial phase fraction */
		cp[i].ss_n_0      		= 0.0;				/* get initial phase fraction */
		cp[i].delta_ss_n    	= 0.0;				/* get initial phase fraction */
		
		for (int ii = 0; ii < gv.len_ox + 1; ii++){
			cp[i].p_em[ii]      = 0.0;
			cp[i].xi_em[ii]     = 0.0;
			cp[i].dguess[ii]    = 0.0;
			cp[i].xeos[ii]      = 0.0;
			cp[i].xeos_0[ii]    = 0.0;
			cp[i].delta_mu[ii]  = 0.0;
			cp[i].dfx[ii]       = 0.0;
			cp[i].mu[ii]        = 0.0;
			cp[i].gbase[ii]     = 0.0;
			cp[i].mu0[ii]       = 0.0;
			cp[i].ss_comp[ii]   = 0.0;
		}
		 
		for (int ii = 0; ii < (gv.len_ox + 1)*2; ii++){
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
									bulk_info 	 z_b,
									SS_ref 				*SS_ref_db
){
	/* reset solution phases */
	for (int iss = 0; iss < gv.len_ss; iss++){

		for (int j = 0; j < gv.n_flags; j++){	
			SS_ref_db[iss].ss_flags[j]   = 0;
		}

		SS_ref_db[iss].min_mode	= 1;
		SS_ref_db[iss].tot_pc 	= 0;
		SS_ref_db[iss].id_pc  	= 0;
		for (int j = 0; j < gv.len_ox; j++){
			SS_ref_db[iss].solvus_id[j] = -1;	
		}

		/* reset levelling pseudocompounds */
		for (int i = 0; i < (SS_ref_db[iss].n_pc); i++){
			SS_ref_db[iss].n_swap[i] = 0;
			SS_ref_db[iss].info[i]   = 0;
			SS_ref_db[iss].G_pc[i]   = 0.0;
			SS_ref_db[iss].DF_pc[i]  = 0.0;
			for (int j = 0; j < gv.len_ox; j++){
				SS_ref_db[iss].comp_pc[i][j]  = 0.0;
			}
			for (int j = 0; j < SS_ref_db[iss].n_em; j++){
				SS_ref_db[iss].p_pc[i][j]  = 0.0;	
				SS_ref_db[iss].mu_pc[i][j] = 0.0;	
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
			SS_ref_db[iss].n_swap_Ppc[i] = 0;
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
			SS_ref_db[iss].factor_Ppc[i] = 0.0;
		}

		/* reset solution phase model parameters */
		for (int j = 0; j < SS_ref_db[iss].n_em; j++){
			SS_ref_db[iss].xi_em[j]      = 0.0;
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
			SS_ref_db[iss].xeos_sf_ok[k] = 0.0;
		}

		for (int j = 0; j < SS_ref_db[iss].n_em; j++){
			SS_ref_db[iss].p[j]     = 0.0;
			SS_ref_db[iss].ape[j]   = 0.0;
		}
		SS_ref_db[iss].forced_stop = 0; 
		SS_ref_db[iss].min_mode    = 1;
		SS_ref_db[iss].nlopt_verb  = 0; // no output by default
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
	d->dG_B_tol	   = -1e-6;
	d->min_F_tol   =  1e6;
	
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
    for (int i = 0; i < gv.len_ox; i++){
		d->gamma_tot[i] 	= 0.0;
		d->gamma_delta[i] 	= 0.0;
	}
	int k;
    for (int i = 0; i < gv.len_ox; i++){
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
	d->n_local_min =  0;
	d->n_filter    =  0;
	d->ph2swp      = -1;
	d->n_swp       =  0;
	d->swp         =  0;
	d->n_Ox        =  z_b.nzEl_val;

	/* initialize arrays */
    for (int i = 0; i < gv.len_ox; i++){
		d->gamma_tot[i] 	= 0.0;
		d->gamma_delta[i] 	= 0.0;
	}
	int k;
    for (int i = 0; i < gv.len_ox; i++){
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

#endif
