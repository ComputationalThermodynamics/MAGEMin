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
#ifndef __TC_init_db_H_
#define __TC_init_db_H_

    #include "../MAGEMin.h"



	/** 
		Metapelite database informations
	**/
	#define n_ox_mpf 12 		// as of now this is n_oxides = 11 + electrochemical potential
	#define n_ss_mpf 18			
	#define n_pp_mpf 26

	typedef struct metapelite_aq_datasets {
		int 	ds_version;
		int 	n_ox;
		int 	n_pp;
		int 	n_ss;
		char    ox[n_ox_mpf][20];
		char    PP[n_pp_mpf][20];
		char    SS[n_ss_mpf][20];

		int 	verifyPC[n_ss_mpf];
		int 	n_SS_PC[n_ss_mpf];
		double 	SS_PC_stp[n_ss_mpf];

		double 	PC_df_add;					/** min value of df under which the PC is added 									*/
		double  solver_switch_T;
		double  min_melt_T;

		double  inner_PGE_ite;				/** number of inner PGE iterations, this has to be made mass or dG dependent 		*/
		double  max_n_phase;				/** maximum mol% phase change during one PGE iteration in wt% 						*/
		double  max_g_phase;				/** maximum delta_G of reference change during PGE 									*/
		double 	max_fac;					/** maximum update factor during PGE under-relax < 0.0, over-relax > 0.0 	 		*/

		double  merge_value;				/** max norm distance between two instances of a solution phase						*/	
		double 	re_in_n;					/** fraction of phase when being reintroduce.  										*/

		double  obj_tol;

	} metapelite_aq_dataset;



	/** 
		Metapelite database informations
	**/
	#define n_ss_mp 18
	#define n_pp_mp 26
	#define n_ox_mp 11
	typedef struct metapelite_datasets {
		int 	ds_version;
		int 	n_ox;
		int 	n_pp;
		int 	n_ss;
		char    ox[n_ox_mp][20];
		char    PP[n_pp_mp][20];
		char    SS[n_ss_mp][20];

		int 	verifyPC[n_ss_mp];
		int 	n_SS_PC[n_ss_mp];
		double 	SS_PC_stp[n_ss_mp];

		double 	PC_df_add;					/** min value of df under which the PC is added 									*/
		double  solver_switch_T;
		double  min_melt_T;

		double  inner_PGE_ite;				/** number of inner PGE iterations, this has to be made mass or dG dependent 		*/
		double  max_n_phase;				/** maximum mol% phase change during one PGE iteration in wt% 						*/
		double  max_g_phase;				/** maximum delta_G of reference change during PGE 									*/
		double 	max_fac;					/** maximum update factor during PGE under-relax < 0.0, over-relax > 0.0 	 		*/

		double  merge_value;				/** max norm distance between two instances of a solution phase						*/	
		double 	re_in_n;					/** fraction of phase when being reintroduce.  										*/

		double  obj_tol;

	} metapelite_dataset;

	/** 
		Metabasite database informations
	**/
		typedef struct metabasite_datasets {
		int 	ds_version;
		int 	n_ox;
		int 	n_pp;
		int 	n_ss;
		char    ox[10][20];
		char    PP[27][20];

		char    SS[17][20];
		int 	verifyPC[17];
		int 	n_SS_PC[17];
		double 	SS_PC_stp[17];

		double 	PC_df_add;					/** min value of df under which the PC is added 									*/
		double  solver_switch_T;
		double  min_melt_T;

		double  inner_PGE_ite;				/** number of inner PGE iterations, this has to be made mass or dG dependent 		*/
		double  max_n_phase;				/** maximum mol% phase change during one PGE iteration in wt% 						*/
		double  max_g_phase;				/** maximum delta_G of reference change during PGE 									*/
		double 	max_fac;					/** maximum update factor during PGE under-relax < 0.0, over-relax > 0.0 	 		*/

		double  merge_value;				/** max norm distance between two instances of a solution phase						*/	
		double 	re_in_n;					/** fraction of phase when being reintroduce.  										*/

		double  obj_tol;

	} metabasite_dataset;

	/** 
		Metabasite ext database informations
	**/
		typedef struct metabasite_ext_datasets {
		int 	ds_version;
		int 	n_ox;
		int 	n_pp;
		int 	n_ss;
		char    ox[10][20];
		char    PP[276][20];

		char    SS[19][20];
		int 	verifyPC[19];
		int 	n_SS_PC[19];
		double 	SS_PC_stp[19];

		double 	PC_df_add;					/** min value of df under which the PC is added 									*/
		double  solver_switch_T;
		double  min_melt_T;

		double  inner_PGE_ite;				/** number of inner PGE iterations, this has to be made mass or dG dependent 		*/
		double  max_n_phase;				/** maximum mol% phase change during one PGE iteration in wt% 						*/
		double  max_g_phase;				/** maximum delta_G of reference change during PGE 									*/
		double 	max_fac;					/** maximum update factor during PGE under-relax < 0.0, over-relax > 0.0 	 		*/

		double  merge_value;				/** max norm distance between two instances of a solution phase						*/	
		double 	re_in_n;					/** fraction of phase when being reintroduce.  										*/

		double  obj_tol;

	} metabasite_ext_dataset;

	/** 
		Igneous database informations 
	**/
	#define n_ox_ig 11 		// as of now this is n_oxides = 11 + electrochemical potential
	#define n_ss_ig 16			
	#define n_pp_ig 26
	typedef struct igneous_datasets {
		int 	ds_version;
		int 	n_ox;
		int 	n_pp;
		int 	n_ss;
		char    ox[n_ox_ig][20];
		char    PP[n_pp_ig][20];
		int 	act_PP[n_pp_ig];
		char    SS[n_ss_ig][20];

		int 	verifyPC[n_ss_ig];
		int 	n_SS_PC[n_ss_ig];
		double 	SS_PC_stp[n_ss_ig];

		double 	PC_df_add;	
		double  solver_switch_T;
		double  min_melt_T;

		double  inner_PGE_ite;				/** number of inner PGE iterations, this has to be made mass or dG dependent 		*/
		double  max_n_phase;				/** maximum mol% phase change during one PGE iteration in wt% 						*/
		double  max_g_phase;				/** maximum delta_G of reference change during PGE 									*/
		double 	max_fac;					/** maximum update factor during PGE under-relax < 0.0, over-relax > 0.0 	 		*/

		double  merge_value;				/** max norm distance between two instances of a solution phase						*/	
		double 	re_in_n;					/** fraction of phase when being reintroduce.  										*/

		double  obj_tol;

	} igneous_dataset;

	/** 
		Igneous Alkali database wet test
	**/
	typedef struct igneous_igad_datasets {
		int 	ds_version;
		int 	n_ox;
		int 	n_pp;
		int 	n_ss;
		char    ox[10][20];
		char    PP[24][20];
		char    SS[12][20];

		int 	verifyPC[12];
		int 	n_SS_PC[12];
		double 	SS_PC_stp[12];

		double 	PC_df_add;	
		double  solver_switch_T;
		double  min_melt_T;

		double  inner_PGE_ite;				/** number of inner PGE iterations, this has to be made mass or dG dependent 		*/
		double  max_n_phase;				/** maximum mol% phase change during one PGE iteration in wt% 						*/
		double  max_g_phase;				/** maximum delta_G of reference change during PGE 									*/
		double 	max_fac;					/** maximum update factor during PGE under-relax < 0.0, over-relax > 0.0 	 		*/

		double  merge_value;				/** max norm distance between two instances of a solution phase						*/	
		double 	re_in_n;					/** fraction of phase when being reintroduce.  										*/

		double  obj_tol;

	} igneous_igad_dataset;


	/** 
		Evans&Frost,2021 database informations
	**/
	typedef struct ultramafic_datasets {
		int 	ds_version;
		int 	n_ox;
		int 	n_pp;
		int 	n_ss;
		char    ox[7][20];
		char    PP[24][20];
		char    SS[12][20];

		int 	verifyPC[12];
		int 	n_SS_PC[12];
		double 	SS_PC_stp[12];

		double 	PC_df_add;					/** min value of df under which the PC is added 									*/
		double  solver_switch_T;
		double  min_melt_T;

		double  inner_PGE_ite;				/** number of inner PGE iterations, this has to be made mass or dG dependent 		*/
		double  max_n_phase;				/** maximum mol% phase change during one PGE iteration in wt% 						*/
		double  max_g_phase;				/** maximum delta_G of reference change during PGE 									*/
		double 	max_fac;					/** maximum update factor during PGE under-relax < 0.0, over-relax > 0.0 	 		*/

		double  merge_value;				/** max norm distance between two instances of a solution phase						*/	
		double 	re_in_n;					/** fraction of phase when being reintroduce.  										*/

		double  obj_tol;

	} ultramafic_dataset;

	typedef struct ultramafic_ext_datasets {
		int 	ds_version;
		int 	n_ox;
		int 	n_pp;
		int 	n_ss;
		char    ox[9][20];
		char    PP[24][20];
		char    SS[15][20];

		int 	verifyPC[15];
		int 	n_SS_PC[15];
		double 	SS_PC_stp[15];

		double 	PC_df_add;					/** min value of df under which the PC is added 									*/
		double  solver_switch_T;
		double  min_melt_T;

		double  inner_PGE_ite;				/** number of inner PGE iterations, this has to be made mass or dG dependent 		*/
		double  max_n_phase;				/** maximum mol% phase change during one PGE iteration in wt% 						*/
		double  max_g_phase;				/** maximum delta_G of reference change during PGE 									*/
		double 	max_fac;					/** maximum update factor during PGE under-relax < 0.0, over-relax > 0.0 	 		*/

		double  merge_value;				/** max norm distance between two instances of a solution phase						*/	
		double 	re_in_n;					/** fraction of phase when being reintroduce.  										*/

		double  obj_tol;

	} ultramafic_ext_dataset;


	/** 
		Metabasite database informations
	**/
	typedef struct mantle_datasets {
		int 	ds_version;
		int 	n_ox;
		int 	n_pp;
		int 	n_ss;
		char    ox[6][20];
		char    PP[8][20];
		char    SS[14][20];

		int 	verifyPC[14];
		int 	n_SS_PC[14];
		double 	SS_PC_stp[14];

		double 	PC_df_add;					/** min value of df under which the PC is added 									*/
		double  solver_switch_T;
		double  min_melt_T;

		double  inner_PGE_ite;				/** number of inner PGE iterations, this has to be made mass or dG dependent 		*/
		double  max_n_phase;				/** maximum mol% phase change during one PGE iteration in wt% 						*/
		double  max_g_phase;				/** maximum delta_G of reference change during PGE 									*/
		double 	max_fac;					/** maximum update factor during PGE under-relax < 0.0, over-relax > 0.0 	 		*/

		double  merge_value;				/** max norm distance between two instances of a solution phase						*/	
		double 	re_in_n;					/** fraction of phase when being reintroduce.  										*/

		double  obj_tol;

	} mantle_dataset;


	/** 
		Metapelite database informations
	**/
	typedef struct metapelite_datasets_ext {
		int 	ds_version;
		int 	n_ox;
		int 	n_pp;
		int 	n_ss;
		char    ox[13][20];
		char    PP[31][20];
		char    SS[24][20];

		int 	verifyPC[24];
		int 	n_SS_PC[24];
		double 	SS_PC_stp[24];

		double 	PC_df_add;					/** min value of df under which the PC is added 									*/
		double  solver_switch_T;
		double  min_melt_T;

		double  inner_PGE_ite;				/** number of inner PGE iterations, this has to be made mass or dG dependent 		*/
		double  max_n_phase;				/** maximum mol% phase change during one PGE iteration in wt% 						*/
		double  max_g_phase;				/** maximum delta_G of reference change during PGE 									*/
		double 	max_fac;					/** maximum update factor during PGE under-relax < 0.0, over-relax > 0.0 	 		*/

		double  merge_value;				/** max norm distance between two instances of a solution phase						*/	
		double 	re_in_n;					/** fraction of phase when being reintroduce.  										*/

		double  obj_tol;

	} metapelite_dataset_ext;

    global_variable global_variable_TC_init( 	global_variable  	 gv,
                                            	bulk_info 			*z_b 	);

    global_variable get_bulk_metapelite( 		global_variable gv);
    global_variable get_bulk_metabasite( 		global_variable gv);
    global_variable get_bulk_igneous( 			global_variable gv);
	global_variable get_bulk_igneous_igad( 		global_variable gv);
    global_variable get_bulk_ultramafic( 		global_variable gv);
    global_variable get_bulk_ultramafic_ext( 	global_variable gv);
    global_variable get_bulk_mantle( 			global_variable gv);
    global_variable get_bulk_metapelite_ext( 	global_variable gv);
#endif