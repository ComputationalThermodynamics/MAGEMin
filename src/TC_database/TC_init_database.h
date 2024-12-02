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
	typedef struct metapelite_datasets {
		int 	ds_version;
		int 	n_ox;
		int 	n_pp;
		int 	n_ss;
		char    ox[11][20];
		char    PP[24][20];
		char    SS[18][20];

		int 	verifyPC[18];
		int 	n_SS_PC[18];
		double 	SS_PC_stp[18];

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
		char    PP[24][20];

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
		Igneous database informations 
	**/
	typedef struct igneous_datasets {
		int 	ds_version;
		int 	n_ox;
		int 	n_pp;
		int 	n_ss;
		char    ox[11][20];
		char    PP[24][20];
		char    SS[15][20];

		int 	verifyPC[15];
		int 	n_SS_PC[15];
		double 	SS_PC_stp[15];

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
		char    PP[23][20];
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
		char    PP[21][20];
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
		char    PP[21][20];
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
		char    PP[26][20];
		char    SS[23][20];

		int 	verifyPC[23];
		int 	n_SS_PC[23];
		double 	SS_PC_stp[23];

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