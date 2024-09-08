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
#ifndef __SB_init_db_H_
#define __SB_init_db_H_

    #include "../MAGEMin.h"

	/** 
		Store oxide informations 
	**/
	typedef struct oxide_datas {
		int 	n_ox;
		char    oxName[15][20];
		double  oxMass[15];
		double  atPerOx[15];
		double  ElEntropy[15]; //standard molar entropy
		double  OPerOx[15];

	} oxide_data;

	/** 
		Igneous database informations 
	**/
	typedef struct stx11_datasets {
		int 	ds_version;
		int 	n_ox;
		int 	n_pp;
		int 	n_ss;
		char    ox[6][20];
		char    PP[10][20];
		char    SS[14][20];

		int 	verifyPC[14];
		int 	n_SS_PC[14];
		double 	SS_PC_stp[14];

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

	} stx11_dataset;

    global_variable global_variable_SB_init( 	global_variable  	 gv,
                                            	bulk_info 			*z_b 	);

    global_variable get_bulk_stx11( global_variable gv);


#endif