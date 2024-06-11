#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <complex.h> 

#include "../uthash.h"
#include "../MAGEMin.h"
#include "../all_endmembers.h"
#include "TC_init_database.h"


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
	gv.db 				= malloc (5 	* sizeof(char)			);
	gv.Phase 			= malloc (50 	* sizeof(char)			);
	gv.sys_in 			= malloc (5 	* sizeof(char)			);
	gv.buffer 			= malloc (10 	* sizeof(char)			);

	gv.arg_bulk 		= malloc (gv.maxlen_ox * sizeof(double)	);
	gv.arg_gamma 		= malloc (gv.maxlen_ox * sizeof(double)	);

	for (i = 0; i < gv.maxlen_ox; i++) {
		gv.arg_bulk[i]  = 0.0;
		gv.arg_gamma[i] = 0.0;
	}

	strcpy(gv.outpath,"./output/");				/** define the outpath to save logs and final results file	 						*/
	strcpy(gv.version,"1.4.9 [05/06/2024]");	/** MAGEMin version 																*/

	/* generate parameters        		*/
	strcpy(gv.buffer,"none");	
	gv.max_n_mSS		= 128;					/** maximum number of metastable pseudocompounds 									*/
	gv.max_n_cp 		= 128;					/** number of considered solution phases 											*/	
	gv.max_ss_size_cp   = 24;					/** maximum size for a solution phase saved in the cp structure                     */
	gv.buffer_n 		= 0.0;					/** factor for QFM buffer 															*/
	gv.limitCaOpx       = 0;					/** limit Ca-bearing  orthopyroxene (add-hoc correction) 							*/
	gv.CaOpxLim         = 1.0;					/** limit Ca-bearing  orthopyroxene (add-hoc correction) 							*/
	gv.mbCpx 			= 0;					/** 0: omphacite LT, 1: augite HT*/
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

	/* levelling parameters 			*/
	gv.em2ss_shift		= 1e-8;					/** small value to shift x-eos of pure endmember from bounds after levelling 		*/
	gv.bnd_filter_pc    = 10.0;					/** value of driving force the pseudocompound is considered 						*/
	gv.bnd_filter_pge   = 2.5;					/** value of driving force the pseudocompound is considered 						*/
	gv.max_G_pc         = 2.5;					/** dG under which PC is considered after their generation		 					*/
	gv.eps_sf_pc		= 1e-10;				/** Minimum value of site fraction under which PC is rejected, 
													don't put it too high as it will conflict with bounds of x-eos					*/

	/* PGE LP pseudocompounds parameters */
	gv.launch_PGE 		= 0;
	gv.n_pc 			= 8192;
	gv.n_Ppc			= 8192;
	gv.max_LP_ite 		= 256;
	gv.save_Ppc_val     = 0.0; 					/** During PGE iterations, if the driving force is < save_Ppc_val, then the 
													pseudocompound is added to the Ppc list 										*/

	/* local minimizer options 	*/
	gv.bnd_val          = 1.0e-10;				/** boundary value for x-eos 										 				*/
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
	gv.EM_database  	=  2; 					
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

	strcpy(gv.File,		"none"); 	/** Filename to be read to have multiple P-T-bulk conditions to solve 	*/
	strcpy(gv.sys_in,	"mol"); 	/** system unit 														*/
	strcpy(gv.db,		"ig"); 		/** database 															*/

	return gv;
}