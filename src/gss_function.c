/**
Function to calculate the reference chemical potential of solid-solutions        

Holland et al., 2018 - Melting of peridotites through to granites                 
Igneous dataset to use with tc-ds633.txt                                         
"bi","cpx","cd","ep","fl","g","hb","ilm","liq","mu", "ol", "opx","fsp","spn" 
*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <complex.h> 
#include "MAGEMin.h"
#include "gem_function.h"
#include "gss_function.h"
#include "NLopt_opt_function.h"
#include "simplex_levelling.h"
#include "toolkit.h"
// #include "gfs_function.h"

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
em_data get_em_data(	int 		 EM_database, 
						int          len_ox,
						bulk_info 	 z_b,
                        double       P,
                        double       T,
						char 		*name, 
						char 		*state		){

	em_data data; 
	PP_ref PP_db   		= G_EM_function(EM_database, len_ox, z_b.id, z_b.bulk_rock, z_b.apo, P, T, name, state);
   	data.ElShearMod  	= PP_db.phase_shearModulus;
   	data.gb  			= PP_db.gbase;

	for (int i = 0; i < len_ox; i++){
		data.C[i] = PP_db.Comp[i];
	}
	return data;
}

em_data get_fs_data(	int             len_ox,
						bulk_info 	    z_b,
                        solvent_prop   *wat,
                        double          P,
                        double          T,
						char 		   *name, 
						char 		   *state		){

	em_data data; 
	PP_ref PP_db   		= G_FS_function(    len_ox,
                                            wat,
                                            z_b.id,
                                            z_b.bulk_rock,
                                            z_b.ElEntropy,
                                            z_b.apo,
                                            P,
                                            T,
                                            name,
                                            state     );

   	data.ElShearMod  	= PP_db.phase_shearModulus;
   	data.gb  			= PP_db.gbase;
    data.charge         = PP_db.charge;

	for (int i = 0; i < len_ox; i++){
		data.C[i] = PP_db.Comp[i];
	}
	return data;
}


/**************************************************************************************/
/* Import text file extracted from Miron et al. (2017) database:                      */
/* aq17-thermofun.json from db.thermohub.org, v. 17.07.2021 16:49:29                  */
/*--------------------------------------------------------------------------          */
/* Miron, G. D., Wagner, T., Kulik, D. A., & Lothenbach, B. (2017). An                */
/* internally consistent thermodynamic dataset for aqueous species in the             */
/* system Ca-Mg-Na-K-Al-Si-OHC-Cl to 800 C and 5 kbar. American Journal               */
/* of Science, 317(7), 755-806.                                                       */
/* DOI: https://doi.org/10.2475/07.2017.01                                            */
/**************************************************************************************/


SS_ref G_SS_aq17_function(SS_ref SS_ref_db, int EM_database, int len_ox, bulk_info z_b, double eps){
    
    int i, j;
    int n_em    = SS_ref_db.n_em;
    double eps2 = 1e-15;
    // char   *EM_tmp[] 		= {"H2O","Al(OH)2+", "Al(OH)3@", "Al(OH)4-", "Al+3", "AlH3SiO4+2", "AlOH+2", "CO3-2", "CO@", "Ca+2", "CaCO3@", "CaCl+", "CaCl2@", "CaHCO3+", "CaHSiO3+", "CaOH+", "CaSiO3@", "Cl-", "HCO3-", "HCl@", "HSiO3-", "K+", "KAlO2@", "KCO3-", "KCl@", "KHCO3@", "KOH@", "Mg+2", "MgCO3@", "MgCl+", "MgCl2@", "MgHCO3+", "MgHSiO3+", "MgOH+", "MgSiO3@", "Na+", "NaAl(OH)4@", "NaCO3-", "NaCl@", "NaHCO3@", "NaHSiO3@", "NaOH@", "OH-", "SiO2@"};
    char   *EM_tmp[] 		= {"H2O","Al(OH)2+", "Al(OH)3@", "Al(OH)4-", "Al+3", "AlH3SiO4+2", "AlOH+2", "Ca+2", "CaHSiO3+", "CaOH+", "CaSiO3@", "HSiO3-", "K+", "KAlO2@", "KOH@", "Mg+2", "MgHSiO3+", "MgOH+", "MgSiO3@", "Na+", "NaAl(OH)4@", "NaHSiO3@", "NaOH@", "OH-", "SiO2@"};
    for (int i = 0; i < n_em; i++){
        strcpy(SS_ref_db.EM_list[i],EM_tmp[i]);
    };

    solvent_prop wat;

    rho_wat_calc(      &wat,
                        SS_ref_db.P*1000.0,
                        SS_ref_db.T,
                        "WP"                           );

    char   solventOpt[] =  {"JN91"};

    if (strcmp( solventOpt, "JN91") == 0){
        propSolvent_JN91_calc(     &wat,
                                    SS_ref_db.T         );
    }
    else if (strcmp( solventOpt, "FE97") == 0){
        propSolvent_FE97_calc(     &wat,
                                    SS_ref_db.P*1000.0,
                                    SS_ref_db.T         );                                
    }
    else if (strcmp( solventOpt, "SV14") == 0){
        propSolvent_SV14_calc(     &wat,
                                    SS_ref_db.P*1000.0,
                                    SS_ref_db.T         );
    }
    SS_ref_db.densityW      = wat.density;
    SS_ref_db.g             = wat.g;
    SS_ref_db.Z             = wat.Z;
    SS_ref_db.epsilon       = wat.epsilon;

    em_data species;

    species 	= get_em_data(	    EM_database, 
                                    len_ox,
                                    z_b,
                                    SS_ref_db.P,
                                    SS_ref_db.T,
                                    "H2O", 
                                    "equilibrium"	);

    SS_ref_db.gbase[0] 		    = species.gb;
    SS_ref_db.z_em[0]           = 1.0;
    SS_ref_db.ElShearMod[0] 	= 0.0;
    SS_ref_db.bounds_ref[0][0]  = 0.0+eps2;  SS_ref_db.bounds_ref[0][1] = 1.0-eps2;
    SS_ref_db.mat_phi[0]        = 0.0; //serve to store charge here

    for (j = 0; j < len_ox; j++){
        SS_ref_db.Comp[0][j] 	= species.C[j];
    }

    /* start at 1, as 0 is water */
    for (i = 1; i < n_em; i++){

        species 	= get_fs_data(	len_ox,
                                    z_b,
                                    &wat,
                                    SS_ref_db.P,
                                    SS_ref_db.T,
                                    SS_ref_db.EM_list[i], 
                                    "equilibrium"	);

        SS_ref_db.gbase[i] 		    = species.gb;
        SS_ref_db.z_em[i]           = 1.0;
        SS_ref_db.ElShearMod[i] 	= 0.0;
        SS_ref_db.bounds_ref[i][0]  = 0.0+eps2;  SS_ref_db.bounds_ref[i][1] = 1.0-eps2;
        SS_ref_db.mat_phi[i]        = species.charge; //serves to store charge here

        for (j = 0; j < len_ox; j++){
            SS_ref_db.Comp[i][j] 	= species.C[j];
        }
    };


    /* copy molar elemental entropy */
    for(int i = 0; i < len_ox; i++){
        SS_ref_db.ElEntropy[i] = z_b.ElEntropy[i];
    }
    SS_ref_db.len_ox           = len_ox;

	// if (z_b.bulk_rock[10] == 0.){ 					
	// 	SS_ref_db.z_em[14]          = 0.0;
	// 	SS_ref_db.bounds_ref[9][0] = eps; 
	// 	SS_ref_db.bounds_ref[9][1] = eps;	
	// }
	// if (z_b.bulk_rock[9] == 0.){ 					
	// 	SS_ref_db.z_em[7]          = 0.0;
    //     SS_ref_db.d_em[7]          = 1.0;
	// 	SS_ref_db.bounds_ref[6][0] = 0.0; 
	// 	SS_ref_db.bounds_ref[6][1] = 0.0;	
	// }
	// if (z_b.bulk_rock[8] == 0.){ 					
	// 	SS_ref_db.z_em[6]          = 0.0;
    //     SS_ref_db.d_em[6]          = 1.0;
	// 	SS_ref_db.bounds_ref[5][0] = 0.0; 
	// 	SS_ref_db.bounds_ref[5][1] = 0.0;	
	// }


    return SS_ref_db;
}

/**************************************************************************************/
/**************************************************************************************/
/**********************METABASITE DATABASE (Gree et al., 2016)*************************/
/**************************************************************************************/
/**************************************************************************************/

/**
   retrieve reference thermodynamic data for mb_liq
*/
SS_ref G_SS_mb_liq_function(SS_ref SS_ref_db, int EM_database, int len_ox, bulk_info z_b, double eps){
    
    int i, j;
    int n_em = SS_ref_db.n_em;
    
    char   *EM_tmp[] 		= {"q4L","abL","kspL","wo1L","sl1L","fa2L","fo2L","watL","anoL"};
    for (int i = 0; i < SS_ref_db.n_em; i++){
        strcpy(SS_ref_db.EM_list[i],EM_tmp[i]);
    };
    int n_xeos = SS_ref_db.n_xeos;
    char   *CV_tmp[] 		= {"q","fsp","na","wo","sil","ol","x","yan"};
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
        strcpy(SS_ref_db.CV_list[i],CV_tmp[i]);
    };
    SS_ref_db.W[0] = 12.0 - 0.4*SS_ref_db.P;
    SS_ref_db.W[1] = -0.5*SS_ref_db.P - 2.0;
    SS_ref_db.W[2] = -5.00;
    SS_ref_db.W[3] = 0.0;
    SS_ref_db.W[4] = 0.0;
    SS_ref_db.W[5] = SS_ref_db.P + 42.0;
    SS_ref_db.W[6] = 18.1 - 0.68*SS_ref_db.P;
    SS_ref_db.W[7] = -0.1*SS_ref_db.P - 29.5;
    SS_ref_db.W[8] = 3.0*SS_ref_db.P - 6.0;
    SS_ref_db.W[9] = -12.0;
    SS_ref_db.W[10] = 10.0;
    SS_ref_db.W[11] = 0.8*SS_ref_db.P - 30.0;
    SS_ref_db.W[12] = 0.3*SS_ref_db.P - 47.3;
    SS_ref_db.W[13] = -0.17*SS_ref_db.P - 4.4;
    SS_ref_db.W[14] = 0.4*SS_ref_db.P + 8.6;
    SS_ref_db.W[15] = -13.0;
    SS_ref_db.W[16] = 0.0;
    SS_ref_db.W[17] = -11.3;
    SS_ref_db.W[18] = 6.80;
    SS_ref_db.W[19] = 10.4 - 0.39*SS_ref_db.P;
    SS_ref_db.W[20] = -0.25*SS_ref_db.P - 16.0;
    SS_ref_db.W[21] = -1.60;
    SS_ref_db.W[22] = 6.50;
    SS_ref_db.W[23] = 4.00;
    SS_ref_db.W[24] = 21.0;
    SS_ref_db.W[25] = 3.50;
    SS_ref_db.W[26] = 12.0;
    SS_ref_db.W[27] = 12.0;
    SS_ref_db.W[28] = 11.0 - 0.5*SS_ref_db.P;
    SS_ref_db.W[29] = 6.40;
    SS_ref_db.W[30] = 18.0;
    SS_ref_db.W[31] = 29.0;
    SS_ref_db.W[32] = -0.95*SS_ref_db.P - 43.5;
    SS_ref_db.W[33] = 29.0 - 0.5*SS_ref_db.P;
    SS_ref_db.W[34] = -0.6*SS_ref_db.P - 26.0;
    SS_ref_db.W[35] = 9.75 - 0.5*SS_ref_db.P;
    
    
    em_data qL_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"qL", 
    										"equilibrium"	);
    
    em_data abL_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"abL", 
    										"equilibrium"	);
    
    em_data kspL_eq 	= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"kspL", 
    										"equilibrium"	);
    
    em_data woL_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"woL", 
    										"equilibrium"	);
    
    em_data silL_eq 	= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"silL", 
    										"equilibrium"	);
    
    em_data faL_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"faL", 
    										"equilibrium"	);
    
    em_data foL_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"foL", 
    										"equilibrium"	);
    
    em_data watL_eq 	= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"h2oL", 
    										"equilibrium"	);
    
    SS_ref_db.gbase[0] 		= 4.0*qL_eq.gb;
    SS_ref_db.gbase[1] 		= abL_eq.gb;
    SS_ref_db.gbase[2] 		= kspL_eq.gb;
    SS_ref_db.gbase[3] 		= woL_eq.gb + 1.3;
    SS_ref_db.gbase[4] 		= silL_eq.gb - 7.8;
    SS_ref_db.gbase[5] 		= -1.4*z_b.P + 2.0*faL_eq.gb - 8.2;
    SS_ref_db.gbase[6] 		= 2.0*foL_eq.gb - 4.0;
    SS_ref_db.gbase[7] 		= watL_eq.gb;
    SS_ref_db.gbase[8] 		= -0.25*z_b.P + silL_eq.gb + woL_eq.gb - 46.5;
    
    SS_ref_db.ElShearMod[0] 	= qL_eq.ElShearMod;
    SS_ref_db.ElShearMod[1] 	= abL_eq.ElShearMod;
    SS_ref_db.ElShearMod[2] 	= kspL_eq.ElShearMod;
    SS_ref_db.ElShearMod[3] 	= woL_eq.ElShearMod;
    SS_ref_db.ElShearMod[4] 	= silL_eq.ElShearMod;
    SS_ref_db.ElShearMod[5] 	= 2.0*faL_eq.ElShearMod;
    SS_ref_db.ElShearMod[6] 	= 2.0*foL_eq.ElShearMod;
    SS_ref_db.ElShearMod[7] 	= watL_eq.ElShearMod;
    SS_ref_db.ElShearMod[8] 	= silL_eq.ElShearMod + woL_eq.ElShearMod;
    
    for (i = 0; i < len_ox; i++){
        SS_ref_db.Comp[0][i] 	= 4.0*qL_eq.C[i];
        SS_ref_db.Comp[1][i] 	= abL_eq.C[i];
        SS_ref_db.Comp[2][i] 	= kspL_eq.C[i];
        SS_ref_db.Comp[3][i] 	= woL_eq.C[i];
        SS_ref_db.Comp[4][i] 	= silL_eq.C[i];
        SS_ref_db.Comp[5][i] 	= 2.0*faL_eq.C[i];
        SS_ref_db.Comp[6][i] 	= 2.0*foL_eq.C[i];
        SS_ref_db.Comp[7][i] 	= watL_eq.C[i];
        SS_ref_db.Comp[8][i] 	= silL_eq.C[i] + woL_eq.C[i];
    }
    
    for (i = 0; i < n_em; i++){
        SS_ref_db.z_em[i] = 1.0;
    };
    
    SS_ref_db.bounds_ref[0][0] = 0.0+eps;  SS_ref_db.bounds_ref[0][1] = 1.0-eps;
    SS_ref_db.bounds_ref[1][0] = 0.0+eps;  SS_ref_db.bounds_ref[1][1] = 1.0-eps;
    SS_ref_db.bounds_ref[2][0] = 0.0+eps;  SS_ref_db.bounds_ref[2][1] = 1.0-eps;
    SS_ref_db.bounds_ref[3][0] = 0.0+eps;  SS_ref_db.bounds_ref[3][1] = 1.0-eps;
    SS_ref_db.bounds_ref[4][0] = 0.0+eps;  SS_ref_db.bounds_ref[4][1] = 1.0-eps;
    SS_ref_db.bounds_ref[5][0] = 0.0+eps;  SS_ref_db.bounds_ref[5][1] = 1.0-eps;
    SS_ref_db.bounds_ref[6][0] = 0.0+eps;  SS_ref_db.bounds_ref[6][1] = 1.0-eps;
    SS_ref_db.bounds_ref[7][0] = 0.0+eps;  SS_ref_db.bounds_ref[7][1] = 1.0-eps;

	// if (z_b.bulk_rock[9] == 0.){ 					
	// 	SS_ref_db.z_em[7]          = 0.0;
	// 	SS_ref_db.bounds_ref[9][0] = eps; 
	// 	SS_ref_db.bounds_ref[9][1] = eps;	
	// }
	// if (z_b.bulk_rock[9] == 0.){ 					
	// 	SS_ref_db.z_em[7]          = 0.0;
    //     SS_ref_db.d_em[7]          = 1.0;
	// 	SS_ref_db.bounds_ref[6][0] = 0.0; 
	// 	SS_ref_db.bounds_ref[6][1] = 0.0;	
	// }
	// if (z_b.bulk_rock[8] == 0.){ 					
	// 	SS_ref_db.z_em[6]          = 0.0;
    //     SS_ref_db.d_em[6]          = 1.0;
	// 	SS_ref_db.bounds_ref[5][0] = 0.0; 
	// 	SS_ref_db.bounds_ref[5][1] = 0.0;	
	// }
    return SS_ref_db;
}

/**
   retrieve reference thermodynamic data for mb_hb
*/
SS_ref G_SS_mb_hb_function(SS_ref SS_ref_db, int EM_database, int len_ox, bulk_info z_b, double eps){
    
    int i, j;
    int n_em = SS_ref_db.n_em;
    
    char   *EM_tmp[] 		= {"tr","tsm","prgm","glm","cumm","grnm","a","b","mrb","kprg","tts"};
    for (int i = 0; i < SS_ref_db.n_em; i++){
        strcpy(SS_ref_db.EM_list[i],EM_tmp[i]);
    };
    int n_xeos = SS_ref_db.n_xeos;
    char   *CV_tmp[] 		= {"x","y","z","a","k","c","f","t","Q1","Q2"};
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
        strcpy(SS_ref_db.CV_list[i],CV_tmp[i]);
    };
    SS_ref_db.W[0] = 20.0;
    SS_ref_db.W[1] = 25.0;
    SS_ref_db.W[2] = 65.0;
    SS_ref_db.W[3] = 45.0;
    SS_ref_db.W[4] = 75.0;
    SS_ref_db.W[5] = 57.0;
    SS_ref_db.W[6] = 63.0;
    SS_ref_db.W[7] = 52.0;
    SS_ref_db.W[8] = 30.0;
    SS_ref_db.W[9] = 85.0;
    SS_ref_db.W[10] = -40.0;
    SS_ref_db.W[11] = 25.0;
    SS_ref_db.W[12] = 70.0;
    SS_ref_db.W[13] = 80.0;
    SS_ref_db.W[14] = 70.0;
    SS_ref_db.W[15] = 72.5;
    SS_ref_db.W[16] = 20.0;
    SS_ref_db.W[17] = -40.0;
    SS_ref_db.W[18] = 35.0;
    SS_ref_db.W[19] = 50.0;
    SS_ref_db.W[20] = 90.0;
    SS_ref_db.W[21] = 106.700000000000;
    SS_ref_db.W[22] = 94.8;
    SS_ref_db.W[23] = 94.8;
    SS_ref_db.W[24] = 40.0;
    SS_ref_db.W[25] = 8.00;
    SS_ref_db.W[26] = 15.0;
    SS_ref_db.W[27] = 100.;
    SS_ref_db.W[28] = 113.500000000000;
    SS_ref_db.W[29] = 100.;
    SS_ref_db.W[30] = 111.200000000000;
    SS_ref_db.W[31] = 0.0;
    SS_ref_db.W[32] = 54.0;
    SS_ref_db.W[33] = 75.0;
    SS_ref_db.W[34] = 33.0;
    SS_ref_db.W[35] = 18.0;
    SS_ref_db.W[36] = 23.0;
    SS_ref_db.W[37] = 80.0;
    SS_ref_db.W[38] = 87.0;
    SS_ref_db.W[39] = 100.;
    SS_ref_db.W[40] = 12.0;
    SS_ref_db.W[41] = 8.00;
    SS_ref_db.W[42] = 91.0;
    SS_ref_db.W[43] = 96.0;
    SS_ref_db.W[44] = 65.0;
    SS_ref_db.W[45] = 20.0;
    SS_ref_db.W[46] = 80.0;
    SS_ref_db.W[47] = 94.0;
    SS_ref_db.W[48] = 95.0;
    SS_ref_db.W[49] = 90.0;
    SS_ref_db.W[50] = 94.0;
    SS_ref_db.W[51] = 95.0;
    SS_ref_db.W[52] = 50.0;
    SS_ref_db.W[53] = 50.0;
    SS_ref_db.W[54] = 35.0;
    
    SS_ref_db.v[0] = 1.00;
    SS_ref_db.v[1] = 1.50;
    SS_ref_db.v[2] = 1.70;
    SS_ref_db.v[3] = 0.800;
    SS_ref_db.v[4] = 1.00;
    SS_ref_db.v[5] = 1.00;
    SS_ref_db.v[6] = 1.00;
    SS_ref_db.v[7] = 1.00;
    SS_ref_db.v[8] = 0.800;
    SS_ref_db.v[9] = 1.70;
    SS_ref_db.v[10] = 1.50;
    
    
    em_data tr_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"tr", 
    										"equilibrium"	);
    
    em_data ts_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"ts", 
    										"equilibrium"	);
    
    em_data parg_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"parg", 
    										"equilibrium"	);
    
    em_data gl_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"gl", 
    										"equilibrium"	);
    
    em_data cumm_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"cumm", 
    										"equilibrium"	);
    
    em_data grun_eq 	= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"grun", 
    										"equilibrium"	);
    
    em_data gr_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"gr", 
    										"equilibrium"	);
    
    em_data andr_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"andr", 
    										"equilibrium"	);
    
    em_data mu_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"mu", 
    										"equilibrium"	);
    
    em_data pa_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"pa", 
    										"equilibrium"	);
    
    em_data dsp_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"dsp", 
    										"equilibrium"	);
    
    em_data ru_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"ru", 
    										"equilibrium"	);
    
    SS_ref_db.gbase[0] 		= tr_eq.gb;
    SS_ref_db.gbase[1] 		= ts_eq.gb + 10.0;
    SS_ref_db.gbase[2] 		= parg_eq.gb - 10.0;
    SS_ref_db.gbase[3] 		= gl_eq.gb - 3.0;
    SS_ref_db.gbase[4] 		= cumm_eq.gb;
    SS_ref_db.gbase[5] 		= grun_eq.gb - 3.0;
    SS_ref_db.gbase[6] 		= 0.428571428571429*cumm_eq.gb + 0.571428571428571*grun_eq.gb - 11.2;
    SS_ref_db.gbase[7] 		= 0.285714285714286*cumm_eq.gb + 0.714285714285714*grun_eq.gb - 13.8;
    SS_ref_db.gbase[8] 		= andr_eq.gb + gl_eq.gb -gr_eq.gb;
    SS_ref_db.gbase[9] 		= 0.02*z_b.T + mu_eq.gb -pa_eq.gb + parg_eq.gb - 7.06;
    SS_ref_db.gbase[10] 		= -2.0*dsp_eq.gb + 2.0*ru_eq.gb + ts_eq.gb + 95.0;
    
    SS_ref_db.ElShearMod[0] 	= tr_eq.ElShearMod;
    SS_ref_db.ElShearMod[1] 	= ts_eq.ElShearMod;
    SS_ref_db.ElShearMod[2] 	= parg_eq.ElShearMod;
    SS_ref_db.ElShearMod[3] 	= gl_eq.ElShearMod;
    SS_ref_db.ElShearMod[4] 	= cumm_eq.ElShearMod;
    SS_ref_db.ElShearMod[5] 	= grun_eq.ElShearMod;
    SS_ref_db.ElShearMod[6] 	= 0.428571428571429*cumm_eq.ElShearMod + 0.571428571428571*grun_eq.ElShearMod;
    SS_ref_db.ElShearMod[7] 	= 0.285714285714286*cumm_eq.ElShearMod + 0.714285714285714*grun_eq.ElShearMod;
    SS_ref_db.ElShearMod[8] 	= andr_eq.ElShearMod + gl_eq.ElShearMod -gr_eq.ElShearMod;
    SS_ref_db.ElShearMod[9] 	= mu_eq.ElShearMod -pa_eq.ElShearMod + parg_eq.ElShearMod;
    SS_ref_db.ElShearMod[10] 	= -2.0*dsp_eq.ElShearMod + 2.0*ru_eq.ElShearMod + ts_eq.ElShearMod;
    
    for (i = 0; i < len_ox; i++){
        SS_ref_db.Comp[0][i] 	= tr_eq.C[i];
        SS_ref_db.Comp[1][i] 	= ts_eq.C[i];
        SS_ref_db.Comp[2][i] 	= parg_eq.C[i];
        SS_ref_db.Comp[3][i] 	= gl_eq.C[i];
        SS_ref_db.Comp[4][i] 	= cumm_eq.C[i];
        SS_ref_db.Comp[5][i] 	= grun_eq.C[i];
        SS_ref_db.Comp[6][i] 	= 0.428571428571429*cumm_eq.C[i] + 0.571428571428571*grun_eq.C[i];
        SS_ref_db.Comp[7][i] 	= 0.285714285714286*cumm_eq.C[i] + 0.714285714285714*grun_eq.C[i];
        SS_ref_db.Comp[8][i] 	= andr_eq.C[i] + gl_eq.C[i] -gr_eq.C[i];
        SS_ref_db.Comp[9][i] 	= mu_eq.C[i] -pa_eq.C[i] + parg_eq.C[i];
        SS_ref_db.Comp[10][i] 	= -2.0*dsp_eq.C[i] + 2.0*ru_eq.C[i] + ts_eq.C[i];
    }
    
    for (i = 0; i < n_em; i++){
        SS_ref_db.z_em[i] = 1.0;
    };
    
    SS_ref_db.bounds_ref[0][0] = 0.0+eps;  SS_ref_db.bounds_ref[0][1] = 1.0-eps;
    SS_ref_db.bounds_ref[1][0] = 0.0+eps;  SS_ref_db.bounds_ref[1][1] = 1.0-eps;
    SS_ref_db.bounds_ref[2][0] = 0.0+eps;  SS_ref_db.bounds_ref[2][1] = 1.0-eps;
    SS_ref_db.bounds_ref[3][0] = 0.0+eps;  SS_ref_db.bounds_ref[3][1] = 1.0-eps;
    SS_ref_db.bounds_ref[4][0] = 0.0+eps;  SS_ref_db.bounds_ref[4][1] = 1.0-eps;
    SS_ref_db.bounds_ref[5][0] = 0.0+eps;  SS_ref_db.bounds_ref[5][1] = 1.0-eps;
    SS_ref_db.bounds_ref[6][0] = 0.0+eps;  SS_ref_db.bounds_ref[6][1] = 1.0-eps;
    SS_ref_db.bounds_ref[7][0] = 0.0+eps;  SS_ref_db.bounds_ref[7][1] = 1.0-eps;
    SS_ref_db.bounds_ref[8][0] = -1.0+eps;  SS_ref_db.bounds_ref[8][1] = 1.0-eps;
    SS_ref_db.bounds_ref[9][0] = -1.0+eps;  SS_ref_db.bounds_ref[9][1] = 1.0-eps;
    


	if (z_b.bulk_rock[8] == 0.){ 	    //O				
		SS_ref_db.z_em[8]          = 0.0;
        SS_ref_db.d_em[8]          = 1.0;
		SS_ref_db.bounds_ref[6][0] = 0.0; 
		SS_ref_db.bounds_ref[6][1] = 0.0;	
	}
	if (z_b.bulk_rock[7] == 0.){ 		//TiO2	
		SS_ref_db.z_em[10]          = 0.0;
        SS_ref_db.d_em[10]          = 1.0;
		SS_ref_db.bounds_ref[7][0]  = 0.0; 
		SS_ref_db.bounds_ref[7][1]  = 0.0;	
	}

    return SS_ref_db;
}

/**
   retrieve reference thermodynamic data for mb_aug
*/
SS_ref G_SS_mb_aug_function(SS_ref SS_ref_db, int EM_database, int len_ox, bulk_info z_b, double eps){
    
    int i, j;
    int n_em = SS_ref_db.n_em;
    
    char   *EM_tmp[] 		= {"di","cenh","cfs","jdm","acmm","ocats","dcats","fmc"};
    for (int i = 0; i < SS_ref_db.n_em; i++){
        strcpy(SS_ref_db.EM_list[i],EM_tmp[i]);
    };
    int n_xeos = SS_ref_db.n_xeos;
    char   *CV_tmp[] 		= {"x","y","f","z","j","Qfm","Qa1"};
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
        strcpy(SS_ref_db.CV_list[i],CV_tmp[i]);
    };
    SS_ref_db.W[0] = 29.8 - 0.03*SS_ref_db.P;
    SS_ref_db.W[1] = 25.8 - 0.03*SS_ref_db.P;
    SS_ref_db.W[2] = 26.0;
    SS_ref_db.W[3] = 21.0;
    SS_ref_db.W[4] = 12.3 - 0.01*SS_ref_db.P;
    SS_ref_db.W[5] = 12.3 - 0.01*SS_ref_db.P;
    SS_ref_db.W[6] = 20.6 - 0.03*SS_ref_db.P;
    SS_ref_db.W[7] = 2.30;
    SS_ref_db.W[8] = 50.0;
    SS_ref_db.W[9] = 62.0;
    SS_ref_db.W[10] = 45.7 - 0.29*SS_ref_db.P;
    SS_ref_db.W[11] = 45.7 - 0.29*SS_ref_db.P;
    SS_ref_db.W[12] = 4.00;
    SS_ref_db.W[13] = 60.0;
    SS_ref_db.W[14] = 58.0;
    SS_ref_db.W[15] = 48.0;
    SS_ref_db.W[16] = 48.0;
    SS_ref_db.W[17] = 3.50;
    SS_ref_db.W[18] = 5.00;
    SS_ref_db.W[19] = 40.0;
    SS_ref_db.W[20] = 40.0;
    SS_ref_db.W[21] = 40.0;
    SS_ref_db.W[22] = 35.0;
    SS_ref_db.W[23] = 35.0;
    SS_ref_db.W[24] = 60.0;
    SS_ref_db.W[25] = 0.01*SS_ref_db.P + 3.8;
    SS_ref_db.W[26] = 50.0;
    SS_ref_db.W[27] = 50.0;
    
    SS_ref_db.v[0] = 1.20;
    SS_ref_db.v[1] = 1.00;
    SS_ref_db.v[2] = 1.00;
    SS_ref_db.v[3] = 1.20;
    SS_ref_db.v[4] = 1.20;
    SS_ref_db.v[5] = 1.90;
    SS_ref_db.v[6] = 1.90;
    SS_ref_db.v[7] = 1.00;
    
    
    em_data di_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"di", 
    										"equilibrium"	);
    
    em_data en_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"en", 
    										"equilibrium"	);
    
    em_data fs_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"fs", 
    										"equilibrium"	);
    
    em_data jd_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"jd", 
    										"equilibrium"	);
    
    em_data acm_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"acm", 
    										"equilibrium"	);
    
    em_data cats_or 	= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"cats", 
    										"ordered"	);
    
    SS_ref_db.gbase[0] 		= di_eq.gb;
    SS_ref_db.gbase[1] 		= 0.048*z_b.P - 0.002*z_b.T + en_eq.gb + 3.5;
    SS_ref_db.gbase[2] 		= 0.045*z_b.P - 0.002*z_b.T + fs_eq.gb + 2.1;
    SS_ref_db.gbase[3] 		= jd_eq.gb + 2.0;
    SS_ref_db.gbase[4] 		= acm_eq.gb - 5.0;
    SS_ref_db.gbase[5] 		= cats_or.gb;
    SS_ref_db.gbase[6] 		= 0.01*z_b.P - 0.002882*z_b.T + cats_or.gb + 3.8;
    SS_ref_db.gbase[7] 		= 0.0465*z_b.P - 0.002*z_b.T + 0.5*en_eq.gb + 0.5*fs_eq.gb - 1.6;
    
    SS_ref_db.ElShearMod[0] 	= di_eq.ElShearMod;
    SS_ref_db.ElShearMod[1] 	= en_eq.ElShearMod;
    SS_ref_db.ElShearMod[2] 	= fs_eq.ElShearMod;
    SS_ref_db.ElShearMod[3] 	= jd_eq.ElShearMod;
    SS_ref_db.ElShearMod[4] 	= acm_eq.ElShearMod;
    SS_ref_db.ElShearMod[5] 	= cats_or.ElShearMod;
    SS_ref_db.ElShearMod[6] 	= cats_or.ElShearMod;
    SS_ref_db.ElShearMod[7] 	= 0.5*en_eq.ElShearMod + 0.5*fs_eq.ElShearMod;
    
    for (i = 0; i < len_ox; i++){
        SS_ref_db.Comp[0][i] 	= di_eq.C[i];
        SS_ref_db.Comp[1][i] 	= en_eq.C[i];
        SS_ref_db.Comp[2][i] 	= fs_eq.C[i];
        SS_ref_db.Comp[3][i] 	= jd_eq.C[i];
        SS_ref_db.Comp[4][i] 	= acm_eq.C[i];
        SS_ref_db.Comp[5][i] 	= cats_or.C[i];
        SS_ref_db.Comp[6][i] 	= cats_or.C[i];
        SS_ref_db.Comp[7][i] 	= 0.5*en_eq.C[i] + 0.5*fs_eq.C[i];
    }
    
    for (i = 0; i < n_em; i++){
        SS_ref_db.z_em[i] = 1.0;
    };
    
    SS_ref_db.bounds_ref[0][0] = 0.0+eps;  SS_ref_db.bounds_ref[0][1] = 1.0-eps;
    SS_ref_db.bounds_ref[1][0] = 0.0+eps;  SS_ref_db.bounds_ref[1][1] = 1.0-eps;
    SS_ref_db.bounds_ref[2][0] = 0.0+eps;  SS_ref_db.bounds_ref[2][1] = 1.0-eps;
    SS_ref_db.bounds_ref[3][0] = 0.0+eps;  SS_ref_db.bounds_ref[3][1] = 1.0-eps;
    SS_ref_db.bounds_ref[4][0] = 0.0+eps;  SS_ref_db.bounds_ref[4][1] = 1.0-eps;
    SS_ref_db.bounds_ref[5][0] = 0.0+eps;  SS_ref_db.bounds_ref[5][1] = 2.0-eps;
    SS_ref_db.bounds_ref[6][0] = 0.0+eps;  SS_ref_db.bounds_ref[6][1] = 1.0-eps;
    
	if (z_b.bulk_rock[8] == 0.){ 	    //O				
		SS_ref_db.z_em[5]          = 0.0;
        SS_ref_db.d_em[4]          = 1.0;
		SS_ref_db.bounds_ref[3][0] = 0.0; 
		SS_ref_db.bounds_ref[3][1] = 0.0;	
	}


    return SS_ref_db;
}

/**
   retrieve reference thermodynamic data for mb_dio
*/
SS_ref G_SS_mb_dio_function(SS_ref SS_ref_db, int EM_database, int len_ox, bulk_info z_b, double eps){
    
    int i, j;
    int n_em = SS_ref_db.n_em;
    
    char   *EM_tmp[] 		= {"jd","di","hed","acmm","om","cfm","jac"};
    for (int i = 0; i < SS_ref_db.n_em; i++){
        strcpy(SS_ref_db.EM_list[i],EM_tmp[i]);
    };
    int n_xeos = SS_ref_db.n_xeos;
    char   *CV_tmp[] 		= {"x","j","t","c","Qaf","Qfm"};
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
        strcpy(SS_ref_db.CV_list[i],CV_tmp[i]);
    };
    SS_ref_db.W[0] = 26.0;
    SS_ref_db.W[1] = 24.0;
    SS_ref_db.W[2] = 5.00;
    SS_ref_db.W[3] = 15.5;
    SS_ref_db.W[4] = 25.2;
    SS_ref_db.W[5] = 3.00;
    SS_ref_db.W[6] = 4.00;
    SS_ref_db.W[7] = 21.0;
    SS_ref_db.W[8] = 15.7500000000000;
    SS_ref_db.W[9] = 2.00;
    SS_ref_db.W[10] = 24.6500000000000;
    SS_ref_db.W[11] = 20.8;
    SS_ref_db.W[12] = 17.2;
    SS_ref_db.W[13] = 2.00;
    SS_ref_db.W[14] = 24.6;
    SS_ref_db.W[15] = 16.4;
    SS_ref_db.W[16] = 22.2;
    SS_ref_db.W[17] = 3.00;
    SS_ref_db.W[18] = 18.4500000000000;
    SS_ref_db.W[19] = 19.5;
    SS_ref_db.W[20] = 24.5500000000000;
    
    
    em_data jd_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"jd", 
    										"equilibrium"	);
    
    em_data di_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"di", 
    										"equilibrium"	);
    
    em_data hed_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"hed", 
    										"equilibrium"	);
    
    em_data acm_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"acm", 
    										"equilibrium"	);
    
    SS_ref_db.gbase[0] 		= jd_eq.gb;
    SS_ref_db.gbase[1] 		= di_eq.gb;
    SS_ref_db.gbase[2] 		= hed_eq.gb;
    SS_ref_db.gbase[3] 		= acm_eq.gb - 7.0;
    SS_ref_db.gbase[4] 		= 0.5*di_eq.gb + 0.5*jd_eq.gb - 2.9;
    SS_ref_db.gbase[5] 		= 0.5*di_eq.gb + 0.5*hed_eq.gb - 1.5;
    SS_ref_db.gbase[6] 		= 0.5*acm_eq.gb + 0.5*jd_eq.gb - 4.5;
    
    SS_ref_db.ElShearMod[0] 	= jd_eq.ElShearMod;
    SS_ref_db.ElShearMod[1] 	= di_eq.ElShearMod;
    SS_ref_db.ElShearMod[2] 	= hed_eq.ElShearMod;
    SS_ref_db.ElShearMod[3] 	= acm_eq.ElShearMod;
    SS_ref_db.ElShearMod[4] 	= 0.5*di_eq.ElShearMod + 0.5*jd_eq.ElShearMod;
    SS_ref_db.ElShearMod[5] 	= 0.5*di_eq.ElShearMod + 0.5*hed_eq.ElShearMod;
    SS_ref_db.ElShearMod[6] 	= 0.5*acm_eq.ElShearMod + 0.5*jd_eq.ElShearMod;
    
    for (i = 0; i < len_ox; i++){
        SS_ref_db.Comp[0][i] 	= jd_eq.C[i];
        SS_ref_db.Comp[1][i] 	= di_eq.C[i];
        SS_ref_db.Comp[2][i] 	= hed_eq.C[i];
        SS_ref_db.Comp[3][i] 	= acm_eq.C[i];
        SS_ref_db.Comp[4][i] 	= 0.5*di_eq.C[i] + 0.5*jd_eq.C[i];
        SS_ref_db.Comp[5][i] 	= 0.5*di_eq.C[i] + 0.5*hed_eq.C[i];
        SS_ref_db.Comp[6][i] 	= 0.5*acm_eq.C[i] + 0.5*jd_eq.C[i];
    }
    
    for (i = 0; i < n_em; i++){
        SS_ref_db.z_em[i] = 1.0;
    };
    
    SS_ref_db.bounds_ref[0][0] = 0.0+eps;  SS_ref_db.bounds_ref[0][1] = 1.0-eps;
    SS_ref_db.bounds_ref[1][0] = 0.0+eps;  SS_ref_db.bounds_ref[1][1] = 1.0-eps;
    SS_ref_db.bounds_ref[2][0] = 0.0+eps;  SS_ref_db.bounds_ref[2][1] = 1.0-eps;
    SS_ref_db.bounds_ref[3][0] = -0.5+eps;  SS_ref_db.bounds_ref[3][1] = 0.5-eps;
    SS_ref_db.bounds_ref[4][0] = -0.5+eps;  SS_ref_db.bounds_ref[4][1] = 0.5-eps;
    SS_ref_db.bounds_ref[5][0] = -0.5+eps;  SS_ref_db.bounds_ref[5][1] = 0.5-eps;

    if (z_b.bulk_rock[8] == 0.){ 	    //O				
		SS_ref_db.z_em[3]          = 0.0;
        SS_ref_db.d_em[3]          = 1.0;
		SS_ref_db.bounds_ref[2][0] = 0.0; 
		SS_ref_db.bounds_ref[2][1] = 0.0;	

		SS_ref_db.z_em[6]          = 0.0;
        SS_ref_db.d_em[6]          = 1.0;
		SS_ref_db.bounds_ref[4][0] = 0.0; 
		SS_ref_db.bounds_ref[4][1] = 0.0;	
	}

    return SS_ref_db;
}

/**
   retrieve reference thermodynamic data for mb_opx
*/
SS_ref G_SS_mb_opx_function(SS_ref SS_ref_db, int EM_database, int len_ox, bulk_info z_b, double eps){
    
    int i, j;
    int n_em = SS_ref_db.n_em;
    
    char   *EM_tmp[] 		= {"en","fs","fm","mgts","fopx","odi"};
    for (int i = 0; i < SS_ref_db.n_em; i++){
        strcpy(SS_ref_db.EM_list[i],EM_tmp[i]);
    };
    int n_xeos = SS_ref_db.n_xeos;
    char   *CV_tmp[] 		= {"x","y","f","c","Q"};
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
        strcpy(SS_ref_db.CV_list[i],CV_tmp[i]);
    };
    SS_ref_db.W[0] = 7.00;
    SS_ref_db.W[1] = 4.00;
    SS_ref_db.W[2] = 13.0 - 0.15*SS_ref_db.P;
    SS_ref_db.W[3] = 11.0 - 0.15*SS_ref_db.P;
    SS_ref_db.W[4] = 0.12*SS_ref_db.P + 32.2;
    SS_ref_db.W[5] = 4.00;
    SS_ref_db.W[6] = 13.0 - 0.15*SS_ref_db.P;
    SS_ref_db.W[7] = 11.6 - 0.15*SS_ref_db.P;
    SS_ref_db.W[8] = 0.084*SS_ref_db.P + 25.54;
    SS_ref_db.W[9] = 17.0 - 0.15*SS_ref_db.P;
    SS_ref_db.W[10] = 15.0 - 0.15*SS_ref_db.P;
    SS_ref_db.W[11] = 0.084*SS_ref_db.P + 22.54;
    SS_ref_db.W[12] = 1.00;
    SS_ref_db.W[13] = 75.4 - 0.94*SS_ref_db.P;
    SS_ref_db.W[14] = 73.4 - 0.94*SS_ref_db.P;
    
    SS_ref_db.v[0] = 1.00;
    SS_ref_db.v[1] = 1.00;
    SS_ref_db.v[2] = 1.00;
    SS_ref_db.v[3] = 1.00;
    SS_ref_db.v[4] = 1.00;
    SS_ref_db.v[5] = 1.20;
    
    
    em_data en_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"en", 
    										"equilibrium"	);
    
    em_data fs_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"fs", 
    										"equilibrium"	);
    
    em_data mgts_eq 	= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"mgts", 
    										"equilibrium"	);
    
    em_data gr_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"gr", 
    										"equilibrium"	);
    
    em_data andr_eq 	= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"andr", 
    										"equilibrium"	);
    
    em_data di_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"di", 
    										"equilibrium"	);
    
    SS_ref_db.gbase[0] 		= en_eq.gb;
    SS_ref_db.gbase[1] 		= fs_eq.gb;
    SS_ref_db.gbase[2] 		= 0.5*en_eq.gb + 0.5*fs_eq.gb - 6.6;
    SS_ref_db.gbase[3] 		= mgts_eq.gb;
    SS_ref_db.gbase[4] 		= 0.5*andr_eq.gb - 0.5*gr_eq.gb + mgts_eq.gb + 2.0;
    SS_ref_db.gbase[5] 		= 0.005*z_b.P + 0.000211*z_b.T + di_eq.gb - 0.1;
    
    SS_ref_db.ElShearMod[0] 	= en_eq.ElShearMod;
    SS_ref_db.ElShearMod[1] 	= fs_eq.ElShearMod;
    SS_ref_db.ElShearMod[2] 	= 0.5*en_eq.ElShearMod + 0.5*fs_eq.ElShearMod;
    SS_ref_db.ElShearMod[3] 	= mgts_eq.ElShearMod;
    SS_ref_db.ElShearMod[4] 	= 0.5*andr_eq.ElShearMod - 0.5*gr_eq.ElShearMod + mgts_eq.ElShearMod;
    SS_ref_db.ElShearMod[5] 	= di_eq.ElShearMod;
    
    for (i = 0; i < len_ox; i++){
        SS_ref_db.Comp[0][i] 	= en_eq.C[i];
        SS_ref_db.Comp[1][i] 	= fs_eq.C[i];
        SS_ref_db.Comp[2][i] 	= 0.5*en_eq.C[i] + 0.5*fs_eq.C[i];
        SS_ref_db.Comp[3][i] 	= mgts_eq.C[i];
        SS_ref_db.Comp[4][i] 	= 0.5*andr_eq.C[i] - 0.5*gr_eq.C[i] + mgts_eq.C[i];
        SS_ref_db.Comp[5][i] 	= di_eq.C[i];
    }
    
    for (i = 0; i < n_em; i++){
        SS_ref_db.z_em[i] = 1.0;
    };
    
    SS_ref_db.bounds_ref[0][0] = 0.0+eps;  SS_ref_db.bounds_ref[0][1] = 1.0-eps;
    SS_ref_db.bounds_ref[1][0] = 0.0+eps;  SS_ref_db.bounds_ref[1][1] = 1.0-eps;
    SS_ref_db.bounds_ref[2][0] = 0.0+eps;  SS_ref_db.bounds_ref[2][1] = 1.0-eps;
    SS_ref_db.bounds_ref[3][0] = 0.0+eps;  SS_ref_db.bounds_ref[3][1] = 1.0-eps;
    SS_ref_db.bounds_ref[4][0] = 0.0+eps;  SS_ref_db.bounds_ref[4][1] = 1.0-eps;
    

    if (z_b.bulk_rock[8] == 0.){ 	    //O				
		SS_ref_db.z_em[4]          = 0.0;
        SS_ref_db.d_em[4]          = 1.0;
		SS_ref_db.bounds_ref[2][0] = 0.0; 
		SS_ref_db.bounds_ref[2][1] = 0.0;	
	}

    return SS_ref_db;
}

/**
   retrieve reference thermodynamic data for mb_g
*/
SS_ref G_SS_mb_g_function(SS_ref SS_ref_db, int EM_database, int len_ox, bulk_info z_b, double eps){
    
    int i, j;
    int n_em = SS_ref_db.n_em;
    
    char   *EM_tmp[] 		= {"py","alm","gr","kho"};
    for (int i = 0; i < SS_ref_db.n_em; i++){
        strcpy(SS_ref_db.EM_list[i],EM_tmp[i]);
    };
    int n_xeos = SS_ref_db.n_xeos;
    char   *CV_tmp[] 		= {"x","z","f"};
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
        strcpy(SS_ref_db.CV_list[i],CV_tmp[i]);
    };
    SS_ref_db.W[0] = 2.50;
    SS_ref_db.W[1] = 31.0;
    SS_ref_db.W[2] = 5.40;
    SS_ref_db.W[3] = 5.00;
    SS_ref_db.W[4] = 22.6;
    SS_ref_db.W[5] = -15.3;
    
    SS_ref_db.v[0] = 1.00;
    SS_ref_db.v[1] = 1.00;
    SS_ref_db.v[2] = 2.70;
    SS_ref_db.v[3] = 1.00;
    
    
    em_data py_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"py", 
    										"equilibrium"	);
    
    em_data alm_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"alm", 
    										"equilibrium"	);
    
    em_data gr_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"gr", 
    										"equilibrium"	);
    
    em_data andr_eq 	= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"andr", 
    										"equilibrium"	);
    
    SS_ref_db.gbase[0] 		= py_eq.gb;
    SS_ref_db.gbase[1] 		= alm_eq.gb;
    SS_ref_db.gbase[2] 		= gr_eq.gb;
    SS_ref_db.gbase[3] 		= andr_eq.gb -gr_eq.gb + py_eq.gb + 27.0;
    
    SS_ref_db.ElShearMod[0] 	= py_eq.ElShearMod;
    SS_ref_db.ElShearMod[1] 	= alm_eq.ElShearMod;
    SS_ref_db.ElShearMod[2] 	= gr_eq.ElShearMod;
    SS_ref_db.ElShearMod[3] 	= andr_eq.ElShearMod -gr_eq.ElShearMod + py_eq.ElShearMod;
    
    for (i = 0; i < len_ox; i++){
        SS_ref_db.Comp[0][i] 	= py_eq.C[i];
        SS_ref_db.Comp[1][i] 	= alm_eq.C[i];
        SS_ref_db.Comp[2][i] 	= gr_eq.C[i];
        SS_ref_db.Comp[3][i] 	= andr_eq.C[i] -gr_eq.C[i] + py_eq.C[i];
    }
    
    for (i = 0; i < n_em; i++){
        SS_ref_db.z_em[i] = 1.0;
    };
    
    SS_ref_db.bounds_ref[0][0] = 0.0+eps;  SS_ref_db.bounds_ref[0][1] = 1.0-eps;
    SS_ref_db.bounds_ref[1][0] = 0.0+eps;  SS_ref_db.bounds_ref[1][1] = 1.0-eps;
    SS_ref_db.bounds_ref[2][0] = 0.0+eps;  SS_ref_db.bounds_ref[2][1] = 1.0-eps;
    
    if (z_b.bulk_rock[8] == 0.){ 	    //O				
		SS_ref_db.z_em[3]          = 0.0;
        SS_ref_db.d_em[3]          = 1.0;
		SS_ref_db.bounds_ref[2][0] = 0.0; 
		SS_ref_db.bounds_ref[2][1] = 0.0;	
	}


    return SS_ref_db;
}

/**
   retrieve reference thermodynamic data for mb_ol
*/
SS_ref G_SS_mb_ol_function(SS_ref SS_ref_db, int EM_database, int len_ox, bulk_info z_b, double eps){
    
    int i, j;
    int n_em = SS_ref_db.n_em;
    
    char   *EM_tmp[] 		= {"fo","fa"};
    for (int i = 0; i < SS_ref_db.n_em; i++){
        strcpy(SS_ref_db.EM_list[i],EM_tmp[i]);
    };
    int n_xeos = SS_ref_db.n_xeos;
    char   *CV_tmp[] 		= {"x"};
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
        strcpy(SS_ref_db.CV_list[i],CV_tmp[i]);
    };
    SS_ref_db.W[0] = 9.00;
    
    
    em_data fo_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"fo", 
    										"equilibrium"	);
    
    em_data fa_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"fa", 
    										"equilibrium"	);
    
    SS_ref_db.gbase[0] 		= fo_eq.gb;
    SS_ref_db.gbase[1] 		= fa_eq.gb;
    
    SS_ref_db.ElShearMod[0] 	= fo_eq.ElShearMod;
    SS_ref_db.ElShearMod[1] 	= fa_eq.ElShearMod;
    
    for (i = 0; i < len_ox; i++){
        SS_ref_db.Comp[0][i] 	= fo_eq.C[i];
        SS_ref_db.Comp[1][i] 	= fa_eq.C[i];
    }
    
    for (i = 0; i < n_em; i++){
        SS_ref_db.z_em[i] = 1.0;
    };
    
    SS_ref_db.bounds_ref[0][0] = 0.0+eps;  SS_ref_db.bounds_ref[0][1] = 1.0-eps;
    
    return SS_ref_db;
}

/**
   retrieve reference thermodynamic data for mb_fsp
*/
SS_ref G_SS_mb_fsp_function(SS_ref SS_ref_db, int EM_database, int len_ox, bulk_info z_b, double eps){
    
    int i, j;
    int n_em = SS_ref_db.n_em;
    
    char   *EM_tmp[] 		= {"ab","an","san"};
    for (int i = 0; i < SS_ref_db.n_em; i++){
        strcpy(SS_ref_db.EM_list[i],EM_tmp[i]);
    };
    int n_xeos = SS_ref_db.n_xeos;
    char   *CV_tmp[] 		= {"ca","k"};
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
        strcpy(SS_ref_db.CV_list[i],CV_tmp[i]);
    };
    SS_ref_db.W[0] = -0.04*SS_ref_db.P - 0.00935*SS_ref_db.T + 14.6;
    SS_ref_db.W[1] = 0.338*SS_ref_db.P - 0.00957*SS_ref_db.T + 24.1;
    SS_ref_db.W[2] = 48.5 - 0.13*SS_ref_db.P;
    
    SS_ref_db.v[0] = 0.674;
    SS_ref_db.v[1] = 0.550;
    SS_ref_db.v[2] = 1.00;
    
    
    em_data ab_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"ab", 
    										"equilibrium"	);
    
    em_data an_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"an", 
    										"equilibrium"	);
    
    em_data san_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"san", 
    										"equilibrium"	);
    
    SS_ref_db.gbase[0] 		= ab_eq.gb;
    SS_ref_db.gbase[1] 		= an_eq.gb;
    SS_ref_db.gbase[2] 		= san_eq.gb;
    
    SS_ref_db.ElShearMod[0] 	= ab_eq.ElShearMod;
    SS_ref_db.ElShearMod[1] 	= an_eq.ElShearMod;
    SS_ref_db.ElShearMod[2] 	= san_eq.ElShearMod;
    
    for (i = 0; i < len_ox; i++){
        SS_ref_db.Comp[0][i] 	= ab_eq.C[i];
        SS_ref_db.Comp[1][i] 	= an_eq.C[i];
        SS_ref_db.Comp[2][i] 	= san_eq.C[i];
    }
    
    for (i = 0; i < n_em; i++){
        SS_ref_db.z_em[i] = 1.0;
    };
    
    SS_ref_db.bounds_ref[0][0] = 0.0+eps;  SS_ref_db.bounds_ref[0][1] = 1.0-eps;
    SS_ref_db.bounds_ref[1][0] = 0.0+eps;  SS_ref_db.bounds_ref[1][1] = 1.0-eps;
    
    return SS_ref_db;
}

/**
   retrieve reference thermodynamic data for mb_abc
*/
SS_ref G_SS_mb_abc_function(SS_ref SS_ref_db, int EM_database, int len_ox, bulk_info z_b, double eps){
    
    int i, j;
    int n_em = SS_ref_db.n_em;
    
    char   *EM_tmp[] 		= {"abm","anm"};
    for (int i = 0; i < SS_ref_db.n_em; i++){
        strcpy(SS_ref_db.EM_list[i],EM_tmp[i]);
    };
    int n_xeos = SS_ref_db.n_xeos;
    char   *CV_tmp[] 		= {"ca"};
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
        strcpy(SS_ref_db.CV_list[i],CV_tmp[i]);
    };
    SS_ref_db.W[0] = 3.40;
    
    SS_ref_db.v[0] = 0.640;
    SS_ref_db.v[1] = 1.00;
    
    
    em_data ab_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"ab", 
    										"equilibrium"	);
    
    em_data an_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"an", 
    										"equilibrium"	);
    
    SS_ref_db.gbase[0] 		= 0.002*z_b.T + ab_eq.gb - 1.746;
    SS_ref_db.gbase[1] 		= an_eq.gb + 10.0;
    
    SS_ref_db.ElShearMod[0] 	= ab_eq.ElShearMod;
    SS_ref_db.ElShearMod[1] 	= an_eq.ElShearMod;
    
    for (i = 0; i < len_ox; i++){
        SS_ref_db.Comp[0][i] 	= ab_eq.C[i];
        SS_ref_db.Comp[1][i] 	= an_eq.C[i];
    }
    
    for (i = 0; i < n_em; i++){
        SS_ref_db.z_em[i] = 1.0;
    };
    
    SS_ref_db.bounds_ref[0][0] = 0.0+eps;  SS_ref_db.bounds_ref[0][1] = 1.0-eps;
    
    return SS_ref_db;
}

/**
   retrieve reference thermodynamic data for mb_k4tr
*/
SS_ref G_SS_mb_k4tr_function(SS_ref SS_ref_db, int EM_database, int len_ox, bulk_info z_b, double eps){
    
    int i, j;
    int n_em = SS_ref_db.n_em;
    
    char   *EM_tmp[] 		= {"ab","an","san"};
    for (int i = 0; i < SS_ref_db.n_em; i++){
        strcpy(SS_ref_db.EM_list[i],EM_tmp[i]);
    };
    int n_xeos = SS_ref_db.n_xeos;
    char   *CV_tmp[] 		= {"na","ca"};
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
        strcpy(SS_ref_db.CV_list[i],CV_tmp[i]);
    };
    SS_ref_db.W[0] = -0.04*SS_ref_db.P - 0.00935*SS_ref_db.T + 14.6;
    SS_ref_db.W[1] = 0.338*SS_ref_db.P - 0.00957*SS_ref_db.T + 24.1;
    SS_ref_db.W[2] = 48.5 - 0.13*SS_ref_db.P;
    
    SS_ref_db.v[0] = 0.674;
    SS_ref_db.v[1] = 0.550;
    SS_ref_db.v[2] = 1.00;
    
    
    em_data ab_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"ab", 
    										"equilibrium"	);
    
    em_data an_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"an", 
    										"equilibrium"	);
    
    em_data san_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"san", 
    										"equilibrium"	);
    
    SS_ref_db.gbase[0] 		= ab_eq.gb;
    SS_ref_db.gbase[1] 		= an_eq.gb;
    SS_ref_db.gbase[2] 		= san_eq.gb;
    
    SS_ref_db.ElShearMod[0] 	= ab_eq.ElShearMod;
    SS_ref_db.ElShearMod[1] 	= an_eq.ElShearMod;
    SS_ref_db.ElShearMod[2] 	= san_eq.ElShearMod;
    
    for (i = 0; i < len_ox; i++){
        SS_ref_db.Comp[0][i] 	= ab_eq.C[i];
        SS_ref_db.Comp[1][i] 	= an_eq.C[i];
        SS_ref_db.Comp[2][i] 	= san_eq.C[i];
    }
    
    for (i = 0; i < n_em; i++){
        SS_ref_db.z_em[i] = 1.0;
    };
    
    SS_ref_db.bounds_ref[0][0] = 0.0+eps;  SS_ref_db.bounds_ref[0][1] = 1.0-eps;
    SS_ref_db.bounds_ref[1][0] = 0.0+eps;  SS_ref_db.bounds_ref[1][1] = 1.0-eps;
    
    return SS_ref_db;
}

/**
   retrieve reference thermodynamic data for mb_sp
*/
SS_ref G_SS_mb_sp_function(SS_ref SS_ref_db, int EM_database, int len_ox, bulk_info z_b, double eps){
    
    int i, j;
    int n_em = SS_ref_db.n_em;
    
    char   *EM_tmp[] 		= {"herc","sp","mt","usp"};
    for (int i = 0; i < SS_ref_db.n_em; i++){
        strcpy(SS_ref_db.EM_list[i],EM_tmp[i]);
    };
    int n_xeos = SS_ref_db.n_xeos;
    char   *CV_tmp[] 		= {"x","y","z"};
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
        strcpy(SS_ref_db.CV_list[i],CV_tmp[i]);
    };
    SS_ref_db.W[0] = 0.0;
    SS_ref_db.W[1] = 18.5;
    SS_ref_db.W[2] = 27.0;
    SS_ref_db.W[3] = 40.0;
    SS_ref_db.W[4] = 30.0;
    SS_ref_db.W[5] = 0.0;
    
    
    em_data herc_eq 	= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"herc", 
    										"equilibrium"	);
    
    em_data sp_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"sp", 
    										"equilibrium"	);
    
    em_data mt_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"mt", 
    										"equilibrium"	);
    
    em_data usp_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"usp", 
    										"equilibrium"	);
    
    SS_ref_db.gbase[0] 		= herc_eq.gb;
    SS_ref_db.gbase[1] 		= sp_eq.gb;
    SS_ref_db.gbase[2] 		= mt_eq.gb;
    SS_ref_db.gbase[3] 		= usp_eq.gb;
    
    SS_ref_db.ElShearMod[0] 	= herc_eq.ElShearMod;
    SS_ref_db.ElShearMod[1] 	= sp_eq.ElShearMod;
    SS_ref_db.ElShearMod[2] 	= mt_eq.ElShearMod;
    SS_ref_db.ElShearMod[3] 	= usp_eq.ElShearMod;
    
    for (i = 0; i < len_ox; i++){
        SS_ref_db.Comp[0][i] 	= herc_eq.C[i];
        SS_ref_db.Comp[1][i] 	= sp_eq.C[i];
        SS_ref_db.Comp[2][i] 	= mt_eq.C[i];
        SS_ref_db.Comp[3][i] 	= usp_eq.C[i];
    }
    
    for (i = 0; i < n_em; i++){
        SS_ref_db.z_em[i] = 1.0;
    };
    
    SS_ref_db.bounds_ref[0][0] = 0.0+eps;  SS_ref_db.bounds_ref[0][1] = 1.0-eps;
    SS_ref_db.bounds_ref[1][0] = 0.0+eps;  SS_ref_db.bounds_ref[1][1] = 1.0-eps;
    SS_ref_db.bounds_ref[2][0] = 0.0+eps;  SS_ref_db.bounds_ref[2][1] = 1.0-eps;
    

	if (z_b.bulk_rock[8] == 0. && z_b.bulk_rock[7] == 0.){ 	    //O	&& TiO2			
		SS_ref_db.z_em[2]          = 0.0;
        SS_ref_db.d_em[2]          = 1.0;
		SS_ref_db.z_em[3]          = 0.0;
        SS_ref_db.d_em[3]          = 1.0;

		SS_ref_db.bounds_ref[2][0] = 0.0; 
		SS_ref_db.bounds_ref[2][1] = 0.0;
		SS_ref_db.bounds_ref[1][0] = 1.0; 
		SS_ref_db.bounds_ref[1][1] = 1.0;	
	}
	if (z_b.bulk_rock[7] == 0. && z_b.bulk_rock[8] != 0.){ 		//TiO2	
		SS_ref_db.z_em[3]          = 0.0;
        SS_ref_db.d_em[3]          = 1.0;
		SS_ref_db.bounds_ref[2][0]  = 0.0; 
		SS_ref_db.bounds_ref[2][1]  = 0.0;	
	}

    return SS_ref_db;
}

/**
   retrieve reference thermodynamic data for mb_ilm
*/
SS_ref G_SS_mb_ilm_function(SS_ref SS_ref_db, int EM_database, int len_ox, bulk_info z_b, double eps){
    
    int i, j;
    int n_em = SS_ref_db.n_em;
    
    char   *EM_tmp[] 		= {"oilm","dilm","dhem"};
    for (int i = 0; i < SS_ref_db.n_em; i++){
        strcpy(SS_ref_db.EM_list[i],EM_tmp[i]);
    };
    int n_xeos = SS_ref_db.n_xeos;
    char   *CV_tmp[] 		= {"x","Q"};
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
        strcpy(SS_ref_db.CV_list[i],CV_tmp[i]);
    };
    SS_ref_db.W[0] = 15.6;
    SS_ref_db.W[1] = 26.6;
    SS_ref_db.W[2] = 11.0;
    
    
    em_data ilm_di 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"ilm", 
    										"disordered"	);
    
    em_data hem_di 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"hem", 
    										"disordered"	);
    
    SS_ref_db.gbase[0] 		= 0.009426*z_b.T + ilm_di.gb - 13.6075;
    SS_ref_db.gbase[1] 		= -0.0021*z_b.T + ilm_di.gb + 1.9928;
    SS_ref_db.gbase[2] 		= hem_di.gb;
    
    SS_ref_db.ElShearMod[0] 	= ilm_di.ElShearMod;
    SS_ref_db.ElShearMod[1] 	= ilm_di.ElShearMod;
    SS_ref_db.ElShearMod[2] 	= hem_di.ElShearMod;
    
    for (i = 0; i < len_ox; i++){
        SS_ref_db.Comp[0][i] 	= ilm_di.C[i];
        SS_ref_db.Comp[1][i] 	= ilm_di.C[i];
        SS_ref_db.Comp[2][i] 	= hem_di.C[i];
    }
    
    for (i = 0; i < n_em; i++){
        SS_ref_db.z_em[i] = 1.0;
    };
    
    SS_ref_db.bounds_ref[0][0] = 0.0+eps;  SS_ref_db.bounds_ref[0][1] = 1.0-eps;
    SS_ref_db.bounds_ref[1][0] = -0.99+eps;  SS_ref_db.bounds_ref[1][1] = 0.99-eps;
    
	if (z_b.bulk_rock[8] == 0.){ 	    //O				
		SS_ref_db.z_em[2]          = 0.0;
        SS_ref_db.d_em[2]          = 1.0;
		SS_ref_db.bounds_ref[0][0] = 1.0; 
		SS_ref_db.bounds_ref[0][1] = 1.0;	
	}
	if (z_b.bulk_rock[7] == 0.){ 		//TiO2	
		SS_ref_db.z_em[0]          = 0.0;
        SS_ref_db.d_em[0]          = 1.0;
		SS_ref_db.bounds_ref[0][0] = 0.0; 
		SS_ref_db.bounds_ref[0][1] = 0.0;	

		SS_ref_db.z_em[1]          = 0.0;
        SS_ref_db.d_em[1]          = 1.0;
		SS_ref_db.bounds_ref[1][0] = 0.0; 
		SS_ref_db.bounds_ref[1][1] = 0.0;	
	}

    /* this lists the index of the order variables */
    SS_ref_db.orderVar      = 1;
    SS_ref_db.idOrderVar[1] = -1.0;
    
    return SS_ref_db;
}

/**
   retrieve reference thermodynamic data for mb_ilmm
*/
SS_ref G_SS_mb_ilmm_function(SS_ref SS_ref_db, int EM_database, int len_ox, bulk_info z_b, double eps){
    
    int i, j;
    int n_em = SS_ref_db.n_em;
    
    char   *EM_tmp[] 		= {"oilm","dilm","dhem","geik"};
    for (int i = 0; i < SS_ref_db.n_em; i++){
        strcpy(SS_ref_db.EM_list[i],EM_tmp[i]);
    };
    int n_xeos = SS_ref_db.n_xeos;
    char   *CV_tmp[] 		= {"c","t","Q"};
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
        strcpy(SS_ref_db.CV_list[i],CV_tmp[i]);
    };
    SS_ref_db.W[0] = 15.6;
    SS_ref_db.W[1] = 26.6;
    SS_ref_db.W[2] = 4.00;
    SS_ref_db.W[3] = 11.0;
    SS_ref_db.W[4] = 4.00;
    SS_ref_db.W[5] = 36.0;
    
    
    em_data ilm_di 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"ilm", 
    										"disordered"	);
    
    em_data hem_di 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"hem", 
    										"disordered"	);
    
    em_data geik_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"geik", 
    										"equilibrium"	);
    
    SS_ref_db.gbase[0] 		= 0.009426*z_b.T + ilm_di.gb - 13.6075;
    SS_ref_db.gbase[1] 		= -0.0021*z_b.T + ilm_di.gb + 1.9928;
    SS_ref_db.gbase[2] 		= hem_di.gb;
    SS_ref_db.gbase[3] 		= geik_eq.gb;
    
    SS_ref_db.ElShearMod[0] 	= ilm_di.ElShearMod;
    SS_ref_db.ElShearMod[1] 	= ilm_di.ElShearMod;
    SS_ref_db.ElShearMod[2] 	= hem_di.ElShearMod;
    SS_ref_db.ElShearMod[3] 	= geik_eq.ElShearMod;
    
    for (i = 0; i < len_ox; i++){
        SS_ref_db.Comp[0][i] 	= ilm_di.C[i];
        SS_ref_db.Comp[1][i] 	= ilm_di.C[i];
        SS_ref_db.Comp[2][i] 	= hem_di.C[i];
        SS_ref_db.Comp[3][i] 	= geik_eq.C[i];
    }
    
    for (i = 0; i < n_em; i++){
        SS_ref_db.z_em[i] = 1.0;
    };
    
    SS_ref_db.bounds_ref[0][0] = 0.0+eps;  SS_ref_db.bounds_ref[0][1] = 1.0-eps;
    SS_ref_db.bounds_ref[1][0] = 0.0+eps;  SS_ref_db.bounds_ref[1][1] = 1.0-eps;
    SS_ref_db.bounds_ref[2][0] = 0.0+eps;  SS_ref_db.bounds_ref[2][1] = 1.0-eps;
    


    return SS_ref_db;
}

/**
   retrieve reference thermodynamic data for mb_ep
*/
SS_ref G_SS_mb_ep_function(SS_ref SS_ref_db, int EM_database, int len_ox, bulk_info z_b, double eps){
    
    int i, j;
    int n_em = SS_ref_db.n_em;
    
    char   *EM_tmp[] 		= {"cz","ep","fep"};
    for (int i = 0; i < SS_ref_db.n_em; i++){
        strcpy(SS_ref_db.EM_list[i],EM_tmp[i]);
    };
    int n_xeos = SS_ref_db.n_xeos;
    char   *CV_tmp[] 		= {"f","Q"};
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
        strcpy(SS_ref_db.CV_list[i],CV_tmp[i]);
    };
    SS_ref_db.W[0] = 1.00;
    SS_ref_db.W[1] = 3.00;
    SS_ref_db.W[2] = 1.00;
    
    
    em_data cz_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"cz", 
    										"equilibrium"	);
    
    em_data ep_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"ep", 
    										"equilibrium"	);
    
    em_data fep_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"fep", 
    										"equilibrium"	);
    
    SS_ref_db.gbase[0] 		= cz_eq.gb;
    SS_ref_db.gbase[1] 		= ep_eq.gb;
    SS_ref_db.gbase[2] 		= fep_eq.gb;
    
    SS_ref_db.ElShearMod[0] 	= cz_eq.ElShearMod;
    SS_ref_db.ElShearMod[1] 	= ep_eq.ElShearMod;
    SS_ref_db.ElShearMod[2] 	= fep_eq.ElShearMod;
    
    for (i = 0; i < len_ox; i++){
        SS_ref_db.Comp[0][i] 	= cz_eq.C[i];
        SS_ref_db.Comp[1][i] 	= ep_eq.C[i];
        SS_ref_db.Comp[2][i] 	= fep_eq.C[i];
    }
    
    for (i = 0; i < n_em; i++){
        SS_ref_db.z_em[i] = 1.0;
    };
    
    SS_ref_db.bounds_ref[0][0] = 0.0+eps;  SS_ref_db.bounds_ref[0][1] = 1.0-eps;
    SS_ref_db.bounds_ref[1][0] = 0.0+eps;  SS_ref_db.bounds_ref[1][1] = 0.5-eps;

	if (z_b.bulk_rock[8] == 0.){ 	    //O			
		SS_ref_db.z_em[1]          = 0.0;
        SS_ref_db.d_em[1]          = 1.0;
		SS_ref_db.z_em[2]          = 0.0;
        SS_ref_db.d_em[2]          = 1.0;
		SS_ref_db.bounds_ref[0][0] = 0.0; 
		SS_ref_db.bounds_ref[0][1] = 0.0;	
		SS_ref_db.bounds_ref[1][0] = 0.0; 
		SS_ref_db.bounds_ref[1][1] = 0.0;	
	}

    return SS_ref_db;
}

/**
   retrieve reference thermodynamic data for mb_bi
*/
SS_ref G_SS_mb_bi_function(SS_ref SS_ref_db, int EM_database, int len_ox, bulk_info z_b, double eps){
    
    int i, j;
    int n_em = SS_ref_db.n_em;
    
    char   *EM_tmp[] 		= {"phl","annm","obi","east","tbi","fbi"};
    for (int i = 0; i < SS_ref_db.n_em; i++){
        strcpy(SS_ref_db.EM_list[i],EM_tmp[i]);
    };
    int n_xeos = SS_ref_db.n_xeos;
    char   *CV_tmp[] 		= {"x","y","f","t","Q"};
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
        strcpy(SS_ref_db.CV_list[i],CV_tmp[i]);
    };
    SS_ref_db.W[0] = 12.0;
    SS_ref_db.W[1] = 4.00;
    SS_ref_db.W[2] = 10.0;
    SS_ref_db.W[3] = 30.0;
    SS_ref_db.W[4] = 8.00;
    SS_ref_db.W[5] = 8.00;
    SS_ref_db.W[6] = 15.0;
    SS_ref_db.W[7] = 32.0;
    SS_ref_db.W[8] = 13.6;
    SS_ref_db.W[9] = 7.00;
    SS_ref_db.W[10] = 24.0;
    SS_ref_db.W[11] = 5.60;
    SS_ref_db.W[12] = 40.0;
    SS_ref_db.W[13] = 1.00;
    SS_ref_db.W[14] = 40.0;
    
    
    em_data phl_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"phl", 
    										"equilibrium"	);
    
    em_data ann_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"ann", 
    										"equilibrium"	);
    
    em_data east_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"east", 
    										"equilibrium"	);
    
    em_data br_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"br", 
    										"equilibrium"	);
    
    em_data ru_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"ru", 
    										"equilibrium"	);
    
    em_data gr_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"gr", 
    										"equilibrium"	);
    
    em_data andr_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"andr", 
    										"equilibrium"	);
    
    SS_ref_db.gbase[0] 		= phl_eq.gb;
    SS_ref_db.gbase[1] 		= ann_eq.gb - 3.0;
    SS_ref_db.gbase[2] 		= 1.0/3.0*ann_eq.gb + 2.0/3.0*phl_eq.gb - 3.0;
    SS_ref_db.gbase[3] 		= east_eq.gb;
    SS_ref_db.gbase[4] 		= -br_eq.gb + phl_eq.gb + ru_eq.gb + 55.0;
    SS_ref_db.gbase[5] 		= 0.5*andr_eq.gb + east_eq.gb - 0.5*gr_eq.gb - 3.0;
    
    SS_ref_db.ElShearMod[0] 	= phl_eq.ElShearMod;
    SS_ref_db.ElShearMod[1] 	= ann_eq.ElShearMod;
    SS_ref_db.ElShearMod[2] 	= 1.0/3.0*ann_eq.ElShearMod + 2.0/3.0*phl_eq.ElShearMod;
    SS_ref_db.ElShearMod[3] 	= east_eq.ElShearMod;
    SS_ref_db.ElShearMod[4] 	= -br_eq.ElShearMod + phl_eq.ElShearMod + ru_eq.ElShearMod;
    SS_ref_db.ElShearMod[5] 	= 0.5*andr_eq.ElShearMod + east_eq.ElShearMod - 0.5*gr_eq.ElShearMod;
    
    for (i = 0; i < len_ox; i++){
        SS_ref_db.Comp[0][i] 	= phl_eq.C[i];
        SS_ref_db.Comp[1][i] 	= ann_eq.C[i];
        SS_ref_db.Comp[2][i] 	= 1.0/3.0*ann_eq.C[i] + 2.0/3.0*phl_eq.C[i];
        SS_ref_db.Comp[3][i] 	= east_eq.C[i];
        SS_ref_db.Comp[4][i] 	= -br_eq.C[i] + phl_eq.C[i] + ru_eq.C[i];
        SS_ref_db.Comp[5][i] 	= 0.5*andr_eq.C[i] + east_eq.C[i] - 0.5*gr_eq.C[i];
    }
    
    for (i = 0; i < n_em; i++){
        SS_ref_db.z_em[i] = 1.0;
    };
    
    SS_ref_db.bounds_ref[0][0] = 0.0+eps;  SS_ref_db.bounds_ref[0][1] = 1.0-eps;
    SS_ref_db.bounds_ref[1][0] = 0.0+eps;  SS_ref_db.bounds_ref[1][1] = 1.0-eps;
    SS_ref_db.bounds_ref[2][0] = 0.0+eps;  SS_ref_db.bounds_ref[2][1] = 1.0-eps;
    SS_ref_db.bounds_ref[3][0] = 0.0+eps;  SS_ref_db.bounds_ref[3][1] = 1.0-eps;
    SS_ref_db.bounds_ref[4][0] = 0.0+eps;  SS_ref_db.bounds_ref[4][1] = 1.0-eps;

	if (z_b.bulk_rock[8] == 0.){ 	    //O				
		SS_ref_db.z_em[5]          = 0.0;
        SS_ref_db.d_em[5]          = 1.0;
		SS_ref_db.bounds_ref[2][0] = 0.0; 
		SS_ref_db.bounds_ref[2][1] = 0.0;	
	}
	if (z_b.bulk_rock[7] == 0.){ 		//TiO2	
		SS_ref_db.z_em[4]          = 0.0;
        SS_ref_db.d_em[4]          = 1.0;
		SS_ref_db.bounds_ref[3][0]  = 0.0; 
		SS_ref_db.bounds_ref[3][1]  = 0.0;	
	}
    return SS_ref_db;
}

/**
   retrieve reference thermodynamic data for mb_mu
*/
SS_ref G_SS_mb_mu_function(SS_ref SS_ref_db, int EM_database, int len_ox, bulk_info z_b, double eps){
    
    int i, j;
    int n_em = SS_ref_db.n_em;
    
    char   *EM_tmp[] 		= {"mu","cel","fcel","pa","mam","fmu"};
    for (int i = 0; i < SS_ref_db.n_em; i++){
        strcpy(SS_ref_db.EM_list[i],EM_tmp[i]);
    };
    int n_xeos = SS_ref_db.n_xeos;
    char   *CV_tmp[] 		= {"x","y","f","n","c"};
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
        strcpy(SS_ref_db.CV_list[i],CV_tmp[i]);
    };
    SS_ref_db.W[0] = 0.2*SS_ref_db.P;
    SS_ref_db.W[1] = 0.2*SS_ref_db.P;
    SS_ref_db.W[2] = 0.353*SS_ref_db.P + 0.0034*SS_ref_db.T + 10.12;
    SS_ref_db.W[3] = 35.0;
    SS_ref_db.W[4] = 0.0;
    SS_ref_db.W[5] = 0.0;
    SS_ref_db.W[6] = 0.25*SS_ref_db.P + 45.0;
    SS_ref_db.W[7] = 50.0;
    SS_ref_db.W[8] = 0.0;
    SS_ref_db.W[9] = 0.25*SS_ref_db.P + 45.0;
    SS_ref_db.W[10] = 50.0;
    SS_ref_db.W[11] = 0.0;
    SS_ref_db.W[12] = 15.0;
    SS_ref_db.W[13] = 30.0;
    SS_ref_db.W[14] = 35.0;
    
    SS_ref_db.v[0] = 0.630;
    SS_ref_db.v[1] = 0.630;
    SS_ref_db.v[2] = 0.630;
    SS_ref_db.v[3] = 0.370;
    SS_ref_db.v[4] = 0.630;
    SS_ref_db.v[5] = 0.630;
    
    
    em_data mu_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"mu", 
    										"equilibrium"	);
    
    em_data cel_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"cel", 
    										"equilibrium"	);
    
    em_data fcel_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"fcel", 
    										"equilibrium"	);
    
    em_data pa_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"pa", 
    										"equilibrium"	);
    
    em_data ma_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"ma", 
    										"equilibrium"	);
    
    em_data gr_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"gr", 
    										"equilibrium"	);
    
    em_data andr_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"andr", 
    										"equilibrium"	);
    
    SS_ref_db.gbase[0] 		= mu_eq.gb;
    SS_ref_db.gbase[1] 		= cel_eq.gb;
    SS_ref_db.gbase[2] 		= fcel_eq.gb;
    SS_ref_db.gbase[3] 		= pa_eq.gb;
    SS_ref_db.gbase[4] 		= ma_eq.gb + 5.0;
    SS_ref_db.gbase[5] 		= 0.5*andr_eq.gb - 0.5*gr_eq.gb + mu_eq.gb + 25.0;
    
    SS_ref_db.ElShearMod[0] 	= mu_eq.ElShearMod;
    SS_ref_db.ElShearMod[1] 	= cel_eq.ElShearMod;
    SS_ref_db.ElShearMod[2] 	= fcel_eq.ElShearMod;
    SS_ref_db.ElShearMod[3] 	= pa_eq.ElShearMod;
    SS_ref_db.ElShearMod[4] 	= ma_eq.ElShearMod;
    SS_ref_db.ElShearMod[5] 	= 0.5*andr_eq.ElShearMod - 0.5*gr_eq.ElShearMod + mu_eq.ElShearMod;
    
    for (i = 0; i < len_ox; i++){
        SS_ref_db.Comp[0][i] 	= mu_eq.C[i];
        SS_ref_db.Comp[1][i] 	= cel_eq.C[i];
        SS_ref_db.Comp[2][i] 	= fcel_eq.C[i];
        SS_ref_db.Comp[3][i] 	= pa_eq.C[i];
        SS_ref_db.Comp[4][i] 	= ma_eq.C[i];
        SS_ref_db.Comp[5][i] 	= 0.5*andr_eq.C[i] - 0.5*gr_eq.C[i] + mu_eq.C[i];
    }
    
    for (i = 0; i < n_em; i++){
        SS_ref_db.z_em[i] = 1.0;
    };
    
    SS_ref_db.bounds_ref[0][0] = 0.0+eps;  SS_ref_db.bounds_ref[0][1] = 1.0-eps;
    SS_ref_db.bounds_ref[1][0] = 0.0+eps;  SS_ref_db.bounds_ref[1][1] = 1.0-eps;
    SS_ref_db.bounds_ref[2][0] = 0.0+eps;  SS_ref_db.bounds_ref[2][1] = 1.0-eps;
    SS_ref_db.bounds_ref[3][0] = 0.0+eps;  SS_ref_db.bounds_ref[3][1] = 1.0-eps;
    SS_ref_db.bounds_ref[4][0] = 0.0+eps;  SS_ref_db.bounds_ref[4][1] = 1.0-eps;

	if (z_b.bulk_rock[8] == 0.){ 	    //O				
		SS_ref_db.z_em[5]          = 0.0;
        SS_ref_db.d_em[5]          = 1.0;
		SS_ref_db.bounds_ref[2][0] = 0.0; 
		SS_ref_db.bounds_ref[2][1] = 0.0;	
	}

    return SS_ref_db;
}

/**
   retrieve reference thermodynamic data for mb_chl
*/
SS_ref G_SS_mb_chl_function(SS_ref SS_ref_db, int EM_database, int len_ox, bulk_info z_b, double eps){
    
    int i, j;
    int n_em = SS_ref_db.n_em;
    
    char   *EM_tmp[] 		= {"clin","afchl","ames","daph","ochl1","ochl4","f3clin"};
    for (int i = 0; i < SS_ref_db.n_em; i++){
        strcpy(SS_ref_db.EM_list[i],EM_tmp[i]);
    };
    int n_xeos = SS_ref_db.n_xeos;
    char   *CV_tmp[] 		= {"x","y","f","QAl","Q1","Q4"};
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
        strcpy(SS_ref_db.CV_list[i],CV_tmp[i]);
    };
    SS_ref_db.W[0] = 17.0;
    SS_ref_db.W[1] = 17.0;
    SS_ref_db.W[2] = 20.0;
    SS_ref_db.W[3] = 30.0;
    SS_ref_db.W[4] = 21.0;
    SS_ref_db.W[5] = 2.00;
    SS_ref_db.W[6] = 16.0;
    SS_ref_db.W[7] = 37.0;
    SS_ref_db.W[8] = 20.0;
    SS_ref_db.W[9] = 4.00;
    SS_ref_db.W[10] = 15.0;
    SS_ref_db.W[11] = 30.0;
    SS_ref_db.W[12] = 29.0;
    SS_ref_db.W[13] = 13.0;
    SS_ref_db.W[14] = 19.0;
    SS_ref_db.W[15] = 18.0;
    SS_ref_db.W[16] = 33.0;
    SS_ref_db.W[17] = 22.0;
    SS_ref_db.W[18] = 24.0;
    SS_ref_db.W[19] = 28.6;
    SS_ref_db.W[20] = 19.0;
    
    
    em_data clin_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"clin", 
    										"equilibrium"	);
    
    em_data afchl_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"afchl", 
    										"equilibrium"	);
    
    em_data ames_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"ames", 
    										"equilibrium"	);
    
    em_data daph_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"daph", 
    										"equilibrium"	);
    
    em_data gr_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"gr", 
    										"equilibrium"	);
    
    em_data andr_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"andr", 
    										"equilibrium"	);
    
    SS_ref_db.gbase[0] 		= clin_eq.gb;
    SS_ref_db.gbase[1] 		= afchl_eq.gb;
    SS_ref_db.gbase[2] 		= ames_eq.gb;
    SS_ref_db.gbase[3] 		= daph_eq.gb;
    SS_ref_db.gbase[4] 		= afchl_eq.gb -clin_eq.gb + daph_eq.gb + 3.0;
    SS_ref_db.gbase[5] 		= afchl_eq.gb - 0.2*clin_eq.gb + 0.2*daph_eq.gb + 2.4;
    SS_ref_db.gbase[6] 		= 0.5*andr_eq.gb + clin_eq.gb - 0.5*gr_eq.gb + 2.0;
    
    SS_ref_db.ElShearMod[0] 	= clin_eq.ElShearMod;
    SS_ref_db.ElShearMod[1] 	= afchl_eq.ElShearMod;
    SS_ref_db.ElShearMod[2] 	= ames_eq.ElShearMod;
    SS_ref_db.ElShearMod[3] 	= daph_eq.ElShearMod;
    SS_ref_db.ElShearMod[4] 	= afchl_eq.ElShearMod -clin_eq.ElShearMod + daph_eq.ElShearMod;
    SS_ref_db.ElShearMod[5] 	= afchl_eq.ElShearMod - 0.2*clin_eq.ElShearMod + 0.2*daph_eq.ElShearMod;
    SS_ref_db.ElShearMod[6] 	= 0.5*andr_eq.ElShearMod + clin_eq.ElShearMod - 0.5*gr_eq.ElShearMod;
    
    for (i = 0; i < len_ox; i++){
        SS_ref_db.Comp[0][i] 	= clin_eq.C[i];
        SS_ref_db.Comp[1][i] 	= afchl_eq.C[i];
        SS_ref_db.Comp[2][i] 	= ames_eq.C[i];
        SS_ref_db.Comp[3][i] 	= daph_eq.C[i];
        SS_ref_db.Comp[4][i] 	= afchl_eq.C[i] -clin_eq.C[i] + daph_eq.C[i];
        SS_ref_db.Comp[5][i] 	= afchl_eq.C[i] - 0.2*clin_eq.C[i] + 0.2*daph_eq.C[i];
        SS_ref_db.Comp[6][i] 	= 0.5*andr_eq.C[i] + clin_eq.C[i] - 0.5*gr_eq.C[i];
    }
    
    for (i = 0; i < n_em; i++){
        SS_ref_db.z_em[i] = 1.0;
    };
    
    SS_ref_db.bounds_ref[0][0] = 0.0+eps;  SS_ref_db.bounds_ref[0][1] = 1.0-eps;
    SS_ref_db.bounds_ref[1][0] = 0.0+eps;  SS_ref_db.bounds_ref[1][1] = 1.0-eps;
    SS_ref_db.bounds_ref[2][0] = 0.0+eps;  SS_ref_db.bounds_ref[2][1] = 1.0-eps;
    SS_ref_db.bounds_ref[3][0] = -1.0+eps;  SS_ref_db.bounds_ref[3][1] = 1.0-eps;
    SS_ref_db.bounds_ref[4][0] = -1.0+eps;  SS_ref_db.bounds_ref[4][1] = 1.0-eps;
    SS_ref_db.bounds_ref[5][0] = -1.0+eps;  SS_ref_db.bounds_ref[5][1] = 1.0-eps;

	if (z_b.bulk_rock[8] == 0.){ 	    //O				
		SS_ref_db.z_em[6]          = 0.0;
        SS_ref_db.d_em[6]          = 1.0;
		SS_ref_db.bounds_ref[2][0] = 0.0; 
		SS_ref_db.bounds_ref[2][1] = 0.0;	
	}
    return SS_ref_db;
}



/**************************************************************************************/
/**************************************************************************************/
/************************METAPELITE DATABASE (White al., 2011)*************************/
/**************************************************************************************/
/**************************************************************************************/

/**
   retrieve reference thermodynamic data for mp_liq
*/
SS_ref G_SS_mp_liq_function(SS_ref SS_ref_db, int EM_database, int len_ox, bulk_info z_b, double eps){
    
    int i, j;
    int n_em = SS_ref_db.n_em;
    
    char   *EM_tmp[] 		= {"q4L","abL","kspL","anL","slL","fo2L","fa2L","h2oL"};
    for (int i = 0; i < SS_ref_db.n_em; i++){
        strcpy(SS_ref_db.EM_list[i],EM_tmp[i]);
    };
    int n_xeos = SS_ref_db.n_xeos;
    char   *CV_tmp[] 		= {"q","fsp","na","an","ol","x","h2o"};
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
        strcpy(SS_ref_db.CV_list[i],CV_tmp[i]);
    };
    SS_ref_db.W[0] = 12.0 - 0.4*SS_ref_db.P;
    SS_ref_db.W[1] = -0.5*SS_ref_db.P - 2.0;
    SS_ref_db.W[2] = 5.0;
    SS_ref_db.W[3] = 12.0;
    SS_ref_db.W[4] = 12.0 - 0.4*SS_ref_db.P;
    SS_ref_db.W[5] = 14.0;
    SS_ref_db.W[6] = 17.0 - 0.5*SS_ref_db.P;
    SS_ref_db.W[7] = 3.0*SS_ref_db.P - 6.0;
    SS_ref_db.W[8] = 0.0;
    SS_ref_db.W[9] = 12.0;
    SS_ref_db.W[10] = 10.0;
    SS_ref_db.W[11] = 2.0;
    SS_ref_db.W[12] = -0.3*SS_ref_db.P - 1.5;
    SS_ref_db.W[13] = -SS_ref_db.P;
    SS_ref_db.W[14] = 12.0;
    SS_ref_db.W[15] = 12.0;
    SS_ref_db.W[16] = 12.0;
    SS_ref_db.W[17] = 9.5 - 0.3*SS_ref_db.P;
    SS_ref_db.W[18] = 0.0;
    SS_ref_db.W[19] = 0.0;
    SS_ref_db.W[20] = 0.0;
    SS_ref_db.W[21] = 7.5 - 0.5*SS_ref_db.P;
    SS_ref_db.W[22] = 12.0;
    SS_ref_db.W[23] = 12.0;
    SS_ref_db.W[24] = 11.0;
    SS_ref_db.W[25] = 18.0;
    SS_ref_db.W[26] = 11.0- 0.5*SS_ref_db.P;
    SS_ref_db.W[27] = 12.0;
    
    
    em_data qL_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"qL", 
    										"equilibrium"	);
    
    em_data abL_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"abL", 
    										"equilibrium"	);
    
    em_data kspL_eq 	= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"kspL", 
    										"equilibrium"	);
    
    em_data anL_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"anL", 
    										"equilibrium"	);
    
    em_data silL_eq 	= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"silL", 
    										"equilibrium"	);
    
    em_data foL_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"foL", 
    										"equilibrium"	);
    
    em_data faL_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"faL", 
    										"equilibrium"	);
    
    em_data h2oL_eq 	= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"h2oL", 
    										"equilibrium"	);
    
    SS_ref_db.gbase[0] 		= 4.0*qL_eq.gb;
    SS_ref_db.gbase[1] 		= abL_eq.gb;
    SS_ref_db.gbase[2] 		= kspL_eq.gb;
    SS_ref_db.gbase[3] 		= anL_eq.gb;
    SS_ref_db.gbase[4] 		= 1.6*silL_eq.gb - 23.0;
    SS_ref_db.gbase[5] 		= 2.0*foL_eq.gb - 10.0;
    SS_ref_db.gbase[6] 		= -1.3*z_b.P + 2.0*faL_eq.gb - 9.0;
    SS_ref_db.gbase[7] 		= h2oL_eq.gb;
    
    SS_ref_db.ElShearMod[0] 	= 4.0*qL_eq.ElShearMod;
    SS_ref_db.ElShearMod[1] 	= abL_eq.ElShearMod;
    SS_ref_db.ElShearMod[2] 	= kspL_eq.ElShearMod;
    SS_ref_db.ElShearMod[3] 	= anL_eq.ElShearMod;
    SS_ref_db.ElShearMod[4] 	= 1.6*silL_eq.ElShearMod;
    SS_ref_db.ElShearMod[5] 	= 2.0*foL_eq.ElShearMod;
    SS_ref_db.ElShearMod[6] 	= 2.0*faL_eq.ElShearMod;
    SS_ref_db.ElShearMod[7] 	= h2oL_eq.ElShearMod;
    
    for (i = 0; i < len_ox; i++){
        SS_ref_db.Comp[0][i] 	= 4.0*qL_eq.C[i];
        SS_ref_db.Comp[1][i] 	= abL_eq.C[i];
        SS_ref_db.Comp[2][i] 	= kspL_eq.C[i];
        SS_ref_db.Comp[3][i] 	= anL_eq.C[i];
        SS_ref_db.Comp[4][i] 	= 1.6*silL_eq.C[i];
        SS_ref_db.Comp[5][i] 	= 2.0*foL_eq.C[i];
        SS_ref_db.Comp[6][i] 	= 2.0*faL_eq.C[i];
        SS_ref_db.Comp[7][i] 	= h2oL_eq.C[i];
    }
    
    for (i = 0; i < n_em; i++){
        SS_ref_db.z_em[i] = 1.0;
    };
    
    SS_ref_db.bounds_ref[0][0] = 0.0+eps;  SS_ref_db.bounds_ref[0][1] = 1.0-eps;
    SS_ref_db.bounds_ref[1][0] = 0.0+eps;  SS_ref_db.bounds_ref[1][1] = 1.0-eps;
    SS_ref_db.bounds_ref[2][0] = 0.0+eps;  SS_ref_db.bounds_ref[2][1] = 1.0-eps;
    SS_ref_db.bounds_ref[3][0] = 0.0+eps;  SS_ref_db.bounds_ref[3][1] = 1.0-eps;
    SS_ref_db.bounds_ref[4][0] = 0.0+eps;  SS_ref_db.bounds_ref[4][1] = 1.0-eps;
    SS_ref_db.bounds_ref[5][0] = 0.0+eps;  SS_ref_db.bounds_ref[5][1] = 1.0-eps;
    SS_ref_db.bounds_ref[6][0] = 0.0+eps;  SS_ref_db.bounds_ref[6][1] = 1.0-eps;
    
	if (z_b.bulk_rock[10] == 0.){ 					
		SS_ref_db.z_em[7]          = 0.0;
		SS_ref_db.bounds_ref[6][0] = eps; 
		SS_ref_db.bounds_ref[6][1] = eps;	
	}

    return SS_ref_db;
}

/**
   retrieve reference thermodynamic data for mp_bi
*/
SS_ref G_SS_mp_bi_function(SS_ref SS_ref_db, int EM_database, int len_ox, bulk_info z_b, double eps){
    
    int i, j;
    int n_em = SS_ref_db.n_em;
    
    char   *EM_tmp[] 		= {"phl","annm","obi","east","tbi","fbi","mmbi"};
    for (int i = 0; i < SS_ref_db.n_em; i++){
        strcpy(SS_ref_db.EM_list[i],EM_tmp[i]);
    };
    int n_xeos = SS_ref_db.n_xeos;
    char   *CV_tmp[] 		= {"x","m","y","f","t","Q"};
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
        strcpy(SS_ref_db.CV_list[i],CV_tmp[i]);
    };
    SS_ref_db.W[0] = 12.0;
    SS_ref_db.W[1] = 4.0;
    SS_ref_db.W[2] = 10.0;
    SS_ref_db.W[3] = 30.0;
    SS_ref_db.W[4] = 8.0;
    SS_ref_db.W[5] = 9.0;
    SS_ref_db.W[6] = 8.0;
    SS_ref_db.W[7] = 15.0;
    SS_ref_db.W[8] = 32.0;
    SS_ref_db.W[9] = 13.6;
    SS_ref_db.W[10] = 6.3;
    SS_ref_db.W[11] = 7.0;
    SS_ref_db.W[12] = 24.0;
    SS_ref_db.W[13] = 5.6;
    SS_ref_db.W[14] = 8.1;
    SS_ref_db.W[15] = 40.0;
    SS_ref_db.W[16] = 1.0;
    SS_ref_db.W[17] = 13.0;
    SS_ref_db.W[18] = 40.0;
    SS_ref_db.W[19] = 30.0;
    SS_ref_db.W[20] = 11.6;
    
    
    em_data phl_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"phl", 
    										"equilibrium"	);
    
    em_data ann_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"ann", 
    										"equilibrium"	);
    
    em_data east_eq 		= get_em_data(	EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"east", 
    										"equilibrium"	);
    
    em_data br_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"br", 
    										"equilibrium"	);
    
    em_data ru_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"ru", 
    										"equilibrium"	);
    
    em_data andr_eq 		= get_em_data(	EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"andr", 
    										"equilibrium"	);
    
    em_data gr_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"gr", 
    										"equilibrium"	);
    
    em_data mnbi_eq 		= get_em_data(	EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"mnbi", 
    										"equilibrium"	);
    
    SS_ref_db.gbase[0] 		= phl_eq.gb;
    SS_ref_db.gbase[1] 		= ann_eq.gb - 3.0;
    SS_ref_db.gbase[2] 		= 1.0/3.0*ann_eq.gb + 2.0/3.0*phl_eq.gb - 3.0;
    SS_ref_db.gbase[3] 		= east_eq.gb;
    SS_ref_db.gbase[4] 		= -br_eq.gb + phl_eq.gb + ru_eq.gb + 55.0;
    SS_ref_db.gbase[5] 		= 0.5*andr_eq.gb + east_eq.gb - 0.5*gr_eq.gb - 3.0;
    SS_ref_db.gbase[6] 		= mnbi_eq.gb - 7.89;
    
    SS_ref_db.ElShearMod[0] 	= phl_eq.ElShearMod;
    SS_ref_db.ElShearMod[1] 	= ann_eq.ElShearMod;
    SS_ref_db.ElShearMod[2] 	= 1.0/3.0*ann_eq.ElShearMod + 2.0/3.0*phl_eq.ElShearMod;
    SS_ref_db.ElShearMod[3] 	= east_eq.ElShearMod;
    SS_ref_db.ElShearMod[4] 	= -br_eq.ElShearMod + phl_eq.ElShearMod + ru_eq.ElShearMod;
    SS_ref_db.ElShearMod[5] 	= 0.5*andr_eq.ElShearMod + east_eq.ElShearMod - 0.5*gr_eq.ElShearMod;
    SS_ref_db.ElShearMod[6] 	= mnbi_eq.ElShearMod;
    
    for (i = 0; i < len_ox; i++){
        SS_ref_db.Comp[0][i] 	= phl_eq.C[i];
        SS_ref_db.Comp[1][i] 	= ann_eq.C[i];
        SS_ref_db.Comp[2][i] 	= 1.0/3.0*ann_eq.C[i] + 2.0/3.0*phl_eq.C[i];
        SS_ref_db.Comp[3][i] 	= east_eq.C[i];
        SS_ref_db.Comp[4][i] 	= -br_eq.C[i] + phl_eq.C[i] + ru_eq.C[i];
        SS_ref_db.Comp[5][i] 	= 0.5*andr_eq.C[i] + east_eq.C[i] - 0.5*gr_eq.C[i];
        SS_ref_db.Comp[6][i] 	= mnbi_eq.C[i];
    }
    
    for (i = 0; i < n_em; i++){
        SS_ref_db.z_em[i] = 1.0;
    };
    
    SS_ref_db.bounds_ref[0][0] = 0.0+eps;  SS_ref_db.bounds_ref[0][1] = 1.0-eps;
    SS_ref_db.bounds_ref[1][0] = 0.0+eps;  SS_ref_db.bounds_ref[1][1] = 1.0-eps;
    SS_ref_db.bounds_ref[2][0] = 0.0+eps;  SS_ref_db.bounds_ref[2][1] = 1.0-eps;
    SS_ref_db.bounds_ref[3][0] = 0.0+eps;  SS_ref_db.bounds_ref[3][1] = 1.0-eps;
    SS_ref_db.bounds_ref[4][0] = 0.0+eps;  SS_ref_db.bounds_ref[4][1] = 1.0-eps;
    SS_ref_db.bounds_ref[5][0] = 0.0+eps;  SS_ref_db.bounds_ref[5][1] = 1.0-eps;
    
	if (z_b.bulk_rock[8] == 0.){ 					
		SS_ref_db.z_em[5]          = 0.0;
        SS_ref_db.d_em[5]          = 1.0;
		SS_ref_db.bounds_ref[3][0] = 0.0; 
		SS_ref_db.bounds_ref[3][1] = 0.0;	
	}
    return SS_ref_db;
}

/**
   retrieve reference thermodynamic data for mp_cd
*/
SS_ref G_SS_mp_cd_function(SS_ref SS_ref_db, int EM_database, int len_ox, bulk_info z_b, double eps){
    
    int i, j;
    int n_em = SS_ref_db.n_em;
    
    char   *EM_tmp[] 		= {"crd","fcrd","hcrd","mncd"};
    for (int i = 0; i < SS_ref_db.n_em; i++){
        strcpy(SS_ref_db.EM_list[i],EM_tmp[i]);
    };
    int n_xeos = SS_ref_db.n_xeos;
    char   *CV_tmp[] 		= {"x","m","h"};
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
        strcpy(SS_ref_db.CV_list[i],CV_tmp[i]);
    };
    SS_ref_db.W[0] = 8.0;
    SS_ref_db.W[1] = 0.0;
    SS_ref_db.W[2] = 6.0;
    SS_ref_db.W[3] = 9.0;
    SS_ref_db.W[4] = 4.0;
    SS_ref_db.W[5] = 6.0;
    
    
    em_data crd_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"crd", 
    										"equilibrium"	);
    
    em_data fcrd_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"fcrd", 
    										"equilibrium"	);
    
    em_data hcrd_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"hcrd", 
    										"equilibrium"	);
    
    em_data mncrd_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"mncrd", 
    										"equilibrium"	);
    
    SS_ref_db.gbase[0] 		= crd_eq.gb;
    SS_ref_db.gbase[1] 		= fcrd_eq.gb;
    SS_ref_db.gbase[2] 		= hcrd_eq.gb;
    SS_ref_db.gbase[3] 		= mncrd_eq.gb - 4.21;
    
    SS_ref_db.ElShearMod[0] 	= crd_eq.ElShearMod;
    SS_ref_db.ElShearMod[1] 	= fcrd_eq.ElShearMod;
    SS_ref_db.ElShearMod[2] 	= hcrd_eq.ElShearMod;
    SS_ref_db.ElShearMod[3] 	= mncrd_eq.ElShearMod;
    
    for (i = 0; i < len_ox; i++){
        SS_ref_db.Comp[0][i] 	= crd_eq.C[i];
        SS_ref_db.Comp[1][i] 	= fcrd_eq.C[i];
        SS_ref_db.Comp[2][i] 	= hcrd_eq.C[i];
        SS_ref_db.Comp[3][i] 	= mncrd_eq.C[i];
    }
    
    for (i = 0; i < n_em; i++){
        SS_ref_db.z_em[i] = 1.0;
    };
    
    SS_ref_db.bounds_ref[0][0] = 0.0+eps;  SS_ref_db.bounds_ref[0][1] = 1.0-eps;
    SS_ref_db.bounds_ref[1][0] = 0.0+eps;  SS_ref_db.bounds_ref[1][1] = 1.0-eps;
    SS_ref_db.bounds_ref[2][0] = 0.0+eps;  SS_ref_db.bounds_ref[2][1] = 1.0-eps;
    
    return SS_ref_db;
}

/**
   retrieve reference thermodynamic data for mp_chl
*/
SS_ref G_SS_mp_chl_function(SS_ref SS_ref_db, int EM_database, int len_ox, bulk_info z_b, double eps){
    
    int i, j;
    int n_em = SS_ref_db.n_em;
    
    char   *EM_tmp[] 		= {"clin","afchl","ames","daph","ochl1","ochl4","f3clin","mmchl"};
    for (int i = 0; i < SS_ref_db.n_em; i++){
        strcpy(SS_ref_db.EM_list[i],EM_tmp[i]);
    };
    int n_xeos = SS_ref_db.n_xeos;
    char   *CV_tmp[] 		= {"x","y","f","m","QAl","Q1","Q4"};
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
        strcpy(SS_ref_db.CV_list[i],CV_tmp[i]);
    };
    SS_ref_db.W[0] = 17.0;
    SS_ref_db.W[1] = 17.0;
    SS_ref_db.W[2] = 20.0;
    SS_ref_db.W[3] = 30.0;
    SS_ref_db.W[4] = 21.0;
    SS_ref_db.W[5] = 2.0;
    SS_ref_db.W[6] = 6.0;
    SS_ref_db.W[7] = 16.0;
    SS_ref_db.W[8] = 37.0;
    SS_ref_db.W[9] = 20.0;
    SS_ref_db.W[10] = 4.0;
    SS_ref_db.W[11] = 15.0;
    SS_ref_db.W[12] = 23.0;
    SS_ref_db.W[13] = 30.0;
    SS_ref_db.W[14] = 29.0;
    SS_ref_db.W[15] = 13.0;
    SS_ref_db.W[16] = 19.0;
    SS_ref_db.W[17] = 17.0;
    SS_ref_db.W[18] = 18.0;
    SS_ref_db.W[19] = 33.0;
    SS_ref_db.W[20] = 22.0;
    SS_ref_db.W[21] = 4.0;
    SS_ref_db.W[22] = 24.0;
    SS_ref_db.W[23] = 28.6;
    SS_ref_db.W[24] = 19.0;
    SS_ref_db.W[25] = 19.0;
    SS_ref_db.W[26] = 22.0;
    SS_ref_db.W[27] = 8.0;
    
    
    em_data clin_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"clin", 
    										"equilibrium"	);
    
    em_data afchl_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"afchl", 
    										"equilibrium"	);
    
    em_data ames_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"ames", 
    										"equilibrium"	);
    
    em_data daph_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"daph", 
    										"equilibrium"	);
    
    em_data gr_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"gr", 
    										"equilibrium"	);
    
    em_data andr_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"andr", 
    										"equilibrium"	);
    
    em_data mnchl_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"mnchl", 
    										"equilibrium"	);
    
    SS_ref_db.gbase[0] 		= clin_eq.gb;
    SS_ref_db.gbase[1] 		= afchl_eq.gb;
    SS_ref_db.gbase[2] 		= ames_eq.gb;
    SS_ref_db.gbase[3] 		= daph_eq.gb;
    SS_ref_db.gbase[4] 		= afchl_eq.gb - clin_eq.gb + daph_eq.gb + 3.0;
    SS_ref_db.gbase[5] 		= afchl_eq.gb - 0.2*clin_eq.gb + 0.2*daph_eq.gb + 2.4;
    SS_ref_db.gbase[6] 		= 0.5*andr_eq.gb + clin_eq.gb - 0.5*gr_eq.gb + 2.0;
    SS_ref_db.gbase[7] 		= mnchl_eq.gb - 5.67;
    
    SS_ref_db.ElShearMod[0] 	= clin_eq.ElShearMod;
    SS_ref_db.ElShearMod[1] 	= afchl_eq.ElShearMod;
    SS_ref_db.ElShearMod[2] 	= ames_eq.ElShearMod;
    SS_ref_db.ElShearMod[3] 	= daph_eq.ElShearMod;
    SS_ref_db.ElShearMod[4] 	= afchl_eq.ElShearMod - clin_eq.ElShearMod + daph_eq.ElShearMod;
    SS_ref_db.ElShearMod[5] 	= afchl_eq.ElShearMod - 0.2*clin_eq.ElShearMod + 0.2*daph_eq.ElShearMod;
    SS_ref_db.ElShearMod[6] 	= 0.5*andr_eq.ElShearMod + clin_eq.ElShearMod - 0.5*gr_eq.ElShearMod;
    SS_ref_db.ElShearMod[7] 	= mnchl_eq.ElShearMod;
    
    for (i = 0; i < len_ox; i++){
        SS_ref_db.Comp[0][i] 	= clin_eq.C[i];
        SS_ref_db.Comp[1][i] 	= afchl_eq.C[i];
        SS_ref_db.Comp[2][i] 	= ames_eq.C[i];
        SS_ref_db.Comp[3][i] 	= daph_eq.C[i];
        SS_ref_db.Comp[4][i] 	= afchl_eq.C[i] - clin_eq.C[i] + daph_eq.C[i];
        SS_ref_db.Comp[5][i] 	= afchl_eq.C[i] - 0.2*clin_eq.C[i] + 0.2*daph_eq.C[i];
        SS_ref_db.Comp[6][i] 	= 0.5*andr_eq.C[i] + clin_eq.C[i] - 0.5*gr_eq.C[i];
        SS_ref_db.Comp[7][i] 	= mnchl_eq.C[i];
    }
    
    for (i = 0; i < n_em; i++){
        SS_ref_db.z_em[i] = 1.0;
    };
    
    SS_ref_db.bounds_ref[0][0] = 0.0+eps;  SS_ref_db.bounds_ref[0][1] = 1.0-eps;
    SS_ref_db.bounds_ref[1][0] = 0.0+eps;  SS_ref_db.bounds_ref[1][1] = 1.0-eps;
    SS_ref_db.bounds_ref[2][0] = 0.0+eps;  SS_ref_db.bounds_ref[2][1] = 1.0-eps;
    SS_ref_db.bounds_ref[3][0] = 0.0+eps;  SS_ref_db.bounds_ref[3][1] = 1.0-eps;
    SS_ref_db.bounds_ref[4][0] = -1.0+eps;  SS_ref_db.bounds_ref[4][1] = 1.0-eps;
    SS_ref_db.bounds_ref[5][0] = -1.0+eps;  SS_ref_db.bounds_ref[5][1] = 1.0-eps;
    SS_ref_db.bounds_ref[6][0] = -1.0+eps;  SS_ref_db.bounds_ref[6][1] = 1.0-eps;
	if (z_b.bulk_rock[8] == 0.){ 					
		SS_ref_db.z_em[6]          = 0.0;
        SS_ref_db.d_em[6]          = 1.0;
		SS_ref_db.bounds_ref[2][0] = 0.0; 
		SS_ref_db.bounds_ref[2][1] = 0.0;	
	}

    /* this lists the index of the order variables */
    // SS_ref_db.orderVar      = 1;
    // SS_ref_db.idOrderVar[3] = -1.0;
    // SS_ref_db.idOrderVar[4] = -1.0;   
    // SS_ref_db.idOrderVar[5] = -1.0;
    return SS_ref_db;
}

/**
   retrieve reference thermodynamic data for mp_ctd
*/
SS_ref G_SS_mp_ctd_function(SS_ref SS_ref_db, int EM_database, int len_ox, bulk_info z_b, double eps){
    
    int i, j;
    int n_em = SS_ref_db.n_em;
    
    char   *EM_tmp[] 		= {"mctd","fctd","mnct","ctdo"};
    for (int i = 0; i < SS_ref_db.n_em; i++){
        strcpy(SS_ref_db.EM_list[i],EM_tmp[i]);
    };
    int n_xeos = SS_ref_db.n_xeos;
    char   *CV_tmp[] 		= {"x","m","f"};
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
        strcpy(SS_ref_db.CV_list[i],CV_tmp[i]);
    };
    SS_ref_db.W[0] = 4.0;
    SS_ref_db.W[1] = 3.0;
    SS_ref_db.W[2] = 1.0;
    SS_ref_db.W[3] = 3.0;
    SS_ref_db.W[4] = 5.0;
    SS_ref_db.W[5] = 4.0;
    
    
    em_data mctd_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"mctd", 
    										"equilibrium"	);
    
    em_data fctd_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"fctd", 
    										"equilibrium"	);
    
    em_data mnctd_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"mnctd", 
    										"equilibrium"	);
    
    em_data andr_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"andr", 
    										"equilibrium"	);
    
    em_data gr_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"gr", 
    										"equilibrium"	);
    
    SS_ref_db.gbase[0] 		= mctd_eq.gb;
    SS_ref_db.gbase[1] 		= fctd_eq.gb;
    SS_ref_db.gbase[2] 		= mnctd_eq.gb + 0.66;
    SS_ref_db.gbase[3] 		= 0.25*andr_eq.gb - 0.25*gr_eq.gb + mctd_eq.gb + 13.5;
    
    SS_ref_db.ElShearMod[0] 	= mctd_eq.ElShearMod;
    SS_ref_db.ElShearMod[1] 	= fctd_eq.ElShearMod;
    SS_ref_db.ElShearMod[2] 	= mnctd_eq.ElShearMod;
    SS_ref_db.ElShearMod[3] 	= 0.25*andr_eq.ElShearMod - 0.25*gr_eq.ElShearMod + mctd_eq.ElShearMod;
    
    for (i = 0; i < len_ox; i++){
        SS_ref_db.Comp[0][i] 	= mctd_eq.C[i];
        SS_ref_db.Comp[1][i] 	= fctd_eq.C[i];
        SS_ref_db.Comp[2][i] 	= mnctd_eq.C[i];
        SS_ref_db.Comp[3][i] 	= 0.25*andr_eq.C[i] - 0.25*gr_eq.C[i] + mctd_eq.C[i];
    }
    
    for (i = 0; i < n_em; i++){
        SS_ref_db.z_em[i] = 1.0;
    };
    
    SS_ref_db.bounds_ref[0][0] = 0.0+eps;  SS_ref_db.bounds_ref[0][1] = 1.0-eps;
    SS_ref_db.bounds_ref[1][0] = 0.0+eps;  SS_ref_db.bounds_ref[1][1] = 1.0-eps;
    SS_ref_db.bounds_ref[2][0] = 0.0+eps;  SS_ref_db.bounds_ref[2][1] = 1.0-eps;
	if (z_b.bulk_rock[8] == 0.){ 					
		SS_ref_db.z_em[3]          = 0.0;
        SS_ref_db.d_em[3]          = 1.0;
		SS_ref_db.bounds_ref[2][0] = 0.0; 
		SS_ref_db.bounds_ref[2][1] = 0.0;	
	}
    return SS_ref_db;
}

/**
   retrieve reference thermodynamic data for mp_ep
*/
SS_ref G_SS_mp_ep_function(SS_ref SS_ref_db, int EM_database, int len_ox, bulk_info z_b, double eps){
    
    int i, j;
    int n_em = SS_ref_db.n_em;
    
    char   *EM_tmp[] 		= {"cz","ep","fep"};
    for (int i = 0; i < SS_ref_db.n_em; i++){
        strcpy(SS_ref_db.EM_list[i],EM_tmp[i]);
    };
    int n_xeos = SS_ref_db.n_xeos;
    char   *CV_tmp[] 		= {"f","Q"};
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
        strcpy(SS_ref_db.CV_list[i],CV_tmp[i]);
    };
    SS_ref_db.W[0] = 1.0;
    SS_ref_db.W[1] = 3.0;
    SS_ref_db.W[2] = 1.0;
    
    
    em_data cz_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"cz", 
    										"equilibrium"	);
    
    em_data ep_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"ep", 
    										"equilibrium"	);
    
    em_data fep_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"fep", 
    										"equilibrium"	);
    
    SS_ref_db.gbase[0] 		= cz_eq.gb;
    SS_ref_db.gbase[1] 		= ep_eq.gb;
    SS_ref_db.gbase[2] 		= fep_eq.gb;
    
    SS_ref_db.ElShearMod[0] 	= cz_eq.ElShearMod;
    SS_ref_db.ElShearMod[1] 	= ep_eq.ElShearMod;
    SS_ref_db.ElShearMod[2] 	= fep_eq.ElShearMod;
    
    for (i = 0; i < len_ox; i++){
        SS_ref_db.Comp[0][i] 	= cz_eq.C[i];
        SS_ref_db.Comp[1][i] 	= ep_eq.C[i];
        SS_ref_db.Comp[2][i] 	= fep_eq.C[i];
    }
    
    for (i = 0; i < n_em; i++){
        SS_ref_db.z_em[i] = 1.0;
    };
    
    SS_ref_db.bounds_ref[0][0] = 0.0+eps;  SS_ref_db.bounds_ref[0][1] = 1.0-eps;
    SS_ref_db.bounds_ref[1][0] = 0.0+eps;  SS_ref_db.bounds_ref[1][1] = 0.5-eps;

	if (z_b.bulk_rock[8] == 0.){ 					
		SS_ref_db.bounds_ref[0][0] = 0.0; 
		SS_ref_db.bounds_ref[0][1] = 0.0;	
		SS_ref_db.bounds_ref[1][0] = 0.0; 
		SS_ref_db.bounds_ref[1][1] = 0.0;	
		SS_ref_db.z_em[1]          = 0.0;
        SS_ref_db.d_em[1]          = 1.0;
		SS_ref_db.z_em[2]          = 0.0;
        SS_ref_db.d_em[2]          = 1.0;
	}
    
    return SS_ref_db;
}

/**
   retrieve reference thermodynamic data for mp_g
*/
SS_ref G_SS_mp_g_function(SS_ref SS_ref_db, int EM_database, int len_ox, bulk_info z_b, double eps){
    
    int i, j;
    int n_em = SS_ref_db.n_em;
    
    char   *EM_tmp[] 		= {"py","alm","spss","gr","kho"};
    for (int i = 0; i < SS_ref_db.n_em; i++){
        strcpy(SS_ref_db.EM_list[i],EM_tmp[i]);
    };
    int n_xeos = SS_ref_db.n_xeos;
    char   *CV_tmp[] 		= {"x","z","m","f"};
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
        strcpy(SS_ref_db.CV_list[i],CV_tmp[i]);
    };
    SS_ref_db.W[0] = 2.5;
    SS_ref_db.W[1] = 2.0;
    SS_ref_db.W[2] = 31.0;
    SS_ref_db.W[3] = 5.4;
    SS_ref_db.W[4] = 2.0;
    SS_ref_db.W[5] = 5.0;
    SS_ref_db.W[6] = 22.6;
    SS_ref_db.W[7] = 0.0;
    SS_ref_db.W[8] = 29.4;
    SS_ref_db.W[9] = -15.3;
    
    SS_ref_db.v[0] = 1.0;
    SS_ref_db.v[1] = 1.0;
    SS_ref_db.v[2] = 1.0;
    SS_ref_db.v[3] = 2.7;
    SS_ref_db.v[4] = 1.0;
    
    
    em_data py_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"py", 
    										"equilibrium"	);
    
    em_data alm_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"alm", 
    										"equilibrium"	);
    
    em_data spss_eq 		= get_em_data(	EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"spss", 
    										"equilibrium"	);
    
    em_data gr_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"gr", 
    										"equilibrium"	);
    
    em_data andr_eq 		= get_em_data(	EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"andr", 
    										"equilibrium"	);
    
    SS_ref_db.gbase[0] 		= py_eq.gb;
    SS_ref_db.gbase[1] 		= alm_eq.gb;
    SS_ref_db.gbase[2] 		= spss_eq.gb;
    SS_ref_db.gbase[3] 		= gr_eq.gb;
    SS_ref_db.gbase[4] 		= andr_eq.gb - gr_eq.gb + py_eq.gb + 27.0;
    
    SS_ref_db.ElShearMod[0] 	= py_eq.ElShearMod;
    SS_ref_db.ElShearMod[1] 	= alm_eq.ElShearMod;
    SS_ref_db.ElShearMod[2] 	= spss_eq.ElShearMod;
    SS_ref_db.ElShearMod[3] 	= gr_eq.ElShearMod;
    SS_ref_db.ElShearMod[4] 	= andr_eq.ElShearMod - gr_eq.ElShearMod + py_eq.ElShearMod;
    
    for (i = 0; i < len_ox; i++){
        SS_ref_db.Comp[0][i] 	= py_eq.C[i];
        SS_ref_db.Comp[1][i] 	= alm_eq.C[i];
        SS_ref_db.Comp[2][i] 	= spss_eq.C[i];
        SS_ref_db.Comp[3][i] 	= gr_eq.C[i];
        SS_ref_db.Comp[4][i] 	= andr_eq.C[i] - gr_eq.C[i] + py_eq.C[i];
    }
    
    for (i = 0; i < n_em; i++){
        SS_ref_db.z_em[i] = 1.0;
    };
    
    SS_ref_db.bounds_ref[0][0] = 0.0+eps;  SS_ref_db.bounds_ref[0][1] = 1.0-eps;
    SS_ref_db.bounds_ref[1][0] = 0.0+eps;  SS_ref_db.bounds_ref[1][1] = 1.0-eps;
    SS_ref_db.bounds_ref[2][0] = 0.0+eps;  SS_ref_db.bounds_ref[2][1] = 1.0-eps;
    SS_ref_db.bounds_ref[3][0] = 0.0+eps;  SS_ref_db.bounds_ref[3][1] = 1.0-eps;

	if (z_b.bulk_rock[8] == 0.){ 					
		SS_ref_db.z_em[4]          = 0.0;
        SS_ref_db.d_em[4]          = 1.0;
		SS_ref_db.bounds_ref[3][0] = 0.0; 
		SS_ref_db.bounds_ref[3][1] = 0.0;	
	}
    return SS_ref_db;
}


/**
   retrieve reference thermodynamic data for mb_ilm
*/
SS_ref G_SS_mp_ilm_function(SS_ref SS_ref_db, int EM_database, int len_ox, bulk_info z_b, double eps){
    
    int i, j;
    int n_em = SS_ref_db.n_em;
    
    char   *EM_tmp[] 		= {"oilm","dilm","dhem"};
    for (int i = 0; i < SS_ref_db.n_em; i++){
        strcpy(SS_ref_db.EM_list[i],EM_tmp[i]);
    };
    int n_xeos = SS_ref_db.n_xeos;
    char   *CV_tmp[] 		= {"x","Q"};
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
        strcpy(SS_ref_db.CV_list[i],CV_tmp[i]);
    };
    SS_ref_db.W[0] = 15.6;
    SS_ref_db.W[1] = 26.6;
    SS_ref_db.W[2] = 11.0;
    
    
    em_data ilm_di 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"ilm", 
    										"disordered"	);
    
    em_data hem_di 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"hem", 
    										"disordered"	);
    
    SS_ref_db.gbase[0] 		= 0.009426*z_b.T + ilm_di.gb - 13.6075;
    SS_ref_db.gbase[1] 		= -0.0021*z_b.T + ilm_di.gb + 1.9928;
    SS_ref_db.gbase[2] 		= hem_di.gb;
    
    SS_ref_db.ElShearMod[0] 	= ilm_di.ElShearMod;
    SS_ref_db.ElShearMod[1] 	= ilm_di.ElShearMod;
    SS_ref_db.ElShearMod[2] 	= hem_di.ElShearMod;
    
    for (i = 0; i < len_ox; i++){
        SS_ref_db.Comp[0][i] 	= ilm_di.C[i];
        SS_ref_db.Comp[1][i] 	= ilm_di.C[i];
        SS_ref_db.Comp[2][i] 	= hem_di.C[i];
    }
    
    for (i = 0; i < n_em; i++){
        SS_ref_db.z_em[i] = 1.0;
    };
    
    SS_ref_db.bounds_ref[0][0] = 0.0+eps;  SS_ref_db.bounds_ref[0][1] = 1.0-eps;
    SS_ref_db.bounds_ref[1][0] = -0.99+eps;  SS_ref_db.bounds_ref[1][1] = 0.99-eps;
    

    /* this lists the index of the order variables */
    SS_ref_db.orderVar      = 1;
    SS_ref_db.idOrderVar[1] = -1.0;


    return SS_ref_db;
}


/**
   retrieve reference thermodynamic data for mp_ilmm
*/
SS_ref G_SS_mp_ilmm_function(SS_ref SS_ref_db, int EM_database, int len_ox, bulk_info z_b, double eps){
    
    int i, j;
    int n_em = SS_ref_db.n_em;
    
    char   *EM_tmp[] 		= {"oilm","dilm","dhem","geik","pnt"};
    for (int i = 0; i < SS_ref_db.n_em; i++){
        strcpy(SS_ref_db.EM_list[i],EM_tmp[i]);
    };
    int n_xeos = SS_ref_db.n_xeos;
    char   *CV_tmp[] 		= {"i","g","m","Q"};
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
        strcpy(SS_ref_db.CV_list[i],CV_tmp[i]);
    };
    SS_ref_db.W[0] = 15.6;
    SS_ref_db.W[1] = 26.6;
    SS_ref_db.W[2] = 4.0;
    SS_ref_db.W[3] = 2.0;
    SS_ref_db.W[4] = 11.0;
    SS_ref_db.W[5] = 4.0;
    SS_ref_db.W[6] = 2.0;
    SS_ref_db.W[7] = 36.0;
    SS_ref_db.W[8] = 25.0;
    SS_ref_db.W[9] = 4.0;
    
    
    em_data ilm_di 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"ilm", 
    										"disordered"	);
    
    em_data hem_di 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"hem", 
    										"disordered"	);
    
    em_data geik_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"geik", 
    										"equilibrium"	);
    
    em_data pnt_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"pnt", 
    										"equilibrium"	);
    
    SS_ref_db.gbase[0] 		= 0.009426*z_b.T + ilm_di.gb - 13.6075;
    SS_ref_db.gbase[1] 		= -0.0021*z_b.T + ilm_di.gb + 1.9928;
    SS_ref_db.gbase[2] 		= hem_di.gb;
    SS_ref_db.gbase[3] 		= geik_eq.gb;
    SS_ref_db.gbase[4] 		= pnt_eq.gb;
    
    SS_ref_db.ElShearMod[0] 	= ilm_di.ElShearMod;
    SS_ref_db.ElShearMod[1] 	= ilm_di.ElShearMod;
    SS_ref_db.ElShearMod[2] 	= hem_di.ElShearMod;
    SS_ref_db.ElShearMod[3] 	= geik_eq.ElShearMod;
    SS_ref_db.ElShearMod[4] 	= pnt_eq.ElShearMod;
    
    for (i = 0; i < len_ox; i++){
        SS_ref_db.Comp[0][i] 	= ilm_di.C[i];
        SS_ref_db.Comp[1][i] 	= ilm_di.C[i];
        SS_ref_db.Comp[2][i] 	= hem_di.C[i];
        SS_ref_db.Comp[3][i] 	= geik_eq.C[i];
        SS_ref_db.Comp[4][i] 	= pnt_eq.C[i];
    }
    
    for (i = 0; i < n_em; i++){
        SS_ref_db.z_em[i] = 1.0;
    };
    
    SS_ref_db.bounds_ref[0][0] = 0.0+eps;  SS_ref_db.bounds_ref[0][1] = 1.0-eps;
    SS_ref_db.bounds_ref[1][0] = 0.0+eps;  SS_ref_db.bounds_ref[1][1] = 1.0-eps;
    SS_ref_db.bounds_ref[2][0] = 0.0+eps;  SS_ref_db.bounds_ref[2][1] = 1.0-eps;
    SS_ref_db.bounds_ref[3][0] = -1.0+eps;  SS_ref_db.bounds_ref[3][1] = 1.0-eps;
	if (z_b.bulk_rock[8] == 0.){ 					
		SS_ref_db.z_em[2]          = 0.0;
        SS_ref_db.d_em[2]          = 1.0;
		SS_ref_db.bounds_ref[0][0] = 1.0; 
		SS_ref_db.bounds_ref[0][1] = 1.0;	
	}


    /* this lists the index of the order variables */
    SS_ref_db.orderVar      = 1;
    SS_ref_db.idOrderVar[3] = -1.0;
    
    return SS_ref_db;
}

/**
   retrieve reference thermodynamic data for mp_ma
*/
SS_ref G_SS_mp_ma_function(SS_ref SS_ref_db, int EM_database, int len_ox, bulk_info z_b, double eps){
    
    int i, j;
    int n_em = SS_ref_db.n_em;
    
    char   *EM_tmp[] 		= {"mut","celt","fcelt","pat","ma","fmu"};
    for (int i = 0; i < SS_ref_db.n_em; i++){
        strcpy(SS_ref_db.EM_list[i],EM_tmp[i]);
    };
    int n_xeos = SS_ref_db.n_xeos;
    char   *CV_tmp[] 		= {"x","y","f","n","c"};
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
        strcpy(SS_ref_db.CV_list[i],CV_tmp[i]);
    };
    SS_ref_db.W[0] = 0.2*SS_ref_db.P;
    SS_ref_db.W[1] = 0.2*SS_ref_db.P;
    SS_ref_db.W[2] = 0.353*SS_ref_db.P + 0.0034*SS_ref_db.T + 10.12;
    SS_ref_db.W[3] = 34.0;
    SS_ref_db.W[4] = 0.0;
    SS_ref_db.W[5] = 0.0;
    SS_ref_db.W[6] = 0.25*SS_ref_db.P + 45.0;
    SS_ref_db.W[7] = 50.0;
    SS_ref_db.W[8] = 0.0;
    SS_ref_db.W[9] = 0.25*SS_ref_db.P + 45.0;
    SS_ref_db.W[10] = 50.0;
    SS_ref_db.W[11] = 0.0;
    SS_ref_db.W[12] = 18.0;
    SS_ref_db.W[13] = 30.0;
    SS_ref_db.W[14] = 35.0;
    
    SS_ref_db.v[0] = 0.63;
    SS_ref_db.v[1] = 0.63;
    SS_ref_db.v[2] = 0.63;
    SS_ref_db.v[3] = 0.37;
    SS_ref_db.v[4] = 0.63;
    SS_ref_db.v[5] = 0.63;
    
    
    em_data mu_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"mu", 
    										"equilibrium"	);
    
    em_data cel_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"cel", 
    										"equilibrium"	);
    
    em_data fcel_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"fcel", 
    										"equilibrium"	);
    
    em_data pa_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"pa", 
    										"equilibrium"	);
    
    em_data ma_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"ma", 
    										"equilibrium"	);
    
    em_data andr_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"andr", 
    										"equilibrium"	);
    
    em_data gr_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"gr", 
    										"equilibrium"	);
    
    SS_ref_db.gbase[0] 		= mu_eq.gb + 1.0;
    SS_ref_db.gbase[1] 		= cel_eq.gb + 5.0;
    SS_ref_db.gbase[2] 		= fcel_eq.gb + 5.0;
    SS_ref_db.gbase[3] 		= pa_eq.gb + 4.0;
    SS_ref_db.gbase[4] 		= ma_eq.gb;
    SS_ref_db.gbase[5] 		= 0.5*andr_eq.gb - 0.5*gr_eq.gb + mu_eq.gb + 25.0;
    
    SS_ref_db.ElShearMod[0] 	= mu_eq.ElShearMod;
    SS_ref_db.ElShearMod[1] 	= cel_eq.ElShearMod;
    SS_ref_db.ElShearMod[2] 	= fcel_eq.ElShearMod;
    SS_ref_db.ElShearMod[3] 	= pa_eq.ElShearMod;
    SS_ref_db.ElShearMod[4] 	= ma_eq.ElShearMod;
    SS_ref_db.ElShearMod[5] 	= 0.5*andr_eq.ElShearMod - 0.5*gr_eq.ElShearMod + mu_eq.ElShearMod;
    
    for (i = 0; i < len_ox; i++){
        SS_ref_db.Comp[0][i] 	= mu_eq.C[i];
        SS_ref_db.Comp[1][i] 	= cel_eq.C[i];
        SS_ref_db.Comp[2][i] 	= fcel_eq.C[i];
        SS_ref_db.Comp[3][i] 	= pa_eq.C[i];
        SS_ref_db.Comp[4][i] 	= ma_eq.C[i];
        SS_ref_db.Comp[5][i] 	= 0.5*andr_eq.C[i] - 0.5*gr_eq.C[i] + mu_eq.C[i];
    }
    
    for (i = 0; i < n_em; i++){
        SS_ref_db.z_em[i] = 1.0;
    };
    
    SS_ref_db.bounds_ref[0][0] = 0.0+eps;  SS_ref_db.bounds_ref[0][1] = 1.0-eps;
    SS_ref_db.bounds_ref[1][0] = 0.0+eps;  SS_ref_db.bounds_ref[1][1] = 1.0-eps;
    SS_ref_db.bounds_ref[2][0] = 0.0+eps;  SS_ref_db.bounds_ref[2][1] = 1.0-eps;
    SS_ref_db.bounds_ref[3][0] = 0.0+eps;  SS_ref_db.bounds_ref[3][1] = 1.0-eps;
    SS_ref_db.bounds_ref[4][0] = 0.0+eps;  SS_ref_db.bounds_ref[4][1] = 1.0-eps;
	if (z_b.bulk_rock[8] == 0.){ 					
		SS_ref_db.z_em[5]          = 0.0;
        SS_ref_db.d_em[5]          = 1.0;
		SS_ref_db.bounds_ref[2][0] = 0.0; 
		SS_ref_db.bounds_ref[2][1] = 0.0;	
	}
    return SS_ref_db;
}

/**
   retrieve reference thermodynamic data for mp_mt
*/
SS_ref G_SS_mp_mt_function(SS_ref SS_ref_db, int EM_database, int len_ox, bulk_info z_b, double eps){
    
    int i, j;
    int n_em = SS_ref_db.n_em;
    
    char   *EM_tmp[] 		= {"imt","dmt","usp"};
    for (int i = 0; i < SS_ref_db.n_em; i++){
        strcpy(SS_ref_db.EM_list[i],EM_tmp[i]);
    };
    int n_xeos = SS_ref_db.n_xeos;
    char   *CV_tmp[] 		= {"x","Q"};
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
        strcpy(SS_ref_db.CV_list[i],CV_tmp[i]);
    };
    SS_ref_db.W[0] = 2.4;
    SS_ref_db.W[1] = 1.0;
    SS_ref_db.W[2] = -5.0;
    
    
    em_data mt_di 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"mt", 
    										"disordered"	);
    
    em_data usp_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"usp", 
    										"equilibrium"	);
    
    SS_ref_db.gbase[0] 		= 0.003166*z_b.T + mt_di.gb - 1.8595;
    SS_ref_db.gbase[1] 		= -0.001184*z_b.T + mt_di.gb + 1.3305;
    SS_ref_db.gbase[2] 		= usp_eq.gb;
    
    SS_ref_db.ElShearMod[0] 	= mt_di.ElShearMod;
    SS_ref_db.ElShearMod[1] 	= mt_di.ElShearMod;
    SS_ref_db.ElShearMod[2] 	= usp_eq.ElShearMod;
    
    for (i = 0; i < len_ox; i++){
        SS_ref_db.Comp[0][i] 	= mt_di.C[i];
        SS_ref_db.Comp[1][i] 	= mt_di.C[i];
        SS_ref_db.Comp[2][i] 	= usp_eq.C[i];
    }
    
    for (i = 0; i < n_em; i++){
        SS_ref_db.z_em[i] = 1.0;
    };
    
    SS_ref_db.bounds_ref[0][0] = 0.0+eps;  SS_ref_db.bounds_ref[0][1] = 1.0-eps;
    SS_ref_db.bounds_ref[1][0] = 0.0+eps;  SS_ref_db.bounds_ref[1][1] = 1.0-eps;

	if (z_b.bulk_rock[8] == 0.){ 					
		SS_ref_db.z_em[0]          = 0.0;
        SS_ref_db.d_em[0]          = 1.0;
		SS_ref_db.z_em[1]          = 0.0;
        SS_ref_db.d_em[1]          = 1.0;
		SS_ref_db.bounds_ref[0][0] = 0.0; 
		SS_ref_db.bounds_ref[0][1] = 0.0;	
		SS_ref_db.bounds_ref[1][0] = 0.0; 
		SS_ref_db.bounds_ref[1][1] = 0.0;	
	}

    return SS_ref_db;
}

/**
   retrieve reference thermodynamic data for mp_mu
*/
SS_ref G_SS_mp_mu_function(SS_ref SS_ref_db, int EM_database, int len_ox, bulk_info z_b, double eps){
    
    int i, j;
    int n_em = SS_ref_db.n_em;
    
    char   *EM_tmp[] 		= {"mut","cel","fcel","pat","ma","fmu"};
    for (int i = 0; i < SS_ref_db.n_em; i++){
        strcpy(SS_ref_db.EM_list[i],EM_tmp[i]);
    };
    int n_xeos = SS_ref_db.n_xeos;
    char   *CV_tmp[] 		= {"x","y","f","n","c"};
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
        strcpy(SS_ref_db.CV_list[i],CV_tmp[i]);
    };
    SS_ref_db.W[0] = 0.2*SS_ref_db.P;
    SS_ref_db.W[1] = 0.2*SS_ref_db.P;
    SS_ref_db.W[2] = 0.353*SS_ref_db.P + 0.0034*SS_ref_db.T + 10.12;
    SS_ref_db.W[3] = 35.0;
    SS_ref_db.W[4] = 0.0;
    SS_ref_db.W[5] = 0.0;
    SS_ref_db.W[6] = 0.25*SS_ref_db.P + 45.0;
    SS_ref_db.W[7] = 50.0;
    SS_ref_db.W[8] = 0.0;
    SS_ref_db.W[9] = 0.25*SS_ref_db.P + 45.0;
    SS_ref_db.W[10] = 50.0;
    SS_ref_db.W[11] = 0.0;
    SS_ref_db.W[12] = 15.0;
    SS_ref_db.W[13] = 30.0;
    SS_ref_db.W[14] = 35.0;
    
    SS_ref_db.v[0] = 0.63;
    SS_ref_db.v[1] = 0.63;
    SS_ref_db.v[2] = 0.63;
    SS_ref_db.v[3] = 0.37;
    SS_ref_db.v[4] = 0.63;
    SS_ref_db.v[5] = 0.63;
    
    
    em_data mu_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"mu", 
    										"equilibrium"	);
    
    em_data cel_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"cel", 
    										"equilibrium"	);
    
    em_data fcel_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"fcel", 
    										"equilibrium"	);
    
    em_data pa_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"pa", 
    										"equilibrium"	);
    
    em_data ma_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"ma", 
    										"equilibrium"	);
    
    em_data andr_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"andr", 
    										"equilibrium"	);
    
    em_data gr_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"gr", 
    										"equilibrium"	);
    
    SS_ref_db.gbase[0] 		= mu_eq.gb;
    SS_ref_db.gbase[1] 		= cel_eq.gb;
    SS_ref_db.gbase[2] 		= fcel_eq.gb;
    SS_ref_db.gbase[3] 		= pa_eq.gb;
    SS_ref_db.gbase[4] 		= ma_eq.gb + 5.0;
    SS_ref_db.gbase[5] 		= 0.5*andr_eq.gb - 0.5*gr_eq.gb + mu_eq.gb + 25.0;
    
    SS_ref_db.ElShearMod[0] 	= mu_eq.ElShearMod;
    SS_ref_db.ElShearMod[1] 	= cel_eq.ElShearMod;
    SS_ref_db.ElShearMod[2] 	= fcel_eq.ElShearMod;
    SS_ref_db.ElShearMod[3] 	= pa_eq.ElShearMod;
    SS_ref_db.ElShearMod[4] 	= ma_eq.ElShearMod;
    SS_ref_db.ElShearMod[5] 	= 0.5*andr_eq.ElShearMod - 0.5*gr_eq.ElShearMod + mu_eq.ElShearMod;
    
    for (i = 0; i < len_ox; i++){
        SS_ref_db.Comp[0][i] 	= mu_eq.C[i];
        SS_ref_db.Comp[1][i] 	= cel_eq.C[i];
        SS_ref_db.Comp[2][i] 	= fcel_eq.C[i];
        SS_ref_db.Comp[3][i] 	= pa_eq.C[i];
        SS_ref_db.Comp[4][i] 	= ma_eq.C[i];
        SS_ref_db.Comp[5][i] 	= 0.5*andr_eq.C[i] - 0.5*gr_eq.C[i] + mu_eq.C[i];
    }
    
    for (i = 0; i < n_em; i++){
        SS_ref_db.z_em[i] = 1.0;
    };
    
    SS_ref_db.bounds_ref[0][0] = 0.0+eps;  SS_ref_db.bounds_ref[0][1] = 1.0-eps;
    SS_ref_db.bounds_ref[1][0] = 0.0+eps;  SS_ref_db.bounds_ref[1][1] = 1.0-eps;
    SS_ref_db.bounds_ref[2][0] = 0.0+eps;  SS_ref_db.bounds_ref[2][1] = 1.0-eps;
    SS_ref_db.bounds_ref[3][0] = 0.0+eps;  SS_ref_db.bounds_ref[3][1] = 1.0-eps;
    SS_ref_db.bounds_ref[4][0] = 0.0+eps;  SS_ref_db.bounds_ref[4][1] = 1.0-eps;
	if (z_b.bulk_rock[8] == 0.){ 					
		SS_ref_db.z_em[5]          = 0.0;
        SS_ref_db.d_em[5]          = 1.0;
		SS_ref_db.bounds_ref[2][0] = 0.0; 
		SS_ref_db.bounds_ref[2][1] = 0.0;	
	}
    return SS_ref_db;
}

/**
   retrieve reference thermodynamic data for mp_opx
*/
SS_ref G_SS_mp_opx_function(SS_ref SS_ref_db, int EM_database, int len_ox, bulk_info z_b, double eps){
    
    int i, j;
    int n_em = SS_ref_db.n_em;
    
    char   *EM_tmp[] 		= {"en","fs","fm","mgts","fopx","mnopx","odi"};
    for (int i = 0; i < SS_ref_db.n_em; i++){
        strcpy(SS_ref_db.EM_list[i],EM_tmp[i]);
    };
    int n_xeos = SS_ref_db.n_xeos;
    char   *CV_tmp[] 		= {"x","m","y","f","c","Q"};
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
        strcpy(SS_ref_db.CV_list[i],CV_tmp[i]);
    };
    SS_ref_db.W[0] = 7.0;
    SS_ref_db.W[1] = 4.0;
    SS_ref_db.W[2] = 13.0 - 0.15*SS_ref_db.P;
    SS_ref_db.W[3] = 11.0 - 0.15*SS_ref_db.P;
    SS_ref_db.W[4] = 5.0;
    SS_ref_db.W[5] = 0.12*SS_ref_db.P + 32.2;
    SS_ref_db.W[6] = 4.0;
    SS_ref_db.W[7] = 13.0 - 0.15*SS_ref_db.P;
    SS_ref_db.W[8] = 11.6 - 0.15*SS_ref_db.P;
    SS_ref_db.W[9] = 4.2;
    SS_ref_db.W[10] = 0.084*SS_ref_db.P + 25.54;
    SS_ref_db.W[11] = 17.0 - 0.15*SS_ref_db.P;
    SS_ref_db.W[12] = 15.0 - 0.15*SS_ref_db.P;
    SS_ref_db.W[13] = 5.1;
    SS_ref_db.W[14] = 0.084*SS_ref_db.P + 22.54;
    SS_ref_db.W[15] = 1.0;
    SS_ref_db.W[16] = 12.0 - 0.15*SS_ref_db.P;
    SS_ref_db.W[17] = 75.4 - 0.94*SS_ref_db.P;
    SS_ref_db.W[18] = 10.6 - 0.15*SS_ref_db.P;
    SS_ref_db.W[19] = 73.4 - 0.94*SS_ref_db.P;
    SS_ref_db.W[20] = 0.084*SS_ref_db.P + 24.54;
    
    SS_ref_db.v[0] = 1.0;
    SS_ref_db.v[1] = 1.0;
    SS_ref_db.v[2] = 1.0;
    SS_ref_db.v[3] = 1.0;
    SS_ref_db.v[4] = 1.0;
    SS_ref_db.v[5] = 1.0;
    SS_ref_db.v[6] = 1.2;
    
    
    em_data en_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"en", 
    										"equilibrium"	);
    
    em_data fs_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"fs", 
    										"equilibrium"	);
    
    em_data mgts_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"mgts", 
    										"equilibrium"	);
    
    em_data andr_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"andr", 
    										"equilibrium"	);
    
    em_data gr_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"gr", 
    										"equilibrium"	);
    
    em_data pxmn_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"pxmn", 
    										"equilibrium"	);
    
    em_data di_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"di", 
    										"equilibrium"	);
    
    SS_ref_db.gbase[0] 		= en_eq.gb;
    SS_ref_db.gbase[1] 		= fs_eq.gb;
    SS_ref_db.gbase[2] 		= 0.5*en_eq.gb + 0.5*fs_eq.gb - 6.6;
    SS_ref_db.gbase[3] 		= mgts_eq.gb;
    SS_ref_db.gbase[4] 		= 0.5*andr_eq.gb - 0.5*gr_eq.gb + mgts_eq.gb + 2.0;
    SS_ref_db.gbase[5] 		= 2.0*pxmn_eq.gb + 6.68;
    SS_ref_db.gbase[6] 		= 0.005*z_b.P + 0.000211*z_b.T + di_eq.gb - 0.1;
    
    SS_ref_db.ElShearMod[0] 	= en_eq.ElShearMod;
    SS_ref_db.ElShearMod[1] 	= fs_eq.ElShearMod;
    SS_ref_db.ElShearMod[2] 	= 0.5*en_eq.ElShearMod + 0.5*fs_eq.ElShearMod;
    SS_ref_db.ElShearMod[3] 	= mgts_eq.ElShearMod;
    SS_ref_db.ElShearMod[4] 	= 0.5*andr_eq.ElShearMod - 0.5*gr_eq.ElShearMod + mgts_eq.ElShearMod;
    SS_ref_db.ElShearMod[5] 	= 2.0*pxmn_eq.ElShearMod;
    SS_ref_db.ElShearMod[6] 	= di_eq.ElShearMod;
    
    for (i = 0; i < len_ox; i++){
        SS_ref_db.Comp[0][i] 	= en_eq.C[i];
        SS_ref_db.Comp[1][i] 	= fs_eq.C[i];
        SS_ref_db.Comp[2][i] 	= 0.5*en_eq.C[i] + 0.5*fs_eq.C[i];
        SS_ref_db.Comp[3][i] 	= mgts_eq.C[i];
        SS_ref_db.Comp[4][i] 	= 0.5*andr_eq.C[i] - 0.5*gr_eq.C[i] + mgts_eq.C[i];
        SS_ref_db.Comp[5][i] 	= 2.0*pxmn_eq.C[i];
        SS_ref_db.Comp[6][i] 	= di_eq.C[i];
    }
    
    for (i = 0; i < n_em; i++){
        SS_ref_db.z_em[i] = 1.0;
    };
    
    SS_ref_db.bounds_ref[0][0] = 0.0+eps;  SS_ref_db.bounds_ref[0][1] = 1.0-eps;
    SS_ref_db.bounds_ref[1][0] = 0.0+eps;  SS_ref_db.bounds_ref[1][1] = 1.0-eps;
    SS_ref_db.bounds_ref[2][0] = 0.0+eps;  SS_ref_db.bounds_ref[2][1] = 2.0-eps;
    SS_ref_db.bounds_ref[3][0] = 0.0+eps;  SS_ref_db.bounds_ref[3][1] = 1.0-eps;
    SS_ref_db.bounds_ref[4][0] = 0.0+eps;  SS_ref_db.bounds_ref[4][1] = 1.0-eps;
    SS_ref_db.bounds_ref[5][0] = 0.0+eps;  SS_ref_db.bounds_ref[5][1] = 1.0-eps;
    
	if (z_b.bulk_rock[8] == 0.){ 					
		SS_ref_db.z_em[4]          = 0.0;
        SS_ref_db.d_em[4]          = 1.0;
		SS_ref_db.bounds_ref[3][0] = 0.0; 
		SS_ref_db.bounds_ref[3][1] = 0.0;	
	}
    return SS_ref_db;
}

/**
   retrieve reference thermodynamic data for mp_fsp
*/
SS_ref G_SS_mp_fsp_function(SS_ref SS_ref_db, int EM_database, int len_ox, bulk_info z_b, double eps){
    
    int i, j;
    int n_em = SS_ref_db.n_em;
    
    char   *EM_tmp[] 		= {"ab","an","san"};
    for (int i = 0; i < SS_ref_db.n_em; i++){
        strcpy(SS_ref_db.EM_list[i],EM_tmp[i]);
    };
    int n_xeos = SS_ref_db.n_xeos;
    char   *CV_tmp[] 		= {"ca","k"};
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
        strcpy(SS_ref_db.CV_list[i],CV_tmp[i]);
    };
    SS_ref_db.W[0] = -0.04*SS_ref_db.P - 0.00935*SS_ref_db.T + 14.6;
    SS_ref_db.W[1] = 0.338*SS_ref_db.P - 0.00957*SS_ref_db.T + 24.1;
    SS_ref_db.W[2] = 48.5 - 0.13*SS_ref_db.P;
    
    SS_ref_db.v[0] = 0.674;
    SS_ref_db.v[1] = 0.55;
    SS_ref_db.v[2] = 1.0;
    
    
    em_data ab_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"ab", 
    										"equilibrium"	);
    
    em_data an_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"an", 
    										"equilibrium"	);
    
    em_data san_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"san", 
    										"equilibrium"	);
    
    SS_ref_db.gbase[0] 		= ab_eq.gb;
    SS_ref_db.gbase[1] 		= an_eq.gb;
    SS_ref_db.gbase[2] 		= san_eq.gb;
    
    SS_ref_db.ElShearMod[0] 	= ab_eq.ElShearMod;
    SS_ref_db.ElShearMod[1] 	= an_eq.ElShearMod;
    SS_ref_db.ElShearMod[2] 	= san_eq.ElShearMod;
    
    for (i = 0; i < len_ox; i++){
        SS_ref_db.Comp[0][i] 	= ab_eq.C[i];
        SS_ref_db.Comp[1][i] 	= an_eq.C[i];
        SS_ref_db.Comp[2][i] 	= san_eq.C[i];
    }
    
    for (i = 0; i < n_em; i++){
        SS_ref_db.z_em[i] = 1.0;
    };
    
    SS_ref_db.bounds_ref[0][0] = 0.0+eps;  SS_ref_db.bounds_ref[0][1] = 1.0-eps;
    SS_ref_db.bounds_ref[1][0] = 0.0+eps;  SS_ref_db.bounds_ref[1][1] = 1.0-eps;
    
    return SS_ref_db;
}

/**
   retrieve reference thermodynamic data for mp_sa
*/
SS_ref G_SS_mp_sa_function(SS_ref SS_ref_db, int EM_database, int len_ox, bulk_info z_b, double eps){
    
    int i, j;
    int n_em = SS_ref_db.n_em;
    
    char   *EM_tmp[] 		= {"spr4","spr5","fspm","spro","ospr"};
    for (int i = 0; i < SS_ref_db.n_em; i++){
        strcpy(SS_ref_db.EM_list[i],EM_tmp[i]);
    };
    int n_xeos = SS_ref_db.n_xeos;
    char   *CV_tmp[] 		= {"x","y","f","Q"};
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
        strcpy(SS_ref_db.CV_list[i],CV_tmp[i]);
    };
    SS_ref_db.W[0] = 10. - 0.02*SS_ref_db.P;
    SS_ref_db.W[1] = 16.;
    SS_ref_db.W[2] = 12.;
    SS_ref_db.W[3] = 8. - 0.02*SS_ref_db.P;
    SS_ref_db.W[4] = 19. - 0.02*SS_ref_db.P;
    SS_ref_db.W[5] = 22. - 0.02*SS_ref_db.P;
    SS_ref_db.W[6] = 1.;
    SS_ref_db.W[7] = 4.;
    SS_ref_db.W[8] = 17.6 - 0.02*SS_ref_db.P;
    SS_ref_db.W[9] = 20. - 0.02*SS_ref_db.P;
    
    
    em_data spr4_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"spr4", 
    										"equilibrium"	);
    
    em_data spr5_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"spr5", 
    										"equilibrium"	);
    
    em_data fspr_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"fspr", 
    										"equilibrium"	);
    
    em_data gr_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"gr", 
    										"equilibrium"	);
    
    em_data andr_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"andr", 
    										"equilibrium"	);
    
    SS_ref_db.gbase[0] 		= spr4_eq.gb;
    SS_ref_db.gbase[1] 		= spr5_eq.gb;
    SS_ref_db.gbase[2] 		= fspr_eq.gb - 2.0;
    SS_ref_db.gbase[3] 		= 0.75*fspr_eq.gb + 0.25*spr4_eq.gb - 3.5;
    SS_ref_db.gbase[4] 		= 0.5*andr_eq.gb - 0.5*gr_eq.gb + spr5_eq.gb - 16.0;
    
    SS_ref_db.ElShearMod[0] 	= spr4_eq.ElShearMod;
    SS_ref_db.ElShearMod[1] 	= spr5_eq.ElShearMod;
    SS_ref_db.ElShearMod[2] 	= fspr_eq.ElShearMod;
    SS_ref_db.ElShearMod[3] 	= 0.75*fspr_eq.ElShearMod + 0.25*spr4_eq.ElShearMod;
    SS_ref_db.ElShearMod[4] 	= 0.5*andr_eq.ElShearMod - 0.5*gr_eq.ElShearMod + spr5_eq.ElShearMod;
    
    for (i = 0; i < len_ox; i++){
        SS_ref_db.Comp[0][i] 	= spr4_eq.C[i];
        SS_ref_db.Comp[1][i] 	= spr5_eq.C[i];
        SS_ref_db.Comp[2][i] 	= fspr_eq.C[i];
        SS_ref_db.Comp[3][i] 	= 0.75*fspr_eq.C[i] + 0.25*spr4_eq.C[i];
        SS_ref_db.Comp[4][i] 	= 0.5*andr_eq.C[i] - 0.5*gr_eq.C[i] + spr5_eq.C[i];
    }
    
    for (i = 0; i < n_em; i++){
        SS_ref_db.z_em[i] = 1.0;
    };
    
    SS_ref_db.bounds_ref[0][0] = 0.0+eps;  SS_ref_db.bounds_ref[0][1] = 1.0-eps;
    SS_ref_db.bounds_ref[1][0] = 0.0+eps;  SS_ref_db.bounds_ref[1][1] = 1.0-eps;
    SS_ref_db.bounds_ref[2][0] = 0.0+eps;  SS_ref_db.bounds_ref[2][1] = 1.0-eps;
    SS_ref_db.bounds_ref[3][0] = -1.0+eps;  SS_ref_db.bounds_ref[3][1] = 1.0-eps;
	if (z_b.bulk_rock[8] == 0.){ 					
		SS_ref_db.z_em[4]          = 0.0;
        SS_ref_db.d_em[4]          = 1.0;
		SS_ref_db.bounds_ref[2][0] = 0.0; 
		SS_ref_db.bounds_ref[2][1] = 0.0;	
	}
    return SS_ref_db;
}

/**
   retrieve reference thermodynamic data for mp_sp
*/
SS_ref G_SS_mp_sp_function(SS_ref SS_ref_db, int EM_database, int len_ox, bulk_info z_b, double eps){
    
    int i, j;
    int n_em = SS_ref_db.n_em;
    
    char   *EM_tmp[] 		= {"herc","sp","mt","usp"};
    for (int i = 0; i < SS_ref_db.n_em; i++){
        strcpy(SS_ref_db.EM_list[i],EM_tmp[i]);
    };
    int n_xeos = SS_ref_db.n_xeos;
    char   *CV_tmp[] 		= {"x","y","z"};
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
        strcpy(SS_ref_db.CV_list[i],CV_tmp[i]);
    };
    SS_ref_db.W[0] = 16.;
    SS_ref_db.W[1] = 2.;
    SS_ref_db.W[2] = 20.;
    SS_ref_db.W[3] = 18.;
    SS_ref_db.W[4] = 36.;
    SS_ref_db.W[5] = 30.;
    
    
    em_data herc_eq 	= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"herc", 
    										"equilibrium"	);
    
    em_data sp_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"sp", 
    										"equilibrium"	);
    
    em_data mt_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"mt", 
    										"equilibrium"	);
    
    em_data usp_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"usp", 
    										"equilibrium"	);
    
    SS_ref_db.gbase[0] 		= herc_eq.gb;
    SS_ref_db.gbase[1] 		= sp_eq.gb;
    SS_ref_db.gbase[2] 		= mt_eq.gb;
    SS_ref_db.gbase[3] 		= usp_eq.gb;
    
    SS_ref_db.ElShearMod[0] 	= herc_eq.ElShearMod;
    SS_ref_db.ElShearMod[1] 	= sp_eq.ElShearMod;
    SS_ref_db.ElShearMod[2] 	= mt_eq.ElShearMod;
    SS_ref_db.ElShearMod[3] 	= usp_eq.ElShearMod;
    
    for (i = 0; i < len_ox; i++){
        SS_ref_db.Comp[0][i] 	= herc_eq.C[i];
        SS_ref_db.Comp[1][i] 	= sp_eq.C[i];
        SS_ref_db.Comp[2][i] 	= mt_eq.C[i];
        SS_ref_db.Comp[3][i] 	= usp_eq.C[i];
    }
    
    for (i = 0; i < n_em; i++){
        SS_ref_db.z_em[i] = 1.0;
    };
    
    SS_ref_db.bounds_ref[0][0] = 0.0+eps;  SS_ref_db.bounds_ref[0][1] = 1.0-eps;
    SS_ref_db.bounds_ref[1][0] = 0.0+eps;  SS_ref_db.bounds_ref[1][1] = 1.0-eps;
    SS_ref_db.bounds_ref[2][0] = 0.0+eps;  SS_ref_db.bounds_ref[2][1] = 1.0-eps;

	if (z_b.bulk_rock[8] == 0.){ 					
		SS_ref_db.z_em[2]          = 0.0;
        SS_ref_db.d_em[2]          = 1.0;
		SS_ref_db.bounds_ref[1][0] = 0.0; 
		SS_ref_db.bounds_ref[1][1] = 0.0;	
	}
    return SS_ref_db;
}

/**
   retrieve reference thermodynamic data for mp_st
*/
SS_ref G_SS_mp_st_function(SS_ref SS_ref_db, int EM_database, int len_ox, bulk_info z_b, double eps){
    
    int i, j;
    int n_em = SS_ref_db.n_em;
    
    char   *EM_tmp[] 		= {"mstm","fst","mnstm","msto","mstt"};
    for (int i = 0; i < SS_ref_db.n_em; i++){
        strcpy(SS_ref_db.EM_list[i],EM_tmp[i]);
    };
    int n_xeos = SS_ref_db.n_xeos;
    char   *CV_tmp[] 		= {"x","m","f","t"};
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
        strcpy(SS_ref_db.CV_list[i],CV_tmp[i]);
    };
    SS_ref_db.W[0] = 16.;
    SS_ref_db.W[1] = 12.;
    SS_ref_db.W[2] = 2.;
    SS_ref_db.W[3] = 20.;
    SS_ref_db.W[4] = 8.;
    SS_ref_db.W[5] = 18.;
    SS_ref_db.W[6] = 36.;
    SS_ref_db.W[7] = 14.;
    SS_ref_db.W[8] = 32.;
    SS_ref_db.W[9] = 30.;
    
    
    em_data mst_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"mst", 
    										"equilibrium"	);
    
    em_data fst_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"fst", 
    										"equilibrium"	);
    
    em_data mnst_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"mnst", 
    										"equilibrium"	);
    
    em_data andr_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"andr", 
    										"equilibrium"	);
    
    em_data gr_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"gr", 
    										"equilibrium"	);
    
    em_data cor_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"cor", 
    										"equilibrium"	);
    
    em_data ru_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"ru", 
    										"equilibrium"	);
    
    SS_ref_db.gbase[0] 		= mst_eq.gb - 8.0;
    SS_ref_db.gbase[1] 		= fst_eq.gb;
    SS_ref_db.gbase[2] 		= mnst_eq.gb - 0.19;
    SS_ref_db.gbase[3] 		= andr_eq.gb - gr_eq.gb + mst_eq.gb + 9.0;
    SS_ref_db.gbase[4] 		= -cor_eq.gb + mst_eq.gb + 1.5*ru_eq.gb + 13.0;
    
    SS_ref_db.ElShearMod[0] 	= mst_eq.ElShearMod;
    SS_ref_db.ElShearMod[1] 	= fst_eq.ElShearMod;
    SS_ref_db.ElShearMod[2] 	= mnst_eq.ElShearMod;
    SS_ref_db.ElShearMod[3] 	= andr_eq.ElShearMod - gr_eq.ElShearMod + mst_eq.ElShearMod;
    SS_ref_db.ElShearMod[4] 	= -cor_eq.ElShearMod + mst_eq.ElShearMod + 1.5*ru_eq.ElShearMod;
    
    for (i = 0; i < len_ox; i++){
        SS_ref_db.Comp[0][i] 	= mst_eq.C[i];
        SS_ref_db.Comp[1][i] 	= fst_eq.C[i];
        SS_ref_db.Comp[2][i] 	= mnst_eq.C[i];
        SS_ref_db.Comp[3][i] 	= andr_eq.C[i] - gr_eq.C[i] + mst_eq.C[i];
        SS_ref_db.Comp[4][i] 	= -cor_eq.C[i] + mst_eq.C[i] + 1.5*ru_eq.C[i];
    }
    
    for (i = 0; i < n_em; i++){
        SS_ref_db.z_em[i] = 1.0;
    };
    
    SS_ref_db.bounds_ref[0][0] = 0.0+eps;  SS_ref_db.bounds_ref[0][1] = 1.0-eps;
    SS_ref_db.bounds_ref[1][0] = 0.0+eps;  SS_ref_db.bounds_ref[1][1] = 1.0-eps;
    SS_ref_db.bounds_ref[2][0] = 0.0+eps;  SS_ref_db.bounds_ref[2][1] = 1.0-eps;
    SS_ref_db.bounds_ref[3][0] = 0.0+eps;  SS_ref_db.bounds_ref[3][1] = 1.0-eps;
	if (z_b.bulk_rock[8] == 0.){ 					
		SS_ref_db.z_em[3]          = 0.0;
        SS_ref_db.d_em[3]          = 1.0;
		SS_ref_db.bounds_ref[2][0] = 0.0; 
		SS_ref_db.bounds_ref[2][1] = 0.0;	
	}
    return SS_ref_db;
}


/**
   retrieve reference thermodynamic data for igp_fper_S11
*/
SS_ref G_SS_ig_fper_function(SS_ref SS_ref_db, int EM_database, int len_ox, bulk_info z_b, double eps){
    
    int i, j;
    int n_em = SS_ref_db.n_em;
    
    char   *EM_tmp[] 		= {"per","wu"};
    for (int i = 0; i < SS_ref_db.n_em; i++){
        strcpy(SS_ref_db.EM_list[i],EM_tmp[i]);
    };
    
    SS_ref_db.W[0] = 13.0;
    
    
    em_data per_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"per", 
    										"equilibrium"	);
    
    em_data wu_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"wu", 
    										"equilibrium"	);
    
    SS_ref_db.gbase[0] 		= per_eq.gb; //tentative correction to fix Stiruxe data
    SS_ref_db.gbase[1] 		= wu_eq.gb;
    
    SS_ref_db.ElShearMod[0] 	= per_eq.ElShearMod;
    SS_ref_db.ElShearMod[1] 	= wu_eq.ElShearMod;
    
    for (i = 0; i < len_ox; i++){
        SS_ref_db.Comp[0][i] 	= per_eq.C[i];
        SS_ref_db.Comp[1][i] 	= wu_eq.C[i];
    }
    
    for (i = 0; i < n_em; i++){
        SS_ref_db.z_em[i] = 1.0;
    };
    
    SS_ref_db.bounds_ref[0][0] = 0.0+eps;  SS_ref_db.bounds_ref[0][1] = 1.0-eps;
    
    return SS_ref_db;
}

/**
  retrieve reference thermodynamic data for biotite 
*/
SS_ref G_SS_ig_bi_function(SS_ref SS_ref_db, int EM_database, int len_ox, bulk_info  z_b, double eps){	
   
    int i, j;
    int n_em = SS_ref_db.n_em;
    
    char   *EM_tmp[] 		= {"phl","annm","obi","eas","tbi","fbi"};
    for (int i = 0; i < SS_ref_db.n_em; i++){
        strcpy(SS_ref_db.EM_list[i],EM_tmp[i]);
    };
    int n_xeos = SS_ref_db.n_xeos;
    char   *CV_tmp[] 		= {"x","y","f","t","Q"};
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
        strcpy(SS_ref_db.CV_list[i],CV_tmp[i]);
    };
    SS_ref_db.W[0]  = 12.;
    SS_ref_db.W[1]  = 4.;
    SS_ref_db.W[2]  = 10.;
    SS_ref_db.W[3]  = 30.;
    SS_ref_db.W[4]  = 8.;
    SS_ref_db.W[5]  = 8.;
    SS_ref_db.W[6]  = 5.;
    SS_ref_db.W[7]  = 32.;
    SS_ref_db.W[8]  = 13.6;
    SS_ref_db.W[9]  = 7.;
    SS_ref_db.W[10] = 24.;
    SS_ref_db.W[11] = 5.6;
    SS_ref_db.W[12] = 40.0;
    SS_ref_db.W[13] = 1.0;
    SS_ref_db.W[14] = 40.0;		
    
    
	em_data phl_eq 		= get_em_data(		EM_database, 
											len_ox,
											z_b,
											SS_ref_db.P,
											SS_ref_db.T,
											"phl", 
											"equilibrium"	);

	em_data ann_eq 		= get_em_data(		EM_database, 
											len_ox,
											z_b,
											SS_ref_db.P,
											SS_ref_db.T,
											"ann", 
											"equilibrium"	);

	em_data east_eq 	= get_em_data(		EM_database, 
											len_ox,
											z_b,
											SS_ref_db.P,
											SS_ref_db.T,
											"east", 
											"equilibrium"	);

	em_data br_eq 		= get_em_data(		EM_database, 
											len_ox,
											z_b,
											SS_ref_db.P,
											SS_ref_db.T,
											"br", 
											"equilibrium"	);

	em_data ru_eq 		= get_em_data(		EM_database, 
											len_ox,
											z_b,
											SS_ref_db.P,
											SS_ref_db.T,
											"ru", 
											"equilibrium"	);

	em_data gr_eq 		= get_em_data(		EM_database, 
											len_ox,
											z_b,
											SS_ref_db.P,
											SS_ref_db.T,
											"gr", 
											"equilibrium"	);

	em_data andr_eq 	= get_em_data(		EM_database, 
											len_ox,
											z_b,
											SS_ref_db.P,
											SS_ref_db.T,
											"andr", 
											"equilibrium"	);
    

	SS_ref_db.gbase[0] 		= phl_eq.gb;
	SS_ref_db.gbase[1] 		= ann_eq.gb - 6.0;
	SS_ref_db.gbase[2] 		= 1./3. * ann_eq.gb + 2./3. * phl_eq.gb - 6.0;
	SS_ref_db.gbase[3] 		= east_eq.gb;
	SS_ref_db.gbase[4] 		= - br_eq.gb + phl_eq.gb + ru_eq.gb + 55.0;
	SS_ref_db.gbase[5] 		= 1./2. * andr_eq.gb + east_eq.gb - 1./2. * gr_eq.gb - 3.0;

	SS_ref_db.ElShearMod[0] = phl_eq.ElShearMod;
	SS_ref_db.ElShearMod[1] = ann_eq.ElShearMod;
	SS_ref_db.ElShearMod[2] = 1./3. * ann_eq.ElShearMod + 2./3. * phl_eq.ElShearMod;
	SS_ref_db.ElShearMod[3] = east_eq.ElShearMod;
	SS_ref_db.ElShearMod[4] = - br_eq.ElShearMod + phl_eq.ElShearMod + ru_eq.ElShearMod;
	SS_ref_db.ElShearMod[5] = 1./2. * andr_eq.ElShearMod + east_eq.ElShearMod - 1./2. * gr_eq.ElShearMod;

	for (i = 0; i < len_ox; i++){
		SS_ref_db.Comp[0][i] = phl_eq.C[i];
		SS_ref_db.Comp[1][i] = ann_eq.C[i];
		SS_ref_db.Comp[2][i] = 1./3. * ann_eq.C[i] + 2./3. * phl_eq.C[i] ;
		SS_ref_db.Comp[3][i] = east_eq.C[i];
		SS_ref_db.Comp[4][i] = - br_eq.C[i] + phl_eq.C[i] + ru_eq.C[i] ;
		SS_ref_db.Comp[5][i] = 1./2. * andr_eq.C[i] + east_eq.C[i] - 1./2. * gr_eq.C[i];
	}
	
    for (i = 0; i < n_em; i++){
        SS_ref_db.z_em[i] = 1.0;
    };
    
    SS_ref_db.bounds_ref[0][0] = 0.0+eps;  SS_ref_db.bounds_ref[0][1] = 1.0-eps;
    SS_ref_db.bounds_ref[1][0] = 0.0+eps;  SS_ref_db.bounds_ref[1][1] = 1.0-eps;
    SS_ref_db.bounds_ref[2][0] = 0.0+eps;  SS_ref_db.bounds_ref[2][1] = 1.0-eps;
    SS_ref_db.bounds_ref[3][0] = 0.0+eps;  SS_ref_db.bounds_ref[3][1] = 1.0-eps;
	SS_ref_db.bounds_ref[4][0] = -1.0+eps;	SS_ref_db.bounds_ref[4][1] = 1.0-eps;	

	if (z_b.bulk_rock[8] == 0.){ 					
		SS_ref_db.z_em[5]          = 0.0;
        SS_ref_db.d_em[5]          = 1.0;
		SS_ref_db.bounds_ref[2][0] = 0.0; 
		SS_ref_db.bounds_ref[2][1] = 0.0;	
	}
	if (z_b.bulk_rock[7] == 0.){ 					
		SS_ref_db.z_em[4]          = 0.0;
        SS_ref_db.d_em[4]          = 1.0;
		SS_ref_db.bounds_ref[3][0] = 0.0; 
		SS_ref_db.bounds_ref[3][1] = 0.0;	
	}


    /* this lists the index of the order variables */
    SS_ref_db.orderVar      = 1;
    SS_ref_db.idOrderVar[4] = -1.0;

    return SS_ref_db;
}

/**
  retrieve reference thermodynamic data for clinopyroxene 
*/
SS_ref G_SS_ig_cpx_function(SS_ref SS_ref_db, int EM_database, int len_ox, bulk_info z_b, double eps){
    
    int i, j;
    int n_em = SS_ref_db.n_em;
    
    char   *EM_tmp[] 		= {"di","cfs","cats","crdi","cess","cbuf","jd","cen","cfm","kjd"};
    for (int i = 0; i < SS_ref_db.n_em; i++){
        strcpy(SS_ref_db.EM_list[i],EM_tmp[i]);
    };
    int n_xeos = SS_ref_db.n_xeos;
    char   *CV_tmp[] 		= {"x","y","o","n","Q","f","cr","t","k"};
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
        strcpy(SS_ref_db.CV_list[i],CV_tmp[i]);
    };
  SS_ref_db.W[0] = 25.8;
    SS_ref_db.W[1] = 13.0 - 0.06*SS_ref_db.P;
    SS_ref_db.W[2] = 8.0;
    SS_ref_db.W[3] = 8.0;
    SS_ref_db.W[4] = 8.0;
    SS_ref_db.W[5] = 26.0;
    SS_ref_db.W[6] = 29.8;
    SS_ref_db.W[7] = 20.6;
    SS_ref_db.W[8] = 26.0;
    SS_ref_db.W[9] = 25.0 - 0.1*SS_ref_db.P;
    SS_ref_db.W[10] = 38.3;
    SS_ref_db.W[11] = 43.3;
    SS_ref_db.W[12] = 24.0;
    SS_ref_db.W[13] = 24.0;
    SS_ref_db.W[14] = 2.3;
    SS_ref_db.W[15] = 3.5;
    SS_ref_db.W[16] = 24.0;
    SS_ref_db.W[17] = 2.0;
    SS_ref_db.W[18] = 2.0;
    SS_ref_db.W[19] = 6.0;
    SS_ref_db.W[20] = 6.0;
    SS_ref_db.W[21] = 45.2 - 0.35*SS_ref_db.P;
    SS_ref_db.W[22] = 27.0 - 0.1*SS_ref_db.P;
    SS_ref_db.W[23] = 6.0;
    SS_ref_db.W[24] = 2.0;
    SS_ref_db.W[25] = 6.0;
    SS_ref_db.W[26] = 3.0;
    SS_ref_db.W[27] = 52.3;
    SS_ref_db.W[28] = 40.3;
    SS_ref_db.W[29] = 3.0;
    SS_ref_db.W[30] = 6.0;
    SS_ref_db.W[31] = 3.0;
    SS_ref_db.W[32] = 57.3;
    SS_ref_db.W[33] = 45.3;
    SS_ref_db.W[34] = 3.0;
    SS_ref_db.W[35] = 16.0;
    SS_ref_db.W[36] = 24.0;
    SS_ref_db.W[37] = 22.0;
    SS_ref_db.W[38] = 16.0;
    SS_ref_db.W[39] = 40.0;
    SS_ref_db.W[40] = 40.0;
    SS_ref_db.W[41] = 28.0;
    SS_ref_db.W[42] = 4.0;
    SS_ref_db.W[43] = 40.0;
    SS_ref_db.W[44] = 40.0;
    
    SS_ref_db.v[0] = 1.2;
    SS_ref_db.v[1] = 1.0;
    SS_ref_db.v[2] = 1.9;
    SS_ref_db.v[3] = 1.9;
    SS_ref_db.v[4] = 1.9;
    SS_ref_db.v[5] = 1.9;
    SS_ref_db.v[6] = 1.2;
    SS_ref_db.v[7] = 1.0;
    SS_ref_db.v[8] = 1.0;
    SS_ref_db.v[9] = 1.2;
  
    em_data di_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
											SS_ref_db.P,
											SS_ref_db.T,
    										"di", 
    										"equilibrium"	);
    
    em_data fs_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
											SS_ref_db.P,
											SS_ref_db.T,
    										"fs", 
    										"equilibrium"	);
    
    em_data cats_eq 	= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
											SS_ref_db.P,
											SS_ref_db.T,
    										"cats", 
    										"equilibrium"	);
    
    em_data kos_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
											SS_ref_db.P,
											SS_ref_db.T,
    										"kos", 
    										"equilibrium"	);
    
    em_data jd_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
											SS_ref_db.P,
											SS_ref_db.T,
    										"jd", 
    										"equilibrium"	);
    
    em_data cats_di 	= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
											SS_ref_db.P,
											SS_ref_db.T,
    										"cats", 
    										"disordered"	);
    
    em_data acm_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
											SS_ref_db.P,
											SS_ref_db.T,
    										"acm", 
    										"equilibrium"	);
    
    em_data ru_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
											SS_ref_db.P,
											SS_ref_db.T,
    										"ru", 
    										"equilibrium"	);
    
    em_data cor_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
											SS_ref_db.P,
											SS_ref_db.T,
    										"cor", 
    										"equilibrium"	);
    
    em_data per_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
											SS_ref_db.P,
											SS_ref_db.T,
    										"per", 
    										"equilibrium"	);
    
    em_data en_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
											SS_ref_db.P,
											SS_ref_db.T,
    										"en", 
    										"equilibrium"	);
    
    em_data abh_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
											SS_ref_db.P,
											SS_ref_db.T,
    										"abh", 
    										"equilibrium"	);
    
    em_data san_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
											SS_ref_db.P,
											SS_ref_db.T,
    										"san", 
    										"equilibrium"	);
    
    SS_ref_db.gbase[0] 		= di_eq.gb;
    SS_ref_db.gbase[1] 		= 0.045*SS_ref_db.P - 0.002*SS_ref_db.T + fs_eq.gb + 2.1;
    SS_ref_db.gbase[2] 		= cats_eq.gb;
    SS_ref_db.gbase[3] 		= cats_di.gb - jd_eq.gb + kos_eq.gb - 4.9;
    SS_ref_db.gbase[4] 		= acm_eq.gb + cats_di.gb - jd_eq.gb - 3.45;
    SS_ref_db.gbase[5] 		= -0.005*SS_ref_db.P - 0.0012*SS_ref_db.T - 0.5*cor_eq.gb + cats_di.gb + 0.5*per_eq.gb + 0.5*ru_eq.gb - 16.2;
    SS_ref_db.gbase[6] 		= jd_eq.gb;
    SS_ref_db.gbase[7] 		= 0.048*SS_ref_db.P - 0.002*SS_ref_db.T + en_eq.gb + 3.5;
    SS_ref_db.gbase[8] 		= 0.0465*SS_ref_db.P - 0.002*SS_ref_db.T + 0.5*en_eq.gb + 0.5*fs_eq.gb - 1.6;
    SS_ref_db.gbase[9] 		= 0.6*SS_ref_db.P - abh_eq.gb + jd_eq.gb + san_eq.gb + 11.7;
    
    SS_ref_db.ElShearMod[0] 	= di_eq.ElShearMod;
    SS_ref_db.ElShearMod[1] 	= fs_eq.ElShearMod;
    SS_ref_db.ElShearMod[2] 	= cats_eq.ElShearMod;
    SS_ref_db.ElShearMod[3] 	= cats_di.ElShearMod - jd_eq.ElShearMod + kos_eq.ElShearMod;
    SS_ref_db.ElShearMod[4] 	= acm_eq.ElShearMod + cats_di.ElShearMod - jd_eq.ElShearMod;
    SS_ref_db.ElShearMod[5] 	= -0.5*cor_eq.ElShearMod + cats_di.ElShearMod + 0.5*per_eq.ElShearMod + 0.5*ru_eq.ElShearMod;
    SS_ref_db.ElShearMod[6] 	= jd_eq.ElShearMod;
    SS_ref_db.ElShearMod[7] 	= en_eq.ElShearMod;
    SS_ref_db.ElShearMod[8] 	= 0.5*en_eq.ElShearMod + 0.5*fs_eq.ElShearMod;
    SS_ref_db.ElShearMod[9] 	= -abh_eq.ElShearMod + jd_eq.ElShearMod + san_eq.ElShearMod;
    
    for (i = 0; i < len_ox; i++){
        SS_ref_db.Comp[0][i] 	= di_eq.C[i];
        SS_ref_db.Comp[1][i] 	= fs_eq.C[i];
        SS_ref_db.Comp[2][i] 	= cats_eq.C[i];
        SS_ref_db.Comp[3][i] 	= cats_di.C[i] - jd_eq.C[i] + kos_eq.C[i];
        SS_ref_db.Comp[4][i] 	= acm_eq.C[i] + cats_di.C[i] - jd_eq.C[i];
        SS_ref_db.Comp[5][i] 	= -0.5*cor_eq.C[i] + cats_di.C[i] + 0.5*per_eq.C[i] + 0.5*ru_eq.C[i];
        SS_ref_db.Comp[6][i] 	= jd_eq.C[i];
        SS_ref_db.Comp[7][i] 	= en_eq.C[i];
        SS_ref_db.Comp[8][i] 	= 0.5*en_eq.C[i] + 0.5*fs_eq.C[i];
        SS_ref_db.Comp[9][i] 	= -abh_eq.C[i] + jd_eq.C[i] + san_eq.C[i];
    }
    
    for (i = 0; i < n_em; i++){
        SS_ref_db.z_em[i] = 1.0;
    };
    
    
    SS_ref_db.bounds_ref[0][0] = 0.0+eps;  SS_ref_db.bounds_ref[0][1] = 1.0-eps;
    SS_ref_db.bounds_ref[1][0] = 0.0+eps;  SS_ref_db.bounds_ref[1][1] = 2.0-eps;
    SS_ref_db.bounds_ref[2][0] = 0.0+eps;  SS_ref_db.bounds_ref[2][1] = 1.0-eps;
    SS_ref_db.bounds_ref[3][0] = 0.0+eps;  SS_ref_db.bounds_ref[3][1] = 1.0-eps;
    SS_ref_db.bounds_ref[4][0] = -1.0+eps;  SS_ref_db.bounds_ref[4][1] = 1.0-eps;
    SS_ref_db.bounds_ref[5][0] = 0.0+eps;  SS_ref_db.bounds_ref[5][1] = 1.0-eps;
    SS_ref_db.bounds_ref[6][0] = 0.0+eps;  SS_ref_db.bounds_ref[6][1] = 1.0-eps;
    SS_ref_db.bounds_ref[7][0] = 0.0+eps;  SS_ref_db.bounds_ref[7][1] = 1.0-eps;
    SS_ref_db.bounds_ref[8][0] = 0.0+eps;  SS_ref_db.bounds_ref[8][1] = 1.0-eps;

	if (z_b.bulk_rock[9] == 0.){ 					
		SS_ref_db.z_em[3]          = 0.0;
        SS_ref_db.d_em[3]          = 1.0;
		SS_ref_db.bounds_ref[6][0] = 0.0; 
		SS_ref_db.bounds_ref[6][1] = 0.0;	
	}

	if (z_b.bulk_rock[8] == 0.){ 					
		SS_ref_db.z_em[4]          = 0.0;
        SS_ref_db.d_em[4]          = 1.0;
		SS_ref_db.bounds_ref[5][0] = 0.0; 
		SS_ref_db.bounds_ref[5][1] = 0.0;	
	}
  	if (z_b.bulk_rock[7] == 0.){ 					
		SS_ref_db.z_em[5]          = 0.0;
        SS_ref_db.d_em[5]          = 1.0;
		SS_ref_db.bounds_ref[7][0] = 0.0; 
		SS_ref_db.bounds_ref[7][1] = 0.0;	
	}  
    return SS_ref_db;
}



/**
  retrieve reference thermodynamic data for cordierite 
*/
SS_ref G_SS_ig_cd_function(SS_ref SS_ref_db, int EM_database, int len_ox, bulk_info z_b, double eps){
    
    int i, j;
    int n_em = SS_ref_db.n_em;
    
    char   *EM_tmp[] 		= {"crd","fcrd","hcrd"};
    for (int i = 0; i < SS_ref_db.n_em; i++){
        strcpy(SS_ref_db.EM_list[i],EM_tmp[i]);
    };
    int n_xeos = SS_ref_db.n_xeos;
    char   *CV_tmp[] 		= {"x","h"};
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
        strcpy(SS_ref_db.CV_list[i],CV_tmp[i]);
    };

    SS_ref_db.W[0] = 6.0;
    SS_ref_db.W[1] = 0.0;
    SS_ref_db.W[2] = 0.0;
    
    
    em_data crd_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
											SS_ref_db.P,
											SS_ref_db.T,
    										"crd", 
    										"equilibrium"	);
    
    em_data fcrd_eq 	= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
											SS_ref_db.P,
											SS_ref_db.T,
    										"fcrd", 
    										"equilibrium"	);
    
    em_data hcrd_eq 	= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
											SS_ref_db.P,
											SS_ref_db.T,
    										"hcrd", 
    										"equilibrium"	);
    
    SS_ref_db.gbase[0] 		= crd_eq.gb;
    SS_ref_db.gbase[1] 		= fcrd_eq.gb;
    SS_ref_db.gbase[2] 		= hcrd_eq.gb;
    
    SS_ref_db.ElShearMod[0] 	= crd_eq.ElShearMod;
    SS_ref_db.ElShearMod[1] 	= fcrd_eq.ElShearMod;
    SS_ref_db.ElShearMod[2] 	= hcrd_eq.ElShearMod;
    
    for (i = 0; i < len_ox; i++){
        SS_ref_db.Comp[0][i] 	= crd_eq.C[i];
        SS_ref_db.Comp[1][i] 	= fcrd_eq.C[i];
        SS_ref_db.Comp[2][i] 	= hcrd_eq.C[i];
    }
    
    
    for (i = 0; i < n_em; i++){
        SS_ref_db.z_em[i] = 1.0;
    };
    
    SS_ref_db.bounds_ref[0][0] = 0.0+eps;  SS_ref_db.bounds_ref[0][1] = 1.0-eps;
    SS_ref_db.bounds_ref[1][0] = 0.0+eps;  SS_ref_db.bounds_ref[1][1] = 1.0-eps;

    return SS_ref_db;
}



/**
  retrieve reference thermodynamic data for epidote 
*/
SS_ref G_SS_ig_ep_function(SS_ref SS_ref_db, int EM_database, int len_ox, bulk_info z_b, double eps){
    
    int i, j;
    int n_em = SS_ref_db.n_em;
    
    char   *EM_tmp[] 		= {"cz","ep","fep"};
    for (int i = 0; i < SS_ref_db.n_em; i++){
        strcpy(SS_ref_db.EM_list[i],EM_tmp[i]);
    };
    int n_xeos = SS_ref_db.n_xeos;
    char   *CV_tmp[] 		= {"f","Q"};
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
        strcpy(SS_ref_db.CV_list[i],CV_tmp[i]);
    };
    
    SS_ref_db.W[0] = 1.0;
    SS_ref_db.W[1] = 3.0;
    SS_ref_db.W[2] = 1.0;
    
    
    em_data cz_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
											SS_ref_db.P,
											SS_ref_db.T,
    										"cz", 
    										"equilibrium"	);
    
    em_data ep_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
											SS_ref_db.P,
											SS_ref_db.T,
    										"ep", 
    										"equilibrium"	);
    
    em_data fep_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
											SS_ref_db.P,
											SS_ref_db.T,
    										"fep", 
    										"equilibrium"	);
    
    SS_ref_db.gbase[0] 		= cz_eq.gb;
    SS_ref_db.gbase[1] 		= ep_eq.gb;
    SS_ref_db.gbase[2] 		= fep_eq.gb;
    
    SS_ref_db.ElShearMod[0] 	= cz_eq.ElShearMod;
    SS_ref_db.ElShearMod[1] 	= ep_eq.ElShearMod;
    SS_ref_db.ElShearMod[2] 	= fep_eq.ElShearMod;
    
    for (i = 0; i < len_ox; i++){
        SS_ref_db.Comp[0][i] 	= cz_eq.C[i];
        SS_ref_db.Comp[1][i] 	= ep_eq.C[i];
        SS_ref_db.Comp[2][i] 	= fep_eq.C[i];
    }
    
    for (i = 0; i < n_em; i++){
        SS_ref_db.z_em[i] = 1.0;
    };
    
    SS_ref_db.bounds_ref[0][0] =  0.0+eps;  SS_ref_db.bounds_ref[0][1] = 1.0-eps;
    SS_ref_db.bounds_ref[1][0] = -0.5+eps;  SS_ref_db.bounds_ref[1][1] = 0.5-eps;

	if (z_b.bulk_rock[8] == 0.){ 					
		SS_ref_db.bounds_ref[0][0] = 0.0; 
		SS_ref_db.bounds_ref[0][1] = 0.0;	
		SS_ref_db.bounds_ref[1][0] = 0.0; 
		SS_ref_db.bounds_ref[1][1] = 0.0;	
		SS_ref_db.z_em[1]          = 0.0;
        SS_ref_db.d_em[1]          = 1.0;
		SS_ref_db.z_em[2]          = 0.0;
        SS_ref_db.d_em[2]          = 1.0;
	}
    /* this lists the index of the order variables */
    SS_ref_db.orderVar      = 1;
    SS_ref_db.idOrderVar[1] = -1.0;

    return SS_ref_db;
}

/**
   retrieve reference thermodynamic data for ig_flH
*/
SS_ref G_SS_ig_fl_function(SS_ref SS_ref_db, int EM_database, int len_ox, bulk_info z_b, double eps){
    
    int i, j;
    int n_em = SS_ref_db.n_em;
    
    char   *EM_tmp[] 		= {"qfL","slfL","wofL","fofL","fafL","jdfL","hmfL","ekfL","tifL","kjfL","H2O"};
    for (int i = 0; i < SS_ref_db.n_em; i++){
        strcpy(SS_ref_db.EM_list[i],EM_tmp[i]);
    };
    int n_xeos = SS_ref_db.n_xeos;
    char   *CV_tmp[] 		= {"wo","sl","fo","fa","jd","hm","ek","ti","kj","h2o"};
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
        strcpy(SS_ref_db.CV_list[i],CV_tmp[i]);
    };
    
    SS_ref_db.W[0] = 0.;
    SS_ref_db.W[1] = 0.;
    SS_ref_db.W[2] = 0.;
    SS_ref_db.W[3] = 0.;
    SS_ref_db.W[4] = 0.;
    SS_ref_db.W[5] = 0.;
    SS_ref_db.W[6] = 0.;
    SS_ref_db.W[7] = 0.;
    SS_ref_db.W[8] = 0.;
    SS_ref_db.W[9] = 59.0 - 0.82*SS_ref_db.P;
    SS_ref_db.W[10] = 0.;
    SS_ref_db.W[11] = 0.;
    SS_ref_db.W[12] = 0.;
    SS_ref_db.W[13] = 0.;
    SS_ref_db.W[14] = 0.;
    SS_ref_db.W[15] = 0.;
    SS_ref_db.W[16] = 0.;
    SS_ref_db.W[17] = 0.;
    SS_ref_db.W[18] = 57.6 - 0.8*SS_ref_db.P;
    SS_ref_db.W[19] = 0.;
    SS_ref_db.W[20] = 0.;
    SS_ref_db.W[21] = 0.;
    SS_ref_db.W[22] = 0.;
    SS_ref_db.W[23] = 0.;
    SS_ref_db.W[24] = 0.;
    SS_ref_db.W[25] = 0.;
    SS_ref_db.W[26] = 72.2 - 0.67*SS_ref_db.P;
    SS_ref_db.W[27] = 0.;
    SS_ref_db.W[28] = 0.;
    SS_ref_db.W[29] = 0.;
    SS_ref_db.W[30] = 0.;
    SS_ref_db.W[31] = 0.;
    SS_ref_db.W[32] = 0.;
    SS_ref_db.W[33] = 71.7 - 1.1*SS_ref_db.P;
    SS_ref_db.W[34] = 0.;
    SS_ref_db.W[35] = 0.;
    SS_ref_db.W[36] = 0.;
    SS_ref_db.W[37] = 0.;
    SS_ref_db.W[38] = 0.;
    SS_ref_db.W[39] = 71.7 - 1.1*SS_ref_db.P;
    SS_ref_db.W[40] = 0.;
    SS_ref_db.W[41] = 0.;
    SS_ref_db.W[42] = 0.;
    SS_ref_db.W[43] = 0.;
    SS_ref_db.W[44] = 57.0 - 0.79*SS_ref_db.P;
    SS_ref_db.W[45] = 0.;
    SS_ref_db.W[46] = 0.;
    SS_ref_db.W[47] = 0.;
    SS_ref_db.W[48] = 73.0 - 0.66*SS_ref_db.P;
    SS_ref_db.W[49] = 0.;
    SS_ref_db.W[50] = 0.;
    SS_ref_db.W[51] = 73.0 - 0.66*SS_ref_db.P;
    SS_ref_db.W[52] = 0.;
    SS_ref_db.W[53] = 75.0 - 0.67*SS_ref_db.P;
    SS_ref_db.W[54] = 44.9 - 1.19*SS_ref_db.P;
    
    
    em_data qL_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
											SS_ref_db.P,
											SS_ref_db.T,
    										"qL", 
    										"equilibrium"	);
    
    em_data silL_eq 	= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
											SS_ref_db.P,
											SS_ref_db.T,
    										"silL", 
    										"equilibrium"	);
    
    em_data woL_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
											SS_ref_db.P,
											SS_ref_db.T,
    										"woL", 
    										"equilibrium"	);
    
    em_data foL_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
											SS_ref_db.P,
											SS_ref_db.T,
    										"foL", 
    										"equilibrium"	);
    
    em_data faL_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
											SS_ref_db.P,
											SS_ref_db.T,
    										"faL", 
    										"equilibrium"	);
    
    em_data abL_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
											SS_ref_db.P,
											SS_ref_db.T,
    										"abL", 
    										"equilibrium"	);
    
    em_data hemL_eq 	= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
											SS_ref_db.P,
											SS_ref_db.T,
    										"hemL", 
    										"equilibrium"	);
    
    em_data eskL_eq 	= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
											SS_ref_db.P,
											SS_ref_db.T,
    										"eskL", 
    										"equilibrium"	);
    
    em_data ruL_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
											SS_ref_db.P,
											SS_ref_db.T,
    										"ruL", 
    										"equilibrium"	);
    
    em_data kspL_eq 	= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
											SS_ref_db.P,
											SS_ref_db.T,
    										"kspL", 
    										"equilibrium"	);
    
    em_data H2O_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
											SS_ref_db.P,
											SS_ref_db.T,
    										"H2O", 
    										"equilibrium"	);
    
    SS_ref_db.gbase[0] 		= -0.051*SS_ref_db.P + 4.0*qL_eq.gb + 2.1;
    SS_ref_db.gbase[1] 		= -0.313*SS_ref_db.P + silL_eq.gb + 6.72;
    SS_ref_db.gbase[2] 		= -0.12*SS_ref_db.P + woL_eq.gb + 0.22;
    SS_ref_db.gbase[3] 		= -0.136*SS_ref_db.P + 2.0*foL_eq.gb + 8.59;
    SS_ref_db.gbase[4] 		= -0.052*SS_ref_db.P + 2.0*faL_eq.gb + 13.56;
    SS_ref_db.gbase[5] 		= -0.099*SS_ref_db.P + abL_eq.gb - qL_eq.gb + 12.32;
    SS_ref_db.gbase[6] 		= -0.077*SS_ref_db.P + 0.5*hemL_eq.gb + 4.05;
    SS_ref_db.gbase[7] 		=  0.245*SS_ref_db.P + 0.5*eskL_eq.gb + 24.75;
    SS_ref_db.gbase[8] 		= -0.489*SS_ref_db.P + ruL_eq.gb + 5.6;
    SS_ref_db.gbase[9] 		= -0.227*SS_ref_db.P + kspL_eq.gb - qL_eq.gb + 12.88;
    SS_ref_db.gbase[10] 		= H2O_eq.gb;
    
    SS_ref_db.ElShearMod[0] 	= 4.0*qL_eq.ElShearMod;
    SS_ref_db.ElShearMod[1] 	= silL_eq.ElShearMod;
    SS_ref_db.ElShearMod[2] 	= woL_eq.ElShearMod;
    SS_ref_db.ElShearMod[3] 	= 2.0*foL_eq.ElShearMod;
    SS_ref_db.ElShearMod[4] 	= 2.0*faL_eq.ElShearMod;
    SS_ref_db.ElShearMod[5] 	= abL_eq.ElShearMod - qL_eq.ElShearMod;
    SS_ref_db.ElShearMod[6] 	= 0.5*hemL_eq.ElShearMod;
    SS_ref_db.ElShearMod[7] 	= 0.5*eskL_eq.ElShearMod;
    SS_ref_db.ElShearMod[8] 	= ruL_eq.ElShearMod;
    SS_ref_db.ElShearMod[9] 	= kspL_eq.ElShearMod - qL_eq.ElShearMod;
    SS_ref_db.ElShearMod[10] 	= H2O_eq.ElShearMod;
    
    for (i = 0; i < len_ox; i++){
        SS_ref_db.Comp[0][i] 	= 4.0*qL_eq.C[i];
        SS_ref_db.Comp[1][i] 	= silL_eq.C[i];
        SS_ref_db.Comp[2][i] 	= woL_eq.C[i];
        SS_ref_db.Comp[3][i] 	= 2.0*foL_eq.C[i];
        SS_ref_db.Comp[4][i] 	= 2.0*faL_eq.C[i];
        SS_ref_db.Comp[5][i] 	= abL_eq.C[i] - qL_eq.C[i];
        SS_ref_db.Comp[6][i] 	= 0.5*hemL_eq.C[i];
        SS_ref_db.Comp[7][i] 	= 0.5*eskL_eq.C[i];
        SS_ref_db.Comp[8][i] 	= ruL_eq.C[i];
        SS_ref_db.Comp[9][i] 	= kspL_eq.C[i] - qL_eq.C[i];
        SS_ref_db.Comp[10][i] 	= H2O_eq.C[i];
    }
    
    
    for (i = 0; i < n_em; i++){
        SS_ref_db.z_em[i] = 1.0;
    };
    
    SS_ref_db.bounds_ref[0][0] = 0.0+eps;  SS_ref_db.bounds_ref[0][1] = 1.0-eps;
    SS_ref_db.bounds_ref[1][0] = 0.0+eps;  SS_ref_db.bounds_ref[1][1] = 1.0-eps;
    SS_ref_db.bounds_ref[2][0] = 0.0+eps;  SS_ref_db.bounds_ref[2][1] = 1.0-eps;
    SS_ref_db.bounds_ref[3][0] = 0.0+eps;  SS_ref_db.bounds_ref[3][1] = 1.0-eps;
    SS_ref_db.bounds_ref[4][0] = 0.0+eps;  SS_ref_db.bounds_ref[4][1] = 1.0-eps;
    SS_ref_db.bounds_ref[5][0] = 0.0+eps;  SS_ref_db.bounds_ref[5][1] = 1.0-eps;
    SS_ref_db.bounds_ref[6][0] = 0.0+eps;  SS_ref_db.bounds_ref[6][1] = 1.0-eps;
    SS_ref_db.bounds_ref[7][0] = 0.0+eps;  SS_ref_db.bounds_ref[7][1] = 1.0-eps;
    SS_ref_db.bounds_ref[8][0] = 0.0+eps;  SS_ref_db.bounds_ref[8][1] = 1.0-eps;
    SS_ref_db.bounds_ref[9][0] = 0.0+eps;  SS_ref_db.bounds_ref[9][1] = 1.0-eps;
    
	if (z_b.bulk_rock[10] == 0.){
		SS_ref_db.z_em[10]         = 0.0;
		SS_ref_db.bounds_ref[9][0] = eps;  
		SS_ref_db.bounds_ref[9][1] = eps;	
	}
    
	if (z_b.bulk_rock[9] == 0.){ 					
		SS_ref_db.z_em[7]          = 0.0;
        SS_ref_db.d_em[7]          = 1.0;
		SS_ref_db.bounds_ref[6][0] = 0.0; 
		SS_ref_db.bounds_ref[6][1] = 0.0;	
	}
	if (z_b.bulk_rock[8] == 0.){ 					
		SS_ref_db.z_em[6]          = 0.0;
        SS_ref_db.d_em[6]          = 1.0;
		SS_ref_db.bounds_ref[5][0] = 0.0; 
		SS_ref_db.bounds_ref[5][1] = 0.0;	
	}
	if (z_b.bulk_rock[7] == 0.){ 					
		SS_ref_db.z_em[8]          = 0.0;
        SS_ref_db.d_em[8]          = 1.0;
		SS_ref_db.bounds_ref[7][0] = 0.0; 
		SS_ref_db.bounds_ref[7][1] = 0.0;	
	}

	return SS_ref_db;
}

/**
  retrieve reference thermodynamic data for garnet 
*/
SS_ref G_SS_ig_g_function(SS_ref SS_ref_db, int EM_database, int len_ox, bulk_info z_b, double eps){
    
    int i, j;
    int n_em = SS_ref_db.n_em;
    
    char   *EM_tmp[] 		= {"py","alm","gr","andr","knom","tig"};
    for (int i = 0; i < SS_ref_db.n_em; i++){
        strcpy(SS_ref_db.EM_list[i],EM_tmp[i]);
    };
    int n_xeos = SS_ref_db.n_xeos;
    char   *CV_tmp[] 		= {"x","c","f","cr","t"};
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
        strcpy(SS_ref_db.CV_list[i],CV_tmp[i]);
    };
    
    SS_ref_db.W[0] = 0.1*SS_ref_db.P + 4.0;
    SS_ref_db.W[1] = 0.04*SS_ref_db.P - 0.01*SS_ref_db.T + 45.4;
    SS_ref_db.W[2] = -0.036*SS_ref_db.P - 0.01*SS_ref_db.T + 107.0;
    SS_ref_db.W[3] = 2.0;
    SS_ref_db.W[4] = 0;
    SS_ref_db.W[5] = 0.1*SS_ref_db.P - 0.01*SS_ref_db.T + 17.0;
    SS_ref_db.W[6] = 0.039*SS_ref_db.P - 0.01*SS_ref_db.T + 65.0;
    SS_ref_db.W[7] = 0.01*SS_ref_db.P + 6.0;
    SS_ref_db.W[8] = 0;
    SS_ref_db.W[9] = 2.0;
    SS_ref_db.W[10] = 0.18*SS_ref_db.P - 0.01*SS_ref_db.T + 1.0;
    SS_ref_db.W[11] = 0;
    SS_ref_db.W[12] = 0.1*SS_ref_db.P - 0.01*SS_ref_db.T + 63.0;
    SS_ref_db.W[13] = 0;
    SS_ref_db.W[14] = 0;
    
    SS_ref_db.v[0] = 1.0;
    SS_ref_db.v[1] = 1.0;
    SS_ref_db.v[2] = 2.5;
    SS_ref_db.v[3] = 2.5;
    SS_ref_db.v[4] = 1.0;
    SS_ref_db.v[5] = 1.0;
    
    
    em_data py_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
											SS_ref_db.P,
											SS_ref_db.T,
    										"py", 
    										"equilibrium"	);
    
    em_data alm_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
											SS_ref_db.P,
											SS_ref_db.T,
    										"alm", 
    										"equilibrium"	);
    
    em_data gr_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
											SS_ref_db.P,
											SS_ref_db.T,
    										"gr", 
    										"equilibrium"	);
    
    em_data andr_eq 	= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
											SS_ref_db.P,
											SS_ref_db.T,
    										"andr", 
    										"equilibrium"	);
    
    em_data knor_eq 	= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
											SS_ref_db.P,
											SS_ref_db.T,
    										"knor", 
    										"equilibrium"	);
    
    em_data ru_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
											SS_ref_db.P,
											SS_ref_db.T,
    										"ru", 
    										"equilibrium"	);
    
    em_data per_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
											SS_ref_db.P,
											SS_ref_db.T,
    										"per", 
    										"equilibrium"	);
    
    em_data cor_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
											SS_ref_db.P,
											SS_ref_db.T,
    										"cor", 
    										"equilibrium"	);
    
    SS_ref_db.gbase[0] 		= py_eq.gb;
    SS_ref_db.gbase[1] 		= alm_eq.gb;
    SS_ref_db.gbase[2] 		= gr_eq.gb;
    SS_ref_db.gbase[3] 		= andr_eq.gb;
    SS_ref_db.gbase[4] 		= knor_eq.gb + 18.2;
    SS_ref_db.gbase[5] 		= -0.0173*SS_ref_db.T - 0.5*cor_eq.gb + 0.5*per_eq.gb + py_eq.gb + 0.5*ru_eq.gb + 46.7;
    
    SS_ref_db.ElShearMod[0] 	= py_eq.ElShearMod;
    SS_ref_db.ElShearMod[1] 	= alm_eq.ElShearMod;
    SS_ref_db.ElShearMod[2] 	= gr_eq.ElShearMod;
    SS_ref_db.ElShearMod[3] 	= andr_eq.ElShearMod;
    SS_ref_db.ElShearMod[4] 	= knor_eq.ElShearMod;
    SS_ref_db.ElShearMod[5] 	= -0.5*cor_eq.ElShearMod + 0.5*per_eq.ElShearMod + py_eq.ElShearMod + 0.5*ru_eq.ElShearMod;
    
    for (i = 0; i < len_ox; i++){
        SS_ref_db.Comp[0][i] 	= py_eq.C[i];
        SS_ref_db.Comp[1][i] 	= alm_eq.C[i];
        SS_ref_db.Comp[2][i] 	= gr_eq.C[i];
        SS_ref_db.Comp[3][i] 	= andr_eq.C[i];
        SS_ref_db.Comp[4][i] 	= knor_eq.C[i];
        SS_ref_db.Comp[5][i] 	= -0.5*cor_eq.C[i] + 0.5*per_eq.C[i] + py_eq.C[i] + 0.5*ru_eq.C[i];
    }
    
    for (i = 0; i < n_em; i++){
        SS_ref_db.z_em[i] = 1.0;
    };
    
    SS_ref_db.bounds_ref[0][0] = 0.0+eps;  SS_ref_db.bounds_ref[0][1] = 1.0-eps;
    SS_ref_db.bounds_ref[1][0] = 0.0+eps;  SS_ref_db.bounds_ref[1][1] = 1.0-eps;
    SS_ref_db.bounds_ref[2][0] = 0.0+eps;  SS_ref_db.bounds_ref[2][1] = 1.0-eps;
    SS_ref_db.bounds_ref[3][0] = 0.0+eps;  SS_ref_db.bounds_ref[3][1] = 1.0-eps;
    SS_ref_db.bounds_ref[4][0] = 0.0+eps;  SS_ref_db.bounds_ref[4][1] = 1.0-eps;

	if (z_b.bulk_rock[9] == 0.){ 					
		SS_ref_db.z_em[4]          = 0.0;
        SS_ref_db.d_em[4]          = 1.0;
		SS_ref_db.bounds_ref[3][0] = 0.0; 
		SS_ref_db.bounds_ref[3][1] = 0.0;	
	}
	if (z_b.bulk_rock[8] == 0.){ 					
		SS_ref_db.z_em[3]          = 0.0;
        SS_ref_db.d_em[3]          = 1.0;
		SS_ref_db.bounds_ref[2][0] = 0.0; 
		SS_ref_db.bounds_ref[2][1] = 0.0;	
	}
	if (z_b.bulk_rock[7] == 0.){ 					
		SS_ref_db.z_em[5]          = 0.0;
        SS_ref_db.d_em[5]          = 1.0;
		SS_ref_db.bounds_ref[4][0] = 0.0; 
		SS_ref_db.bounds_ref[4][1] = 0.0;	
	}
	return SS_ref_db;
}

/**
  retrieve reference thermodynamic data for horblende 
*/
SS_ref G_SS_ig_hb_function(SS_ref SS_ref_db, int EM_database, int len_ox, bulk_info z_b, double eps){
    
    int i, j;
    int n_em = SS_ref_db.n_em;
    
    char   *EM_tmp[] 		= {"tr","tsm","prgm","glm","cumm","grnm","a","b","mrb","kprg","tts"};
    for (int i = 0; i < SS_ref_db.n_em; i++){
        strcpy(SS_ref_db.EM_list[i],EM_tmp[i]);
    };
    int n_xeos = SS_ref_db.n_xeos;
    char   *CV_tmp[] 		= {"x","y","z","a","k","c","f","t","Q1","Q2"};
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
        strcpy(SS_ref_db.CV_list[i],CV_tmp[i]);
    };
    
    SS_ref_db.W[0] = 20.0;
    SS_ref_db.W[1] = 25.0;
    SS_ref_db.W[2] = 65.0;
    SS_ref_db.W[3] = 45.0;
    SS_ref_db.W[4] = 75.0;
    SS_ref_db.W[5] = 57.0;
    SS_ref_db.W[6] = 63.0;
    SS_ref_db.W[7] = 52.0;
    SS_ref_db.W[8] = 30.0;
    SS_ref_db.W[9] = 85.0;
    SS_ref_db.W[10] = -40.0;
    SS_ref_db.W[11] = 25.0;
    SS_ref_db.W[12] = 70.0;
    SS_ref_db.W[13] = 80.0;
    SS_ref_db.W[14] = 70.0;
    SS_ref_db.W[15] = 72.5;
    SS_ref_db.W[16] = 20.0;
    SS_ref_db.W[17] = -40.0;
    SS_ref_db.W[18] = 35.0;
    SS_ref_db.W[19] = 50.0;
    SS_ref_db.W[20] = 90.0;
    SS_ref_db.W[21] = 106.7;
    SS_ref_db.W[22] = 94.8;
    SS_ref_db.W[23] = 94.8;
    SS_ref_db.W[24] = 40.0;
    SS_ref_db.W[25] = 8.0;
    SS_ref_db.W[26] = 15.0;
    SS_ref_db.W[27] = 100.0;
    SS_ref_db.W[28] = 113.5;
    SS_ref_db.W[29] = 100.0;
    SS_ref_db.W[30] = 111.2;
    SS_ref_db.W[31] = 0.0;
    SS_ref_db.W[32] = 54.0;
    SS_ref_db.W[33] = 75.0;
    SS_ref_db.W[34] = 33.0;
    SS_ref_db.W[35] = 18.0;
    SS_ref_db.W[36] = 23.0;
    SS_ref_db.W[37] = 80.0;
    SS_ref_db.W[38] = 87.0;
    SS_ref_db.W[39] = 100.0;
    SS_ref_db.W[40] = 12.0;
    SS_ref_db.W[41] = 8.0;
    SS_ref_db.W[42] = 91.0;
    SS_ref_db.W[43] = 96.0;
    SS_ref_db.W[44] = 65.0;
    SS_ref_db.W[45] = 20.0;
    SS_ref_db.W[46] = 80.0;
    SS_ref_db.W[47] = 94.0;
    SS_ref_db.W[48] = 95.0;
    SS_ref_db.W[49] = 90.0;
    SS_ref_db.W[50] = 94.0;
    SS_ref_db.W[51] = 95.0;
    SS_ref_db.W[52] = 50.0;
    SS_ref_db.W[53] = 50.0;
    SS_ref_db.W[54] = 35.0;
    
    SS_ref_db.v[0] = 1.0;
    SS_ref_db.v[1] = 1.5;
    SS_ref_db.v[2] = 1.7;
    SS_ref_db.v[3] = 0.8;
    SS_ref_db.v[4] = 1.0;
    SS_ref_db.v[5] = 1.0;
    SS_ref_db.v[6] = 1.0;
    SS_ref_db.v[7] = 1.0;
    SS_ref_db.v[8] = 0.8;
    SS_ref_db.v[9] = 1.7;
    SS_ref_db.v[10] = 1.5;
    
    
    em_data tr_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
											SS_ref_db.P,
											SS_ref_db.T,
    										"tr", 
    										"equilibrium"	);
    
    em_data ts_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
											SS_ref_db.P,
											SS_ref_db.T,
    										"ts", 
    										"equilibrium"	);
    
    em_data parg_eq 	= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
											SS_ref_db.P,
											SS_ref_db.T,
    										"parg", 
    										"equilibrium"	);
    
    em_data gl_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
											SS_ref_db.P,
											SS_ref_db.T,
    										"gl", 
    										"equilibrium"	);
    
    em_data cumm_eq 	= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
											SS_ref_db.P,
											SS_ref_db.T,
    										"cumm", 
    										"equilibrium"	);
    
    em_data grun_eq 	= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
											SS_ref_db.P,
											SS_ref_db.T,
    										"grun", 
    										"equilibrium"	);
    
    em_data gr_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
											SS_ref_db.P,
											SS_ref_db.T,
    										"gr", 
    										"equilibrium"	);
    
    em_data andr_eq 	= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
											SS_ref_db.P,
											SS_ref_db.T,
    										"andr", 
    										"equilibrium"	);
    
    em_data pa_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
											SS_ref_db.P,
											SS_ref_db.T,
    										"pa", 
    										"equilibrium"	);
    
    em_data mu_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
											SS_ref_db.P,
											SS_ref_db.T,
    										"mu", 
    										"equilibrium"	);
    
    em_data ru_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
											SS_ref_db.P,
											SS_ref_db.T,
    										"ru", 
    										"equilibrium"	);
    
    em_data dsp_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
											SS_ref_db.P,
											SS_ref_db.T,
    										"dsp", 
    										"equilibrium"	);
    
    SS_ref_db.gbase[0] 		= tr_eq.gb;
    SS_ref_db.gbase[1] 		= ts_eq.gb + 10.0;
    SS_ref_db.gbase[2] 		= parg_eq.gb - 10.0;
    SS_ref_db.gbase[3] 		= gl_eq.gb - 3.0;
    SS_ref_db.gbase[4] 		= cumm_eq.gb;
    SS_ref_db.gbase[5] 		= grun_eq.gb - 3.0;
    SS_ref_db.gbase[6] 		= 3.0/7.0*cumm_eq.gb + 4.0/7.0*grun_eq.gb - 11.2;
    SS_ref_db.gbase[7] 		= 2.0/7.0*cumm_eq.gb + 5.0/7.0*grun_eq.gb - 13.8;
    SS_ref_db.gbase[8] 		= andr_eq.gb + gl_eq.gb - gr_eq.gb;
    SS_ref_db.gbase[9] 		= 0.02*SS_ref_db.T + mu_eq.gb - pa_eq.gb + parg_eq.gb - 7.06;
    SS_ref_db.gbase[10] 		= -2.0*dsp_eq.gb + 2.0*ru_eq.gb + ts_eq.gb + 95.0;
    
    SS_ref_db.ElShearMod[0] 	= tr_eq.ElShearMod;
    SS_ref_db.ElShearMod[1] 	= ts_eq.ElShearMod;
    SS_ref_db.ElShearMod[2] 	= parg_eq.ElShearMod;
    SS_ref_db.ElShearMod[3] 	= gl_eq.ElShearMod;
    SS_ref_db.ElShearMod[4] 	= cumm_eq.ElShearMod;
    SS_ref_db.ElShearMod[5] 	= grun_eq.ElShearMod;
    SS_ref_db.ElShearMod[6] 	= 3.0/7.0*cumm_eq.ElShearMod + 4.0/7.0*grun_eq.ElShearMod;
    SS_ref_db.ElShearMod[7] 	= 2.0/7.0*cumm_eq.ElShearMod + 5.0/7.0*grun_eq.ElShearMod;
    SS_ref_db.ElShearMod[8] 	= andr_eq.ElShearMod + gl_eq.ElShearMod - gr_eq.ElShearMod;
    SS_ref_db.ElShearMod[9] 	= mu_eq.ElShearMod - pa_eq.ElShearMod + parg_eq.ElShearMod;
    SS_ref_db.ElShearMod[10] 	= -2.0*dsp_eq.ElShearMod + 2.0*ru_eq.ElShearMod + ts_eq.ElShearMod;
    
    for (i = 0; i < len_ox; i++){
        SS_ref_db.Comp[0][i] 	= tr_eq.C[i];
        SS_ref_db.Comp[1][i] 	= ts_eq.C[i];
        SS_ref_db.Comp[2][i] 	= parg_eq.C[i];
        SS_ref_db.Comp[3][i] 	= gl_eq.C[i];
        SS_ref_db.Comp[4][i] 	= cumm_eq.C[i];
        SS_ref_db.Comp[5][i] 	= grun_eq.C[i];
        SS_ref_db.Comp[6][i] 	= 3.0/7.0*cumm_eq.C[i] + 4.0/7.0*grun_eq.C[i];
        SS_ref_db.Comp[7][i] 	= 2.0/7.0*cumm_eq.C[i] + 5.0/7.0*grun_eq.C[i];
        SS_ref_db.Comp[8][i] 	= andr_eq.C[i] + gl_eq.C[i] - gr_eq.C[i];
        SS_ref_db.Comp[9][i] 	= mu_eq.C[i] - pa_eq.C[i] + parg_eq.C[i];
        SS_ref_db.Comp[10][i] 	= -2.0*dsp_eq.C[i] + 2.0*ru_eq.C[i] + ts_eq.C[i];
    }
    
    for (i = 0; i < n_em; i++){
        SS_ref_db.z_em[i] = 1.0;
    };
    
    SS_ref_db.bounds_ref[0][0] = 0.0+eps;  SS_ref_db.bounds_ref[0][1] = 1.0-eps;
    SS_ref_db.bounds_ref[1][0] = 0.0+eps;  SS_ref_db.bounds_ref[1][1] = 1.0-eps;
    SS_ref_db.bounds_ref[2][0] = 0.0+eps;  SS_ref_db.bounds_ref[2][1] = 1.0-eps;
    SS_ref_db.bounds_ref[3][0] = 0.0+eps;  SS_ref_db.bounds_ref[3][1] = 1.0-eps;
    SS_ref_db.bounds_ref[4][0] = 0.0+eps;  SS_ref_db.bounds_ref[4][1] = 1.0-eps;
    SS_ref_db.bounds_ref[5][0] = 0.0+eps;  SS_ref_db.bounds_ref[5][1] = 1.0-eps;
    SS_ref_db.bounds_ref[6][0] = 0.0+eps;  SS_ref_db.bounds_ref[6][1] = 1.0-eps;
    SS_ref_db.bounds_ref[7][0] = 0.0+eps;  SS_ref_db.bounds_ref[7][1] = 1.0-eps;
    SS_ref_db.bounds_ref[8][0] = -1.0+eps;  SS_ref_db.bounds_ref[8][1] = 1.0-eps;
    SS_ref_db.bounds_ref[9][0] = -1.0+eps;  SS_ref_db.bounds_ref[9][1] = 1.0-eps;
    

	if (z_b.bulk_rock[8] == 0.){ 					
		SS_ref_db.z_em[8]          = 0.0;
        SS_ref_db.d_em[8]          = 1.0;
		SS_ref_db.bounds_ref[6][0] = 0.0; 
		SS_ref_db.bounds_ref[6][1] = 0.0;	
	}
	if (z_b.bulk_rock[7] == 0.){ 					
		SS_ref_db.z_em[10]          = 0.0;
        SS_ref_db.d_em[10]          = 1.0;
		SS_ref_db.bounds_ref[7][0] = 0.0; 
		SS_ref_db.bounds_ref[7][1] = 0.0;	
	}
    

    /* this lists the index of the order variables */
    // SS_ref_db.orderVar      = 1;
    // SS_ref_db.idOrderVar[8] = -1.0;
    // SS_ref_db.idOrderVar[9] = -1.0;


	return SS_ref_db;	
}


/**
   retrieve reference thermodynamic data for ig_ilm
*/
SS_ref G_SS_ig_ilm_function(SS_ref SS_ref_db, int EM_database, int len_ox, bulk_info z_b, double eps){
    
    int i, j;
    int n_em = SS_ref_db.n_em;
    
    char   *EM_tmp[] 		= {"oilm","dilm","dhem"};
    for (int i = 0; i < SS_ref_db.n_em; i++){
        strcpy(SS_ref_db.EM_list[i],EM_tmp[i]);
    };

    int n_xeos = SS_ref_db.n_xeos;
    char   *CV_tmp[] 		= {"x","Q"};
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
        strcpy(SS_ref_db.CV_list[i],CV_tmp[i]);
    };

 
    SS_ref_db.W[0] = 7.05;
    SS_ref_db.W[1] = 14.3;
    SS_ref_db.W[2] = 7.25;
    
    
    em_data ilm_or 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
											SS_ref_db.P,
											SS_ref_db.T,
    										"ilm", 
    										"ordered"	);
    
    em_data ilm_di 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
											SS_ref_db.P,
											SS_ref_db.T,
    										"ilm", 
    										"disordered"	);
    
    em_data hem_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
											SS_ref_db.P,
											SS_ref_db.T,
    										"hem", 
    										"equilibrium"	);
    
    SS_ref_db.gbase[0] 		= ilm_or.gb;
    SS_ref_db.gbase[1] 		= ilm_di.gb;
    SS_ref_db.gbase[2] 		= hem_eq.gb;
    
    SS_ref_db.ElShearMod[0] 	= ilm_or.ElShearMod;
    SS_ref_db.ElShearMod[1] 	= ilm_di.ElShearMod;
    SS_ref_db.ElShearMod[2] 	= hem_eq.ElShearMod;
    
    for (i = 0; i < len_ox; i++){
        SS_ref_db.Comp[0][i] 	= ilm_or.C[i];
        SS_ref_db.Comp[1][i] 	= ilm_di.C[i];
        SS_ref_db.Comp[2][i] 	= hem_eq.C[i];
    }
    
    for (i = 0; i < n_em; i++){
        SS_ref_db.z_em[i] = 1.0;
    };
    
    SS_ref_db.bounds_ref[0][0] = 0.0+eps;  SS_ref_db.bounds_ref[0][1] = 1.0-eps;
    SS_ref_db.bounds_ref[1][0] = -0.99+eps;  SS_ref_db.bounds_ref[1][1] = 0.99-eps;

	if (z_b.bulk_rock[8] == 0.){ 					
		SS_ref_db.z_em[2]          = 0.0;
        SS_ref_db.d_em[2]          = 1.0;
		SS_ref_db.bounds_ref[0][0] = 1.0; 
		SS_ref_db.bounds_ref[0][1] = 1.0;	
	}

    return SS_ref_db;
}
/**
   retrieve reference thermodynamic data for ig_liqHw
*/
SS_ref G_SS_ig_liq_function(SS_ref SS_ref_db, int EM_database, int len_ox, bulk_info z_b, double eps){
        
    int i, j;
    int n_em = SS_ref_db.n_em;
    
    char   *EM_tmp[] 		= {"q4L","slL","wo1L","fo2L","fa2L","jdL","hmL","ekL","tiL","kjL","ctL","wat1L"};
    for (int i = 0; i < SS_ref_db.n_em; i++){
        strcpy(SS_ref_db.EM_list[i],EM_tmp[i]);
    };
    int n_xeos = SS_ref_db.n_xeos;
    char   *CV_tmp[] 		= {"wo","sl","fo","fa","jd","hm","ek","ti","kj","yct","h2o"};
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
        strcpy(SS_ref_db.CV_list[i],CV_tmp[i]);
    };
    SS_ref_db.W[0] = 9.5 - 0.1*SS_ref_db.P;
    SS_ref_db.W[1] = -10.3;
    SS_ref_db.W[2] = -3.12*SS_ref_db.P - 26.5;
    SS_ref_db.W[3] = -0.55*SS_ref_db.P - 12.0;
    SS_ref_db.W[4] = -0.13*SS_ref_db.P - 15.1;
    SS_ref_db.W[5] = 20.0;
    SS_ref_db.W[6] = 0;
    SS_ref_db.W[7] = 24.6;
    SS_ref_db.W[8] = -0.05*SS_ref_db.P - 17.8;
    SS_ref_db.W[9] = -14.6;
    SS_ref_db.W[10] = 17.8 - 0.61*SS_ref_db.P;
    SS_ref_db.W[11] = 0.85*SS_ref_db.P - 26.5;
    SS_ref_db.W[12] = 2.2;
    SS_ref_db.W[13] = 2.5;
    SS_ref_db.W[14] = 16.8;
    SS_ref_db.W[15] = -5.0;
    SS_ref_db.W[16] = 0;
    SS_ref_db.W[17] = 15.2 - 0.04*SS_ref_db.P;
    SS_ref_db.W[18] = 7.0;
    SS_ref_db.W[19] = 4.0;
    SS_ref_db.W[20] = 23.7 - 0.94*SS_ref_db.P;
    SS_ref_db.W[21] = 0.11*SS_ref_db.P + 25.5;
    SS_ref_db.W[22] = 14.0;
    SS_ref_db.W[23] = -1.2;
    SS_ref_db.W[24] = 0;
    SS_ref_db.W[25] = 0;
    SS_ref_db.W[26] = 18.0;
    SS_ref_db.W[27] = -1.1;
    SS_ref_db.W[28] = 9.5;
    SS_ref_db.W[29] = 40.3 - 0.86*SS_ref_db.P;
    SS_ref_db.W[30] = 18.0;
    SS_ref_db.W[31] = 1.5;
    SS_ref_db.W[32] = 0;
    SS_ref_db.W[33] = 0;
    SS_ref_db.W[34] = 7.5;
    SS_ref_db.W[35] = 3.0;
    SS_ref_db.W[36] = -5.6;
    SS_ref_db.W[37] = 9.4 - 1.58*SS_ref_db.P;
    SS_ref_db.W[38] = 7.5 - 0.05*SS_ref_db.P;
    SS_ref_db.W[39] = -30.0;
    SS_ref_db.W[40] = 0;
    SS_ref_db.W[41] = 6.7;
    SS_ref_db.W[42] = 10.0;
    SS_ref_db.W[43] = -6.5;
    SS_ref_db.W[44] = 9.2 - 1.58*SS_ref_db.P;
    SS_ref_db.W[45] = 10.0;
    SS_ref_db.W[46] = 0;
    SS_ref_db.W[47] = 0.14*SS_ref_db.P + 16.5;
    SS_ref_db.W[48] = -5.9;
    SS_ref_db.W[49] = 7.6;
    SS_ref_db.W[50] = -0.06*SS_ref_db.P - 8.3;
    SS_ref_db.W[51] = 0;
    SS_ref_db.W[52] = 0;
    SS_ref_db.W[53] = 10.0;
    SS_ref_db.W[54] = 0;
    SS_ref_db.W[55] = 60.0 - 0.66*SS_ref_db.P;
    SS_ref_db.W[56] = 0;
    SS_ref_db.W[57] = 0;
    SS_ref_db.W[58] = 0;
    SS_ref_db.W[59] = 30.0 - 0.66*SS_ref_db.P;
    SS_ref_db.W[60] = 9.0;
    SS_ref_db.W[61] = 0;
    SS_ref_db.W[62] = 30.0 - 0.6*SS_ref_db.P;
    SS_ref_db.W[63] = -5.6;
    SS_ref_db.W[64] = 0.22*SS_ref_db.P - 0.1;
    SS_ref_db.W[65] = 0.05*SS_ref_db.P + 17.3;
    
    SS_ref_db.v[0] = 100.0;
    SS_ref_db.v[1] = 120.0;
    SS_ref_db.v[2] = 140.0;
    SS_ref_db.v[3] = 240.0;
    SS_ref_db.v[4] = 100.0;
    SS_ref_db.v[5] = 120.0;
    SS_ref_db.v[6] = 100.0;
    SS_ref_db.v[7] = 100.0;
    SS_ref_db.v[8] = 100.0;
    SS_ref_db.v[9] = 100.0;
    SS_ref_db.v[10] = 100.0;
    SS_ref_db.v[11] = 100.0;
    
    
    em_data qL_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"qL", 
    										"equilibrium"	);
    
    em_data silL_eq 	= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"silL", 
    										"equilibrium"	);
    
    em_data woL_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"woL", 
    										"equilibrium"	);
    
    em_data foL_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"foL", 
    										"equilibrium"	);
    
    em_data faL_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"faL", 
    										"equilibrium"	);
    
    em_data abL_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"abL", 
    										"equilibrium"	);
    
    em_data hemL_eq 	= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"hemL", 
    										"equilibrium"	);
    
    em_data eskL_eq 	= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"eskL", 
    										"equilibrium"	);
    
    em_data ruL_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"ruL", 
    										"equilibrium"	);
    
    em_data kspL_eq     = get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"kspL", 
    										"equilibrium"	);
    
    em_data h2oL_eq 	= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"h2oL", 
    										"equilibrium"	);
    
    SS_ref_db.gbase[0] 		= -0.059*SS_ref_db.P + 4.0*qL_eq.gb + 0.22;
    SS_ref_db.gbase[1] 		= -0.318*SS_ref_db.P + silL_eq.gb + 6.2;
    SS_ref_db.gbase[2] 		= -0.114*SS_ref_db.P + woL_eq.gb - 0.45;
    SS_ref_db.gbase[3] 		= -0.131*SS_ref_db.P + 2.0*foL_eq.gb + 8.67;
    SS_ref_db.gbase[4] 		= -0.055*SS_ref_db.P + 2.0*faL_eq.gb + 13.7;
    SS_ref_db.gbase[5] 		= -0.089*SS_ref_db.P + abL_eq.gb -qL_eq.gb + 12.19;
    SS_ref_db.gbase[6] 		= -0.032*SS_ref_db.P + 0.5*hemL_eq.gb + 3.3;
    SS_ref_db.gbase[7] 		= 0.245*SS_ref_db.P + 0.5*eskL_eq.gb + 24.85;
    SS_ref_db.gbase[8] 		= -0.489*SS_ref_db.P + ruL_eq.gb + 5.58;
    SS_ref_db.gbase[9] 		= -0.21*SS_ref_db.P + kspL_eq.gb -qL_eq.gb + 11.98;
    SS_ref_db.gbase[10] 		= 0.053*SS_ref_db.P + 0.055*SS_ref_db.T -qL_eq.gb + silL_eq.gb + woL_eq.gb - 108.3;
    SS_ref_db.gbase[11] 		= 0.00087*SS_ref_db.P - 0.0039*SS_ref_db.T + h2oL_eq.gb + 3.2;
    
    SS_ref_db.ElShearMod[0] 	= 4.0*qL_eq.ElShearMod;
    SS_ref_db.ElShearMod[1] 	= silL_eq.ElShearMod;
    SS_ref_db.ElShearMod[2] 	= woL_eq.ElShearMod;
    SS_ref_db.ElShearMod[3] 	= 2.0*foL_eq.ElShearMod;
    SS_ref_db.ElShearMod[4] 	= 2.0*faL_eq.ElShearMod;
    SS_ref_db.ElShearMod[5] 	= abL_eq.ElShearMod -qL_eq.ElShearMod;
    SS_ref_db.ElShearMod[6] 	= 0.5*hemL_eq.ElShearMod;
    SS_ref_db.ElShearMod[7] 	= 0.5*eskL_eq.ElShearMod;
    SS_ref_db.ElShearMod[8] 	= ruL_eq.ElShearMod;
    SS_ref_db.ElShearMod[9] 	= kspL_eq.ElShearMod -qL_eq.ElShearMod;
    SS_ref_db.ElShearMod[10] 	= -qL_eq.ElShearMod + silL_eq.ElShearMod + woL_eq.ElShearMod;
    SS_ref_db.ElShearMod[11] 	= h2oL_eq.ElShearMod;
    
    for (i = 0; i < len_ox; i++){
        SS_ref_db.Comp[0][i] 	= 4.0*qL_eq.C[i];
        SS_ref_db.Comp[1][i] 	= silL_eq.C[i];
        SS_ref_db.Comp[2][i] 	= woL_eq.C[i];
        SS_ref_db.Comp[3][i] 	= 2.0*foL_eq.C[i];
        SS_ref_db.Comp[4][i] 	= 2.0*faL_eq.C[i];
        SS_ref_db.Comp[5][i] 	= abL_eq.C[i] -qL_eq.C[i];
        SS_ref_db.Comp[6][i] 	= 0.5*hemL_eq.C[i];
        SS_ref_db.Comp[7][i] 	= 0.5*eskL_eq.C[i];
        SS_ref_db.Comp[8][i] 	= ruL_eq.C[i];
        SS_ref_db.Comp[9][i] 	= kspL_eq.C[i] -qL_eq.C[i];
        SS_ref_db.Comp[10][i] 	= -qL_eq.C[i] + silL_eq.C[i] + woL_eq.C[i];
        SS_ref_db.Comp[11][i] 	= h2oL_eq.C[i];
    }
    
    for (i = 0; i < n_em; i++){
        SS_ref_db.z_em[i] = 1.0;
    };
    
    SS_ref_db.bounds_ref[0][0] = 0.0+eps;  SS_ref_db.bounds_ref[0][1] = 1.0-eps;
    SS_ref_db.bounds_ref[1][0] = 0.0+eps;  SS_ref_db.bounds_ref[1][1] = 1.0-eps;
    SS_ref_db.bounds_ref[2][0] = 0.0+eps;  SS_ref_db.bounds_ref[2][1] = 1.0-eps;
    SS_ref_db.bounds_ref[3][0] = 0.0+eps;  SS_ref_db.bounds_ref[3][1] = 1.0-eps;
    SS_ref_db.bounds_ref[4][0] = 0.0+eps;  SS_ref_db.bounds_ref[4][1] = 1.0-eps;
    SS_ref_db.bounds_ref[5][0] = 0.0+eps;  SS_ref_db.bounds_ref[5][1] = 1.0-eps;
    SS_ref_db.bounds_ref[6][0] = 0.0+eps;  SS_ref_db.bounds_ref[6][1] = 1.0-eps;
    SS_ref_db.bounds_ref[7][0] = 0.0+eps;  SS_ref_db.bounds_ref[7][1] = 1.0-eps;
    SS_ref_db.bounds_ref[8][0] = 0.0+eps;  SS_ref_db.bounds_ref[8][1] = 1.0-eps;
    SS_ref_db.bounds_ref[9][0] = 0.0+eps;  SS_ref_db.bounds_ref[9][1] = 1.0-eps;
    SS_ref_db.bounds_ref[10][0] = 0.0+eps;  SS_ref_db.bounds_ref[10][1] = 1.0-eps;
    
	if (z_b.bulk_rock[10] == 0.){ 					
		SS_ref_db.z_em[11]          = 0.0;
        SS_ref_db.d_em[11]           = 1.0;
		SS_ref_db.bounds_ref[10][0] = 0.0; 
		SS_ref_db.bounds_ref[10][1] = 0.0;	
	}
	if (z_b.bulk_rock[9] == 0.){ 					
		SS_ref_db.z_em[7]          = 0.0;
        SS_ref_db.d_em[7]          = 1.0;
		SS_ref_db.bounds_ref[6][0] = 0.0; 
		SS_ref_db.bounds_ref[6][1] = 0.0;	
	}
	if (z_b.bulk_rock[8] == 0.){ 					
		SS_ref_db.z_em[6]          = 0.0;
        SS_ref_db.d_em[6]          = 1.0;
		SS_ref_db.bounds_ref[5][0] = 0.0; 
		SS_ref_db.bounds_ref[5][1] = 0.0;	
	}
	if (z_b.bulk_rock[7] == 0.){ 					
		SS_ref_db.z_em[8]          = 0.0;
        SS_ref_db.d_em[8]          = 1.0;
		SS_ref_db.bounds_ref[7][0] = 0.0; 
		SS_ref_db.bounds_ref[7][1] = 0.0;	
	}
    return SS_ref_db;
}



/**
  retrieve reference thermodynamic data for muscovite
*/
SS_ref G_SS_ig_mu_function(SS_ref SS_ref_db, int EM_database, int len_ox, bulk_info z_b, double eps){
    
    int i, j;
    int n_em = SS_ref_db.n_em;
    
    char   *EM_tmp[] 		= {"mu","cel","fcel","pa","mam","fmu"};
    for (int i = 0; i < SS_ref_db.n_em; i++){
        strcpy(SS_ref_db.EM_list[i],EM_tmp[i]);
    };
    int n_xeos = SS_ref_db.n_xeos;
    char   *CV_tmp[] 		= {"x","y","f","n","c"};
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
        strcpy(SS_ref_db.CV_list[i],CV_tmp[i]);
    };
    SS_ref_db.W[0] = 0.2*SS_ref_db.P;
    SS_ref_db.W[1] = 0.2*SS_ref_db.P;
    SS_ref_db.W[2] = 0.353*SS_ref_db.P + 0.0034*SS_ref_db.T + 10.12;
    SS_ref_db.W[3] = 35.0;
    SS_ref_db.W[4] = 0.0;
    SS_ref_db.W[5] = 0.0;
    SS_ref_db.W[6] = 0.25*SS_ref_db.P + 45.0;
    SS_ref_db.W[7] = 50.0;
    SS_ref_db.W[8] = 0.0;
    SS_ref_db.W[9] = 0.25*SS_ref_db.P + 45.0;
    SS_ref_db.W[10] = 50.0;
    SS_ref_db.W[11] = 0.0;
    SS_ref_db.W[12] = 15.0;
    SS_ref_db.W[13] = 30.0;
    SS_ref_db.W[14] = 35.0;
    
    SS_ref_db.v[0] = 0.63;
    SS_ref_db.v[1] = 0.63;
    SS_ref_db.v[2] = 0.63;
    SS_ref_db.v[3] = 0.37;
    SS_ref_db.v[4] = 0.63;
    SS_ref_db.v[5] = 0.63;
    
    
    em_data mu_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
											SS_ref_db.P,
											SS_ref_db.T,
    										"mu", 
    										"equilibrium"	);
    
    em_data cel_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
											SS_ref_db.P,
											SS_ref_db.T,
    										"cel", 
    										"equilibrium"	);
    
    em_data fcel_eq 	= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
											SS_ref_db.P,
											SS_ref_db.T,
    										"fcel", 
    										"equilibrium"	);
    
    em_data pa_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
											SS_ref_db.P,
											SS_ref_db.T,
    										"pa", 
    										"equilibrium"	);
    
    em_data ma_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
											SS_ref_db.P,
											SS_ref_db.T,
    										"ma", 
    										"equilibrium"	);
    
    em_data gr_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
											SS_ref_db.P,
											SS_ref_db.T,
    										"gr", 
    										"equilibrium"	);
    
    em_data andr_eq 	= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
											SS_ref_db.P,
											SS_ref_db.T,
    										"andr", 
    										"equilibrium"	);
    
    SS_ref_db.gbase[0] 		= mu_eq.gb;
    SS_ref_db.gbase[1] 		= cel_eq.gb;
    SS_ref_db.gbase[2] 		= fcel_eq.gb;
    SS_ref_db.gbase[3] 		= pa_eq.gb;
    SS_ref_db.gbase[4] 		= ma_eq.gb + 5.0; //6.5 from Eleanor?
    SS_ref_db.gbase[5] 		= 0.5*andr_eq.gb - 0.5*gr_eq.gb + mu_eq.gb + 25.0;
    
    SS_ref_db.ElShearMod[0] 	= mu_eq.ElShearMod;
    SS_ref_db.ElShearMod[1] 	= cel_eq.ElShearMod;
    SS_ref_db.ElShearMod[2] 	= fcel_eq.ElShearMod;
    SS_ref_db.ElShearMod[3] 	= pa_eq.ElShearMod;
    SS_ref_db.ElShearMod[4] 	= ma_eq.ElShearMod;
    SS_ref_db.ElShearMod[5] 	= 0.5*andr_eq.ElShearMod - 0.5*gr_eq.ElShearMod + mu_eq.ElShearMod;
    
    for (i = 0; i < len_ox; i++){
        SS_ref_db.Comp[0][i] 	= mu_eq.C[i];
        SS_ref_db.Comp[1][i] 	= cel_eq.C[i];
        SS_ref_db.Comp[2][i] 	= fcel_eq.C[i];
        SS_ref_db.Comp[3][i] 	= pa_eq.C[i];
        SS_ref_db.Comp[4][i] 	= ma_eq.C[i];
        SS_ref_db.Comp[5][i] 	= 0.5*andr_eq.C[i] - 0.5*gr_eq.C[i] + mu_eq.C[i];
    }
    
    for (i = 0; i < n_em; i++){
        SS_ref_db.z_em[i] = 1.0;
    };
    
    SS_ref_db.bounds_ref[0][0] = 0.0+eps;  SS_ref_db.bounds_ref[0][1] = 1.0-eps;
    SS_ref_db.bounds_ref[1][0] = 0.0+eps;  SS_ref_db.bounds_ref[1][1] = 1.0-eps;
    SS_ref_db.bounds_ref[2][0] = 0.0+eps;  SS_ref_db.bounds_ref[2][1] = 1.0-eps;
    SS_ref_db.bounds_ref[3][0] = 0.0+eps;  SS_ref_db.bounds_ref[3][1] = 1.0-eps;
    SS_ref_db.bounds_ref[4][0] = 0.0+eps;  SS_ref_db.bounds_ref[4][1] = 1.0-eps;
    
	/* define box bounds according to bulk-rock */
	if (z_b.bulk_rock[8] == 0.){
		SS_ref_db.z_em[5]          = 0.0;
		SS_ref_db.bounds_ref[2][0] = eps;
		SS_ref_db.bounds_ref[2][1] = eps;
	}

	return SS_ref_db;
}

/**
  retrieve reference thermodynamic data for olivine
*/
SS_ref G_SS_ig_ol_function(SS_ref SS_ref_db, int EM_database, int len_ox, bulk_info z_b, double eps){
    
    int i, j;
    int n_em = SS_ref_db.n_em;
    
    char   *EM_tmp[] 		= {"mont","fa","fo","cfm"};
    for (int i = 0; i < SS_ref_db.n_em; i++){
        strcpy(SS_ref_db.EM_list[i],EM_tmp[i]);
    };
    int n_xeos = SS_ref_db.n_xeos;
    char   *CV_tmp[] 		= {"x","c","Q"};
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
        strcpy(SS_ref_db.CV_list[i],CV_tmp[i]);
    };
    
    SS_ref_db.W[0] = 24.0;
    SS_ref_db.W[1] = 38.0;
    SS_ref_db.W[2] = 24.0;
    SS_ref_db.W[3] = 9.0;
    SS_ref_db.W[4] = 4.5;
    SS_ref_db.W[5] = 4.5;
    
    
    em_data mont_eq 	= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
											SS_ref_db.P,
											SS_ref_db.T,
    										"mont", 
    										"equilibrium"	);
    
    em_data fa_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
											SS_ref_db.P,
											SS_ref_db.T,
    										"fa", 
    										"equilibrium"	);
    
    em_data fo_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
											SS_ref_db.P,
											SS_ref_db.T,
    										"fo", 
    										"equilibrium"	);
    
    SS_ref_db.gbase[0] 		= mont_eq.gb;
    SS_ref_db.gbase[1] 		= fa_eq.gb;
    SS_ref_db.gbase[2] 		= fo_eq.gb;
    SS_ref_db.gbase[3] 		= 0.5*fa_eq.gb + 0.5*fo_eq.gb;
    
    SS_ref_db.ElShearMod[0] 	= mont_eq.ElShearMod;
    SS_ref_db.ElShearMod[1] 	= fa_eq.ElShearMod;
    SS_ref_db.ElShearMod[2] 	= fo_eq.ElShearMod;
    SS_ref_db.ElShearMod[3] 	= 0.5*fa_eq.ElShearMod + 0.5*fo_eq.ElShearMod;
    
    for (i = 0; i < len_ox; i++){
        SS_ref_db.Comp[0][i] 	= mont_eq.C[i];
        SS_ref_db.Comp[1][i] 	= fa_eq.C[i];
        SS_ref_db.Comp[2][i] 	= fo_eq.C[i];
        SS_ref_db.Comp[3][i] 	= 0.5*fa_eq.C[i] + 0.5*fo_eq.C[i];
    }
    
    for (i = 0; i < n_em; i++){
        SS_ref_db.z_em[i] = 1.0;
    };
    
    SS_ref_db.bounds_ref[0][0] = 0.0+eps;  SS_ref_db.bounds_ref[0][1] = 1.0-eps;
    SS_ref_db.bounds_ref[1][0] = 0.0+eps;  SS_ref_db.bounds_ref[1][1] = 1.0-eps;
    SS_ref_db.bounds_ref[2][0] = -0.5+eps;  SS_ref_db.bounds_ref[2][1] = 0.5-eps;
    
    return SS_ref_db;
}

/**
  retrieve reference thermodynamic data for orthopyroxene
*/
SS_ref G_SS_ig_opx_function(SS_ref SS_ref_db, int EM_database, int len_ox, bulk_info z_b, double eps){
    
    int i, j;
    int n_em = SS_ref_db.n_em;
    
    char   *EM_tmp[] 		= {"en","fs","fm","odi","mgts","cren","obuf","mess","ojd"};
    for (int i = 0; i < SS_ref_db.n_em; i++){
        strcpy(SS_ref_db.EM_list[i],EM_tmp[i]);
    };
    int n_xeos = SS_ref_db.n_xeos;
    char   *CV_tmp[] 		= {"x","y","c","Q","f","t","cr","j"};
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
        strcpy(SS_ref_db.CV_list[i],CV_tmp[i]);
    };
    
    SS_ref_db.W[0] = 7.0;
    SS_ref_db.W[1] = 4.0;
    SS_ref_db.W[2] = 29.4;
    SS_ref_db.W[3] = 12.5 - 0.04*SS_ref_db.P;
    SS_ref_db.W[4] = 8.0;
    SS_ref_db.W[5] = 6.0;
    SS_ref_db.W[6] = 8.0;
    SS_ref_db.W[7] = 35.0;
    SS_ref_db.W[8] = 4.0;
    SS_ref_db.W[9] = 0.08*SS_ref_db.P + 21.5;
    SS_ref_db.W[10] = 11.0 - 0.15*SS_ref_db.P;
    SS_ref_db.W[11] = 10.0;
    SS_ref_db.W[12] = 7.0;
    SS_ref_db.W[13] = 10.0;
    SS_ref_db.W[14] = 35.0;
    SS_ref_db.W[15] = 0.08*SS_ref_db.P + 18.0;
    SS_ref_db.W[16] = 15.0 - 0.15*SS_ref_db.P;
    SS_ref_db.W[17] = 12.0;
    SS_ref_db.W[18] = 8.0;
    SS_ref_db.W[19] = 12.0;
    SS_ref_db.W[20] = 35.0;
    SS_ref_db.W[21] = 75.5 - 0.84*SS_ref_db.P;
    SS_ref_db.W[22] = 20.0;
    SS_ref_db.W[23] = 40.0;
    SS_ref_db.W[24] = 20.0;
    SS_ref_db.W[25] = 35.0;
    SS_ref_db.W[26] = 2.0;
    SS_ref_db.W[27] = 10.0;
    SS_ref_db.W[28] = 2.0;
    SS_ref_db.W[29] = 7.0;
    SS_ref_db.W[30] = 6.0;
    SS_ref_db.W[31] = 2.0;
    SS_ref_db.W[32] = -11.0;
    SS_ref_db.W[33] = 6.0;
    SS_ref_db.W[34] = 20.0;
    SS_ref_db.W[35] = -11.0;
    
    SS_ref_db.v[0] = 1.0;
    SS_ref_db.v[1] = 1.0;
    SS_ref_db.v[2] = 1.0;
    SS_ref_db.v[3] = 1.2;
    SS_ref_db.v[4] = 1.0;
    SS_ref_db.v[5] = 1.0;
    SS_ref_db.v[6] = 1.0;
    SS_ref_db.v[7] = 1.0;
    SS_ref_db.v[8] = 1.2;
    
    
    em_data en_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
											SS_ref_db.P,
											SS_ref_db.T,
    										"en", 
    										"equilibrium"	);
    
    em_data fs_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
											SS_ref_db.P,
											SS_ref_db.T,
    										"fs", 
    										"equilibrium"	);
    
    em_data di_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
											SS_ref_db.P,
											SS_ref_db.T,
    										"di", 
    										"equilibrium"	);
    
    em_data mgts_eq 	= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
											SS_ref_db.P,
											SS_ref_db.T,
    										"mgts", 
    										"equilibrium"	);
    
    em_data kos_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
											SS_ref_db.P,
											SS_ref_db.T,
    										"kos", 
    										"equilibrium"	);
    
    em_data jd_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
											SS_ref_db.P,
											SS_ref_db.T,
    										"jd", 
    										"equilibrium"	);
    
    em_data ru_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
											SS_ref_db.P,
											SS_ref_db.T,
    										"ru", 
    										"equilibrium"	);
    
    em_data cor_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
											SS_ref_db.P,
											SS_ref_db.T,
    										"cor", 
    										"equilibrium"	);
    
    em_data per_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
											SS_ref_db.P,
											SS_ref_db.T,
    										"per", 
    										"equilibrium"	);
    
    em_data acm_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
											SS_ref_db.P,
											SS_ref_db.T,
    										"acm", 
    										"equilibrium"	);
    
    SS_ref_db.gbase[0] 		= en_eq.gb;
    SS_ref_db.gbase[1] 		= fs_eq.gb;
    SS_ref_db.gbase[2] 		= 0.5*en_eq.gb + 0.5*fs_eq.gb - 6.6;
    SS_ref_db.gbase[3] 		= 0.005*SS_ref_db.P + di_eq.gb + 2.8;
    SS_ref_db.gbase[4] 		= mgts_eq.gb;
    SS_ref_db.gbase[5] 		= 0.05*SS_ref_db.P + 0.0155*SS_ref_db.T - jd_eq.gb + kos_eq.gb + mgts_eq.gb - 25.9;
    SS_ref_db.gbase[6] 		= -0.0061*SS_ref_db.P - 0.0051*SS_ref_db.T - 0.5*cor_eq.gb + 0.5*per_eq.gb + mgts_eq.gb + 0.5*ru_eq.gb - 5.0;
    SS_ref_db.gbase[7] 		= -0.089*SS_ref_db.P + acm_eq.gb - jd_eq.gb + mgts_eq.gb + 4.8;
    SS_ref_db.gbase[8] 		= jd_eq.gb + 18.8;
    
    SS_ref_db.ElShearMod[0] 	= en_eq.ElShearMod;
    SS_ref_db.ElShearMod[1] 	= fs_eq.ElShearMod;
    SS_ref_db.ElShearMod[2] 	= 0.5*en_eq.ElShearMod + 0.5*fs_eq.ElShearMod;
    SS_ref_db.ElShearMod[3] 	= di_eq.ElShearMod;
    SS_ref_db.ElShearMod[4] 	= mgts_eq.ElShearMod;
    SS_ref_db.ElShearMod[5] 	= -jd_eq.ElShearMod + kos_eq.ElShearMod + mgts_eq.ElShearMod;
    SS_ref_db.ElShearMod[6] 	= -0.5*cor_eq.ElShearMod + 0.5*per_eq.ElShearMod + mgts_eq.ElShearMod + 0.5*ru_eq.ElShearMod;
    SS_ref_db.ElShearMod[7] 	= acm_eq.ElShearMod - jd_eq.ElShearMod + mgts_eq.ElShearMod;
    SS_ref_db.ElShearMod[8] 	= jd_eq.ElShearMod;
    
    for (i = 0; i < len_ox; i++){
        SS_ref_db.Comp[0][i] 	= en_eq.C[i];
        SS_ref_db.Comp[1][i] 	= fs_eq.C[i];
        SS_ref_db.Comp[2][i] 	= 0.5*en_eq.C[i] + 0.5*fs_eq.C[i];
        SS_ref_db.Comp[3][i] 	= di_eq.C[i];
        SS_ref_db.Comp[4][i] 	= mgts_eq.C[i];
        SS_ref_db.Comp[5][i] 	= -jd_eq.C[i] + kos_eq.C[i] + mgts_eq.C[i];
        SS_ref_db.Comp[6][i] 	= -0.5*cor_eq.C[i] + 0.5*per_eq.C[i] + mgts_eq.C[i] + 0.5*ru_eq.C[i];
        SS_ref_db.Comp[7][i] 	= acm_eq.C[i] - jd_eq.C[i] + mgts_eq.C[i];
        SS_ref_db.Comp[8][i] 	= jd_eq.C[i];
    }
    
    
    for (i = 0; i < n_em; i++){
        SS_ref_db.z_em[i] = 1.0;
    };
    
    SS_ref_db.bounds_ref[0][0] = 0.0+eps;  SS_ref_db.bounds_ref[0][1] = 1.0-eps;
    SS_ref_db.bounds_ref[1][0] = 0.0+eps;  SS_ref_db.bounds_ref[1][1] = 2.0-eps;
    SS_ref_db.bounds_ref[2][0] = 0.0+eps;  SS_ref_db.bounds_ref[2][1] = 1.0-eps;
    SS_ref_db.bounds_ref[3][0] = -1.0+eps;  SS_ref_db.bounds_ref[3][1] = 1.0-eps;
    SS_ref_db.bounds_ref[4][0] = 0.0+eps;  SS_ref_db.bounds_ref[4][1] = 1.0-eps;
    SS_ref_db.bounds_ref[5][0] = 0.0+eps;  SS_ref_db.bounds_ref[5][1] = 1.0-eps;
    SS_ref_db.bounds_ref[6][0] = 0.0+eps;  SS_ref_db.bounds_ref[6][1] = 1.0-eps;
    SS_ref_db.bounds_ref[7][0] = 0.0+eps;  SS_ref_db.bounds_ref[7][1] = 1.0-eps;

	if (z_b.bulk_rock[9] == 0.){ 					
		SS_ref_db.z_em[5]          = 0.0;
        SS_ref_db.d_em[5]          = 1.0;
		SS_ref_db.bounds_ref[6][0] = 0.0; 
		SS_ref_db.bounds_ref[6][1] = 0.0;	
	}
	if (z_b.bulk_rock[8] == 0.){ 					
		SS_ref_db.z_em[7]          = 0.0;
        SS_ref_db.d_em[7]          = 1.0;
		SS_ref_db.bounds_ref[4][0] = 0.0; 
		SS_ref_db.bounds_ref[4][1] = 0.0;	
	}
	if (z_b.bulk_rock[7] == 0.){ 					
		SS_ref_db.z_em[6]          = 0.0;
        SS_ref_db.d_em[6]          = 1.0;
		SS_ref_db.bounds_ref[5][0] = 0.0; 
		SS_ref_db.bounds_ref[5][1] = 0.0;	
	}
    
    return SS_ref_db;
}

/**
  retrieve reference thermodynamic data for plagioclase4T (late 2021 update of TC, given by Eleanor)
*/
SS_ref G_SS_ig_fsp_function(SS_ref SS_ref_db, int EM_database, int len_ox, bulk_info z_b, double eps){
    
    int i, j;
    int n_em = SS_ref_db.n_em;
    
    char   *EM_tmp[] 		= {"ab","an","san"};
    for (int i = 0; i < SS_ref_db.n_em; i++){
        strcpy(SS_ref_db.EM_list[i],EM_tmp[i]);
    };
    int n_xeos = SS_ref_db.n_xeos;
    char   *CV_tmp[] 		= {"ca","k"};
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
        strcpy(SS_ref_db.CV_list[i],CV_tmp[i]);
    };
    SS_ref_db.W[0] = -0.04*SS_ref_db.P - 0.00935*SS_ref_db.T + 14.6;
    SS_ref_db.W[1] = 0.338*SS_ref_db.P - 0.00957*SS_ref_db.T + 24.1;
    SS_ref_db.W[2] = 48.5 - 0.13*SS_ref_db.P;
    
    SS_ref_db.v[0] = 0.674;
    SS_ref_db.v[1] = 0.55;
    SS_ref_db.v[2] = 1.0;
    
    
    em_data ab_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
											SS_ref_db.P,
											SS_ref_db.T,
    										"ab", 
    										"equilibrium"	);
    
    em_data an_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
											SS_ref_db.P,
											SS_ref_db.T,
    										"an", 
    										"equilibrium"	);
    
    em_data san_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
											SS_ref_db.P,
											SS_ref_db.T,
    										"san", 
    										"equilibrium"	);
    
    SS_ref_db.gbase[0] 		= ab_eq.gb;
    SS_ref_db.gbase[1] 		= an_eq.gb;
    SS_ref_db.gbase[2] 		= san_eq.gb;
    
    SS_ref_db.ElShearMod[0] 	= ab_eq.ElShearMod;
    SS_ref_db.ElShearMod[1] 	= an_eq.ElShearMod;
    SS_ref_db.ElShearMod[2] 	= san_eq.ElShearMod;
    
    for (i = 0; i < len_ox; i++){
        SS_ref_db.Comp[0][i] 	= ab_eq.C[i];
        SS_ref_db.Comp[1][i] 	= an_eq.C[i];
        SS_ref_db.Comp[2][i] 	= san_eq.C[i];
    }
    
    for (i = 0; i < n_em; i++){
        SS_ref_db.z_em[i] = 1.0;
    };
    
    SS_ref_db.bounds_ref[0][0] = 0.0+eps;  SS_ref_db.bounds_ref[0][1] = 1.0-eps;
    SS_ref_db.bounds_ref[1][0] = 0.0+eps;  SS_ref_db.bounds_ref[1][1] = 1.0-eps;
    
    return SS_ref_db;
}


/**
  retrieve reference thermodynamic data for spinel
*/
SS_ref G_SS_ig_spn_function(SS_ref SS_ref_db, int EM_database, int len_ox, bulk_info z_b, double eps){
    
    int i, j;
    int n_em = SS_ref_db.n_em;
    
    char   *EM_tmp[] 		= {"nsp","isp","nhc","ihc","nmt","imt","pcr","qndm"};
    for (int i = 0; i < SS_ref_db.n_em; i++){
        strcpy(SS_ref_db.EM_list[i],EM_tmp[i]);
    };
    int n_xeos = SS_ref_db.n_xeos;
    char   *CV_tmp[] 		= {"x","y","c","t","Q1","Q2","Q3"};
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
        strcpy(SS_ref_db.CV_list[i],CV_tmp[i]);
    };
    
    SS_ref_db.W[0] = -8.2;
    SS_ref_db.W[1] = 3.5;
    SS_ref_db.W[2] = -13.0;
    SS_ref_db.W[3] = 43.2;
    SS_ref_db.W[4] = 49.1;
    SS_ref_db.W[5] = -5.0;
    SS_ref_db.W[6] = 22.5;
    SS_ref_db.W[7] = 4.4;
    SS_ref_db.W[8] = -6.0;
    SS_ref_db.W[9] = 36.8;
    SS_ref_db.W[10] = 20.0;
    SS_ref_db.W[11] = 14.0;
    SS_ref_db.W[12] = 21.5;
    SS_ref_db.W[13] = -8.2;
    SS_ref_db.W[14] = 18.1;
    SS_ref_db.W[15] = 49.0;
    SS_ref_db.W[16] = -19.0;
    SS_ref_db.W[17] = 35.1;
    SS_ref_db.W[18] = -4.0;
    SS_ref_db.W[19] = 7.6;
    SS_ref_db.W[20] = -11.0;
    SS_ref_db.W[21] = 9.0;
    SS_ref_db.W[22] = 18.1;
    SS_ref_db.W[23] = 11.9;
    SS_ref_db.W[24] = 62.2;
    SS_ref_db.W[25] = -6.4;
    SS_ref_db.W[26] = 24.3;
    SS_ref_db.W[27] = 60.0;
    
    
    em_data sp_or 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
											SS_ref_db.P,
											SS_ref_db.T,
    										"sp", 
    										"ordered"	);
    
    em_data herc_or 	= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
											SS_ref_db.P,
											SS_ref_db.T,
    										"herc", 
    										"ordered"	);
    
    em_data mt_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
											SS_ref_db.P,
											SS_ref_db.T,
    										"mt", 
    										"equilibrium"	);
    
    em_data picr_eq 	= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
											SS_ref_db.P,
											SS_ref_db.T,
    										"picr", 
    										"equilibrium"	);
    
    em_data qnd_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
											SS_ref_db.P,
											SS_ref_db.T,
    										"qnd", 
    										"equilibrium"	);
    
    SS_ref_db.gbase[0] 		= sp_or.gb;
    SS_ref_db.gbase[1] 		= -0.005763*SS_ref_db.T + sp_or.gb + 23.6;
    SS_ref_db.gbase[2] 		= herc_or.gb;
    SS_ref_db.gbase[3] 		= -0.005763*SS_ref_db.T + herc_or.gb + 23.6;
    SS_ref_db.gbase[4] 		= 0.005763*SS_ref_db.T + mt_eq.gb;
    SS_ref_db.gbase[5] 		= mt_eq.gb + 0.3;
    SS_ref_db.gbase[6] 		= picr_eq.gb;
    SS_ref_db.gbase[7] 		= qnd_eq.gb - 30.0;
    
    SS_ref_db.ElShearMod[0] 	= sp_or.ElShearMod;
    SS_ref_db.ElShearMod[1] 	= sp_or.ElShearMod;
    SS_ref_db.ElShearMod[2] 	= herc_or.ElShearMod;
    SS_ref_db.ElShearMod[3] 	= herc_or.ElShearMod;
    SS_ref_db.ElShearMod[4] 	= mt_eq.ElShearMod;
    SS_ref_db.ElShearMod[5] 	= mt_eq.ElShearMod;
    SS_ref_db.ElShearMod[6] 	= picr_eq.ElShearMod;
    SS_ref_db.ElShearMod[7] 	= qnd_eq.ElShearMod;
    
    for (i = 0; i < len_ox; i++){
        SS_ref_db.Comp[0][i] 	= sp_or.C[i];
        SS_ref_db.Comp[1][i] 	= sp_or.C[i];
        SS_ref_db.Comp[2][i] 	= herc_or.C[i];
        SS_ref_db.Comp[3][i] 	= herc_or.C[i];
        SS_ref_db.Comp[4][i] 	= mt_eq.C[i];
        SS_ref_db.Comp[5][i] 	= mt_eq.C[i];
        SS_ref_db.Comp[6][i] 	= picr_eq.C[i];
        SS_ref_db.Comp[7][i] 	= qnd_eq.C[i];
    }
    
    for (i = 0; i < n_em; i++){
        SS_ref_db.z_em[i] = 1.0;
    };
    
    SS_ref_db.bounds_ref[0][0] = 0.0+eps;  SS_ref_db.bounds_ref[0][1] = 1.0-eps;
    SS_ref_db.bounds_ref[1][0] = 0.0+eps;  SS_ref_db.bounds_ref[1][1] = 1.0-eps;
    SS_ref_db.bounds_ref[2][0] = 0.0+eps;  SS_ref_db.bounds_ref[2][1] = 1.0-eps;
    SS_ref_db.bounds_ref[3][0] = 0.0+eps;  SS_ref_db.bounds_ref[3][1] = 1.0-eps;
    SS_ref_db.bounds_ref[4][0] = -1.0+eps;  SS_ref_db.bounds_ref[4][1] = 1.0-eps;
    SS_ref_db.bounds_ref[5][0] = -1.0+eps;  SS_ref_db.bounds_ref[5][1] = 1.0-eps;
    SS_ref_db.bounds_ref[6][0] = -1.0+eps;  SS_ref_db.bounds_ref[6][1] = 1.0-eps;

	if (z_b.bulk_rock[9] == 0.){ 
		SS_ref_db.z_em[6]          = 0.0;
		SS_ref_db.d_em[6]          = 1.0;
		SS_ref_db.bounds_ref[2][0] = 0.0; 
		SS_ref_db.bounds_ref[2][1] = 0.0;	
	}
	if (z_b.bulk_rock[8] == 0.){ 					
		SS_ref_db.z_em[4]          = 0.0;
        SS_ref_db.d_em[4]          = 1.0;
		SS_ref_db.z_em[5]          = 0.0;
        SS_ref_db.d_em[5]          = 1.0;
		SS_ref_db.bounds_ref[6][0] = 0.0; 
		SS_ref_db.bounds_ref[6][1] = 0.0;	
		SS_ref_db.bounds_ref[1][0] = 0.0; 
		SS_ref_db.bounds_ref[1][1] = 0.0;	
	}
	if (z_b.bulk_rock[7] == 0.){ 					
		SS_ref_db.z_em[7]          = 0.0;
        SS_ref_db.d_em[7]          = 1.0;
		SS_ref_db.bounds_ref[3][0] = 0.0; 
		SS_ref_db.bounds_ref[3][1] = 0.0;	
	}
    /* this lists the index of the order variables */
    SS_ref_db.orderVar      = 1;
    SS_ref_db.idOrderVar[4] = -1.0;
    SS_ref_db.idOrderVar[5] = -1.0;
    SS_ref_db.idOrderVar[6] = -1.0;


	return SS_ref_db;
}



/**************************************************************************************/
/**************************************************************************************/
/*********************Evan&Frost DATABASE (Evans&Frost , 2021)*************************/
/**************************************************************************************/
/**************************************************************************************/

/**
   retrieve reference thermodynamic data for ev_fluid
*/
SS_ref G_SS_um_fluid_function(SS_ref SS_ref_db, int EM_database, int len_ox, bulk_info z_b, double eps){
    
    int i, j;
    int n_em = SS_ref_db.n_em;

    char   *EM_tmp[] 		= {"H2","H2O"};
    for (int i = 0; i < SS_ref_db.n_em; i++){
        strcpy(SS_ref_db.EM_list[i],EM_tmp[i]);
    };
    int n_xeos = SS_ref_db.n_xeos;
    char   *CV_tmp[] 		= {"x"};
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
        strcpy(SS_ref_db.CV_list[i],CV_tmp[i]);
    };
    em_data H2_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"H2", 
    										"equilibrium"	);
    
    em_data H2O_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"H2O", 
    										"equilibrium"	);
    
    SS_ref_db.gbase[0] 		= H2_eq.gb;
    SS_ref_db.gbase[1] 		= H2O_eq.gb;

    SS_ref_db.ElShearMod[0] 	= H2_eq.ElShearMod;
    SS_ref_db.ElShearMod[1] 	= H2O_eq.ElShearMod;

    for (i = 0; i < len_ox; i++){
        SS_ref_db.Comp[0][i] 	= H2_eq.C[i];
        SS_ref_db.Comp[1][i] 	= H2O_eq.C[i];
    }

    for (i = 0; i < n_em; i++){
        SS_ref_db.z_em[i] = 1.0;
    };

    SS_ref_db.bounds_ref[0][0] = 0.0+eps;  SS_ref_db.bounds_ref[0][1] = 0.01-eps;

    return SS_ref_db;
}

/**
   retrieve reference thermodynamic data for ev_ol
*/
SS_ref G_SS_um_ol_function(SS_ref SS_ref_db, int EM_database, int len_ox, bulk_info z_b, double eps){
    
    int i, j;
    int n_em = SS_ref_db.n_em;
    
    char   *EM_tmp[] 		= {"fo","fa"};
    for (int i = 0; i < SS_ref_db.n_em; i++){
        strcpy(SS_ref_db.EM_list[i],EM_tmp[i]);
    };
    int n_xeos = SS_ref_db.n_xeos;
    char   *CV_tmp[] 		= {"x"};
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
        strcpy(SS_ref_db.CV_list[i],CV_tmp[i]);
    };
    SS_ref_db.W[0] = 8.0;
    
    
    em_data fo_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"fo", 
    										"equilibrium"	);
    
    em_data fa_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"fa", 
    										"equilibrium"	);
    
    SS_ref_db.gbase[0] 		= fo_eq.gb;
    SS_ref_db.gbase[1] 		= fa_eq.gb;
    
    SS_ref_db.ElShearMod[0] 	= fo_eq.ElShearMod;
    SS_ref_db.ElShearMod[1] 	= fa_eq.ElShearMod;
    
    for (i = 0; i < len_ox; i++){
        SS_ref_db.Comp[0][i] 	= fo_eq.C[i];
        SS_ref_db.Comp[1][i] 	= fa_eq.C[i];
    }
    
    for (i = 0; i < n_em; i++){
        SS_ref_db.z_em[i] = 1.0;
    };
    
    SS_ref_db.bounds_ref[0][0] = 0.0+eps;  SS_ref_db.bounds_ref[0][1] = 1.0-eps;

    return SS_ref_db;
}

/**
   retrieve reference thermodynamic data for ev_br
*/
SS_ref G_SS_um_br_function(SS_ref SS_ref_db, int EM_database, int len_ox, bulk_info z_b, double eps){
    
    int i, j;
    int n_em = SS_ref_db.n_em;
    
    char   *EM_tmp[] 		= {"br","fbr"};
    for (int i = 0; i < SS_ref_db.n_em; i++){
        strcpy(SS_ref_db.EM_list[i],EM_tmp[i]);
    };
    int n_xeos = SS_ref_db.n_xeos;
    char   *CV_tmp[] 		= {"x"};
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
        strcpy(SS_ref_db.CV_list[i],CV_tmp[i]);
    };
    em_data br_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"br", 
    										"equilibrium"	);
    
    em_data fo_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"fo", 
    										"equilibrium"	);
    
    em_data fa_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"fa", 
    										"equilibrium"	);
    
    SS_ref_db.gbase[0] 		= br_eq.gb;
    SS_ref_db.gbase[1] 		= br_eq.gb + fa_eq.gb/2.0 - fo_eq.gb/2.0 + 2.0;
    
    SS_ref_db.ElShearMod[0] 	= br_eq.ElShearMod;
    SS_ref_db.ElShearMod[1] 	= br_eq.ElShearMod + fa_eq.ElShearMod/2.0 - fo_eq.ElShearMod/2.0;
    
    for (i = 0; i < len_ox; i++){
        SS_ref_db.Comp[0][i] 	= br_eq.C[i];
        SS_ref_db.Comp[1][i] 	= br_eq.C[i] + fa_eq.C[i]/2.0 - fo_eq.C[i]/2.0;
    }
    
    for (i = 0; i < n_em; i++){
        SS_ref_db.z_em[i] = 1.0;
    };
    
    SS_ref_db.bounds_ref[0][0] = 0.0+eps;  SS_ref_db.bounds_ref[0][1] = 1.0-eps;

    return SS_ref_db;
}

/**
   retrieve reference thermodynamic data for ev_ch
*/
SS_ref G_SS_um_ch_function(SS_ref SS_ref_db, int EM_database, int len_ox, bulk_info z_b, double eps){
    
    int i, j;
    int n_em = SS_ref_db.n_em;
    
    char   *EM_tmp[] 		= {"chum","chuf"};
    for (int i = 0; i < SS_ref_db.n_em; i++){
        strcpy(SS_ref_db.EM_list[i],EM_tmp[i]);
    };
    int n_xeos = SS_ref_db.n_xeos;
    char   *CV_tmp[] 		= {"x"};
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
        strcpy(SS_ref_db.CV_list[i],CV_tmp[i]);
    };
    SS_ref_db.W[0] = 36.000;
    
    
    em_data chum_eq 		= get_em_data(	EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"chum", 
    										"equilibrium"	);
    
    em_data fo_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"fo", 
    										"equilibrium"	);
    
    em_data fa_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"fa", 
    										"equilibrium"	);
    
    SS_ref_db.gbase[0] 		= chum_eq.gb;
    SS_ref_db.gbase[1] 		= chum_eq.gb + 9.0*fa_eq.gb/2.0 - 9.0*fo_eq.gb/2.0 - 5.0;
    
    SS_ref_db.ElShearMod[0] 	= chum_eq.ElShearMod;
    SS_ref_db.ElShearMod[1] 	= chum_eq.ElShearMod + 9.0*fa_eq.ElShearMod/2.0 - 9.0*fo_eq.ElShearMod/2.0;
    
    for (i = 0; i < len_ox; i++){
        SS_ref_db.Comp[0][i] 	= chum_eq.C[i];
        SS_ref_db.Comp[1][i] 	= chum_eq.C[i] + 9.0*fa_eq.C[i]/2.0 - 9.0*fo_eq.C[i]/2.0;
    }
    
    for (i = 0; i < n_em; i++){
        SS_ref_db.z_em[i] = 1.0;
    };
    
    SS_ref_db.bounds_ref[0][0] = 0.0+eps;  SS_ref_db.bounds_ref[0][1] = 1.0-eps;

    return SS_ref_db;
}

/**
   retrieve reference thermodynamic data for ev_atg
*/
SS_ref G_SS_um_atg_function(SS_ref SS_ref_db, int EM_database, int len_ox, bulk_info z_b, double eps){
    
    int i, j;
    int n_em = SS_ref_db.n_em;
    
    char   *EM_tmp[] 		= {"atgf","fatg","atgo","aatg","oatg"};
    for (int i = 0; i < SS_ref_db.n_em; i++){
        strcpy(SS_ref_db.EM_list[i],EM_tmp[i]);
    };
    int n_xeos = SS_ref_db.n_xeos;
    char   *CV_tmp[] 		= {"x","y","f","t"};
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
        strcpy(SS_ref_db.CV_list[i],CV_tmp[i]);
    };
    SS_ref_db.W[0] = 12.000;
    SS_ref_db.W[1] = 4.0;
    SS_ref_db.W[2] = 10.000;
    SS_ref_db.W[3] = 8.0;
    SS_ref_db.W[4] = 8.0;
    SS_ref_db.W[5] = 15.000;
    SS_ref_db.W[6] = 13.600;
    SS_ref_db.W[7] = 7.0;
    SS_ref_db.W[8] = 5.6000;
    SS_ref_db.W[9] = 2.0;
    
    
    em_data atg_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"atg", 
    										"equilibrium"	);
    
    em_data fta_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"fta", 
    										"equilibrium"	);
    
    em_data ta_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"ta", 
    										"equilibrium"	);
    
    em_data tats_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"tats", 
    										"equilibrium"	);
    
    em_data andr_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"andr", 
    										"equilibrium"	);
    
    em_data gr_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"gr", 
    										"equilibrium"	);
    
    SS_ref_db.gbase[0] 		= atg_eq.gb/16.0;
    SS_ref_db.gbase[1] 		= atg_eq.gb/16.0 + fta_eq.gb - ta_eq.gb - 28.0;
    SS_ref_db.gbase[2] 		= atg_eq.gb/16.0 + fta_eq.gb/3.0 - ta_eq.gb/3.0 - 11.33;
    SS_ref_db.gbase[3] 		= atg_eq.gb/16.0 - ta_eq.gb + tats_eq.gb - 36.0;
    SS_ref_db.gbase[4] 		= andr_eq.gb/2.0 + atg_eq.gb/16.0 - gr_eq.gb/2.0 - ta_eq.gb + tats_eq.gb - 5.0;
    
    SS_ref_db.ElShearMod[0] 	= atg_eq.ElShearMod/16.0;
    SS_ref_db.ElShearMod[1] 	= atg_eq.ElShearMod/16.0 + fta_eq.ElShearMod - ta_eq.ElShearMod;
    SS_ref_db.ElShearMod[2] 	= atg_eq.ElShearMod/16.0 + fta_eq.ElShearMod/3.0 - ta_eq.ElShearMod/3.0;
    SS_ref_db.ElShearMod[3] 	= atg_eq.ElShearMod/16.0 - ta_eq.ElShearMod + tats_eq.ElShearMod;
    SS_ref_db.ElShearMod[4] 	= andr_eq.ElShearMod/2.0 + atg_eq.ElShearMod/16.0 - gr_eq.ElShearMod/2.0 - ta_eq.ElShearMod + tats_eq.ElShearMod;
    
    for (i = 0; i < len_ox; i++){
        SS_ref_db.Comp[0][i] 	= atg_eq.C[i]/16.0;
        SS_ref_db.Comp[1][i] 	= atg_eq.C[i]/16.0 + fta_eq.C[i] - ta_eq.C[i];
        SS_ref_db.Comp[2][i] 	= atg_eq.C[i]/16.0 + fta_eq.C[i]/3.0 - ta_eq.C[i]/3.0;
        SS_ref_db.Comp[3][i] 	= atg_eq.C[i]/16.0 - ta_eq.C[i] + tats_eq.C[i];
        SS_ref_db.Comp[4][i] 	= andr_eq.C[i]/2.0 + atg_eq.C[i]/16.0 - gr_eq.C[i]/2.0 - ta_eq.C[i] + tats_eq.C[i];
    }
    
    for (i = 0; i < n_em; i++){
        SS_ref_db.z_em[i] = 1.0;
    };
    
    SS_ref_db.bounds_ref[0][0] = 0.0+eps;  SS_ref_db.bounds_ref[0][1] = 1.0-eps;
    SS_ref_db.bounds_ref[1][0] = 0.0+eps;  SS_ref_db.bounds_ref[1][1] = 1.0-eps;
    SS_ref_db.bounds_ref[2][0] = 0.0+eps;  SS_ref_db.bounds_ref[2][1] = 1.0-eps;
    SS_ref_db.bounds_ref[3][0] = -1.0+eps;  SS_ref_db.bounds_ref[3][1] = 1.0-eps;

    return SS_ref_db;
}

/**
   retrieve reference thermodynamic data for ev_g
*/
SS_ref G_SS_um_g_function(SS_ref SS_ref_db, int EM_database, int len_ox, bulk_info z_b, double eps){
    
    int i, j;
    int n_em = SS_ref_db.n_em;
    
    char   *EM_tmp[] 		= {"py","alm"};
    for (int i = 0; i < SS_ref_db.n_em; i++){
        strcpy(SS_ref_db.EM_list[i],EM_tmp[i]);
    };
    int n_xeos = SS_ref_db.n_xeos;
    char   *CV_tmp[] 		= {"x"};
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
        strcpy(SS_ref_db.CV_list[i],CV_tmp[i]);
    };
    SS_ref_db.W[0] = 0.1*SS_ref_db.P + 4.0;
    
    
    em_data py_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"py", 
    										"equilibrium"	);
    
    em_data alm_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"alm", 
    										"equilibrium"	);
    
    SS_ref_db.gbase[0] 		= py_eq.gb;
    SS_ref_db.gbase[1] 		= alm_eq.gb;
    
    SS_ref_db.ElShearMod[0] 	= py_eq.ElShearMod;
    SS_ref_db.ElShearMod[1] 	= alm_eq.ElShearMod;
    
    for (i = 0; i < len_ox; i++){
        SS_ref_db.Comp[0][i] 	= py_eq.C[i];
        SS_ref_db.Comp[1][i] 	= alm_eq.C[i];
    }
    
    for (i = 0; i < n_em; i++){
        SS_ref_db.z_em[i] = 1.0;
    };
    
    //!!! RANGE MISSING, WAIT FOR UPDATE !!!
    SS_ref_db.bounds_ref[0][0] = 0.0+eps;  SS_ref_db.bounds_ref[0][1] = 1.0-eps;

    return SS_ref_db;
}

/**
   retrieve reference thermodynamic data for ev_ta
*/
SS_ref G_SS_um_ta_function(SS_ref SS_ref_db, int EM_database, int len_ox, bulk_info z_b, double eps){
    
    int i, j;
    int n_em = SS_ref_db.n_em;
    
    char   *EM_tmp[] 		= {"ta","fta","tao","tats","ota","tap"};
    for (int i = 0; i < SS_ref_db.n_em; i++){
        strcpy(SS_ref_db.EM_list[i],EM_tmp[i]);
    };
    int n_xeos = SS_ref_db.n_xeos;
    char   *CV_tmp[] 		= {"x","y","f","v","Q"};
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
        strcpy(SS_ref_db.CV_list[i],CV_tmp[i]);
    };
    SS_ref_db.W[0] = 12.000;
    SS_ref_db.W[1] = 8.0;
    SS_ref_db.W[2] = 10.000;
    SS_ref_db.W[3] = 9.5000;
    SS_ref_db.W[4] = 55.000;
    SS_ref_db.W[5] = 4.0;
    SS_ref_db.W[6] = 16.500;
    SS_ref_db.W[7] = 16.300;
    SS_ref_db.W[8] = 43.000;
    SS_ref_db.W[9] = 12.500;
    SS_ref_db.W[10] = 12.300;
    SS_ref_db.W[11] = 52.000;
    SS_ref_db.W[12] = 0.50;
    SS_ref_db.W[13] = 65.000;
    SS_ref_db.W[14] = 66.500;
    
    
    em_data ta_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"ta", 
    										"equilibrium"	);
    
    em_data fta_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"fta", 
    										"equilibrium"	);
    
    em_data tats_eq 		= get_em_data(	EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"tats", 
    										"equilibrium"	);
    
    em_data andr_eq 		= get_em_data(	EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"andr", 
    										"equilibrium"	);
    
    em_data gr_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"gr", 
    										"equilibrium"	);
    
    em_data tap_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"tap", 
    										"equilibrium"	);
    
    SS_ref_db.gbase[0] 		= ta_eq.gb;
    SS_ref_db.gbase[1] 		= fta_eq.gb;
    SS_ref_db.gbase[2] 		= 2.0*fta_eq.gb/3.0 + ta_eq.gb/3.0 - 2.0;
    SS_ref_db.gbase[3] 		= tats_eq.gb;
    SS_ref_db.gbase[4] 		= andr_eq.gb/2.0 - gr_eq.gb/2.0 + tats_eq.gb + 4.0;
    SS_ref_db.gbase[5] 		= tap_eq.gb;
    
    SS_ref_db.ElShearMod[0] 	= ta_eq.ElShearMod;
    SS_ref_db.ElShearMod[1] 	= fta_eq.ElShearMod;
    SS_ref_db.ElShearMod[2] 	= 2.0*fta_eq.ElShearMod/3.0 + ta_eq.ElShearMod/3.0;
    SS_ref_db.ElShearMod[3] 	= tats_eq.ElShearMod;
    SS_ref_db.ElShearMod[4] 	= andr_eq.ElShearMod/2.0 - gr_eq.ElShearMod/2.0 + tats_eq.ElShearMod;
    SS_ref_db.ElShearMod[5] 	= tap_eq.ElShearMod;
    
    for (i = 0; i < len_ox; i++){
        SS_ref_db.Comp[0][i] 	= ta_eq.C[i];
        SS_ref_db.Comp[1][i] 	= fta_eq.C[i];
        SS_ref_db.Comp[2][i] 	= 2.0*fta_eq.C[i]/3.0 + ta_eq.C[i]/3.0;
        SS_ref_db.Comp[3][i] 	= tats_eq.C[i];
        SS_ref_db.Comp[4][i] 	= andr_eq.C[i]/2.0 - gr_eq.C[i]/2.0 + tats_eq.C[i];
        SS_ref_db.Comp[5][i] 	= tap_eq.C[i];
    }
    
    for (i = 0; i < n_em; i++){
        SS_ref_db.z_em[i] = 1.0;
    };
    
    SS_ref_db.bounds_ref[0][0] = 0.0+eps;  SS_ref_db.bounds_ref[0][1] = 1.0-eps;
    SS_ref_db.bounds_ref[1][0] = 0.0+eps;  SS_ref_db.bounds_ref[1][1] = 1.0-eps;
    SS_ref_db.bounds_ref[2][0] = 0.0+eps;  SS_ref_db.bounds_ref[2][1] = 1.0-eps;
    SS_ref_db.bounds_ref[3][0] = 0.0+eps;  SS_ref_db.bounds_ref[3][1] = 1.0-eps;
    SS_ref_db.bounds_ref[4][0] = -1.0+eps; SS_ref_db.bounds_ref[4][1] = 1.0-eps;

    return SS_ref_db;
}

/**
   retrieve reference thermodynamic data for ev_chl
*/
SS_ref G_SS_um_chl_function(SS_ref SS_ref_db, int EM_database, int len_ox, bulk_info z_b, double eps){
    
    int i, j;
    int n_em = SS_ref_db.n_em;
    
    char   *EM_tmp[] 		= {"clin","afchl","ames","daph","ochl1","ochl4","f3clin"};
    for (int i = 0; i < SS_ref_db.n_em; i++){
        strcpy(SS_ref_db.EM_list[i],EM_tmp[i]);
    };
    int n_xeos = SS_ref_db.n_xeos;
    char   *CV_tmp[] 		= {"x","y","f","m","t","QA1","Q1","Q4"};
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
        strcpy(SS_ref_db.CV_list[i],CV_tmp[i]);
    };
    SS_ref_db.W[0] = 17.00;
    SS_ref_db.W[1] = 17.00;
    SS_ref_db.W[2] = 20.00;
    SS_ref_db.W[3] = 30.00;
    SS_ref_db.W[4] = 21.00;
    SS_ref_db.W[5] = 2.000;
    SS_ref_db.W[6] = 16.00;
    SS_ref_db.W[7] = 37.00;
    SS_ref_db.W[8] = 20.00;
    SS_ref_db.W[9] = 4.000;
    SS_ref_db.W[10] = 15.00;
    SS_ref_db.W[11] = 30.00;
    SS_ref_db.W[12] = 29.00;
    SS_ref_db.W[13] = 13.00;
    SS_ref_db.W[14] = 19.00;
    SS_ref_db.W[15] = 18.00;
    SS_ref_db.W[16] = 33.00;
    SS_ref_db.W[17] = 22.00;
    SS_ref_db.W[18] = 24.00;
    SS_ref_db.W[19] = 28.60;
    SS_ref_db.W[20] = 19.00;
    
    
    em_data clin_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"clin", 
    										"equilibrium"	);
    
    em_data afchl_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"afchl", 
    										"equilibrium"	);
    
    em_data ames_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"ames", 
    										"equilibrium"	);
    
    em_data daph_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"daph", 
    										"equilibrium"	);
    
    em_data gr_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"gr", 
    										"equilibrium"	);
    
    em_data andr_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"andr", 
    										"equilibrium"	);
    
    SS_ref_db.gbase[0] 		= clin_eq.gb;
    SS_ref_db.gbase[1] 		= afchl_eq.gb;
    SS_ref_db.gbase[2] 		= ames_eq.gb;
    SS_ref_db.gbase[3] 		= daph_eq.gb;
    SS_ref_db.gbase[4] 		= afchl_eq.gb - clin_eq.gb + daph_eq.gb + 3.0;
    SS_ref_db.gbase[5] 		= afchl_eq.gb - clin_eq.gb/5.0 + daph_eq.gb/5 + 2.4;
    SS_ref_db.gbase[6] 		= andr_eq.gb/2.0 + clin_eq.gb - gr_eq.gb/2.0 + 40;
    
    SS_ref_db.ElShearMod[0] 	= clin_eq.ElShearMod;
    SS_ref_db.ElShearMod[1] 	= afchl_eq.ElShearMod;
    SS_ref_db.ElShearMod[2] 	= ames_eq.ElShearMod;
    SS_ref_db.ElShearMod[3] 	= daph_eq.ElShearMod;
    SS_ref_db.ElShearMod[4] 	= afchl_eq.ElShearMod - clin_eq.ElShearMod + daph_eq.ElShearMod;
    SS_ref_db.ElShearMod[5] 	= afchl_eq.ElShearMod - clin_eq.ElShearMod/5.0 + daph_eq.ElShearMod/5.0;
    SS_ref_db.ElShearMod[6] 	= andr_eq.ElShearMod/2.0 + clin_eq.ElShearMod - gr_eq.ElShearMod/2.0;
    
    for (i = 0; i < len_ox; i++){
        SS_ref_db.Comp[0][i] 	= clin_eq.C[i];
        SS_ref_db.Comp[1][i] 	= afchl_eq.C[i];
        SS_ref_db.Comp[2][i] 	= ames_eq.C[i];
        SS_ref_db.Comp[3][i] 	= daph_eq.C[i];
        SS_ref_db.Comp[4][i] 	= afchl_eq.C[i] - clin_eq.C[i] + daph_eq.C[i];
        SS_ref_db.Comp[5][i] 	= afchl_eq.C[i] - clin_eq.C[i]/5.0 + daph_eq.C[i]/5.0;
        SS_ref_db.Comp[6][i] 	= andr_eq.C[i]/2.0 + clin_eq.C[i] - gr_eq.C[i]/2.0;
    }
    
    for (i = 0; i < n_em; i++){
        SS_ref_db.z_em[i] = 1.0;
    };
    
    SS_ref_db.bounds_ref[0][0] = 0.0+eps;  SS_ref_db.bounds_ref[0][1] = 1.0-eps;
    SS_ref_db.bounds_ref[1][0] = 0.0+eps;  SS_ref_db.bounds_ref[1][1] = 1.0-eps;
    SS_ref_db.bounds_ref[2][0] = 0.0+eps;  SS_ref_db.bounds_ref[2][1] = 1.0-eps;
    SS_ref_db.bounds_ref[3][0] = 0.0+eps;  SS_ref_db.bounds_ref[3][1] = 1.0-eps;
    SS_ref_db.bounds_ref[4][0] = -1.0+eps;  SS_ref_db.bounds_ref[4][1] = 1.0-eps;
    SS_ref_db.bounds_ref[5][0] = -1.0+eps;  SS_ref_db.bounds_ref[5][1] = 1.0-eps;
    
    return SS_ref_db;
}


/**
   retrieve reference thermodynamic data for ev_anth
*/
SS_ref G_SS_um_anth_function(SS_ref SS_ref_db, int EM_database, int len_ox, bulk_info z_b, double eps){
    
    int i, j;
    int n_em = SS_ref_db.n_em;
    
    char   *EM_tmp[] 		= {"anth","gedf","fant","a","b"};
    for (int i = 0; i < SS_ref_db.n_em; i++){
        strcpy(SS_ref_db.EM_list[i],EM_tmp[i]);
    };
    int n_xeos = SS_ref_db.n_xeos;
    char   *CV_tmp[] 		= {"x","y","z","a","c","Q1","Q2"};
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
        strcpy(SS_ref_db.CV_list[i],CV_tmp[i]);
    };
    SS_ref_db.W[0] = 25.00;
    SS_ref_db.W[1] = 33.00;
    SS_ref_db.W[2] = 18.00;
    SS_ref_db.W[3] = 23.00;
    SS_ref_db.W[4] = 39.50;
    SS_ref_db.W[5] = 29.00;
    SS_ref_db.W[6] = 34.60;
    SS_ref_db.W[7] = 12.00;
    SS_ref_db.W[8] = 8.000;
    SS_ref_db.W[9] = 20.00;
    
    SS_ref_db.v[0] = 1.000;
    SS_ref_db.v[1] = 1.500;
    SS_ref_db.v[2] = 1.000;
    SS_ref_db.v[3] = 1.000;
    SS_ref_db.v[4] = 1.000;
    
    
    em_data anth_eq 		= get_em_data(	EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"anth", 
    										"equilibrium"	);
    
    em_data ged_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"ged", 
    										"equilibrium"	);
    
    em_data fanth_eq 		= get_em_data(	EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"fanth", 
    										"equilibrium"	);
    
    SS_ref_db.gbase[0] 		= anth_eq.gb;
    SS_ref_db.gbase[1] 		= ged_eq.gb + 22.0;
    SS_ref_db.gbase[2] 		= fanth_eq.gb + 7.0;
    SS_ref_db.gbase[3] 		= 3.0*anth_eq.gb/7.0 + 4.0*fanth_eq.gb/7.0 - 5.5;
    SS_ref_db.gbase[4] 		= 2.0*anth_eq.gb/7.0 + 5.0*fanth_eq.gb/7.0 - 6.7;
    
    SS_ref_db.ElShearMod[0] 	= anth_eq.ElShearMod;
    SS_ref_db.ElShearMod[1] 	= ged_eq.ElShearMod;
    SS_ref_db.ElShearMod[2] 	= fanth_eq.ElShearMod;
    SS_ref_db.ElShearMod[3] 	= 3.0*anth_eq.ElShearMod/7.0 + 4.0*fanth_eq.ElShearMod/7.0;
    SS_ref_db.ElShearMod[4] 	= 2.0*anth_eq.ElShearMod/7.0 + 5.0*fanth_eq.ElShearMod/7.0;
    
    for (i = 0; i < len_ox; i++){
        SS_ref_db.Comp[0][i] 	= anth_eq.C[i];
        SS_ref_db.Comp[1][i] 	= ged_eq.C[i];
        SS_ref_db.Comp[2][i] 	= fanth_eq.C[i];
        SS_ref_db.Comp[3][i] 	= 3.0*anth_eq.C[i]/7.0 + 4.0*fanth_eq.C[i]/7.0;
        SS_ref_db.Comp[4][i] 	= 2.0*anth_eq.C[i]/7.0 + 5.0*fanth_eq.C[i]/7.0;
    }
    
    for (i = 0; i < n_em; i++){
        SS_ref_db.z_em[i] = 1.0;
    };
    
    SS_ref_db.bounds_ref[0][0] = 0.0+eps;  SS_ref_db.bounds_ref[0][1] = 1.0-eps;
    SS_ref_db.bounds_ref[1][0] = 0.0+eps;  SS_ref_db.bounds_ref[1][1] = 1.0-eps;
    SS_ref_db.bounds_ref[2][0] = -1.0+eps;  SS_ref_db.bounds_ref[2][1] = 1.0-eps;
    SS_ref_db.bounds_ref[3][0] = -1.0+eps;  SS_ref_db.bounds_ref[3][1] = 1.0-eps;
    
    return SS_ref_db;
}

/**
   retrieve reference thermodynamic data for ev_spi
*/
SS_ref G_SS_um_spi_function(SS_ref SS_ref_db, int EM_database, int len_ox, bulk_info z_b, double eps){
    
    int i, j;
    int n_em = SS_ref_db.n_em;
    
    char   *EM_tmp[] 		= {"herc","sp","mt"};
    for (int i = 0; i < SS_ref_db.n_em; i++){
        strcpy(SS_ref_db.EM_list[i],EM_tmp[i]);
    };
    int n_xeos = SS_ref_db.n_xeos;
    char   *CV_tmp[] 		= {"x","y"};
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
        strcpy(SS_ref_db.CV_list[i],CV_tmp[i]);
    };
    SS_ref_db.W[0] = 0.0;
    SS_ref_db.W[1] = 18.500;
    SS_ref_db.W[2] = 40.000;
    
    
    em_data herc_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"herc", 
    										"equilibrium"	);
    
    em_data sp_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"sp", 
    										"equilibrium"	);
    
    em_data mt_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"mt", 
    										"equilibrium"	);
    
    SS_ref_db.gbase[0] 		= herc_eq.gb;
    SS_ref_db.gbase[1] 		= sp_eq.gb;
    SS_ref_db.gbase[2] 		= mt_eq.gb;
    
    SS_ref_db.ElShearMod[0] 	= herc_eq.ElShearMod;
    SS_ref_db.ElShearMod[1] 	= sp_eq.ElShearMod;
    SS_ref_db.ElShearMod[2] 	= mt_eq.ElShearMod;
    
    for (i = 0; i < len_ox; i++){
        SS_ref_db.Comp[0][i] 	= herc_eq.C[i];
        SS_ref_db.Comp[1][i] 	= sp_eq.C[i];
        SS_ref_db.Comp[2][i] 	= mt_eq.C[i];
    }
    
    for (i = 0; i < n_em; i++){
        SS_ref_db.z_em[i] = 1.0;
    };
    
    SS_ref_db.bounds_ref[0][0] = 0.0+eps;  SS_ref_db.bounds_ref[0][1] = 1.0-eps;
    SS_ref_db.bounds_ref[1][0] = 0.0+eps;  SS_ref_db.bounds_ref[1][1] = 1.0-eps;

    return SS_ref_db;
}

/**
   retrieve reference thermodynamic data for ev_opx
*/
SS_ref G_SS_um_opx_function(SS_ref SS_ref_db, int EM_database, int len_ox, bulk_info z_b, double eps){
    
    int i, j;
    int n_em = SS_ref_db.n_em;
    
    char   *EM_tmp[] 		= {"en","fs","fm","mgts","fopx"};
    for (int i = 0; i < SS_ref_db.n_em; i++){
        strcpy(SS_ref_db.EM_list[i],EM_tmp[i]);
    };
    int n_xeos = SS_ref_db.n_xeos;
    char   *CV_tmp[] 		= {"x","y","f","Q"};
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
        strcpy(SS_ref_db.CV_list[i],CV_tmp[i]);
    };
    SS_ref_db.W[0] = 7.0;
    SS_ref_db.W[1] = 4.0;
    SS_ref_db.W[2] = 13.0 - 0.15*SS_ref_db.P;
    SS_ref_db.W[3] = 11.0 - 0.15*SS_ref_db.P;
    SS_ref_db.W[4] = 4.0;
    SS_ref_db.W[5] = 13.0 - 0.15*SS_ref_db.P;
    SS_ref_db.W[6] = 11.6 - 0.15*SS_ref_db.P;
    SS_ref_db.W[7] = 17.0 - 0.15*SS_ref_db.P;
    SS_ref_db.W[8] = 15.0 - 0.15*SS_ref_db.P;
    SS_ref_db.W[9] = 1.0;
    
    
    em_data en_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"en", 
    										"equilibrium"	);
    
    em_data fs_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"fs", 
    										"equilibrium"	);
    
    em_data mgts_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"mgts", 
    										"equilibrium"	);
    
    em_data andr_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"andr", 
    										"equilibrium"	);
    
    em_data gr_eq 		= get_em_data(		EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"gr", 
    										"equilibrium"	);
    
    SS_ref_db.gbase[0] 		= en_eq.gb;
    SS_ref_db.gbase[1] 		= fs_eq.gb;
    SS_ref_db.gbase[2] 		= en_eq.gb/2.0 + fs_eq.gb/2.0 - 6.6;
    SS_ref_db.gbase[3] 		= mgts_eq.gb;
    SS_ref_db.gbase[4] 		= andr_eq.gb/2.0 - gr_eq.gb/2.0 + mgts_eq.gb + 2.0;
    
    SS_ref_db.ElShearMod[0] 	= en_eq.ElShearMod;
    SS_ref_db.ElShearMod[1] 	= fs_eq.ElShearMod;
    SS_ref_db.ElShearMod[2] 	= en_eq.ElShearMod/2.0 + fs_eq.ElShearMod/2.0;
    SS_ref_db.ElShearMod[3] 	= mgts_eq.ElShearMod;
    SS_ref_db.ElShearMod[4] 	= andr_eq.ElShearMod/2.0 - gr_eq.ElShearMod/2.0 + mgts_eq.ElShearMod;
    
    for (i = 0; i < len_ox; i++){
        SS_ref_db.Comp[0][i] 	= en_eq.C[i];
        SS_ref_db.Comp[1][i] 	= fs_eq.C[i];
        SS_ref_db.Comp[2][i] 	= en_eq.C[i]/2.0 + fs_eq.C[i]/2.0;
        SS_ref_db.Comp[3][i] 	= mgts_eq.C[i];
        SS_ref_db.Comp[4][i] 	= andr_eq.C[i]/2.0 - gr_eq.C[i]/2.0 + mgts_eq.C[i];
    }
    
    for (i = 0; i < n_em; i++){
        SS_ref_db.z_em[i] = 1.0;
    };
    
    SS_ref_db.bounds_ref[0][0] = 0.0+eps;  SS_ref_db.bounds_ref[0][1] = 1.0-eps;
    SS_ref_db.bounds_ref[1][0] = 0.0+eps;  SS_ref_db.bounds_ref[1][1] = 1.0-eps;
    SS_ref_db.bounds_ref[2][0] = 0.0+eps;  SS_ref_db.bounds_ref[2][1] = 1.0-eps;
    SS_ref_db.bounds_ref[3][0] = 0.0+eps;  SS_ref_db.bounds_ref[3][1] = 1.0-eps;

    return SS_ref_db;
}

/**
   retrieve reference thermodynamic data for ev_po
*/
SS_ref G_SS_um_po_function(SS_ref SS_ref_db, int EM_database, int len_ox, bulk_info z_b, double eps){
    
    int i, j;
    int n_em = SS_ref_db.n_em;
    
    char   *EM_tmp[] 		= {"trov","trot"};
    for (int i = 0; i < SS_ref_db.n_em; i++){
        strcpy(SS_ref_db.EM_list[i],EM_tmp[i]);
    };
    int n_xeos = SS_ref_db.n_xeos;
    char   *CV_tmp[] 		= {"y"};
    for (int i = 0; i < SS_ref_db.n_xeos; i++){
        strcpy(SS_ref_db.CV_list[i],CV_tmp[i]);
    };
    SS_ref_db.W[0] = -3.190;
    
    
    em_data trov_eq 		= get_em_data(	EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"trov", 
    										"equilibrium"	);
    
    em_data trot_eq 		= get_em_data(	EM_database, 
    										len_ox,
    										z_b,
    										SS_ref_db.P,
    										SS_ref_db.T,
    										"trot", 
    										"equilibrium"	);
    
    SS_ref_db.gbase[0] 		= trov_eq.gb;
    SS_ref_db.gbase[1] 		= trot_eq.gb;
    
    SS_ref_db.ElShearMod[0] 	= trov_eq.ElShearMod;
    SS_ref_db.ElShearMod[1] 	= trot_eq.ElShearMod;
    
    for (i = 0; i < len_ox; i++){
        SS_ref_db.Comp[0][i] 	= trov_eq.C[i];
        SS_ref_db.Comp[1][i] 	= trot_eq.C[i];
    }
    
    for (i = 0; i < n_em; i++){
        SS_ref_db.z_em[i] = 1.0;
    };
    
    SS_ref_db.bounds_ref[0][0] = 0.0+eps;  SS_ref_db.bounds_ref[0][1] = 1.0-eps;
    
    return SS_ref_db;
}

SS_ref G_SS_mb_EM_function(		global_variable 	 gv,
								SS_ref 				 SS_ref_db, 
								int 				 EM_database, 
								bulk_info 	 		 z_b, 
								char   				*name				){
									  
	double eps 		   	= gv.bnd_val;
	double P 			= SS_ref_db.P;
	double T 			= SS_ref_db.T;	
					   
	SS_ref_db.ss_flags[0]  = 1;

	/* Associate the right solid-solution data */
	for (int FD = 0; FD < gv.n_Diff; FD++){				/* cycle twice in order to get gb_P_eps to calculate densities later on */
		
		if (FD == 8 || FD == 9){				// dG/dP0 to get Volume at P = 1bar
			SS_ref_db.P = 1.+ gv.gb_P_eps*gv.pdev[0][FD];
			SS_ref_db.T = T + gv.gb_T_eps*gv.pdev[1][FD];
		}
		else{
			SS_ref_db.P = P + gv.gb_P_eps*gv.pdev[0][FD];
			SS_ref_db.T = T + gv.gb_T_eps*gv.pdev[1][FD];
		}

        if (strcmp( name, "liq") == 0 ){
			/* turn of liquid when T < 600°C) */
			if ( T < gv.min_melt_T){
				SS_ref_db.ss_flags[0]  = 0;
			}
            SS_ref_db  = G_SS_mb_liq_function(SS_ref_db, EM_database, gv.len_ox, z_b, eps);	}
        else if (strcmp( name, "hb") == 0 ){
            SS_ref_db  = G_SS_mb_hb_function(SS_ref_db, EM_database, gv.len_ox, z_b, eps);	}
        else if (strcmp( name, "aug") == 0 ){
            SS_ref_db  = G_SS_mb_aug_function(SS_ref_db, EM_database, gv.len_ox, z_b, eps);	}
        else if (strcmp( name, "dio") == 0 ){
            SS_ref_db  = G_SS_mb_dio_function(SS_ref_db, EM_database, gv.len_ox, z_b, eps);	}
        else if (strcmp( name, "opx") == 0 ){
            SS_ref_db  = G_SS_mb_opx_function(SS_ref_db, EM_database, gv.len_ox, z_b, eps);	}
        else if (strcmp( name, "g") == 0 ){
            SS_ref_db  = G_SS_mb_g_function(SS_ref_db, EM_database, gv.len_ox, z_b, eps);	}
        else if (strcmp( name, "ol") == 0 ){
            SS_ref_db  = G_SS_mb_ol_function(SS_ref_db, EM_database, gv.len_ox, z_b, eps);	}
        else if (strcmp( name, "fsp") == 0 ){
            SS_ref_db  = G_SS_mb_fsp_function(SS_ref_db, EM_database, gv.len_ox, z_b, eps);	}
        else if (strcmp( name, "abc") == 0 ){
            SS_ref_db  = G_SS_mb_abc_function(SS_ref_db, EM_database, gv.len_ox, z_b, eps);	}
        else if (strcmp( name, "k4tr") == 0 ){
            SS_ref_db  = G_SS_mb_k4tr_function(SS_ref_db, EM_database, gv.len_ox, z_b, eps);	}
        else if (strcmp( name, "sp") == 0 ){
            SS_ref_db  = G_SS_mb_sp_function(SS_ref_db, EM_database, gv.len_ox, z_b, eps);	}
        else if (strcmp( name, "ilm") == 0 ){
			if (z_b.bulk_rock[gv.TiO2_id] == 0.){
				SS_ref_db.ss_flags[0]  = 0;
			}
            SS_ref_db  = G_SS_mb_ilm_function(SS_ref_db, EM_database, gv.len_ox, z_b, eps);	}
        else if (strcmp( name, "ilmm") == 0 ){
			if (z_b.bulk_rock[gv.TiO2_id] == 0.){
				SS_ref_db.ss_flags[0]  = 0;
			}
            SS_ref_db  = G_SS_mb_ilmm_function(SS_ref_db, EM_database, gv.len_ox, z_b, eps);	}
        else if (strcmp( name, "ep") == 0 ){
            SS_ref_db  = G_SS_mb_ep_function(SS_ref_db, EM_database, gv.len_ox, z_b, eps);	}
        else if (strcmp( name, "bi") == 0 ){
            SS_ref_db  = G_SS_mb_bi_function(SS_ref_db, EM_database, gv.len_ox, z_b, eps);	}
        else if (strcmp( name, "mu") == 0 ){
            SS_ref_db  = G_SS_mb_mu_function(SS_ref_db, EM_database, gv.len_ox, z_b, eps);	}
        else if (strcmp( name, "chl") == 0 ){
            SS_ref_db  = G_SS_mb_chl_function(SS_ref_db, EM_database, gv.len_ox, z_b, eps);	}
        else{
            printf("\nsolid solution '%s' is not in the database\n",name);	}

        for (int j = 0; j < SS_ref_db.n_em; j++){
            SS_ref_db.mu_array[FD][j] = SS_ref_db.gbase[j];
        }
    }

	for (int j = 0; j < SS_ref_db.n_xeos; j++){
		SS_ref_db.bounds[j][0] = SS_ref_db.bounds_ref[j][0];
		SS_ref_db.bounds[j][1] = SS_ref_db.bounds_ref[j][1];
	}

	/* Calculate the number of atoms in the bulk-rock composition */
	double fbc     = 0.0;
	for (int i = 0; i < gv.len_ox; i++){
		fbc += z_b.bulk_rock[i]*z_b.apo[i];
	}

	/* get the numer of atoms per endmember, needed to update normalization factor for liquid */
	for (int i = 0; i < SS_ref_db.n_em; i++){
		SS_ref_db.ape[i] = 0.0;
		for (int j = 0; j < gv.len_ox; j++){
			SS_ref_db.ape[i] += SS_ref_db.Comp[i][j]*z_b.apo[j];
		}
	}
	
	SS_ref_db.fbc = z_b.fbc;	
	
	if (gv.verbose == 1){
		printf(" %4s:",name);

		/* display Gibbs free energy of reference? */
		for (int j = 0; j < SS_ref_db.n_em; j++){
			printf(" %+12.5f",SS_ref_db.gbase[j]);
		}
		printf("\n");

		if (1 == 1){
			/* display molar composition */
            printf("\n S   A   C   M   F   K   N   T   O   H\n");
			for (int i = 0; i < SS_ref_db.n_em; i++){
				for (int j = 0; j < gv.len_ox; j++){
					printf(" %.1f",SS_ref_db.Comp[i][j]);
				}
				printf("\n");
			}
			printf("\n");
		}
	}

	return SS_ref_db;
};

/**
  checks if it can satisfy the mass constraint
*/
SS_ref G_SS_ig_EM_function(		global_variable 	 gv,
								SS_ref 				 SS_ref_db, 
								int 				 EM_database, 
								bulk_info 	 		 z_b, 
								char   				*name				){
									  
	double eps 		   	= gv.bnd_val;
	double P 			= SS_ref_db.P;
	double T 			= SS_ref_db.T;	
					   
	SS_ref_db.ss_flags[0]  = 1;

	/* Associate the right solid-solution data */
	for (int FD = 0; FD < gv.n_Diff; FD++){				/* cycle twice in order to get gb_P_eps to calculate densities later on */
		
		if (FD == 8 || FD == 9){				// dG/dP0 to get Volume at P = 1bar
			SS_ref_db.P = 1.+ gv.gb_P_eps*gv.pdev[0][FD];
			SS_ref_db.T = T + gv.gb_T_eps*gv.pdev[1][FD];
		}
		else{
			SS_ref_db.P = P + gv.gb_P_eps*gv.pdev[0][FD];
			SS_ref_db.T = T + gv.gb_T_eps*gv.pdev[1][FD];
		}

		if (strcmp( name, "bi") == 0 ){
			// if no H2O, deactivate
			if (z_b.bulk_rock[gv.H2O_id] == 0.){
				SS_ref_db.ss_flags[0]  = 0;
			}
			SS_ref_db  = G_SS_ig_bi_function(SS_ref_db, EM_database, gv.len_ox, z_b, eps);	}
		else if (strcmp( name, "cd") == 0){
			// if no H2O, deactivate
			if (z_b.bulk_rock[gv.H2O_id] == 0.){
				SS_ref_db.ss_flags[0]  = 0;
			}
			SS_ref_db  = G_SS_ig_cd_function(SS_ref_db, EM_database, gv.len_ox, z_b, eps);	}
		else if (strcmp( name, "cpx") == 0){
			SS_ref_db  = G_SS_ig_cpx_function(SS_ref_db, EM_database, gv.len_ox, z_b, eps);	}	
		else if (strcmp( name, "ep") == 0){
			// if no h2O, deactivate
			if (z_b.bulk_rock[gv.H2O_id] == 0.){
				SS_ref_db.ss_flags[0]  = 0;
			}
			SS_ref_db  = G_SS_ig_ep_function(SS_ref_db, EM_database, gv.len_ox, z_b, eps);	}
		else if (strcmp( name, "fl") == 0){
			// if no H2O, deactivate
			if (z_b.bulk_rock[gv.H2O_id] == 0.){
				SS_ref_db.ss_flags[0]  = 0;
			}
			SS_ref_db  = G_SS_ig_fl_function(SS_ref_db, EM_database, gv.len_ox, z_b, eps);	}		
		else if (strcmp( name, "g") == 0){
			SS_ref_db  = G_SS_ig_g_function(SS_ref_db, EM_database, gv.len_ox, z_b, eps);		}
		else if (strcmp( name, "hb") == 0){
			// if no H2O, deactivate
			if (z_b.bulk_rock[gv.H2O_id] == 0.){
				SS_ref_db.ss_flags[0]  = 0;
			}
			SS_ref_db  = G_SS_ig_hb_function(SS_ref_db, EM_database, gv.len_ox, z_b, eps);	}	
		else if (strcmp( name, "ilm") == 0){
			if (z_b.bulk_rock[gv.TiO2_id] == 0.){
				SS_ref_db.ss_flags[0]  = 0;
			}
			SS_ref_db  = G_SS_ig_ilm_function(SS_ref_db, EM_database, gv.len_ox, z_b, eps);	}
		else if (strcmp( name, "liq") == 0){
			/* turn of liquid when T < 600°C) */
			if ( T < gv.min_melt_T){
				SS_ref_db.ss_flags[0]  = 0;
			}
			SS_ref_db = G_SS_ig_liq_function(SS_ref_db, EM_database, gv.len_ox, z_b, eps);	}
		else if (strcmp( name, "mu") == 0){
			// if no H2O, deactivate
			if (z_b.bulk_rock[gv.H2O_id] == 0.){
				SS_ref_db.ss_flags[0]  = 0;
			}
			SS_ref_db  = G_SS_ig_mu_function(SS_ref_db, EM_database, gv.len_ox, z_b, eps);	}	
		else if (strcmp( name, "ol") == 0){
			SS_ref_db  = G_SS_ig_ol_function(SS_ref_db, EM_database, gv.len_ox, z_b, eps);	}
		else if (strcmp( name, "opx") == 0){
			SS_ref_db  = G_SS_ig_opx_function(SS_ref_db, EM_database, gv.len_ox, z_b, eps);	
        if (gv.limitCaOpx == 1){ SS_ref_db.bounds_ref[2][1] =  gv.CaOpxLim - eps;          }}
        else if (strcmp( name, "fper") == 0 ){
            SS_ref_db  = G_SS_ig_fper_function(SS_ref_db, EM_database, gv.len_ox, z_b, eps);}
		else if (strcmp( name, "fsp") == 0){
			SS_ref_db  = G_SS_ig_fsp_function(SS_ref_db, EM_database, gv.len_ox, z_b, eps);	}	
		else if (strcmp( name, "spn") == 0){
			SS_ref_db  = G_SS_ig_spn_function(SS_ref_db, EM_database, gv.len_ox, z_b, eps);	}
		else{
			printf("\nsolid solution '%s' is not in the database\n",name);	}	
		
		for (int j = 0; j < SS_ref_db.n_em; j++){
			SS_ref_db.mu_array[FD][j] = SS_ref_db.gbase[j];
		}
	}

	for (int j = 0; j < SS_ref_db.n_xeos; j++){
		SS_ref_db.bounds[j][0] = SS_ref_db.bounds_ref[j][0];
		SS_ref_db.bounds[j][1] = SS_ref_db.bounds_ref[j][1];
	}

	/* Calculate the number of atoms in the bulk-rock composition */
	double fbc     = 0.0;
	for (int i = 0; i < gv.len_ox; i++){
		fbc += z_b.bulk_rock[i]*z_b.apo[i];
	}

	/* get the numer of atoms per endmember, needed to update normalization factor for liquid */
	for (int i = 0; i < SS_ref_db.n_em; i++){
		SS_ref_db.ape[i] = 0.0;
		for (int j = 0; j < gv.len_ox; j++){
			SS_ref_db.ape[i] += SS_ref_db.Comp[i][j]*z_b.apo[j];
		}
	}
	
	SS_ref_db.fbc = z_b.fbc;	
	
	if (gv.verbose == 1){
		printf(" %4s:",name);

		/* display Gibbs free energy of reference? */
		for (int j = 0; j < SS_ref_db.n_em; j++){
			printf(" %+12.5f",SS_ref_db.gbase[j]);
		}
		printf("\n");

		if (1 == 1){
			/* display molar composition */
            printf("\n S   A   C   M   F   K   N   T   O   Cr  H\n");
			for (int i = 0; i < SS_ref_db.n_em; i++){
				for (int j = 0; j < gv.len_ox; j++){
					printf(" %.1f",SS_ref_db.Comp[i][j]);
				}
				printf("\n");
			}
			printf("\n");
		}
	}

	return SS_ref_db;
};

/**
  checks if it can satisfy the mass constraint
*/
SS_ref G_SS_mp_EM_function(		global_variable 	 gv,
								SS_ref 				 SS_ref_db, 
								int 				 EM_database, 
								bulk_info 	 		 z_b, 
								char   				*name				){
									  
	double eps 		   	= gv.bnd_val;
	double P 			= SS_ref_db.P;
	double T 			= SS_ref_db.T;	
					   
	SS_ref_db.ss_flags[0]  = 1;

	/* Associate the right solid-solution data */
	for (int FD = 0; FD < gv.n_Diff; FD++){				/* cycle twice in order to get gb_P_eps to calculate densities later on */
		
		if (FD == 8 || FD == 9){				// dG/dP0 to get Volume at P = 1bar
			SS_ref_db.P = 1.+ gv.gb_P_eps*gv.pdev[0][FD];
			SS_ref_db.T = T + gv.gb_T_eps*gv.pdev[1][FD];
		}
		else{
			SS_ref_db.P = P + gv.gb_P_eps*gv.pdev[0][FD];
			SS_ref_db.T = T + gv.gb_T_eps*gv.pdev[1][FD];
		}

		if (strcmp( name, "liq") == 0){
			/* turn of liquid when T < 600°C) */
			if ( T < gv.min_melt_T){
				SS_ref_db.ss_flags[0]  = 0;
			}
			SS_ref_db = G_SS_mp_liq_function(SS_ref_db, EM_database, gv.len_ox, z_b, eps);	    }
        else if (strcmp( name, "bi") == 0 ){
			// if no H2O, deactivate
			if (z_b.bulk_rock[gv.H2O_id] == 0.){
				SS_ref_db.ss_flags[0]  = 0;
			}
			SS_ref_db  = G_SS_mp_bi_function(SS_ref_db, EM_database, gv.len_ox, z_b, eps);	    }
		else if (strcmp( name, "fsp") == 0){
			SS_ref_db  = G_SS_mp_fsp_function(SS_ref_db, EM_database, gv.len_ox, z_b, eps);	}
		else if (strcmp( name, "g") == 0){
			SS_ref_db  = G_SS_mp_g_function(SS_ref_db, EM_database, gv.len_ox, z_b, eps);	}
        else if (strcmp( name, "ep") == 0 ){
			// if no H2O, deactivate
			if (z_b.bulk_rock[gv.H2O_id] == 0.){
				SS_ref_db.ss_flags[0]  = 0;
			}
			SS_ref_db  = G_SS_mp_ep_function(SS_ref_db, EM_database, gv.len_ox, z_b, eps);	    }
        else if (strcmp( name, "ma") == 0 ){
			// if no H2O, deactivate
			if (z_b.bulk_rock[gv.H2O_id] == 0.){
				SS_ref_db.ss_flags[0]  = 0;
			}
			SS_ref_db  = G_SS_mp_ma_function(SS_ref_db, EM_database, gv.len_ox, z_b, eps);	    }
        else if (strcmp( name, "mu") == 0 ){
			// if no H2O, deactivate
			if (z_b.bulk_rock[gv.H2O_id] == 0.){
				SS_ref_db.ss_flags[0]  = 0;
			}
			SS_ref_db  = G_SS_mp_mu_function(SS_ref_db, EM_database, gv.len_ox, z_b, eps);	    }
		else if (strcmp( name, "opx") == 0){
			SS_ref_db  = G_SS_mp_opx_function(SS_ref_db, EM_database, gv.len_ox, z_b, eps);	    }
		else if (strcmp( name, "sa") == 0){
			SS_ref_db  = G_SS_mp_sa_function(SS_ref_db, EM_database, gv.len_ox, z_b, eps);	    }
		else if (strcmp( name, "cd") == 0){
			// if no H2O, deactivate
			if (z_b.bulk_rock[gv.H2O_id] == 0.){
				SS_ref_db.ss_flags[0]  = 0;
			}
			SS_ref_db  = G_SS_mp_cd_function(SS_ref_db, EM_database, gv.len_ox, z_b, eps);	    }
        else if (strcmp( name, "st") == 0 ){
			// if no H2O, deactivate
			if (z_b.bulk_rock[gv.H2O_id] == 0.){
				SS_ref_db.ss_flags[0]  = 0;
			}
			SS_ref_db  = G_SS_mp_st_function(SS_ref_db, EM_database, gv.len_ox, z_b, eps);	    }
        else if (strcmp( name, "chl") == 0 ){
			// if no H2O, deactivate
			if (z_b.bulk_rock[gv.H2O_id] == 0.){
				SS_ref_db.ss_flags[0]  = 0;
			}
			SS_ref_db  = G_SS_mp_chl_function(SS_ref_db, EM_database, gv.len_ox, z_b, eps);	    }
        else if (strcmp( name, "ctd") == 0 ){
			// if no H2O, deactivate
			if (z_b.bulk_rock[gv.H2O_id] == 0.){
				SS_ref_db.ss_flags[0]  = 0;
			}
			SS_ref_db  = G_SS_mp_ctd_function(SS_ref_db, EM_database, gv.len_ox, z_b, eps);	    }
		else if (strcmp( name, "sp") == 0){
			SS_ref_db  = G_SS_mp_sp_function(SS_ref_db, EM_database, gv.len_ox, z_b, eps);	    }
		else if (strcmp( name, "ilm") == 0){
			if (z_b.bulk_rock[gv.TiO2_id] == 0.){
				SS_ref_db.ss_flags[0]  = 0;
			}
			SS_ref_db  = G_SS_mp_ilm_function(SS_ref_db, EM_database, gv.len_ox, z_b, eps);    }
		else if (strcmp( name, "ilmm") == 0){
			if (z_b.bulk_rock[gv.TiO2_id] == 0.){
				SS_ref_db.ss_flags[0]  = 0;
			}
			SS_ref_db  = G_SS_mp_ilmm_function(SS_ref_db, EM_database, gv.len_ox, z_b, eps);    }
		else if (strcmp( name, "mt") == 0){
			SS_ref_db  = G_SS_mp_mt_function(SS_ref_db, EM_database, gv.len_ox, z_b, eps);	    }
		else if (strcmp( name, "aq17") == 0){
			SS_ref_db  = G_SS_aq17_function(SS_ref_db, EM_database, gv.len_ox, z_b, eps);	    }
		else{
			printf("\nsolid solution '%s' is not in the database\n",name);	                    }	
		for (int j = 0; j < SS_ref_db.n_em; j++){
			SS_ref_db.mu_array[FD][j] = SS_ref_db.gbase[j];
			// printf(" %+10.10f",SS_ref_db.gbase[j]);
		}
		// printf("\n");
	}

	for (int j = 0; j < SS_ref_db.n_xeos; j++){
		SS_ref_db.bounds[j][0] = SS_ref_db.bounds_ref[j][0];
		SS_ref_db.bounds[j][1] = SS_ref_db.bounds_ref[j][1];
	}

	/* Calculate the number of atoms in the bulk-rock composition */
	double fbc     = 0.0;
	for (int i = 0; i < gv.len_ox; i++){
		fbc += z_b.bulk_rock[i]*z_b.apo[i];
	}

	/* get the numer of atoms per endmember, needed to update normalization factor for liquid */
	for (int i = 0; i < SS_ref_db.n_em; i++){
		SS_ref_db.ape[i] = 0.0;
		for (int j = 0; j < gv.len_ox; j++){
			SS_ref_db.ape[i] += SS_ref_db.Comp[i][j]*z_b.apo[j];
		}
	}
	
	SS_ref_db.fbc = z_b.fbc;	
	
	if (gv.verbose == 1){
		printf(" %4s:",name);

		/* display Gibbs free energy of reference? */
		for (int j = 0; j < SS_ref_db.n_em; j++){
			printf(" %+12.5f",SS_ref_db.gbase[j]);
		}
		printf("\n");

		if (1 == 1){
			/* display molar composition */
            printf("\n S   A   C   M   F   K   N   T   O   M   H  \n");
			for (int i = 0; i < SS_ref_db.n_em; i++){
				for (int j = 0; j < gv.len_ox; j++){
					printf(" %.1f",SS_ref_db.Comp[i][j]);
				}
				printf("\n");
			}
			printf("\n");
		}
	}

	return SS_ref_db;
};


/**
  checks if it can satisfy the mass constraint
*/
SS_ref G_SS_um_EM_function(		global_variable 	 gv,
								SS_ref 				 SS_ref_db, 
								int 				 EM_database, 
								bulk_info 	 		 z_b, 
								char   				*name				){
									  
	double eps 		   	= gv.bnd_val;
	double P 			= SS_ref_db.P;
	double T 			= SS_ref_db.T;	
					   
	SS_ref_db.ss_flags[0]  = 1;

	/* Associate the right solid-solution data */
	for (int FD = 0; FD < gv.n_Diff; FD++){				/* cycle twice in order to get gb_P_eps to calculate densities later on */
		//printf("%d %d %s...........\n",FD,EM_database,name);	
		if (FD == 8 || FD == 9){				// dG/dP0 to get Volume at P = 1bar
			SS_ref_db.P = 1.+ gv.gb_P_eps*gv.pdev[0][FD];
			SS_ref_db.T = T + gv.gb_T_eps*gv.pdev[1][FD];
		}
		else{
			SS_ref_db.P = P + gv.gb_P_eps*gv.pdev[0][FD];
			SS_ref_db.T = T + gv.gb_T_eps*gv.pdev[1][FD];
		}

		if (strcmp( name, "fl") == 0 ){
			// if no H2O, deactivate
			if (z_b.bulk_rock[gv.H2O_id] == 0.){
				SS_ref_db.ss_flags[0]  = 0;
			}
			SS_ref_db  = G_SS_um_fluid_function(SS_ref_db, EM_database, gv.len_ox, z_b, eps);}
		else if (strcmp( name, "br") == 0){
			// if no H2O, deactivate
			if (z_b.bulk_rock[gv.H2O_id] == 0.){
				SS_ref_db.ss_flags[0]  = 0;
			}
			SS_ref_db  = G_SS_um_br_function(SS_ref_db, EM_database, gv.len_ox, z_b, eps);	}
		else if (strcmp( name, "ol") == 0){
			SS_ref_db  = G_SS_um_ol_function(SS_ref_db, EM_database, gv.len_ox, z_b, eps);	}	
		else if (strcmp( name, "ch") == 0){
			// if no h2O, deactivate
			if (z_b.bulk_rock[gv.H2O_id] == 0.){
				SS_ref_db.ss_flags[0]  = 0;
			}
			SS_ref_db  = G_SS_um_ch_function(SS_ref_db, EM_database, gv.len_ox, z_b, eps);	}
		else if (strcmp( name, "atg") == 0){
			// if no H2O, deactivate
			if (z_b.bulk_rock[gv.H2O_id] == 0.){
				SS_ref_db.ss_flags[0]  = 0;
			}
			SS_ref_db  = G_SS_um_atg_function(SS_ref_db, EM_database, gv.len_ox, z_b, eps);	}		
		else if (strcmp( name, "g") == 0){
			SS_ref_db  = G_SS_um_g_function(SS_ref_db, EM_database, gv.len_ox, z_b, eps);		}
		else if (strcmp( name, "ta") == 0){
			// if no H2O, deactivate
			if (z_b.bulk_rock[gv.H2O_id] == 0.){
				SS_ref_db.ss_flags[0]  = 0;
			}
			SS_ref_db  = G_SS_um_ta_function(SS_ref_db, EM_database, gv.len_ox, z_b, eps);	}	
		else if (strcmp( name, "po") == 0){
			// if no S, deactivate
			SS_ref_db  = G_SS_um_po_function(SS_ref_db, EM_database, gv.len_ox, z_b, eps);	}
		else if (strcmp( name, "chl") == 0){
			// if no H2O, deactivate
			if (z_b.bulk_rock[gv.H2O_id] == 0.){
				SS_ref_db.ss_flags[0]  = 0;
			}
			SS_ref_db  = G_SS_um_chl_function(SS_ref_db, EM_database, gv.len_ox, z_b, eps);	}	
		else if (strcmp( name, "anth") == 0){
			// if no H2O, deactivate
			if (z_b.bulk_rock[gv.H2O_id] == 0.){
				SS_ref_db.ss_flags[0]  = 0;
			}
			SS_ref_db  = G_SS_um_anth_function(SS_ref_db, EM_database, gv.len_ox, z_b, eps);}
		else if (strcmp( name, "opx") == 0){
			SS_ref_db  = G_SS_um_opx_function(SS_ref_db, EM_database, gv.len_ox, z_b, eps);	}
		else if (strcmp( name, "spi") == 0){
			SS_ref_db  = G_SS_um_spi_function(SS_ref_db, EM_database, gv.len_ox, z_b, eps);	}
		else{
			printf("\nsolid solution '%s' is not in the database\n",name);	}	
		
		for (int j = 0; j < SS_ref_db.n_em; j++){
			SS_ref_db.mu_array[FD][j] = SS_ref_db.gbase[j];
			// printf(" %+10.10f",SS_ref_db.gbase[j]);
		}
		// printf("\n");
	}

	for (int j = 0; j < SS_ref_db.n_xeos; j++){
		SS_ref_db.bounds[j][0] = SS_ref_db.bounds_ref[j][0];
		SS_ref_db.bounds[j][1] = SS_ref_db.bounds_ref[j][1];
	}

	/* Calculate the number of atoms in the bulk-rock composition */
	double fbc     = 0.0;
	for (int i = 0; i < gv.len_ox; i++){
		fbc += z_b.bulk_rock[i]*z_b.apo[i];
	}

	/* get the numer of atoms per endmember, needed to update normalization factor for liquid */
	for (int i = 0; i < SS_ref_db.n_em; i++){
		SS_ref_db.ape[i] = 0.0;
		for (int j = 0; j < gv.len_ox; j++){
			SS_ref_db.ape[i] += SS_ref_db.Comp[i][j]*z_b.apo[j];
		}
	}

	SS_ref_db.fbc = z_b.fbc;	
	
	if (gv.verbose == 1){
		printf(" %4s:",name);

		/* display Gibbs free energy of reference? */
		for (int j = 0; j < SS_ref_db.n_em; j++){
			printf(" %+12.5f",SS_ref_db.gbase[j]);
		}
		printf("\n");

		if (1 == 1){
			/* display molar composition */
            if (EM_database == 4){
                printf("\n S   A   M   F   O   H   S\n");
            }
			for (int i = 0; i < SS_ref_db.n_em; i++){
				for (int j = 0; j < gv.len_ox; j++){
					printf(" %.1f",SS_ref_db.Comp[i][j]);
				}
				printf("\n");
			}
			printf("\n");
		}
	}

	return SS_ref_db;
};
