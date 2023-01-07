/**
Function to calculate the reference chemical potential of solid-solutions        

Holland et al., 2018 - Melting of peridotites through to granites                 
Igneous dataset to use with tc-ds633.txt                                         
"bi","cpx","cd","ep","fl","g","hb","ilm","liq","mu", "ol", "opx","pl4T","spn" 
*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <complex.h> 
#include <lapacke.h> 

#include "MAGEMin.h"
#include "gem_function.h"
#include "gss_function.h"
#include "NLopt_opt_function.h"
#include "simplex_levelling.h"
#include "toolkit.h"


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
typedef struct em_datas{
	double C[15];
	double ElShearMod;
	double gb;		
} em_data;

em_data get_em_data(	int 		 EM_database, 
						int          len_ox,
						bulk_info 	 z_b,
                        double       P,
                        double       T,
						char 		*name, 
						char 		*state		){

	em_data data; 
	PP_ref PP_db   		= G_EM_function(EM_database, len_ox, z_b.bulk_rock, z_b.apo, P, T, name, state);
   	data.ElShearMod  	= PP_db.phase_shearModulus;
   	data.gb  			= PP_db.gbase;

	for (int i = 0; i < len_ox; i++){
		data.C[i] = PP_db.Comp[i];
	}
	return data;
}


/**
  retrieve reference thermodynamic data for biotite 
*/
SS_ref G_SS_ig_bi_function(SS_ref SS_ref_db, int EM_database, int len_ox, bulk_info  z_b, double eps){	
								
	char   *EM_tmp[]   = {"phl","annm","obi","east","tbi","fbi"};	
	for (int i = 0; i < SS_ref_db.n_em; i++){ 
		strcpy(SS_ref_db.EM_list[i],EM_tmp[i]);			
	}
									  
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

    int i,j;
	
	int n_em = SS_ref_db.n_em;
	
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
	
	/* define endmembers to be deactivated when bulk-rock composition equal zero */
	for (i = 0; i < n_em; i++){
		SS_ref_db.z_em[i] = 1.0;
	}

	SS_ref_db.bounds_ref[0][0] = 0.+eps;	SS_ref_db.bounds_ref[0][1] = 1.-eps;
	SS_ref_db.bounds_ref[1][0] = 0.+eps;	SS_ref_db.bounds_ref[1][1] = 1.-eps;
	SS_ref_db.bounds_ref[2][0] = 0.+eps; 	SS_ref_db.bounds_ref[2][1] = 1.-eps;
	SS_ref_db.bounds_ref[3][0] = 0.+eps; 	SS_ref_db.bounds_ref[3][1] = 1.-eps;
	SS_ref_db.bounds_ref[4][0] = 0.+eps;	SS_ref_db.bounds_ref[4][1] = 1.-eps;	

	if (z_b.bulk_rock[7] == 0.){
		SS_ref_db.z_em[4]          = 0.0;
		SS_ref_db.bounds_ref[3][0] = eps; 
		SS_ref_db.bounds_ref[3][1] = eps;
	}
	if (z_b.bulk_rock[7] == 0.){
		SS_ref_db.z_em[5]          = 0.0;
		SS_ref_db.bounds_ref[2][0] = eps; 
		SS_ref_db.bounds_ref[2][1] = eps;
	}
	
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
    SS_ref_db.gbase[3] 		= cats_di.gb - 1.0*jd_eq.gb + kos_eq.gb - 4.9;
    SS_ref_db.gbase[4] 		= acm_eq.gb + cats_di.gb - 1.0*jd_eq.gb - 3.45;
    SS_ref_db.gbase[5] 		= -0.005*SS_ref_db.P - 0.0012*SS_ref_db.T - 0.5*cor_eq.gb + cats_di.gb + 0.5*per_eq.gb + 0.5*ru_eq.gb - 16.2;
    SS_ref_db.gbase[6] 		= jd_eq.gb;
    SS_ref_db.gbase[7] 		= 0.048*SS_ref_db.P - 0.002*SS_ref_db.T + en_eq.gb + 3.5;
    SS_ref_db.gbase[8] 		= 0.0465*SS_ref_db.P - 0.002*SS_ref_db.T + 0.5*en_eq.gb + 0.5*fs_eq.gb - 1.6;
    SS_ref_db.gbase[9] 		= 0.6*SS_ref_db.P - 1.0*abh_eq.gb + jd_eq.gb + san_eq.gb + 11.7;
    
    SS_ref_db.ElShearMod[0] 	= di_eq.ElShearMod;
    SS_ref_db.ElShearMod[1] 	= fs_eq.ElShearMod;
    SS_ref_db.ElShearMod[2] 	= cats_eq.ElShearMod;
    SS_ref_db.ElShearMod[3] 	= cats_di.ElShearMod - 1.0*jd_eq.ElShearMod + kos_eq.ElShearMod;
    SS_ref_db.ElShearMod[4] 	= acm_eq.ElShearMod + cats_di.ElShearMod - 1.0*jd_eq.ElShearMod;
    SS_ref_db.ElShearMod[5] 	= -0.5*cor_eq.ElShearMod + cats_di.ElShearMod + 0.5*per_eq.ElShearMod + 0.5*ru_eq.ElShearMod;
    SS_ref_db.ElShearMod[6] 	= jd_eq.ElShearMod;
    SS_ref_db.ElShearMod[7] 	= en_eq.ElShearMod;
    SS_ref_db.ElShearMod[8] 	= 0.5*en_eq.ElShearMod + 0.5*fs_eq.ElShearMod;
    SS_ref_db.ElShearMod[9] 	= -1.0*abh_eq.ElShearMod + jd_eq.ElShearMod + san_eq.ElShearMod;
    
    for (i = 0; i < len_ox; i++){
        SS_ref_db.Comp[0][i] 	= di_eq.C[i];
        SS_ref_db.Comp[1][i] 	= fs_eq.C[i];
        SS_ref_db.Comp[2][i] 	= cats_eq.C[i];
        SS_ref_db.Comp[3][i] 	= cats_di.C[i] - 1.0*jd_eq.C[i] + kos_eq.C[i];
        SS_ref_db.Comp[4][i] 	= acm_eq.C[i] + cats_di.C[i] - 1.0*jd_eq.C[i];
        SS_ref_db.Comp[5][i] 	= -0.5*cor_eq.C[i] + cats_di.C[i] + 0.5*per_eq.C[i] + 0.5*ru_eq.C[i];
        SS_ref_db.Comp[6][i] 	= jd_eq.C[i];
        SS_ref_db.Comp[7][i] 	= en_eq.C[i];
        SS_ref_db.Comp[8][i] 	= 0.5*en_eq.C[i] + 0.5*fs_eq.C[i];
        SS_ref_db.Comp[9][i] 	= -1.0*abh_eq.C[i] + jd_eq.C[i] + san_eq.C[i];
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
    
	if (z_b.bulk_rock[7] == 0.){	//TiO2
		SS_ref_db.z_em[5]          = 0.0;
		SS_ref_db.bounds_ref[7][0] = eps;  
		SS_ref_db.bounds_ref[7][1] = eps;
	}
	if (z_b.bulk_rock[8] == 0.){ 		//O
		SS_ref_db.z_em[4]          = 0.0;
		SS_ref_db.bounds_ref[5][0] = eps; 
		SS_ref_db.bounds_ref[5][1] = eps;
	}
	if (z_b.bulk_rock[9] == 0.){ 	//Cr2O3
		SS_ref_db.z_em[3]          = 0.0;
		SS_ref_db.bounds_ref[6][0] = eps;  
		SS_ref_db.bounds_ref[6][1] = eps;
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

	if (z_b.bulk_rock[10] == 0.){
		SS_ref_db.z_em[2]          = 0.0;
		SS_ref_db.bounds_ref[1][0] = eps;	
		SS_ref_db.bounds_ref[1][1] = eps;
	}

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
    
    return SS_ref_db;
}

/**
  retrieve reference thermodynamic data for fluid 
*/
SS_ref G_SS_ig_fl_function(SS_ref SS_ref_db, int EM_database, int len_ox, bulk_info z_b, double eps){
    
    int i, j;
    int n_em = SS_ref_db.n_em;
    
    char   *EM_tmp[] 		= {"qfL","slfL","wofL","fofL","fafL","jdfL","hmfL","ekfL","tifL","kjfL","h2o"};
    for (int i = 0; i < SS_ref_db.n_em; i++){
        strcpy(SS_ref_db.EM_list[i],EM_tmp[i]);
    };
    
    SS_ref_db.W[0] = 0;
    SS_ref_db.W[1] = 0;
    SS_ref_db.W[2] = 0;
    SS_ref_db.W[3] = 0;
    SS_ref_db.W[4] = 0;
    SS_ref_db.W[5] = 0;
    SS_ref_db.W[6] = 0;
    SS_ref_db.W[7] = 0;
    SS_ref_db.W[8] = 0;
    SS_ref_db.W[9] = 59.0 - 0.82*SS_ref_db.P;
    SS_ref_db.W[10] = 0;
    SS_ref_db.W[11] = 0;
    SS_ref_db.W[12] = 0;
    SS_ref_db.W[13] = 0;
    SS_ref_db.W[14] = 0;
    SS_ref_db.W[15] = 0;
    SS_ref_db.W[16] = 0;
    SS_ref_db.W[17] = 0;
    SS_ref_db.W[18] = 57.6 - 0.8*SS_ref_db.P;
    SS_ref_db.W[19] = 0;
    SS_ref_db.W[20] = 0;
    SS_ref_db.W[21] = 0;
    SS_ref_db.W[22] = 0;
    SS_ref_db.W[23] = 0;
    SS_ref_db.W[24] = 0;
    SS_ref_db.W[25] = 0;
    SS_ref_db.W[26] = 72.2 - 0.67*SS_ref_db.P;
    SS_ref_db.W[27] = 0;
    SS_ref_db.W[28] = 0;
    SS_ref_db.W[29] = 0;
    SS_ref_db.W[30] = 0;
    SS_ref_db.W[31] = 0;
    SS_ref_db.W[32] = 0;
    SS_ref_db.W[33] = 71.7 - 1.1*SS_ref_db.P;
    SS_ref_db.W[34] = 0;
    SS_ref_db.W[35] = 0;
    SS_ref_db.W[36] = 0;
    SS_ref_db.W[37] = 0;
    SS_ref_db.W[38] = 0;
    SS_ref_db.W[39] = 71.7 - 1.1*SS_ref_db.P;
    SS_ref_db.W[40] = 0;
    SS_ref_db.W[41] = 0;
    SS_ref_db.W[42] = 0;
    SS_ref_db.W[43] = 0;
    SS_ref_db.W[44] = 57.0 - 0.79*SS_ref_db.P;
    SS_ref_db.W[45] = 0;
    SS_ref_db.W[46] = 0;
    SS_ref_db.W[47] = 0;
    SS_ref_db.W[48] = 73.0 - 0.66*SS_ref_db.P;
    SS_ref_db.W[49] = 0;
    SS_ref_db.W[50] = 0;
    SS_ref_db.W[51] = 73.0 - 0.66*SS_ref_db.P;
    SS_ref_db.W[52] = 0;
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
    SS_ref_db.gbase[5] 		= -0.099*SS_ref_db.P + abL_eq.gb - 1.0*qL_eq.gb + 12.32;
    SS_ref_db.gbase[6] 		= -0.077*SS_ref_db.P + 0.5*hemL_eq.gb + 4.05;
    SS_ref_db.gbase[7] 		=  0.245*SS_ref_db.P + 0.5*eskL_eq.gb + 24.75;
    SS_ref_db.gbase[8] 		= -0.489*SS_ref_db.P + ruL_eq.gb + 5.6;
    SS_ref_db.gbase[9] 		= -0.227*SS_ref_db.P + kspL_eq.gb - 1.0*qL_eq.gb + 12.88;
    SS_ref_db.gbase[10] 		= H2O_eq.gb;
    
    SS_ref_db.ElShearMod[0] 	= 4.0*qL_eq.ElShearMod;
    SS_ref_db.ElShearMod[1] 	= silL_eq.ElShearMod;
    SS_ref_db.ElShearMod[2] 	= woL_eq.ElShearMod;
    SS_ref_db.ElShearMod[3] 	= 2.0*foL_eq.ElShearMod;
    SS_ref_db.ElShearMod[4] 	= 2.0*faL_eq.ElShearMod;
    SS_ref_db.ElShearMod[5] 	= abL_eq.ElShearMod - 1.0*qL_eq.ElShearMod;
    SS_ref_db.ElShearMod[6] 	= 0.5*hemL_eq.ElShearMod;
    SS_ref_db.ElShearMod[7] 	= 0.5*eskL_eq.ElShearMod;
    SS_ref_db.ElShearMod[8] 	= ruL_eq.ElShearMod;
    SS_ref_db.ElShearMod[9] 	= kspL_eq.ElShearMod - 1.0*qL_eq.ElShearMod;
    SS_ref_db.ElShearMod[10] 	= H2O_eq.ElShearMod;
    
    for (i = 0; i < len_ox; i++){
        SS_ref_db.Comp[0][i] 	= 4.0*qL_eq.C[i];
        SS_ref_db.Comp[1][i] 	= silL_eq.C[i];
        SS_ref_db.Comp[2][i] 	= woL_eq.C[i];
        SS_ref_db.Comp[3][i] 	= 2.0*foL_eq.C[i];
        SS_ref_db.Comp[4][i] 	= 2.0*faL_eq.C[i];
        SS_ref_db.Comp[5][i] 	= abL_eq.C[i] - 1.0*qL_eq.C[i];
        SS_ref_db.Comp[6][i] 	= 0.5*hemL_eq.C[i];
        SS_ref_db.Comp[7][i] 	= 0.5*eskL_eq.C[i];
        SS_ref_db.Comp[8][i] 	= ruL_eq.C[i];
        SS_ref_db.Comp[9][i] 	= kspL_eq.C[i] - 1.0*qL_eq.C[i];
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
    
	if (z_b.bulk_rock[7] == 0.){
		SS_ref_db.z_em[8]          = 0.0;
		SS_ref_db.bounds_ref[7][0] = eps;  
		SS_ref_db.bounds_ref[7][1] = eps;
	}
	if (z_b.bulk_rock[8] == 0.){
		SS_ref_db.z_em[6]          = 0.0;
		SS_ref_db.bounds_ref[5][0] = eps;  
		SS_ref_db.bounds_ref[5][1] = eps;
	}
	if (z_b.bulk_rock[9] == 0.){
		SS_ref_db.z_em[7]          = 0.0;
		SS_ref_db.bounds_ref[6][0] = eps;  
		SS_ref_db.bounds_ref[6][1] = eps;
	}
	if (z_b.bulk_rock[10] == 0.){
		SS_ref_db.z_em[10]         = 0.0;
		SS_ref_db.bounds_ref[9][0] = eps;  
		SS_ref_db.bounds_ref[9][1] = eps;	
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
    
	if (z_b.bulk_rock[8] == 0.){
		SS_ref_db.z_em[3]          = 0.0;
		SS_ref_db.bounds_ref[2][0] = eps;  
		SS_ref_db.bounds_ref[2][1] = eps;
	}
	if (z_b.bulk_rock[9] == 0.){
		SS_ref_db.z_em[4]          = 0.0;
		SS_ref_db.bounds_ref[3][0] = eps;  
		SS_ref_db.bounds_ref[3][1] = eps;
	}
	if (z_b.bulk_rock[7] == 0.){
		SS_ref_db.z_em[5]          = 0.0;
		SS_ref_db.bounds_ref[4][0] = eps;	
		SS_ref_db.bounds_ref[4][1] = eps;	
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
    SS_ref_db.gbase[6] 		= 0.428571428571429*cumm_eq.gb + 0.571428571428571*grun_eq.gb - 11.2;
    SS_ref_db.gbase[7] 		= 0.285714285714286*cumm_eq.gb + 0.714285714285714*grun_eq.gb - 13.8;
    SS_ref_db.gbase[8] 		= andr_eq.gb + gl_eq.gb - 1.0*gr_eq.gb;
    SS_ref_db.gbase[9] 		= 0.02*SS_ref_db.T + mu_eq.gb - 1.0*pa_eq.gb + parg_eq.gb - 7.06;
    SS_ref_db.gbase[10] 		= -2.0*dsp_eq.gb + 2.0*ru_eq.gb + ts_eq.gb + 95.0;
    
    SS_ref_db.ElShearMod[0] 	= tr_eq.ElShearMod;
    SS_ref_db.ElShearMod[1] 	= ts_eq.ElShearMod;
    SS_ref_db.ElShearMod[2] 	= parg_eq.ElShearMod;
    SS_ref_db.ElShearMod[3] 	= gl_eq.ElShearMod;
    SS_ref_db.ElShearMod[4] 	= cumm_eq.ElShearMod;
    SS_ref_db.ElShearMod[5] 	= grun_eq.ElShearMod;
    SS_ref_db.ElShearMod[6] 	= 0.428571428571429*cumm_eq.ElShearMod + 0.571428571428571*grun_eq.ElShearMod;
    SS_ref_db.ElShearMod[7] 	= 0.285714285714286*cumm_eq.ElShearMod + 0.714285714285714*grun_eq.ElShearMod;
    SS_ref_db.ElShearMod[8] 	= andr_eq.ElShearMod + gl_eq.ElShearMod - 1.0*gr_eq.ElShearMod;
    SS_ref_db.ElShearMod[9] 	= mu_eq.ElShearMod - 1.0*pa_eq.ElShearMod + parg_eq.ElShearMod;
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
        SS_ref_db.Comp[8][i] 	= andr_eq.C[i] + gl_eq.C[i] - 1.0*gr_eq.C[i];
        SS_ref_db.Comp[9][i] 	= mu_eq.C[i] - 1.0*pa_eq.C[i] + parg_eq.C[i];
        SS_ref_db.Comp[10][i] 	= -2.0*dsp_eq.C[i] + 2.0*ru_eq.C[i] + ts_eq.C[i];
    }
    
    for (i = 0; i < n_em; i++){
        SS_ref_db.z_em[i] = 1.0;
    };
    
    SS_ref_db.bounds_ref[0][0] = 0.0+eps;  SS_ref_db.bounds_ref[0][1] = 1.0-eps;
    SS_ref_db.bounds_ref[1][0] = 0.0+eps;  SS_ref_db.bounds_ref[1][1] = 2.0-eps;
    SS_ref_db.bounds_ref[2][0] = 0.0+eps;  SS_ref_db.bounds_ref[2][1] = 1.0-eps;
    SS_ref_db.bounds_ref[3][0] = 0.0+eps;  SS_ref_db.bounds_ref[3][1] = 1.0-eps;
    SS_ref_db.bounds_ref[4][0] = 0.0+eps;  SS_ref_db.bounds_ref[4][1] = 1.0-eps;
    SS_ref_db.bounds_ref[5][0] = 0.0+eps;  SS_ref_db.bounds_ref[5][1] = 1.0-eps;
    SS_ref_db.bounds_ref[6][0] = 0.0+eps;  SS_ref_db.bounds_ref[6][1] = 1.0-eps;
    SS_ref_db.bounds_ref[7][0] = 0.0+eps;  SS_ref_db.bounds_ref[7][1] = 1.0-eps;
    SS_ref_db.bounds_ref[8][0] = -1.0+eps;  SS_ref_db.bounds_ref[8][1] = 1.0-eps;
    SS_ref_db.bounds_ref[9][0] = -1.0+eps;  SS_ref_db.bounds_ref[9][1] = 1.0-eps;
    

	if (z_b.bulk_rock[7] == 0.){
		SS_ref_db.z_em[10]         = 0.0;
		SS_ref_db.bounds_ref[7][0] = eps;   
		SS_ref_db.bounds_ref[7][1] = eps;
	}
	if (z_b.bulk_rock[8] == 0.){
		SS_ref_db.z_em[8]          = 0.0;
		SS_ref_db.bounds_ref[6][0] = eps;   
		SS_ref_db.bounds_ref[6][1] = eps;
	}

	return SS_ref_db;	
}

/**
  retrieve reference thermodynamic data for ilmenite 
*/
SS_ref G_SS_ig_ilm_function(SS_ref SS_ref_db, int EM_database, int len_ox, bulk_info z_b, double eps){
    
    int i, j;
    int n_em = SS_ref_db.n_em;
    
    char   *EM_tmp[] 		= {"oilm","dilm","dhem"};
    for (int i = 0; i < SS_ref_db.n_em; i++){
        strcpy(SS_ref_db.EM_list[i],EM_tmp[i]);
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
    
    return SS_ref_db;
}

/**
  retrieve reference thermodynamic data for liquid (melt) 
*/
SS_ref G_SS_ig_liq_function(SS_ref SS_ref_db, int EM_database, int len_ox, bulk_info z_b, double eps){
    
    int i, j;
    int n_em = SS_ref_db.n_em;
    
    char   *EM_tmp[] 		= {"q4L","sl1L","wo1L","fo2L","fa2L","jdL","hmL","ekL","tiL","kjL","ctL","h2o1L"};
    for (int i = 0; i < SS_ref_db.n_em; i++){
        strcpy(SS_ref_db.EM_list[i],EM_tmp[i]);
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
    
    em_data kspL_eq 	= get_em_data(		EM_database, 
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
    SS_ref_db.gbase[5] 		= -0.089*SS_ref_db.P + abL_eq.gb - 1.0*qL_eq.gb + 12.19;
    SS_ref_db.gbase[6] 		= -0.032*SS_ref_db.P + 0.5*hemL_eq.gb + 3.3;
    SS_ref_db.gbase[7] 		= 0.245*SS_ref_db.P + 0.5*eskL_eq.gb + 24.85;
    SS_ref_db.gbase[8] 		= -0.489*SS_ref_db.P + ruL_eq.gb + 5.58;
    SS_ref_db.gbase[9] 		= -0.21*SS_ref_db.P + kspL_eq.gb - 1.0*qL_eq.gb + 11.98;
    SS_ref_db.gbase[10] 		= 0.053*SS_ref_db.P + 0.055*SS_ref_db.T - 1.0*qL_eq.gb + silL_eq.gb + woL_eq.gb - 108.3;
    SS_ref_db.gbase[11] 		= 0.00087*SS_ref_db.P - 0.0039*SS_ref_db.T + h2oL_eq.gb + 3.2;
    
    SS_ref_db.ElShearMod[0] 	= 4.0*qL_eq.ElShearMod;
    SS_ref_db.ElShearMod[1] 	= silL_eq.ElShearMod;
    SS_ref_db.ElShearMod[2] 	= woL_eq.ElShearMod;
    SS_ref_db.ElShearMod[3] 	= 2.0*foL_eq.ElShearMod;
    SS_ref_db.ElShearMod[4] 	= 2.0*faL_eq.ElShearMod;
    SS_ref_db.ElShearMod[5] 	= abL_eq.ElShearMod - 1.0*qL_eq.ElShearMod;
    SS_ref_db.ElShearMod[6] 	= 0.5*hemL_eq.ElShearMod;
    SS_ref_db.ElShearMod[7] 	= 0.5*eskL_eq.ElShearMod;
    SS_ref_db.ElShearMod[8] 	= ruL_eq.ElShearMod;
    SS_ref_db.ElShearMod[9] 	= kspL_eq.ElShearMod - 1.0*qL_eq.ElShearMod;
    SS_ref_db.ElShearMod[10] 	= -1.0*qL_eq.ElShearMod + silL_eq.ElShearMod + woL_eq.ElShearMod;
    SS_ref_db.ElShearMod[11] 	= h2oL_eq.ElShearMod;
    
    for (i = 0; i < len_ox; i++){
        SS_ref_db.Comp[0][i] 	= 4.0*qL_eq.C[i];
        SS_ref_db.Comp[1][i] 	= silL_eq.C[i];
        SS_ref_db.Comp[2][i] 	= woL_eq.C[i];
        SS_ref_db.Comp[3][i] 	= 2.0*foL_eq.C[i];
        SS_ref_db.Comp[4][i] 	= 2.0*faL_eq.C[i];
        SS_ref_db.Comp[5][i] 	= abL_eq.C[i] - 1.0*qL_eq.C[i];
        SS_ref_db.Comp[6][i] 	= 0.5*hemL_eq.C[i];
        SS_ref_db.Comp[7][i] 	= 0.5*eskL_eq.C[i];
        SS_ref_db.Comp[8][i] 	= ruL_eq.C[i];
        SS_ref_db.Comp[9][i] 	= kspL_eq.C[i] - 1.0*qL_eq.C[i];
        SS_ref_db.Comp[10][i] 	= -1.0*qL_eq.C[i] + silL_eq.C[i] + woL_eq.C[i];
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
    
	if (z_b.bulk_rock[7] == 0.){
		SS_ref_db.z_em[8]          = 0.0;
		SS_ref_db.bounds_ref[7][0] = eps; 
		SS_ref_db.bounds_ref[7][1] = eps;	
	}
	if (z_b.bulk_rock[8] == 0.){					// eps seems to give bad results, uses 0.0 instead
		SS_ref_db.z_em[6]          = 0.0;
		SS_ref_db.bounds_ref[5][0] = eps; 
		SS_ref_db.bounds_ref[5][1] = eps;	
	}
	if (z_b.bulk_rock[9] == 0.){
		SS_ref_db.z_em[7]          = 0.0;
		SS_ref_db.bounds_ref[6][0] = eps; 
		SS_ref_db.bounds_ref[6][1] = eps;	
	}
	if (z_b.bulk_rock[10] == 0.){ 					// no h2o, cannot be 0 for this xeos
		SS_ref_db.z_em[11]          = 0.0;
		SS_ref_db.bounds_ref[10][0] = eps; 
		SS_ref_db.bounds_ref[10][1] = eps;	
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
    SS_ref_db.gbase[4] 		= ma_eq.gb + 6.5;
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
    SS_ref_db.bounds_ref[2][0] = 0.0+eps;  SS_ref_db.bounds_ref[2][1] = 1.0-eps;
    
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
    SS_ref_db.gbase[5] 		= 0.05*SS_ref_db.P + 0.0155*SS_ref_db.T - 1.0*jd_eq.gb + kos_eq.gb + mgts_eq.gb - 25.9;
    SS_ref_db.gbase[6] 		= -0.0061*SS_ref_db.P - 0.0051*SS_ref_db.T - 0.5*cor_eq.gb + 0.5*per_eq.gb + mgts_eq.gb + 0.5*ru_eq.gb - 5.0;
    SS_ref_db.gbase[7] 		= -0.089*SS_ref_db.P + acm_eq.gb - 1.0*jd_eq.gb + mgts_eq.gb + 4.8;
    SS_ref_db.gbase[8] 		= jd_eq.gb + 18.8;
    
    SS_ref_db.ElShearMod[0] 	= en_eq.ElShearMod;
    SS_ref_db.ElShearMod[1] 	= fs_eq.ElShearMod;
    SS_ref_db.ElShearMod[2] 	= 0.5*en_eq.ElShearMod + 0.5*fs_eq.ElShearMod;
    SS_ref_db.ElShearMod[3] 	= di_eq.ElShearMod;
    SS_ref_db.ElShearMod[4] 	= mgts_eq.ElShearMod;
    SS_ref_db.ElShearMod[5] 	= -1.0*jd_eq.ElShearMod + kos_eq.ElShearMod + mgts_eq.ElShearMod;
    SS_ref_db.ElShearMod[6] 	= -0.5*cor_eq.ElShearMod + 0.5*per_eq.ElShearMod + mgts_eq.ElShearMod + 0.5*ru_eq.ElShearMod;
    SS_ref_db.ElShearMod[7] 	= acm_eq.ElShearMod - 1.0*jd_eq.ElShearMod + mgts_eq.ElShearMod;
    SS_ref_db.ElShearMod[8] 	= jd_eq.ElShearMod;
    
    for (i = 0; i < len_ox; i++){
        SS_ref_db.Comp[0][i] 	= en_eq.C[i];
        SS_ref_db.Comp[1][i] 	= fs_eq.C[i];
        SS_ref_db.Comp[2][i] 	= 0.5*en_eq.C[i] + 0.5*fs_eq.C[i];
        SS_ref_db.Comp[3][i] 	= di_eq.C[i];
        SS_ref_db.Comp[4][i] 	= mgts_eq.C[i];
        SS_ref_db.Comp[5][i] 	= -1.0*jd_eq.C[i] + kos_eq.C[i] + mgts_eq.C[i];
        SS_ref_db.Comp[6][i] 	= -0.5*cor_eq.C[i] + 0.5*per_eq.C[i] + mgts_eq.C[i] + 0.5*ru_eq.C[i];
        SS_ref_db.Comp[7][i] 	= acm_eq.C[i] - 1.0*jd_eq.C[i] + mgts_eq.C[i];
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
    
	if (z_b.bulk_rock[7] == 0.){ //tio2
		SS_ref_db.z_em[6]          = 0.0;
		SS_ref_db.bounds_ref[5][1] = eps;
		SS_ref_db.bounds_ref[5][0] = eps;
	}
	if (z_b.bulk_rock[8] == 0.){ //o
		SS_ref_db.z_em[7]          = 0.0;
		SS_ref_db.bounds_ref[4][0] = eps;  
		SS_ref_db.bounds_ref[4][1] = eps;	
	}
	if (z_b.bulk_rock[9] == 0.){ //cr2o3
		SS_ref_db.z_em[5]          = 0.0;
		SS_ref_db.bounds_ref[6][0] = eps;  
		SS_ref_db.bounds_ref[6][1] = eps;
	}

	return SS_ref_db;
}

/**
  retrieve reference thermodynamic data for plagioclase4T (late 2021 update of TC, given by Eleanor)
*/
SS_ref G_SS_ig_pl4T_function(SS_ref SS_ref_db, int EM_database, int len_ox, bulk_info z_b, double eps){
    
    int i, j;
    int n_em = SS_ref_db.n_em;
    
    char   *EM_tmp[] 		= {"ab","an","san"};
    for (int i = 0; i < SS_ref_db.n_em; i++){
        strcpy(SS_ref_db.EM_list[i],EM_tmp[i]);
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
    
	if (z_b.bulk_rock[7] == 0.){ //tio2
		SS_ref_db.z_em[7]          = 0.0;
		SS_ref_db.bounds_ref[3][0] = eps; 
		SS_ref_db.bounds_ref[3][1] = eps;
	}
	if (z_b.bulk_rock[8] == 0.){ //fe2o3
		SS_ref_db.bounds_ref[1][0] = eps;   
		SS_ref_db.bounds_ref[1][1] = eps;
		SS_ref_db.bounds_ref[6][0] = eps;   
		SS_ref_db.bounds_ref[6][1] = eps;
	}
	if (z_b.bulk_rock[9] == 0.){ //cr2o3
		SS_ref_db.z_em[6]          = 0.0;
		SS_ref_db.bounds_ref[2][0] = eps;  
		SS_ref_db.bounds_ref[2][1] = eps;
	}

	return SS_ref_db;
}

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
			if (z_b.bulk_rock[10] == 0.){
				SS_ref_db.ss_flags[0]  = 0;
			}
			SS_ref_db  = G_SS_ig_bi_function(SS_ref_db, EM_database, gv.len_ox, z_b, eps);	}
		else if (strcmp( name, "cd") == 0){
			// if no H2O, deactivate
			if (z_b.bulk_rock[10] == 0.){
				SS_ref_db.ss_flags[0]  = 0;
			}
			SS_ref_db  = G_SS_ig_cd_function(SS_ref_db, EM_database, gv.len_ox, z_b, eps);	}
		else if (strcmp( name, "cpx") == 0){
			SS_ref_db  = G_SS_ig_cpx_function(SS_ref_db, EM_database, gv.len_ox, z_b, eps);	}	
		else if (strcmp( name, "ep") == 0){
			// if no h2O, deactivate
			if (z_b.bulk_rock[8] == 0. || z_b.bulk_rock[10] == 0.){
				SS_ref_db.ss_flags[0]  = 0;
			}
			SS_ref_db  = G_SS_ig_ep_function(SS_ref_db, EM_database, gv.len_ox, z_b, eps);	}
		else if (strcmp( name, "fl") == 0){
			// if no H2O, deactivate
			if (z_b.bulk_rock[10] == 0.){
				SS_ref_db.ss_flags[0]  = 0;
			}
			SS_ref_db  = G_SS_ig_fl_function(SS_ref_db, EM_database, gv.len_ox, z_b, eps);	}		
		else if (strcmp( name, "g") == 0){
			SS_ref_db  = G_SS_ig_g_function(SS_ref_db, EM_database, gv.len_ox, z_b, eps);		}
		else if (strcmp( name, "hb") == 0){
			// if no H2O, deactivate
			if (z_b.bulk_rock[10] == 0.){
				SS_ref_db.ss_flags[0]  = 0;
			}
			SS_ref_db  = G_SS_ig_hb_function(SS_ref_db, EM_database, gv.len_ox, z_b, eps);	}	
		else if (strcmp( name, "ilm") == 0){
			// if no TiO2, deactivate
			if (z_b.bulk_rock[7] == 0.){
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
			if (z_b.bulk_rock[10] == 0.){
				SS_ref_db.ss_flags[0]  = 0;
			}
			SS_ref_db  = G_SS_ig_mu_function(SS_ref_db, EM_database, gv.len_ox, z_b, eps);	}	
		else if (strcmp( name, "ol") == 0){
			SS_ref_db  = G_SS_ig_ol_function(SS_ref_db, EM_database, gv.len_ox, z_b, eps);	}
		else if (strcmp( name, "opx") == 0){
			SS_ref_db  = G_SS_ig_opx_function(SS_ref_db, EM_database, gv.len_ox, z_b, eps);	}
		else if (strcmp( name, "pl4T") == 0){
			SS_ref_db  = G_SS_ig_pl4T_function(SS_ref_db, EM_database, gv.len_ox, z_b, eps);	}	
		else if (strcmp( name, "spn") == 0){
			SS_ref_db  = G_SS_ig_spn_function(SS_ref_db, EM_database, gv.len_ox, z_b, eps);	}
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
		for (int j = SS_ref_db.n_em; j < gv.len_ox+1; j++){
			printf("%13s","-");
		}
		printf("\n");

		if (1 == 1){
			/* display molar composition */
			for (int i = 0; i < SS_ref_db.n_em; i++){
				for (int j = 0; j < gv.len_ox; j++){
					printf(" %+10f",SS_ref_db.Comp[i][j]);
				}
				printf("\n");
			}
			printf("\n");
		}

	}

	return SS_ref_db;
};

