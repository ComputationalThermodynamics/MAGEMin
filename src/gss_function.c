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

#define nEl 11

typedef struct get_datas{
	double comp[nEl];
} get_data;

/** 
  function to easely get gb and comp in order to define solid solutions
*/
get_data get_gb_comp(	double 		*density, 
						double 		*gb_tmp,
						PP_ref 		 PP_db,
						get_data 	 data, 
						int 		 EM_database, 
						double 		*bulk_rock, 
						double 		 P, 
						double 		 T, 
						char 		*name, 
						char 		*state		){
					 
	PP_db  = G_EM_function(EM_database, bulk_rock, P, T, name, state);
   *gb_tmp = PP_db.gbase;

	for (int i = 0; i < nEl; i++){
		data.comp[i] = PP_db.Comp[i];
	}
	return data;
}

/**
  retrieve reference thermodynamic data for biotite 
*/
SS_ref G_SS_bi_function(SS_ref SS_ref_bi_db, int EM_database, double *bulk_rock, double P, 
								  double T, double eps){	
								
	char   *EM_tmp[] 		= {"phl","annm","obi","east","tbi","fbi"};	
	for (int i = 0; i < SS_ref_bi_db.n_em; i++){ 
		strcpy(SS_ref_bi_db.EM_list[i],EM_tmp[i]);			
	}
									  
    SS_ref_bi_db.W[0] = 12.;
    SS_ref_bi_db.W[1] = 4.;
    SS_ref_bi_db.W[2] = 10.;
    SS_ref_bi_db.W[3] = 30.;
    SS_ref_bi_db.W[4] = 8.;
    SS_ref_bi_db.W[5] = 8.;
    SS_ref_bi_db.W[6] = 5.;
    SS_ref_bi_db.W[7] = 32.;
    SS_ref_bi_db.W[8] = 13.6;
    SS_ref_bi_db.W[9] = 7.;
    SS_ref_bi_db.W[10] = 24.;
    SS_ref_bi_db.W[11] = 5.6;
    SS_ref_bi_db.W[12] = 40.0;
    SS_ref_bi_db.W[13] = 1.0;
    SS_ref_bi_db.W[14] = 40.0;			  	
									  								  			  
    PP_ref PP_db;

	double gb_tmp;
    double density;

    int i,j;
	
	int n_em = SS_ref_bi_db.n_em;
	
	get_data chem_comp1       = get_gb_comp(&density, &gb_tmp, PP_db, chem_comp1, EM_database, bulk_rock, P, T, "phl", "equilibrium");
	double gb1       = gb_tmp;
	SS_ref_bi_db.density[0] = density;
	
	get_data chem_comp2       = get_gb_comp(&density, &gb_tmp, PP_db, chem_comp2, EM_database, bulk_rock, P, T, "ann", "equilibrium");
	double gb_ann    = gb_tmp;	
	double gb2       = gb_ann - 6.;
	SS_ref_bi_db.density[1] = density;
	
	get_data chem_comp_phl    = get_gb_comp(&density, &gb_tmp, PP_db, chem_comp_phl, EM_database, bulk_rock, P, T, "phl", "equilibrium");
	double gb_phl_od  = gb_tmp;	
	double rho_phl_od = density;
	
	get_data chem_comp_ann    = get_gb_comp(&density, &gb_tmp, PP_db, chem_comp_ann, EM_database, bulk_rock, P, T, "ann", "equilibrium");
	double gb_ann_od  = gb_tmp;		
	double rho_ann_od = density;
	
	double gb3       = (gb_ann_od + gb_phl_od*2.)/3. -6.0;
	SS_ref_bi_db.density[2] = (rho_ann_od + rho_phl_od*2.)/3.;

	double chem_comp3[nEl];
	for (i = 0; i < nEl; i++){
		chem_comp3[i] = (chem_comp2.comp[i] + 2.0*chem_comp1.comp[i])/3.0;
	}
	
	get_data chem_comp4       = get_gb_comp(&density, &gb_tmp, PP_db, chem_comp4, EM_database, bulk_rock, P, T, "east", "equilibrium");
	double gb4       = gb_tmp;
	SS_ref_bi_db.density[3] = density;
	
	get_data chem_comp_br       = get_gb_comp(&density, &gb_tmp, PP_db, chem_comp_br, EM_database, bulk_rock, P, T, "br", "equilibrium");
	double gb_br       = gb_tmp;	
	double rho_br      = density;
	
	get_data chem_comp_ru       = get_gb_comp(&density, &gb_tmp, PP_db, chem_comp_ru, EM_database, bulk_rock, P, T, "ru", "equilibrium");
	double gb_ru       = gb_tmp;	
	double rho_ru      = density;	
	double gb5         = -gb_br + gb1 + gb_ru + 55.0;
	SS_ref_bi_db.density[4] = -rho_br + SS_ref_bi_db.density[0] + rho_ru;	
	
	double chem_comp5[nEl];
	for (i = 0; i < nEl; i++){
		chem_comp5[i] = - chem_comp_br.comp[i] + chem_comp1.comp[i] + chem_comp_ru.comp[i];
	}
	
	get_data chem_comp_gr       = get_gb_comp(&density, &gb_tmp, PP_db, chem_comp_gr, EM_database, bulk_rock, P, T, "gr", "equilibrium");
	double gb_gr       = gb_tmp;
	double rho_gr      = density;		
	
	get_data chem_comp_andr       = get_gb_comp(&density, &gb_tmp, PP_db, chem_comp_andr, EM_database, bulk_rock, P, T, "andr", "equilibrium");
	double gb_andr       = gb_tmp;
	double rho_andr      = density;
	
	double gb6         = (gb_andr + 2.0*gb4 - gb_gr)/2. - 3.0;
	SS_ref_bi_db.density[5] = (rho_andr + 2.0*SS_ref_bi_db.density[3] - rho_gr)/2.;
	double chem_comp6[nEl];
	for (i = 0; i < nEl; i++){
		chem_comp6[i] = (chem_comp_andr.comp[i] + 2.*chem_comp4.comp[i] - chem_comp_gr.comp[i])/2.0;
	}
	
	SS_ref_bi_db.gbase[0] = gb1;
	SS_ref_bi_db.gbase[1] = gb2;
	SS_ref_bi_db.gbase[2] = gb3;
	SS_ref_bi_db.gbase[3] = gb4;
	SS_ref_bi_db.gbase[4] = gb5;
	SS_ref_bi_db.gbase[5] = gb6;
	
	for (i = 0; i < nEl; i++){
		SS_ref_bi_db.Comp[0][i] = chem_comp1.comp[i];
		SS_ref_bi_db.Comp[1][i] = chem_comp2.comp[i];
		SS_ref_bi_db.Comp[2][i] = chem_comp3[i];
		SS_ref_bi_db.Comp[3][i] = chem_comp4.comp[i];
		SS_ref_bi_db.Comp[4][i] = chem_comp5[i];
		SS_ref_bi_db.Comp[5][i] = chem_comp6[i];
	}
	
	/* define endmembers to be deactivated when bulk-rock composition equal zero */
	for (i = 0; i < n_em; i++){
		SS_ref_bi_db.z_em[i] = 1.0;
	}

	SS_ref_bi_db.box_bounds_default[0][0] = 0.+eps;	SS_ref_bi_db.box_bounds_default[0][1] = 1.-eps;
	SS_ref_bi_db.box_bounds_default[1][0] = 0.+eps;	SS_ref_bi_db.box_bounds_default[1][1] = 1.-eps;
	SS_ref_bi_db.box_bounds_default[2][0] = 0.+eps; SS_ref_bi_db.box_bounds_default[2][1] = 1.-eps;
	SS_ref_bi_db.box_bounds_default[3][0] = 0.+eps; SS_ref_bi_db.box_bounds_default[3][1] = 1.-eps;
	SS_ref_bi_db.box_bounds_default[4][0] = 0.+eps;	SS_ref_bi_db.box_bounds_default[4][1] = 1.-eps;	

	if (bulk_rock[7] == 0.){
		SS_ref_bi_db.z_em[4]          = 0.0;
		SS_ref_bi_db.box_bounds_default[3][0] = eps; 
		SS_ref_bi_db.box_bounds_default[3][1] = eps;
	}
	if (bulk_rock[7] == 0.){
		SS_ref_bi_db.z_em[5]          = 0.0;
		SS_ref_bi_db.box_bounds_default[2][0] = eps; 
		SS_ref_bi_db.box_bounds_default[2][1] = eps;
	}
	
	return SS_ref_bi_db;
}

/**
  retrieve reference thermodynamic data for clinopyroxene 
*/
SS_ref G_SS_cpx_function(SS_ref SS_ref_cpx_db, int EM_database, double *bulk_rock, double P, 
								  double T, double eps){	
									  	
	char   *EM_tmp[] 		= {"di","cfs","cats","crdi","cess","cbuf","jd","cen","cfm","kjd"};	
	for (int i = 0; i < SS_ref_cpx_db.n_em; i++){ 
		strcpy(SS_ref_cpx_db.EM_list[i],EM_tmp[i]);			
	}
		
    SS_ref_cpx_db.W[0] = 25.8;
    SS_ref_cpx_db.W[1] = 13.0 - 0.06*P;
    SS_ref_cpx_db.W[2] = 8.0;
    SS_ref_cpx_db.W[3] = 8.0;
    SS_ref_cpx_db.W[4] = 8.0;
    SS_ref_cpx_db.W[5] = 26.0;
    SS_ref_cpx_db.W[6] = 29.8;
    SS_ref_cpx_db.W[7] = 20.6;
    SS_ref_cpx_db.W[8] = 26.0;
    SS_ref_cpx_db.W[9] = 25.0 - 0.10*P;
    SS_ref_cpx_db.W[10] = 38.3;
    SS_ref_cpx_db.W[11] = 43.3;
    SS_ref_cpx_db.W[12] = 24.0;
    SS_ref_cpx_db.W[13] = 24.0;
    SS_ref_cpx_db.W[14] = 2.3;
    SS_ref_cpx_db.W[15] = 3.5;
    SS_ref_cpx_db.W[16] = 24.0;
    SS_ref_cpx_db.W[17] = 2.0;
    SS_ref_cpx_db.W[18] = 2.0;
    SS_ref_cpx_db.W[19] = 6.0;
    SS_ref_cpx_db.W[20] = 6.0;
    SS_ref_cpx_db.W[21] = 45.2 - 0.35*P;
    SS_ref_cpx_db.W[22] = 27.0 - 0.10*P;
    SS_ref_cpx_db.W[23] = 6.0;
    SS_ref_cpx_db.W[24] = 2.0;
    SS_ref_cpx_db.W[25] = 6.0;
    SS_ref_cpx_db.W[26] = 3.0;
    SS_ref_cpx_db.W[27] = 52.3;
    SS_ref_cpx_db.W[28] = 40.3;
    SS_ref_cpx_db.W[29] = 3.0;
    SS_ref_cpx_db.W[30] = 6.0;
    SS_ref_cpx_db.W[31] = 3.0;
    SS_ref_cpx_db.W[32] = 57.3;
    SS_ref_cpx_db.W[33] = 45.3;
    SS_ref_cpx_db.W[34] = 3.0;
    SS_ref_cpx_db.W[35] = 16.0;
    SS_ref_cpx_db.W[36] = 24.0;
    SS_ref_cpx_db.W[37] = 22.0;
    SS_ref_cpx_db.W[38] = 16.0;
    SS_ref_cpx_db.W[39] = 40.0;
    SS_ref_cpx_db.W[40] = 40.0;
    SS_ref_cpx_db.W[41] = 28.0;
    SS_ref_cpx_db.W[42] = 4.0;
    SS_ref_cpx_db.W[43] = 40.0;
    SS_ref_cpx_db.W[44] = 40.0;
    SS_ref_cpx_db.v[0] = 1.20;
    SS_ref_cpx_db.v[1] = 1.00;
    SS_ref_cpx_db.v[2] = 1.90;
    SS_ref_cpx_db.v[3] = 1.90;
    SS_ref_cpx_db.v[4] = 1.90;
    SS_ref_cpx_db.v[5] = 1.90;
    SS_ref_cpx_db.v[6] = 1.20;
    SS_ref_cpx_db.v[7] = 1.00;
    SS_ref_cpx_db.v[8] = 1.00;
    SS_ref_cpx_db.v[9] = 1.20;

    PP_ref PP_db;

	double gb_tmp;
    double density;
    int i,j;
	
	int n_em  = SS_ref_cpx_db.n_em;

	get_data chem_comp1       = get_gb_comp(&density, &gb_tmp, PP_db, chem_comp1, EM_database, bulk_rock, P, T, "di", "equilibrium");
	double gb1       = gb_tmp;	
	SS_ref_cpx_db.density[0] = density;
	
	get_data chem_comp2       = get_gb_comp(&density, &gb_tmp, PP_db, chem_comp2, EM_database, bulk_rock, P, T, "fs", "equilibrium");
	double gb_fs     = gb_tmp;
	double gb2       = gb_fs + 2.1 - 0.002*T + 0.045*P;
	SS_ref_cpx_db.density[1] = density;

	get_data chem_comp3       = get_gb_comp(&density, &gb_tmp, PP_db, chem_comp3, EM_database, bulk_rock, P, T, "cats", "equilibrium");
	double gb3       = gb_tmp;	
	SS_ref_cpx_db.density[2] = density;	
	
	get_data chem_comp_cats_d = get_gb_comp(&density, &gb_tmp, PP_db, chem_comp_cats_d, EM_database, bulk_rock, P, T, "cats", "disordered");
	double gb_cats_d = gb_tmp;	
	double rho_cats_d = density;	
	
	get_data chem_comp7       = get_gb_comp(&density, &gb_tmp, PP_db, chem_comp7, EM_database, bulk_rock, P, T, "jd", "equilibrium");
	double gb7       = gb_tmp;	
	SS_ref_cpx_db.density[6] = density;

	get_data chem_comp_kos    = get_gb_comp(&density, &gb_tmp, PP_db, chem_comp_kos, EM_database, bulk_rock, P, T, "kos", "equilibrium");
	double gb_kos    = gb_tmp;	
	double rho_kos   = density;
	double gb4       = gb_cats_d + gb_kos - gb7 - 4.9;
	SS_ref_cpx_db.density[3] = rho_cats_d + rho_kos- SS_ref_cpx_db.density[6];
	
	double chem_comp4[nEl];	
	for (i = 0; i < nEl; i++){
		chem_comp4[i] = chem_comp_cats_d.comp[i] + chem_comp_kos.comp[i] - chem_comp7.comp[i];
	}

	get_data chem_comp_acm    = get_gb_comp(&density, &gb_tmp, PP_db, chem_comp_acm, EM_database, bulk_rock, P, T, "acm", "equilibrium");
	double gb_acm    = gb_tmp;	
	double rho_acm   = density;
	double gb5       = gb_cats_d + gb_acm - gb7 - 3.45;
	SS_ref_cpx_db.density[4] = rho_cats_d + rho_acm - SS_ref_cpx_db.density[6];

	double chem_comp5[nEl];	
	for (i = 0; i < nEl; i++){
		chem_comp5[i] = chem_comp_cats_d.comp[i] + chem_comp_acm.comp[i] - chem_comp7.comp[i];
	}

	get_data chem_comp_per     = get_gb_comp(&density, &gb_tmp, PP_db, chem_comp_per, EM_database, bulk_rock, P, T, "per", "equilibrium");
	double gb_per     = gb_tmp;	
	double rho_per    = density;
	
	get_data chem_comp_ru     = get_gb_comp(&density, &gb_tmp, PP_db, chem_comp_ru, EM_database, bulk_rock, P, T, "ru", "equilibrium");
	double gb_ru     = gb_tmp;	
	double rho_ru    = density;

	get_data chem_comp_cor    = get_gb_comp(&density, &gb_tmp, PP_db, chem_comp_cor, EM_database, bulk_rock, P, T, "cor", "equilibrium");
	double gb_cor    = gb_tmp;	
	double rho_cor   = density;
	
	double gb6       = gb_cats_d + (gb_per+gb_ru-gb_cor)/2. -16.2 -0.0012*T -0.005*P;
	SS_ref_cpx_db.density[5] = rho_cats_d + (rho_per+rho_ru-rho_cor)/2.;
	double chem_comp6[nEl];	
	for (i = 0; i < nEl; i++){
		chem_comp6[i] = chem_comp_cats_d.comp[i] + (chem_comp_per.comp[i] + chem_comp_ru.comp[i] - chem_comp_cor.comp[i])/2.0;
	}

	get_data chem_comp8    = get_gb_comp(&density, &gb_tmp, PP_db, chem_comp8, EM_database, bulk_rock, P, T, "en", "equilibrium");
	double gb_en  = gb_tmp;	
	double rho_en = density;
	double gb8    = gb_en + 3.5 - 0.002*T + 0.048*P;
	SS_ref_cpx_db.density[7] = rho_en;
	
	double gb9    = (gb_en + gb_fs)/2. - 1.60 - 0.002*T + 0.0465*P;
	SS_ref_cpx_db.density[8] = (rho_en + SS_ref_cpx_db.density[1])/2.;
	
	double chem_comp9[nEl];	
	for (i = 0; i < nEl; i++){
		chem_comp9[i] = (chem_comp8.comp[i] + chem_comp2.comp[i])/2.0;
	}
	
	get_data chem_comp_abh    = get_gb_comp(&density, &gb_tmp, PP_db, chem_comp_abh, EM_database, bulk_rock, P, T, "abh", "equilibrium");
	double gb_abh    = gb_tmp;	
	double rho_abh   = density;
	
	get_data chem_comp_san    = get_gb_comp(&density, &gb_tmp, PP_db, chem_comp_san, EM_database, bulk_rock, P, T, "san", "equilibrium");
	double gb_san    = gb_tmp;		
	double rho_san   = density;
	double gb10      = gb7 - gb_abh + gb_san + 11.7 + 0.6*P;
	SS_ref_cpx_db.density[9] = SS_ref_cpx_db.density[6] - rho_abh + rho_san;
	
	double chem_comp10[nEl];	
	for (i = 0; i < nEl; i++){
		chem_comp10[i] = chem_comp7.comp[i] - chem_comp_abh.comp[i] + chem_comp_san.comp[i];
	}

	SS_ref_cpx_db.gbase[0] = gb1;
	SS_ref_cpx_db.gbase[1] = gb2;
	SS_ref_cpx_db.gbase[2] = gb3;
	SS_ref_cpx_db.gbase[3] = gb4;
	SS_ref_cpx_db.gbase[4] = gb5;
	SS_ref_cpx_db.gbase[5] = gb6;
	SS_ref_cpx_db.gbase[6] = gb7;
	SS_ref_cpx_db.gbase[7] = gb8;
	SS_ref_cpx_db.gbase[8] = gb9;
	SS_ref_cpx_db.gbase[9] = gb10;
	
	for (i = 0; i < nEl; i++){
		SS_ref_cpx_db.Comp[0][i] = chem_comp1.comp[i];
		SS_ref_cpx_db.Comp[1][i] = chem_comp2.comp[i];
		SS_ref_cpx_db.Comp[2][i] = chem_comp3.comp[i];
		SS_ref_cpx_db.Comp[3][i] = chem_comp4[i];
		SS_ref_cpx_db.Comp[4][i] = chem_comp5[i];
		SS_ref_cpx_db.Comp[5][i] = chem_comp6[i];
		SS_ref_cpx_db.Comp[6][i] = chem_comp7.comp[i];
		SS_ref_cpx_db.Comp[7][i] = chem_comp8.comp[i];
		SS_ref_cpx_db.Comp[8][i] = chem_comp9[i];
		SS_ref_cpx_db.Comp[9][i] = chem_comp10[i];
	}
	
	/* define endmembers to be deactivated when bulk-rock composition equal zero */
	for (i = 0; i < n_em; i++){
		SS_ref_cpx_db.z_em[i] = 1.0;
	}

	SS_ref_cpx_db.box_bounds_default[0][0] = 0.+eps;  SS_ref_cpx_db.box_bounds_default[0][1] = 1.-eps;
	SS_ref_cpx_db.box_bounds_default[1][0] = 0.+eps;  SS_ref_cpx_db.box_bounds_default[1][1] = 2.-eps;
	SS_ref_cpx_db.box_bounds_default[2][0] = 0.+eps;  SS_ref_cpx_db.box_bounds_default[2][1] = 1.-eps;
	SS_ref_cpx_db.box_bounds_default[3][0] = 0.+eps;  SS_ref_cpx_db.box_bounds_default[3][1] = 1.-eps;
	SS_ref_cpx_db.box_bounds_default[4][0] = -1.+eps; SS_ref_cpx_db.box_bounds_default[4][1] = 1.-eps;	
	SS_ref_cpx_db.box_bounds_default[5][0] = 0.+eps;  SS_ref_cpx_db.box_bounds_default[5][1] = 1.-eps;
	SS_ref_cpx_db.box_bounds_default[6][0] = 0.+eps;  SS_ref_cpx_db.box_bounds_default[6][1] = 1.-eps;
	SS_ref_cpx_db.box_bounds_default[7][0] = 0.+eps;  SS_ref_cpx_db.box_bounds_default[7][1] = 1.-eps;
	SS_ref_cpx_db.box_bounds_default[8][0] = 0.+eps;  SS_ref_cpx_db.box_bounds_default[8][1] = 1.-eps;	
	
	if (bulk_rock[7] == 0.){	//TiO2
		SS_ref_cpx_db.z_em[5]          = 0.0;
		SS_ref_cpx_db.box_bounds_default[7][0] = eps;  
		SS_ref_cpx_db.box_bounds_default[7][1] = eps;
	}
	if (bulk_rock[8] == 0.){ 		//O
		SS_ref_cpx_db.z_em[4]          = 0.0;
		SS_ref_cpx_db.box_bounds_default[5][0] = eps; 
		SS_ref_cpx_db.box_bounds_default[5][1] = eps;
	}
	if (bulk_rock[9] == 0.){ 	//Cr2O3
		SS_ref_cpx_db.z_em[3]          = 0.0;
		SS_ref_cpx_db.box_bounds_default[6][0] = eps;  
		SS_ref_cpx_db.box_bounds_default[6][1] = eps;
	}

	return SS_ref_cpx_db;
}


/**
  retrieve reference thermodynamic data for cordierite 
*/
SS_ref G_SS_cd_function(SS_ref SS_ref_cd_db, int EM_database, double *bulk_rock, double P, 
								  double T, double eps){
		
	char   *EM_tmp[] 		= {"crd","fcrd","hcrd"};	
	for (int i = 0; i < SS_ref_cd_db.n_em; i++){ 
		strcpy(SS_ref_cd_db.EM_list[i],EM_tmp[i]);			
	}
	
    SS_ref_cd_db.W[0] = 6.;
    SS_ref_cd_db.W[1] = 0.;
    SS_ref_cd_db.W[2] = 0.;  
    
    PP_ref PP_db;

	double gb_tmp;
    double density;
    int i,j;
	
	int n_em = SS_ref_cd_db.n_em;
	
	get_data chem_comp1       = get_gb_comp(&density, &gb_tmp, PP_db, chem_comp1, EM_database, bulk_rock, P, T, "crd", "equilibrium");
	double gb1       = gb_tmp;	
	SS_ref_cd_db.density[0] = density;
	
	get_data chem_comp2       = get_gb_comp(&density, &gb_tmp, PP_db,chem_comp2,  EM_database, bulk_rock, P, T, "fcrd", "equilibrium");
	double gb2       = gb_tmp;		
	SS_ref_cd_db.density[1] = density;
	
	get_data chem_comp3       = get_gb_comp(&density, &gb_tmp, PP_db, chem_comp3, EM_database, bulk_rock, P, T, "hcrd", "equilibrium");
	double gb3       = gb_tmp;	
	SS_ref_cd_db.density[2] = density;
	
	SS_ref_cd_db.gbase[0] = gb1;
	SS_ref_cd_db.gbase[1] = gb2;
	SS_ref_cd_db.gbase[2] = gb3;
	
	for (i = 0; i < nEl; i++){
		SS_ref_cd_db.Comp[0][i] = chem_comp1.comp[i];
		SS_ref_cd_db.Comp[1][i] = chem_comp2.comp[i];
		SS_ref_cd_db.Comp[2][i] = chem_comp3.comp[i];
	}
	
	/* define endmembers to be deactivated when bulk-rock composition equal zero */
	for (i = 0; i < n_em; i++){
		SS_ref_cd_db.z_em[i] = 1.0;
	}

	/* define box bounds */
	SS_ref_cd_db.box_bounds_default[0][0] = 0.+eps;	SS_ref_cd_db.box_bounds_default[0][1] = 1.-eps;
	SS_ref_cd_db.box_bounds_default[1][0] = 0.+eps;	SS_ref_cd_db.box_bounds_default[1][1] = 1.-eps;

	if (bulk_rock[10] == 0.){
		SS_ref_cd_db.z_em[2]          = 0.0;
		SS_ref_cd_db.box_bounds_default[1][0] = eps;	
		SS_ref_cd_db.box_bounds_default[1][1] = eps;
	}
	
	return SS_ref_cd_db;	
}

/**
  retrieve reference thermodynamic data for epidote 
*/
SS_ref G_SS_ep_function(SS_ref SS_ref_ep_db, int EM_database, double *bulk_rock, double P, 
								  double T, double eps){
									  	
	char   *EM_tmp[] 		= {"cz","ep","fep"};	
	for (int i = 0; i < SS_ref_ep_db.n_em; i++){ 
		strcpy(SS_ref_ep_db.EM_list[i],EM_tmp[i]);			
	}
		
    SS_ref_ep_db.W[0] = 1.;
    SS_ref_ep_db.W[1] = 3.;
    SS_ref_ep_db.W[2] = 1.; 
		
    PP_ref PP_db;

	double gb_tmp;
    double density;
    int i,j;
	
	int n_em = SS_ref_ep_db.n_em;
	
	get_data chem_comp1       = get_gb_comp(&density, &gb_tmp, PP_db, chem_comp1, EM_database, bulk_rock, P, T, "cz", "equilibrium");
	double gb1       = gb_tmp;	
	SS_ref_ep_db.density[0] = density;
	
	get_data chem_comp2       = get_gb_comp(&density, &gb_tmp, PP_db, chem_comp2, EM_database, bulk_rock, P, T, "ep", "equilibrium"); //ordered?
	double gb2       = gb_tmp;		
	SS_ref_ep_db.density[1] = density;
	
	get_data chem_comp3       = get_gb_comp(&density, &gb_tmp, PP_db, chem_comp3, EM_database, bulk_rock, P, T, "fep", "equilibrium");
	double gb3       = gb_tmp;	
	SS_ref_ep_db.density[2] = density;
	
	SS_ref_ep_db.gbase[0] = gb1;
	SS_ref_ep_db.gbase[1] = gb2;
	SS_ref_ep_db.gbase[2] = gb3;
	
	for (i = 0; i < nEl; i++){
		SS_ref_ep_db.Comp[0][i] = chem_comp1.comp[i];
		SS_ref_ep_db.Comp[1][i] = chem_comp2.comp[i];
		SS_ref_ep_db.Comp[2][i] = chem_comp3.comp[i];
	}
	
	/* define endmembers to be deactivated when bulk-rock composition equal zero */
	for (i = 0; i < n_em; i++){
		SS_ref_ep_db.z_em[i] = 1.0;
	}

	/* define box bounds */
	SS_ref_ep_db.box_bounds_default[0][0] = 0.+eps;	SS_ref_ep_db.box_bounds_default[0][1] = 1.-eps;
	SS_ref_ep_db.box_bounds_default[1][0] = 0.+eps;	SS_ref_ep_db.box_bounds_default[1][1] = 0.5-eps;


	return SS_ref_ep_db;	
}

/**
  retrieve reference thermodynamic data for fluid 
*/
SS_ref G_SS_fl_function(SS_ref SS_ref_fl_db, int EM_database, double *bulk_rock, double P, 
								  double T, double eps){
									  	
	char   *EM_tmp[] 		= {"qfL","slfL","wofL","fofL","fafL","jdfL","hmfL","ekfL","tifL","kjfL","h2o"};	
	for (int i = 0; i < SS_ref_fl_db.n_em; i++){ 
		strcpy(SS_ref_fl_db.EM_list[i],EM_tmp[i]);			
	}
							  
    SS_ref_fl_db.W[0] = 0.0;
    SS_ref_fl_db.W[1] = 0.0;
    SS_ref_fl_db.W[2] = 0.0;
    SS_ref_fl_db.W[3] = 0.0;
    SS_ref_fl_db.W[4] = 0.0;
    SS_ref_fl_db.W[5] = 0.0;
    SS_ref_fl_db.W[6] = 0.0;
    SS_ref_fl_db.W[7] = 0.0;
    SS_ref_fl_db.W[8] = 0.0;
    SS_ref_fl_db.W[9] = 59.0 - 0.82*P;
    SS_ref_fl_db.W[10] = 0.0;
    SS_ref_fl_db.W[11] = 0.0;
    SS_ref_fl_db.W[12] = 0.0;
    SS_ref_fl_db.W[13] = 0.0;
    SS_ref_fl_db.W[14] = 0.0;
    SS_ref_fl_db.W[15] = 0.0;
    SS_ref_fl_db.W[16] = 0.0;
    SS_ref_fl_db.W[17] = 0.0;
    SS_ref_fl_db.W[18] = 57.6 - 0.80*P;
    SS_ref_fl_db.W[19] = 0.0;
    SS_ref_fl_db.W[20] = 0.0;
    SS_ref_fl_db.W[21] = 0.0;
    SS_ref_fl_db.W[22] = 0.0;
    SS_ref_fl_db.W[23] = 0.0;
    SS_ref_fl_db.W[24] = 0.0;
    SS_ref_fl_db.W[25] = 0.0;
    SS_ref_fl_db.W[26] = 72.2 - 0.67*P;
    SS_ref_fl_db.W[27] = 0.0;
    SS_ref_fl_db.W[28] = 0.0;
    SS_ref_fl_db.W[29] = 0.0;
    SS_ref_fl_db.W[30] = 0.0;
    SS_ref_fl_db.W[31] = 0.0;
    SS_ref_fl_db.W[32] = 0.0;
    SS_ref_fl_db.W[33] = 71.7 - 1.10*P;
    SS_ref_fl_db.W[34] = 0.0;
    SS_ref_fl_db.W[35] = 0.0;
    SS_ref_fl_db.W[36] = 0.0;
    SS_ref_fl_db.W[37] = 0.0;
    SS_ref_fl_db.W[38] = 0.0;
    SS_ref_fl_db.W[39] = 71.7 - 1.10*P;
    SS_ref_fl_db.W[40] = 0.0;
    SS_ref_fl_db.W[41] = 0.0;
    SS_ref_fl_db.W[42] = 0.0;
    SS_ref_fl_db.W[43] = 0.0;
    SS_ref_fl_db.W[44] = 57.0 - 0.79*P;
    SS_ref_fl_db.W[45] = 0.0;
    SS_ref_fl_db.W[46] = 0.0;
    SS_ref_fl_db.W[47] = 0.0;
    SS_ref_fl_db.W[48] = 73.0 - 0.66*P;
    SS_ref_fl_db.W[49] = 0.0;
    SS_ref_fl_db.W[50] = 0.0;
    SS_ref_fl_db.W[51] = 73.0 - 0.66*P;
    SS_ref_fl_db.W[52] = 0.0;
    SS_ref_fl_db.W[53] = 75.0 - 0.67*P;
    SS_ref_fl_db.W[54] = 44.9 - 1.19*P;

		
    PP_ref PP_db;

	double gb_tmp;
    double density;
    int i,j;
	
	int n_em = SS_ref_fl_db.n_em;

	get_data chem_comp_qL     = get_gb_comp(&density, &gb_tmp, PP_db, chem_comp_qL, EM_database, bulk_rock, P, T, "qL", "equilibrium");
	double gb_qL     = gb_tmp;
	double rho_qL    = density;
	double gb1       = 4.*gb_qL + 2.10 - 0.051*P;	
	SS_ref_fl_db.density[0] = rho_qL;
	
	double chem_comp1[nEl];	
	for (i = 0; i < nEl; i++){
		chem_comp1[i] = 4.*chem_comp_qL.comp[i];
	}	
	
	get_data chem_comp2       = get_gb_comp(&density, &gb_tmp, PP_db, chem_comp2, EM_database, bulk_rock, P, T, "silL", "equilibrium");
	double gb_silL   = gb_tmp;
	double gb2       = gb_silL + 6.72 - 0.313*P;
	SS_ref_fl_db.density[1] = density;
	
	get_data chem_comp3       = get_gb_comp(&density, &gb_tmp, PP_db, chem_comp3, EM_database, bulk_rock, P, T, "woL", "equilibrium");
	double gb_woL    = gb_tmp;
	double gb3       = gb_woL + 0.22 - 0.120*P;
	SS_ref_fl_db.density[2] = density;
	
	get_data chem_comp_foL    = get_gb_comp(&density, &gb_tmp, PP_db, chem_comp_foL, EM_database, bulk_rock, P, T, "foL", "equilibrium");
	double gb_foL    = gb_tmp;	
	double gb4       = 2.*gb_foL + 8.59 - 0.136*P;
	SS_ref_fl_db.density[3] = density;
	double chem_comp4[nEl];	
	for (i = 0; i < nEl; i++){
		chem_comp4[i] = 4.*chem_comp_foL.comp[i];
	}	

	get_data chem_comp_faL    = get_gb_comp(&density, &gb_tmp, PP_db, chem_comp_faL, EM_database, bulk_rock, P, T, "faL", "equilibrium");
	double gb_faL    = gb_tmp;	
	double gb5       = 2.*gb_faL + 13.56 - 0.052*P;
	SS_ref_fl_db.density[4] = density;
	double chem_comp5[nEl];	
	for (i = 0; i < nEl; i++){
		chem_comp5[i] = 4.*chem_comp_faL.comp[i];
	}	
	
	get_data chem_comp_abL    = get_gb_comp(&density, &gb_tmp, PP_db, chem_comp_abL, EM_database, bulk_rock, P, T, "abL", "equilibrium");
	double gb_abL    = gb_tmp;
	double rho_abL   = density;	
	double gb6       = gb_abL - gb_qL + 12.32 - 0.099*P;
	SS_ref_fl_db.density[5] = rho_abL - rho_qL;
	double chem_comp6[nEl];	
	for (i = 0; i < nEl; i++){
		chem_comp6[i] = chem_comp_abL.comp[i] - chem_comp_qL.comp[i];
	}	
	
	get_data chem_comp_hemL   = get_gb_comp(&density, &gb_tmp, PP_db, chem_comp_hemL, EM_database, bulk_rock, P, T, "hemL", "equilibrium");
	double gb_hemL   = gb_tmp;	
	double gb7       = gb_hemL/2. + 4.05 - 0.077*P;
	SS_ref_fl_db.density[6] = density;
	double chem_comp7[nEl];	
	for (i = 0; i < nEl; i++){
		chem_comp7[i] = chem_comp_hemL.comp[i]/2.;
	}	
	
	get_data chem_comp_eskL   = get_gb_comp(&density, &gb_tmp, PP_db, chem_comp_eskL, EM_database, bulk_rock, P, T, "eskL", "equilibrium");
	double gb_eskL   = gb_tmp;	
	double gb8       = gb_eskL/2. + 24.75 + 0.245*P;
	SS_ref_fl_db.density[7] = density;
	double chem_comp8[nEl];	
	for (i = 0; i < nEl; i++){
		chem_comp8[i] = chem_comp_eskL.comp[i]/2.;
	}	
	
	get_data chem_comp9       = get_gb_comp(&density, &gb_tmp, PP_db, chem_comp9, EM_database, bulk_rock, P, T, "ruL", "equilibrium");
	double gb_tiL    = gb_tmp;	
	double gb9       = gb_tiL + 5.60 - 0.489*P;
	SS_ref_fl_db.density[8] = density;
	
	get_data chem_comp_kspL   = get_gb_comp(&density, &gb_tmp, PP_db, chem_comp_kspL, EM_database, bulk_rock, P, T, "kspL", "equilibrium");
	double gb_kspL   = gb_tmp;	
	double rho_kspL  = density;
	double gb10      = gb_kspL - gb_qL + 12.88 - 0.227*P;
	SS_ref_fl_db.density[9] = rho_kspL - rho_qL;
	double chem_comp10[nEl];	
	for (i = 0; i < nEl; i++){
		chem_comp10[i] = chem_comp_kspL.comp[i] - chem_comp_qL.comp[i];
	}	
	
	get_data chem_comp11      = get_gb_comp(&density, &gb_tmp, PP_db, chem_comp11, EM_database, bulk_rock, P, T, "H2O", "equilibrium");
	double gb11      = gb_tmp;
	SS_ref_fl_db.density[10] = density;

	SS_ref_fl_db.gbase[0] = gb1;
	SS_ref_fl_db.gbase[1] = gb2;
	SS_ref_fl_db.gbase[2] = gb3;
	SS_ref_fl_db.gbase[3] = gb4;
	SS_ref_fl_db.gbase[4] = gb5;
	SS_ref_fl_db.gbase[5] = gb6;
	SS_ref_fl_db.gbase[6] = gb7;
	SS_ref_fl_db.gbase[7] = gb8;
	SS_ref_fl_db.gbase[8] = gb9;
	SS_ref_fl_db.gbase[9] = gb10;
	SS_ref_fl_db.gbase[10] = gb11;
		
	for (i = 0; i < nEl; i++){
		SS_ref_fl_db.Comp[0][i] = chem_comp1[i];
		SS_ref_fl_db.Comp[1][i] = chem_comp2.comp[i];
		SS_ref_fl_db.Comp[2][i] = chem_comp3.comp[i];
		SS_ref_fl_db.Comp[3][i] = chem_comp4[i];
		SS_ref_fl_db.Comp[4][i] = chem_comp5[i];
		SS_ref_fl_db.Comp[5][i] = chem_comp6[i];
		SS_ref_fl_db.Comp[6][i] = chem_comp7[i];
		SS_ref_fl_db.Comp[7][i] = chem_comp8[i];
		SS_ref_fl_db.Comp[8][i] = chem_comp9.comp[i];
		SS_ref_fl_db.Comp[9][i] = chem_comp10[i];
		SS_ref_fl_db.Comp[10][i] = chem_comp11.comp[i];
	}

	/* define endmembers to be deactivated when bulk-rock composition equal zero */
	for (i = 0; i < n_em; i++){
		SS_ref_fl_db.z_em[i] = 1.0;
	}		

	SS_ref_fl_db.box_bounds_default[0][0] = 0.+eps;  SS_ref_fl_db.box_bounds_default[0][1] = 1.-eps;
	SS_ref_fl_db.box_bounds_default[1][0] = 0.+eps;  SS_ref_fl_db.box_bounds_default[1][1] = 1.-eps;
	SS_ref_fl_db.box_bounds_default[2][0] = 0.+eps;  SS_ref_fl_db.box_bounds_default[2][1] = 1.-eps;
	SS_ref_fl_db.box_bounds_default[3][0] = 0.+eps;  SS_ref_fl_db.box_bounds_default[3][1] = 1.-eps;
	SS_ref_fl_db.box_bounds_default[4][0] = 0.+eps;  SS_ref_fl_db.box_bounds_default[4][1] = 1.-eps;	
	SS_ref_fl_db.box_bounds_default[5][0] = 0.+eps;  SS_ref_fl_db.box_bounds_default[5][1] = 1.-eps;
	SS_ref_fl_db.box_bounds_default[6][0] = 0.+eps;  SS_ref_fl_db.box_bounds_default[6][1] = 1.-eps;
	SS_ref_fl_db.box_bounds_default[7][0] = 0.+eps;  SS_ref_fl_db.box_bounds_default[7][1] = 1.-eps;
	SS_ref_fl_db.box_bounds_default[8][0] = 0.+eps;  SS_ref_fl_db.box_bounds_default[8][1] = 1.-eps;	
	SS_ref_fl_db.box_bounds_default[9][0] = 0.+eps;  SS_ref_fl_db.box_bounds_default[9][1] = 1.-eps;	

	if (bulk_rock[7] == 0.){
		SS_ref_fl_db.z_em[8]          = 0.0;
		SS_ref_fl_db.box_bounds_default[7][0] = eps;  
		SS_ref_fl_db.box_bounds_default[7][1] = eps;
	}
	if (bulk_rock[8] == 0.){
		SS_ref_fl_db.z_em[6]          = 0.0;
		SS_ref_fl_db.box_bounds_default[5][0] = eps;  
		SS_ref_fl_db.box_bounds_default[5][1] = eps;
	}
	if (bulk_rock[9] == 0.){
		SS_ref_fl_db.z_em[7]          = 0.0;
		SS_ref_fl_db.box_bounds_default[6][0] = eps;  
		SS_ref_fl_db.box_bounds_default[6][1] = eps;
	}
	if (bulk_rock[10] == 0.){
		SS_ref_fl_db.z_em[10]         = 0.0;
		SS_ref_fl_db.box_bounds_default[9][0] = eps;  
		SS_ref_fl_db.box_bounds_default[9][1] = eps;	
	}

	return SS_ref_fl_db;
}

/**
  retrieve reference thermodynamic data for garnet 
*/
SS_ref G_SS_g_function(SS_ref SS_ref_g_db, int EM_database, double *bulk_rock, double P, 
								  double T, double eps){
									  	
	char   *EM_tmp[] 		= {"py","alm","gr","andr","knom","tig"};	
	for (int i = 0; i < SS_ref_g_db.n_em; i++){ 
		strcpy(SS_ref_g_db.EM_list[i],EM_tmp[i]);			
	}
		
    SS_ref_g_db.W[0] = 4.0 + 0.10*P;
    SS_ref_g_db.W[1] = 45.4 - 0.010*T + 0.04*P;
    SS_ref_g_db.W[2] = 107.0 - 0.010*T - 0.036*P;
    SS_ref_g_db.W[3] = 2.0;
    SS_ref_g_db.W[4] = 0.0;
    SS_ref_g_db.W[5] = 17.0 - 0.010*T + 0.10*P;
    SS_ref_g_db.W[6] = 65.0 - 0.010*T + 0.039*P;
    SS_ref_g_db.W[7] = 6.0 + 0.01*P;
    SS_ref_g_db.W[8] = 0.0;
    SS_ref_g_db.W[9] = 2.0;
    SS_ref_g_db.W[10] = 1.0 - 0.010*T + 0.18*P;
    SS_ref_g_db.W[11] = 0.0;
    SS_ref_g_db.W[12] = 63.0 - 0.010*T + 0.10*P;
    SS_ref_g_db.W[13] = 0.0;
    SS_ref_g_db.W[14] = 0.0;
    
    SS_ref_g_db.v[0] = 1.00;
    SS_ref_g_db.v[1] = 1.00;
    SS_ref_g_db.v[2] = 2.50;
    SS_ref_g_db.v[3] = 2.50;
    SS_ref_g_db.v[4] = 1.00;
    SS_ref_g_db.v[5] = 1.00;
				
    PP_ref PP_db;

	double gb_tmp;
    double density;
    int i,j;
	
	int n_em = SS_ref_g_db.n_em;

	get_data chem_comp1       = get_gb_comp(&density, &gb_tmp, PP_db, chem_comp1, EM_database, bulk_rock, P, T, "py", "equilibrium");
	double gb1       = gb_tmp;
	SS_ref_g_db.density[0] = density;
	
	get_data chem_comp2       = get_gb_comp(&density, &gb_tmp, PP_db, chem_comp2, EM_database, bulk_rock, P, T, "alm", "equilibrium");
	double gb2       = gb_tmp;
	SS_ref_g_db.density[1] = density;
	
	get_data chem_comp3       = get_gb_comp(&density, &gb_tmp, PP_db, chem_comp3, EM_database, bulk_rock, P, T, "gr", "equilibrium");
	double gb3       = gb_tmp;
	SS_ref_g_db.density[2] = density;
	
	get_data chem_comp4       = get_gb_comp(&density, &gb_tmp, PP_db, chem_comp4, EM_database, bulk_rock, P, T, "andr", "equilibrium");
	double gb4       = gb_tmp;
	SS_ref_g_db.density[3] = density;
	
	get_data chem_comp5       = get_gb_comp(&density, &gb_tmp, PP_db, chem_comp5, EM_database, bulk_rock, P, T, "knor", "equilibrium");
	double gb5       = gb_tmp;
	gb5             += 18.2;
	SS_ref_g_db.density[4] = density;
	
	get_data chem_comp_per     = get_gb_comp(&density, &gb_tmp, PP_db, chem_comp_per, EM_database, bulk_rock, P, T, "per", "equilibrium");
	double gb_per     = gb_tmp;	
	double rho_per    = density;
	
	get_data chem_comp_ru     = get_gb_comp(&density, &gb_tmp, PP_db, chem_comp_ru, EM_database, bulk_rock, P, T, "ru", "equilibrium");
	double gb_ru     = gb_tmp;	
	double rho_ru    = density;
	
	get_data chem_comp_cor    = get_gb_comp(&density, &gb_tmp, PP_db, chem_comp_cor, EM_database, bulk_rock, P, T, "cor", "equilibrium");
	double gb_cor    = gb_tmp;	
	double rho_cor   = density;
	double gb6       = gb1 + (gb_per+gb_ru-gb_cor)/2. + 46.70 - 0.0173*T;
	SS_ref_g_db.density[5] = SS_ref_g_db.density[0] + (rho_per+rho_ru-rho_cor)/2.;

	double chem_comp6[nEl];	
	for (i = 0; i < nEl; i++){
		chem_comp6[i] = chem_comp1.comp[i] + (chem_comp_per.comp[i]+chem_comp_ru.comp[i]-chem_comp_cor.comp[i])/2.;
	}	
	
	SS_ref_g_db.gbase[0] = gb1;
	SS_ref_g_db.gbase[1] = gb2;
	SS_ref_g_db.gbase[2] = gb3;
	SS_ref_g_db.gbase[3] = gb4;
	SS_ref_g_db.gbase[4] = gb5;
	SS_ref_g_db.gbase[5] = gb6;
	
	for (i = 0; i < nEl; i++){
		SS_ref_g_db.Comp[0][i] = chem_comp1.comp[i];
		SS_ref_g_db.Comp[1][i] = chem_comp2.comp[i];
		SS_ref_g_db.Comp[2][i] = chem_comp3.comp[i];
		SS_ref_g_db.Comp[3][i] = chem_comp4.comp[i];
		SS_ref_g_db.Comp[4][i] = chem_comp5.comp[i];
		SS_ref_g_db.Comp[5][i] = chem_comp6[i];
	}
	
	/* define endmembers to be deactivated when bulk-rock composition equal zero */
	for (i = 0; i < n_em; i++){
		SS_ref_g_db.z_em[i] = 1.0;
	}
	
	SS_ref_g_db.box_bounds_default[0][0] = 0.+eps;	SS_ref_g_db.box_bounds_default[0][1] = 1.-eps;
	SS_ref_g_db.box_bounds_default[1][0] = 0.+eps;	SS_ref_g_db.box_bounds_default[1][1] = 1.-eps;
	SS_ref_g_db.box_bounds_default[2][0] = 0.+eps;  SS_ref_g_db.box_bounds_default[2][1] = 1.-eps;
	SS_ref_g_db.box_bounds_default[3][0] = 0.+eps;  SS_ref_g_db.box_bounds_default[3][1] = 1.-eps;
	SS_ref_g_db.box_bounds_default[4][0] = 0.+eps;	SS_ref_g_db.box_bounds_default[4][1] = 1.-eps;	

	if (bulk_rock[8] == 0.){
		SS_ref_g_db.z_em[3]          = 0.0;
		SS_ref_g_db.box_bounds_default[2][0] = eps;  
		SS_ref_g_db.box_bounds_default[2][1] = eps;
	}
	if (bulk_rock[9] == 0.){
		SS_ref_g_db.z_em[4]          = 0.0;
		SS_ref_g_db.box_bounds_default[3][0] = eps;  
		SS_ref_g_db.box_bounds_default[3][1] = eps;
	}
	if (bulk_rock[7] == 0.){
		SS_ref_g_db.z_em[5]          = 0.0;
		SS_ref_g_db.box_bounds_default[4][0] = eps;	
		SS_ref_g_db.box_bounds_default[4][1] = eps;	
	}

	return SS_ref_g_db;
}

/**
  retrieve reference thermodynamic data for horblende 
*/
SS_ref G_SS_hb_function(SS_ref SS_ref_hb_db, int EM_database, double *bulk_rock, double P, 
								  double T, double eps){
									  	
	char   *EM_tmp[] 		= {"tr","tsm","prgm","glm","cumm","grnm","a","b","mrb","kprg","tts"};	
	for (int i = 0; i < SS_ref_hb_db.n_em; i++){ 
		strcpy(SS_ref_hb_db.EM_list[i],EM_tmp[i]);			
	}
	
    SS_ref_hb_db.W[0] = 20.0;
    SS_ref_hb_db.W[1] = 25.0;
    SS_ref_hb_db.W[2] = 65.0;
    SS_ref_hb_db.W[3] = 45.0;
    SS_ref_hb_db.W[4] = 75.0;
    SS_ref_hb_db.W[5] = 57.0;
    SS_ref_hb_db.W[6] = 63.0;
    SS_ref_hb_db.W[7] = 52.0;
    SS_ref_hb_db.W[8] = 30.0;
    SS_ref_hb_db.W[9] = 85.0;
    SS_ref_hb_db.W[10] = -40.0;
    SS_ref_hb_db.W[11] = 25.0;
    SS_ref_hb_db.W[12] = 70.0;
    SS_ref_hb_db.W[13] = 80.0;
    SS_ref_hb_db.W[14] = 70.0;
    SS_ref_hb_db.W[15] = 72.5;
    SS_ref_hb_db.W[16] = 20.0;
    SS_ref_hb_db.W[17] = -40.0;
    SS_ref_hb_db.W[18] = 35.0;
    SS_ref_hb_db.W[19] = 50.0;
    SS_ref_hb_db.W[20] = 90.0;
    SS_ref_hb_db.W[21] = 106.7;
    SS_ref_hb_db.W[22] = 94.8;
    SS_ref_hb_db.W[23] = 94.8;
    SS_ref_hb_db.W[24] = 40.0;
    SS_ref_hb_db.W[25] = 8.0;
    SS_ref_hb_db.W[26] = 15.0;
    SS_ref_hb_db.W[27] = 100.0;
    SS_ref_hb_db.W[28] = 113.5;
    SS_ref_hb_db.W[29] = 100.0;
    SS_ref_hb_db.W[30] = 111.2;
    SS_ref_hb_db.W[31] = 0.0;
    SS_ref_hb_db.W[32] = 54.0;
    SS_ref_hb_db.W[33] = 75.0;
    SS_ref_hb_db.W[34] = 33.0;
    SS_ref_hb_db.W[35] = 18.0;
    SS_ref_hb_db.W[36] = 23.0;
    SS_ref_hb_db.W[37] = 80.0;
    SS_ref_hb_db.W[38] = 87.0;
    SS_ref_hb_db.W[39] = 100.0;
    SS_ref_hb_db.W[40] = 12.0;
    SS_ref_hb_db.W[41] = 8.0;
    SS_ref_hb_db.W[42] = 91.0;
    SS_ref_hb_db.W[43] = 96.0;
    SS_ref_hb_db.W[44] = 65.0;
    SS_ref_hb_db.W[45] = 20.0;
    SS_ref_hb_db.W[46] = 80.0;
    SS_ref_hb_db.W[47] = 94.0;
    SS_ref_hb_db.W[48] = 95.0;
    SS_ref_hb_db.W[49] = 90.0;
    SS_ref_hb_db.W[50] = 94.0;
    SS_ref_hb_db.W[51] = 95.0;
    SS_ref_hb_db.W[52] = 50.0;
    SS_ref_hb_db.W[53] = 50.0;
    SS_ref_hb_db.W[54] = 35.0;
    
    SS_ref_hb_db.v[0] = 1.00;
    SS_ref_hb_db.v[1] = 1.50;
    SS_ref_hb_db.v[2] = 1.70;
    SS_ref_hb_db.v[3] = 0.80;
    SS_ref_hb_db.v[4] = 1.00;
    SS_ref_hb_db.v[5] = 1.00;
    SS_ref_hb_db.v[6] = 1.00;
    SS_ref_hb_db.v[7] = 1.00;
    SS_ref_hb_db.v[8] = 0.80;
    SS_ref_hb_db.v[9] = 1.70;
    SS_ref_hb_db.v[10] = 1.50;

    PP_ref PP_db;

	double gb_tmp;
    double density;
    int i,j;
	
	int n_em = SS_ref_hb_db.n_em;
	
	get_data chem_comp1       = get_gb_comp(&density, &gb_tmp, PP_db, chem_comp1, EM_database, bulk_rock, P, T, "tr", "equilibrium");
	double gb1       = gb_tmp;
	SS_ref_hb_db.density[0] = density;
	
	get_data chem_comp2       = get_gb_comp(&density, &gb_tmp, PP_db, chem_comp2, EM_database, bulk_rock, P, T, "ts", "equilibrium");
	double gb_ts     = gb_tmp;
	double gb2       = gb_ts + 10.0;
	SS_ref_hb_db.density[1] = density;
	
	get_data chem_comp3       = get_gb_comp(&density, &gb_tmp, PP_db, chem_comp3, EM_database, bulk_rock, P, T, "parg", "equilibrium");
	double gb_parg   = gb_tmp;
	double gb3       = gb_parg - 10.0;
	SS_ref_hb_db.density[2] = density;

	get_data chem_comp4       = get_gb_comp(&density, &gb_tmp, PP_db, chem_comp4, EM_database, bulk_rock, P, T, "gl", "equilibrium");
	double gb_gl   = gb_tmp;
	double gb4       = gb_gl - 3.0;	
	SS_ref_hb_db.density[3] = density;

	get_data chem_comp5       = get_gb_comp(&density, &gb_tmp, PP_db, chem_comp5, EM_database, bulk_rock, P, T, "cumm", "equilibrium");
	double gb5       = gb_tmp;
	SS_ref_hb_db.density[4] = density;

	get_data chem_comp6       = get_gb_comp(&density, &gb_tmp, PP_db, chem_comp6, EM_database, bulk_rock, P, T, "grun", "equilibrium");
	double gb_grun   = gb_tmp;
	double gb6       = gb_grun - 3.0;	
	SS_ref_hb_db.density[5] = density;

	double gb7       = (3.0*gb5+4.0*gb_grun)/7.0 - 11.2;
	SS_ref_hb_db.density[6] = (3.0*SS_ref_hb_db.density[4]+4.0*SS_ref_hb_db.density[5])/7.0;

	
	double chem_comp7[nEl];	
	for (i = 0; i < nEl; i++){
		chem_comp7[i] = (3.0*chem_comp5.comp[i] + 4.0*chem_comp6.comp[i])/7.0;
	}	
	
	double gb8       = (2.0*gb5+5.0*gb_grun)/7.0 - 13.8;
	SS_ref_hb_db.density[7]   = (2.0*SS_ref_hb_db.density[4]+5.0*SS_ref_hb_db.density[5])/7.0 - 13.8;
	
	double chem_comp8[nEl];	
	for (i = 0; i < nEl; i++){
		chem_comp8[i] = (2.0*chem_comp5.comp[i] + 5.0*chem_comp6.comp[i])/7.0;
	}	
	
	get_data chem_comp_gr       = get_gb_comp(&density, &gb_tmp, PP_db, chem_comp_gr, EM_database, bulk_rock, P, T, "gr", "equilibrium");
	double gb_gr       = gb_tmp;
	double rho_gr      = density;
	get_data chem_comp_andr     = get_gb_comp(&density, &gb_tmp, PP_db, chem_comp_andr, EM_database, bulk_rock, P, T, "andr", "equilibrium");
	double gb_andr     = gb_tmp;
	double rho_andr    = density;
	
	double gb9         = gb_gl - gb_gr + gb_andr;
	SS_ref_hb_db.density[8] = SS_ref_hb_db.density[3] - rho_gr + rho_andr;
	double chem_comp9[nEl];	
	for (i = 0; i < nEl; i++){
		chem_comp9[i] = chem_comp4.comp[i] - chem_comp_gr.comp[i] + chem_comp_andr.comp[i];
	}	
	
	get_data chem_comp_mu       = get_gb_comp(&density, &gb_tmp, PP_db, chem_comp_mu, EM_database, bulk_rock, P, T, "mu", "equilibrium");
	double gb_mu       = gb_tmp;
	double rho_mu      = density;
	
	get_data chem_comp_pa     = get_gb_comp(&density, &gb_tmp, PP_db, chem_comp_pa, EM_database, bulk_rock, P, T, "pa", "equilibrium");
	double gb_pa     = gb_tmp;
	double rho_pa    = density;
	
	double gb10      = gb_mu - gb_pa + gb_parg - 7.06 + 0.020*T;
	SS_ref_hb_db.density[9]  = rho_mu - rho_pa + SS_ref_hb_db.density[2];
	double chem_comp10[nEl];	
	for (i = 0; i < nEl; i++){
		chem_comp10[i] = chem_comp_mu.comp[i] - chem_comp_pa.comp[i] + chem_comp3.comp[i];
	}	
	
	get_data chem_comp_dsp       = get_gb_comp(&density, &gb_tmp, PP_db, chem_comp_dsp, EM_database, bulk_rock, P, T, "dsp", "equilibrium");
	double gb_dsp       = gb_tmp;
	double rho_dsp      = density;
	
	get_data chem_comp_ru       = get_gb_comp(&density, &gb_tmp, PP_db, chem_comp_ru, EM_database, bulk_rock, P, T, "ru", "equilibrium");
	double gb_ru       = gb_tmp;
	double rho_ru      = density;
	double gb11        = -2.*gb_dsp + 2.0*gb_ru + gb_ts + 95.00;
	SS_ref_hb_db.density[10]  = -2.0*rho_dsp + 2.0*rho_ru + SS_ref_hb_db.density[1];
	double chem_comp11[nEl];	
	for (i = 0; i < nEl; i++){
		chem_comp11[i] = -2.0*chem_comp_dsp.comp[i] + 2.0*chem_comp_ru.comp[i] + chem_comp2.comp[i];
	}
	
	SS_ref_hb_db.gbase[0] = gb1;
	SS_ref_hb_db.gbase[1] = gb2;
	SS_ref_hb_db.gbase[2] = gb3;
	SS_ref_hb_db.gbase[3] = gb4;
	SS_ref_hb_db.gbase[4] = gb5;
	SS_ref_hb_db.gbase[5] = gb6;
	SS_ref_hb_db.gbase[6] = gb7;
	SS_ref_hb_db.gbase[7] = gb8;
	SS_ref_hb_db.gbase[8] = gb9;
	SS_ref_hb_db.gbase[9] = gb10;
	SS_ref_hb_db.gbase[10] = gb11;
		
	for (i = 0; i < nEl; i++){
		SS_ref_hb_db.Comp[0][i] = chem_comp1.comp[i];
		SS_ref_hb_db.Comp[1][i] = chem_comp2.comp[i];
		SS_ref_hb_db.Comp[2][i] = chem_comp3.comp[i];
		SS_ref_hb_db.Comp[3][i] = chem_comp4.comp[i];
		SS_ref_hb_db.Comp[4][i] = chem_comp5.comp[i];
		SS_ref_hb_db.Comp[5][i] = chem_comp6.comp[i];
		SS_ref_hb_db.Comp[6][i] = chem_comp7[i];
		SS_ref_hb_db.Comp[7][i] = chem_comp8[i];
		SS_ref_hb_db.Comp[8][i] = chem_comp9[i];
		SS_ref_hb_db.Comp[9][i] = chem_comp10[i];
		SS_ref_hb_db.Comp[10][i] = chem_comp11[i];
	}

	/* define endmembers to be deactivated when bulk-rock composition equal zero */
	for (i = 0; i < n_em; i++){
		SS_ref_hb_db.z_em[i] = 1.0;
	}		

	/* define box bounds according to bulk-rock */
	SS_ref_hb_db.box_bounds_default[0][0] = 0.+eps;   SS_ref_hb_db.box_bounds_default[0][1] = 1.-eps;
	SS_ref_hb_db.box_bounds_default[1][0] = 0.+eps;   SS_ref_hb_db.box_bounds_default[1][1] = 1.-eps;
	SS_ref_hb_db.box_bounds_default[2][0] = 0.+eps;   SS_ref_hb_db.box_bounds_default[2][1] = 1.-eps;
	SS_ref_hb_db.box_bounds_default[3][0] = 0.+eps;   SS_ref_hb_db.box_bounds_default[3][1] = 1.-eps;
	SS_ref_hb_db.box_bounds_default[4][0] = 0.+eps;   SS_ref_hb_db.box_bounds_default[4][1] = 1.-eps;	
	SS_ref_hb_db.box_bounds_default[5][0] = 0.+eps;   SS_ref_hb_db.box_bounds_default[5][1] = 1.-eps;
	SS_ref_hb_db.box_bounds_default[6][0] = 0.+eps;   SS_ref_hb_db.box_bounds_default[6][1] = 1.-eps;
	SS_ref_hb_db.box_bounds_default[7][0] = 0.+eps;   SS_ref_hb_db.box_bounds_default[7][1] = 1.-eps;
	SS_ref_hb_db.box_bounds_default[8][0] = -1.+eps;  SS_ref_hb_db.box_bounds_default[8][1] = 1.-eps;	
	SS_ref_hb_db.box_bounds_default[9][0] = -1.+eps;  SS_ref_hb_db.box_bounds_default[9][1] = 1.-eps;

	if (bulk_rock[7] == 0.){
		SS_ref_hb_db.z_em[10]         = 0.0;
		SS_ref_hb_db.box_bounds_default[7][0] = eps;   
		SS_ref_hb_db.box_bounds_default[7][1] = eps;
	}
	if (bulk_rock[8] == 0.){
		SS_ref_hb_db.z_em[8]          = 0.0;
		SS_ref_hb_db.box_bounds_default[6][0] = eps;   
		SS_ref_hb_db.box_bounds_default[6][1] = eps;
	}

	return SS_ref_hb_db;	
}

/**
  retrieve reference thermodynamic data for ilmenite 
*/
SS_ref G_SS_ilm_function(SS_ref SS_ref_ilm_db, int EM_database, double *bulk_rock, double P, 
								  double T, double eps){
									  	
	char   *EM_tmp[] 		= {"oilm","dilm","dhem"};	
	for (int i = 0; i < SS_ref_ilm_db.n_em; i++){ 
		strcpy(SS_ref_ilm_db.EM_list[i],EM_tmp[i]);			
	}
							  
    SS_ref_ilm_db.W[0] = 7.05;
    SS_ref_ilm_db.W[1] = 14.3;
    SS_ref_ilm_db.W[2] = 7.25;    						  
					
    PP_ref PP_db;

	double gb_tmp;
    double density;
    int i,j;
	
	int n_em = SS_ref_ilm_db.n_em;
	
	get_data chem_comp1       = get_gb_comp(&density, &gb_tmp, PP_db, chem_comp1, EM_database, bulk_rock, P, T, "ilm", "ordered");
	double gb1       = gb_tmp;	
	SS_ref_ilm_db.density[0] = density;			  
	
	get_data chem_comp2       = get_gb_comp(&density, &gb_tmp, PP_db, chem_comp2, EM_database, bulk_rock, P, T, "ilm", "disordered");
	double gb2       = gb_tmp;	
	SS_ref_ilm_db.density[1] = density;			  

	get_data chem_comp3       = get_gb_comp(&density, &gb_tmp, PP_db, chem_comp3, EM_database, bulk_rock, P, T, "hem", "equilibrium");
	double gb3       = gb_tmp;	
	SS_ref_ilm_db.density[2] = density;			  
	
	SS_ref_ilm_db.gbase[0] = gb1;
	SS_ref_ilm_db.gbase[1] = gb2;
	SS_ref_ilm_db.gbase[2] = gb3;
	
	for (i = 0; i < nEl; i++){
		SS_ref_ilm_db.Comp[0][i] = chem_comp1.comp[i];
		SS_ref_ilm_db.Comp[1][i] = chem_comp2.comp[i];
		SS_ref_ilm_db.Comp[2][i] = chem_comp3.comp[i];
	}

	/* define endmembers to be deactivated when bulk-rock composition equal zero */
	for (i = 0; i < n_em; i++){
		SS_ref_ilm_db.z_em[i] = 1.0;
	}
	
	SS_ref_ilm_db.box_bounds_default[0][0] =  0. 	+eps;	SS_ref_ilm_db.box_bounds_default[0][1] = 1.0  - eps;
	SS_ref_ilm_db.box_bounds_default[1][0] =  -0.99 +eps;	SS_ref_ilm_db.box_bounds_default[1][1] = 0.99  - eps;

	return SS_ref_ilm_db;										  
}

/**
  retrieve reference thermodynamic data for liquid (melt) 
*/
SS_ref G_SS_liq_function(SS_ref SS_ref_liq_db, int EM_database, double *bulk_rock, double P, 
								  double T, double eps){
									  	
	char   *EM_tmp[] 		= {"q4L","sl1L","wo1L","fo2L","fa2L","jdL","hmL","ekL","tiL","kjL","ctL","h2o1L"};	
	for (int i = 0; i < SS_ref_liq_db.n_em; i++){ 
		strcpy(SS_ref_liq_db.EM_list[i],EM_tmp[i]);			
	}
		  
    SS_ref_liq_db.W[0] = 9.5 - 0.10*P;
    SS_ref_liq_db.W[1] = -10.3;
    SS_ref_liq_db.W[2] = -26.5 - 3.12*P;
    SS_ref_liq_db.W[3] = -12.0 - 0.55*P;
    SS_ref_liq_db.W[4] = -15.1 - 0.13*P;
    SS_ref_liq_db.W[5] = 20.0;
    SS_ref_liq_db.W[6] = 0.0;
    SS_ref_liq_db.W[7] = 24.6;
    SS_ref_liq_db.W[8] = -17.8 - 0.05*P;
    SS_ref_liq_db.W[9] = -14.6;
    SS_ref_liq_db.W[10] = 17.8 - 0.61*P;
    SS_ref_liq_db.W[11] = -26.5 + 0.85*P;
    SS_ref_liq_db.W[12] = 2.2;
    SS_ref_liq_db.W[13] = 2.5;
    SS_ref_liq_db.W[14] = 16.8;
    SS_ref_liq_db.W[15] = -5.0;
    SS_ref_liq_db.W[16] = 0.0;
    SS_ref_liq_db.W[17] = 15.2 - 0.04*P;
    SS_ref_liq_db.W[18] = 7.0;
    SS_ref_liq_db.W[19] = 4.0;
    SS_ref_liq_db.W[20] = 23.7 - 0.94*P;
    SS_ref_liq_db.W[21] = 25.5 + 0.11*P;
    SS_ref_liq_db.W[22] = 14.0;
    SS_ref_liq_db.W[23] = -1.2;
    SS_ref_liq_db.W[24] = 0.0;
    SS_ref_liq_db.W[25] = 0.0;
    SS_ref_liq_db.W[26] = 18.0;
    SS_ref_liq_db.W[27] = -1.1;
    SS_ref_liq_db.W[28] = 9.5;
    SS_ref_liq_db.W[29] = 40.3 - 0.86*P;
    SS_ref_liq_db.W[30] = 18.0;
    SS_ref_liq_db.W[31] = 1.5;
    SS_ref_liq_db.W[32] = 0.0;
    SS_ref_liq_db.W[33] = 0.0;
    SS_ref_liq_db.W[34] = 7.5;
    SS_ref_liq_db.W[35] = 3.0;
    SS_ref_liq_db.W[36] = -5.6;
    SS_ref_liq_db.W[37] = 9.4 - 1.58*P;
    SS_ref_liq_db.W[38] = 7.5 - 0.05*P;
    SS_ref_liq_db.W[39] = -30.0;
    SS_ref_liq_db.W[40] = 0.0;
    SS_ref_liq_db.W[41] = 6.7;
    SS_ref_liq_db.W[42] = 10.0;
    SS_ref_liq_db.W[43] = -6.5;
    SS_ref_liq_db.W[44] = 9.2 - 1.58*P;
    SS_ref_liq_db.W[45] = 10.0;
    SS_ref_liq_db.W[46] = 0.0;
    SS_ref_liq_db.W[47] = 16.5 + 0.14*P;
    SS_ref_liq_db.W[48] = -5.9;
    SS_ref_liq_db.W[49] = 7.6;
    SS_ref_liq_db.W[50] = -8.3 - 0.06*P; 
    SS_ref_liq_db.W[51] = 0.0;
    SS_ref_liq_db.W[52] = 0.0;
    SS_ref_liq_db.W[53] = 10.0;
    SS_ref_liq_db.W[54] = 0.0;
    SS_ref_liq_db.W[55] = 60.0 - 0.66*P;
    SS_ref_liq_db.W[56] = 0.0;
    SS_ref_liq_db.W[57] = 0.0;
    SS_ref_liq_db.W[58] = 0.0;
    SS_ref_liq_db.W[59] = 30.0 - 0.66*P; 
    SS_ref_liq_db.W[60] = 9.0;
    SS_ref_liq_db.W[61] = 0.0;
    SS_ref_liq_db.W[62] = 30.0 - 0.60*P;
    SS_ref_liq_db.W[63] = -5.6;
    SS_ref_liq_db.W[64] = -0.1 + 0.22*P;
    SS_ref_liq_db.W[65] = 17.3 + 0.05*P;

    SS_ref_liq_db.v[0] = 100.00;
    SS_ref_liq_db.v[1] = 120.00;
    SS_ref_liq_db.v[2] = 140.00;
    SS_ref_liq_db.v[3] = 240.00;
    SS_ref_liq_db.v[4] = 100.00;
    SS_ref_liq_db.v[5] = 120.00;
    SS_ref_liq_db.v[6] = 100.00;
    SS_ref_liq_db.v[7] = 100.00;
    SS_ref_liq_db.v[8] = 100.00;
    SS_ref_liq_db.v[9] = 100.00;
    SS_ref_liq_db.v[10] = 100.00;
    SS_ref_liq_db.v[11] = 100.00;  
		  
    PP_ref PP_db;

	double gb_tmp;
    double density;
    int i,j;
	
	int n_em = SS_ref_liq_db.n_em;
	
	get_data chem_comp_qL       = get_gb_comp(&density, &gb_tmp, PP_db, chem_comp_qL, EM_database, bulk_rock, P, T, "qL", "equilibrium");
	double gb_qL       = gb_tmp;
	double rho_qL      = density;
	double gb1         = 4.*gb_qL + 0.22 - 0.059*P;
	SS_ref_liq_db.density[0] = rho_qL;
	
	double chem_comp1[nEl];	
	for (i = 0; i < nEl; i++){
		chem_comp1[i] = 4.0*chem_comp_qL.comp[i];
	}
	get_data chem_comp2       = get_gb_comp(&density, &gb_tmp, PP_db, chem_comp2, EM_database, bulk_rock, P, T, "silL", "equilibrium");
	double gb_silL   = gb_tmp;
	double rho_silL  = density;
	double gb2       = gb_silL + 6.20 - 0.318*P;
	SS_ref_liq_db.density[1] = rho_silL;

	
	get_data chem_comp3       = get_gb_comp(&density, &gb_tmp, PP_db, chem_comp3, EM_database, bulk_rock, P, T, "woL", "equilibrium");
	double gb_woL   = gb_tmp;
	double rho_woL  = density;
	double gb3      = gb_woL - 0.45 - 0.114*P;
	SS_ref_liq_db.density[2] = rho_woL;


	get_data chem_comp_foL       = get_gb_comp(&density, &gb_tmp, PP_db, chem_comp_foL, EM_database, bulk_rock, P, T, "foL", "equilibrium");
	double gb_foL    = gb_tmp;
	double rho_foL   = density;
	double gb4       = 2.0*gb_foL + 8.67 - 0.131*P;	
	SS_ref_liq_db.density[3] = rho_foL;
	double chem_comp4[nEl];	
	for (i = 0; i < nEl; i++){
		chem_comp4[i] = 2.0*chem_comp_foL.comp[i];
	}
	
	get_data chem_comp_faL       = get_gb_comp(&density, &gb_tmp, PP_db, chem_comp_faL, EM_database, bulk_rock, P, T, "faL", "equilibrium");
	double gb_faL    = gb_tmp;
	double rho_faL   = density;
	double gb5       = 2.0*gb_faL + 13.70 - 0.055*P;	
	SS_ref_liq_db.density[4] = rho_faL;
	double chem_comp5[nEl];	
	for (i = 0; i < nEl; i++){
		chem_comp5[i] = 2.0*chem_comp_faL.comp[i];
	}
	
	get_data chem_comp_abL    = get_gb_comp(&density, &gb_tmp, PP_db, chem_comp_abL, EM_database, bulk_rock, P, T, "abL", "equilibrium");
	double gb_abL    = gb_tmp;
	double rho_abL   = density;
	double gb6       = gb_abL - gb_qL + 12.19 -0.089*P;
	SS_ref_liq_db.density[5] = rho_abL - rho_qL;
	double chem_comp6[nEl];	
	for (i = 0; i < nEl; i++){
		chem_comp6[i] = chem_comp_abL.comp[i] - chem_comp_qL.comp[i];
	}

	get_data chem_comp_hemL    = get_gb_comp(&density, &gb_tmp, PP_db, chem_comp_hemL, EM_database, bulk_rock, P, T, "hemL", "equilibrium");
	double gb_hemL    = gb_tmp;
	double rho_hemL   = density;	
	double gb7        = gb_hemL/2.0 + 3.30 - 0.032*P;
	SS_ref_liq_db.density[6] = rho_hemL;
	double chem_comp7[nEl];	
	for (i = 0; i < nEl; i++){
		chem_comp7[i] = chem_comp_hemL.comp[i]/2.0;
	}
	
	get_data chem_comp_eskL    = get_gb_comp(&density, &gb_tmp, PP_db, chem_comp_eskL, EM_database, bulk_rock, P, T, "eskL", "equilibrium");
	double gb_eskL    = gb_tmp;		
	double rho_eskL   = density;
	double gb8        = gb_eskL/2.0 + 24.85 + 0.245*P;
	SS_ref_liq_db.density[7] = rho_eskL;
	double chem_comp8[nEl];	
	for (i = 0; i < nEl; i++){
		chem_comp8[i] = chem_comp_eskL.comp[i]/2.0;
	}
	
	get_data chem_comp9        = get_gb_comp(&density, &gb_tmp, PP_db, chem_comp9, EM_database, bulk_rock, P, T, "ruL", "equilibrium");
	double gb_tiL     = gb_tmp;	
	double gb9        = gb_tiL + 5.58 - 0.489*P;
	SS_ref_liq_db.density[8] = density;

	get_data chem_comp_kspL    = get_gb_comp(&density, &gb_tmp, PP_db, chem_comp_kspL, EM_database, bulk_rock, P, T, "kspL", "equilibrium");
	double gb_kspL    = gb_tmp;	
	double rho_kspL   = density;
	double gb10       = gb_kspL - gb_qL + 11.98 - 0.210*P;
	SS_ref_liq_db.density[9] = rho_kspL - rho_qL;
	double chem_comp10[nEl];	
	for (i = 0; i < nEl; i++){
		chem_comp10[i] = chem_comp_kspL.comp[i] - chem_comp_qL.comp[i];
	}
	
	double gb11       = gb_woL + gb_silL - gb_qL - 108.30 + 0.055*T + 0.053*P;
	SS_ref_liq_db.density[10] = rho_woL + rho_silL - rho_qL;
	double chem_comp11[nEl];	
	for (i = 0; i < nEl; i++){
		chem_comp11[i] = chem_comp3.comp[i] + chem_comp2.comp[i] - chem_comp_qL.comp[i];
	}
	
	get_data chem_comp12    = get_gb_comp(&density, &gb_tmp, PP_db, chem_comp12, EM_database, bulk_rock, P, T, "h2oL", "equilibrium");
	double gb_h2oL    = gb_tmp;		
	double gb12       = gb_h2oL + 3.20 - 0.0039*T + 0.00087*P;
	SS_ref_liq_db.density[11] = density;
	
	SS_ref_liq_db.gbase[0]  = gb1;
	SS_ref_liq_db.gbase[1]  = gb2;
	SS_ref_liq_db.gbase[2]  = gb3;
	SS_ref_liq_db.gbase[3]  = gb4;
	SS_ref_liq_db.gbase[4]  = gb5;
	SS_ref_liq_db.gbase[5]  = gb6;
	SS_ref_liq_db.gbase[6]  = gb7;
	SS_ref_liq_db.gbase[7]  = gb8;
	SS_ref_liq_db.gbase[8]  = gb9;
	SS_ref_liq_db.gbase[9]  = gb10;
	SS_ref_liq_db.gbase[10] = gb11;
	SS_ref_liq_db.gbase[11] = gb12;
		
	for (i = 0; i < nEl; i++){
		SS_ref_liq_db.Comp[0][i] = chem_comp1[i];
		SS_ref_liq_db.Comp[1][i] = chem_comp2.comp[i];
		SS_ref_liq_db.Comp[2][i] = chem_comp3.comp[i];
		SS_ref_liq_db.Comp[3][i] = chem_comp4[i];
		SS_ref_liq_db.Comp[4][i] = chem_comp5[i];
		SS_ref_liq_db.Comp[5][i] = chem_comp6[i];
		SS_ref_liq_db.Comp[6][i] = chem_comp7[i];
		SS_ref_liq_db.Comp[7][i] = chem_comp8[i];
		SS_ref_liq_db.Comp[8][i] = chem_comp9.comp[i];
		SS_ref_liq_db.Comp[9][i] = chem_comp10[i];
		SS_ref_liq_db.Comp[10][i] = chem_comp11[i];
		SS_ref_liq_db.Comp[11][i] = chem_comp12.comp[i];
	}

	/* define endmembers to be deactivated when bulk-rock composition equal zero */
	for (i = 0; i < n_em; i++){
		SS_ref_liq_db.z_em[i] = 1.0;
	}		

	SS_ref_liq_db.box_bounds_default[0][0]  = 0.+eps;   SS_ref_liq_db.box_bounds_default[0][1]  = 1.-eps;
	SS_ref_liq_db.box_bounds_default[1][0]  = 0.+eps;   SS_ref_liq_db.box_bounds_default[1][1]  = 1.-eps;
	SS_ref_liq_db.box_bounds_default[2][0]  = 0.+eps;   SS_ref_liq_db.box_bounds_default[2][1]  = 1.-eps;
	SS_ref_liq_db.box_bounds_default[3][0]  = 0.+eps;   SS_ref_liq_db.box_bounds_default[3][1]  = 1.-eps;
	SS_ref_liq_db.box_bounds_default[4][0]  = 0.+eps;   SS_ref_liq_db.box_bounds_default[4][1]  = 1.-eps;	
	SS_ref_liq_db.box_bounds_default[5][0]  = 0.+eps;   SS_ref_liq_db.box_bounds_default[5][1]  = 1.-eps;
	SS_ref_liq_db.box_bounds_default[6][0]  = 0.+eps;   SS_ref_liq_db.box_bounds_default[6][1]  = 1.-eps;
	SS_ref_liq_db.box_bounds_default[7][0]  = 0.+eps;   SS_ref_liq_db.box_bounds_default[7][1]  = 1.-eps;
	SS_ref_liq_db.box_bounds_default[8][0]  = 0.+eps;   SS_ref_liq_db.box_bounds_default[8][1]  = 1.-eps;	
	SS_ref_liq_db.box_bounds_default[9][0]  = 0.+eps;   SS_ref_liq_db.box_bounds_default[9][1]  = 1.-eps;	
	SS_ref_liq_db.box_bounds_default[10][0] = 0.+eps;   SS_ref_liq_db.box_bounds_default[10][1] = 1.-eps;	

	if (bulk_rock[7] == 0.){
		SS_ref_liq_db.z_em[8]          = 0.0;
		SS_ref_liq_db.box_bounds_default[7][0] = eps; 
		SS_ref_liq_db.box_bounds_default[7][1] = eps;	
	}
	if (bulk_rock[8] == 0.){					// eps seems to give bad results, uses 0.0 instead
		SS_ref_liq_db.z_em[6]          = 0.0;
		SS_ref_liq_db.box_bounds_default[5][0] = eps; 
		SS_ref_liq_db.box_bounds_default[5][1] = eps;	
	}
	if (bulk_rock[9] == 0.){
		SS_ref_liq_db.z_em[7]          = 0.0;
		SS_ref_liq_db.box_bounds_default[6][0] = eps; 
		SS_ref_liq_db.box_bounds_default[6][1] = eps;	
	}
	if (bulk_rock[10] == 0.){ 					// no h2o, cannot be 0 for this xeos
		SS_ref_liq_db.z_em[11]          = 0.0;
		SS_ref_liq_db.box_bounds_default[10][0] = eps; 
		SS_ref_liq_db.box_bounds_default[10][1] = eps;	
	}

	return SS_ref_liq_db;	
}

/**
  retrieve reference thermodynamic data for muscovite
*/
SS_ref G_SS_mu_function(SS_ref SS_ref_mu_db, int EM_database, double *bulk_rock, double P, 
								  double T, double eps){		
									  	
	char   *EM_tmp[] 		= {"mu","cel","fcel","pa","mam","fmu"};	
	for (int i = 0; i < SS_ref_mu_db.n_em; i++){ 
		strcpy(SS_ref_mu_db.EM_list[i],EM_tmp[i]);			
	}
						
    SS_ref_mu_db.W[0] = 0. + 0.20*P;
    SS_ref_mu_db.W[1] = 0. + 0.20*P;
    SS_ref_mu_db.W[2] = 10.12 + 0.0034*T + 0.353*P;
    SS_ref_mu_db.W[3] = 35.0;
    SS_ref_mu_db.W[4] = 0.;
    SS_ref_mu_db.W[5] = 0.;
    SS_ref_mu_db.W[6] = 45.0 + 0.25*P;
    SS_ref_mu_db.W[7] = 50.0;
    SS_ref_mu_db.W[8] = 0.;
    SS_ref_mu_db.W[9] = 45.0 + 0.25*P;
    SS_ref_mu_db.W[10] = 50.0;
    SS_ref_mu_db.W[11] = 0.;
    SS_ref_mu_db.W[12] = 15.0;
    SS_ref_mu_db.W[13] = 30.0;
    SS_ref_mu_db.W[14] = 35.0;
    SS_ref_mu_db.v[0] = 0.63;
    SS_ref_mu_db.v[1] = 0.63;
    SS_ref_mu_db.v[2] = 0.63;
    SS_ref_mu_db.v[3] = 0.37;
    SS_ref_mu_db.v[4] = 0.63;
    SS_ref_mu_db.v[5] = 0.63;
					
    PP_ref PP_db;

	double gb_tmp;
    double density;
    int i,j;
	
	int n_em = SS_ref_mu_db.n_em;
	
	get_data chem_comp1       = get_gb_comp(&density, &gb_tmp, PP_db, chem_comp1, EM_database, bulk_rock, P, T, "mu", "equilibrium");
	double gb1       = gb_tmp;
	SS_ref_mu_db.density[0] = density;
	
	get_data chem_comp2       = get_gb_comp(&density, &gb_tmp, PP_db, chem_comp2, EM_database, bulk_rock, P, T, "cel", "equilibrium");
	double gb2       = gb_tmp;	
	SS_ref_mu_db.density[1] = density;
	
	get_data chem_comp3       = get_gb_comp(&density, &gb_tmp, PP_db, chem_comp3, EM_database, bulk_rock, P, T, "fcel", "equilibrium");
	double gb3       = gb_tmp;
	SS_ref_mu_db.density[2] = density;
	
	get_data chem_comp4       = get_gb_comp(&density, &gb_tmp, PP_db, chem_comp4, EM_database, bulk_rock, P, T, "pa", "equilibrium");
	double gb4       = gb_tmp;
	SS_ref_mu_db.density[3] = density;
	
	get_data chem_comp5       = get_gb_comp(&density, &gb_tmp, PP_db, chem_comp5, EM_database, bulk_rock, P, T, "ma", "equilibrium");
	double gb5       = gb_tmp + 6.5;
	SS_ref_mu_db.density[4] = density;
	
	get_data chem_comp_gr       = get_gb_comp(&density, &gb_tmp, PP_db, chem_comp_gr, EM_database, bulk_rock, P, T, "gr", "equilibrium");
	double gb_gr       = gb_tmp;
	double rho_gr      = density;
	
	get_data chem_comp_andr     = get_gb_comp(&density, &gb_tmp, PP_db, chem_comp_andr, EM_database, bulk_rock, P, T, "andr", "equilibrium");
	double gb_andr     = gb_tmp;	
	double rho_andr    = density;
	double gb6         = (gb_andr - gb_gr)/2.0 + gb1 + 25.0;
	SS_ref_mu_db.density[5] = (rho_andr - rho_gr)/2.0;
	
	double chem_comp6[nEl];
	for (i = 0; i < nEl; i++){
		chem_comp6[i] = (chem_comp_andr.comp[i] - chem_comp_gr.comp[i])/2.0 + chem_comp1.comp[i];
	}

	SS_ref_mu_db.gbase[0] = gb1;
	SS_ref_mu_db.gbase[1] = gb2;
	SS_ref_mu_db.gbase[2] = gb3;
	SS_ref_mu_db.gbase[3] = gb4;
	SS_ref_mu_db.gbase[4] = gb5;
	SS_ref_mu_db.gbase[5] = gb6;
	
	for (i = 0; i < nEl; i++){
		SS_ref_mu_db.Comp[0][i] = chem_comp1.comp[i];
		SS_ref_mu_db.Comp[1][i] = chem_comp2.comp[i];
		SS_ref_mu_db.Comp[2][i] = chem_comp3.comp[i];
		SS_ref_mu_db.Comp[3][i] = chem_comp4.comp[i];
		SS_ref_mu_db.Comp[4][i] = chem_comp5.comp[i];
		SS_ref_mu_db.Comp[5][i] = chem_comp6[i];
	}
	
	SS_ref_mu_db.box_bounds_default[0][0] = 0.+eps;	SS_ref_mu_db.box_bounds_default[0][1] = 1.0 - eps;
	SS_ref_mu_db.box_bounds_default[1][0] = 0.+eps;	SS_ref_mu_db.box_bounds_default[1][1] = 1.0 - eps;
	SS_ref_mu_db.box_bounds_default[2][0] = 0.+eps; SS_ref_mu_db.box_bounds_default[2][1] = 1.0 - eps;
	SS_ref_mu_db.box_bounds_default[3][0] = 0.+eps; SS_ref_mu_db.box_bounds_default[3][1] = 1.0 - eps;
	SS_ref_mu_db.box_bounds_default[4][0] = 0.+eps;	SS_ref_mu_db.box_bounds_default[4][1] = 1.0 - eps;	

	/* define endmembers to be deactivated when bulk-rock composition equal zero */
	for (i = 0; i < n_em; i++){
		SS_ref_mu_db.z_em[i] = 1.0;
	}

	/* define box bounds according to bulk-rock */
	if (bulk_rock[8] == 0.){
		SS_ref_mu_db.z_em[5]          = 0.0;
		SS_ref_mu_db.box_bounds_default[2][0] = eps;
		SS_ref_mu_db.box_bounds_default[2][1] = eps;
	}

	return SS_ref_mu_db;
}

/**
  retrieve reference thermodynamic data for olivine
*/
SS_ref G_SS_ol_function(SS_ref SS_ref_ol_db, int EM_database, double *bulk_rock, double P, 
								  double T, double eps){		
									  	
	char   *EM_tmp[] 		= {"mont","fa","fo","cfm"};	
	for (int i = 0; i < SS_ref_ol_db.n_em; i++){ 
		strcpy(SS_ref_ol_db.EM_list[i],EM_tmp[i]);			
	}
								  
    SS_ref_ol_db.W[0] = 24.;
    SS_ref_ol_db.W[1] = 38.;
    SS_ref_ol_db.W[2] = 24.;
    SS_ref_ol_db.W[3] = 9.;
    SS_ref_ol_db.W[4] = 4.5;
    SS_ref_ol_db.W[5] = 4.5;    
			
    PP_ref PP_db;

	double gb_tmp;
    double density;
    int i,j;
	
	int n_em = SS_ref_ol_db.n_em;
	
	get_data chem_comp1       = get_gb_comp(&density, &gb_tmp, PP_db, chem_comp1, EM_database, bulk_rock, P, T, "mont", "equilibrium");
	double gb1       = gb_tmp;
	SS_ref_ol_db.density[0] = density;
	
	get_data chem_comp2       = get_gb_comp(&density, &gb_tmp, PP_db, chem_comp2, EM_database, bulk_rock, P, T, "fa", "equilibrium");
	double gb2       = gb_tmp;
	SS_ref_ol_db.density[1] = density;
	
	get_data chem_comp3       = get_gb_comp(&density, &gb_tmp, PP_db, chem_comp3, EM_database, bulk_rock, P, T, "fo", "equilibrium");
	double gb3       = gb_tmp;
	SS_ref_ol_db.density[2] = density;
	
	double gb4       = (gb2+gb3)/2.;
	SS_ref_ol_db.density[3] = (SS_ref_ol_db.density[1] + SS_ref_ol_db.density[2])/2.0;
	
	double chem_comp4[nEl];
	for (i = 0; i < nEl; i++){
		chem_comp4[i] = (chem_comp2.comp[i] + chem_comp3.comp[i])/2.;
	}

	SS_ref_ol_db.gbase[0] = gb1;
	SS_ref_ol_db.gbase[1] = gb2;
	SS_ref_ol_db.gbase[2] = gb3;
	SS_ref_ol_db.gbase[3] = gb4;

	for (i = 0; i < nEl; i++){
		SS_ref_ol_db.Comp[0][i] = chem_comp1.comp[i];
		SS_ref_ol_db.Comp[1][i] = chem_comp2.comp[i];
		SS_ref_ol_db.Comp[2][i] = chem_comp3.comp[i];
		SS_ref_ol_db.Comp[3][i] = chem_comp4[i];
	}

	/* define box bounds according to bulk-rock */
	SS_ref_ol_db.box_bounds_default[0][0] = 0. + eps;	SS_ref_ol_db.box_bounds_default[0][1] = 1. - eps;
	SS_ref_ol_db.box_bounds_default[1][0] = 0. + eps;	SS_ref_ol_db.box_bounds_default[1][1] = 1. - eps;
	SS_ref_ol_db.box_bounds_default[2][0] = 0. + eps; 	SS_ref_ol_db.box_bounds_default[2][1] = 1. - eps;

	/* define endmembers to be deactivated when bulk-rock composition equal zero */
	for (i = 0; i < n_em; i++){
		SS_ref_ol_db.z_em[i] = 1.0;
	}

	return SS_ref_ol_db;
}

/**
  retrieve reference thermodynamic data for orthopyroxene
*/
SS_ref G_SS_opx_function(SS_ref SS_ref_opx_db, int EM_database, double *bulk_rock, double P, 
								  double T, double eps){		
									  	
	char   *EM_tmp[] 		= {"en","fs","fm","odi","mgts","cren","obuf","mess","ojd"};	
	for (int i = 0; i < SS_ref_opx_db.n_em; i++){ 
		strcpy(SS_ref_opx_db.EM_list[i],EM_tmp[i]);			
	}
		
    SS_ref_opx_db.W[0] = 7.0;
    SS_ref_opx_db.W[1] = 4.0;
    SS_ref_opx_db.W[2] = 29.4;
    SS_ref_opx_db.W[3] = 12.5 - 0.04*P;
    SS_ref_opx_db.W[4] = 8.0;
    SS_ref_opx_db.W[5] = 6.0;
    SS_ref_opx_db.W[6] = 8.0;
    SS_ref_opx_db.W[7] = 35.0;
    SS_ref_opx_db.W[8] = 4.0;
    SS_ref_opx_db.W[9] = 21.5 + 0.08*P;
    SS_ref_opx_db.W[10] = 11.0 - 0.15*P;
    SS_ref_opx_db.W[11] = 10.0;
    SS_ref_opx_db.W[12] = 7.0;
    SS_ref_opx_db.W[13] = 10.0;
    SS_ref_opx_db.W[14] = 35.0;
    SS_ref_opx_db.W[15] = 18.0 + 0.08*P;
    SS_ref_opx_db.W[16] = 15.0 - 0.15*P;
    SS_ref_opx_db.W[17] = 12.0;
    SS_ref_opx_db.W[18] = 8.0;
    SS_ref_opx_db.W[19] = 12.0;
    SS_ref_opx_db.W[20] = 35.0;
    SS_ref_opx_db.W[21] = 75.5 - 0.84*P;
    SS_ref_opx_db.W[22] = 20.0;
    SS_ref_opx_db.W[23] = 40.0;
    SS_ref_opx_db.W[24] = 20.0;
    SS_ref_opx_db.W[25] = 35.0;
    SS_ref_opx_db.W[26] = 2.0;
    SS_ref_opx_db.W[27] = 10.0;
    SS_ref_opx_db.W[28] = 2.0;
    SS_ref_opx_db.W[29] = 7.0;
    SS_ref_opx_db.W[30] = 6.0;
    SS_ref_opx_db.W[31] = 2.0;
    SS_ref_opx_db.W[32] = -11.0;
    SS_ref_opx_db.W[33] = 6.0;
    SS_ref_opx_db.W[34] = 20.0;
    SS_ref_opx_db.W[35] = -11.0;
    
    SS_ref_opx_db.v[0] = 1.00;
    SS_ref_opx_db.v[1] = 1.00;
    SS_ref_opx_db.v[2] = 1.00;
    SS_ref_opx_db.v[3] = 1.20;
    SS_ref_opx_db.v[4] = 1.00;
    SS_ref_opx_db.v[5] = 1.00;
    SS_ref_opx_db.v[6] = 1.00;
    SS_ref_opx_db.v[7] = 1.00;
    SS_ref_opx_db.v[8] = 1.20;
					
    PP_ref PP_db;

	double gb_tmp;
    double density;
    int i,j;
	
	int n_em = SS_ref_opx_db.n_em;
	
	get_data chem_comp1       = get_gb_comp(&density, &gb_tmp, PP_db, chem_comp1, EM_database, bulk_rock, P, T, "en", "equilibrium");
	double gb1       = gb_tmp;	
	SS_ref_opx_db.density[0] = density;

	get_data chem_comp2       = get_gb_comp(&density, &gb_tmp, PP_db, chem_comp2, EM_database, bulk_rock, P, T, "fs", "equilibrium");
	double gb2       = gb_tmp;
	SS_ref_opx_db.density[1] = density;
	
	double gb3       = (gb1 + gb2)/2. - 6.6;
	SS_ref_opx_db.density[2] = (SS_ref_opx_db.density[0] + SS_ref_opx_db.density[1])/2.;
	
	double chem_comp3[nEl];	
	for (i = 0; i < nEl; i++){
		chem_comp3[i] = (chem_comp1.comp[i] + chem_comp2.comp[i])/2.0;
	}
	
	get_data chem_comp4       = get_gb_comp(&density, &gb_tmp, PP_db, chem_comp4, EM_database, bulk_rock, P, T, "di", "equilibrium");
	double gb4       = gb_tmp + 2.8 + 0.005*P;	
	SS_ref_opx_db.density[3] = density;
	
	get_data chem_comp5       = get_gb_comp(&density, &gb_tmp, PP_db, chem_comp5, EM_database, bulk_rock, P, T, "mgts", "equilibrium");
	double gb5       = gb_tmp;	
	SS_ref_opx_db.density[4] = density;
	
	get_data chem_comp_kos    = get_gb_comp(&density, &gb_tmp, PP_db, chem_comp_kos, EM_database, bulk_rock, P, T, "kos", "equilibrium");
	double gb_kos    = gb_tmp;	
	double rho_kos   = density;
	
	get_data chem_comp_jd     = get_gb_comp(&density, &gb_tmp, PP_db, chem_comp_jd, EM_database, bulk_rock, P, T, "jd", "equilibrium");
	double gb_jd     = gb_tmp;	
	double rho_jd    = density;
	
	double gb6       = gb5 + gb_kos - gb_jd + (-25.9 + 0.0155*T + 0.05*P);
	SS_ref_opx_db.density[5] = SS_ref_opx_db.density[4] + rho_kos - rho_jd;
	double chem_comp6[nEl];	
	for (i = 0; i < nEl; i++){
		chem_comp6[i] = chem_comp5.comp[i] + chem_comp_kos.comp[i] - chem_comp_jd.comp[i];
	}
	
	get_data chem_comp_per     = get_gb_comp(&density, &gb_tmp, PP_db, chem_comp_per, EM_database, bulk_rock, P, T, "per", "equilibrium");
	double gb_per     = gb_tmp;	
	double rho_per    = density;
	
	get_data chem_comp_ru     = get_gb_comp(&density, &gb_tmp, PP_db, chem_comp_ru, EM_database, bulk_rock, P, T, "ru", "equilibrium");
	double gb_ru     = gb_tmp;	
	double rho_ru    = density;
	
	get_data chem_comp_cor    = get_gb_comp(&density, &gb_tmp, PP_db, chem_comp_cor, EM_database, bulk_rock, P, T, "cor", "equilibrium");
	double gb_cor    = gb_tmp;	
	double rho_cor   = density;
	
	double gb7       = gb5 + (gb_per + gb_ru - gb_cor)/2.0 + (-5.0 - 0.0051*T - 0.0061*P);
	SS_ref_opx_db.density[6] = SS_ref_opx_db.density[4] + (rho_per + rho_ru - rho_cor)/2.0;
	
	double chem_comp7[nEl];	
	for (i = 0; i < nEl; i++){
		chem_comp7[i] = chem_comp5.comp[i] + (chem_comp_per.comp[i] + chem_comp_ru.comp[i] - chem_comp_cor.comp[i])/2.0;
	}
	
	get_data chem_comp_acm    = get_gb_comp(&density, &gb_tmp, PP_db, chem_comp_acm, EM_database, bulk_rock, P, T, "acm", "equilibrium");
	double gb_acm    = gb_tmp;
	double rho_acm   = density;	
	
	double gb8       = gb5 + gb_acm - gb_jd + (4.80 - 0.089*P);
	SS_ref_opx_db.density[7] = SS_ref_opx_db.density[4] + rho_acm - rho_jd;
	
	double chem_comp8[nEl];	
	for (i = 0; i < nEl; i++){
		chem_comp8[i] = chem_comp5.comp[i] + chem_comp_acm.comp[i] - chem_comp_jd.comp[i];
	}
	
	double gb9       = gb_jd + 18.80;
	SS_ref_opx_db.density[8] = rho_jd;
	double chem_comp9[nEl];	
	for (i = 0; i < nEl; i++){
		chem_comp9[i] = chem_comp_jd.comp[i];
	}
	
	SS_ref_opx_db.gbase[0] = gb1;
	SS_ref_opx_db.gbase[1] = gb2;
	SS_ref_opx_db.gbase[2] = gb3;
	SS_ref_opx_db.gbase[3] = gb4;
	SS_ref_opx_db.gbase[4] = gb5;
	SS_ref_opx_db.gbase[5] = gb6;
	SS_ref_opx_db.gbase[6] = gb7;
	SS_ref_opx_db.gbase[7] = gb8;
	SS_ref_opx_db.gbase[8] = gb9;
	
	for (i = 0; i < nEl; i++){
		SS_ref_opx_db.Comp[0][i] = chem_comp1.comp[i];
		SS_ref_opx_db.Comp[1][i] = chem_comp2.comp[i];
		SS_ref_opx_db.Comp[2][i] = chem_comp3[i];
		SS_ref_opx_db.Comp[3][i] = chem_comp4.comp[i];
		SS_ref_opx_db.Comp[4][i] = chem_comp5.comp[i];
		SS_ref_opx_db.Comp[5][i] = chem_comp6[i];
		SS_ref_opx_db.Comp[6][i] = chem_comp7[i];
		SS_ref_opx_db.Comp[7][i] = chem_comp8[i];
		SS_ref_opx_db.Comp[8][i] = chem_comp9[i];
	}

	/* define endmembers to be deactivated when bulk-rock composition equal zero */
	for (i = 0; i < n_em; i++){
		SS_ref_opx_db.z_em[i] = 1.0;
	}
	
	SS_ref_opx_db.box_bounds_default[0][0] = 0.+eps;  SS_ref_opx_db.box_bounds_default[0][1] = 1.-eps;
	SS_ref_opx_db.box_bounds_default[1][0] = 0.+eps;  SS_ref_opx_db.box_bounds_default[1][1] = 2.-eps;
	SS_ref_opx_db.box_bounds_default[2][0] = 0.+eps;  SS_ref_opx_db.box_bounds_default[2][1] = 1.-eps;
	SS_ref_opx_db.box_bounds_default[3][0] = -1.+eps; SS_ref_opx_db.box_bounds_default[3][1] = 1.-eps;
	SS_ref_opx_db.box_bounds_default[4][0] = 0.+eps;  SS_ref_opx_db.box_bounds_default[4][1] = 1.-eps;	
	SS_ref_opx_db.box_bounds_default[5][0] = 0.+eps;  SS_ref_opx_db.box_bounds_default[5][1] = 1.-eps;
	SS_ref_opx_db.box_bounds_default[6][0] = 0.+eps;  SS_ref_opx_db.box_bounds_default[6][1] = 1.-eps;
	SS_ref_opx_db.box_bounds_default[7][0] = 0.+eps;  SS_ref_opx_db.box_bounds_default[7][1] = 1.-eps;

	if (bulk_rock[7] == 0.){ //tio2
		SS_ref_opx_db.z_em[6]          = 0.0;
		SS_ref_opx_db.box_bounds_default[5][1] = eps;
		SS_ref_opx_db.box_bounds_default[5][0] = eps;
	}
	if (bulk_rock[8] == 0.){ //o
		SS_ref_opx_db.z_em[7]          = 0.0;
		SS_ref_opx_db.box_bounds_default[4][0] = eps;  
		SS_ref_opx_db.box_bounds_default[4][1] = eps;	
	}
	if (bulk_rock[9] == 0.){ //cr2o3
		SS_ref_opx_db.z_em[5]          = 0.0;
		SS_ref_opx_db.box_bounds_default[6][0] = eps;  
		SS_ref_opx_db.box_bounds_default[6][1] = eps;
	}

	return SS_ref_opx_db;
}

/**
  retrieve reference thermodynamic data for plagioclase4T (late 2021 update of TC, given by Eleanor)
*/
SS_ref G_SS_pl4T_function(SS_ref SS_ref_pl4T_db, int EM_database, double *bulk_rock, double P, 
								  double T, double eps){
									  	
	char   *EM_tmp[] 		= {"ab","an","san"};	
	for (int i = 0; i < SS_ref_pl4T_db.n_em; i++){ 
		strcpy(SS_ref_pl4T_db.EM_list[i],EM_tmp[i]);			
	}
			
    SS_ref_pl4T_db.W[0]    = 14.6 - 0.00935*T - 0.04*P;
    SS_ref_pl4T_db.W[1]    = 24.1 - 0.00957*T + 0.338*P;
    SS_ref_pl4T_db.W[2]    = 48.5 - 0.13*P;
    
    SS_ref_pl4T_db.v[0]    = 0.674;
    SS_ref_pl4T_db.v[1]    = 0.55;
    SS_ref_pl4T_db.v[2]    = 1.0;
							  
    PP_ref PP_db;

	double gb_tmp;
    double density;
    int i,j;
	
	int n_em = SS_ref_pl4T_db.n_em;

	get_data chem_comp1       = get_gb_comp(&density, &gb_tmp, PP_db, chem_comp1, EM_database, bulk_rock, P, T, "ab", "equilibrium");
	double gb1       = gb_tmp;	
	SS_ref_pl4T_db.density[0] = density;

	get_data chem_comp2       = get_gb_comp(&density, &gb_tmp, PP_db, chem_comp2, EM_database, bulk_rock, P, T, "an", "equilibrium");
	double gb2       = gb_tmp;	
	SS_ref_pl4T_db.density[1] = density;
						  
	get_data chem_comp3       = get_gb_comp(&density, &gb_tmp, PP_db, chem_comp3, EM_database, bulk_rock, P, T, "san", "equilibrium");
	double gb3       = gb_tmp;	
	SS_ref_pl4T_db.density[2] = density;

	SS_ref_pl4T_db.gbase[0] = gb1;
	SS_ref_pl4T_db.gbase[1] = gb2;
	SS_ref_pl4T_db.gbase[2] = gb3;
	
	for (i = 0; i < nEl; i++){
		SS_ref_pl4T_db.Comp[0][i] = chem_comp1.comp[i];
		SS_ref_pl4T_db.Comp[1][i] = chem_comp2.comp[i];
		SS_ref_pl4T_db.Comp[2][i] = chem_comp3.comp[i];
	}
	
	/* define box bounds */
	SS_ref_pl4T_db.box_bounds_default[0][0] = 0.+eps;	SS_ref_pl4T_db.box_bounds_default[0][1] = 1.-eps;
	SS_ref_pl4T_db.box_bounds_default[1][0] = 0.+eps;	SS_ref_pl4T_db.box_bounds_default[1][1] = 1.-eps;

	/* define endmembers to be deactivated when bulk-rock composition equal zero */
	for (i = 0; i < n_em; i++){
		SS_ref_pl4T_db.z_em[i] = 1.0;
	}

	return SS_ref_pl4T_db;										  
}

/**
  retrieve reference thermodynamic data for spinel
*/
SS_ref G_SS_spn_function(SS_ref SS_ref_spn_db, int EM_database, double *bulk_rock, double P, 
								  double T, double eps){	
									  	
	char   *EM_tmp[] 		= {"nsp","isp","nhc","ihc","nmt","imt","pcr","qndm"};	
	for (int i = 0; i < SS_ref_spn_db.n_em; i++){ 
		strcpy(SS_ref_spn_db.EM_list[i],EM_tmp[i]);			
	}
									  
    SS_ref_spn_db.W[0] = -8.2;
    SS_ref_spn_db.W[1] = 3.5;
    SS_ref_spn_db.W[2] = -13.0;
    SS_ref_spn_db.W[3] = 43.2;
    SS_ref_spn_db.W[4] = 49.1;
    SS_ref_spn_db.W[5] = -5.0;
    SS_ref_spn_db.W[6] = 22.5;
    SS_ref_spn_db.W[7] = 4.4;
    SS_ref_spn_db.W[8] = -6.0;
    SS_ref_spn_db.W[9] = 36.8;
    SS_ref_spn_db.W[10] = 20.0;
    SS_ref_spn_db.W[11] = 14.0;
    SS_ref_spn_db.W[12] = 21.5;
    SS_ref_spn_db.W[13] = -8.2;
    SS_ref_spn_db.W[14] = 18.1;
    SS_ref_spn_db.W[15] = 49.0;
    SS_ref_spn_db.W[16] = -19.0;
    SS_ref_spn_db.W[17] = 35.1;
    SS_ref_spn_db.W[18] = -4.0;
    SS_ref_spn_db.W[19] = 7.6;
    SS_ref_spn_db.W[20] = -11.0;
    SS_ref_spn_db.W[21] = 9.0;
    SS_ref_spn_db.W[22] = 18.1;
    SS_ref_spn_db.W[23] = 11.9;
    SS_ref_spn_db.W[24] = 62.2;
    SS_ref_spn_db.W[25] = -6.4;
    SS_ref_spn_db.W[26] = 24.3;
    SS_ref_spn_db.W[27] = 60.0;					  	
									  
    PP_ref PP_db;

	double gb_tmp;
    double density;
    int i,j;
	
	int n_em = SS_ref_spn_db.n_em;
	
	get_data chem_comp1       = get_gb_comp(&density, &gb_tmp, PP_db, chem_comp1, EM_database, bulk_rock, P, T, "sp", "ordered");
	double gb1       = gb_tmp;	
	SS_ref_spn_db.density[0] = density;
	
	get_data chem_comp2       = get_gb_comp(&density, &gb_tmp, PP_db, chem_comp2, EM_database, bulk_rock, P, T, "sp", "ordered");
	double gb2       = gb_tmp + 23.6 - 0.00576303*T;
	SS_ref_spn_db.density[1] = density;
	
	get_data chem_comp3       = get_gb_comp(&density, &gb_tmp, PP_db, chem_comp3, EM_database, bulk_rock, P, T, "herc", "ordered");
	double gb3       = gb_tmp;
	SS_ref_spn_db.density[2] = density;
		
	get_data chem_comp4       = get_gb_comp(&density, &gb_tmp, PP_db, chem_comp4, EM_database, bulk_rock, P, T, "herc", "ordered");
	double gb4       = gb_tmp + 23.60 - 0.00576303*T;		
	SS_ref_spn_db.density[3] = density;
		
	get_data chem_comp5       = get_gb_comp(&density, &gb_tmp, PP_db, chem_comp5, EM_database, bulk_rock, P, T, "mt", "equilibrium");
	double gb5       = gb_tmp + 0.00576303*T;	
	SS_ref_spn_db.density[4] = density;	
	
	get_data chem_comp6       = get_gb_comp(&density, &gb_tmp, PP_db, chem_comp6, EM_database, bulk_rock, P, T, "mt", "equilibrium");
	double gb6       = gb_tmp + 0.3;
	SS_ref_spn_db.density[5] = density;
	
	get_data chem_comp7       = get_gb_comp(&density, &gb_tmp, PP_db, chem_comp7, EM_database, bulk_rock, P, T, "picr", "equilibrium"); //ordered
	double gb7       = gb_tmp;	
	SS_ref_spn_db.density[6] = density;
	
	get_data chem_comp8       = get_gb_comp(&density, &gb_tmp, PP_db, chem_comp8, EM_database, bulk_rock, P, T, "qnd", "equilibrium");
	double gb8       = gb_tmp -30.0;			
	SS_ref_spn_db.density[7] = density;
	
	SS_ref_spn_db.gbase[0] = gb1;
	SS_ref_spn_db.gbase[1] = gb2;
	SS_ref_spn_db.gbase[2] = gb3;
	SS_ref_spn_db.gbase[3] = gb4;
	SS_ref_spn_db.gbase[4] = gb5;
	SS_ref_spn_db.gbase[5] = gb6;
	SS_ref_spn_db.gbase[6] = gb7;
	SS_ref_spn_db.gbase[7] = gb8;
	
	for (i = 0; i < nEl; i++){
		SS_ref_spn_db.Comp[0][i] = chem_comp1.comp[i];
		SS_ref_spn_db.Comp[1][i] = chem_comp2.comp[i];
		SS_ref_spn_db.Comp[2][i] = chem_comp3.comp[i];
		SS_ref_spn_db.Comp[3][i] = chem_comp4.comp[i];
		SS_ref_spn_db.Comp[4][i] = chem_comp5.comp[i];
		SS_ref_spn_db.Comp[5][i] = chem_comp6.comp[i];
		SS_ref_spn_db.Comp[6][i] = chem_comp7.comp[i];
		SS_ref_spn_db.Comp[7][i] = chem_comp8.comp[i];
	}

	/* define endmembers to be deactivated when bulk-rock composition equal zero */
	for (i = 0; i < n_em; i++){
		SS_ref_spn_db.z_em[i] = 1.0;
	}
	SS_ref_spn_db.box_bounds_default[0][0] =  0.0 +eps;   SS_ref_spn_db.box_bounds_default[0][1] = 1.0 -eps;
	SS_ref_spn_db.box_bounds_default[1][0] =  0.0 +eps;   SS_ref_spn_db.box_bounds_default[1][1] = 1.0 -eps;
	SS_ref_spn_db.box_bounds_default[2][0] =  0.0 +eps;   SS_ref_spn_db.box_bounds_default[2][1] = 1.0 -eps;
	SS_ref_spn_db.box_bounds_default[3][0] =  0.0 +eps;   SS_ref_spn_db.box_bounds_default[3][1] = 1.0 -eps;
	SS_ref_spn_db.box_bounds_default[4][0] = -1.0 +eps;   SS_ref_spn_db.box_bounds_default[4][1] = 1.0 -eps;	
	SS_ref_spn_db.box_bounds_default[5][0] = -1.0 +eps;   SS_ref_spn_db.box_bounds_default[5][1] = 1.0 -eps;
	SS_ref_spn_db.box_bounds_default[6][0] = -1.0 +eps;   SS_ref_spn_db.box_bounds_default[6][1] = 1.0 -eps;

	if (bulk_rock[7] == 0.){ //tio2
		SS_ref_spn_db.z_em[7]          = 0.0;
		SS_ref_spn_db.box_bounds_default[3][0] = eps; 
		SS_ref_spn_db.box_bounds_default[3][1] = eps;
	}
	if (bulk_rock[8] == 0.){ //fe2o3
		SS_ref_spn_db.box_bounds_default[1][0] = eps;   
		SS_ref_spn_db.box_bounds_default[1][1] = eps;
		SS_ref_spn_db.box_bounds_default[6][0] = eps;   
		SS_ref_spn_db.box_bounds_default[6][1] = eps;
		SS_ref_spn_db.z_em[4]          = 0.0;
		SS_ref_spn_db.z_em[5]          = 0.0;
	}
	if (bulk_rock[9] == 0.){ //cr2o3
		SS_ref_spn_db.z_em[6]          = 0.0;
		SS_ref_spn_db.box_bounds_default[2][0] = eps;  
		SS_ref_spn_db.box_bounds_default[2][1] = eps;
	}

	return SS_ref_spn_db;
}

/**
  checks if it can satisfy the mass constraint
*/
SS_ref G_SS_EM_function(		global_variable 	 gv,
								SS_ref 				 SS_ref_db, 
								int 				 EM_database, 
								struct 	bulk_info 	 z_b, 
								char   				*name				){
									  
	double eps 		   	= gv.bnd_val;
	double P 			= z_b.P;
	double T 			= z_b.T;	
					   
	SS_ref_db.ss_flags[0]  = 1;

	/* Associate the right solid-solution data */
	for (int FD = 0; FD < gv.n_Diff; FD++){				/* cycle twice in order to get gb_P_eps to calculate densities later on */
		
		P = z_b.P + gv.gb_P_eps*gv.numDiff[0][FD];
		T = z_b.T + gv.gb_T_eps*gv.numDiff[1][FD];

		if (strcmp( name, "bi") == 0 ){
			// if no H2O, deactivate
			if (z_b.bulk_rock[10] == 0.){
				SS_ref_db.ss_flags[0]  = 0;
			}
			SS_ref_db  = G_SS_bi_function(SS_ref_db, EM_database, z_b.bulk_rock, P, T, eps);	}
		else if (strcmp( name, "cd") == 0){
			// if no H2O, deactivate
			if (z_b.bulk_rock[10] == 0.){
				SS_ref_db.ss_flags[0]  = 0;
			}
			SS_ref_db  = G_SS_cd_function(SS_ref_db, EM_database, z_b.bulk_rock, P, T, eps);	}
		else if (strcmp( name, "cpx") == 0){
			SS_ref_db  = G_SS_cpx_function(SS_ref_db, EM_database, z_b.bulk_rock, P, T, eps);	}	
		else if (strcmp( name, "ep") == 0){
			// if no h2O, deactivate
			if (z_b.bulk_rock[8] == 0. || z_b.bulk_rock[10] == 0.){
				SS_ref_db.ss_flags[0]  = 0;
			}
			SS_ref_db  = G_SS_ep_function(SS_ref_db, EM_database, z_b.bulk_rock, P, T, eps);	}
		else if (strcmp( name, "fl") == 0){
			// if no H2O, deactivate
			if (z_b.bulk_rock[10] == 0.){
				SS_ref_db.ss_flags[0]  = 0;
			}
			SS_ref_db  = G_SS_fl_function(SS_ref_db, EM_database, z_b.bulk_rock, P, T, eps);	}		
		else if (strcmp( name, "g") == 0){
			SS_ref_db  = G_SS_g_function(SS_ref_db, EM_database, z_b.bulk_rock, P, T, eps);		}
		else if (strcmp( name, "hb") == 0){
			// if no H2O, deactivate
			if (z_b.bulk_rock[10] == 0.){
				SS_ref_db.ss_flags[0]  = 0;
			}
			SS_ref_db  = G_SS_hb_function(SS_ref_db, EM_database, z_b.bulk_rock, P, T, eps);	}	
		else if (strcmp( name, "ilm") == 0){
			// if no TiO2, deactivate
			if (z_b.bulk_rock[7] == 0.){
				SS_ref_db.ss_flags[0]  = 0;
			}
			SS_ref_db  = G_SS_ilm_function(SS_ref_db, EM_database, z_b.bulk_rock, P, T, eps);	}
		else if (strcmp( name, "liq") == 0){
			/* turn of liquid when T < 500°C) */
			if ( T < 773.0){
				SS_ref_db.ss_flags[0]  = 0;
			}
			SS_ref_db = G_SS_liq_function(SS_ref_db, EM_database, z_b.bulk_rock, P, T, eps);	}
		else if (strcmp( name, "mu") == 0){
			// if no H2O, deactivate
			if (z_b.bulk_rock[10] == 0.){
				SS_ref_db.ss_flags[0]  = 0;
			}
			SS_ref_db  = G_SS_mu_function(SS_ref_db, EM_database, z_b.bulk_rock, P, T, eps);	}	
		else if (strcmp( name, "ol") == 0){
			SS_ref_db  = G_SS_ol_function(SS_ref_db, EM_database, z_b.bulk_rock, P, T, eps);	}
		else if (strcmp( name, "opx") == 0){
			SS_ref_db  = G_SS_opx_function(SS_ref_db, EM_database, z_b.bulk_rock, P, T, eps);	}
		else if (strcmp( name, "pl4T") == 0){
			SS_ref_db  = G_SS_pl4T_function(SS_ref_db, EM_database, z_b.bulk_rock, P, T, eps);	}	
		else if (strcmp( name, "spn") == 0){
			SS_ref_db  = G_SS_spn_function(SS_ref_db, EM_database, z_b.bulk_rock, P, T, eps);	}
		else{
			printf("\nsolid solution '%s' is not in the database\n",name);	}	
		
		for (int j = 0; j < SS_ref_db.n_em; j++){
			SS_ref_db.mu_array[FD][j] = SS_ref_db.gbase[j];
		}
	}

	for (int j = 0; j < SS_ref_db.n_xeos; j++){
		SS_ref_db.box_bounds[j][0] = SS_ref_db.box_bounds_default[j][0];
		SS_ref_db.box_bounds[j][1] = SS_ref_db.box_bounds_default[j][1];
	}

	/* Calculate the number of atoms in the bulk-rock composition */
	double fbc     = 0.0;
	for (int i = 0; i < gv.len_ox; i++){
		fbc += z_b.bulk_rock[i]*z_b.apo[i];
	}

	/* get the numer of atoms per endmember, needed to update normalization factor for liquid */
	for (int i = 0; i < SS_ref_db.n_em; i++){
		SS_ref_db.ape[i] = 0.0;
		for (int j = 0; j < nEl; j++){
			SS_ref_db.ape[i] += SS_ref_db.Comp[i][j]*z_b.apo[j];
		}
	}
	
	SS_ref_db.fbc = z_b.fbc;	
	
	if (gv.verbose == 1){
		printf(" %4s:",name);
		for (int j = 0; j < SS_ref_db.n_em; j++){
			printf(" %+12.5f",SS_ref_db.gbase[j]);
		}
		for (int j = SS_ref_db.n_em; j < gv.len_ox; j++){
			printf("%13s","-");
		}
		printf("\n");
	}
				
		
	return SS_ref_db;
};

