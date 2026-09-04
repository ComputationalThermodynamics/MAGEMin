/*@ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 **
 **   Project      : MAGEMin
 **   License      : GNU GENERAL PUBLIC LICENSE Version 3, 29 June 2007
 **   Developers   : Nicolas Riel, Boris Kaus
 **   Contributors : Moccetti, N. B., Dominguez, H., Assunção J., Green E., Dolejš, D., Berlie N., and Rummel L.
 **   Organization : Institute of Geosciences, Johannes-Gutenberg University, Mainz
 **   Contact      : nriel[at]uni-mainz.de, kaus[at]uni-mainz.de
 **
 ** ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ @*/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <complex.h>

#include "../uthash.h"
#include "../MAGEMin.h"
#include "../initialize.h"
#include "../all_endmembers.h"
#include "GH_init_database.h"

/* Same 16-oxide axis table as TC_database/TC_init_database.c's oxide_info,
   kept as a research-group-local copy (same convention SB_database uses
   for oxide_info_sb) so "gh" doesn't depend on tc internals. */
oxide_data oxide_info_gh = {
    16,
    {"SiO2","Al2O3","CaO","MgO","FeO","K2O","Na2O","TiO2","O","MnO","Cr2O3","H2O","CO2","S","Cl","ecp"},
    {"Si","Al","Ca","Mg","Fe","K","Na","Ti","O","Mn","Cr","H","C","S","Cl","ecp"},
    {60.08,101.96,56.08,40.30,71.85,94.2,61.98,79.88,16.0,70.94,151.99,18.015,44.01,32.06,35.453,0.0},
    {3.0,5.0,2.0,2.0,2.0,3.0,3.0,3.0,1.0,2.0,5.0,3.0,3.0,1.0,1.0,0.0},
    {66.7736,108.653,42.9947,40.3262,38.7162,69.1514,61.1729,70.3246,30.5827,40.1891,106.9795,69.5449,62.8768,9.5557,33.2556,0.0},
    {2,3,1,1,1,1,1,2,1,1,3,1,2,0,0,0},
    {1,2,1,1,1,2,2,1,1,1,2,2,1,1,1,0}
};

gh_dataset gh_db = {
	1,							/* ds_version (single Stage-A dataset)				*/
	13,							/* number of oxides									*/
	25,							/* number of pure phases							*/
	16,							/* number of solution phases						*/
	{"SiO2"	,"Al2O3","CaO"	,"MgO"	,"FeO"	,"K2O"	,"Na2O"	,"TiO2"	,"O"	,"MnO"	,"Cr2O3","H2O"	,"CO2"				            },
	{"aTiO2","aAl2O3","qfm" ,"hm"   ,"q"	,"crst"	,"trd"	,"cor"	,"sill"	,"and"	,"ky"	,"ru"	,"sph"	,"perov","cc"	,"arag"	,"mgs"	,"sid"	,
	 "dol"	,"spu"	,"til"	,"mu"	,"aeg"	,"aen"	,"O2"	/* ,"H2O"*/															        },
	{"liq"	,"ol"	,"fsp"	,"bi"	,"g"	,"hb"	,"lc"	,"mel"	,"cum"	,"spn"	,"cpx"	,"opx"	,"fl"	,"rhm"	,"nph"	,"kls"	},

	{1		,1		,1		,1		,1		,1		,1		,1		,1		,1		,1		,1		,1		,1		,1		,1		}, // allow solvus?
	{1820	,41		,231	,41		,231	,231	,231	,165	,41		,1001	,3003	,3003	,41		,1820	,165	,165	}, // # of pseudocompound
	{0.15	,0.025	,0.05	,0.05	,0.05	,0.05	,0.05	,0.10	,0.025	,0.05	,0.10	,0.10	,0.025	,0.05	,0.125	,0.125	}, // discretization step

	6.0, 						/** max dG under which a phase is considered to be reintroduced  					*/
	673.15,						/** max temperature above which PGE solver is active 								*/
	773.15,						/** minimum temperature above which melt is considered 								*/

	4,							/** number of inner PGE iterations, this has to be made mass or dG dependent 		*/
	0.025,						/** maximum mol% phase change during one PGE iteration in wt% 						*/
	2.5,						/** maximum delta_G of reference change during PGE 									*/
	1.0,						/** maximum update factor during PGE under-relax < 0.0, over-relax > 0.0 	 		*/

	1e-1,						/** merge instances of solution phase if norm < val 								*/
	1e-4,						/** fraction of solution phase when re-introduced 									*/
	1e-6						/** objective function tolerance 				 									*/
};

/** pMELTS (gv.EM_database==2): pMELTS minus CO2 entirely  */
gh_dataset gh_db_pmelts_dataset = {
	1,							/* ds_version (single Stage-A dataset)				*/
	12,							/* number of oxides									*/
	19,							/* number of pure phases							*/
	15,							/* number of solution phases						*/
	{"SiO2"	,"Al2O3","CaO"	,"MgO"	,"FeO"	,"K2O"	,"Na2O"	,"TiO2"	,"O"	,"MnO"	,"Cr2O3","H2O"						            },
	{"aTiO2","aAl2O3","qfm" ,"hm"   ,"q"	,"crst"	,"trd"	,"cor"	,"sill"	,"and"	,"ky"	,"ru"	,"sph"	,"perov","mu"	,"aeg"	,"aen"	,"O2"	,"H2O"	},
	{"liq"	,"ol"	,"fsp"	,"bi"	,"g"	,"hb"	,"lc"	,"mel"	,"cum"	,"spn"	,"cpx"	,"opx"	,"rhm"	,"nph"	,"kls"	},

	{1		,1		,1		,1		,1		,1		,1		,1		,1		,1		,1		,1		,1		,1		,1		}, // allow solvus?
	{4368	,41		,231	,41		,231	,231	,231	,165	,41		,1001	,3003	,3003	,1820	,165	,165	}, // # of pseudocompound
	{0.15	,0.025	,0.05	,0.05	,0.05	,0.05	,0.05	,0.10	,0.025	,0.05	,0.10	,0.10	,0.05	,0.125	,0.125	}, // discretization step

	6.0, 						/** max dG under which a phase is considered to be reintroduced  					*/
	673.15,						/** max temperature above which PGE solver is active 								*/
	773.15,						/** minimum temperature above which melt is considered 								*/

	4,							/** number of inner PGE iterations, this has to be made mass or dG dependent 		*/
	0.025,						/** maximum mol% phase change during one PGE iteration in wt% 						*/
	2.5,						/** maximum delta_G of reference change during PGE 									*/
	1.0,						/** maximum update factor during PGE under-relax < 0.0, over-relax > 0.0 	 		*/

	1e-1,						/** merge instances of solution phase if norm < val 								*/
	1e-4,						/** fraction of solution phase when re-introduced 									*/
	1e-6						/** objective function tolerance 				 									*/
};

global_variable global_variable_GH_init(   global_variable      gv,
                                            bulk_info           *z_b    ){
    int i, j;

    gh_dataset db       = (gv.EM_database == 2) ? gh_db_pmelts_dataset : gh_db;
    gv.EM_dataset       = db.ds_version;
    gv.len_pp           = db.n_pp + gv.n_mu_fix;
    gv.len_ss           = db.n_ss;
    gv.len_ox           = db.n_ox;

    gv.PC_df_add        = db.PC_df_add;
    gv.solver_switch_T  = db.solver_switch_T;
    gv.min_melt_T       = db.min_melt_T;

    gv.inner_PGE_ite    = db.inner_PGE_ite;
    gv.max_n_phase      = db.max_n_phase;
    gv.max_g_phase      = db.max_g_phase;
    gv.max_fac          = db.max_fac;

    gv.merge_value      = db.merge_value;
    gv.re_in_n          = db.re_in_n;
    gv.obj_tol          = db.obj_tol;

    gv.ox               = malloc (gv.len_ox * sizeof(char*)   );
    for (i = 0; i < gv.len_ox; i++){
        gv.ox[i]        = malloc(20 * sizeof(char));
        strcpy(gv.ox[i],db.ox[i]);
    }

    gv.PP_list          = malloc ((gv.len_pp>0?gv.len_pp:1) * sizeof(char*)    );
    for (i = 0; i < db.n_pp; i++){
        gv.PP_list[i]   = malloc(20 * sizeof(char));
        strcpy(gv.PP_list[i],db.PP[i]);
    }
    for (i = 0; i < gv.n_mu_fix; i++){
        gv.PP_list[db.n_pp+i] = malloc(20 * sizeof(char));
        sprintf(gv.PP_list[db.n_pp+i], "mu%d", i);
    }

    gv.SS_list          = malloc ((gv.len_ss) * sizeof (char*) );
    gv.n_SS_PC          = malloc ((gv.len_ss) * sizeof (int)   );
    gv.verifyPC         = malloc ((gv.len_ss) * sizeof (int)   );
    gv.SS_PC_stp        = malloc ((gv.len_ss) * sizeof (double));
    for (i = 0; i < gv.len_ss; i++){
        gv.SS_list[i]   = malloc(20 * sizeof(char)             );
        strcpy(gv.SS_list[i],db.SS[i]);
        gv.verifyPC[i]  = db.verifyPC[i];
        gv.n_SS_PC[i]   = db.n_SS_PC[i];
        gv.SS_PC_stp[i] = db.SS_PC_stp[i];
    }
    gv.act_PP           = malloc ((gv.len_pp>0?gv.len_pp:1) * sizeof (int) );
    for (i = 0; i < gv.len_pp; i++){ gv.act_PP[i] = 1; };

    /**
        ALLOCATE MEMORY OF OTHER GLOBAL VARIABLES
        (mirrors global_variable_SB_init's generic PGE/LP working-memory
        footprint, identical regardless of which oxide/phase list is active)
    */
    gv.n_min            = malloc ((gv.len_ss) * sizeof (int)  );
    gv.n_ss_ph          = malloc ((gv.len_ss) * sizeof (int)  );
    gv.bulk_rock        = malloc (gv.len_ox * sizeof(double)  );
    gv.PGE_mass_norm    = malloc (gv.it_f*2 * sizeof (double) );
    gv.Alg              = malloc (gv.it_f*2 * sizeof (int)    );
    gv.gamma_norm       = malloc (gv.it_f*2 * sizeof (double) );
    gv.gibbs_ev         = malloc (gv.it_f*2 * sizeof (double) );
    gv.ite_time         = malloc (gv.it_f*2 * sizeof (double) );

    gv.n_Diff = 8;
    gv.pdev = malloc (2 * sizeof(double*));
    for (i = 0; i < 2; i++){
        gv.pdev[i] = malloc (gv.n_Diff * sizeof(double));
    }
    gv.pdev[0][0]  =  0.0;  gv.pdev[1][0]  =  1.0;
    gv.pdev[0][1]  =  0.0;  gv.pdev[1][1]  = -1.0;
    gv.pdev[0][2]  =  1.0;  gv.pdev[1][2]  =  1.0;
    gv.pdev[0][3]  =  1.0;  gv.pdev[1][3]  = -1.0;
    gv.pdev[0][4]  =  2.0;  gv.pdev[1][4]  =  0.0;
    gv.pdev[0][5]  =  1.0;  gv.pdev[1][5]  =  0.0;
    gv.pdev[0][6]  =  3.0;  gv.pdev[1][6]  =  0.0;
    gv.pdev[0][7]  =  0.0;  gv.pdev[1][7]  =  0.0;

    gv.V_cor = malloc (2 * sizeof(double));

    gv.dGamma           = malloc (gv.len_ox * sizeof(double)   );
    gv.gam_tot          = malloc (gv.len_ox * sizeof (double)  );
    gv.gam_tot_0        = malloc (gv.len_ox * sizeof (double)  );
    gv.delta_gam_tot    = malloc (gv.len_ox * sizeof (double)  );
    gv.mass_residual    = malloc (gv.len_ox * sizeof(double)   );

    gv.lwork            = 64;
    gv.ipiv             = malloc ((gv.len_ox*3) * sizeof (int)    );
    gv.work             = malloc ((gv.len_ox*gv.lwork) * sizeof (double)  );
    gv.n_solvi          = malloc ((gv.len_ss) * sizeof (int)   );

    gv.n_flags     = 5;

    gv.pp_n             = malloc ((gv.len_pp>0?gv.len_pp:1) * sizeof(double)   );
    gv.pp_n_mol         = malloc ((gv.len_pp>0?gv.len_pp:1) * sizeof(double)   );
    gv.pp_n_wt          = malloc ((gv.len_pp>0?gv.len_pp:1) * sizeof(double)   );
    gv.pp_n_vol         = malloc ((gv.len_pp>0?gv.len_pp:1) * sizeof(double)   );
    gv.pp_xi            = malloc ((gv.len_pp>0?gv.len_pp:1) * sizeof(double)   );
    gv.delta_pp_n       = malloc ((gv.len_pp>0?gv.len_pp:1) * sizeof(double)   );
    gv.delta_pp_xi      = malloc ((gv.len_pp>0?gv.len_pp:1) * sizeof(double)   );
    gv.pp_flags         = malloc ((gv.len_pp>0?gv.len_pp:1) * sizeof(int*)     );

    for (i = 0; i < (gv.len_pp>0?gv.len_pp:1); i++){
        gv.pp_flags[i]      = malloc (gv.n_flags  * sizeof(int));
    }

    gv.A_PGE  = malloc ((gv.len_ox*gv.len_ox*9)    * sizeof(double));
    gv.A0_PGE = malloc ((gv.len_ox*gv.len_ox*9)    * sizeof(double));
    gv.b_PGE  = malloc ((gv.len_ox*3)              * sizeof(double));

    gv.cp_id  = malloc ((gv.len_ox)                * sizeof(int)  );
    gv.pp_id  = malloc ((gv.len_ox)                * sizeof(int)  );

    gv.dn_cp  = malloc ((gv.len_ox)                * sizeof(double));
    gv.dn_pp  = malloc ((gv.len_ox)                * sizeof(double));

    gv.A  = malloc ((gv.len_ox) * sizeof(double*));
    gv.A2 = malloc ((gv.len_ox) * sizeof(double*));
    for (i = 0; i < (gv.len_ox); i++){
        gv.A[i]  = malloc ((gv.len_ox) * sizeof(double));
        gv.A2[i] = malloc ((gv.len_ox) * sizeof(double));

        gv.pc_id = malloc (gv.len_ox * sizeof(int));
    }
    gv.b    = malloc (gv.len_ox * sizeof(double));
    gv.b1   = malloc (gv.len_ox * sizeof(double));
    gv.tmp1 = malloc (gv.len_ox * sizeof(double));
    gv.tmp2 = malloc (gv.len_ox * sizeof(double));
    gv.tmp3 = malloc (gv.len_ox * sizeof(double));
    gv.n_ss_array = malloc (gv.len_ss* sizeof(double));

    z_b->apo            = malloc (gv.len_ox * sizeof (double) );
    z_b->masspo         = malloc (gv.len_ox * sizeof (double) );
    z_b->opo            = malloc (gv.len_ox * sizeof (double) );
    z_b->cpo            = malloc (gv.len_ox * sizeof (double) );
    z_b->ElEntropy      = malloc (gv.len_ox * sizeof (double) );
    z_b->id             = malloc (gv.len_ox * sizeof (int)    );
    z_b->elName         = malloc (gv.len_ox * sizeof (char*)  );
    for (i = 0; i < (gv.len_ox); i++){
        z_b->elName[i]  = malloc(20 * sizeof(char));
    }

    gv.Al2O3_id = -1;
    gv.TiO2_id  = -1;
    gv.CaO_id   = -1;
    gv.Na2O_id  = -1;
    gv.FeO_id   = -1;
    gv.Fe_id    = -1;
    gv.MgO_id   = -1;
    gv.Cr2O3_id = -1;
    gv.O_id     = -1;
    gv.MnO_id   = -1;
    gv.CO2_id   = -1;

    /* z_b mirrors gv's oxide ids so G_SS_*_function() (which only receives z_b, not gv)
       can test a specific oxide's abundance without hardcoding its positional index */
    z_b->Al2O3_id = -1;
    z_b->TiO2_id  = -1;
    z_b->CaO_id   = -1;
    z_b->Na2O_id  = -1;
    z_b->FeO_id   = -1;
    z_b->MgO_id   = -1;
    z_b->Cr2O3_id = -1;
    z_b->O_id     = -1;
    z_b->MnO_id   = -1;
    z_b->CO2_id   = -1;

    oxide_data ox_in = oxide_info_gh;
    for (i = 0; i < gv.len_ox; i++){
        for (j = 0; j < ox_in.n_ox; j++){
            if (strcmp( gv.ox[i], ox_in.oxName[j]) == 0){
                if      (strcmp( gv.ox[i], "Al2O3") == 0){ gv.Al2O3_id = i; z_b->Al2O3_id = i; }
                else if (strcmp( gv.ox[i], "TiO2")  == 0){ gv.TiO2_id  = i; z_b->TiO2_id  = i; }
                else if (strcmp( gv.ox[i], "O")     == 0){ gv.O_id     = i; z_b->O_id     = i; }
                else if (strcmp( gv.ox[i], "CaO")   == 0){ gv.CaO_id   = i; z_b->CaO_id   = i; }
                else if (strcmp( gv.ox[i], "Na2O")  == 0){ gv.Na2O_id  = i; z_b->Na2O_id  = i; }
                else if (strcmp( gv.ox[i], "MgO")   == 0){ gv.MgO_id   = i; z_b->MgO_id   = i; }
                else if (strcmp( gv.ox[i], "FeO")   == 0){ gv.FeO_id   = i; z_b->FeO_id   = i; }
                else if (strcmp( gv.ox[i], "Cr2O3") == 0){ gv.Cr2O3_id = i; z_b->Cr2O3_id = i; }
                else if (strcmp( gv.ox[i], "MnO")   == 0){ gv.MnO_id   = i; z_b->MnO_id   = i; }
                else if (strcmp( gv.ox[i], "CO2")   == 0){ gv.CO2_id   = i; z_b->CO2_id   = i; }
                z_b->apo[i]         = ox_in.atPerOx[j];
                z_b->masspo[i]      = ox_in.oxMass[j];
                z_b->opo[i]         = ox_in.OPerOx[j];
                z_b->cpo[i]         = ox_in.catPerOx[j];
                z_b->ElEntropy[i]   = ox_in.ElEntropy[j];
                strcpy(z_b->elName[i],ox_in.elName[j]);
                z_b->id[i]          = j;
                break;
            }
        }
    }
    z_b->bulk_rock_cat  = malloc (gv.len_ox * sizeof (double) );
    z_b->bulk_rock      = malloc (gv.len_ox * sizeof (double) );
    z_b->nzEl_array     = malloc (gv.len_ox * sizeof (int) );
    z_b->zEl_array      = malloc (gv.len_ox * sizeof (int) );

    gv.n_em_db = (gv.EM_database == 2) ? 12 : 13;

    return gv;
}

/* Simple test bulk-rock compositions for the "gh" beta test
   - not meant to be geologically precise, only to exercise the liquid
   model's minimization pipeline. Values are oxide mol% (sum to ~100),
   matching the convention tc's get_bulk_igneous test bulks (KLB1, RE46,
   ...) use - MAGEMin's default --sys_in is "mol", so raw gv.bulk_rock
   entries are used as-is (no wt%->mol conversion applied, see
   toolkit.c's read_in_bulk); wt%-style figures would need dividing
   through by oxide molar mass first. */
global_variable get_bulk_gh( global_variable gv) {
    if (gv.test != -1){
        if (gv.verbose == 1){
            printf("\n");
            printf("   - Minimization using in-built bulk-rock  : test %2d\n",gv.test);
        }
    }
    else{
        gv.test = 0;
    }

    if (gv.test == 0){ /* rough basalt */
        /* SiO2  Al2O3  CaO   MgO   FeO  K2O  Na2O  TiO2  O    MnO   Cr2O3 H2O  CO2 */
        gv.bulk_rock[0]  = 49.90;   /** SiO2  */
        gv.bulk_rock[1]  = 8.82;    /** Al2O3 */
        gv.bulk_rock[2]  = 10.69;   /** CaO   */
        gv.bulk_rock[3]  = 11.90;   /** MgO   */
        gv.bulk_rock[4]  = 7.51;    /** FeO   */
        gv.bulk_rock[5]  = 0.32;    /** K2O   */
        gv.bulk_rock[6]  = 2.42;    /** Na2O  */
        gv.bulk_rock[7]  = 1.13;    /** TiO2  */
        gv.bulk_rock[8]  = 0.37;    /** O     */
        gv.bulk_rock[9]  = 0.13;    /** MnO   */
        gv.bulk_rock[10] = 0.02;    /** Cr2O3 */
        gv.bulk_rock[11] = 6.66;    /** H2O   */
        gv.bulk_rock[12] = 0.14;    /** CO2   */
    }
    else if (gv.test == 1){ /* rough rhyolite */
        /* SiO2  Al2O3  CaO   MgO   FeO  K2O  Na2O  TiO2  O    MnO   Cr2O3 H2O  CO2 */
        gv.bulk_rock[0]  = 76.58;   /** SiO2  */
        gv.bulk_rock[1]  = 7.93;    /** Al2O3 */
        gv.bulk_rock[2]  = 1.33;    /** CaO   */
        gv.bulk_rock[3]  = 0.46;    /** MgO   */
        gv.bulk_rock[4]  = 1.30;    /** FeO   */
        gv.bulk_rock[5]  = 2.97;    /** K2O   */
        gv.bulk_rock[6]  = 3.51;    /** Na2O  */
        gv.bulk_rock[7]  = 0.16;    /** TiO2  */
        gv.bulk_rock[8]  = 0.39;    /** O     */
        gv.bulk_rock[9]  = 0.04;    /** MnO   */
        gv.bulk_rock[10] = 0.4;   /** Cr2O3 */
        gv.bulk_rock[11] = 50.0;    /** H2O   */
        gv.bulk_rock[12] = 2.0;    /** CO2   */
    }
    else if (gv.test == 2){ /* rough rhyolite */
        /* SiO2  Al2O3  CaO   MgO   FeO  K2O  Na2O  TiO2  O    MnO   Cr2O3 H2O  CO2 */
        gv.bulk_rock[0]  = 76.58;   /** SiO2  */
        gv.bulk_rock[1]  = 7.93;    /** Al2O3 */
        gv.bulk_rock[2]  = 1.33;    /** CaO   */
        gv.bulk_rock[3]  = 0.46;    /** MgO   */
        gv.bulk_rock[4]  = 1.30;    /** FeO   */
        gv.bulk_rock[5]  = 2.97;    /** K2O   */
        gv.bulk_rock[6]  = 3.51;    /** Na2O  */
        gv.bulk_rock[7]  = 0.16;    /** TiO2  */
        gv.bulk_rock[8]  = 0.0;    /** O     */
        gv.bulk_rock[9]  = 0.0;    /** MnO   */
        gv.bulk_rock[10] = 0.00;   /** Cr2O3 */
        gv.bulk_rock[11] = 30.18;    /** H2O   */
        gv.bulk_rock[12] = 0.0;    /** CO2   */
    }
    else{
        printf("Unknown test %i for 'gh' research group - please specify a different test! \n", gv.test);
        exit(EXIT_FAILURE);
    }
    return gv;
}

/* pMELTS (gv.EM_database==2) test bulk-rock compositions: identical to
   get_bulk_gh's above, minus the CO2 entry (index 12) - pMELTS drops CO2
   entirely as an oxide axis (see gh_db_pmelts_dataset), so this dataset's
   gv.bulk_rock only has 12 slots. */
global_variable get_bulk_pmelts_dataset( global_variable gv) {
    if (gv.test != -1){
        if (gv.verbose == 1){
            printf("\n");
            printf("   - Minimization using in-built bulk-rock  : test %2d\n",gv.test);
        }
    }
    else{
        gv.test = 0;
    }

    if (gv.test == 0){ /* rough basalt */
        /* SiO2  Al2O3  CaO   MgO   FeO  K2O  Na2O  TiO2  O    MnO   Cr2O3 H2O */
        gv.bulk_rock[0]  = 49.90;   /** SiO2  */
        gv.bulk_rock[1]  = 8.82;    /** Al2O3 */
        gv.bulk_rock[2]  = 10.69;   /** CaO   */
        gv.bulk_rock[3]  = 11.90;   /** MgO   */
        gv.bulk_rock[4]  = 7.51;    /** FeO   */
        gv.bulk_rock[5]  = 0.32;    /** K2O   */
        gv.bulk_rock[6]  = 2.42;    /** Na2O  */
        gv.bulk_rock[7]  = 1.13;    /** TiO2  */
        gv.bulk_rock[8]  = 0.37;    /** O     */
        gv.bulk_rock[9]  = 0.13;    /** MnO   */
        gv.bulk_rock[10] = 0.02;    /** Cr2O3 */
        gv.bulk_rock[11] = 6.66;    /** H2O   */
    }
    else if (gv.test == 1){ /* rough rhyolite */
        /* SiO2  Al2O3  CaO   MgO   FeO  K2O  Na2O  TiO2  O    MnO   Cr2O3 H2O */
        gv.bulk_rock[0]  = 76.58;   /** SiO2  */
        gv.bulk_rock[1]  = 7.93;    /** Al2O3 */
        gv.bulk_rock[2]  = 1.33;    /** CaO   */
        gv.bulk_rock[3]  = 0.46;    /** MgO   */
        gv.bulk_rock[4]  = 1.30;    /** FeO   */
        gv.bulk_rock[5]  = 2.97;    /** K2O   */
        gv.bulk_rock[6]  = 3.51;    /** Na2O  */
        gv.bulk_rock[7]  = 0.16;    /** TiO2  */
        gv.bulk_rock[8]  = 0.39;    /** O     */
        gv.bulk_rock[9]  = 0.04;    /** MnO   */
        gv.bulk_rock[10] = 0.004;   /** Cr2O3 */
        gv.bulk_rock[11] = 5.18;    /** H2O   */
    }
    else if (gv.test == 2){ /* rough rhyolite */
        /* SiO2  Al2O3  CaO   MgO   FeO  K2O  Na2O  TiO2  O    MnO   Cr2O3 H2O */
        gv.bulk_rock[0]  = 76.58;   /** SiO2  */
        gv.bulk_rock[1]  = 7.93;    /** Al2O3 */
        gv.bulk_rock[2]  = 1.33;    /** CaO   */
        gv.bulk_rock[3]  = 0.46;    /** MgO   */
        gv.bulk_rock[4]  = 1.30;    /** FeO   */
        gv.bulk_rock[5]  = 2.97;    /** K2O   */
        gv.bulk_rock[6]  = 3.51;    /** Na2O  */
        gv.bulk_rock[7]  = 0.16;    /** TiO2  */
        gv.bulk_rock[8]  = 0.0;     /** O     */
        gv.bulk_rock[9]  = 0.0;     /** MnO   */
        gv.bulk_rock[10] = 0.00;    /** Cr2O3 */
        gv.bulk_rock[11] = 5.18;    /** H2O   */
    }
    else{
        printf("Unknown test %i for 'gh' research group - please specify a different test! \n", gv.test);
        exit(EXIT_FAILURE);
    }
    return gv;
}
