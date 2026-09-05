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
/* Source: Pourteau et al. (2014), Contrib Mineral Petrol; Berman (1988)
   EOS formalism via the THERIAK JUN92.bs lineage. */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../MAGEMin.h"
#include "../initialize.h"
#include "BR_init_database.h"
#include "BR_PP_endmembers.h"

oxide_data oxide_info_br = {
    16,
    {"SiO2","Al2O3","MgO","FeO","K2O","Na2O","H2O","CaO","TiO2","O","MnO","Cr2O3","CO2","S","Cl","ecp"},
    {"Si","Al","Mg","Fe","K","Na","H","Ca","Ti","O","Mn","Cr","C","S","Cl","ecp"},
    {60.08,101.96,40.30,71.85,94.2,61.98,18.015,56.08,79.88,16.0,70.94,151.99,44.01,32.06,35.453,0.0},
    {3.0,5.0,2.0,2.0,3.0,3.0,3.0,2.0,3.0,1.0,2.0,5.0,3.0,1.0,1.0,0.0},
    {66.7736,108.653,40.3262,38.7162,69.1514,61.1729,69.5449,42.9947,70.3246,30.5827,40.1891,106.9795,62.8768,9.5557,33.2556,0.0},
    {2,3,1,1,1,1,1,1,2,1,1,3,2,0,0,0},
    {1,2,1,1,2,2,2,1,1,1,1,2,1,1,1,0}
};

/* Active subset of the full 85-entry BR_PP_endmembers.c table, curated for "po". */
br_dataset br_db = {
    1,
    10,
    14,
    17,
    {"SiO2","Al2O3","MgO","FeO","K2O","Na2O","H2O","CaO","TiO2","O"},
    {"cor","coe","q","ky","and","sill","mtlc","ftlc","law","glc","dsp","h2o","mt","hem"},
    {"ctd","car","chl","mica","talc","ilm","bt","ol","ep","opx","spl","stau","crd","grt","omph","amphx","fsp"},

    {1,1,1,1, 1,1,1,1,1,1,1,1,1,1,1,1,1, 0,0,0,0,0,0,0},
    {201,201,1001,1001, 201,201,201,201,201,201,201,201,201,201,861,861,861, 0,0,0,0,0,0,0},
    {0.001,0.001,0.01,0.01, 0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.005,0.005,0.005, 0.0,0.0,0.0,0.0,0.0,0.0,0.0},

    6.0,
    673.15,
    773.15,

    4,
    0.025,
    2.5,
    1.0,

    1e-1,
    1e-4,
    1e-6
};

global_variable global_variable_BR_init(   global_variable      gv,
                                            bulk_info           *z_b    ){
    int i, j;

    br_dataset db       = br_db;
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
        if (gv.mu_fix_idx[i] >= 0 && gv.mu_fix_idx[i] < gv.len_ox){
            sprintf(gv.PP_list[db.n_pp+i], "m%s", gv.ox[gv.mu_fix_idx[i]]);
        }
        else{
            sprintf(gv.PP_list[db.n_pp+i], "mu%d", i);
        }
    }

    gv.SS_list          = malloc ((gv.len_ss>0?gv.len_ss:1) * sizeof (char*) );
    gv.n_SS_PC          = malloc ((gv.len_ss>0?gv.len_ss:1) * sizeof (int)   );
    gv.verifyPC         = malloc ((gv.len_ss>0?gv.len_ss:1) * sizeof (int)   );
    gv.SS_PC_stp        = malloc ((gv.len_ss>0?gv.len_ss:1) * sizeof (double));
    for (i = 0; i < gv.len_ss; i++){
        gv.SS_list[i]   = malloc(20 * sizeof(char)             );
        strcpy(gv.SS_list[i],db.SS[i]);
        gv.verifyPC[i]  = db.verifyPC[i];
        gv.n_SS_PC[i]   = db.n_SS_PC[i];
        gv.SS_PC_stp[i] = db.SS_PC_stp[i];
    }
    gv.act_PP           = malloc ((gv.len_pp>0?gv.len_pp:1) * sizeof (int) );
    for (i = 0; i < gv.len_pp; i++){ gv.act_PP[i] = 1; };

    gv.n_min            = malloc ((gv.len_ss>0?gv.len_ss:1) * sizeof (int)  );
    gv.n_ss_ph          = malloc ((gv.len_ss>0?gv.len_ss:1) * sizeof (int)  );
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
    gv.n_solvi          = malloc ((gv.len_ss>0?gv.len_ss:1) * sizeof (int)   );

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
    gv.n_ss_array = malloc ((gv.len_ss>0?gv.len_ss:1) * sizeof(double));

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

    oxide_data ox_in = oxide_info_br;
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

    gv.n_em_db = 0;

    return gv;
}

global_variable get_bulk_br( global_variable gv) {
    if (gv.test != -1){
        if (gv.verbose == 1){
            printf("\n");
            printf("   - Minimization using in-built bulk-rock  : test %2d\n",gv.test);
        }
    }
    else{
        gv.test = 0;
    }

    if (gv.test == 0){
        /* SiO2  Al2O3  MgO   FeO  K2O  Na2O  H2O  CaO  TiO2  O */
        gv.bulk_rock[0]  = 64.578;  /** SiO2  */
        gv.bulk_rock[1]  = 13.651;  /** Al2O3 */
        gv.bulk_rock[2]  = 5.529;   /** MgO   */
        gv.bulk_rock[3]  = 8.025;   /** FeO   */
        gv.bulk_rock[4]  = 2.943;   /** K2O   */
        gv.bulk_rock[5]  = 2.0;     /** Na2O  */
        gv.bulk_rock[6]  = 13.5;    /** H2O   */
        gv.bulk_rock[7]  = 1.586;   /** CaO   */
        gv.bulk_rock[8]  = 0.907;   /** TiO2  */
        gv.bulk_rock[9]  = 0.5;     /** O     */
    }
    else if (gv.test == 1){
        /* SiO2  Al2O3  MgO   FeO  K2O  Na2O  H2O  CaO  TiO2  O */
        gv.bulk_rock[0]  = 64.578;  /** SiO2  */
        gv.bulk_rock[1]  = 13.651;  /** Al2O3 */
        gv.bulk_rock[2]  = 5.529;   /** MgO   */
        gv.bulk_rock[3]  = 8.025;   /** FeO   */
        gv.bulk_rock[4]  = 2.943;   /** K2O   */
        gv.bulk_rock[5]  = 2.0;     /** Na2O  */
        gv.bulk_rock[6]  = 13.5;    /** H2O   */
        gv.bulk_rock[7]  = 1.586;   /** CaO   */
        gv.bulk_rock[8]  = 0.907;   /** TiO2  */
        gv.bulk_rock[9]  = 0.0;     /** O     */
    }
    else{
        printf("Unknown test %i for 'br' research group - please specify a different test! \n", gv.test);
        exit(EXIT_FAILURE);
    }
    return gv;
}
