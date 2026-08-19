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
#ifndef __GH_init_db_H_
#define __GH_init_db_H_

    #include "../MAGEMin.h"

    typedef struct gh_datasets {
        int     ds_version;
        int     n_ox;
        int     n_pp;
        int     n_ss;
        char    ox[13][20];
        char    PP[25][20];
        char    SS[16][20];

        int     verifyPC[16];
        int     n_SS_PC[16];
        double  SS_PC_stp[16];

        double  PC_df_add;
        double  solver_switch_T;
        double  min_melt_T;

        double  inner_PGE_ite;
        double  max_n_phase;
        double  max_g_phase;
        double  max_fac;

        double  merge_value;
        double  re_in_n;

        double  obj_tol;

    } gh_dataset;

    global_variable global_variable_GH_init(   global_variable      gv,
                                                bulk_info           *z_b    );

    global_variable get_bulk_gh( global_variable gv );
    global_variable get_bulk_pmelts_dataset( global_variable gv );

#endif
